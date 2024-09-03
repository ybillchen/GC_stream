from collections.abc import Iterable

import numpy as np
from scipy.optimize import minimize

import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates.matrix_utilities import matrix_transpose, rotation_matrix

class StreamDataset(object):
    '''StreamDataset'''
    def __init__(self, ref, star_id, ra, dec, plx, mu_ra, mu_dec, mag, color, d, vr, evr, rvr, stream_id, prog):
        super().__init__()
        self.ref = ref
        self.star_id = star_id
        self.ra = ra
        self.dec = dec
        self.plx = plx
        self.mu_ra = mu_ra
        self.mu_dec = mu_dec
        self.mag = mag
        self.color = color
        self.d = d
        self.vr = vr
        self.evr = evr
        self.rvr = rvr
        self.prog = prog
        if isinstance(stream_id, Iterable):
            assert all(i == stream_id[0] for i in stream_id)
            self.stream_id = stream_id[0]
        else:
            self.stream_id = stream_id

        self.constuct_coord()

        self.c_GreatCircle = None
        self.phi1 = None
        self.phi2 = None
        self.phi2hat = None

    def constuct_coord(self):
        if self.ref == 'Ibata+21' or self.ref == 'mock':
            self.c = coord.SkyCoord(self.ra*u.deg, self.dec*u.deg, distance=self.d*u.kpc,
                pm_ra_cosdec=self.mu_ra*u.mas/u.yr, pm_dec=self.mu_dec*u.mas/u.yr, frame='icrs')
        elif self.ref == 'Ibata+24':
            self.c = coord.SkyCoord(self.ra*u.deg, self.dec*u.deg, 
                pm_ra_cosdec=self.mu_ra*u.mas/u.yr, pm_dec=self.mu_dec*u.mas/u.yr, frame='icrs')
        elif self.ref == 'Ibata+19':
            self.c = coord.SkyCoord(self.ra, self.dec, 
                frame='icrs', unit=(u.hourangle, u.deg))
        else:
            self.c = None
        self.c_prog = coord.SkyCoord(self.prog['RA']*u.deg, self.prog['DEC']*u.deg, distance=self.prog['Rsun']*u.kpc,
            pm_ra_cosdec=self.prog['mu_alpha']*u.mas/u.yr, pm_dec=self.prog['mu_delta']*u.mas/u.yr,
            radial_velocity=self.prog['<RV>']*u.km/u.s)

    def apply_mask(self, mask):
        assert len(mask) == len(self.ra)
        self.star_id = self.star_id[mask] if not self.star_id is None else None
        self.ra = self.ra[mask] if not self.ra is None else None
        self.dec = self.dec[mask] if not self.dec is None else None
        self.plx = self.plx[mask] if not self.plx is None else None
        self.mu_ra = self.mu_ra[mask] if not self.mu_ra is None else None
        self.mu_dec = self.mu_dec[mask] if not self.mu_dec is None else None
        self.mag = self.mag[mask] if not self.mag is None else None
        self.color = self.color[mask] if not self.color is None else None
        self.d = self.d[mask] if not self.d is None else None
        self.vr = self.vr[mask] if not self.vr is None else None
        self.evr = self.evr[mask] if not self.evr is None else None
        self.rvr = self.rvr[mask] if not self.rvr is None else None
        self.c = self.c[mask] if not self.c is None else None
        self.phi1 = self.phi1[mask] if not self.phi1 is None else None
        self.phi2 = self.phi2[mask] if not self.phi2 is None else None
        self.phi2hat = self.phi2hat[mask] if not self.phi2hat is None else None
        self.c_GreatCircle = self.c_GreatCircle[mask] if not self.c_GreatCircle is None else None

    def cut_mag(self, mag_limit):
        self.apply_mask(self.mag <= mag_limit)

    def construct_greatcircle(self):

        record = np.inf
        for guess1 in range(-180,181,30):
            for guess2 in range(-180,181,30):
                o = minimize(rotate1, (guess1,guess2), args=(self.c.ra.to_value(u.rad), self.c.dec.to_value(u.rad)))
                if o.fun < record:
                    alpha = o.x[0]
                    beta = o.x[1]
                    record = o.fun
        # print(np.sqrt(record)*180/np.pi)
        o = minimize(rotate2, (0), args=(alpha, beta, self.c.ra.to_value(u.rad), self.c.dec.to_value(u.rad)))
        # print(np.sqrt(o.fun)*180/np.pi)
        gamma = o.x[0]
        # print(alpha, beta, gamma)
        self.mat = rotation_matrix(gamma*u.deg, 'z') @ rotation_matrix(beta*u.deg, 'x') @ rotation_matrix(alpha*u.deg, 'z')

        class GreatCircle(coord.BaseCoordinateFrame):
            default_representation = coord.SphericalRepresentation
            default_differential = coord.SphericalCosLatDifferential

            frame_specific_representation_info = {
                coord.SphericalRepresentation: [
                    coord.RepresentationMapping('lon', 'phi1'),
                    coord.RepresentationMapping('lat', 'phi2'),
                    coord.RepresentationMapping('distance', 'distance'),
                ]
            }

        @coord.frame_transform_graph.transform(coord.StaticMatrixTransform, coord.ICRS, GreatCircle)
        def icrs_to_GreatCircle():
            return self.mat
        
        @coord.frame_transform_graph.transform(coord.StaticMatrixTransform, GreatCircle, coord.ICRS)
        def GreatCircle_to_icrs():
            return matrix_transpose(self.mat)
        
        self.GreatCircle = GreatCircle

    def construct_greatcircle_coord(self, GreatCircle=None):
        if GreatCircle is None:
            GreatCircle = self.GreatCircle
        self.c_GreatCircle = self.c.transform_to(GreatCircle())
        self.phi1 = self.c_GreatCircle.phi1.wrap_at(180*u.deg).to_value(u.deg)
        self.phi2 = self.c_GreatCircle.phi2.to_value(u.deg)
        self.c_prog_GreatCircle = self.c_prog.transform_to(GreatCircle())
        self.prog_phi1 = self.c_prog_GreatCircle.phi1.wrap_at(180*u.deg).to_value(u.deg)
        self.prog_phi2 = self.c_prog_GreatCircle.phi2.to_value(u.deg)

    def construct_polynomial_coord(self, p=None, deg=2):
        if p is None:
            p = np.polyfit(self.phi1, self.phi2, deg=deg)
        self.p = p
        self.fit_poly = np.poly1d(p)
        self.phi2hat = self.phi2 - self.fit_poly(self.phi1)
        self.prog_phi2hat = self.prog_phi2 - self.fit_poly(self.prog_phi1)

    def set_release_time(self, release_time):
        assert self.ref == 'mock'
        self.release_time = release_time


def rotate1(angles, ra, dec):
    alpha = angles[0]
    beta = angles[1]
    mat = rotation_matrix(beta*u.deg, 'x') @ rotation_matrix(alpha*u.deg, 'z')
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    pos = np.row_stack((x,y,z))
    pos_new = mat @ pos
    x = pos_new[0]
    y = pos_new[1]
    z = pos_new[2]

    phi1 = np.arctan2(y, x)
    phi2 = np.arcsin(z)

    return np.mean(phi2**2)

def rotate2(angles, alpha, beta, ra, dec):
    gamma = angles[0]
    mat = rotation_matrix(gamma*u.deg, 'z') @ rotation_matrix(beta*u.deg, 'x') @ rotation_matrix(alpha*u.deg, 'z')
    x = np.cos(dec) * np.cos(ra)
    y = np.cos(dec) * np.sin(ra)
    z = np.sin(dec)

    pos = np.row_stack((x,y,z))
    pos_new = mat @ pos
    x = pos_new[0]
    y = pos_new[1]
    z = pos_new[2]

    phi1 = np.arctan2(y, x)
    phi2 = np.arcsin(z)

    return np.mean(phi1**2)
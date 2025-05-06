import numpy as np
import agama
from scipy.signal import argrelmin, argrelmax
import astropy.units as u
import astropy.coordinates as coord

agama.setUnits(mass=1, length=1, velocity=1)
timeunitgyr = agama.getUnits()['time'].value/1000

class StreamOutput:
    
    def __init__(self,time_orbit,orbit_sat,xv_stream,n_release):
        self.time_orbit = time_orbit
        self.orbit_sat = orbit_sat
        self.xv_stream = xv_stream
        self.n_release = n_release


    def shift_stream_coordinates(self, xv_stream, sun_velocity=(11.1, 232.24, 7.25), return_R=False):
        """
        Shifts the stream coordinates to a phi_1, phi_2 frame using Neil's routine. Also calculates the 
        proper motions (note that this assumes a solar position and velocity).
        """


        
        #get positions
        pos = xv_stream[:,:3]*np.array([1,-1,-1])

        # Make covariance matrix
        covMat = pos.T.dot(pos)

        # diagonalize the covariance matrix
        VarOnPrincipalAxis, principalAxes = np.linalg.eig(covMat) 

        # sort the principal axes in decreasing order
        ix = np.argsort(-VarOnPrincipalAxis)
        principalAxes = principalAxes[:,ix]
        VarOnPrincipalAxis = VarOnPrincipalAxis[ix]

        # Make sure the eigenvectors form a right-handed coordinate system
        if np.linalg.det(principalAxes) < 0:
            principalAxes[:,0] *= -1

        # Check that the first eigenvector is in the correct direction
        if np.sum(np.sum(principalAxes[:,0]*pos,axis=1))<0:
            principalAxes[:,0] *= -1
            principalAxes[:,1] *= -1

        # Make rotation matrix
        R = principalAxes.T

        # Rotate data
        pos_rot = R.dot(pos.T).T

        # Get PHI1 and PHI2
        phi1_temp = np.arctan2(pos_rot[:,1],pos_rot[:,0])
        phi2_temp = np.arctan2(pos_rot[:,2],np.sqrt(pos_rot[:,0]**2 + pos_rot[:,1]**2))

        phi1 = phi1_temp ## in radians
        phi2 = phi2_temp ## in radians

        ## Now, calculate the proper motions
        # Get velocities
        vel = xv_stream[:, 3:6] - np.array(sun_velocity)  # Subtract solar motion

        # Rotate velocities
        vel_rot = R @ vel.T
        vx_rot, vy_rot, vz_rot = vel_rot

        # Convert to proper motions
        dist = np.linalg.norm(pos, axis=1)  # Distance in kpc
        coords_gal = coord.SkyCoord(x=pos[:,0]*u.kpc, y=pos[:,1]*u.kpc, z=pos[:,2]*u.kpc, frame='galactocentric')
        coords_icrs = coords_gal.transform_to(coord.ICRS)
        dist = coords_icrs.distance.kpc
        #pm_phi1 = (-vy_rot * np.cos(phi1) + vx_rot * np.sin(phi1)) / (4.74047 * dist)  # mas/yr
        pm_phi1 = vx_rot*np.cos(phi2) / (4.74047 * dist)  # mas/yr
        pm_phi2 = vz_rot / (4.74047 * dist)  # mas/yr
   
    
    

       
        if return_R:
            return (phi1*180/np.pi)*u.deg, (phi2*180/np.pi)*u.deg, pm_phi1, pm_phi2, R
        else:
            return (phi1*180/np.pi)*u.deg, (phi2*180/np.pi)*u.deg, pm_phi1, pm_phi2





class ProgenitorOrbit:

    def __init__(self, time_orbit, orbit_sat):
        self.t = time_orbit
        self.orbit_sat = orbit_sat

        self.x = self.orbit_sat.T[0]
        self.y = self.orbit_sat.T[1]
        self.z = self.orbit_sat.T[2]
        self.r = np.sqrt(self.x**2 + self.y**2 + self.z**2)

        self.vx = self.orbit_sat.T[3]
        self.vy = self.orbit_sat.T[4]
        self.vz = self.orbit_sat.T[5]

    
    def pericenter(self):
        pericenters = self.r[argrelmin(self.r)]
        pericenter_times = self.t[argrelmin(self.r)]

        return pericenters, pericenter_times
    

    def apocenter(self):
        apocenters = self.r[argrelmax(self.r)]
        apocenter_times = self.t[argrelmax(self.r)]

        return apocenters, apocenter_times
    

    def angular_momentum(self):
        r_vec = np.array([self.x, self.y, self.z]).T
        v_vec = np.array([self.vx, self.vy, self.vz]).T

        L = np.cross(r_vec, v_vec)*u.kpc*u.km/u.s
        return L.to(u.kpc**2/u.Myr)
    
    def energy(self, potential):
        vx = self.vx*u.km/u.s
        vy = self.vy*u.km/u.s
        vz = self.vz*u.km/u.s
        v2 = (vx**2 + vy**2 + vz**2).to(u.kpc**2/u.Myr**2)

        phi = (potential.potential(np.array([self.x, self.y, self.z]).T)/timeunitgyr**2)*(u.kpc**2/u.Gyr**2)

        E = (1/2)*(v2) + phi.to(u.kpc**2/u.Myr**2)
        return E

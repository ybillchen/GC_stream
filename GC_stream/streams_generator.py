"""
BSD 3-Clause License
Copyright (c) 2024 Yingtian Chen
All rights reserved.
"""

import numpy as np

gravG = 4.30091727003628e-6

class Orbit(object):
    """Orbit class"""
    def __init__(self, posvel, times):
        """
        posvel: (N,6) ndarray
        """
        super(Orbit, self).__init__()
        sorted_arg = np.argsort(times)
        self.posvel = posvel[sorted_arg]
        self.times = times[sorted_arg]

class MassHistory(object):
    """Mass history class"""
    def __init__(self, masses, times):
        super(MassHistory, self).__init__()
        sorted_arg = np.argsort(times)
        self.masses = masses[sorted_arg]
        self.times = times[sorted_arg]

class BaseStreamsGenerator(object):
    """Base class for streams generator"""

    def __init__(self, pot):
        super(BaseStreamsGenerator, self).__init__()
        self.pot = pot
        self.name = 'Base'

    def get_f(self, rsat, time):
        return 3 - self.pot.dlnMdlnr(rsat, time)

    def rtid(self, rsat, msat, time):
        f = self.get_f(rsat, time)
        return rsat*(msat/(f*self.pot.Menclose(rsat, time)))**(1/3)

    def sample(self, psat, vsat, msat, time, rng):
        """
        Sample one spray particle
        """
        raise Exception("Undefined method! \
            Must define in actual generator class")

    def sample_stream(self, t_begin, t_end, 
        orbit, mass_history, dm, rng):
        """
        Sample entire stream from t_begin to t_end,
        using orbit and mass history
        """

        t_begin = np.max((t_begin, mass_history.times[0]))
        t_end = np.min((t_end, mass_history.times[-1]))

        # This may not be super accurate, since the last particle is 
        # actually less than dm. But doesn't matter if dm is tiny.
        mass_sample = np.arange(mass_history.masses[0]-dm, 
            mass_history.masses[-1]-dm, -dm)

        # Mass is always decreasing, so should be flipped
        time_sample = np.interp(mass_sample, 
            mass_history.masses[::-1], mass_history.times[::-1])

        mask = (time_sample>=t_begin) & (time_sample<=t_end)
        mass_sample = mass_sample[mask]
        time_sample = time_sample[mask]

        posvel_sample = []
        for i in range(6):
            posvel_sample.append(np.interp(time_sample, 
                orbit.times, orbit.posvel[:,i]))
        posvel_sample = np.column_stack((posvel_sample))

        posvel_ej = np.full_like(posvel_sample, fill_value=np.nan)

        for i in range(len(time_sample)):
            pi, vi = self.sample(
                posvel_sample[i,:3], posvel_sample[i,3:], 
                mass_sample[i], time_sample[i], rng)
            posvel_ej[i,:3] = pi
            posvel_ej[i,3:] = vi

        return time_sample, posvel_ej

class F15StreamsGenerator(BaseStreamsGenerator):
    """Fardal et al. (2015)"""

    def __init__(self, pot, Rapo, Rperi, ft=1, gala=False):
        super(F15StreamsGenerator, self).__init__(pot)
        self.name = 'F15'
        self.Rapo = Rapo
        self.Rperi = Rperi
        self.ft = ft
        self.gala = bool=(gala)

    def sample(self, psat, vsat, msat, time, rng):
        # calculate the rotation matrix
        # after rotation, pnew is on x-axis, vnew has no z-component
        dir_x = psat
        dir_x = dir_x / np.sqrt(np.sum(dir_x**2))
        dir_z = np.cross(psat, vsat)
        dir_z = dir_z / np.sqrt(np.sum(dir_z**2))
        dir_y = np.cross(dir_z, dir_x)
        dir_y = dir_y / np.sqrt(np.sum(dir_y**2))

        # [x']   [dir_x]   [x]
        # [y'] = [dir_y] @ [y]
        # [z']   [dir_z]   [z]
        R = np.array([dir_x, dir_y, dir_z]) # rotation matrix
        Rinv = R.T # inverse of an orthogonal matrix is its transpose

        # rotation
        pnew = R @ psat
        vnew = R @ vsat

        direction = 1 - 2*rng.integers(2) # 1:trailing -1:leading 
        rsat = pnew[0]
        rtid = self.rtid(rsat, msat, time) * direction

        # calculate the ejection position and velocity
        # Fardal et al. (2015)
        if self.gala:
            sig = 0.5
        else:
            Omega_peri = self.pot.Vcirc(self.Rperi, time) / self.Rperi # units doesn't matter
            Omega_apo = self.pot.Vcirc(self.Rapo, time) / self.Rapo # units doesn't matter
            ga_peri = Omega_peri**2 * self.get_f(self.Rperi, time)
            ga_apo = Omega_apo**2 * self.get_f(self.Rapo, time)
            Racc = ga_peri / ga_apo
            assert Racc >= 1
            sig = np.min((0.4, 0.15*self.ft**2*Racc**(2/3)))
        kr = rng.normal(2, sig)
        kphi = 0
        kvr = 0
        kvt = rng.normal(0.3, sig)
        kz = rng.normal(0, 0.5)
        kvz = rng.normal(0, 0.5)

        rej = rsat + kr*rtid
        phiej = kphi*rtid/rsat
        vrej = (1+kvr)*vnew[0]
        vtej = vnew[1] + kvt*self.pot.Vcirc(rsat, time)*rtid/rsat
        zej = kz*rtid/rsat
        vzej = kvz*self.pot.Vcirc(rsat, time)*rtid/rsat

        pejnew = np.array([
            rej*np.cos(phiej), 
            rej*np.sin(phiej), 
            zej])

        vejnew = np.array([
            vrej*np.cos(phiej)-vtej*np.sin(phiej), 
            vrej*np.sin(phiej)+vtej*np.cos(phiej),
            vzej])
        
        # rotation back to the original coordinates
        pej = Rinv @ pejnew
        vej = Rinv @ vejnew

        return pej, vej


class R24StreamsGenerator(BaseStreamsGenerator):
    """Roberts et al. (2024, QSG-PS)"""

    def __init__(self, pot, fe=1.5, eps=0.57):
        super(R24StreamsGenerator, self).__init__(pot)
        self.name = 'R24'
        self.fe = fe
        self.eps = eps

    def sample(self, psat, vsat, msat, time, rng):
        # calculate the rotation matrix
        # after rotation, pnew is on x-axis, vnew has no z-component
        dir_x = psat
        dir_x = dir_x / np.sqrt(np.sum(dir_x**2))
        dir_z = np.cross(psat, vsat)
        dir_z = dir_z / np.sqrt(np.sum(dir_z**2))
        dir_y = np.cross(dir_z, dir_x)
        dir_y = dir_y / np.sqrt(np.sum(dir_y**2))

        # [x']   [dir_x]   [x]
        # [y'] = [dir_y] @ [y]
        # [z']   [dir_z]   [z]
        R = np.array([dir_x, dir_y, dir_z]) # rotation matrix
        Rinv = R.T # inverse of an orthogonal matrix is its transpose

        # rotation
        pnew = R @ psat
        vnew = R @ vsat

        direction = 1 - 2*rng.integers(2) # 1:trailing -1:leading 
        rJ = self.rtid(pnew[0], msat, time)
        pejnew = np.array([pnew[0] + self.fe*rJ*direction, 0, 0])

        a = 1.305 * 0.15 * rJ
        disp = np.sqrt(gravG*msat/(6*a)) * (1 + (self.fe*rJ/a)**2)**(-0.25) # Eq. (10)

        vejnew = np.copy(vnew)
        vejnew[0] += rng.normal(0.0, disp)
        vejnew[1] += rng.normal(
            self.eps*self.pot.Vcirc(pnew[0], time)*self.fe*rJ*direction/pnew[0], disp)
        vejnew[2] += rng.normal(0.0, disp)
        
        # rotation back to the original coordinates
        pej = Rinv @ pejnew
        vej = Rinv @ vejnew

        return pej, vej

class SphStreamsGenerator(BaseStreamsGenerator):
    """Streams generator using spherical coordinates"""

    def __init__(self, pot, mean, cov, orbit_dependent):
        super(SphStreamsGenerator, self).__init__(pot)
        self.name = 'Sph'

        # variables in the multivariate normal distributions:
        # 1. Dr_rtid : Dr/rtid
        # 2. phi     : position azimuth (arcdeg)
        # 3. theta   : position latitude (arcdeg)
        # 4. Dv_vesc : Dv/vesc
        # 5. alpha   : velocity azimuth (arcdeg)
        # 6. beta    : velocity latitude (arcdeg)
        self.mean = mean # array of mean values
        self.cov = cov # covariance matrix
        self.orbit_dependent = bool(orbit_dependent)

    def sample(self, psat, vsat, msat, time, rng):
        # calculate the rotation matrix
        # after rotation, pnew is on x-axis, vnew has no z-component
        dir_x = psat
        dir_x = dir_x / np.sqrt(np.sum(dir_x**2))
        dir_z = np.cross(psat, vsat)
        dir_z = dir_z / np.sqrt(np.sum(dir_z**2))
        dir_y = np.cross(dir_z, dir_x)
        dir_y = dir_y / np.sqrt(np.sum(dir_y**2))

        # [x']   [dir_x]   [x]
        # [y'] = [dir_y] @ [y]
        # [z']   [dir_z]   [z]
        R = np.array([dir_x, dir_y, dir_z]) # rotation matrix
        Rinv = R.T # inverse of an orthogonal matrix is its transpose

        # rotation
        pnew = R @ psat
        vnew = R @ vsat

        direction = 1 - 2*rng.integers(2) # 1:trailing -1:leading 
        rsat = pnew[0]
        rtid = self.rtid(rsat, msat, time)

        # calculate the ejection position and velocity
        [Dr_rtid, phi, theta, Dv_vesc, alpha, beta] = \
            rng.multivariate_normal(self.mean, self.cov)

        if self.orbit_dependent:
            vnew_r = vnew[0]
            vnew_t = vnew[1]
            vc = self.pot.Vcirc(rsat, time)
            Dr_rtid += -0.5*vnew_r/vc
            phi += -30*(-1+vnew_t/vc)

        Dr = Dr_rtid * rtid
        vesc = np.sqrt(2*gravG*msat/Dr) # escape velocity
        Dv = Dv_vesc * vesc

        if direction < 0: # leading
            phi += 180
            alpha += 180

        # convert degrees to radians
        phi *= (np.pi/180)
        theta *= (np.pi/180)
        alpha *= (np.pi/180)
        beta *= (np.pi/180)

        pejnew = pnew + np.array([
            Dr*np.cos(theta)*np.cos(phi),
            Dr*np.cos(theta)*np.sin(phi),
            Dr*np.sin(theta)])

        vejnew = vnew + np.array([
            Dv*np.cos(beta)*np.cos(alpha),
            Dv*np.cos(beta)*np.sin(alpha),
            Dv*np.sin(beta)])
        
        # rotation back to the original coordinates
        pej = Rinv @ pejnew
        vej = Rinv @ vejnew

        return pej, vej


class SphModStreamsGenerator(BaseStreamsGenerator):
    """Streams generator using spherical coordinates (modified)"""

    def __init__(self, pot, mean, cov, orbit_dependent):
        super(SphModStreamsGenerator, self).__init__(pot)
        self.name = 'SphMod'

        # variables in the multivariate normal distributions:
        # 1. Dr_rtid : Dr/rtid
        # 2. phi     : position azimuth (arcdeg)
        # 3. theta   : position latitude (arcdeg)
        # 4. Dv_vesc : Dv/vesc
        # 5. alpha   : velocity azimuth (arcdeg)
        # 6. beta    : velocity latitude (arcdeg)
        self.mean = mean # array of mean values
        self.cov = cov # covariance matrix
        self.orbit_dependent = bool(orbit_dependent)

    def sample(self, psat, vsat, msat, time, rng):
        # calculate the rotation matrix
        # after rotation, pnew is on x-axis, vnew has no z-component
        dir_x = psat
        dir_x = dir_x / np.sqrt(np.sum(dir_x**2))
        dir_z = np.cross(psat, vsat)
        dir_z = dir_z / np.sqrt(np.sum(dir_z**2))
        dir_y = np.cross(dir_z, dir_x)
        dir_y = dir_y / np.sqrt(np.sum(dir_y**2))

        # [x']   [dir_x]   [x]
        # [y'] = [dir_y] @ [y]
        # [z']   [dir_z]   [z]
        R = np.array([dir_x, dir_y, dir_z]) # rotation matrix
        Rinv = R.T # inverse of an orthogonal matrix is its transpose

        # rotation
        pnew = R @ psat
        vnew = R @ vsat

        direction = 1 - 2*rng.integers(2) # 1:trailing -1:leading 
        rsat = pnew[0]
        rtid = self.rtid(rsat, msat, time)

        # calculate the ejection position and velocity
        [Dr_rtid, phi, theta, Dv_vesc, alpha, beta] = \
            rng.multivariate_normal(self.mean, self.cov)

        if self.orbit_dependent:
            vnew_r = vnew[0]
            vnew_t = vnew[1]
            vc = self.pot.Vcirc(rsat, time)
            Dr_rtid += -0.5*vnew_r/vc
            phi += -30*(-1+vnew_t/vc)

        Dr = Dr_rtid * rtid
        vesc = np.sqrt(2*gravG*msat/Dr) # escape velocity
        Dv = Dv_vesc * vesc

        if direction < 0: # leading
            phi += 180
            alpha += 180

        theta /= np.sqrt(self.cov[2,2]) 
        theta *= 12 if Dr_rtid > 1 else 22

        # convert degrees to radians
        phi *= (np.pi/180)
        theta *= (np.pi/180)
        alpha *= (np.pi/180)
        beta *= (np.pi/180)

        pejnew = pnew + np.array([
            Dr*np.cos(theta)*np.cos(phi),
            Dr*np.cos(theta)*np.sin(phi),
            Dr*np.sin(theta)])

        vejnew = vnew + np.array([
            Dv*np.cos(beta)*np.cos(alpha),
            Dv*np.cos(beta)*np.sin(alpha),
            Dv*np.sin(beta)])
        
        # rotation back to the original coordinates
        pej = Rinv @ pejnew
        vej = Rinv @ vejnew

        return pej, vej

class C24StreamsGenerator(SphStreamsGenerator):
    """Chen et al. (2024)"""
    def __init__(self, pot,
        mean = np.array([1.5, -30.0, 0.0, 1.0, 20.0, 0.0]),
        std = np.array([0.4, 23.0, 12.0, 0.0, 20.0, 22.0]),
        R_04 = -0.625,
        orbit_dependent = True):
        
        cov = std.T @ std
        cov_04 = R_04*std[0]*std[4]
        cov[0,4] = cov_04
        cov[4,0] = cov_04
        
        super(C24StreamsGenerator, self).__init__(pot, mean, cov, orbit_dependent)
            

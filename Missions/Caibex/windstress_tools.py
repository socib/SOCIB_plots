import numpy as np
import numexpr as ne

def cdntc(sp, Ta):
        '''
        cd, u10 = cdntc(sp, z , Ta) computes the neutral drag coefficient and
        wind speed at 10m given the wind speed and air temperature at height z
        following Smith (1988), J. Geophys. Res., 93, 311-326.
    
        INPUT:   sp - wind speed  [m/s]
             self.Z - measurement height [m]
            self.Ta - air temperature (optional)  [C]
    
        OUTPUT:  cd - neutral drag coefficient at 10m
                  u10 - wind speed at 10m  [m/s]
        ''' 
        g = 9.8        # acceleration due to gravity [m/s^2]
        sigmaSB     = 5.6697e-8  # Stefan-Boltzmann constant [W/m^2/K^4]
        eps_air     = 0.62197    # molecular weight ratio (water/air)
        gas_const_R = 287.04     # gas constant for dry air [J/kg/K]
        CtoK        = 273.16     # conversion factor for [C] to [K]

        kappa          = 0.4     # von Karman's constant
        Charnock_alpha = 0.011   # Charnock constant (for determining roughness length
        R_roughness   = 0.11     # limiting roughness Reynolds # for aerodynamically

                              
        # ------ defaults suitable for boundary-layer studies
        cp            = 1004.7    # heat capacity of air [J/kg/K]
        rho_air       = 1.22      # air density (when required as constant) [kg/m^2]
        Ta_default    = 20        # default air temperature [C]

        z = 10.
        
        # Iteration endpoint
        TOL=.00001
        
        # Get air viscosity
        visc = viscair(Ta)
        #print '----------visc',visc
        
        # Remove any sp==0 to prevent division by zero
        np.place(sp, sp == 0, 0.1)
        
        # Initial guess
        ustaro = np.zeros(sp.shape)
        ustarn = np.copy(sp)
        ustarn = ne.evaluate('ustarn * 0.036')
        #print 'ustarn1', ustarn
        
        # Iterate to find z0 and ustar
        while np.any(ne.evaluate('abs(ustarn - ustaro)') > TOL):
            ustaro = np.copy(ustarn)
            z0 = ne.evaluate('Charnock_alpha * ustaro**2 / g + \
                              R_roughness * visc / ustaro')
            ustarn = ne.evaluate('sp * (kappa / log(z / z0))')
     
        #print 'ustarn2', ustarn
        sqrcd = ne.evaluate('kappa / log(10. / z0)')
        cd = ne.evaluate('sqrcd**2') # np.power(sqrcd, 2)
        #u10 = ustarn / sqrcd
        u10 = ne.evaluate('ustarn / sqrcd')
        #print cd, u10
        return cd, u10
        
        
def stresstc(u, v, Ta, rho_air):
        '''
        taux, tauy = stresstc(u, v) computes the neutral wind stress given the
        wind speed and air temperature at height z following Smith (1988),
        J. Geophys. Res., 93, 311-326

        INPUT:  u, v  - wind speed components   [m/s]
                z     - measurement height  [m]
                Ta    - air temperature [C]
                rho_air  - air density  [kg/m^3]

        OUTPUT: taux, tauy - wind stress components  [N/m^2]
        Adapted from http://woodshole.er.usgs.gov/operations/sea-mat/air_sea-html/as_consts.html
        '''
        # NOTE: windspd is calculated in __main__ before call to self.stresstc, however if
        # windspd_fill_mask == True there can be issues at the coast in taux, tauy
        # Hence, we calculate it again without mask filling
        sp = np.hypot(u, v) # this line commented Nov2014 cos redundant from _get_wspd
        cd, u10 = cdntc(sp, Ta)
       
        taux = cd*rho_air*u*u
        tauy = cd*rho_air*v*v
        return taux, tauy



def air_dens(Ta, RH, Pa, Q):
        '''
        rhoa=AIR_DENS(Ta,RH,Pa) computes the density of moist air.
        Air pressure is optional.

        INPUT:   Ta  - air temperature Ta  [C]
                 RH  - relative humidity  [%]
                 Pa  - air pressure (optional) [mb]
                 Q   - Q air
        OUTPUT:  rhoa - air density  [kg/m^3]
        '''
        #if Pa is not None:
            #self.Pa = Pa # pressure in mb
        
        #if Q is None:
            #Q = self.qsat(Ta, Pa)
            #Q *= 0.01
            #Q *= RH # specific humidity of air [kg/kg]
        
        o_six_one = 1. / self.EPS_AIR - 1. # 0.61 (moisture correction for temp.)
        o_six_one *= Q
        o_six_one += 1.
        Tv = np.copy(Ta)
        Tv += self.C2K
        Tv *= o_six_one # air virtual temperature
        Tv *= self.GAS_CONST_R
        rhoa = np.copy(Pa)
        rhoa *= 100.
        rhoa /= Tv # air density [kg/m^3]
        return rhoa

def viscair(Ta):
        '''
        vis=VISCAIR(Ta) computes the kinematic viscosity of dry air as a
        function of air temperature following Andreas (1989), CRREL Report
        89-11.

        INPUT: Ta - air temperature [C]
            1.326e-5*(1 + 6.542e-3*Ta + 8.301e-6*Ta.^2 - 4.84e-9*Ta.^3)
        OUTPUT: vis - air viscosity [m^2/s]
        '''
        Ta = np.asarray(Ta)
        #va = 6.542e-3 * Ta
        #va += 8.301e-6 * Ta**2
        #va -= 4.84e-9 * Ta**3
        #va += 1
        #va *= 1.326e-5
        va = ne.evaluate('1.326e-5 * (1 + 6.542e-3 * Ta + 8.301e-6 * \
                          Ta**2 - 4.84e-9 * Ta**3)')
        return va

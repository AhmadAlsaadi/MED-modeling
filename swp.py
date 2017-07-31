# -*- coding: utf-8 -*-
import numpy as np
# seawater density
def Density(Temp=25, Salinity=0):
    '''this function calculate the density in Kg/m3 of water at different temperature (C)
    and different salinity in weight percent. the function is valid in the following ranges
    Temperature range 10 to 180 C
    Salinity range 0 to 16 wt%
    reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-526
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 27/5/2012.'''
    Sal=np.array(Salinity)
    T=np.array(Temp)
    if Sal.ndim==0:
        Sal=np.array([Sal])
    WD = [10. ** 3. * ((4.032219 * (0.5) + 0.115313 * ((2 * S * 10. - 150.) / 150.) + 3.26 * 10. ** -4. * \
        (2. * ((2. * S * 10. - 150.) / 150.) ** 2. - 1.)) * (0.5) + (-0.108199 * (0.5) + 1.571 * 10. ** -3. * \
        ((2. * S * 10. - 150.) / 150.) - 4.23 * 10. ** -4. * (2. * ((2. * S * 10. - 150.) / 150.) ** 2 - 1)) \
        * ((2. * T - 200.) / 160.) + (-0.012247 * (0.5) + 1.74 * 10. ** -3. * ((2. * S * 10. - 150.) / 150.) \
        - 9. * 10. ** -6. * (2. * ((2. * S * 10. - 150.) / 150.) ** 2. - 1.)) * (2. * ((2. * T - 200.) / 160.) **\
        2. - 1.)+ (6.92 * 10. ** -4. * (0.5) - 8.7 * 10. ** -5. * ((2. * S * 10. - 150.) / 150.) - 5.3 * 10. ** -5. \
        * (2. * ((2. * S * 10. - 150.) / 150.) ** 2. - 1.)) * (4. * ((2. * T - 200.) / 160.) ** 3. - 3. * \
        ((2. * T - 200.) / 160.))) for S in Sal]
    return WD

# specific heat
def SpecificHeat(Temp=25, Salinity=0):
    '''this function calculate the specific heat in J/kg.C of water at different temperature (C)
    and different salinity in weight percent. the function is valid in the following ranges
    Temperature range 20 to 180 C
    Salinity range 2 to 16 wt%
    reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-528
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 27/5/2012.'''
    
    Sal=np.array(Salinity)
    T=np.array(Temp)
    if Sal.ndim==0:
        Sal=np.array([Sal])
        
    SpHeat = [(((4206.8 - 0.0066197 * S * 10000. + 1.2288e-06 * 10. ** -2. * (S * 10000.) ** 2.) \
                + (-1.1262 + 0.0054178 * 10. ** -2. * S * 10000. - 2.2719e-06 * 10. ** -4. * (S * 10000.) ** 2.) \
                * T + (1.2026 * 10. ** -2. - 0.0053566 * 10. ** -4. * S * 10000. + 1.8906e-06 * 10. ** -6. * \
                (S * 10000.) ** 2.) * T ** 2. + (6.8777 * 10. ** -7. + 0.001517 * 10. ** -6. * S * \
                10000. - 4.4268e-06 * 10. ** -9. * (S * 10000.) ** 2.) * T ** 3.) * 10. ** -3.) * 1000. for S in Sal]
    return SpHeat

# Seawater Dynamic viscosity
def DynamicViscosity(Temp=25, Salinity=0):
    '''this function calculate the dynamic viscosity in Kg/m.s of water at different temperature (C)
    and different salinity in weight percent. the function is valid in the following ranges
    Temperature range 10 to 180 C
    Salinity range 0 to 130 g/kg
    reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-530
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 27/5/2012.'''
    
    Sal=np.array(Salinity)
    T=np.array(Temp)
    if Sal.ndim==0:
        Sal=np.array([Sal])
    
    DynVisc = [(np.e**(-3.79418 + 604.129 / (139.18 + T))) * (1. + (1.474 * 10. ** -3. + 1.5 * 10. ** -5. * T - 3.927 * 10. ** -8. * T ** 2.) * \
                    (S * 10.) + (1.0734 * 10. ** -5. - 8.5 * 10. ** -8. * T + 2.23 * 10. ** -10. * T ** 2.) * (S * 10.) ** 2.) * 0.001 for S in Sal]
    return DynVisc

# Saturated water vapor viscosity
def SatVapViscosity(Temp=25):
    '''this function calculate the dynamic viscosity of Saturated water vapor in Kg/m.s at different temperature (C)
    'and different salinity in weight percent. the function is valid in the following ranges
    'Temperature range 10 to 180 C
    'reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-554
    'creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 23/3/2013.'''
    SatVapVis = np.power(np.e,(-3.609417664 + 275.928958 / (-227.0446083 - 0.896081232 * Temp - 0.002291383 * Temp ** 2))) * 0.001
    return SatVapVis

# Specific volume of water vapor
def SpecificVolume(Temp=25):
    '''this function calculate the specific volume of water vapor in M3/kg at different temperature in C .
    'reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-548
    'creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 4/22/2013.'''
    SpVol = 0.003172222 * (647.286 / (273.14 + Temp) - 1) * np.power(np.e,(83.63213098 - 0.668265339 * (Temp + 273.15) ** 1. \
                    + 0.002495964 * (Temp + 273.15) ** 2. - 5.04185e-06 * (Temp + 273.15) ** 3. + 5.34205e-09 * \
                    (Temp + 273.15) ** 4. - 2.3279e-12 * (Temp + 273.15) ** 5.))
    return SpVol

# Enthalpy of saturated liquid water
def SatLiquidEnthalpy(Temp=25):
    '''this function calculate the saturated water liquid enthalpy in kj/kg of water at different temperature (C)
    the function is valid in the following range
    Temperature range 5 to 200 C
    reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-538
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 1/6/2012.'''
    T=np.array(Temp)
    SatLiqEnthalpy = - 0.033635409 + 4.207557011*T - 6.200339*10.**-4.*T**2 + 4.459374*10.**-6.*T**3
    return SatLiqEnthalpy
    
#*****************************************************

# Seawater thermal conductivity
def ThermalConductivity(Temp=25, Salinity=0):
    '''this function calculate the thermal conductivity in W/m.C of water at different temperature (C)
    and different salinity in weight percent. the function is valid in the following ranges
    Temperature range 20 to 180 C and Salinity range 0 to 160 g/kg
    reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-532
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 27/5/2012.'''
    Sal=np.array(Salinity)
    T=np.array(Temp)
    if Sal.ndim==0:
        Sal=np.array([Sal])
    ThermalCond = [10 ** (np.log10(240 + 0.0002 * S * 10) + 0.434 * (2.3 - \
                  (343.5 + 0.037 * S * 10) / (T + 273.15)) * (1 - (T + 273.15) / (647.3 \
                  + 0.03 * S * 10.)) ** (1. / 3.)) / 1000. for S in Sal]

    return ThermalCond
    
#*******************************************************

# Enthalpy of saturated water vapor
def SatVaporEnthalpy(Temp=0):
    '''this function calculate the saturated vapor enthalpy in kj/kg of water at different temperature (C)
    the function is valid in the following range
    Temperature range 0.01 to 200 C
    reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-536
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 1/6/2012.'''
    T=np.array(Temp)
    SatVapEnth = (2501689.845 + 1806.916015 * T + 0.5087717 * T ** 2. - 0.011221 * T ** 3.) / 1000.
    return SatVapEnth
    
#********************************************************

# Water latent heat of evaporation
def LatentHeat(Temp=25):
    '''this function calculate the latent heat of evaporation in kj/kg of water at different temperature (C)
    the function is valid in the following range
    Temperature range 5 to 200 C
    reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-538
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 1/6/2012.'''
    T=np.array(Temp)
    LatHeat = 2501.897149 - 2.407064037 * T + 1.192217 * 10. ** -3. * T ** 2. - 1.5863 * 10 ** -5. * Temp ** 3.
    return LatHeat
    
#********************************************************

# saturation pressure of water vapor
def vap_pressure(Temp=25, Salinity=0):
    '''this function calculate the partial pressure in Pascal (kPa) of water at different temperature (C)
    and different salinity in weight percent. the function is valid in the following ranges
    Temperature range 0 to 200 C. Salinity range 1 to 160 g/kg
    reference: D. Winter*, J. Koschikowski, M. Wieghaus, Desalination using
    membranedistillation: Experimental studies on full scale spiral wound
    modules. Journal of Membrane Science 375 (2011) 104–112.
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 27/5/2012.
    '''
    Sal=np.array(Salinity)
    T=np.array(Temp)
    if Sal.ndim==0:
        Sal=np.array([Sal])
        
    p_result = [100 * 220.93 * np.power(np.e,((647.25 / (T + 273.15)) * (-7.8889166 * (1 - ((T + 273.15) / 647.25))+ \
                2.5514255 * (1 - ((T + 273.15) / 647.25)) ** 1.5 - 6.716169 * (1 - ((T + 273.15) / 647.25)) ** 2 + 33.239495 * \
                (1 - ((T + 273.15) / 647.25)) ** 2.5 - 105.38479 * (1 - ((T + 273.15) / 647.25)) ** 3 + 174.35319 * \
                (1 - ((T + 273.15) / 647.25)) ** 3.5 - 148.39348 * (1 - ((T + 273.15) / 647.25)) ** 4 + 48.631602 * \
                (1 - ((273.15 + T) / 647.25)) ** 4.5)))* (1 / (1 + 0.57357 * (10. * S / (1000. - 10. * S)))) for S in Sal]
    
    return p_result

#*****************************************************************
def BPE(Temp=25,Salinity=0):
    '''this function calculate the boiling point elevation of seawater in C at different temperatures (C)
    and different salinity in weight percent (wt%). the function is valid in the following ranges
    Temperature range 10 to 180 C. Salinity range 0 to 16 wt%
    reference: H.T. El-Dessouky and H.M. Ettouney, Fundamentals of salt water desalination. Elsevier 2002,p-566
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 27/5/2012.
    '''
    Sal=np.array(Salinity)
    T=np.array(Temp)
    if Sal.ndim==0:
        Sal=np.array([Sal])
    BPE_result=[(8.325e-2 + 1.883e-4*T + 4.02e-6*T**2)*S + (- 7.625e-4 + 9.02e-5*T - 5.2e-7*T**2)*S**2 \
                +(1.522e-4 - 3e-6*T - 3e-8*T**2)*S**3 for S in Sal]
    return BPE_result
#******************************************************************************
def SurfaceTension(Temp=25,Salinity=0):
    '''this function calculate the surface tension of seawater in mN/m at different temperature (C)
    and different salinity in weight percent. the function is valid in the following ranges
    Temperature range 0 to 90 C. Salinity range 0 to 121 g/kg
    reference: An Experimental Investigation of the surface tension of seawater, Kishor Govind Nayar, thesis, 19 May 2014.
    creater: Ahmad alsaadi, KAUST Saudi Arabia Thuwal. 27/5/2012.
    '''
    Sal=np.array(Salinity)
    T=np.array(Temp)
    if Sal.ndim==0:
        Sal=np.array([Sal])
    SurfTen = [(235.8*(1-(T+273.15)/647.096)**1.256*(1-0.625*(1-(T+273.15)/647.096)))*(1+3.766e-4*S*10+2.347e-6*S*10*T) for S in Sal]
    return SurfTen
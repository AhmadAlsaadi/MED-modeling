import swp
import numpy as np
import math
def Re(temp,sal,V, dH):
    
    '''this function calculate Reynolds number of water at different temperature
    and different salinity.
    the function require the hydrodynamic diameter in meters and fluid velocity 
    in meter per second to determine Re value 
    the function return Re value.'''
    # creating a vector of different feed temperatures from 5C to 95C in step of 5 C. 
    Dens=np.array(swp.Density(temp,sal))
    Dvis=np.array(swp.DynamicViscosity(temp,sal))  
    # calculating Re  
    Reynolds=Dens*dH*V/Dvis
    return Reynolds

def FluidRegime(Re):
    if Re.mean() <2000:
        FR='Laminar'
    else:
        FR='Turbulent'
    return FR
        
def Pr(Temp,Salinity):
    '''this function calculate Prandtl number at different temperature and different salinity'''
    Prandtl = np.array(swp.SpecificHeat(Temp,Salinity))*np.array(swp.DynamicViscosity(Temp,Salinity))/np.array(swp.ThermalConductivity(Temp,Salinity))
    return Prandtl

def Nu(Temp,Salinity,V,dH,D=1.,L=10.):
    ReValue=Re(Temp,Salinity,V,dH)
    if FluidRegime(ReValue)=='Laminar':
    # for laminar flow regime
        Nusselt=0.13*np.power(ReValue,0.64)*np.power(Pr(Temp,Salinity),.38)
    else:
        Nusselt=0.036*np.power(ReValue,0.8)*np.power(Pr(Temp,Salinity),0.33)*np.power((D/L),.055)
    return Nusselt
def HTC(Temp,Salinity,V,dH):
    '''This function the convective heat transfer coefficient in W/m2.C based on Nusselt
    number at different temperatures and different salinity.
    it takes fluid velocity (m/sec) and channel hydraulic diameter (m) as parameters.'''
    CHTC=Nu(Temp,Salinity,V,dH)*np.array(swp.ThermalConductivity(Temp,Salinity))/dH
    return CHTC
    
    #checked
def SpacerCorrectionFactor(porosity,thickness,filamentDiameter,angle):
    SCF=1.904*(filamentDiameter/thickness)**-0.039*porosity**0.75*(math.sin(math.pi*angle/(360)))**0.086
    return SCF
    
def ChannelHydrolicDiameter(SpacerPorosity,SpacerThickness,FilamentDiameter):
    CHD=4.0*SpacerPorosity*SpacerThickness*FilamentDiameter/(2.0*FilamentDiameter + \
        4.0*(1.0-SpacerPorosity)*SpacerThickness)
    return CHD
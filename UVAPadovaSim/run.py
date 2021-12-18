import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import pandas as pd
import csv
from sklearn.metrics import mean_absolute_error
from basal_val import Basal
from ODEs import *  


def model(z,t, GS, IS, EG, RE, GU, GR, IR, basal, bolus, D, IIRb):

    

    Gp, Gt, Qsto1, Qsto2, Qgut, XL, I1, X, Isc1, Isc2, Ip, Il =z
   

    if t>10:
        D=0
    else:
        D=D/10
    
    
    if t<1:
        IIR=bolus + basal.IIRb
       
    else:
        IIR=basal.IIRb
        
        
     
    EGP = EG.EGP(Gp,XL)
    Uii = GU.Uii()
    E   = RE.E(Gp)
    Uid = GU.Uid(X, Gp, Gt)
    I   = IS.I(Ip)
    
    Qsto= GR.Qsto(Qsto1, Qsto2)
    kempt = GR.kempt(D, Qsto)
    Ra    = GR.Ra(Qgut)
    Ria   = IR.Ria(Isc1, Isc2)
    
    
    dXL_dt = EG.dXL_dt(XL, I1)
    dI1_dt = EG.dI1_dt(I1, I)
    dX_dt  = GU.dX_dt(X, I, basal)

    
    dQsto1_dt = GR.dQsto1_dt(Qsto1, D)
    dQsto2_dt = GR.dQsto2_dt(Qsto1, Qsto2, kempt)
    dQgut_dt  = GR.dQgut_dt(Qgut, Qsto2, kempt)
    
    
        
    dIsc1_dt = IR.dIsc1_dt(Isc1, IIR)
    dIsc2_dt = IR.dIsc2_dt(Isc1, Isc2)


    dIp_dt = IS.dIp_dt(Ip, Il, Ria)
    dIl_dt = IS.dIl_dt(Il, Ip)
    
    dGp_dt = GS.dGp_dt(Gp, EGP, Ra, Uii, Gt, E)
    dGt_dt = GS.dGt_dt(Gt, Gp, Uid)


    return dGp_dt, dGt_dt, dQsto1_dt,dQsto2_dt,dQgut_dt,dXL_dt, dI1_dt, dX_dt, dIsc1_dt, dIsc2_dt, dIp_dt, dIl_dt
    


def simulator(Gb, IIRb, bolus, D):



    loc = 17
    GS = GlucoseSubsystem(loc)
    IS = InsulinSubsystem(loc)
    EG = EndogenousGlucoseProduction(loc)
    RE = RenalExcretion(loc)
    GU = GlucoseUtilisation(loc)
    GR = GlucoseRateOfAppearance(loc)
    IR = InsulinRateOfAppearance(loc)
    p = Params(loc)

    
   
    D = D*1000
    basal = Basal(loc, Gb, IIRb)
    bolus = bolus/ p.params.BW / 60 * 6000


    z0 = [basal.Gpb, basal.Gtb, basal.Qsto1b, basal.Qsto2b, basal.Qgutb, basal.XLb, basal.I1b, basal.Xb, basal.Isc1ss, basal.Isc2ss, basal.Ipb, basal.Ilb]

   

    t = np.linspace(1,60,60)
    y = odeint(model,z0,t, args=(GS, IS, EG, RE, GU, GR, IR, basal, bolus, D, IIRb))

    Vg = p.params.Vg
    simG = np.array([(y[i][0]/Vg) for i in range(len(y))])


    return(simG)






    


    

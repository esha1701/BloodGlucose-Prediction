#from patient_params import Params
import pandas as pd
class Basal():

    def __init__(self, loc, Gb, IIRb):
        patient_params = pd.read_csv("vpatient_params.csv")
        self.params = patient_params.iloc[loc, :]
  
        self.Gpb=Gb * self.params.Vg
        self.IIRb=IIRb/ self.params.BW / 60 * 6000
       
        self.Gtb= (self.params.k1/self.params.k2) * self.Gpb
        self.Ipb=self.IIRb/(self.params.m2 + self.params.m4 - ((self.params.m1 * self.params.m2)/(self.params.m1 + self.params.m2)))
        self.Ilb=self.Ipb * (self.params.m2/(self.params.m1 + self.params.m3))
        self.Ib =self.Ipb/self.params.Vi
        self.I1b=self.Ib
        self.XLb=self.Ib

        self.Isc1ss=self.IIRb/(self.params.kd + self.params.ka1)
        self.Isc2ss=self.Isc1ss * (self.params.kd/self.params.ka2)
 
    Qsto1b=0
    Qsto2b=0
    Qgutb=0
    Xb=0


from patient_params import Params
import numpy as np




class GlucoseSubsystem(Params):
    
    
    def dGp_dt(self, Gp, EGP, Ra, Uii, Gt, E):
        return EGP + Ra - Uii - E - self.params.k1 * Gp + self.params.k2 * Gt
      
 
    def dGt_dt(self, Gt, Gp, Uid):
        return - Uid + self.params.k1 * Gp - self.params.k2 * Gt
    
    def G(self, Gp):
        return Gp/self.params.Vg


class InsulinSubsystem(Params):

    def dIp_dt(self, Ip, Il, Ria):
        return (-(self.params.m2 + self.params.m4) * Ip + self.params.m1 * Il + Ria)

    def dIl_dt(self, Il, Ip):
        return (-(self.params.m1 + self.params.m3) * Il + self.params.m2 * Ip)

    def I(self, Ip):
        return Ip/self.params.Vi


class EndogenousGlucoseProduction(Params):

    def EGP(self, Gp, XL):
        return max(0,self.params.kp1 - self.params.kp2 * Gp - self.params.kp3 * XL)
        
        

    def dXL_dt(self, XL, I1):
        return -self.params.ki * (XL-I1)

    def dI1_dt(self, I1, I):
        return -self.params.ki * (I1-I)

    

class GlucoseUtilisation(Params):

    def Uii(self):
        return self.params.Fcns
        

    def Uid(self, X, Gp, Gt):
        return ((self.params.Vm0 + self.params.Vmx * X) * Gt)/ (self.params.Km0 + Gt)

    def dX_dt(self, X, I, basal):
        return - self.params.p2u * X + (self.params.p2u * (I - basal.Ib))



class GlucoseRateOfAppearance(Params):
   
        
    def Qsto(self, Qsto1, Qsto2):
        return Qsto1 + Qsto2


    def dQsto1_dt(self, Qsto1, D):
        return - self.params.kmax * Qsto1 + D

    def dQsto2_dt(self, Qsto1, Qsto2, kempt):
     
        return -1 * kempt *  Qsto2 + self.params.kmax * Qsto1
        

    def dQgut_dt(self, Qgut, Qsto2, kempt):
        return -self.params.kabs * Qgut + kempt * Qsto2

    def Ra(self, Qgut):
        return (self.params.f * self.params.kabs * Qgut)/self.params.BW

    def kempt(self, D, Qsto):

        if D>0:
            alpha = 5 / (2 * (1 - self.params.b) * D)
            beta  = 5 / (2 * self.params.c * D)
            return self.params.kmin + ((self.params.kmax - self.params.kmin)/2 * (np.tanh(alpha * (Qsto - self.params.b * D)) - np.tanh(beta * (Qsto - self.params.c*D))+2))
       
        else:
            return self.params.kmax






class RenalExcretion(Params):

    def E(self, Gp):
        if Gp> self.params.ke2:
            return self.params.ke1 * (Gp - self.params.ke2)
        else:
            return 0





class InsulinRateOfAppearance(Params):
    
    def Ria(self, Isc1, Isc2):
        return self.params.ka1 * Isc1 + self.params.ka2 * Isc2

    def dIsc1_dt(self, Isc1, IIR):
        return -(self.params.kd + self.params.ka1) * Isc1 + IIR

    def dIsc2_dt(self,Isc1, Isc2):
        return self.params.kd * Isc1 - self.params.ka2 * Isc2



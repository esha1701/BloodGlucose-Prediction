import pandas as pd

class Params:

    def __init__(self,loc):
        patient_params = pd.read_csv("vpatient_params.csv")
        self.params = patient_params.iloc[loc, :]
        
        

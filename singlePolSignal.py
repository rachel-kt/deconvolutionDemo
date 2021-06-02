###### Code written by Ovidiu Radulescu, University of Montpellier, June 2019

####### Required libraries######

import numpy as np

####### signal from one polymerase###########

def Signal_par(ypos,Intensity_for_1_Polym,TailleSeqMarq):
    S = np.ones(len(ypos))*Intensity_for_1_Polym
    ind2 = np.where(ypos < TailleSeqMarq)[0]
    S[ind2] = (1+ypos[ind2])/TailleSeqMarq*Intensity_for_1_Polym      
    return S

# -*- coding: utf-8 -*-
"""
Created on Sat Jan  9 17:37:12 2021

@author: mdouaihy
"""

#### Code written by Ovidiu Radulescu, University of Montpellier, June 2019
##  computes signal given polymerase positions
# input: transcription start positions(Trans_position) in min spacings, Parameters
# output: sum of signals
# call function: getSignal()
# Parameters = {FreqEchSimu, FreqEchImg,DureeSimu,NbrSondeFluo,...
#            ProbeByIntensitie_nb,TaillePreMarq,TailleSeqMarq, TaillePostMarq, Polym_speed};


import numpy as np
from singlePolSignal import Signal_par


def sumSignal1_par(Trans_positions,FreqEchSimu, FreqEchImg, TaillePreMarq, TailleSeqMarq, TaillePostMarq,
                   Polym_speed, frame_num, Intensity_for_1_Polym):
    
### compute signal from positions
    Taille = (TaillePreMarq+TailleSeqMarq+TaillePostMarq)
    Sum_signals_matrix = np.zeros((frame_num,len(Trans_positions)))
    ximage = (np.transpose(np.tile(1+np.arange(frame_num), (len(Trans_positions),1))))/FreqEchImg*Polym_speed #### frame positions in bp
    xpos = np.multiply(np.divide((Trans_positions+1),FreqEchSimu),Polym_speed)-Taille
    t1=np.tile(xpos+TaillePreMarq, (frame_num,1))
    ypos = np.subtract(ximage, t1)
    ind = np.logical_and((ypos > 0),(ypos < (TailleSeqMarq + TaillePostMarq))) #
    Sum_signals_matrix[ind] = Sum_signals_matrix[ind] + Signal_par(ypos[ind]-1,Intensity_for_1_Polym,TailleSeqMarq)
    Sum_signals=np.sum(Sum_signals_matrix, axis = 1)
    return Sum_signals

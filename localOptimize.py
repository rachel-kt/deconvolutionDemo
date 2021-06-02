###########################################################################
##### Code written by Ovidiu Radulescu, University of Montpellier, June 2019
##### this program implements the local optimisation algorithm for polymerase positions %%%
###########################################################################
import numpy as np
from sumSignalDroso import sumSignal1_par

def optimize_local1_par(target,guess,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq,
                        TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym):
    def GD_y_fitness(x):
        return sum((sumSignal1_par(x,FreqEchSimu, FreqEchImg, TaillePreMarq,
                                    TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)-target)**2)  ## positions of poly

    positions = np.where(guess==1)[0]
    Nbr_poly_estimate = len(positions)
    
    shift_window = round(num_possible_poly/Nbr_poly_estimate)+20 ## one position can move [-s_w,s_w]

    Min_fit = GD_y_fitness(positions)
  
    for posi_i in range(len(positions)):
        new_pos = positions.copy()
        for j in range(-shift_window,shift_window+1):
            displaced = positions[posi_i] + j ### displaced position
            if displaced < 0 or displaced > (num_possible_poly -1) or len(np.where( (positions<=displaced) & (positions >= displaced))[0]) !=0:
                continue
            new_pos[posi_i]=displaced
            fitness= GD_y_fitness(new_pos)
            if fitness < Min_fit:
                positions[posi_i] = displaced.copy()
                Min_fit = fitness

    return positions, Min_fit

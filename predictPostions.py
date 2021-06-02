import numpy as np
from sumSignalDroso import sumSignal1_par
from scipy.io import savemat
from scipy.io import loadmat
import timeit
from localOptimize import optimize_local1_par
from gaForDeconvolution import optimize_ga1_par
import matplotlib.pyplot as plt
import timeit


def predictPositions(cellnumber, DataExp, generations, fParam):
        start_time = timeit.default_timer()
        iexp = cellnumber
        print('\nProcessing cell number '+str(iexp+1))
        TolFun=1e-12
        if 'Parameters.npz' in fParam:
            deconParameters=np.load(fParam)

        elif '.mat' in fname:
            deconParameters=loadmat(fParam)

        Polym_speed = deconParameters['Polym_speed'] 
        TaillePreMarq = deconParameters['TaillePreMarq']
        TailleSeqMarq = deconParameters['TailleSeqMarq']
        TaillePostMarq = deconParameters['TaillePostMarq']
        EspaceInterPolyMin = deconParameters['EspaceInterPolyMin']
        FrameLen = deconParameters['FrameLen']
        Intensity_for_1_Polym = deconParameters['Intensity_for_1_Polym']
        FreqEchImg = deconParameters['FreqEchImg']
        DureeSignal = deconParameters['DureeSignal']
        FreqEchSimu = deconParameters['FreqEchSimu']

        sd=DataExp.shape
        frame_num=len(DataExp) ### number of frames
        DureeSimu = frame_num*FrameLen  ### film duration in s
        DureeAnalysee = DureeSignal + DureeSimu ###(s)
        num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed)) # maximal number of polymerase positions

        DataExpSmooth=DataExp[:,iexp] # taking the data cell by cell
        area = FreqEchImg/Polym_speed*Intensity_for_1_Polym*(TaillePostMarq+TailleSeqMarq/2)

        Nbr_poly_estimate = min([round( sum(DataExpSmooth) / area),num_possible_poly]) # first estimation of number of polymerase
                
            
        x_GA_art_list = optimize_ga1_par(TolFun, DataExpSmooth,Nbr_poly_estimate,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq, TailleSeqMarq, TaillePostMarq, 
                                          Polym_speed, frame_num, Intensity_for_1_Polym, generations)
        x_GA_art=  x_GA_art_list[0]


        result_local = optimize_local1_par(DataExpSmooth,x_GA_art,num_possible_poly,FreqEchSimu, FreqEchImg, TaillePreMarq,
                            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)
        print('\nDone processing cell number ',str(iexp+1))
        positions_fit = result_local[0]
        Min_Fit = result_local[1]
        prediction = sumSignal1_par(positions_fit,FreqEchSimu, FreqEchImg, TaillePreMarq,
                                            TailleSeqMarq, TaillePostMarq,  Polym_speed, frame_num, Intensity_for_1_Polym)


        ## store prediction
        #Fit[iexp] = Min_Fit

        #for i in range(len(positions_fit)):
        #    PosPred[positions_fit[i],iexp]=1 # fill Positions of polymerases with 1
        #DataPred[:,iexp] = prediction
        elapsed_time = timeit.default_timer() - start_time
        print("Total time to process cell {} is {}".format(iexp+1, elapsed_time))    
        return [iexp, Min_Fit, prediction, DataExp, positions_fit]

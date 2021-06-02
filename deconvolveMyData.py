# deconvolveMyData

from predictPostions import *
import multiprocessing as mp
import numpy as np
import shutil
import os


def deconvolveMyData(DataFilePath, number_of_workers, fileFormatData='data.npz',):
    outputFolder = DataFilePath   
    fParam = './drosoParameters.npz'
    for DataFileName in os.listdir(DataFilePath):
        if DataFileName.endswith(fileFormatData):
        # we load the result from read_data.m
            print(DataFileName)

            #fname = 'data_'+names[iii]+'CalibratedTraces.npy'
            #fname="../../artificial_data_short/results_artificial_2state/datafiles/data_carola_dataset1_artificial.mat"
            fname = os.path.join(DataFilePath,DataFileName)
            if '.npz' in fname:
                matcontent=np.load(fname)

            elif '.mat' in fname:
                matcontent=loadmat(fname)

            DataExp=matcontent['DataExp']
            if 'Parameters.npz' in fParam:
                deconParameters=np.load(fParam)

            ### calculate data specific parameters

            sd=DataExp.shape
            nloops = sd[1]
            #nloops = 2
            frame_num=len(DataExp) ### number of frames
            FrameLen = deconParameters['FrameLen']
            DureeSignal = deconParameters['DureeSignal']
            DureeSimu = frame_num*FrameLen  ### film duration in s
            DureeAnalysee = DureeSignal + DureeSimu ###(s)
            EspaceInterPolyMin = deconParameters['EspaceInterPolyMin']
            Polym_speed = deconParameters['Polym_speed']
            num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed)) # maximal number of polymerase positions

            PosPred=np.zeros((num_possible_poly,nloops)) # np.zeros(num_possible_poly,(len(DataExp[0]))) # short for positions predictions
            DataPred =np.zeros((sd[0],sd[1])) #signal prediction
            Fit=np.zeros((nloops))


            #### to save the data

            generations=400

            # ------ Parallel pool this part ------ #

            # Step 1: Init multiprocessing.Pool()
            pool = mp.Pool(number_of_workers)

            # Step 2: `pool.starmap` the `predictPositions()`
            # predictPositions(cellnumber, DataExp, generations, fParam):
            result = pool.starmap(predictPositions, [(iexp, DataExp, generations, fParam) for iexp in range(nloops)])

            # Step 3: Don't forget to close
            pool.close()  
            print("pools closed")


            # rearrange the returned results [iexp, Min_Fit, prediction, DataExp, positions_fit]
            for ll in range(len(result)):
                iexp = result[ll][0]
                Fit[iexp] = result[ll][1]
                DataPred[:,iexp] = result[ll][2]
                DataExp = result[ll][3]
                positions_fit = result[ll][4]
                for i in range(len(positions_fit)):
                    PosPred[positions_fit[i],iexp]=1 # fill Positions of polymerases with 1

            #fname = 'python_artificial_data_'+'test_3'+'_00.npz'
            #fname = 'dataset1_deconvolution_2_states.npz'
            outFolder = os.path.join(outputFolder,'resultDec')
            if os.path.exists(outFolder): 
                shutil.rmtree(outFolder, ignore_errors = True)  
            os.mkdir(outFolder)  
            fname = os.path.join(outFolder,DataFileName.replace('_M2_artificial_data.npz','predictions'))
            print("Results saved in {}".format(fname))
            np.savez(fname, Fit=Fit, DataPred=DataPred, DataExp=DataExp, PosPred=PosPred)

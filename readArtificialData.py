#------------------------ Read data ------------------------#

import os
import numpy as np
import pandas as pd

DataFilePath = './demo_data/'
outputFolder = DataFilePath
class readArtificialData:
     def __init__(self,DataFilePath, outputFolder):
        for DataFileName in os.listdir(DataFilePath):
            if DataFileName.endswith(".xlsx"):
                DataPath = os.path.join(DataFilePath, DataFileName)
                DataExp =pd.read_excel(DataPath,header=None).to_numpy()
                DataSaveName=os.path.join(outputFolder,'data_carola_'+DataFileName.replace('.xlsx','')) 
                Frames=np.arange(0,len(DataExp))
                Samples=np.arange(1,len(DataExp)+1)
                np.savez(DataSaveName,DataExp=DataExp,Samples=Samples,Frames=Frames)

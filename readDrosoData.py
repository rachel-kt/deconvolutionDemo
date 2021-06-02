import os
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
import shutil
import timeit

#from joblib import Parallel, delayed

### imput: index of the file
class readDataInFile:
    def __init__(self,DataFilePath):
        PyFilePath = os.path.join(DataFilePath,'npzdatafiles/')
        ### where the images will be stored
        DataFilePath0 = os.path.join(DataFilePath,'images')
        ### creating the folder
        if os.path.exists(PyFilePath): 
            shutil.rmtree(PyFilePath, ignore_errors = True)  
        os.mkdir(PyFilePath) 
        if os.path.exists(DataFilePath0): 
            shutil.rmtree(DataFilePath0, ignore_errors = True)  
        os.mkdir(DataFilePath0) 
        self.PyFilePath = PyFilePath
        self.DataFilePath0=DataFilePath0
        self.DataFilePath=DataFilePath
        self.lists, self.nexp = self.ListingFiles(self.DataFilePath)    
        self.file_name_list=self.lists
        
    def read_data(self,data_i):
        nexp=self.nexp
        file_name_list=self.file_name_list
        DataFileName = self.DataFilePath + file_name_list[data_i] #path to the specified file
        DataExp =pd.read_excel(DataFileName,usecols=np.arange(4,300),skiprows=1).to_numpy() #extracting the data from the excel sheet
        DataExp=DataExp.transpose()[~np.all(DataExp.transpose() == 0, axis=1)].transpose() #avoid zero cells
        Spot_names=pd.read_excel(DataFileName,usecols=np.arange(4,len(DataExp[0])+4),skipfooter =len(DataExp)+1).columns.tolist()
        ## creating a folder to store the images
        dirname=file_name_list[data_i].replace('_CalibratedTraces.xls','')

        WriteTo=self.DataFilePath0+dirname+'_image'
        if os.path.exists(WriteTo): 
            shutil.rmtree(WriteTo, ignore_errors = True)  
        os.mkdir(WriteTo)    

        ## file name and size of dataexp
        n=DataExp.shape
        fn=self.file_name_list[data_i].replace('_','').replace('.xls','')

        ## ploting the intensity spot for each nuclei in a 6*6 plot
        for iii in range(n[1]):

            ifig= math.floor((iii)/36)+1 
            h=plt.figure(ifig+data_i*10)

            plt.subplot(6,6, (iii%36+1))
            plt.plot(np.arange(1,len(DataExp)+1), DataExp[:,iii], color='black', linewidth=0.1)
            h.subplots_adjust(hspace = 0.5, wspace=0.5)
            plt.title(Spot_names[iii].replace('_','-'))

            if (data_i+1)%36==0 or (data_i+1)== nexp:
                plt.savefig(WriteTo+'fig'+str(ifig)+'.pdf')

        #saving the data
        sd=DataExp.shape      
        Frames=np.arange(0,len(DataExp))
        Samples=np.arange(1,len(DataExp)+1)
        fname=self.PyFilePath+fn+'_data.npz'
        self.DataExp=DataExp
        self.Frames=Frames
        self.Samples=Samples
        np.savez(fname,DataExp=DataExp,Samples=Samples,Frames=Frames)
        self.fname=fname
        return self.DataFilePath0
    
    def ListingFiles(self,DataFilePath):
        self.DataFilePath=DataFilePath
        file_name_list = np.array(os.listdir(self.DataFilePath)) # list of subdirectories containing data from different genotypes
        nexp = len(file_name_list) # number of files inside he folder
        #  DataFilePath0='../images/' # place to store images
        #just to avoid erros in the opening of the files
        for i in range(nexp):
                if '._' in file_name_list[i]:
                    file_name_list[i]=file_name_list[i].replace('._','',1)
        file_name_list=np.unique(file_name_list)
        nexp=len(file_name_list)
        return file_name_list, nexp
    

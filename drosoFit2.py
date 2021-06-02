import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
import math
from ecdfEstimate import *
from scipy.optimize import least_squares
import pandas as pd
import timeit
from extentForPlot import *

class fit2:
    def __init__(self,path,parameter):
        self.path=path
        self.parameterpath=parameter
        
        ### parameters used
        filecontent=np.load(self.parameterpath)
        Polym_speed=filecontent['Polym_speed']
        TaillePreMarq=filecontent['TaillePreMarq']
        TailleSeqMarq=filecontent['TailleSeqMarq']
        TaillePostMarq=filecontent['TaillePostMarq']
        EspaceInterPolyMin=filecontent['EspaceInterPolyMin']
        FrameLen=filecontent['FrameLen']
        Intensity_for_1_Polym=filecontent['Intensity_for_1_Polym']
        FreqEchImg=filecontent['FreqEchImg']
        DureeSignal=filecontent['DureeSignal']
        
        FreqEchSimu = 1/(EspaceInterPolyMin/Polym_speed) # how many interval(possible poly start position) in 1s
        self.FreqEchSimu =FreqEchSimu 
        
        ####### parameters for the plots
        fsz=16 #figure size
        lw = 2
        msz=10
        
         ## function needed to set the the parameters for the color map
        cm_jet= plt.cm.get_cmap('jet') # set the colormap to jet array
        
        DataFilePath0 = 'python_Results2'
        if os.path.exists(DataFilePath0):
            shutil.rmtree(DataFilePath0, ignore_errors = True)

        os.mkdir(DataFilePath0)

        xlsfilename = DataFilePath0 + '/notebook_python_results.xlsx'

        # Setting the Names of the Data output in the excel sheet
        xls_cont = pd.DataFrame({'Data': [],'k1p': [], 'k1m': [], 'k2': [], 'p1': [],
                                 'p2': [],'Obj': [],'Nuclei': [],'Frames': []})   

        # Create a Pandas Excel writer using XlsxWriter as the engine.
        writer = pd.ExcelWriter(xlsfilename, engine='xlsxwriter')

        # Convert the dataframe to an XlsxWriter Excel object.
        xls_cont.to_excel(writer, sheet_name='Sheet1', index=False)
        
        ## loading the result of the deconvolution
        NPZFilePath = path;  
        file_name_list = np.array(os.listdir(NPZFilePath)) # list of the data
        self.file_name_list=file_name_list
        nexp = len(file_name_list) # length of the list
        nfiles=nexp
        
        
        #######################################################################
        pooled =0; #### if this is 1, pool all the result files from NPZfilePath
        ############### if not, use each file separately########################
    
        if not pooled:
            nfiles=nexp
        else:
            nfiles=1
            
        ### starting the fit for each file
        for ifile in range(nfiles):
            
            if pooled: 
                
            ### first compute max dimension
                nmax=1
                nmaxpos=1
                
                for iii in range(nexp):
                    fname=file_name_list[iii]
                    ffname=NPZFilePath+fname
                    fnameContent=np.load(ffname)
                    DataPred = fnameContent['DataPred']
                    DataExp=fnameContent['DataExp']
                    Fit = fnameContent['Fit']
                    PosPred=fnameContent['PosPred']
                    n2 =DataExp.shape
                    n3= PosPred.shape
            
                    if n2[0] >nmax:
                        nmax = n2[0]
            
                    if n3[0]> nmaxpos:
                        nmaxpos = n3[0]
               
                #### lump files, Procustes method 
                dataExp=np.empty((nmax, 0), int)
                dataPred=np.empty((nmax,0), int)
                posPred= np.empty((nmaxpos,0), int)
                tmax = np.empty((0, n2[1]), int)

        
            
                for iii in range(nexp): 
                    fname=file_name_list[iii]
                    ffname=NPZFilePath+fname
                    fnameContent=np.load(ffname)
                    DataPred = fnameContent['DataPred']
                    DataExp=fnameContent['DataExp']
                    Fit = fnameContent['Fit']
                    PosPred=fnameContent['PosPred']
                    n2 =DataExp.shape
                    n3= PosPred.shape
      
                    DataExp=np.append(DataExp,np.zeros((nmax-n2[0],n2[1])),axis=0)
                    DataPred=np.append(DataPred,np.zeros((nmax-n2[0],n2[1])),axis=0)
                    PosPred=np.append(PosPred,np.zeros((nmaxpos-n3[0],n3[1])),axis=0)
            
                    # we are adding all the data from different files together
                    dataExp = np.append(dataExp, DataExp, axis=1) 
                    dataPred = np.append(dataPred, DataPred, axis=1)
                    posPred = np.append(posPred, PosPred, axis=1)

                    tmax=np.append(tmax, n2[0]/FreqEchImg*np.ones(n2[1]))
            
        
                DataExp = dataExp.copy()
                DataPred=dataPred.copy()
                PosPred = posPred.copy()
            else:
                fname =   file_name_list[ifile]
                #### full path file name
                ffname = NPZFilePath+ fname
                fnameContent=np.load(ffname)
                DataPred = fnameContent['DataPred']
                DataExp=fnameContent['DataExp']
                Fit = fnameContent['Fit']
                PosPred=fnameContent['PosPred']
                n2=DataExp.shape
                tmax=n2[0]/FreqEchImg*np.ones(n2[1]) #### movie length, the same for all nuclei in a data sets 
  
            ### extract short name from result file name
            iend=fname.index('predictions')
            name=fname[0:iend]   
            self.name=name
            
            ### where to write figure files 
            dirwrite = DataFilePath0+'/'+name+'_result'
            if os.path.exists(dirwrite):
                shutil.rmtree(dirwrite, ignore_errors = True)

            os.mkdir(dirwrite)

            n = DataExp.shape
            nexp = n[1]
            
            
            ## parameters
            DureeSimu = n[0]*FrameLen #in s
            frame_num = n[0]
            DureeAnalysee = DureeSignal + DureeSimu # (s)
            num_possible_poly = round(DureeAnalysee/(EspaceInterPolyMin/Polym_speed))
            
   
            MIntensity = np.array([])
            T0 = np.array([])

            ## find first hit which is 1/5th of the max intensity
            for data_i in range(nexp):
                ifig= math.floor((data_i)/36)+1

                max_intensity=max(DataPred[:,data_i])
                MIntensity=np.append(MIntensity,max_intensity)

                ihit=np.where(DataPred[:,data_i]> max_intensity/5)[0]

                if len(ihit)== 0:
                    ihit=n[0]
                    t0o1 = (ihit)/FreqEchImg 
                else:
                    ihit=min(np.where(DataPred[:,data_i]> max_intensity/5)[0])
                    t0o1 = (ihit+1)/FreqEchImg 

                t0=t0o1   
                T0 = np.append(T0, t0) #stores the  first hit for each nuclei
                
                
                h = plt.figure(ifig+ ifile)
                plt.subplot(6,6, (data_i%36+1))    
                plt.subplots_adjust(hspace=.5)
                plt.fill_between(np.array([t0/60, tmax[data_i]/60]), np.array([150,150]), facecolors =np.array([0.9, 0.9, 0.9])) #this what is really analyzed
                plt.plot(np.arange(0, frame_num)/FreqEchImg/60, DataExp[:,data_i].T, color = 'k', linewidth = 0.1)
                plt.plot(np.arange(0, frame_num)/FreqEchImg/60, DataPred[:,data_i].T, color = 'r', linewidth = 0.1)
                plt.xlim(0, 40)
                plt.ylim(0, 150)
                plt.title(str(round(data_i+1, 3))+':'+str(round(max_intensity,3))  +':' +str(round((tmax[data_i]-t0),3)))

                if data_i%36==35 or (data_i+1)==nexp:
                    figfile=dirwrite+'/figure'+str(ifig)+'.pdf'
                    h.savefig(figfile)
                    plt.close()
        
            ### Figure showing Data Signal Prediction
            h=plt.figure(40)
            sz= DataPred.shape
            Y_normal = np.arange(1,sz[1]+1)
            Y=Y_normal[::-1]
            X = np.arange(0, sz[0])/FreqEchImg/60
            plt.imshow(DataPred.T, cmap=cm_jet, extent=extentForPlot(X).result + extentForPlot(Y).result, aspect='auto', origin='upper')
            plt.xlabel('Time [min]', fontsize=12)
            plt.ylabel('Transcription site', fontsize=12)
            cb= plt.colorbar()
            cb.ax.tick_params(labelsize=fsz)
            figfile=dirwrite+'/DataPred_'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()

            ### Figure showing Data Signal Experimental
            h = plt.figure(50)
            plt.imshow(DataExp.T, cmap=cm_jet, extent=extentForPlot(X).result + extentForPlot(Y).result, aspect='auto', origin='upper')
            plt.xlabel('Time [min]', fontsize=12)
            plt.ylabel('Transcription site', fontsize=12)
            plt.colorbar()
            figfile=dirwrite+'/DataExp_'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()

            ### Figure showing Data Position Prediction
            h=plt.figure(60)
            Y_normal=np.arange(1, len(PosPred[0])+1)
            Y=Y_normal[::-1]
            X=np.arange(0,len(PosPred))*EspaceInterPolyMin/Polym_speed/60  -(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed/60 ### time
            plt.imshow(PosPred.T, cmap='gray', extent=extentForPlot(X).result + extentForPlot(Y).result, aspect='auto', origin='upper')
            plt.xlabel('Time [min]', fontsize=12)
            plt.ylabel('Transcription site', fontsize=12)
            figfile=dirwrite+'/PosPred'+name+'.pdf'
            h.savefig(figfile, dpi=800) 
            plt.close()
            
            ### compute distribution of spacings
            nn=PosPred.shape

            dt=np.array([])
            dtc=np.array([])
            
            figfile=dirwrite +'/PosPred'+name+'.txt' ###  text file for pol positions 
            fid = open(figfile,'w+')

            
            for i in range(nn[1]): #for all cells
                times = (np.where(PosPred[:,i]==1)[0]+1) / FreqEchSimu -(TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed 
                fid.writelines([' \n'+ str(times/60)])

                if len(times) !=0:
                    dtimes = np.diff(times)

                    ### find first index larger than T0

                    istart= min(np.where(times > T0[i] - (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed )[0])
                    dt=np.append(dt, dtimes[istart:])
                    if tmax[i]-times[-1]>0:
                        dtc = np.append(dtc, tmax[i]-times[-1]) ### the very last time

            fid.close()

            #Store: define matrices that store parameters and objective functions 
            # for each iteration of the least_square_function
            # to have better results of least square we ran the function for 100 iterations with random initial values for each iteration
            store = np.empty((0,5))

            if len(dt)!=0: 
                for cens in range(2):
                    if cens:
                        xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc),np.append( np.zeros(len(dt)), np.ones(len(dtc)) ) )
                    else:
                
                        xs,fs,flo,fup =ecdf_bounds(dt )

                    ## fit distribution of spacings using combination of two exponentials

                    xs = xs[:-1]
                    fs = fs[:-1]
                    flo = flo[:-1]
                    fup = fup[:-1]

                    ##############

                    sN = np.sqrt(len(xs))

                    def exp_fitness(k): #objective function
                        return np.abs(np.log(k[2]*np.exp(k[0]*xs)+(1-k[2])*np.exp(k[1]*xs) ) -np.log(1-fs)) /sN # k: parameters

                    k00 = np.array([-0.01, -0.001, 1])
                    amp = np.array([np.log(100), np.log(100)])
                    NbIterationinFit = 100

                    test=0
                    for mc in range(NbIterationinFit):
                        print('iteration nbr ',mc)
                        while test==0: #test is just to re-do the iteration until we encounter no error
                            ## change k00 which is the initial value
                            factor = np.exp(amp* (2* np.random.uniform(size=2)-1))
                            k0 = k00.copy()
                            k0[:2] = k0[0:2]*factor
                            k0[2]=2*np.random.uniform(size=1)-1

                            ## sort k0(1:2)
                            k0[0:2] = np.sort(k0[0:2])

                            ## impose constraints
                            if not (k0[0]*k0[2]+k0[1]*(1-k0[2])<0):
                                while not (k0[0]*k0[2]+k0[1]*(1-k0[2])<0):
                                    k0[2] = 2*np.random.uniform(size=1)-1 # A1, A2 value

                            # Use the fcn lsqnonlin which is the least square function

                            try:
                                k = least_squares(exp_fitness, k0,bounds=(-np.inf,[0,0,np.inf]) ,ftol = (1e-8),max_nfev= 1e6, xtol= (1e-10)).x
                                obj = sum(exp_fitness(k)**2)
                                test = 1
                            except:
                                pass
                        test=0

                        # write down results


                        ## sort k
                        A = np.array([k[2],1-k[2]])  # A values before sorting 
                        kk = np.sort(k[0:2])
                        IX = np.argsort(k[0:2])
                        k[:2]= kk
                        A=A[IX]
                        k[2]=A[0]
                        k_obj=np.append(k,np.array([obj]))
                        to_store=np.append(k_obj,cens)
                        to_store=to_store.reshape(1,len(to_store))
                        store = np.append(store,to_store,axis=0)

                        
                # select optimal 
        
                # ind is index of the least square with real numbers
                ind = np.where(np.max(np.abs(store[:,0:3].imag )  ,axis=1) <1e-10)[0]
                # objmin is minimun of the above results (least square value of ind)
                objmin = np.min(store[ind,3])

                #ind is the index of the minimum value of the least square with real value
                indmin=np.argmin(store[ind,3])
                imin = ind[indmin]

                #overflow help us set the Uncertainty interval
                overflow = 1

                # ind index where the least square are real  and less than < (1+overflow)*objmin)
                ind = np.where( (store[:,3] < (1+overflow)*objmin) & (np.max(np.abs(store[:,0:3].imag   ),axis=1) <1e-10)) [0]
                ksel = store[ind,0:3].real
                kmin = store[imin,0:3].real
                censmin=store[imin,4].real

                if censmin:
                    xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc),np.append( np.zeros(len(dt)), np.ones(len(dtc)) ) )
                else:
                    xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc))


                ## Survival function
                h = plt.figure(70)
                plt.semilogy(xs, 1-fs, 'o', color='r', mfc = 'none', markersize=9, linestyle='') #empirical function
                plt.semilogy(xs, 1-flo, '--r') # lower confidence
                plt.semilogy(xs, 1-fup, '--r') # upper confidence
                pred = kmin[2]*np.exp(kmin[0]*xs)+(1-kmin[2])*np.exp(kmin[1]*xs)
                plt.semilogy(xs, pred, 'k', linewidth = 3) # predicted 2 exp
                plt.xlim(0,250)
                plt.ylim(1e-6, 1)
                plt.xlabel('Time [s]',fontsize=fsz)
                plt.ylabel('Survival function',fontsize = fsz)

                # compute 3 rates k1p,m k2 from the 3 parameters
                S1 = kmin[2]*kmin[0]+(1-kmin[2])*kmin[1]
                k2 = -S1
                S2 = kmin[2]*(kmin[0])**2+(1-kmin[2])*(kmin[1])**2
                S3 = kmin[2]*(kmin[0])**3+(1-kmin[2])*(kmin[1])**3
                k1m = S1-S2/S1
                k1p = (S3*S1-S2**2)/S1/(S1**2-S2)
                plt.title('k_1^-='+str("%.1g" % k1m)+'k_1^+='+str("%.1g" % k1p)+'k_2='+str("%.1g" % k2))

                figfile = dirwrite+'/Fit_'+name+'.pdf'
                h.savefig(figfile)
                plt.close()


                # compute interval
                S1 = ksel[:,2]*ksel[:,0]+(1-ksel[:,2])*ksel[:,1]
                K2 = -S1
                S2 = ksel[:,2]*ksel[:,0]**2+(1-ksel[:,2])*ksel[:,1]**2
                S3 = ksel[:,2]*ksel[:,0]**3+(1-ksel[:,2])*ksel[:,1]**3
                K1m = S1-S2/S1
                K1p = (S3*S1-S2**2)/S1/(S1**2-S2)
                #################

                ## optimal
                res = np.array([k1p,k1m,k2,k1m/(k1m+k1p),k1p/(k1m+k1p),objmin,len(DataExp[0]),len(DataExp)])

                P1=K1m/(K1m+K1p)
                P2=K1p/(K1m+K1p)

                resl= np.array([np.min(K1p),np.min(K1m),np.min(K2),np.min(P1),np.min(P2)])

                resl = np.max(np.vstack([resl, np.zeros(5)]), axis=0)


                resh= np.array([np.max(K1p),np.max(K1m),np.max(K2),np.max(P1),np.max(P2)])

                df1 = pd.DataFrame([res.tolist(), #best result
                                    resl.tolist(), # low 
                                    resh.tolist()]) # high


                df1.to_excel(writer,sheet_name='Sheet1', startrow=4*(ifile+1)-2, startcol=1, header=False, index=False)

                df2 = pd.DataFrame([name.replace('result_','')]) #filename
                df2.to_excel(writer,sheet_name='Sheet1', startrow=4*(1+ifile)-2, startcol=0, header=False, index=False)


            else:
                df2 = pd.DataFrame([name.replace('result_','')]) #filename
                df2.to_excel(writer,sheet_name='Sheet1', startrow=4*(1+ifile)-2, startcol=1, header=False, index=False)
                res = [0, 0, 0, 0, 0, 0, len(DataExp[0]), len(DataExp)]
                df1 = pd.DataFrame(res)
                df1.to_excel(writer,sheet_name='Sheet1', startrow=4*(ifile+1)-2, startcol=2, header=False, index=False)
        writer.save()


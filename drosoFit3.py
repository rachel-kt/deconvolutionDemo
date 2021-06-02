import numpy as np
import os
import shutil
import matplotlib.pyplot as plt
import math
from scipy.optimize import least_squares
import pandas as pd
import openpyxl as pxl
import timeit
from ecdfEstimate import *
from extentForPlot import *

class fit3:
    def __init__(self,path,parameterpath):
        self.path=path
        
        self.parameterpath=parameterpath

        ### parameters used
        filecontent=np.load(parameterpath)
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
        
        ### folder to write  output files
        DataFilePath0 = 'python_Results3'
        if os.path.exists(DataFilePath0):
            shutil.rmtree(DataFilePath0, ignore_errors = True)
        os.mkdir(DataFilePath0)
        
        ### creating the excel sheet for the results
        xlsfilename = DataFilePath0 + '/notebook_python_results_03.xlsx'
         
            # Setting the Names of the Data output in the excel sheet
        xls_cont = pd.DataFrame({'Data': [],'M1:k1p' : [],'M1:k1m' : [],'M1:k2p' : [],
                                 'M1:k2m' : [],'M1:k3' : [],'M1:p1' : [],'M1:p2': [],
                                 'M1:p3' : [],'M2:k1p': [],'M2:k1m': [],'M2:k2p': [],
                                 'M2:k2m': [],'M2:k3': [],'M2:p1': [],'M2:p2': [],
                                 'M2:p3': [],'lambda1': [],'lambda2': [],
                                 'lambda3': [],'A1': [],'A2': [],'A3': [],'S1': [],
                                 'Obj': [],'Nuclei': [],'Frames': []})  
        
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
            dirwrite = DataFilePath0+'/'+name+'_result3'
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
            
            ## store the first hit for each nuclei
            T0 = np.array([])
            MIntensity = np.array([]
                                 )
            ## find first hit which is 1/5th of the max intensity
            for data_i in range(nexp):
                ifig= math.floor((data_i)/36)+1

                max_intensity=max(DataPred[:,data_i])

                ###  find first hit
                ihit=np.where(DataPred[:,data_i]> max_intensity/5)[0]

                if len(ihit)== 0:
                    ihit=n[0]
                    t0 = (ihit)/FreqEchImg # t0 in s Option 1 : 1/5th of the max intensity
                else:
                    ihit=min(np.where(DataPred[:,data_i]> max_intensity/5)[0])
                    t0 = (ihit+1)/FreqEchImg # t0 in s Option 1 : 1/5th of the max intensity

                T0 = np.append(T0, t0) # stores the  first hit for each nuclei

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

                if len(times) !=0 and  T0[i] < tmax[i]:
                    dtimes = np.diff(times)

                    ### find first index larger than T0

                    istart= min(np.where(times > T0[i] - (TaillePreMarq+TailleSeqMarq+TaillePostMarq)/Polym_speed )[0])
                    dt=np.append(dt, dtimes[istart:])
                    dtc = np.append(dtc, tmax[i]-times[-1]) ### the very last time

            fid.close()

            #Store: define matrices that store parameters and objective functions 
            # for each iteration of the least_square_function
            # to have better results of least square we ran the function for 100 iterations with random initial values for each iteration
            store = np.empty((0,6))

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
                        return (np.log(k[3]*np.exp(k[0]*xs)+k[4]*np.exp(k[1]*xs)+(1-k[3]-k[4])*np.exp(k[2]*xs) ) -np.log(1-fs)) /sN # k: parameters

                    k00 = np.array([-0.01, -0.01,-0.001, 0.25, 0.25])
                    amp = np.array([np.log(100), np.log(100), np.log(100), np.log(1), np.log(1)])

                    NbIterationinFit = 100

                    test=0
                    for mc in range(NbIterationinFit):
                        print('iteration nbr ', mc)
                        while test==0: #test is just to re-do the iteration until we encounter no error

                            ## change k00 which is the initial value
                            factor = np.exp(amp* (2* np.random.uniform(size=5)-1))
                            k0 = k00*factor
                            k0[3:5] = 2*np.random.uniform(size=2)-1

                            ## sort k0(1:2)
                            k0[0:3] = np.sort(k0[0:3])

                            ## impose constraints
                            if not (sum(k0[3:5])<1 and k0[0]*k0[3]+ k0[1] * k0[4] +k0[2]*(1-sum(k0[3:5])) <0):
                                while not (sum(k0[3:5])<1 and k0[0]*k0[3]+ k0[1] * k0[4] +k0[2]*(1-sum(k0[3:5])) <0):
                                    k0[3:5] = 2*np.random.uniform(size=2)-1 # A1, A2 value

                            # Use the fcn lsqnonlin

                            try:
                                k = least_squares(exp_fitness, k0,bounds=(-np.inf,[0,0,0,np.inf,np.inf]), ftol = (1e-8), xtol= (1e-10), method='trf').x
                                obj = sum(exp_fitness(k)**2)
                                test = 1
                            except:
                                pass
                        test=0

                        # write down results


                        ## sort k
                        A = np.array([k[3],k[4], 1-k[3]-k[4]])  # A values before sorting 
                        kk = np.sort(k[0:3])
                        IX = np.argsort(k[0:3])
                        k[:3]= kk
                        A=A[IX]
                        k[3:5]=A[0:2]

                        store = np.append(store,np.append(k,np.array([obj])).reshape(1,len(np.append(k,np.array([obj])))),axis=0)

                # select optimal and suboptimal obj < 2 objmin
                ind = np.where(np.max(np.abs(store[:,0:3].imag )  ,axis=1) <1e-10)[0] #takes the part where we have the imaginary part of lambda i almost 0
                objmin = np.min(store[ind,-1]).real #minimum of the objectives for all the previous lambdas
                indmin=np.argmin(store[ind,-1]) #finding the positions of the lambda where we have the lowest objective
                imin = ind[indmin] #finding the index of the positions of the lambda where we have the lowest objective
                overflow = 1 #help us set the Uncertainty interval

                ind = np.where( (store[:,-1] < (1+overflow)*objmin) & (np.max(np.abs(store[:,0:3].imag   ),axis=1) <1e-10)) [0] #taking suboptimal objectives such that they are <2*optimal objective and they satisfy the fact that the imaginary part of lambda is almost 0
                ksel = store[ind,0:5].real #taking an array of lambda i's and Ai's where we have the the objective function less than <2*times the minimum objective function
                kmin = store[imin,0:5].real# taking the lambda i's and A i's that are the fittest
                censmin=store[imin,-1]

                if censmin:
                    xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc),np.append( np.zeros(len(dt)), np.ones(len(dtc)) ) )
                else:
                    xs,fs,flo,fup =ecdf_bounds(np.append(dt,dtc))


                h = plt.figure(70)
                plt.semilogy(xs, 1-fs, 'o', color='r', mfc = 'none', markersize=9, linestyle='') #empirical function
                plt.semilogy(xs, 1-flo, '--r') # lower confidence
                plt.semilogy(xs, 1-fup, '--r') # upper confidence
                pred = kmin[3]*np.exp(kmin[0]*xs)+kmin[4]*np.exp(kmin[1]*xs)+(1-kmin[3]-kmin[4])*np.exp(kmin[2]*xs) 
                plt.semilogy(xs, pred, 'k', linewidth = 2) # predicted 2 exp
                plt.xlim(0,250)
                plt.ylim(1e-6, 1)
                plt.xlabel('Time [s]',fontsize=fsz)
                plt.ylabel('Survival function',fontsize = fsz)

                # compute 5 rates k1p,m k2p,m k3 from the 5 parameters
                l1 = kmin[0]
                l2 = kmin[1]
                l3 = kmin[2]
                A1 = kmin[3]
                A2 = kmin[4]
                A3 = 1-A1-A2
                L1=l1+l2+l3
                L2=l1*l2+l1*l3+l2*l3
                L3=l1*l2*l3
                S1=A1*l1+A2*l2+A3*l3
                S2=A1*l1**2+A2*l2**2+A3*l3**2
                S3=A1*l1**3+A2*l2**3+A3*l3**3

                ### Model M2
                kk3= -S1 
                kk1p = 1/2 * ( -L1+S2/S1 +  np.sqrt((S1*L1-S2) **2-4*L3*S1)/S1 ) 
                kk2p = 1/2 * ( -L1+S2/S1 -  np.sqrt((S1*L1-S2) **2-4*L3*S1)/S1 ) 
                kk1m = 1/2 * (S1-S2/S1 - (-S1 **2*L1+S1*S2+S1*L2-L3+S2 **2/S1-S3)/ np.sqrt((S1*L1-S2) **2-4*L3*S1)) 
                kk2m = 1/2 * (S1-S2/S1 + (-S1 **2*L1+S1*S2+S1*L2-L3+S2 **2/S1-S3)/ np.sqrt((S1*L1-S2) **2-4*L3*S1)) 
                pp1=kk1m*kk2p/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m) 
                pp2=kk1p*kk2m/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m) 
                pp3=kk1p*kk2p/(kk1p*kk2p+kk1m*kk2p+kk1p*kk2m) 
                # model M1
                k1p=-L3 *(S1  **2-S2) /(S2  **2-S1 *S3)   # k1p
                k2p=-(S2  **2-S1 *S3) /S1 /(S1  **2-S2)   # k2p
                k3=-S1   # k3
                k2m=(S1  **2-S2) /S1   # k2m
                k1m=-A1 *A2 *A3 *(l1-l2)  **2 *(l1-l3)  **2 *(l2-l3)  **2 *S1 /(S1  **2-S2) /(S2  **2-S1 *S3)   # k1m
                p1=k1m*k2m/(k1p*k2m+k1m*k2m+k1p*k2p) 
                p2=k1p*k2m/(k1p*k2m+k1m*k2m+k1p*k2p) 
                p3=k1p*k2p/(k1p*k2m+k1m*k2m+k1p*k2p) 

                plt.title('k_1^-='+str("%.1g" % k1m)+'k_1^+='+str("%.1g" % k1p) +'k_2^-='+str("%.1g" % k2m)+'k_2^+='+str("%.1g" % k2p)+'k_3='+str("%.1g" % k3))

                figfile = dirwrite+'/Fit_'+name+'.pdf' 
                h.savefig(figfile)
                plt.close()

                nn = DataExp.shape

                #optimal 

                res= np.array([k1p,k1m,k2p,k2m,k3,p1,p2,p3,kk1p,kk1m,kk2p,kk2m,kk3,pp1,pp2,pp3,l1,l2,l3,A1,A2,A3,S1,objmin,nn[1],nn[0]]) # optimum

                # compute interval
                ##################################################
                l1=ksel[:,0] 
                l1sel=l1.copy()
                l2=ksel[:,1] 
                l2sel=l2.copy()
                l3=ksel[:,2] 
                l3sel=l3.copy()
                A1=ksel[:,3] 
                A1sel=A1.copy()
                A2=ksel[:,4] 
                A2sel=A2.copy()
                A3=1-A1-A2 
                A3sel=A3.copy()
                L1=l1+l2+l3 
                L2=l1 *l2+l1 *l3+l2 *l3 
                L3=l1 *l2 *l3 
                S1=A1 *l1+A2 *l2+A3 *l3 
                S1sel=S1.copy()
                S2=A1 *l1  **2+A2 *l2  **2+A3 *l3  **2 
                S3=A1 *l1 **3+A2 *l2 **3+A3 *l3 **3 

                #### model M2
                KK1p = 1/2 * ( -L1+S2 /S1 + np.sqrt((S1 *L1-S2) **2-4*L3 *S1) /S1 ) 
                KK2p = 1/2 * ( -L1+S2 /S1 - np.sqrt((S1 *L1-S2) **2-4*L3 *S1) /S1 )  
                KK1m = 1/2 * (S1-S2 /S1 - (-S1 **2 *L1+S1 *S2+S1 *L2-L3+S2 **2 /S1-S3) /np.sqrt((S1 *L1-S2) **2-4*L3 *S1))  
                KK2m = 1/2 * (S1-S2 /S1 + (-S1 **2 *L1+S1 *S2+S1 *L2-L3+S2 **2 /S1-S3) /np.sqrt((S1 *L1-S2) **2-4*L3 *S1))  
                KK3=-S1  #### k3
                PP1=KK1m *KK2p /(KK1p *KK2p+KK1m *KK2p+KK1p *KK2m) 
                PP2=KK1p *KK2m /(KK1p *KK2p+KK1m *KK2p+KK1p *KK2m) 
                PP3=KK1p *KK2p /(KK1p *KK2p+KK1m *KK2p+KK1p *KK2m) 
                #### model M1
                K1p=-L3 *(S1 **2-S2) /(S2 **2-S1 *S3)  ### k1p
                K2p=-(S2 **2-S1 *S3) /S1 /(S1 **2-S2)  ### k2p
                K3=-S1.copy()  #### k3
                K2m=(S1 **2-S2) /S1  ### k2m
                K1m=-A1 *A2 *A3 *(l1-l2) **2 *(l1-l3) **2 *(l2-l3) **2 *S1 /(S1 **2-S2) /(S2 **2-S1 *S3)  ### k1m
                P1=K1m *K2m /(K1p *K2m+K1m *K2m+K1p *K2p) 
                P2=K1p *K2m /(K1p *K2m+K1m *K2m+K1p *K2p) 
                P3=K1p *K2p /(K1p *K2m+K1m *K2m+K1p *K2p) 


                resl= np.array([np.nanmin(K1p),np.nanmin(K1m),np.nanmin(K2p),np.nanmin(K2m),np.nanmin(K3),np.nanmin(P1),np.nanmin(P2),np.nanmin(P3),np.nanmin(KK1p),
                                np.nanmin(KK1m),np.nanmin(KK2p),np.nanmin(KK2m),np.nanmin(KK3),np.nanmin(PP1),np.nanmin(PP2),np.nanmin(PP3),np.nanmin(l1sel),np.nanmin(l2sel),np.nanmin(l3sel),np.nanmin(A1sel),np.nanmin(A2sel),np.nanmin(A3sel),np.nanmin(S1sel)])

                resl[0:16] = np.max(np.vstack([resl[0:16], np.zeros(16)]), axis=0)


                resh= np.array([np.nanmax(K1p),np.nanmax(K1m),np.nanmax(K2p),np.nanmax(K2m),np.nanmax(K3),np.nanmax(P1),np.nanmax(P2),np.nanmax(P3),np.nanmax(KK1p),np.nanmax(KK1m),np.nanmax(KK2p),np.nanmax(KK2m),np.nanmax(KK3),np.nanmax(PP1),np.nanmax(PP2),np.nanmax(PP3),np.nanmax(l1sel),np.nanmax(l2sel),np.nanmax(l3sel),np.nanmax(A1sel),np.nanmax(A2sel),np.nanmax(A3sel),np.nanmax(S1sel)])


                df1 = pd.DataFrame([res.tolist(), #best result
                                    resl.tolist(), # low 
                                    resh.tolist()]) # high

                df1.to_excel(writer,sheet_name='Sheet1', startrow=4*(ifile+1)-2, startcol=1, header=False, index=False)
                df2 = pd.DataFrame([name.replace('result_','')]) #filename
                df2.to_excel(writer,sheet_name='Sheet1', startrow=4*(ifile+1)-2, startcol=0, header=False, index=False)


            else:
                df2 = pd.DataFrame([name.replace('result_','')]) #filename
                df2.to_excel(writer,sheet_name='Sheet1', startrow=4*(ifile+1)-2, startcol=0, header=False, index=False)
                res = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, nn[1], nn[0]]
                self.df1 = pd.DataFrame(res)
                df1.to_excel(writer,sheet_name='Sheet1', startrow=4*(ifile+1)-2, startcol=1, header=False, index=False)
                
        writer.save()

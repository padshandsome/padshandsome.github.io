# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 09:04:20 2020

@author: compu
"""


# -*- coding: utf-8 -*-
"""
Created on Fri Sep  4 11:36:31 2020


@author: compu
"""

import sys, os
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplot
import scipy.signal
from scipy.signal import detrend, cheby2, filtfilt, butter, argrelextrema, hilbert, detrend, find_peaks, lfilter,savgol_filter, blackmanharris,resample
from scipy.fftpack import dct,fft
from biosppy.signals import ecg
import hrvanalysis
import PPG_RR_INTERVALT
def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data,lowcut, highcut, fs, order=5):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def hrv_analysis(rrintervals):
    #rrintervals=rrintervals*1000
    return hrvanalysis.extract_features.get_time_domain_features(rrintervals)
sigma=1

#唯一需要修改的地方，哦，后面保存数据和图片也可能要根据文件来修改
file='sch_opensignals_0007804d2e7c_2020-07-29_10-40-01.h5'#ecg文件名称
threshold_max=220;threshold_min=30;drawpicture,savefile=0,1
runecg=True;
ecghour=file[-11:-9];ecgminute=file[-8:-6];ecgsecond=file[-5:-3]
#ecghour=9;ecgminute=18;ecgsecond=19;
start=int(ecghour)*3600+int(ecgminute)*60+int(ecgsecond)#改ecg文件时间，就是取上面文件的时间
userid="724"
fstart,fend=564,918#ppg文件的位置 600,893 46,387;0,320;534,922;532,878;512,875；529,882;526,865;538,838



#filepath按照文件保存的路径来改
date=file[-22:-12]
#date=file[-14:-4]
path="D:/work/"
filepath=path+"/ppg+mot_"+userid+"/ppg/"+date+"/"+userid+"/"
motionpath=path+"/ppg+mot_"+userid+"/mot/"+date+"/"+userid+"/"
ecgpath='D:/ECG data/data/ECGsignal_63267/ECGsignal/'+date+'/'+userid+'/'+file

'''
files=os.listdir(filepath)
N_ecg = pd.read_csv(ecgpath)
ecgtesttime=N_ecg['datetime']
N_ecg = N_ecg['ecgraw']
'''
files=os.listdir(filepath)
ecgfs=300

f = h5py.File(ecgpath,'r')
f.keys()  
if(ecgpath[-25]=='a' or ecgpath[-25]=='A'):
    ecg_channel_1 = f['00:07:80:46:F2:A8/raw/channel_1'][:]
else:
    ecg_channel_1 = f['00:07:80:4D:2E:7C/raw/channel_1'][:]
    
ecg_channel_1 = ecg_channel_1.astype(np.float64)


if(userid=='724' or userid=='667'):
    N_ecg = (ecg_channel_1)
else:
    N_ecg = -(ecg_channel_1)


#N_ecg=(ecg_channel_1)
ecg1 = np.squeeze(N_ecg)

fs=25
j=0
hrvresults=pd.DataFrame();ppgresults=pd.DataFrame();ecgresults=pd.DataFrame();ecghr_min=pd.DataFrame();ecghr_max=pd.DataFrame();ecghr_mean=pd.DataFrame();rr_test=pd.DataFrame()

timelist=[]
index=[]
rrinterval=pd.DataFrame()
fposition=fstart
for file in files[fstart:fend]:
    results=[]
    #hrvresults=[]
    parameters=[]
    
    ppg0=np.loadtxt(filepath+file,delimiter=',')
    
    '''
    修改ppg窗口
    '''
    
    ppg0=ppg0 #修改窗
    winlen="1 minute"
    pics=1#第一段
    #ppg0=-ppg0
    ppg1=ppg0
    
    hours = int(file[-12:-10])# 时
    mins = int(file[-9:-7])# 分
    seconds = int(file[-6:-4])
    
    if True:
      all_sec = 3600*hours+60*mins+seconds
    
      start_index = (all_sec - start)*ecgfs
    
      end_index = start_index+60*ecgfs
      print(str(start_index)+"-"+str(end_index))
      if(start_index<0):
          print("ecg time is invalid!")
          print(start_index,end_index)
          print(str(hours)+"-"+str(mins)+"-"+str(seconds))
          break
      '''
      datedate="2020/10/20"+" "+file[-12:-10]+":"+file[-9:-7]+":"+file[-6:-4]#根据ppg文件找ecg信号的位置，datedate是fitbit ecg自带的时间戳
      for newecgtimeindex in range(len(ecgtesttime)):
          if(ecgtesttime[newecgtimeindex]==datedate):
              break;
            
      start_index = newecgtimeindex
      end_index = start_index + 60*ecgfs
      if(end_index>len(ecgtesttime)):
          continue;
     '''
     
      ecg2 = ecg1[start_index:end_index]
      
      ecgsave=pd.DataFrame(ecg2)
    
    
      time=str(hours)+str(mins)+str(seconds)
      print(str(hours)+":"+str(mins)+":"+str(seconds))
      '''
       修改ecg窗口，与ppg窗口对应
      '''
     #ecg2=ecg2[0:int((end_index-start_index)/6)] 
    
      ecgheartbeat_min=[]
      ecgheartbeat_max=[]
      ecgheartbeat_mean=[]
      
      rpeakk = ecg.hamilton_segmenter(ecg2, sampling_rate=ecgfs)[0]
      if drawpicture:
          
          plt.figure(figsize=(100,10))
          #plt.title(date+" "+str(len(rpeakk))
          plt.title(date+" "+str(len(rpeakk)))
          plt.plot(ecg2)
          plt.plot(rpeakk,ecg2[rpeakk],'x')
          #plt.plot(rpeakk[2],ecg2[rpeakk[2]],'o')
          plt.savefig("C:/Users/compu/Desktop/"+"jj.png")
         
        
      rr_secgwhole=[]
      for i in range(0,10):
            
          ecg_win=ecg2[i*6*ecgfs:(i+1)*6*ecgfs]
          r_peaks_1 = ecg.hamilton_segmenter(ecg_win, sampling_rate=ecgfs)[0]
          '''
          plt.figure(figsize=(20,10))
          plt.title(date+" "+time)
          plt.plot(ecg_win)
          plt.plot(r_peaks_1,ecg_win[r_peaks_1],'x')
          '''
          rr_secg = np.diff(r_peaks_1)/(ecgfs)
          
          rr_secg=rr_secg*1000
          rr_secgwhole.append(rr_secg)
          ecgresults0=hrv_analysis(rr_secg)
            
          ecgheartbeat_min.append(ecgresults0['min_hr'])
          ecgheartbeat_max.append(ecgresults0['max_hr'])
          ecgheartbeat_mean.append(ecgresults0['mean_hr'])
       
    
      
      #rr_secgwhole=np.diff(r_peaks)/(ecgfs)
      #rr_secgwhole=rr_secgwhole*1000
      
      rr_secgwhole1=[]
      for i in range(0,10):
          for j in rr_secgwhole[i]:
              rr_secgwhole1.append(j)
      ecgresults0=hrv_analysis(rr_secgwhole1)
        
    
        
     
      #ecgrealtime=r_peaks
        
      diffindex=0; 
      mean_diff=np.diff(ecgheartbeat_mean);
       
      print(ecgresults0['min_hr'])
      if(ecgresults0['range_nni']>500) or True:
          for maxhr,minhr,diff,diffindex in zip(ecgheartbeat_max,ecgheartbeat_min,mean_diff,range(len(mean_diff))):
              if(maxhr>threshold_max) or (minhr<threshold_min) or diff>10 or diff<-10:
                   
                  ecgheartbeat_min[diffindex]=0;ecgheartbeat_min[diffindex+1]=0;
                  ecgheartbeat_max[diffindex]=0;ecgheartbeat_max[diffindex+1]=0;
                  ecgheartbeat_mean[diffindex]=0;ecgheartbeat_mean[diffindex+1]=0;
          if ecgheartbeat_max[9]>threshold_max or ecgheartbeat_min[9]<threshold_min:
              ecgheartbeat_max[9]=0;ecgheartbeat_min[9]=0;ecgheartbeat_mean[9]=0;
          final_rr=[];final_ecg_time=[]
          for i in range(0,10):
            time0=i*6*ecgfs;time1=(i+1)*6*ecgfs
            if(ecgheartbeat_min[i]!=0):
                final_rr.append(rr_secgwhole[i]); 
      else:
          final_rr = rr_secgwhole
          
      if(len(np.where(np.array(ecgheartbeat_mean) ==0)[0])>=6):
          continue;
      
      final_rr1=[]
      for i in range(0,len(final_rr)):
          for j in final_rr[i]:
              final_rr1.append(j)
      
      
      diff_rr=np.diff(final_rr1)
      reject_rr=[]
      
      '''
      for rrindex in range(0,len(diff_rr)-2):
          sumr=sum(diff_rr[rrindex:rrindex+3])
          abssumr=sum(abs(diff_rr[rrindex:rrindex+3]))
          
          if(sumr!=0 and abssumr/abs(sumr)>25):
              reject_rr.append(final_rr[rrindex])
              reject_rr.append(final_rr[rrindex+1])
              reject_rr.append(final_rr[rrindex+2])
              reject_rr.append(final_rr[rrindex+3])
          elif(sumr==0 and abssumr>500):
              #reject_rr.append(final_rr[rrindex])
              reject_rr.append(final_rr[rrindex+1])
              reject_rr.append(final_rr[rrindex+2])
              reject_rr.append(final_rr[rrindex+3])
      '''
     
        
      ecgresults0=hrv_analysis(final_rr1)
        
      rr_test=pd.concat([rr_test,pd.DataFrame(rr_secgwhole).T],axis=0)
      ecgresults0=pd.DataFrame([ecgresults0])
      ecgresults=pd.concat([ecgresults,ecgresults0],axis=0)
        
      ecgheartbeat_min=pd.DataFrame(ecgheartbeat_min) 
      ecgheartbeat_max=pd.DataFrame(ecgheartbeat_max)
      ecgheartbeat_mean=pd.DataFrame(ecgheartbeat_mean)
    
      ecgheartbeat_min = ecgheartbeat_min.T
      ecgheartbeat_max = ecgheartbeat_max.T
      ecgheartbeat_mean=ecgheartbeat_mean.T
    
      ecghr_min=pd.concat([ecghr_min,ecgheartbeat_min],axis=0)
      ecghr_max=pd.concat([ecghr_max,ecgheartbeat_max],axis=0)
      ecghr_mean=pd.concat([ecghr_mean,ecgheartbeat_mean],axis=0)
          
    
    ppg0=ppg0[8:]
    ppg0=ppg0[:-8]
    
    
    ppg0=scipy.signal.resample(ppg0,60*ecgfs)
    #ppg0=np.sign(np.diff(ppg0))*(1-1/(1+np.square(np.diff(ppg0)/sigma))) # probabilistic filter ***
    bprob25,aprob25 = cheby2(4,40,[30/60/(ecgfs/2),220/60/(ecgfs/2)],'band')#15-200bpm #40,250
    #bprob25,aprob25=cheby2(5,40,[.01,.25],'band')
    ppg0f=filtfilt(bprob25,aprob25,ppg0)
    
    #peaks=find_peaks(ppg0f,prominence=1)[0]
    #ppg1f=ppg1f/np.median(ppg1f)*np.median(ppg0f)
    
    
    
    
    acc=np.loadtxt(motionpath+file,delimiter=',',skiprows=1)
    accn=np.linalg.norm(acc,axis=1)
    baccn=butter_bandpass_filter(detrend(accn),0.4,4,25,order=4)

   
    #rr_s=np.diff(peaks)/(fs)
    _,_,_,_,_,_,_,_,rr_s=PPG_RR_INTERVALT.ppg_process(1,ppg0f,ppg0f[:60*ecgfs],baccn,sample_rate=ecgfs,Filtered_Show=False,Raw_Data_Show=False,R_Peak_Show=False,RR_interval_Show=False,label1="ppg",cheby2f=-1)  
    #plt.savefig("D:/work/ppg+mot_"+userid+"/"+"/ppg/"+date+"/hrvanalysis/rejection pro pics/"+date+"_"+time2+".jpg")
    #plt.close()
    #rr_s1,_,_,_,_,_,_,_=PPG_RR_INTERVALT.ppg_process(1,ppg1f,ppg1f[:1500],sample_rate=25,Filtered_Show=True,Raw_Data_Show=False,R_Peak_Show=True,RR_interval_Show=False,label1="ppg",cheby2f=0)  
    
    
    #rr_s=rr_s/1000
    time1=str(hours)+":"+str(mins)+":"+str(seconds)
    if runecg:
        if(len(str(hours))==1):
          hours="0"+str(hours)
        if(len(str(mins))==1):
          mins="0"+str(mins)
        if(len(str(seconds))==1):
          seconds="0"+str(seconds)
        time2=str(hours)+"-"+str(mins)+"-"+str(seconds)
        timelist.append(time1)#save time
        index.append(fposition)#sasve index
    
    if(len(rr_s)<=15):
        pass
    else:
        ppgresults0=hrv_analysis(rr_s)
        ppgresults0=pd.DataFrame([ppgresults0])
        hrvresults=pd.concat([hrvresults,ppgresults0],axis=0)
        #hrvresults.to_csv(path+"/ppg+mot_"+userid+"/ppg/"+date+"/"+"hrvanalysis"+"/"+"upsampled300_"+date+"_rejectionmodified.csv")
        #timelist.append(time1)
        #index.append(fposition)
    fposition=fposition+1 
    
  

    #parameters.append("ecg")
    #parameters=pd.DataFrame(parameters)    
    #hrvresults=pd.concat([hrvresults,ppgresults],axis=0)
    #ecgresults=pd.DataFrame([ecgresults])
    
    #ecgresults=pd.concat([ecgresults,ecgresults0],axis=0)
    #ecgresults.to_csv(path+"/ppg+mot_"+userid+"/ppg/"+date+"/"+"hrvanalysis"+"/"+"ecg_"+date+".csv")   
    

timelist=pd.DataFrame(timelist)
index=pd.DataFrame(index)
timelist=pd.concat([timelist,index],axis=1)

if runecg and savefile:
  ecghr_min.to_csv("D:/ECG data/data/ECGsignal_63267/ECGsignal/"+date+"/"+userid+"_"+date+"_ecg_6sec_minhr_withrejecitonrules.csv")
  ecghr_max.to_csv("D:/ECG data/data/ECGsignal_63267/ECGsignal/"+date+"/"+userid+"_"+date+"_ecg_6sec_maxhr_withrejectionrules.csv")
  ecghr_mean.to_csv("D:/ECG data/data/ECGsignal_63267/ECGsignal/"+date+"/"+userid+"_"+date+"_ecg_6sec_meanhr_withrejectionrules.csv")
  timelist.to_csv("D:/ECG data/data/ECGsignal_63267/ECGsignal/"+date+"/"+userid+"_"+date+"_ecg_time.csv")
  ecgresults.to_csv("D:/ECG data/data/ECGsignal_63267/ECGsignal/"+date+"/"+userid+"_"+date+"_hrvanalysis_withrejectionrules.csv");

#timelist.to_csv(path+"/ppg+mot_"+userid+"/ppg/"+date+"/"+"hrvanalysis"+"/"+'xxtime'+".csv")

#timelist.to_csv("D:/ECG data/data/ECGsignal_63267/ECGsignal/"+date+"/"+userid+"/"+userid+"_"+date+"_ecgtime_remove_ectopicbeats.csv")
#hrvresults=pd.concat([hrvresults,timelist],axis=0
    
    #hrvresults=pd.concat([hrvresults,parameters],axis=1) 
    
    #results=pd.DataFrame(results)
    #results=pd.concat([parameters,results],axis=1)
    #results.to_csv(path+"ppg+mot_"+userid+"/ppg/"+date+"/"+winlen+"/1/"+date+"_"+time+"_"+".csv")     #保存ppg interval的数据

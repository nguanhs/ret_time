import numpy as np
from numpy.linalg import inv
import pandas as pd
from scipy.signal import lfilter,savgol_filter
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt
import src.ret_tim_func as ret
from src.raw_xml import rawdata
from src.mzml_xml import rawmzml
from src.read_datb_class import database
from src.lodes2 import lodes2
import src.use_func as usf
import os
import sys

'''
## Ladder
ladder='./YTL/ladder.mzML'
ladder_format='mzML'
if ladder_format=='mzML':
    ladd=rawmzml(ladder)
#raw.rawtomzml()
allraw=ladd.readxml()
print(list(allraw[0].keys()))
print(len(allraw))
cyc=ladd.check_cyc(allraw)
target_mzs=[527.5, 689.5, 851.5, 1013.5, 1175.5, 1337.5, 1499.5]
time1,intent1=ladd.ladder_ret(allraw,target_mzs)
ptime=[]
for i in range(len(target_mzs)):
#    plt.plot(time1[i],intent1[i])
    max1=ret.in_peaks(time1[i],intent1[i])
    ptime.append(max1[0]['time'])
    ptime.append(max1[1]['time'])

gu_index=sorted(ptime)
print(gu_index)
#plt.xticks(ptime1)
#print(ptime1)
#plt.show()
#exit()

gu_index=[13.329506666667, 14.961555, 19.092311666667, 21.58939, 25.81464, 28.503386666667, 31.962108333333, 34.555823333333, 36.764713333333, 39.165001666667, 41.278675, 43.29519, 45.217025, 47.042606666667]    

print('a')
'''
gu_index=[13.329506666667, 14.961555, 19.092311666667, 21.58939, 25.81464, 28.503386666667, 31.962108333333, 34.555823333333, 36.764713333333, 39.165001666667, 41.278675, 43.29519, 45.217025, 47.042606666667]

high_mannose={609:'MAN1',771:'MAN2',933:'MAN3',1095:'MAN4',1257:'MAN5',1419:'MAN6'
              ,1581:'MAN7',1743:'MAN8',1905:'MAN9',1045:'MAN10'}
#mass_tol=0.5

## Raw or mzML
def raw_input(fraw):
    if fraw[len(fraw)-3:]=='raw':
        convertsh='convert_rawtomzml.sh'
        fmzml=fraw[:len(fraw)-3]+'mzML'
        raw=rawdata(fraw,fmzml,convertsh)
    #    print(fmzml)
        raw.rawtomzml()
    if fraw[len(fraw)-4:]=='mzML':
        fmzml=fraw
        convertsh='convert_rawtomzml.sh'
        raw=rawdata(fraw,fmzml,convertsh)   
    allraw=raw.readxml()
    return allraw

## Determine high mannose masses
def raw_hman_info(allraw):
    i=0
    hman=[]
    for a in allraw:
        scan_number=a.get('scan')
        level=a.get('ms')
        stime=a.get('scan start time')
        tandem=a.get('Tandem_mass')
#    print([scan_number,level,stime,tandem])
        if level=='2':
            hman.append([int(tandem[0]),high_mannose.get(int(tandem[0]))])
        if i!=0 and level=='1':
            n_spect=i
            break
        i=i+1
    return n_spect,hman

##Determine matching and similarity scores
def match_cal(allraw,hman,method,mass_tol):
    amatch=[];all_results=[]
    for a in hman:
        #read database
    #    a=[1743,'MAN8']
    #    a=[1045, 'MAN10']
    #    a=[1581, 'MAN7']
        mass=a[0]
        nglycan=a[1]
        inputexcel=usf.database_excel(nglycan)
        dat=database(inputexcel,nglycan)
        datab=dat.data()

        
        ## retention time according to mass
    #    print(mass)
    #    a,b=ret.ms2ret_mass(allraw,mass)
    #    plt.plot(a,b)
    #    plt.show()
    #    exit()
        match,match1=ret.match_data_input2(allraw,mass,datab)
        #print(match1)
        # compare with database
        results1=[]
        match=match1
        amatch.append(match)
        if len(match)==0:
            print('Does not found any match in the retention time.')
        test=lodes2(nglycan)    
        for i in range(len(match)):
            peak_index=[match[i][3],match[i][4]]
            for ret_time_num in peak_index:
                rawms_spectra=test.cal_score2(allraw,mass,datab,ret_time_num)
    #            print(rawms_spectra[0]['ms_level'],rawms_spectra[2]['ms_level'])
    #            for c in rawms_spectra:
    #                plt.plot(c['spectra'][0],c['spectra'][1])
    #                plt.title(c['Ret_peak'])
    #                print(c['Ret_peak'])
    #                plt.show()
    #            plt.plot(rawms_spectra[0]['spectra'][0],rawms_spectra[0]['spectra'][1])
    #            plt.show()
    #            exit()
                #print(list(rawms_spectra[0].keys()))
                #print(match[i][0])
                #rawms_spectra
                iso_list=match[i][0]
                scores1=test.ms_sim4(iso_list,rawms_spectra,datab,method,mass_tol)
                scores1['Ret_peak']=rawms_spectra[0].get('Ret_peak')
                results1.append(scores1)
        all_results.append(results1)
    return amatch,all_results

def present_results(hman,amatch,all_results,allraw,method,raw_file):
    i=0
    for a in hman:
        match=amatch[i]
        results1=all_results[i]
        i=i+1
        if len(match)!=0:
            mass=a[0]
            nglycan=a[1]
            test=lodes2(nglycan)
            pict=test.isomer_picture()
            prefix=raw_file.split('.')[0]
            fig_name=prefix+'_'+nglycan+'_'+str(mass)+'_'+str(method)+'.png'
            table_name=prefix+'_'+nglycan+'_'+str(mass)+'_table_1.png'
            usf.plot_res3(match,results1,allraw,mass,pict,gu_index,fig_name,method)
        #usf.spect_table(match,results1,table_name)






#allraw=raw_input('Blackbean_M7.raw')
'''
exit()

#ms_level,targeted_mass=raw.ms_lev(allraw)
#tim1,int1=raw.ms2ret(allraw)

##Check raw's information
# Determine tandem mass cycle
high_mannose={609:'MAN1',771:'MAN2',933:'MAN3',1095:'MAN4',1257:'MAN5',1419:'MAN6'
              ,1581:'MAN7',1743:'MAN8',1905:'MAN9',1045:'MAN10'}
i=0
hman=[]
for a in allraw:
    scan_number=a.get('scan')
    level=a.get('ms')
    stime=a.get('scan start time')
    tandem=a.get('Tandem_mass')
#    print([scan_number,level,stime,tandem])
    if level=='2':
        hman.append([int(tandem[0]),high_mannose.get(int(tandem[0]))])
    if i!=0 and level=='1':
        n_spect=i
        break
    i=i+1
    
print(n_spect,hman)

for a in hman:
    #read database
#    a=[1743,'MAN8']
#    a=[1045, 'MAN10']
#    a=[1581, 'MAN7']
    mass=a[0]
    nglycan=a[1]
    inputexcel=usf.database_excel(nglycan)
    dat=database(inputexcel,nglycan)
    datab=dat.data()   

    
    ## retention time according to mass
#    print(mass)
#    a,b=ret.ms2ret_mass(allraw,mass)
#    plt.plot(a,b)
#    plt.show()
#    exit()
    match,match1=ret.match_data_input2(allraw,mass,datab)
    print(match1)
    # compare with database
    results1=[]
    match=match1
    if len(match)==0:
        print('Does not found any match in the retention time.')
        break
    for i in range(len(match)):
        peak_index=[match[i][3],match[i][4]]
        for ret_time_num in peak_index:
            test=lodes2(nglycan)
            rawms_spectra=test.cal_score2(allraw,mass,datab,ret_time_num)
#            print(rawms_spectra[0]['ms_level'],rawms_spectra[2]['ms_level'])
#            for c in rawms_spectra:
#                plt.plot(c['spectra'][0],c['spectra'][1])
#                plt.title(c['Ret_peak'])
#                print(c['Ret_peak'])
#                plt.show()
#            plt.plot(rawms_spectra[0]['spectra'][0],rawms_spectra[0]['spectra'][1])
#            plt.show()
#            exit()
            #print(list(rawms_spectra[0].keys()))
            #print(match[i][0])
            #rawms_spectra
            iso_list=match[i][0]
            scores1=test.ms_sim4(iso_list,rawms_spectra,datab,3)
            scores1['Ret_peak']=rawms_spectra[0].get('Ret_peak')
            results1.append(scores1)
#    for a in results1:
#        print(a)
    pict=test.isomer_picture()
    fig_name=nglycan+'_'+str(mass)+'_3.png'
    table_name=nglycan+'_'+str(mass)+'_table_1.png'
    usf.plot_res3(match,results1,allraw,mass,pict,gu_index,fig_name)
    usf.spect_table(match,results1,table_name)
#    plt.show()
#    exit()
'''

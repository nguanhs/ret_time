import numpy as np
import pandas as pd
import csv
import base64
import zlib
import os
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
import src.ret_tim_func as ret

    
def isotpemass():
    atm={'C':12.000000,'H':1.007825,'O':15.994915,'N':14.003074,'Na':22.989769}
    hex=[('C',6),('H',12),('O',6)]
    hexnac=[('C',8),('H',15),('O',6),('N',1)] 
    water=atm['O']+2*atm['H']
    sod=atm['Na']
    mass=0
    for a in hex:
        mass=mass+atm[a[0]]*a[1]
    mhex=mass
    print(mass,mass+atm['Na'])
    mass=0
    for a in hexnac:
        mass=mass+atm[a[0]]*a[1]
    print(mass,mass+atm['Na'])
    mhexnac=mass
    hman=['Man1','Man2','Man3','Man4','Man5','Man6','Man7','Man8','Man9','Man10']
    ihman_iso={}
    for i in range(10):
        j=i+1
        print(j*mhex+2*mhexnac+sod-(j+1)*water,j*mhex+2*mhexnac-(j+1)*water)
        if j!=10:
            ihman_iso[hman[i]]=j*mhex+2*mhexnac+sod-(j+1)*water
        elif j==10:
            ihman_iso[hman[i]]=(j*mhex+2*mhexnac+2*sod-(j+1)*water)/2
    de1=[18,60,90,120,101,162,203,221,324,424,486] 
    de1_com={18:[('H',2),('O',1)],
             46:[('C',2),('H',6),('O',1)], #H3C-CH2(OH)
             72:[('C',3),('H',4),('O',2)], 
             60:[('C',2),('H',4),('O',2)],
             90:[('C',3),('H',6),('O',3)],
             120:[('C',4),('H',8),('O',4)],
             101:[('C',4),('H',7),('O',2),('N',1)],
             203:[('C',8),('H',13),('O',5),('N',1)],
             144:[('C',6),('H',8),('O',4)], #dehyd_dehyd_hex
             162:[('C',6),('H',10),('O',5)], #dehyd_hex
             180:[('C',6),('H',12),('O',6)], #hex
             221:[('C',8),('H',15),('O',6),('N',1)], # hexnac
             252:[('C',6),('H',10),('O',5),('C',3),('H',6),('O',3)], # hex+hex_cr90
             324:[('C',12),('H',20),('O',10)], #hex+hex_dehyd
             396:[('C',12),('H',20),('O',10),('C',3),('H',4),('O',2)], #hex+hex_dehyd+72
             342:[('C',12),('H',22),('O',11)], #2hex
             414:[('C',12),('H',22),('O',11),('C',3),('H',4),('O',2)],#2hex+hex_cr72
             486:[('C',18),('H',30),('O',15)], #hex+hex+hex_dehyd
             528:[('C',18),('H',30),('O',15),('C',2),('H',2),('O',1)], #hex+hex+hex_dehyd+C2OH2 (+42)
             558:[('C',18),('H',30),('O',15),('C',3),('H',4),('O',2)], #hex+hex+hex_dehyd+72
             588:[('C',18),('H',30),('O',15),('C',4),('H',6),('O',3)], #hex+hex+hex_dehyd+102
             504:[('C',18),('H',32),('O',16)], #hex+hex+hex
             576:[('C',18),('H',32),('O',16),('C',3),('H',4),('O',2)], #hex+hex+hex+cr90
             424:[('C',16),('H',28),('O',11),('N',2)], #hexnac+hexnac
             406:[('C',16),('H',26),('O',10),('N',2)], #hexnac+hexnac_dehyd
             545:[('C',6),('H',10),('O',5),('C',6),('H',10),('O',5),('C',8),('H',15),('O',6),('N',1)], #2*hex_dehyd+hexnac
             707:[('C',18),('H',30),('O',15),('C',8),('H',15),('O',6),('N',1)], #3*hex+hexnac
             630:[('C',24),('H',38),('O',19)], #hex+hex+hex+hex_dehyd_dehyd
             648:[('C',24),('H',40),('O',20)], #hex+hex+hex+hex_dehyd
             738:[('C',24),('H',40),('O',20),('C',3),('H',6),('O',3)], #hex+hex+hex+hex_dehyd+cr90
             720:[('C',24),('H',40),('O',20),('C',3),('H',4),('O',2)], #hex+hex+hex+hex_dehyd+72
             568:[('C',6),('H',10),('O',5),('C',16),('H',26),('O',10),('N',2)], # hex+hexnac+hexnac_dehyd 
             586:[('C',6),('H',10),('O',5),('C',16),('H',28),('O',11),('N',2)], # hex+hexnac+hexnac
             730:[('C',12),('H',20),('O',10),('C',16),('H',26),('O',10),('N',2)], # 2hex+hexnac+hexnac_dehyd
             748:[('C',12),('H',20),('O',10),('C',16),('H',28),('O',11),('N',2)], # 2hex+hexnac+hexnac
             810:[('C',30),('H',50),('O',25)], #5hex_dehyd
             869:[('C',24),('H',40),('O',20),('C',8),('H',15),('O',6),('N',1)], #4*hex+hexnac
             892:[('C',18),('H',30),('O',15),('C',16),('H',26),('O',10),('N',2)], # 3hex+hexnac+hexnac_dehyd 
             
             
             
             0:[0]}
             
    frag=[203,226,244,347,365,509,527,671,689]	
    frag_com={203:[('C',6),('H',12),('O',6),('Na',1)], #hex+Na
              226:[('C',8),('H',13),('O',5),('N',1),('Na',1)], #dehyd_hexnac+Na
              244:[('C',8),('H',15),('O',6),('N',1),('Na',1)], #hexnac+Na
              347:[('C',12),('H',20),('O',10),('Na',1)], #dehyd_hex-hex+Na
              365:[('C',12),('H',22),('O',11),('Na',1)], #hex-hex+Na
              437:[('C',15),('H',26),('O',13),('Na',1)], #hex+hex+hex_cr+Na
              447:[('C',16),('H',28),('O',11),('N',2),('Na',1)], #hexnac+hexnac+Na
              509:[('C',18),('H',30),('O',15),('Na',1)], #hex+hex+hex_dehyd+Na
              527:[('C',18),('H',32),('O',16),('Na',1)], #hex+hex+hex+Na
              671:[('C',24),('H',40),('O',20),('Na',1)], #hex+hex+hex+hex_dehyd+Na
              689:[('C',24),('H',42),('O',21),('Na',1)], #hex+hex+hex+hex+Na
              712:[('C',18),('H',30),('O',15),('C',8),('H',13),('O',5),('N',1),('Na',1)], #3hex+hexnac_dehyd+Na
              771:[('C',12),('H',20),('O',10),('C',16),('H',28),('O',11),('N',2),('Na',1)], # 2hex+hexnac+hexnac+Na
              833:[('C',30),('H',50),('O',25),('Na',1)], #4hex+hex_dehyd+Na
              851:[('C',30),('H',52),('O',26),('Na',1)], #5hex+Na
              
              874:[('C',24),('H',40),('O',20),('C',8),('H',13),('O',5),('N',1),('Na',1)], #4hex+hexnac_dehyd+Na
              995:[('C',36),('H',60),('O',30),('Na',1)], #6hex_dehyd+Na
              1013:[('C',36),('H',62),('O',31),('Na',1)], #6hex+Na
              1157:[('C',42),('H',70),('O',35),('Na',1)], #7hex_dehyd+Na
              1175:[('C',42),('H',72),('O',36),('Na',1)], #7hex+Na
              1319:[('C',48),('H',80),('O',40),('Na',1)], #8hex_dehyd+Na
              1337:[('C',48),('H',82),('O',41),('Na',1)], #8hex+Na
              1481:[('C',54),('H',90),('O',45),('Na',1)], #9hex_dehyd+Na
              1499:[('C',54),('H',92),('O',46),('Na',1)], #9hex+Na
              1643:[('C',60),('H',100),('O',50),('Na',1)], #10hex_dehyd+Na
              1661:[('C',60),('H',102),('O',51),('Na',1)], #10hex+Na

              0:[]}

    print(cal_isomass(de1_com[545]))
    print(cal_isomass(de1_com[648]))
    print(cal_isomass(de1_com[46]))
    print(cal_isomass(de1_com[72]))
    print(cal_isomass(frag_com[509]))
    print(cal_isomass(frag_com[1013]))


    im=ihman_iso['Man10']
    im1=cal_isomass(frag_com[1661])
    im2=cal_isomass(frag_com[1643])
#    im3=cal_isomass(frag_com[1319])
#    print([round(im,3),round(im1,3),round(im2,3),round(im3,3)])
    print([round(im,3),round(im1,3),round(im2,3)])  


    a=[]
    nm=1045
    im=ihman_iso['Man10']
    sel=[1036,994.5,964,943.5,934.5,883,802,833]
    for b in sel:
        de2=nm-b
        #a.append(round(im-cal_isomass(de1_com[de2]),3))
#        print(de2,cal_isomass(de1_com[de2*2])/2,cal_isomass(de1_com[de2*2]))
        a.append(round(im-cal_isomass(de1_com[de2*2])/2,3))
    print(a)

    a=[]
    nm=1045
    im=ihman_iso['Man10']
    sel=[802,721,691.5,671,640,610.5]
    for b in sel:
        de2=nm-b
        #a.append(round(im-cal_isomass(de1_com[de2]),3))
        a.append(round(im-cal_isomass(de1_com[de2*2])/2,3)) 
    print(a)

    a=[]
    nm=1045
    im=ihman_iso['Man10']
    sel=[964,883,802,994.5]
    for b in sel:
        de2=nm-b
        #a.append(round(im-cal_isomass(de1_com[de2]),3))
        a.append(round(im-cal_isomass(de1_com[de2*2])/2,3))
    print(a)



    a=[]
    nm=1661
    im=cal_isomass(frag_com[nm])
    sel=[1643,1499,1337,1247,1175,1085,1013,923]
    for b in sel:
        de2=int(nm-b)
        a.append(round(im-cal_isomass(de1_com[de2]),3))
    print(a)


    a=[]
    nm=1643
    im=cal_isomass(frag_com[nm])
    sel=[1481,1319,1157,995]
    for b in sel:
        de2=int(nm-b)
        a.append(round(im-cal_isomass(de1_com[de2]),3))
    print(a)

    exit()
    a=[]
    nm=527
    im=cal_isomass(frag_com[nm])
    sel=[509,467,437,407,365,347,275]
    for b in sel:
        de2=int(nm-b)
        a.append(round(im-cal_isomass(de1_com[de2]),3))
    print(a)

##[1581,1157,1175,527]

 
 
def cal_isomass(at):
    atm={'C':12.000000,'H':1.007825,'O':15.994915,'N':14.003074,'Na':22.989769}
    mass=0
    for a in at:
        mass=mass+atm[a[0]]*a[1]
    return mass
    

def get_msspectra(dat,isomer,pref):
    for a in dat:
        if a['isomer']==isomer:
            ma=a[pref+'_mass']
            inten=a[pref+'_intensity']
            break
    return ma,inten

def area(x,y,x0):
    x1=x0-0.5;x2=x0+0.5
    xx1=[];yy1=[]
    for xx in x:
        if xx >=x1 and xx <=x2:
            xx1.append(xx)
            i=list(x).index(xx)
            yy1.append(y[i])
    area=0
    for i in range(len(xx1)-1):
        area=area+(xx1[i+1]-xx1[i])*(yy1[i+1]+yy1[i])/2.0
    return area
    
def t_inter(retention_peaks,ret_int1,ret_time_num):
    t1=retention_peaks[ret_time_num].get('time')
    #h1=test1[i].get('max_intens')/2.0
    h1=0.05
    ih=int(retention_peaks[ret_time_num].get('time_i'))
    r=[]
    t=[]
    for j in range(100):
        i1=ih-j
        if ret_int1[1][i1] < h1:
            t_ini=ret_int1[0][i1+1]
            break
    for j in range(100):
        i1=ih+j+1
        if  ret_int1[1][i1] < h1:
            t_end=ret_int1[0][i1-1]
            break
    return t_ini,t_end

def spectra_profile(dat,t_ini,t_end,ms_sel):
    mi=[];inte=[];mrange=[];intensity=[]
    for a in dat:
        tim=float(a['scan start time'])
        if tim >= t_ini and tim <= t_end:
            try:
                mass=int(float(a['selected mass']))
            except:
                mass='no'
            if mass==ms_sel:
                b1=a['masses']
                decoded_source = base64.b64decode(b1.encode('ascii'))
                decoded_source = zlib.decompress(decoded_source)
                output = np.frombuffer(bytearray(decoded_source), dtype=np.float64)
                mi.append(output)

                b1=a['intensity']
                decoded_source = base64.b64decode(b1.encode('ascii'))
                decoded_source = zlib.decompress(decoded_source)
                output = np.frombuffer(bytearray(decoded_source), dtype=np.float64)
                inte.append(output)
    intent_com=[]
    mass_com=[]
    arr=0
    for i in range(len(mi)):
        arr=np.add(arr,inte[i])
    mass_com=mi[0]
    intent_com=arr
    maxi=max(intent_com)
    for i in range(len(intent_com)):
        intent_com[i]=intent_com[i]/maxi
    return mass_com,intent_com,maxi

def spectra_profile1(dat,mass_ms2,t_ini,t_end,ms_sel):
    mi=[];inte=[];mrange=[];intensity=[]
    for a in dat:
        tim=float(a['scan start time'])
        
        if tim >= t_ini and tim <= t_end:
            try:
                mass=int(float(a['selected mass']))
            except:
                mass='no'
            #print(mass,a['Tandem_mass'][0])
            if mass==ms_sel and int(a['Tandem_mass'][0])==mass_ms2:
                b1=a['masses']
                decoded_source = base64.b64decode(b1.encode('ascii'))
                decoded_source = zlib.decompress(decoded_source)
                output = np.frombuffer(bytearray(decoded_source), dtype=np.float64)
                mi.append(output)

                b1=a['intensity']
                decoded_source = base64.b64decode(b1.encode('ascii'))
                decoded_source = zlib.decompress(decoded_source)
                output = np.frombuffer(bytearray(decoded_source), dtype=np.float64)
                inte.append(output)
    intent_com=[]
    mass_com=[]
    arr=0
    for i in range(len(mi)):
        arr=np.add(arr,inte[i])
    mass_com=mi[0]
    intent_com=arr
    maxi=max(intent_com)
    for i in range(len(intent_com)):
        intent_com[i]=intent_com[i]/maxi
    return mass_com,intent_com,maxi


def plot_res(match,results,tim1,int1,max_tint1,ms_level,targeted_mass,pict):
    fig, axs = plt.subplots(2,figsize=(10,7))
    max_tint1.sort(key=lambda x: x.get('time'),reverse=False)
#    fig, axs = plt.subplots(2)
    axs[1].axis('off')
    linew=0.5
    ti1=[]
    for i in range(len(max_tint1)):
        ti1.append(max_tint1[i].get('time'))
    axs[0].set_xticks(ti1)
    axs[0].plot(tim1,int1)
    axs[0].set_xlim([ti1[0]-3, ti1[len(ti1)-1]+3])
    axs[0].set_ylim([0,1.8])
    iso_m=[];ms_score1=[]
#    max_tint1.sort(key=lambda x: x.get('time'),reverse=False)
    k1=0
    for i in range(len(max_tint1)):
        iso=[]
        sco=[]
        axs[1].text(0, 0.95-k1*linew, 'Time:'+str(round(max_tint1[i].get('time'),1)))
        k1=k1+1
        j1=0
        for a in results:
            if a['Ret_peak']==max_tint1[i].get('time'):
                j1=j1+1
                iso.append(a['isomer'])
                sco.append(a['MSn_score'])
                l1=' '
                j=0
                for b in a['MSn_score']:
                    l1=l1+',MS'+str(ms_level[j])+'_'+str(targeted_mass[j])+':'+str(round(b,2))+' '
                    j=j+1
                l1=a['isomer']+l1
                axs[1].text(0, 0.95-k1*linew, l1)
                k1=k1+1
        #im = axs[0].imread('bird.jpg')
        
        axs[0].text(max_tint1[i].get('time')-0.5,max_tint1[i].get('max_intens')+0.2,iso )
        k1=k1+1 
    j1=0
    for a in match:
        pic=pict[a[0]]
  
#        pic=d1+a['isomer']+'.png'
        img = plt.imread(pic)
        im = OffsetImage(img,zoom=0.08)
        ab = AnnotationBbox(im, (0.8,0.6-j1*0.5), xycoords='axes fraction', box_alignment=(1.1,-0.1))
        axs[1].add_artist(ab)
        j1=j1+1
    plt.show

def plot_res1(match,results,allraw,pict):

    tim1,int1=ret.ms2ret(allraw)
    max_tint1=ret.in_peaks2(tim1,int1,0.1)
    ms_level,targeted_mass=ret.ms_lev(allraw)
    
    fig, axs = plt.subplots(2,figsize=(10,14))
#    ax1 = plt.subplot2grid((4, 5), (0, 0), colspan=5)
#    ax2 = plt.subplot2grid((4, 5), (1, 0), colspan=2, rowspan=3)
#    ax3=[]

    ax1 = plt.subplot2grid((7, 5), (0, 0), colspan=5)
    ax2 = plt.subplot2grid((7, 5), (1, 0), colspan=5, rowspan=2)
    ax3=[]
    
    ax1.set_xlabel('Time(Minutes)',fontsize=15)
    ax1.set_title('LC Retention time',fontsize=15)
    ax2.set_title('MSn similarity score',fontsize=15)
    l3=0
    for l1 in range(4):
        for l2 in range(5):
            ax3a = plt.subplot2grid((7, 5), (3+l1, 0+l2))
            ax3a.axis('off')
            ax3.append(ax3a)

    max_tint1.sort(key=lambda x: x.get('time'),reverse=False)
#    fig, axs = plt.subplots(2)
    #axs[1].axis('off')
    ax2.axis('off')
    linew=0.05
    ti1=[]
    d1='./image_man/'
#    pict=[d1+'7E1.png', d1+'7D3.png', d1+'7G1.png', d1+'7E2.png', d1+'7D2.png', d1+'7D1.png', d1+'7F1.png']
    for i in range(len(max_tint1)):
        ti1.append(max_tint1[i].get('time'))
    ''''
    axs[0].set_xticks(ti1)
    axs[0].plot(tim1,int1)
    axs[0].set_xlim([ti1[0]-3, ti1[len(ti1)-1]+3])
    axs[0].set_ylim([0,1.8])
    '''
    ax1.set_xticks(ti1)
    ax1.set_xlim([ti1[0]-3, ti1[len(ti1)-1]+3])
    ax1.set_ylim([0,1.8])
    ax1.plot(tim1,int1)
    iso_m=[];ms_score1=[]
#    max_tint1.sort(key=lambda x: x.get('time'),reverse=False)
    #score_name=[]
    #j=0
    #for b in results[]['MSn_score']:
    #    l1=l1+',MS'+str(ms_level[j])+'_'+str(targeted_mass[j])+': '+str(round(b,2))+' '
    #                j=j+1
    #for a in
    k1=0
    for i in range(len(max_tint1)):
        iso=[]
        sco=[]
        ax2.text(-0.3, 0.95-k1*linew, 'Peak at time = '+str(round(max_tint1[i].get('time'),1))+'min')
        k1=k1+1
        j1=0
        for a in results:
            if a['Ret_peak']==max_tint1[i].get('time'):
                j1=j1+1
                iso.append(a['isomer'])
                sco.append(a['MSn_score'])
                l1=' '
                j=0
                for b in a['MSn_score']:
                    l1=l1+','+b[0]+': '+str(round(b[1],2))+' '

                l1=a['isomer']+l1
                ax2.text(-0.2, 0.95-k1*linew, l1,)
                k1=k1+1
        #im = axs[0].imread('bird.jpg')
        
        ax1.text(max_tint1[i].get('time')-0.5,max_tint1[i].get('max_intens')+0.2,iso )
        k1=k1+1
    j1=0
    for a in match:
        pic=pict[a[0]]
        ax3[j1].axis('off')
#        pic=d1+a['isomer']+'.png'
        img = plt.imread(pic)
        im = OffsetImage(img,zoom=0.05)
        #ab = AnnotationBbox(im, (0.8,0.6-j1*0.5), xycoords='axes fraction', box_alignment=(1.1,-0.1))
         #ab = AnnotationBbox(im, (0.0,0.0), xycoords='axes fraction')
        #ax3[j1].add_artist(ab)
        ax3[j1].imshow(img)
        j1=j1+1
    plt.tight_layout()
    plt.show

def plot_res2(match,results,allraw,pict):

    tim1,int1=ret.ms2ret(allraw)
    max_tint1=ret.in_peaks2(tim1,int1,0.1)
    ms_level,targeted_mass=ret.ms_lev(allraw)
    
    fig, axs = plt.subplots(2,figsize=(10,14))

    ax1 = plt.subplot2grid((7, 5), (0, 0), colspan=5)
    ax2 = plt.subplot2grid((7, 5), (1, 0), colspan=5, rowspan=2)
    ax3=[]
    
    ax1.set_xlabel('Time(Minutes)',fontsize=15)
    ax1.set_title('LC Retention time',fontsize=15)
    ax2.set_title('MSn similarity score',fontsize=15)
    l3=0
    for l1 in range(4):
        for l2 in range(5):
            ax3a = plt.subplot2grid((7, 5), (3+l1, 0+l2))
            ax3a.axis('off')
            ax3.append(ax3a)

    max_tint1.sort(key=lambda x: x.get('time'),reverse=False)
    ax2.axis('off')
    linew=0.05
    ti1=[]
    d1='./image_man/'
    for i in range(len(max_tint1)):
        ti1.append(max_tint1[i].get('time'))

    ax1.set_xticks(ti1)
    ax1.set_xlim([ti1[0]-3, ti1[len(ti1)-1]+3])
    ax1.set_ylim([0,1.8])
    ax1.plot(tim1,int1)
    iso_m=[];ms_score1=[]

    k1=0
    for i in range(len(max_tint1)):
        iso=[]
        sco=[]
        ax2.text(-0.3, 0.95-k1*linew, 'Peak at time = '+str(round(max_tint1[i].get('time'),1))+'min')
        k1=k1+1
        j1=0
        iso_name=[]
        for a in results:
            if a['Ret_peak']==max_tint1[i].get('time'):
                j1=j1+1
                go1=0
                
                for is1 in iso_name:
                    if is1==a['isomer']:
                        go1=1
                        break
                
                if go1==0:
                    iso_name.append(a['isomer'])
                    iso.append(a['isomer'])
                    sco.append(a['MSn_score'])
                    l1=' '
                    j=0
                    for b in a['MSn_score']:
                        l1=l1+','+b[0]+': '+str(round(b[1],2))+' '

                    l1=a['isomer']+l1
                    ax2.text(-0.2, 0.95-k1*linew, l1,)
                    k1=k1+1
                
        #im = axs[0].imread('bird.jpg')
        
        ax1.text(max_tint1[i].get('time')-0.5,max_tint1[i].get('max_intens')+0.2,iso )
        k1=k1+1
    j1=0
    for a in match:
        pic=pict[a[0]]
        ax3[j1].axis('off')
#        pic=d1+a['isomer']+'.png'
        img = plt.imread(pic)
        im = OffsetImage(img,zoom=0.05)
        #ab = AnnotationBbox(im, (0.8,0.6-j1*0.5), xycoords='axes fraction', box_alignment=(1.1,-0.1))
         #ab = AnnotationBbox(im, (0.0,0.0), xycoords='axes fraction')
        #ax3[j1].add_artist(ab)
        ax3[j1].imshow(img)
        j1=j1+1
    plt.tight_layout()
    plt.show()

def plot_res3(match,results,allraw,mass_ms2,pict,ladder_ret,fig_name,method):
   
    cmap1 = mcolors.LinearSegmentedColormap.from_list("my_cmap1", ["red", "orange", "green","blue","purple"])

    tim1,int1=ret.ms2ret_mass(allraw,mass_ms2)
    max_tint1=ret.in_peaks2(tim1,int1,0.1)
    ms_level,targeted_mass=ret.ms_lev_mass(allraw,mass_ms2)
    
    fig, axs = plt.subplots(2,figsize=(10,14))

    ax1 = plt.subplot2grid((7, 5), (0, 0), colspan=5)
    ax2 = plt.subplot2grid((7, 5), (1, 0), colspan=5, rowspan=3)
    ax3=[]
    #ax1.set_position([0.0, ax1.get_position().y0, ax1.get_position().width, ax1.get_position().height])
    #ax1.set_position([0.0, 0.0, 0.8, 0.8])
    #ax1.set_xlabel('Time(Minutes)',fontsize=15)
    ax1.set_xlabel('Gu unit',fontsize=15)
    ax1.set_title('LC Retention time',fontsize=15)
    if method==1:
        ax2.set_title('MSn similarity score (integrated difference)',fontsize=15,loc='left')
    if method==2:
        ax2.set_title('MSn similarity score (dot product)',fontsize=15,loc='left')
    if method==3:
        ax2.set_title('MSn similarity score (KL divergence)',fontsize=15,loc='left')
    l3=0
    for l1 in range(3):
        for l2 in range(5):
            ax3a = plt.subplot2grid((7, 5), (4+l1, 0+l2))
            ax3a.axis('off')
            ax3.append(ax3a)

    max_tint1.sort(key=lambda x: x.get('time'),reverse=False)
    ax2.axis('off')
    linew=0.03
    ti1=[]
    d1='./image_man/'
    for i in range(len(max_tint1)):
        ti2=max_tint1[i].get('time')
        #print(ti2)
        if ti2 > ladder_ret[0] and ti2 < ladder_ret[len(ladder_ret)-1]:
            ti3=interpolate(ti2, ladder_ret)
        if ti2 > ladder_ret[len(ladder_ret)-1]:
            ti3=len(ladder_ret)+(ti2-ladder_ret[len(ladder_ret)-1])
        if ti2 < ladder_ret[0]:
            ti3=ti2/ladder_ret[0]
        ti3=round(ti3,1)
        #ti1.append(max_tint1[i].get('time'))
        ti1.append(ti3)
    
    larg=[]
    for a in match:
        ti2=float(a[5])
        if ti2 > ladder_ret[0] and ti2 < ladder_ret[len(ladder_ret)-1]:
            ti3=interpolate(ti2, ladder_ret)
        if ti2 > ladder_ret[len(ladder_ret)-1]:
            ti3=len(ladder_ret)+(ti2-ladder_ret[len(ladder_ret)-1])
        if ti2 < ladder_ret[0]:
            ti3=ti2/ladder_ret[0]
        larg.append(ti3)

        ti2=float(a[6])
        if ti2 > ladder_ret[0] and ti2 < ladder_ret[len(ladder_ret)-1]:
            ti3=interpolate(ti2, ladder_ret)
        if ti2 > ladder_ret[len(ladder_ret)-1]:
            ti3=len(ladder_ret)+(ti2-ladder_ret[len(ladder_ret)-1])
        if ti2 < ladder_ret[0]:
            ti3=ti2/ladder_ret[0]
        larg.append(ti3)
        
    large=max(larg)
    small=min(larg)
    xtic=np.arange(int(ti1[0]-1), int(ti1[len(ti1)-1]+1))
    ax1.set_xticks(xtic)
    #ax1.set_xlim([ti1[0]-1, ti1[len(ti1)-1]+1])
    print(larg)
    ax1.set_xlim(small-1, large+1)
    ax1.set_ylim([0,1.8])
    tim2=[]
    for a in tim1:
        if a > ladder_ret[0] and a < ladder_ret[len(ladder_ret)-1]:
            b=interpolate(a, ladder_ret)
        if a < ladder_ret[0]:
            b=a/ladder_ret[0]
        if a > ladder_ret[len(ladder_ret)-1]:
            b=len(ladder_ret)+(a-ladder_ret[len(ladder_ret)-1])
        tim2.append(b)
    ax1.plot(tim2,int1)
    iso_m=[];ms_score1=[]

    k1=0
    mlist=[]
    for i in range(len(max_tint1)):
        iso=[]
        sco=[]
        my=0
        for a in results:
            if a['Ret_peak']==max_tint1[i].get('time'):
                ax2.text(0.0, 0.95-k1*linew, 'Peak at time = '+str(round(max_tint1[i].get('time'),1))+'min, Gu_index='+str(round(ti1[i],1)))
                k1=k1+1
                my=1
                mlist.append(i)
                break
        j1=0
        iso_name=[]
        sc_list=[]
        sc_text=[]
        for a in results:
            if a['Ret_peak']==max_tint1[i].get('time'):
                j1=j1+1
                go1=0
                
                for is1 in iso_name:
                    if is1==a['isomer']:
                        go1=1
                        break
                
                if go1==0:
                    iso_name.append(a['isomer'])
                    iso.append(a['isomer'])
                    sco.append(a['MSn_score'])
                    l1=' '
                    j=0
                    aver=0
                    for b in a['MSn_score']:
                        l1=l1+','+b[0]+': '+str(round(b[1],2))+' '
                        aver=aver+b[1]/len(b)
                    l1=a['isomer']+l1
                    sc_list.append((aver,l1))
                    sc_text.append(l1)
                    #ax2.text(0.05, 0.95-k1*linew, l1,)
                    #k1=k1+1
        if method!=3:
            sc_list_sort = sorted(sc_list, key=lambda x: x[0], reverse=True)
        elif method==3:
            sc_list_sort = sorted(sc_list, key=lambda x: x[0])
        for b in sc_list_sort:
            ax2.text(0.05, 0.95-k1*linew, b[1])
            k1=k1+1
        #im = axs[0].imread('bird.jpg')
        
        #ax1.text(max_tint1[i].get('time')-0.5,max_tint1[i].get('max_intens')+0.2,iso )
        if len(mlist) > 1 and my==1:
            n1=mlist[len(mlist)-1]
            y1=max_tint1[n1].get('max_intens')
            y2=max_tint1[i].get('max_intens')
            if (ti1[i]-ti1[n1])< 1 and abs(y1-y2)<0.1:
                y_loc1=y_loc1+0.3
            else:
                
                y_loc1=max_tint1[i].get('max_intens')+0.4
        else:
            y_loc1=max_tint1[i].get('max_intens')+0.4
        #ax1.text(ti1[i]-0.5,max_tint1[i].get('max_intens')+0.2,iso )
        if my==1:
            color=i*1.0/len(ti1)
            ax1.text(ti1[i]-0.5,y_loc1,iso,color=cmap1(color))
            ax1.annotate(iso, xy=(ti1[i],max_tint1[i].get('max_intens')), xytext=(ti1[i]-0.5,y_loc1),arrowprops=dict(facecolor=cmap1(color), edgecolor=cmap1(color), arrowstyle='->'),color=cmap1(color))
            k1=k1+1
    j1=0
    for a in match:
        pic=pict[a[0]]
        ax3[j1].axis('off')
#        pic=d1+a['isomer']+'.png'
        img = plt.imread(pic)
        im = OffsetImage(img,zoom=0.05)
        #ab = AnnotationBbox(im, (0.8,0.6-j1*0.5), xycoords='axes fraction', box_alignment=(1.1,-0.1))
         #ab = AnnotationBbox(im, (0.0,0.0), xycoords='axes fraction')
        #ax3[j1].add_artist(ab)
        ax3[j1].imshow(img)
        j1=j1+1

    plt.tight_layout()
    plt.savefig(fig_name)
    #plt.show()
    

def interpolate(x, index_rt):
    sorted_index_rt = sorted(enumerate(index_rt, start=1), key=lambda x: x[1])
    for i in range(len(sorted_index_rt) - 1):
        if sorted_index_rt[i][1] <= x <= sorted_index_rt[i+1][1]:
            lower_rank = sorted_index_rt[i][0]
            lower_value = sorted_index_rt[i][1]
            upper_rank = sorted_index_rt[i+1][0]
            upper_value = sorted_index_rt[i+1][1]
            break
    interpolation_rank = lower_rank + (x - lower_value) / (upper_value - lower_value) * (upper_rank - lower_rank)
    return interpolation_rank

def area1(x,y,xrange):
    x1=[];y1=[]
    for i in range(len(x)):
        if x[i] >= xrange[0] and x[i] <= xrange[1]:
            x1.append(x[i])
            y1.append(y[i])
    ar=0
    for i in range(len(x1)-1):
        ar=ar+(y1[i]+y1[i+1])*(x1[i+1]-x1[i])/2.0
    return ar

def scr_1(spec_in,spec_dat,max1):
    arc=0
    are1=[]
    are2=[]
    for a in max1:
        #print(a['mass'])
        x1=a['mass']-1
        x2=a['mass']+1
        range=[x1,x2]
        
        ar1=area1(spec_in[0],spec_in[1],range)
        ar2=area1(spec_dat[0],spec_dat[1],range)
        are1.append(ar1)
        are2.append(ar2)
    ma1=max(are1)
    ma2=max(are2)
    #print(ma1,ma2)
    #arc=arc+abs(are1[i]/ma1-are2[i]/ma2)/len(max1)
    i=0
    for a in max1:
        arc=arc+abs(are1[i]/ma1-are2[i]/ma2)/len(max1)
        i=i+1
    tabl=[max1,are1,are2,ma1,ma2]
    return  1-arc,tabl
'''
def scr_2(spec_in,spec_dat,max1):
    arc=0;s1=0;s2=0;arc1=0;arc2=0
    for a in max1:
        #print(a['mass'])
        x1=a['mass']-1
        x2=a['mass']+1
        range=[x1,x2] 
        ar1=area1(spec_in[0],spec_in[1],range)
        ar2=area1(spec_dat[0],spec_dat[1],range) 
        s1=s1+ar1;s2=s2+ar2
        #arc=arc+abs(ar1-ar2)/len(max1)
        print(a['mass'],ar1,ar2)
        arc=arc+ar2*np.log(ar2/ar1)
   
        if ar1==0:
            arc1=arc1+0
        else:
            arc1=arc1+ar1*np.log(ar1/ar2)
        arc2=arc2+np.log(1.0/ar2)
    print(arc/s2+np.log(s1/s2),arc1/s1+np.log(s2/s1),arc2/(len(max1))+np.log(s2/len(max1)))
'''

def scr_2(spec_in,spec_dat,max1):
    arc=0
    are1=[]
    are2=[]
    for a in max1:
        #print(a['mass'])
        x1=a['mass']-1
        x2=a['mass']+1
        range=[x1,x2]

        ar1=area1(spec_in[0],spec_in[1],range)
        ar2=area1(spec_dat[0],spec_dat[1],range)
        are1.append(ar1)
        are2.append(ar2)
    ma1=max(are1)
    ma2=max(are2)
    i=0
    tabl=[max1,are1,are2,ma1,ma2]
    for a in are1:
        are1[i]=are1[i]/ma1
        are2[i]=are2[i]/ma2
        i=i+1
    #arc=arc+abs(are1[i]/ma1-are2[i]/ma2)/len(max1)
        #tabl=[max1,are1,are2,ma1,ma2]
    i=0
    le1=0
    le2=0
    for a in max1:
        #arc=arc+abs(are1[i]/ma1-are2[i]/ma2)/len(max1)
        arc=arc+abs(are1[i]*are2[i])
        le1=are1[i]*are1[i]+le1
        le2=are2[i]*are2[i]+le2
        i=i+1
    return  (arc**0.5/(le1**0.5*le2**0.5)),tabl


def scr_3(spec_in,spec_dat,max1):
    arc=0
    are1=[]
    are2=[]
    for a in max1:
        #print(a['mass'])
        x1=a['mass']-1
        x2=a['mass']+1
        range=[x1,x2]

        ar1=area1(spec_in[0],spec_in[1],range)
        ar2=area1(spec_dat[0],spec_dat[1],range)
        are1.append(ar1)
        are2.append(ar2)
    ma1=max(are1)
    ma2=max(are2)

    i=0
    for a in are1:
        are1[i]=are1[i]/ma1
        are2[i]=are2[i]/ma2
        if are1[i]==0:
            are1[i]=0.000001
        if are2[i]==0:
            are2[i]=0.000001
        i=i+1
        
    sum_ma1=sum(are1)
    sum_ma2=sum(are2)
    i=0
    for a in are1:
        are1[i]=are1[i]/sum_ma1
        are2[i]=are2[i]/sum_ma2
        i=i+1

    i=0
    le1=0
    le2=0
    for a in max1:
        #arc=arc+abs(are1[i]/ma1-are2[i]/ma2)/len(max1)
        la1=are2[i]*np.log(are2[i]/are1[i])
        arc=arc+la1
        i=i+1
    tabl=[max1,are1,are2,sum_ma1,sum_ma2]
    return  arc,tabl

def spect_table(match,results1,name):
    iso=[]
    ret=[]
    for a in results1:
        iso.append(a['isomer'])
        ret.append(a['Ret_peak'])
    
    i=0
    iso1=[]
    for i in range(len(iso)):
        if i==0:
            iso1.append(iso[i])
        if i!=0:
            rep=0
            for j in range(i):
                if iso[i]==iso[j]:
                    rep=1
                    break
            if rep==0:
                iso1.append(iso[i])
    print(iso)
    print(iso1)
    data1=[[iso1[i]] for i in range(len(iso1))]
    
    
    
    for b in results1:
        for i in range(len(iso1)):
            if b['isomer']==iso1[i]:
                inum=i
                mass1=['mass']
 
        for i in range(len(b['spectra_table'])):
            dat2=b['spectra_table'][i][2]
            dat1=b['spectra_table'][i][1]
            dat1.insert(0,b['spectra_table'][i][5])
            dat1.insert(0,b['Ret_peak'])
            data1[inum].append(dat1)
            data1[inum].append(dat2)
            #print(dat1)
            #xprint(len(a['spectra_table']))

    print(data1[0])
    
    title=['Isomer','m/z type','retention time','selected m/z','scaled intensity']

   
        
    

    
    '''
    for a in results1:
        nn=len(a['spectra_table'])
        print(nn)
        fig, axs = plt.subplots(nn, 1, figsize=(10, 4))
        for i in range(len(axs)):
            axs[i].set_title('Table'+'_'+a['spectra_table'][i][5])
            axs[i].axis('off')
            mass1=['mass']
            for b in a['spectra_table'][i][0]:
                mass1.append(b['mass'])
            for j in range(len(a['spectra_table'][i][1])):
                a['spectra_table'][i][1][j]=round(a['spectra_table'][i][1][j],4)
                a['spectra_table'][i][2][j]=round(a['spectra_table'][i][2][j],4)
            a['spectra_table'][i][1].insert(0,'input')
            a['spectra_table'][i][2].insert(0,'database')
            data=[mass1,a['spectra_table'][i][1],a['spectra_table'][i][2]]
            table = axs[i].table(cellText=data, loc='center')
            table.set_fontsize(28)
            #fig.set_size_inches(8, 4)
            #table.scale(1, 2)
        plt.show()
#        print(a['spectra_table'][1][5])
        
        break
    '''

def database_excel(nglycan):
    up_nglycan=nglycan.upper()
    match up_nglycan:
        case 'MAN1':
            inputexcel='./database/Man-1_database with chromatogram.xlsx'
        case 'MAN2':
            inputexcel='./database/Man-2_database with chromatogram.xlsx'
        case 'MAN3':
            inputexcel='./database/Man-3_database with chromatogram.xlsx'
        case 'MAN4':
            inputexcel='./database/Man-4_database with chromatogram.xlsx'
        case 'MAN5':
            inputexcel='./database/Man-5_database with chromatogram.xlsx'
        case 'MAN6':
            inputexcel='./database/Man-6_database with chromatogram.xlsx'
        case 'MAN7':
            inputexcel='./database/Man-7_database with chromatogram.xlsx'
        case 'MAN8':
            inputexcel='./database/Man-8_database with chromatogram.xlsx'
        case 'MAN9':
            inputexcel='./database/Man-9_database with chromatogram.xlsx'
        case 'MAN10':
            inputexcel='./database/Man-10_database with chromatogram.xlsx'

    return inputexcel


    
    

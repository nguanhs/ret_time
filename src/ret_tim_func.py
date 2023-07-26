import numpy as np
from numpy.linalg import inv
import pandas as pd
from scipy.signal import lfilter,savgol_filter
from scipy.interpolate import UnivariateSpline

def in_line(in1):
    ## in1 is the csv readed by panda
    for i in range(1000):
        a=in1.loc[i,0]
        b=in1.loc[i,1]
#        print(b)
        try:
            a=float(a)
            b=float(b)
            line=i
        except:
            line=2000
        if line==i:
            break
    return line

def data_line(in1,col):
    ## in1 is the csv readed by panda
    line=len(in1.loc[:,col])
    for i in range(len(in1.loc[:,col])):
        a=in1.loc[i,col]
        b=in1.loc[i,col+1]
#        print(b)
        try:
            a=float(a)
            b=float(b)
            line=i
        except:
            line=2000
        if line==i:
            break
    endline=len(in1.loc[:,col])
    for i in range(line,len(in1.loc[:,0]),1):
        a=in1.loc[i,col]
        b=in1.loc[i,col+1]
        if pd.isnull(in1.loc[i,col])==True:
            endline=i
            break
    endline=endline-1
    return line,endline


def read_input(in_file):
    df1 = pd.read_csv('input.csv',header=None)
    start=in_line(df1)
    #print(start,df1.loc[start,0],df1.loc[start-1,0])
    lstart=start
    lend=df1.index.stop-1
#    print(lstart,lend,df1.loc[lend,0])
    tim1=[] ## retention time
    int1=[] ## intensity
    for i in range(start,lend,1):
        tim1.append(float(df1.loc[i,0]))
        int1.append(float(df1.loc[i,1]))
    return tim1,int1

def read_database(in_file):
    df2 = pd.read_csv(in_file,header=None)
    nisomers=int(len(df2.loc[20])/2)

#    tim2=[[] for i in range(nisomers)] ## retention time
#    int2=[[] for i in range(nisomers)]
    tim2=[]
    int2=[]
    for i in range(nisomers):
        i1=2*i
        i2=2*i+1
        sline,eline=data_line(df2,i1)
        t1=[]
        in1=[]
        for j in range(sline,eline,1):
            t1.append(float(df2.loc[j,i1]))
            in1.append(float(df2.loc[j,i2]))
        tim2.append(t1)
        int2.append(in1)
    return nisomers,tim2,int2

def in_peaks(tim,int):
    gra1=grad1(tim,int)
    max_tint1=[] ## Dict. for time and intensity maximum.
    for i in range(len(gra1)-1):
        if gra1[i] > 0 and gra1[i+1] < 0:
            it1={'time':tim[i+1],'max_intens':int[i+1]}
            max_tint1.append(it1)
    max_tint1.sort(key=lambda x: x.get('max_intens'),reverse=True)
    return max_tint1

def in_peaks1(tim,int):
    gra1=grad1(tim,int)
    max_tint1=[] ## Dict. for time and intensity maximum.
    for i in range(len(gra1)-1):
        if gra1[i] > 0 and gra1[i+1] < 0:
            it1={'time':tim[i+1],'max_intens':int[i+1],'time_i':i+1}
            max_tint1.append(it1)
    max_tint1.sort(key=lambda x: x.get('max_intens'),reverse=True)
    return max_tint1

def in_peaks2(tim,int,theshold):
    gra1=grad1(tim,int)
    max_tint1=[] ## Dict. for time and intensity maximum.
    for i in range(len(gra1)-1):
        if gra1[i] > 0 and gra1[i+1] < 0 and int[i+1] > theshold:
            it1={'time':tim[i+1],'max_intens':int[i+1],'time_i':i+1}
            max_tint1.append(it1)
#    max_tint1.sort(key=lambda x: x.get('max_intens'),reverse=True)
    return max_tint1

def in_peaks3(tim,int,theshold):
    gra1=grad1(tim,int)
    max_tint1=[] ## Dict. for time and intensity maximum.
    for i in range(len(gra1)-1):
        if gra1[i] > 0 and gra1[i+1] < 0 and int[i+1] > theshold:
            it1={'mass':tim[i+1],'max_intens':int[i+1],'time_i':i+1}
            max_tint1.append(it1)
#    max_tint1.sort(key=lambda x: x.get('max_intens'),reverse=True)
    return max_tint1

def in_peaks4(tim,int,seleted_peaks):
    gra1=grad1(tim,int)
    max_tint1=[] ## Dict. for time and intensity maximum.
    for i in range(len(gra1)-1):
        dm=[]
        for a in seleted_peaks:
            dd=abs(tim[i+1]-a)
            dm.append(dd)   
        if gra1[i] > 0 and gra1[i+1] < 0 and min(dm) < 0.6:
            it1={'mass':tim[i+1],'max_intens':int[i+1],'time_i':i+1}
            max_tint1.append(it1)
#    max_tint1.sort(key=lambda x: x.get('max_intens'),reverse=True)
    return max_tint1

def database_peaks1(nisomers,tim2,int2):
    max_tint2=[]
    for i in range(nisomers):
        t1=tim2[i][:]
        in1=int2[i][:]
        col=in_peaks(t1,in1)
        max_tint2.append(col)
    return max_tint2


def database_peaks(nisomers,tim2,int2):
    max_tint2=[]
    for i in range(nisomers):
        t1=tim2[i][:]
        in1=int2[i][:]
        col=in_peaks(t1,in1)
        max_tint2.append(col)
    return max_tint2

def get_isomer(dat):
    iso_type=[]
    for a in dat:
        iso_type.append(a['isomer'])
    return iso_type

def get_ret1(dat):
    tim=[];inten=[];max_tint=[]
    for a in dat:
        tim.append(a['ret_time'])
        inten.append(a['ret_intensity'])
        max1=[]
        for b in a['ret_peaks']:
            dt=[]
            for j in range(len(a['ret_time'])):
                dt.append(abs(a['ret_time'][j]-b))
            mdt=min(dt)
            k=dt.index(mdt)
            it1={'time':a['ret_time'][k],'max_intents':a['ret_intensity'][k]}
            max1.append(it1)
        max_tint.append(max1)
    return tim,inten,max_tint

def get_ret(dat):
    tim=[];inten=[]
    for a in dat:
        tim.append(a['ret_time'])
        inten.append(a['ret_intensity'])

    return tim,inten


def get_msspectra(dat,isomer,pref):
    for a in dat:
        if a['isomer']==isomer:
            ma=a[pref+'_mass']
            inten=a[pref+'_intensity']
            break
    return ma,inten

def match_data_input1(allraw,datab):
    tim_in,int_in=ms2ret(allraw)
    max_tint_in=in_peaks2(tim_in,int_in,0.1)
    numb_of_peaks=len(max_tint_in)

    isomers=get_isomer(datab)
    nisomers=len(isomers)
#    tim_d,int_d=get_ret(datab)
    tim_d1,int_d1,max_tint_d1=get_ret1(datab)
#    max_tint_d=database_peaks(nisomers,tim_d,int_d)
#    print(max_tint_d1)
#    print(max_tint_d)
    match,match1=match_data_input(max_tint_in,max_tint_d1,nisomers,numb_of_peaks)
    mat1=match_data_input_1(max_tint_in,max_tint_d1,nisomers,numb_of_peaks)
    return match, mat1 

def match_data_input2(allraw,mass,datab):
    tim_in,int_in=ms2ret_mass(allraw,mass)
    max_tint_in=in_peaks2(tim_in,int_in,0.1)
    numb_of_peaks=len(max_tint_in)

    isomers=get_isomer(datab)
    nisomers=len(isomers)
#    tim_d,int_d=get_ret(datab)
    tim_d1,int_d1,max_tint_d1=get_ret1(datab)
#    max_tint_d=database_peaks(nisomers,tim_d,int_d)
#    print(max_tint_d1)
#    print(max_tint_d)
    match,match1=match_data_input(max_tint_in,max_tint_d1,nisomers,numb_of_peaks)
    mat1=match_data_input_1(max_tint_in,max_tint_d1,nisomers,numb_of_peaks)
    return match, mat1



def match_data_input(max_tint1,max_tint2,nisomers,numb_of_peaks):
    match=[];match1=[]
    for i in range(nisomers):
        dt1=[]
        dt2=[]
        t1=max_tint2[i][0].get('time')
        t2=max_tint2[i][1].get('time')
        for j in range(numb_of_peaks):
            t3=max_tint1[j].get('time')
            dt1.append(abs(t3-t1))
            dt2.append(abs(t3-t2))
        mdt1=min(dt1)
        idt1=dt1.index(mdt1)
        mdt2=min(dt2)
        idt2=dt2.index(mdt2)
        if mdt1 < 1 and mdt2 < 1:
#            m=[i,mdt1,mdt2,max_tint1[idt1].get('time'),max_tint1[idt2].get('time'),t1,t2]
            m=[i,mdt1,mdt2,idt1,idt2,t1,t2]
            m1={'iso_num':i,'peak_num1':idt1,'peak_num2':idt2,'time_dat1':t1,'time_dat2':t2}
            match.append(m)
            match1.append(m1)
    return match,match1

def match_data_input_1(max_tint1,max_tint2,nisomers,numb_of_peaks):
    match=[];match1=[]
#    for i in range(nisomers):
    for i in range(numb_of_peaks):
#        dt1=[]
#        dt2=[]
        t3=max_tint1[i].get('time')
        for j in range(nisomers):
            t1=max_tint2[j][0].get('time')
            t2=max_tint2[j][1].get('time')
            dt1=abs(t3-t1)
            if dt1 < 1:
                #if i==1:
                #print(j,t1,t3)
                #    print('ok')
                for k in range(numb_of_peaks):
                    k1=k
                    t4=max_tint1[k1].get('time')
                    dt2=abs(t4-t2)
                    #if i==1:
                    #print(k1,dt2,t4,'ok')
                    if dt2 < 1:
                    # found match.
                        m1={'iso_num':j,'peak_num1':i,'peak_num2':k1,'time_dat1':t1,'time_dat2':t2}
                        m2=[j,dt1,dt2,i,k1,t1,t2] 
                        match.append(m2)
                '''
                if i+2==numb_of_peaks:
                    k1=i+1
                    t3=max_tint1[k1].get('time')
                    dt2=abs(t3-t2)
                    if dt2 < 1:
                    # found match.
                        m1={'iso_num':j,'peak_num1':i,'peak_num2':k1,'time_dat1':t1,'time_dat2':t2}
                        m2=[j,dt1,dt2,i,k1,t1,t2]
                        match.append(m2)
                '''
    return match


def grad1(x,y):
    dy=[]
    for i in range(len(x)-1):
        dy1=(y[i+1]-y[i])/(x[i+1]-x[i])
        dy.append(dy1)
    return dy

def contri(max_tint1,max_tint2,match):
##Calculate int(Gg)
    gg1=[]
    for a in match:
        j=a[3]
        g1=max_tint1[j].get('max_intens')
        g1a=max_tint2[a[0]][0].get('max_intens')
        j=a[4]
        g2=max_tint1[j].get('max_intens')
        g2a=max_tint2[a[0]][1].get('max_intens')
#        gg1.append((g1*g1a+g2*g2a)/100)
        gg1.append((g1*g1a+g2*g2a))


    ##Calculate int(gg)
    gg2=[]
    for a in match:
        aa=[]
        for b in match:
            if a[0]==b[0]:
                g1=max_tint2[a[0]][0].get('max_intens')
                g2=max_tint2[a[0]][1].get('max_intens')
                #g=(g1*g1+g2*g2)/10000
                g=(g1*g1+g2*g2)
            if a[0] != b[0]:
                g=0
                for i in range(2):
                    t1=max_tint2[a[0]][i].get('time')
                    g1=max_tint2[a[0]][i].get('max_intens')
                    for j in range(2):
                        t2=max_tint2[b[0]][j].get('time')
                        g2=max_tint2[b[0]][j].get('max_intens')
                        if abs(t1-t2) < 0.5:
                            #g=g+g1*g2/10000
                            g=g+g1*g2
            aa.append(g)
        gg2.append(aa)


    gg2a=np.array(gg2)
    gg1a=np.array(gg1)
    igg2a=inv(gg2a)
    aa=np.matmul(igg2a, gg1a)
#    aa=aa/100
    aa=aa
    return aa

def contri1(allraw,datab,match):
    tim_in,int_in=ms2ret(allraw)
    max_tint1=in_peaks2(tim_in,int_in,0.1)
#    numb_of_peaks=len(max_tint_in)

    isomers=get_isomer(datab)
    nisomers=len(isomers)
    tim_d,int_d=get_ret(datab)
    max_tint2=database_peaks(nisomers,tim_d,int_d)


##Calculate int(Gg)
    gg1=[]
    for a in match:
        j=a[3]
        g1=max_tint1[j].get('max_intens')
        g1a=max_tint2[a[0]][0].get('max_intens')
        j=a[4]
        g2=max_tint1[j].get('max_intens')
        g2a=max_tint2[a[0]][1].get('max_intens')
#        gg1.append((g1*g1a+g2*g2a)/100)
        gg1.append((g1*g1a+g2*g2a))


    ##Calculate int(gg)
    gg2=[]
    for a in match:
        aa=[]
        for b in match:
            if a[0]==b[0]:
                g1=max_tint2[a[0]][0].get('max_intens')
                g2=max_tint2[a[0]][1].get('max_intens')
                #g=(g1*g1+g2*g2)/10000
                g=(g1*g1+g2*g2)
            if a[0] != b[0]:
                g=0
                for i in range(2):
                    t1=max_tint2[a[0]][i].get('time')
                    g1=max_tint2[a[0]][i].get('max_intens')
                    for j in range(2):
                        t2=max_tint2[b[0]][j].get('time')
                        g2=max_tint2[b[0]][j].get('max_intens')
                        if abs(t1-t2) < 0.5:
                            #g=g+g1*g2/10000
                            g=g+g1*g2
            aa.append(g)
        gg2.append(aa)


    gg2a=np.array(gg2)
    gg1a=np.array(gg1)
    igg2a=inv(gg2a)
    aa=np.matmul(igg2a, gg1a)
#    aa=aa/100
    aa=aa
    return aa



#def compare(allraw,datab,a[0],a[3],a[4]):
#    a=1

## get ms2 total ion count
def ms2ret(dat):
    time=[];intensity=[]
    for a in dat:
#         print(a['ms'])
        if a['ms']=='2':
            time.append(float(a['scan start time']))
            intensity.append(float(a['total ion current']))
    max_tint1=in_peaks(time,intensity)
    for i in range(len(intensity)):
        intensity[i]=intensity[i]/max_tint1[0]['max_intens']

    return time,intensity

def ms2ret_mass(dat,mass):
    time=[];intensity=[]
    ptim=0
    for a in dat:
#         print(a['ms'])
        #print(a['ms'],a['Tandem_mass'][0])
        
        if a['ms']=='2' and int(a['Tandem_mass'][0])==mass:
            if ptim==float(a['scan start time']):
                #print(ptim)
                break
            time.append(float(a['scan start time']))
            intensity.append(float(a['total ion current']))
            ptim=float(a['scan start time'])
            #print(a['scan start time'],a['total ion current'])
            #print(a)
    max_tint1=in_peaks(time,intensity)
    for i in range(len(intensity)):
        intensity[i]=intensity[i]/max_tint1[0]['max_intens']

    return time,intensity

##get ms_level and selected  mass
def ms_lev(dat):
    ms_level=[]
    selected_mass=[]
    i=0
    for a in dat:
        l=a['ms']
        mass=a['selected mass']
        if i > 1 and l=='1':
            break
        if i != 0:
            ms_level.append(l)
            selected_mass.append(int(float(mass)))
        i=i+1
    return ms_level,selected_mass
  
def ms_lev_mass(dat,mass_ms2):
    ms_level=[]
    selected_mass=[]
    i=0
    for a in dat:
        l=a['ms']
        mass=a['selected mass']
        tandem=a['Tandem_mass']
        
        if i > 1 and l=='1':
            break
        if i != 0 and int(tandem[0])==mass_ms2:
            ms_level.append(l)
            selected_mass.append(int(float(mass)))
        i=i+1
    return ms_level,selected_mass

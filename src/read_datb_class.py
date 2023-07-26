import numpy as np
import pandas as pd


class database():
    '''
    Read database from excel
    '''

    def __init__(self, inputexcel, nglycan):
        self.inputexcel=inputexcel       #/pathtofile/excelfile
        self.nglycan=nglycan             #nglycan type: Man7,Man8 ....

    def data(self):
        df1=pd.read_excel(self.inputexcel,sheet_name=None)
        aa=list(df1.keys())
        a1=list(df1[aa[0]].axes[1])
        nl=int(len(a1)/2)
        ngly=self.nglycan
        match ngly.upper():
            case "MAN1":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['Man1']=[20.1,26.7]

            case "MAN2":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['2F1']=[28.6,36.7]
               pt['2E1']=[27.3,34.6]

            case "MAN3":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['3D1']=[36.2,44.5]
               pt['3F1']=[17.8,23.5]
               pt['3E2']=[23.4,30.0]
               pt['3E1']=[30.6,38.0]

            case "MAN4":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['4D2']=[31.4,38.7]
               pt['4D1']=[38.9,47.4]
               pt['4D3']=[24.2,30.6]
               pt['4E1']=[29.2,36.0]
               pt['4E2']=[28.9,36.0]
               pt['4E3']=[25.7,32.8]
               pt['4F1']=[18.6,24.6]

            case "MAN5":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['5F2']=[30.8,38.9]
               pt['5D2']=[22.9,28.9]
               pt['5D1']=[28.5,34.9]
               pt['5E4']=[24.9,30.9]
               pt['5E2']=[36.7,45.2]
               pt['5E3']=[33.7,41.6]
               pt['5E1']=[37.0,44.4]
               pt['5F1']=[24.7,30.2]

            case "MAN6":
               pt={}
               pt['6H1']=[28.7,35.2]
               pt['6F1']=[29.8,36.0]
               pt['6F2']=[39.2,47.5]
               pt['6D3']=[28.5,34.8]
               pt['6D1']=[27.5,33.4]
               pt['6D2']=[23.6,30.3]
               pt['6G1']=[25.4,31.1]
               pt['6E2']=[29.0,35.1]
               pt['6E1']=[23.2,28.8]

            case "MAN7":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['7F1']=[28.0,34.2]
               pt['7E1']=[27.1,31.7]
               pt['7D3']=[28.6,34.4]
               pt['7G1']=[31.2,37.8]
               pt['7E2']=[29.5,36.1]
               pt['7D2']=[23.9,29.9]
               pt['7D1']=[28.7,34.7]

            case "MAN8":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['8D3']=[32.8,39.2]
               pt['8E2']=[29.1,35.5]
               pt['8G1']=[25.3,30.8]
               pt['8E1']=[26.7,31.4]

            case "MAN9":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['9E1']=[25.7,31.5]
               pt['9D1']=[29.0,34.0]
               pt['9D2']=[34.2,41.2]

            case "MAN10":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['10D1']=[29.1,35.1]



        dat=[]
        for i in range(len(a1)):
            j=i%2
            if j==1:
                a2=a1[i].split()[0]
                s1={'isomer':a2,'ret_peaks':pt[a2]}
                dat.append(s1)

        for i in range(len(aa)):
            a1=list(df1[aa[i]].axes[1])
            a2=df1[aa[i]].to_numpy()
            a2=np.transpose(a2)
            j2=0
            for j in range(len(a1)):
                j1=j%2
                j2=int(j/2)
                
                a3=[]
                for l1 in a2[j]:
                    try:
                        a4=float(l1)
                    except:
                        l1=np.nan
                    if np.isnan(l1):
                        break
                    else:
                        a3.append(float(l1))
                maxa=max(a3)
                if j1==0 and i==0:
                    dat[j2]["ret_time"]=a3
                if j1==1 and i==0:
                    for k1 in range(len(a3)):
                        a3[k1]=a3[k1]/maxa
                    dat[j2]["ret_intensity"]=a3
           
                if j1==0 and i > 0:
                    b=a1[j].split()[0]
                    for k in range(len(dat)):
                        if dat[k]['isomer']==b:
                            j3=k 
                    dat[j3][aa[i]+"_mass"]=a3
                if j1==1 and i > 0:
                    b=a1[j].split()[0]
                    for k1 in range(len(a3)):
                        a3[k1]=a3[k1]/maxa
                    for k in range(len(dat)):
                        if dat[k]['isomer']==b:
                            j3=k
                    dat[j3][aa[i]+"_intensity"]=a3
        return dat   
        
    def get_isomer(self,dat):
        iso_type=[]
        for a in dat:
            iso_type.append(a['isomer'])
        return iso_type

    def get_ret(self,dat):
        tim=[];inten=[]
        for a in dat:
            tim.append(a['ret_time'])
            inten.append(a['ret_intensity']) 
        return tim,inten

    def get_retpeak(self,dat):
        ngly=self.nglycan
        tim=[];inten=[];ptim=[]
        match ngly.upper():

            case "MAN1":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['Man1']=[20.1,26.7]

            case "MAN2":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['2F1']=[28.6,36.7]
               pt['2E1']=[27.3,34.6]

            case "MAN3":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['3D1']=[36.2,44.5]
               pt['3F1']=[17.8,23.5]
               pt['3E2']=[23.4,30.0]
               pt['3E1']=[30.6,38.0]

            case "MAN4":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['4D2']=[31.4,38.7]
               pt['4D1']=[38.9,47.4]
               pt['4D3']=[24.2,30.6]
               pt['4E1']=[29.2,36.0]
               pt['4E2']=[28.9,36.0]
               pt['4E3']=[25.7,32.8]
               pt['4F1']=[18.6,24.6]

            case "MAN5":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['5F2']=[30.8,38.9]
               pt['5D2']=[22.9,28.9]
               pt['5D1']=[28.5,34.9]
               pt['5E4']=[24.9,30.9]
               pt['5E2']=[36.7,45.2]
               pt['5E3']=[33.7,41.6]
               pt['5E1']=[37.0,44.4]
               pt['5F1']=[24.7,30.2]

            case "MAN6":
               pt={}
               pt['6H1']=[28.7,35.2]
               pt['6F1']=[29.8,36.0]
               pt['6F2']=[39.2,47.5]
               pt['6D3']=[28.5,34.8]
               pt['6D1']=[27.5,33.4]
               pt['6D2']=[23.6,30.3]
               pt['6G1']=[25.4,31.1]
               pt['6E2']=[29.0,35.1]
               pt['6E1']=[23.2,28.8]

            case "MAN7":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['7F1']=[28.0,34.2]
               pt['7E1']=[27.1,31.7]
               pt['7D3']=[28.6,34.4]
               pt['7G1']=[31.2,37.8]
               pt['7E2']=[29.5,36.1]
               pt['7D2']=[23.9,29.9]
               pt['7D1']=[28.7,34.7]

            case "MAN8":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['8D3']=[32.8,39.2]
               pt['8E2']=[29.1,35.5]
               pt['8G1']=[25.3,30.8]
               pt['8E1']=[26.7,31.4]

            case "MAN9":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['9E1']=[25.7,31.5]
               pt['9D1']=[29.0,34.0]
               pt['9D2']=[34.2,41.2]

            case "MAN10":
               #pt=[{} for _ in range(nl)]
               pt={}
               pt['10D1']=[29.1,35.1]



        for a in dat:
            tim.append(a['ret_time'])
            inten.append(a['ret_intensity'])           
            ptim.append(pt[a['isomer']])


        return tim,inten,ptim





    def get_msspectra(self,dat,isomer,pref):
#        ma=[];inten=[]
        for a in dat:
            if a['isomer']==isomer:
                ma=a[pref+'_mass']
                inten=a[pref+'_intensity']
                break
        return ma,inten

def dataread(infile):
    a=1

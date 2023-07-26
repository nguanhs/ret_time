import numpy as np
import xml.etree.ElementTree as et
import base64
import zlib
import os
import src.ret_tim_func as ret


class rawdata():
    '''
    Convert *raw to mzml and extract the data
    ''' 

    def __init__(self, input_raw, output_mzml,thermoscript):
        self.input_raw=input_raw        #/pathtofile/*raw
        self.output_mzml=output_mzml    #/pathtofile/*mzml
        self.thermoscript=thermoscript  #shell script to convert raw to mzml
    
    def rawtomzml(self):
        '''
        convert .raw file to .mzml file
        '''
        cwd=os.getcwd()
        com="sh "+cwd+"/src/"+self.thermoscript+" "+self.input_raw 
        os.system(com)

    def readxml(self):
        ##check mxml existance     
        '''
        read data into list of dictionary of each scan
        ['scan', 'ms', 'total ion current', 'scan start time', 'experiment', 'selected mass','Tandem_mass', 'masses', 'intensity']
        '''
        checkexist=os.path.exists(self.output_mzml)
        if checkexist==False:
            print(self.output_mzml+' is not found!')
            exit()
        txml = et.parse(self.output_mzml)
        root=txml.getroot()

        i=0
        dat=[]
        a2=0
        #for neighbor in root[6][0].iter('{http://psi.hupo.org/ms/mzml}cvParam'):
        #    print(neighbor.attrib,type(neighbor.attrib))
        for neighbor in root.iter():
            a1=neighbor.attrib
            if neighbor.tag=='{http://psi.hupo.org/ms/mzml}cvParam':
                if a1['name']=="ms level":
                    i=i+1
                    is1={"scan":i,"ms":a1['value']}    
                if a1['name']=="total ion current":
                    is1["total ion current"]=a1['value']
                if a1['name']=="scan start time":
                    is1["scan start time"]=a1['value']
                if a1['name']=="filter string":
                    is1["experiment"]=a1['value']
                    c1=a1['value'].split()
                    if c1[5]=='ms':
                        is1["selected mass"]=c1[6]
                        is1["Tandem_mass"]='no'
                    elif c1[5]=='ms2':
                        is1["selected mass"]=c1[6].split('@')[0]
                        is1["Tandem_mass"]=[float(c1[6].split('@')[0])]
                    elif c1[5]=='ms3':
                        is1["selected mass"]=c1[7].split('@')[0]
                        is1["Tandem_mass"]=[float(c1[6].split('@')[0]),float(c1[7].split('@')[0])]
                    elif c1[5]=='ms4':
                        is1["selected mass"]=c1[8].split('@')[0]
                        is1["Tandem_mass"]=[float(c1[6].split('@')[0]),float(c1[7].split('@')[0]),float(c1[8].split('@')[0])]
 
                if a1['name']=="m/z array":
                    a2=1
                if a1['name']=="intensity array":
                    a2=2
            if neighbor.tag=='{http://psi.hupo.org/ms/mzml}binary' and a2==1:
                is1["masses"]=neighbor.text
            if neighbor.tag=='{http://psi.hupo.org/ms/mzml}binary' and a2==2:
                is1["intensity"]=neighbor.text
                a2=0
                dat.append(is1)
        return dat       
    '''
    def ms2ret(self,dat):
        time=[];intensity=[]
        for a in dat:
#            print(a['ms'])
            if a['ms']=='2':
                time.append(float(a['scan start time']))
                intensity.append(float(a['total ion current']))
        return time,intensity
    '''
    def ms2ret(self,dat):
        time=[];intensity=[]
        for a in dat:
#            print(a['ms'])
            if a['ms']=='2':
                time.append(float(a['scan start time']))
                intensity.append(float(a['total ion current']))
        max_tint1=ret.in_peaks(time,intensity)
        for i in range(len(intensity)):
            intensity[i]=intensity[i]/max_tint1[0]['max_intens']

        return time,intensity

    def ms_lev(self,dat):
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


    def spectra_centroid(self,dat,t_ini,t_end,ms_sel):
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

        ind1=[]
        dm=0.0000001
        dn=7
        for arr in mi:
            i2=[]
            for m1 in arr:
                b=round(m1,dn)
                i1=(b-400.0)/dm
                i2.append(int(i1))
            ind1.append(i2)
        sum1=0
        a1=ind1[0][0]
        ind2=[a1]
        for i1 in range(len(ind1)):
            for a in ind1[i1]:
                for b in ind2:
                    if a != b:
                        ind2.append(a)
                        break
        
        ind2.sort() 
        intent_com=[]
        mass_com=[]
        for i1 in ind2:
            in1=0
            j=0
            for arr in ind1:
                for i in range(len(arr)):
                    if arr[i]==i1:
                        in1=in1+inte[j][i]
                        break
                j=j+1
            intent_com.append(in1)
            mass1=400.0+dm*i1
            mass_com.append(mass1)
        return mass_com,intent_com,mi,inte


    def spectra_profile(self,dat,t_ini,t_end,ms_sel):
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
#        return mass_com,intent_com,mi,inte
        return mass_com,intent_com,maxi

    def t_inter(self,retention_peaks,ret_int1,ret_time_num):
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


#import lxml
#from lxml import etree
import xml.etree.ElementTree as et
import base64
import zlib

txml = et.parse('./blackbeans_Man7a.mzML')
root=txml.getroot()
#print(root.attrib)

#for child in root:
#    print(type(child.tag))


#for child in root.iter('{http://psi.hupo.org/ms/mzml}run'):
#    print(child.tag, child.attrib)



#print(type(root))
#for child in root:
#for neighbor in root[6].iter('{http://psi.hupo.org/ms/mzml}binary'):   
#    print(neighbor.attrib)
'''
for neighbor in root[0].iter():
    print(neighbor.tag)

for neighbor in root[1].iter():
    print(neighbor.tag)
'''
for neighbor in root.iter():
    print(neighbor.tag)

print('---')
#exit()
#for n in root:
#    print(n.tag)

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
        if a1['name']=="total ion current":
            is1["total ion current"]=a1['value']
#    if a1['name']=="isolation window target m/z":
#        is1["selected mass"]=a1['value']
#        a2=1
        if a1['name']=="filter string":
            is1["experiment"]=a1['value']
        if a1['name']=="m/z array":
            a2=1
        if a1['name']=="intensity array":
            a2=2
#            dat.append(is1)
    if neighbor.tag=='{http://psi.hupo.org/ms/mzml}binary' and a2==1:
        is1["masses"]=neighbor.text
        print(neighbor.text)
        exit()
    if neighbor.tag=='{http://psi.hupo.org/ms/mzml}binary' and a2==2:
        is1["intensity"]=neighbor.text
        a2=0
        dat.append(is1)
       

#        if a1['name']=="m/z array":
#            a2=1
#            dat.append(is1)

#for neighbor in root[6][0].iter('{http://psi.hupo.org/ms/mzml}'):

print(dat)
#        exit()
    


#for neighbor in root.findall('./run/spectrumList/spectrum/[cvParam]'):
#    print(neighbor.attrib)
#

a1='{http://psi.hupo.org/ms/mzml}'
a2='spectrum'

#for a in root.iter(a1+a2):
#    print(a.tag,a.attrib)

#for a in root[6].iter('{http://psi.hupo.org/ms/mzml}binary'): 
#    print(a.text)
 

import numpy as np
import base64
import zlib
import os
import src.ret_tim_func as ret
import src.use_func as usf

class lodes2():
    '''
    Nglycan structure determination by LODES
    '''

    def __init__(self, nglycan):
#        self.rawdata=rawdata                    #Raw spectra in python list of dictionaries 
#        self.database=datab                  #Database in python dictionary 
#        self.ret_int1=ret_int1                 #Retention time and intensity in python list of lists
#        self.retention_peaks=retention_peaks    #retention_peaks in python dictionary
        self.nglycan=nglycan                    #kind of nglycan
#        self.raw=self.raw()

    def mass_select(self):
        ngly=self.nglycan
        ngly.upper()
        met=[]
        match ngly.upper():
            case "MAN7":
                targeted_mass=[1581,1175,1157]
                ms_level=[2,3,3]
        return targeted_mass,ms_level

    def isomer_picture(self):
        ngly=self.nglycan
        ngly.upper()
        met=[]
        d1='./image_man/'
        match ngly.upper():

            case "MAN1":
                isomer=['Man1']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

            case "MAN2":
                isomer=['2E1','2F1']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

            case "MAN3":
                isomer=['3D1','3F1','3E2','3E1']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

            case "MAN4":
                isomer=['4D3','4D1','4D2','4F1','4E3','4E2','4E1']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

            case "MAN5":
                isomer=['5F1', '5E1', '5E3', '5E2', '5D2', '5E4', '5D1','5F2']
                # pref='prefix'
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)
            
            case "MAN6":
                isomer=['6H1', '6E1', '6E2', '6G1', '6D2', '6D1', '6D3','6F2','6F1']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

            case "MAN7":
                isomer=['7E1', '7D3', '7G1', '7E2', '7D2', '7D1', '7F1']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

            case "MAN8":
                isomer=['8E1', '8G1', '8E2', '8D3']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

            case "MAN9":
                isomer=['9E1', '9D1', '9D2']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

            case "MAN10":
                isomer=['10D1']
                pict=[]
                for a in isomer:
                    b=d1+a+'.png'
                    pict.append(b)

        return pict

    def method(self):
        ngly=self.nglycan
        ngly.upper()
        met=[]
        de1=[18,101,162,203,221,324,424,486]
        de1_com={18:[('H',2),('O',1)],101:[('C',4),('H',7),('O',2),('N',1)],203:[('C',8),('H',13),('O',5),('N',1)],
            162:[('C',6),('H',10),('O',5)],221:[('C',8),('H',15),('O',6),('N',1)],324:[('C',12),('H',20),('O',10)],
            486:[('C',18),('H',30),('O',15)],424:[('C',16),('H',28),('O',11),('N',2)]}
        frag=[203,226,244]
        frag_com={203:[('C',6),('H',12),('O',6),('Na',1)],226:[('C',8),('H',13),('O',5),('N',1),('Na',1)],
             244:[('C',8),('H',15),('O',6),('N',1),('Na',1)]}
        match ngly.upper():
            case "MAN1":
                isomer=['Man1']
                pref=['Man-1_MS2-609']
                m_sel=[609]
                m_sel_iso=[609.211907]
                ms_level=[2]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[609],'m_sel_iso':[609.211907]})

                #minmass=[1100,500,500,200]
                minmass=[100,100,100,100]
                ms2a=[591,508,406,447,388]
                ms2a_iso=[591.201, 508.164, 406.133, 447.159, 388.122]
                ms2b=[203,226,244]
                ms2b_iso=[203.053, 226.069, 244.08]
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso}

            case "MAN2":
                isomer=['2E1','2F1']
                pref=['Man-2_MS2-771','Man-2_MS3-365']
                m_sel=[771,365]
                ms_level=[2,3]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[771,365,347],'m_sel_iso':[771.265, 365.106, 347.095]})
                met.append({1:isomer[1],'m_sel':[771,365,347],'m_sel_iso':[771.265, 365.106, 347.095]})
                m_sel_iso=[771.265,365.106, 347.095]                

                #minmass=[1100,500,500,200]
                minmass=[100,100,100,100]
                ms2a=[753,670,550,447,609,568]
                ms2a_iso=[753.254, 670.217, 550.175, 447.159, 609.212, 568.185]
                ms2b=[347,226,244,365]
                ms2b_iso=[347.095, 226.069, 244.08, 365.106]
                ms2c=[]
                ms3_365=[347,305,275,245,203,185] ##Cion
                ms3_365_iso=[347.095, 305.085, 275.074, 245.064, 203.053, 185.043]
                ms3_347=[329,287,301,275,227,203,185] ##B ion
                ms3_347_iso=[329.085, 287.074, 301.054, 275.074, 227.053, 203.053, 185.043]
#                peaks_sel={'ms2a':ms2a}
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms3_365':ms3_365,
                           'ms3_347':ms3_347,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms3_365_iso':ms3_365_iso,
                           'ms3_347_iso':ms3_347_iso}


            case "MAN3":
                isomer=['3D1','3F1','3E2','3E1']
                pref=['Man-3_MS2-933','Man-3_MS3-509','Man-3_MS3-527','Man-3_MS4-365']
                m_sel=[933,509,527,365]
                m_sel_iso=[933.318, 509.148, 527.159, 365.106]
                ms_level=[2,3,3,4]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[933,527,509],'m_sel_iso':[933.318, 527.159, 509.148]})
                met.append({1:isomer[1],'m_sel':[933,527,509],'m_sel_iso':[933.318, 527.159, 509.148]})
                met.append({2:isomer[2],'m_sel':[933,527,509,365],'m_sel_iso':[933.318, 527.159, 509.148, 365.106]})
                met.append({3:isomer[3],'m_sel':[933,527,509,365],'m_sel_iso':[933.318, 527.159, 509.148, 365.106]})

                #minmass=[1100,500,500,200]
                minmass=[100,100,100,100]
                ms2a=[915,832,771,730,712]
                ms2a_iso=[915.307, 832.27, 771.265, 730.238, 712.228]
                ms2b=[509,447,609]
                ms2b_iso=[509.148, 447.159, 609.212]
                ms2c=[]
                ms3_527=[509,467,437,407,365,347,275]
                ms3_527_iso=[509.148, 467.138, 437.127, 407.117, 365.106, 347.095, 275.074]
                ms3_509=[449,437,365,347,329]
                ms3_509_iso=[449.127, 437.127, 365.106, 347.095, 329.085]
                ms4_365=[347,305,275,245,203]
                ms4_365_iso=[347.095, 305.085, 275.074, 245.064, 203.053]
#                peaks_sel={'ms2a':ms2a}
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms3_527':ms3_527,
                           'ms3_509':ms3_509,'ms4_365':ms4_365,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms3_527_iso':ms3_527_iso,
                           'ms3_509_iso':ms3_509_iso,'ms4_365_iso':ms4_365_iso}


            case "MAN4":
                isomer=['4D3','4D1','4D2','4F1','4E3','4E2','4E1']
                pref=['Man-4_MS2-1095','Man-4_MS3-689','Man-4_MS3-671','Man-4_MS4-365','Man-4_MS4-527']
                m_sel=[1095,689,671,365,527]
                m_sel_iso=[1095.37, 689.212, 671.201, 365.106, 527.159]
                ms_level=[2,3,3,4,4]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[1095,689,671,365],'m_sel_iso':[1095.37, 689.212, 671.201, 365.106]})
                met.append({1:isomer[1],'m_sel':[1095,689,671,365],'m_sel_iso':[1095.37, 689.212, 671.201, 365.106]})
                met.append({2:isomer[2],'m_sel':[1095,689,671,365],'m_sel_iso':[1095.37, 689.212, 671.201, 365.106]})
                met.append({3:isomer[3],'m_sel':[1095,689,671,527],'m_sel_iso':[1095.37, 689.212, 671.201, 527.159]})
                met.append({4:isomer[4],'m_sel':[1095,689,671,527],'m_sel_iso':[1095.37, 689.212, 671.201, 527.159]})
                met.append({5:isomer[5],'m_sel':[1095,689,671,527],'m_sel_iso':[1095.37, 689.212, 671.201, 527.159]})
                met.append({6:isomer[6],'m_sel':[1095,689,671,527],'m_sel_iso':[1095.37, 689.212, 671.201, 527.159]})
                
                #minmass=[1100,500,500,200]
                minmass=[100,100,100,100]
                ms2a=[1077,994,933,892,874,671]
                ms2a_iso=[1077.36, 994.323, 933.318, 892.291, 874.28, 671.201]
                ms2b=[609,550,509,447,771]
                ms2b_iso=[609.212, 550.175, 509.148, 447.159, 771.265]
                ms2c=[]
                ms3_689=[671,629,599,569,527,509,437]
                ms3_689_iso=[671.201, 629.191, 599.18, 569.169, 527.159, 509.148, 437.127]
                ms3_671=[599,527,509,347]
                ms3_671_iso=[599.18, 527.159, 509.148, 347.095]
                ms4_365=[347,305,275,245,203]
                ms4_365_iso=[347.095, 305.085, 275.074, 245.064, 203.053]
                ms4_527=[509,467,437,407,365,347,275]
                ms4_527_iso=[509.148, 467.138, 437.127, 407.117, 365.106, 347.095, 275.074]
#                peaks_sel={'ms2a':ms2a}
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms3_689':ms3_689,
                           'ms3_671':ms3_671,'ms4_365':ms4_365,'ms4_527':ms4_527,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms3_689_iso':ms3_689_iso,
                           'ms3_671_iso':ms3_671_iso,'ms4_365_iso':ms4_365_iso,'ms4_527_iso':ms4_527_iso}


            case "MAN5":
                isomer=['5F1', '5E1', '5E3', '5E2', '5D2', '5E4', '5D1','5F2']
                pref=['Man-5_MS2-1257','Man-5_MS3-851','Man-5_MS3-833','Man-5_MS4-527']
                m_sel=[1257,851,833,527]
                m_sel_iso=[1257.423, 851.264, 833.254, 527.159]
                ms_level=[2,3,3,4]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[1257,851,833],'m_sel_iso':[1257.423, 851.264, 833.254]})
                met.append({1:isomer[1],'m_sel':[1257,851,833,527],'m_sel_iso':[1257.423, 851.264, 833.254,527.159]})
                met.append({2:isomer[2],'m_sel':[1257,851,833],'m_sel_iso':[1257.423, 851.264, 833.254]})
                met.append({3:isomer[3],'m_sel':[1257,851,833,527],'m_sel_iso':[1257.423, 851.264, 833.254,527.159]})
                met.append({4:isomer[4],'m_sel':[1257,851,833],'m_sel_iso':[1257.423, 851.264, 833.254]})
                met.append({5:isomer[5],'m_sel':[1257,851,833],'m_sel_iso':[1257.423, 851.264, 833.254]})
                met.append({6:isomer[6],'m_sel':[1257,851,833],'m_sel_iso':[1257.423, 851.264, 833.254]})
                met.append({7:isomer[7],'m_sel':[1257,851,833],'m_sel_iso':[1257.423, 851.264, 833.254]})
                #minmass=[1100,500,500,200]
                minmass=[100,100,100,100]
                ms2a=[1239,1156,1095,1054,1036,833]
                ms2a_iso=[1239.413,1156.376,1095.37,1054.344,1036.333,833.254]
                ms2b=[933,771,712,671,609]
                ms2b_iso=[933.318, 771.265, 712.228, 671.201, 609.212]
                ms2c=[]
                ms3_851=[833,791,761,731,689,599,527,509,437]
                ms3_851_iso=[833.254, 791.243, 761.233, 731.222, 689.212, 599.18, 527.159, 509.148, 437.127]
                ms3_833=[761,689,671,509,347]
                ms3_833_iso=[761.233, 689.212, 671.201, 509.148, 347.095]
                ms4_527=[509,467,437,407,365,347,275]
                ms4_527_iso=[509.148, 467.138, 437.127, 407.117, 365.106, 347.095, 275.074]
#                peaks_sel={'ms2a':ms2a}
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms3_851':ms3_851,
                           'ms3_833':ms3_833,'ms4_527':ms4_527,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms3_851_iso':ms3_851_iso,
                           'ms3_833_iso':ms3_833_iso,'ms4_527_iso':ms4_527_iso}
                
            case "MAN6":
                isomer=['6H1', '6E1', '6E2', '6G1', '6D2', '6D1', '6D3','6F2','6F1']
                pref=['Man-6_MS2-1419','Man-6_MS3-1013','Man-6_MS3-995','Man-6_MS4-527']
                m_sel=[1419,1013,995,527]
                m_sel_iso=[1419.476,1013.317, 995.307, 527.159]
                ms_level=[2,3,3,4]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[1419,1013,995],'m_sel_iso':[1419.476,1013.317, 995.307]})
                met.append({1:isomer[1],'m_sel':[1419,1013,995],'m_sel_iso':[1419.476,1013.317, 995.307]})
                met.append({2:isomer[2],'m_sel':[1419,1013,995,527],'m_sel_iso':[1419.476,1013.317, 995.307,527.159]})
                met.append({3:isomer[3],'m_sel':[1419,1013,995],'m_sel_iso':[1419.476,1013.317, 995.307]})
                met.append({4:isomer[4],'m_sel':[1419,1013,995],'m_sel_iso':[1419.476,1013.317, 995.307]})
                met.append({5:isomer[5],'m_sel':[1419,1013,995],'m_sel_iso':[1419.476,1013.317, 995.307]})
                met.append({6:isomer[6],'m_sel':[1419,1013,995,527],'m_sel_iso':[1419.476,1013.317, 995.307,527.159]})
                met.append({7:isomer[7],'m_sel':[1419,1013,995],'m_sel_iso':[1419.476,1013.317, 995.307]})
                met.append({8:isomer[8],'m_sel':[1419,1013,995],'m_sel_iso':[1419.476,1013.317, 995.307]})
                #minmass=[1100,500,500,200]
                minmass=[100,100,100,100]
                ms2a=[1401,1318,1257,1216,1198,1095,995]
                ms2a_iso=[1401.465, 1318.428, 1257.423, 1216.397, 1198.386, 1095.37, 995.307]
                ms2b=[1095,933,874,833,771,712,671]
                ms2b_iso=[1095.37, 933.318, 874.28, 833.254, 771.265, 712.228, 671.201]
                ms2c=[1095,1257,1318] ##check 1095,1257,1419 vs 1480
                ms2c_iso=[1095.37, 1257.423, 1318.428]
                ms3_1013=[995,953,923,893,851,761,689,599,527,509,437]
                ms3_1013_iso=[995.307, 953.296, 923.286, 893.275, 851.264, 761.233, 689.212, 599.18, 527.159, 509.148, 437.127]
                ms3_995=[923,851,833,671,599,509]
                ms3_995_iso=[923.286, 851.264, 833.254, 671.201, 599.18, 509.148]
                ms4_527=[509,467,437,407,365,347,275]
                ms4_527_iso=[509.148, 467.138, 437.127, 407.117, 365.106, 347.095, 275.074]
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms2c':ms2c,'ms3_1013':ms3_1013,
                           'ms3_995':ms3_995,'ms4_527':ms4_527,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms2c_iso':ms2c_iso,'ms3_1013_iso':ms3_1013_iso,
                           'ms3_995_iso':ms3_995_iso,'ms4_527_iso':ms4_527_iso}

            case "MAN7":
                isomer=['7E1', '7D3', '7G1', '7E2', '7D2', '7D1', '7F1']
                pref=['Man-7_MS2-1581','Man-7_MS3-1157','Man-7_MS3-1175','Man-7_MS4-527']
                m_sel=[1581,1157,1175,527]
                m_sel_iso=[1581.529, 1157.36, 1175.37, 527.159]
                ms_level=[2,3,3,4]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[1581,1157,1175],'m_sel_iso':[1581.529, 1157.36, 1175.37]})
                met.append({1:isomer[1],'m_sel':[1581,1157,1175,527],'m_sel_iso':[1581.529, 1157.36, 1175.37,527.159]})
                met.append({2:isomer[2],'m_sel':[1581,1157,1175],'m_sel_iso':[1581.529, 1157.36, 1175.37]})
                met.append({3:isomer[3],'m_sel':[1581,1157,1175],'m_sel_iso':[1581.529, 1157.36, 1175.37]})
                met.append({4:isomer[4],'m_sel':[1581,1157,1175],'m_sel_iso':[1581.529, 1157.36, 1175.37]})
                met.append({5:isomer[5],'m_sel':[1581,1157,1175,527],'m_sel_iso':[1581.529, 1157.36, 1175.37,527.159]})
                met.append({6:isomer[6],'m_sel':[1581,1157,1175],'m_sel_iso':[1581.529, 1157.36, 1175.37]})
                minmass=[100,100,100,100]
                ms2a=[1563,1480,1419,1378,1360,1257,1095,1157]
                ms2a_iso=[1563.518, 1480.481, 1419.476, 1378.449, 1360.439, 1257.423, 1095.37, 1157.36]
                ms2b=[1257,1095,1036,995,933,874,833,771]
                ms2b_iso=[1257.423, 1095.37, 1036.333, 995.307, 933.318, 874.28, 833.254, 771.265]
                ms2c=[1095,1257,1419,1480] ##check 1095,1257,1419 vs 1480
                ms2c_iso=[1095.37, 1257.423, 1419.476, 1480.481]
                ms3_1175=[1157,1013,923,851,761,689,599,527]
                ms3_1175_iso=[1157.36, 1013.317, 923.286, 851.264, 761.233, 689.212, 599.18, 527.159]
                ms3_1157=[995,833,671,509]
                ms3_1157_iso=[995.307, 833.254, 671.201, 509.148]
                ms4_527=[509,467,437,407,365,347,275]
                ms4_527_iso=[509.148, 467.138, 437.127, 407.117, 365.106, 347.095, 275.074]
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms2c':ms2c,'ms3_1175':ms3_1175,
                           'ms3_1157':ms3_1157,'ms4_527':ms4_527,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms2c_iso':ms2c_iso,'ms3_1175_iso':ms3_1175_iso,
                           'ms3_1157_iso':ms3_1157_iso,'ms4_527_iso':ms4_527_iso}
                #minmass=[400,400,400,200]

            case "MAN8":
                isomer=['8E1', '8G1', '8E2', '8D3']
                pref=['Man-8_MS2-1743','Man-8_MS3-1337','Man-8_MS3-1319']
                m_sel=[1743,1337,1319]
                m_sel_iso=[1743.582, 1337.423, 1319.412]
                ms_level=[2,3,3]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[1743,1337,1319],'m_sel_iso':[1743.582, 1337.423, 1319.412]})
                met.append({1:isomer[1],'m_sel':[1743,1337,1319],'m_sel_iso':[1743.582, 1337.423, 1319.412]})
                met.append({2:isomer[2],'m_sel':[1743,1337,1319],'m_sel_iso':[1743.582, 1337.423, 1319.412]})
                met.append({3:isomer[3],'m_sel':[1743,1337,1319],'m_sel_iso':[1743.582, 1337.423, 1319.412]})
                minmass=[100,100,100,100]
                ms2a=[1725,1642,1581,1540,1522,1319,1419,1257]
                ms2b=[1257,1095,1036,995,933,874]
                ms2c=[1581,1419,1257,1642] ##check 1095,1257,1419 vs 1480
                ms3_1337=[1319,1175,1013,923,851,761,689,599]
                ms3_1319=[1157,995,833,671]
                
                ms2a_iso=[1725.571, 1642.534, 1581.529, 1540.502, 1522.492, 1319.412, 1419.476, 1257.423]
                ms2b_iso=[1257.423, 1095.37, 1036.333, 995.307, 933.318, 874.28]
                ms2c_iso=[1581.529, 1419.476, 1257.423, 1642.534] ##check 1095,1257,1419 vs 1480
                ms3_1337_iso=[1319.412, 1175.37, 1013.317, 923.286, 851.264, 761.233, 689.212,  599.18]
                ms3_1319_iso=[1157.36, 995.307, 833.254, 671.201]
            
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms2c':ms2c,'ms3_1337':ms3_1337,
                           'ms3_1319':ms3_1319,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms2c_iso':ms2c_iso,'ms3_1337_iso':ms3_1337_iso,
                           'ms3_1319_iso':ms3_1319_iso}


            case "MAN9":
                isomer=['9E1', '9D1', '9D2']
                pref=['Man-9_MS2-1905','Man-9_MS3-1499','Man-9_MS3-1481']
                m_sel=[1905,1499,1481]
                m_sel_iso=[1905.635, 1499.476, 1481.465]
                ms_level=[2,3,3]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[1905,1499,1481],'m_sel_iso':[1905.635, 1499.476, 1481.465]})
                met.append({1:isomer[1],'m_sel':[1905,1499,1481],'m_sel_iso':[1905.635, 1499.476, 1481.465]})
                met.append({2:isomer[2],'m_sel':[1905,1499,1481],'m_sel_iso':[1905.635, 1499.476, 1481.465]})
                minmass=[100,100,100]
                ms2a=[1887,1804,1743,1702,1684,1581,1419,1481]
                ms2b=[1419,1257,1198,1157,1095,1036]
                ms2c=[1743,1581,1419,1804] ##check 1095,1257,1419 vs 1480
                ms3_1499=[1481,1337,1175,1085,1013,923,851,761]
                ms3_1481=[1319,1157,995,833]
                
                ms2a_iso=[1887.624, 1804.587, 1743.582, 1702.555, 1684.545, 1581.529, 1419.476, 1481.465]
                ms2b_iso=[1419.476, 1257.423, 1198.386, 1157.36, 1095.37, 1036.333]
                ms2c_iso=[1743.582, 1581.529, 1419.476, 1804.587] ##check 1095,1257,1419 vs 1480
                ms3_1499_iso=[1481.465, 1337.423, 1175.37, 1085.338, 1013.317, 923.286, 851.264, 761.233]
                ms3_1481_iso=[1319.412, 1157.36, 995.307, 833.254]

                
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms2c':ms2c,'ms3_1499':ms3_1499,
                           'ms3_1481':ms3_1481,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms2c_iso':ms2c_iso,'ms3_1499_iso':ms3_1499_iso,
                           'ms3_1481_iso':ms3_1481_iso}


            case "MAN10":
                isomer=['10D1']
                pref=['Man-10_MS2-1045','Man-10_MS3-1661','Man-10_MS3-1643']
                m_sel=[1045,1661,1643]
                m_sel_iso=[1045.339, 1661.529, 1643.518]
                ms_level=[2,3,3]
                ion_sel=[]
                met.append({0:isomer[0],'m_sel':[1045,1661,1643],'m_sel_iso':[1045.339, 1661.529, 1643.518]})
                minmass=[100,100,100]
                ms2a=[1036,994.5,964,943.5,934.5,883,802,833]
                ms2b=[802,721,691.5,671,640,610.5]
                ms2c=[964,883,802,994.5] ##check 1095,1257,1419 vs 1480
                ms3_1661=[1643,1499,1337,1247,1175,1085,1013,923]
                ms3_1643=[1481,1319,1157,995]
                
                
                ms2a_iso=[1036.333, 994.815, 964.312, 943.799, 934.794, 883.286, 802.259, 833.254]
                ms2b_iso=[802.259, 721.233, 691.714, 671.201, 640.206, 610.688]
                ms2c_iso=[964.312, 883.286, 802.259, 994.815] ##check 1095,1257,1419 vs 1480
                ms3_1661_iso=[1643.518, 1499.476, 1337.423, 1247.391, 1175.37, 1085.338, 1013.317, 923.286]
                ms3_1643_iso=[1481.465, 1319.412, 1157.36, 995.307]
                
                peaks_sel={'ms2a':ms2a,'ms2b':ms2b,'ms2c':ms2c,'ms3_1661':ms3_1661,
                           'ms3_1643':ms3_1643,'ms2a_iso':ms2a_iso,'ms2b_iso':ms2b_iso,'ms2c_iso':ms2c_iso,'ms3_1661_iso':ms3_1661_iso,
                           'ms3_1643_iso':ms3_1643_iso}


        return met,peaks_sel

    def database(self):
        ngly=self.nglycan
        ngly.upper()
        met=[]
        # pref='prefix'
        match ngly.upper():


            case "MAN1":
                isomer=['Man1']
                pref=['Man-1_MS2-609']
                m_sel=[609]
                ms_level=[2]

            case "MAN2":
                isomer=['2E1','2F1']
                pref=['Man-2_MS2-771','Man-2_MS3-365']
                m_sel=[771,365]
                ms_level=[2,3]

            case "MAN3":
                isomer=['3D1','3F1','3E2','3E1']
                pref=['Man-3_MS2-933','Man-3_MS3-509','Man-3_MS3-527','Man-3_MS4-365']
                m_sel=[933,509,527,365]
                ms_level=[2,3,3,4]

            case "MAN4":
                isomer=['4D3','4D1','4D2','4F1','4E3','4E2','4E1']
                pref=['Man-4_MS2-1095','Man-4_MS3-689','Man-4_MS3-671','Man-4_MS4-365','Man-4_MS4-527']
                m_sel=[1095,689,671,365,527]
                ms_level=[2,3,3,4,4]


            case "MAN5":
                isomer=['5F1', '5E1', '5E3', '5E2', '5D2', '5E4', '5D1','5F2']
                # pref='prefix'
                pref=['Man-5_MS2-1257','Man-5_MS3-851','Man-5_MS3-833','Man-5_MS4-527']
                m_sel=[1257,851,833,527]
                ms_level=[2,3,3,4]
            
            case "MAN6":
                isomer=['6H1', '6E1', '6E2', '6G1', '6D2', '6D1', '6D3','6F2','6F1']
                pref=['Man-6_MS2-1419','Man-6_MS3-1013','Man-6_MS3-995','Man-6_MS4-527']
                m_sel=[1419,1013,995,527]
                ms_level=[2,3,3,4]

            case "MAN7":
                isomer=['7E1', '7D3', '7G1', '7E2', '7D2', '7D1', '7F1']
                pref=['Man-7_MS2-1581','Man-7_MS3-1157','Man-7_MS3-1175','Man-7_MS4-527']
                m_sel=[1581,1157,1175,527]
                ms_level=[2,3,3,4]

            case "MAN8":
                isomer=['8E1', '8G1', '8E2', '8D3']
                pref=['Man-8_MS2-1743','Man-8_MS3-1337','Man-8_MS3-1319']
                m_sel=[1743,1337,1319]
                ms_level=[2,3,3]

            case "MAN9":
                isomer=['9E1', '9D1', '9D2']
                pref=['Man-9_MS2-1905','Man-9_MS3-1499','Man-9_MS3-1481']
                m_sel=[1905,1499,1481]
                ms_level=[2,3,3]

            case "MAN10":
                isomer=['10D1']
                pref=['Man-10_MS2-1045','Man-10_MS3-1661','Man-10_MS3-1643']
                m_sel=[1045,1661,1643]
                ms_level=[2,3,3]



        return pref


    def cal_score1(self,allraw,datab,ret_time_num):

        tim_in,int_in=ret.ms2ret(allraw)
        max_tint1=ret.in_peaks2(tim_in,int_in,0.1)
        ms_level,mass_select=ret.ms_lev(allraw)

        isomers=ret.get_isomer(datab)
        nisomers=len(isomers)
        tim_d,int_d=ret.get_ret(datab)
        max_tint2=ret.database_peaks(nisomers,tim_d,int_d)

        ret_int1=[tim_in,int_in] 
        t_ini,t_end=usf.t_inter(max_tint1,ret_int1,ret_time_num)
#        mass_select,ms_level=self.mass_select()
        rawms_spectra=[]
        i1=0
        for ms_sel in mass_select:
            a,b,c=usf.spectra_profile(allraw,t_ini,t_end,ms_sel)
            ms_peaks=ret.in_peaks2(a,b,0.01)
            time1=max_tint1[ret_time_num].get('time')
            s1={'ms_level':ms_level[i1],'selected_mass':ms_sel,'spectra':[a,b,c],'ms_peaks':ms_peaks,'Ret_peak':time1}
            i1=i1+1
            rawms_spectra.append(s1)
            #plt.figure(i+1)
            #plt.plot(a,b)
        return rawms_spectra

    def cal_score2(self,allraw,mass,datab,ret_time_num):

        tim_in,int_in=ret.ms2ret_mass(allraw,mass)
        max_tint1=ret.in_peaks2(tim_in,int_in,0.1)
        ms_level,mass_select=ret.ms_lev_mass(allraw,mass)
        #print(ms_level,mass_select)
        isomers=ret.get_isomer(datab)
        nisomers=len(isomers)
        tim_d,int_d=ret.get_ret(datab)
        max_tint2=ret.database_peaks(nisomers,tim_d,int_d)

        ret_int1=[tim_in,int_in]
        t_ini,t_end=usf.t_inter(max_tint1,ret_int1,ret_time_num)
#        mass_select,ms_level=self.mass_select()
        rawms_spectra=[]
        i1=0
        for ms_sel in mass_select:
            a,b,c=usf.spectra_profile1(allraw,mass,t_ini,t_end,ms_sel)
            ms_peaks=ret.in_peaks2(a,b,0.01)
            time1=max_tint1[ret_time_num].get('time')
            s1={'ms_level':ms_level[i1],'selected_mass':ms_sel,'spectra':[a,b,c],'ms_peaks':ms_peaks,'Ret_peak':time1}
            i1=i1+1
            rawms_spectra.append(s1)
            #plt.figure(i+1)
            #plt.plot(a,b)
        return rawms_spectra


    def ms_sim3(self,iso_list,in_spec,datab):
        #initialize database,methods
        pref=self.database()
        met,peaks_sel=self.method()
        #print(len(in_spec))
        #['ms_level', 'selected_mass', 'spectra', 'ms_peaks', 'Ret_peak']
        ## check input whether they got
        score=[]
        method=met[iso_list]
        isomer=method.get(iso_list)
        sc=[]
        for i in range(len(in_spec)):
            ml=in_spec[i].get('ms_level')
            ma=in_spec[i].get('selected_mass')
            mspect=in_spec[i].get('spectra')

            m_sel_d=method.get('m_sel')
            index=m_sel_d.index(ma)
            a,b=usf.get_msspectra(datab,isomer,pref[index])
            spec_dat=[a,b]
     
            max1=ret.in_peaks3(a,b,0.05)
            #print(peaks_sel.get('ms2a'))
            if ml=='2':
                max2=ret.in_peaks4(a,b,peaks_sel.get('ms2a'))
                scr1=usf.scr_1(mspect,spec_dat,max2)
                sc.append(['ms2_part1',scr1])
                max2=ret.in_peaks4(a,b,peaks_sel.get('ms2b'))
                scr1=usf.scr_1(mspect,spec_dat,max2)
                sc.append(['ms2_part2',scr1])

            else:    
                str_peak='ms'+ml+'_'+str(ma)
                #print(str_peak,ml,ma)
                max2=ret.in_peaks4(a,b,peaks_sel.get(str_peak))
                scr1=usf.scr_1(mspect,spec_dat,max2)
                sc.append([str_peak,scr1])

        ssc=1
        for a in sc:
            ssc=ssc*a[1]
        sc1={'isomer':isomer,'MSn_score':sc,'multi_score':ssc}
        return sc1

    def ms_sim4(self,iso_list,in_spec,datab,score_calc,mass_tol):
        #initialize database,methods
        pref=self.database()
        met,peaks_sel=self.method()
        #print(len(in_spec))
        #['ms_level', 'selected_mass', 'spectra', 'ms_peaks', 'Ret_peak']
        ## check input whether they got
        score=[]
        method=met[iso_list]
        isomer=method.get(iso_list)
        sc=[]
        sc_table=[]
        for i in range(len(in_spec)):
            ml=in_spec[i].get('ms_level')
            ma=in_spec[i].get('selected_mass')
            mspect=in_spec[i].get('spectra')

            m_sel_d=method.get('m_sel')
            index=m_sel_d.index(ma)
            a,b=usf.get_msspectra(datab,isomer,pref[index])
            spec_dat=[a,b]
     
            max1=ret.in_peaks3(a,b,0.05)
            #print(peaks_sel.get('ms2a'))
            if ml=='2':
                max2=peaks_sel.get('ms2a_iso') 
                if score_calc==1:
                    scr1,tabl=usf.scr_1(mspect,spec_dat,max2,mass_tol)
                if score_calc==2:
                    scr1,tabl=usf.scr_2(mspect,spec_dat,max2,mass_tol)
                if score_calc==3:
                    scr1,tabl=usf.scr_3(mspect,spec_dat,max2,mass_tol)
                sc.append(['ms2_part1',scr1])
                tabl.append('ms2_part1')
                sc_table.append(tabl)
                max2=peaks_sel.get('ms2b_iso')
                if score_calc==1:
                    scr1,tabl=usf.scr_1(mspect,spec_dat,max2,mass_tol)
                if score_calc==2:
                    scr1,tabl=usf.scr_2(mspect,spec_dat,max2,mass_tol)
                if score_calc==3:
                    scr1,tabl=usf.scr_3(mspect,spec_dat,max2,mass_tol)
                sc.append(['ms2_part2',scr1])
                tabl.append('ms2_part2')
                sc_table.append(tabl)
                
                try:
                    max2=peaks_sel.get('ms2c_iso')
                    if score_calc==1:
                        scr1,tabl=usf.scr_1(mspect,spec_dat,max2,mass_tol)
                    if score_calc==2:
                        scr1,tabl=usf.scr_2(mspect,spec_dat,max2,mass_tol)
                    if score_calc==3:
                        scr1,tabl=usf.scr_3(mspect,spec_dat,max2,mass_tol)
                    sc.append(['ms2_part3',scr1])
                    tabl.append('ms2_part3')
                    sc_table.append(tabl)
                except:
                    a=1

            else:
                str_peak='ms'+ml+'_'+str(ma)+'_iso'
                max2=peaks_sel.get(str_peak)
                if score_calc==1:
                    scr1,tabl=usf.scr_1(mspect,spec_dat,max2,mass_tol)
                if score_calc==2:
                    scr1,tabl=usf.scr_2(mspect,spec_dat,max2,mass_tol)
                if score_calc==3:
                    scr1,tabl=usf.scr_3(mspect,spec_dat,max2,mass_tol)
                sc.append([str_peak,scr1])
                tabl.append(str_peak)
                sc_table.append(tabl)
        ssc=1
        for a in sc:
            ssc=ssc*a[1]
        sc1={'isomer':isomer,'MSn_score':sc,'multi_score':ssc,'spectra_table':sc_table}
        return sc1


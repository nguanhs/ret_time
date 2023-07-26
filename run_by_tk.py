import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
import subprocess
import os
from src import main_progcheck3 as pro1
import sys


if len(sys.argv) < 1:
    print("No input file provided.")
    exit()

raw_file=sys.argv[1]
mass_tol=float(sys.argv[2])
met1=[];met1.append(sys.argv[3]);met1.append(sys.argv[4]);met1.append(sys.argv[5])


#f1=open(finput,'r')

'''
    if len(a.split())>1:
        print("Error in file's name format.")
        print("Please make sure:\n 1. One file name per line.\n 2. Not white space in the file's name")
        break
    raw_file=a.split()[0]
'''
    
allraw=pro1.raw_input(raw_file)
n_spect,hman=pro1.raw_hman_info(allraw)
print(n_spect,hman)
#scoring_scheme=1 for "integrated_difference"
#scoring_scheme=2 for "Dot_product"
#scoring_scheme=3 for "DK_divergence"
i=1
for a in met1:
    if a=='1':
        match,results1=pro1.match_cal(allraw,hman,i,mass_tol)
        pro1.present_results(hman,match,results1,allraw,i,raw_file)
    i=i+1





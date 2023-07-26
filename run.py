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

finput=sys.argv[1]

f1=open(finput,'r')

for a in f1:
    if len(a.split())>1:
        print("Error in file's name format.")
        print("Please make sure:\n 1. One file name per line.\n 2. Not white space in the file's name")
        break
    raw_file=a.split()[0]
    
    allraw=pro1.raw_input(raw_file)
    n_spect,hman=pro1.raw_hman_info(allraw)
    print(n_spect,hman)
    #scoring_scheme=1 for "integrated_difference"
    #scoring_scheme=2 for "Dot_product"
    #scoring_scheme=3 for "DK_divergence"

    match,results1=pro1.match_cal(allraw,hman,1)
    match,results2=pro1.match_cal(allraw,hman,2)
    match,results3=pro1.match_cal(allraw,hman,3)

    pro1.present_results(hman,match,results1,allraw,1,raw_file)
#    pro1.present_results(hman,match,results2,allraw,2,raw_file)
#    pro1.present_results(hman,match,results3,allraw,3,raw_file)




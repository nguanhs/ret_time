from tkinter.constants import CENTER
from tkinter.messagebox import *
import tkinter as tk
from tkinter import filedialog
import tkinter.colorchooser as cc
from PIL import Image,ImageDraw,ImageFont
from tkinter import ttk

import subprocess
import os
from src import main_progcheck3 as pro1
import sys


def locate():
#    file_path = filedialog.askopenfilename(filetypes=[("COM and Log Files", "*.com *.log")])
    file_path = filedialog.askopenfilename()
    if file_path:
        with open(file_path, "r") as file:
            for a in file:
                print(a)
            #file_contents = file.read()
            #text_box.insert(tk.END, file_contents)
            #print(a)

def loadFile():
    if loadFile_en.get() is None:
        file_path = filedialog.askopenfilename(filetypes = [("mzml and raw Files", "*.mzML *.raw")])
        loadFile_en.insert(0,file_path)
    else:
        file_path = filedialog.askopenfilename(filetypes =[("mzml and raw Files", "*.mzML *.raw")])
        loadFile_en.delete(0,'end')
        loadFile_en.insert(0,file_path)

def handle_selection(event):
    selected_item = dropdown.get()
    print("Selected item:", selected_item)


def on_checkbox_select():
    selected_choices = [choice_var.get() for choice_var in choice_vars]
    print("Selected choices:", selected_choices)

def get_entry_value():
    value = mass_tol.get()
    print("Entered value:", value)

def run():
    com1=' '
    com2=[]
    filename=loadFile_en.get()
    if filename:
        com1=com1+' '+filename
        com2.append(filename)
    else:
        print("Input raw not provided. Please enter a raw file.")

    mt = mass_tol.get()
    if mt:
        com1=com1+' '+str(mt)
        com2.append(mt)
    else:
        print("Machine mass accurcy not provided. Please enter the value.")

    selected_choices = [choice_var.get() for choice_var in choice_vars]
    if selected_choices:
        com1=com1+' '+ str(selected_choices[0])+' '+ str(selected_choices[1])+' '+ str(selected_choices[2])
        com1a=str(selected_choices[0])+' '+ str(selected_choices[1])+' '+ str(selected_choices[2])
        com2.append(str(selected_choices[0]));com2.append(str(selected_choices[1]));com2.append(str(selected_choices[2]))
    else:
        print("Similarity score's scheme not provided. Please select at least one scheme.")
 
    print(com2)

    if filename and mt and selected_choices:
#        subprocess.run(["python3", com1]) 
        subprocess.run(["python3", 'run_by_tk.py']+com2)

# Create the Tkinter window

window = tk.Tk()
window.title('LODES program')
window.geometry('600x400+50+50')

'''
frame1 = tk.Frame(window)
frame1.pack(side="top")


frame2 = tk.Frame(window)
frame2.pack(side="top")

frame3 = tk.Frame(window)
frame3.pack(side="top")
'''



# Create a text box widget
#text_box = tk.Text(window)
#text_box.pack()

# Create the Tkinter window

# Create a label
label1 = tk.Label(text="Select file:",height=1)
label1.place(x=0 ,y=0)
#label1.pack(side='left')


label2 = tk.Label(text="Machine mass accuracy(in m/z):",height=1)
label2.place(x=0 ,y=50)
#label2.pack(side='left')


label3 = tk.Label( text="Select Similarity Score Scheme(or Schemes):",height=1)
label3.place(x=0 ,y=130)
#label3.pack(side='left')


# Select raw file:
loadFile_en = tk.Entry(width=40)
loadFile_en.place(x=70 ,y=0)

loadFile_btn = tk.Button(text="search",height=1,command=loadFile)
loadFile_btn.place(x=400 ,y=0)


# Entry widget to capture user input
mass_tol = tk.Entry(width=10)
mass_tol.place(x=20 ,y=100)



choices = ["Integrated Diff.", "Dot product", "KL Divergence"]
choice_vars = []

i=0
for choice in choices:
    var = tk.IntVar()
    choice_vars.append(var)
    checkbox = tk.Checkbutton(window, text=choice, variable=var, command=on_checkbox_select)
    x1=[0,100,200]
    checkbox.place(x=x1[i] ,y=160)
    i=i+1

# Button to retrieve the entered value
button = tk.Button(window, text="RUN", command=run)
button.place(x=200 ,y=180)



window.mainloop()


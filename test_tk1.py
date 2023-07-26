from tkinter.constants import CENTER
from tkinter.messagebox import *
import tkinter as tk
from tkinter import filedialog
import tkinter.colorchooser as cc
from PIL import Image,ImageDraw,ImageFont

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
        file_path = filedialog.askopenfilename(filetypes = (("png files","*.png"),("all files","*.*")))
        loadFile_en.insert(0,file_path)
    else:
        file_path = filedialog.askopenfilename(filetypes = (("png files","*.png"),("all files","*.*")))
        loadFile_en.delete(0,'end')
        loadFile_en.insert(0,file_path)

def handle_selection(event):
    selected_item = dropdown.get()
    print("Selected item:", selected_item)

# Create the Tkinter window

window = tk.Tk()
window.title('LODES program')
window.geometry('600x400+50+50')

frame1 = tk.Frame(window)
frame1.pack(side="top")


frame2 = tk.Frame(window)
frame2.pack(side="top")

frame3 = tk.Frame(window)
frame3.pack(side="top")




# Create a text box widget
#text_box = tk.Text(window)
#text_box.pack()

# Create the Tkinter window

# Create a label
label1 = tk.Label(frame1, text="Select file:")
label1.pack(side='left')

label2 = tk.Label(frame2, text="Machine mass accuracy")
label2.pack(side='left')


label3 = tk.Label(frame3, text="Select Similarity Score Scheme(s)")
label3.pack(side='left')


# Select raw file:
loadFile_en = tk.Entry(width=40)
loadFile_en.place(x=70 ,y=0)




# Create a dropdown selection widget
dropdown = ttk.Combobox(frame1, values=["Item 1", "Item 2", "Item 3"])
dropdown.pack(side='left')
#dropdown.grid(row=0, column=1)

# Set a default value for the dropdown
dropdown.current(0)

# Create a button to open the file
open_button = tk.Button(frame1, text="Loactes file", command=locate)
open_button.pack(side='left',ipadx=10,ipady=10)
#open_button.grid(row=0, column=0)


# exit button
exit_button = tk.Button(
    frame2,
    text='Exit',
    command=lambda: window.quit()
)

exit_button.pack(side='left')
#exit_button.grid(row=1, column=0,columnspan=2)

# Start the Tkinter event loop
window.mainloop()


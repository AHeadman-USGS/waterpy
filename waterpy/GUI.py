
import os
import subprocess
from waterpy.main import waterpy
import PIL
from PIL import Image, ImageTk
from waterpy import geospatial

#import everything from tkinter
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import tkinter.messagebox


#Create Window object
window=Tk()
window.title("WATERpy")
window.iconbitmap('usgs-logo-green.ico')
w = 525 # width for the Tk root
h = 500 # height for the Tk root

# get screen width and height
ws = window.winfo_screenwidth() # width of the screen
hs = window.winfo_screenheight() # height of the screen

# calculate x and y coordinates for the Tk root window
x = (ws/2) - (w/2)
y = (hs/2) - (h/2)

# set the dimensions of the screen 
# and where it is placed
window.geometry('%dx%d+%d+%d' % (w, h, x, y))
    


frame = tkinter.Frame(window)
frame.configure(bg ="#A4C1A4")
frame.grid(row=0, column=0, sticky='nswe')
frame.grid_rowconfigure(0, weight=1)
frame.grid_columnconfigure(1, weight=1)
window.grid_rowconfigure(0, weight=1)
window.grid_columnconfigure(0, weight=1)

frame.grid_rowconfigure(9, minsize = 20)
frame.grid_rowconfigure(0, minsize = 10)
frame.grid_rowconfigure(2, minsize =20)
frame.grid_columnconfigure(0, minsize = 20)
frame.grid_columnconfigure(3, minsize = 20)
frame.grid_columnconfigure(5, minsize = 20)

#Tabs

tabControl = ttk.Notebook(frame, height = 400, width = 450)
tab1 = ttk.Frame(tabControl)
tab2 = ttk.Frame(tabControl)
tabControl.add(tab1, text = 'Geospatial')
tabControl.add(tab2, text = 'WATERpy')
tabControl.grid(row = 2, column = 2)



filename = ""  
inifilename = ""
folderpath2 = ""
folderpath = ""


L1 = tkinter.Label(frame, text = "Water Availability Tool for Environmental Resources", width = 40, height = 1, relief= 'raised')
L1.place(x=275, y=25, anchor = 'center')

    

#GEOSPATIAL TAB
CheckVar1 = tkinter.IntVar()
CheckVar2 = tkinter.IntVar()
    
def run_model():
    print(filename)
    print(folderpath)
    if CheckVar1.get() == 1:
        geospatial.geospatial(str(filename), str(folderpath), ts = True)
        tkinter.messagebox.showinfo('Waterpy', 'Geospatial Compiling Sucessful. If you would like to run another basin, please close and relaunch the application.')
    else:
        geospatial.geospatial(str(filename), str(folderpath), ts = False)
        tkinter.messagebox.showinfo('Waterpy', 'Geospatial Compiling Sucessful. If you would like to run another basin, please close and relaunch the application.')
    

Intro = tkinter.Label(tab1, text = "Geospatial.py requires a user-supplied shapefile of a basin. The program will use this to compile three files: Basin Characteristics, Time Series, and TWI.", height=2, width = 100, wraplength = 420, justify = "center")
Intro.place(x = -130, y = 25, anchor = "w")   

def getFilePath():
    global filename
    filename = filedialog.askopenfilename()
    shpfilePath.set(filename)  
    if filename:
        filepath= os.path.abspath(filename)
    return filepath

def FilePath():
    getFilePath()
    
def getFolderPath():
    global folderpath
    dbfolder_selected = filedialog.askdirectory()
    dbFolderPath.set(dbfolder_selected)
    if dbfolder_selected:
        folderpathtemp = os.path.abspath(dbfolder_selected)
        folderpath = os.path.dirname(folderpathtemp)
    return folderpath
    
def FolderPath():
    getFolderPath()

  
    
shpfilePath = StringVar()  
E = Entry(tab1,textvariable=shpfilePath, width = 35)
E.place(x = 70, y = 95)
btnFind = ttk.Button(tab1, text="Browse",command=FilePath)
btnFind.place(x = 300, y = 92)

T = tkinter.Label(tab1, text = "Shapefile for geospatial sampling", height=1, bg = "#A4C1A4")
T.place(x=135, y =65, anchor = 'w')

ts = tkinter.Label(tab1, text = "Create a time series?", height = 1, bg = "#A4C1A4")
ts.place(x = 170, y = 145, anchor = 'w')

C1 = tkinter.Checkbutton(tab1, text = "Yes, I need temperature and precipitation data ", variable = CheckVar1, anchor = 'w', onvalue= 1, state = 'active', offvalue = 0, height=1)
C1.place(x=220, y=175, anchor = 'center')   

C2 = tkinter.Checkbutton(tab1, text = "No, I already have temperature (Celsius) and precipitation (mm/day) data", variable = CheckVar2, anchor = 'w', onvalue= 1, state = 'active', offvalue = 0, height=1)
C2.place(x=220, y=200, anchor = 'center')  

p = tkinter.Label(tab1, text = "*If I include daily potential-evapotranspiration data, I understand that it will be used in the simulation.  Otherwise, it will be estimated using the Hamon equation", height = 3, wraplength = 380)
p.place(x = 220, y = 235, anchor = 'center')

B1=tkinter.Button(tab1, text = "Run", width = 12, height=1, command = run_model)
B1.place(x=220, y=380, anchor = 'center')

D = tkinter.Label(tab1, text = "Path to Database Folder", height=1, bg = "#A4C1A4")
D.place(x=160, y =285, anchor = 'w')

dbFolderPath = StringVar()  
E3 = Entry(tab1,textvariable=dbFolderPath, width = 35)
E3.place(x = 70, y = 310)
btnFind3 = ttk.Button(tab1, text="Browse", command = FolderPath)
btnFind3.place(x = 300, y = 307)



# WATERpy Tab

def getiniPath():
    global inifilename
    inifilename = filedialog.askopenfilename()
    inifilepath.set(inifilename)
    if inifilename:
        filepath2 = os.path.abspath(inifilename)
    return filepath2

def FilePath2():
    getiniPath()
    
def getFolderPath2():
    global folderpath2
    outfolder_selected = filedialog.askdirectory()
    outfolderPath.set(outfolder_selected)
    if outfolder_selected:
        folderpath2= os.path.abspath(outfolder_selected)
    return folderpath2
    
def FolderPath2():
    getFolderPath2()

def run_water():
    # print(inifilename)
    ini = os.path.dirname(inifilename)
    os.chdir(ini)
    subprocess.call("waterpy run modelconfig.ini")
    tkinter.messagebox.showinfo('Waterpy', 'Model Run Sucessful')  
    os.startfile(folderpath2)
    
    
    
Intro2 = tkinter.Label(tab2, text = "Before running WATERpy, the modelconfig.ini file must be updated with the paths to your Input, Output, and Database folders", height=2, width = 100, wraplength = 420, justify = "center")
Intro2.place(x = -130, y = 25, anchor = "w") 

m = tkinter.Label(tab2, text = "Path to Modelconfig.ini file", height=1, bg = "#A4C1A4")
m.place(x=140, y =85, anchor = 'w')
inifilepath = StringVar()   
E2 = Entry(tab2,textvariable=inifilepath, width = 35)
E2.place(x = 70, y = 105)
btn2Find = ttk.Button(tab2, text="Browse",command=FilePath2)
btn2Find.place(x = 300, y = 102)


j = tkinter.Label(tab2, text = "Path to Outputs Folder", height=1, bg = "#A4C1A4")
j.place(x=155, y =165, anchor = 'w')
outfolderPath = StringVar()   
E4 = Entry(tab2,textvariable=outfolderPath, width = 35)
E4.place(x = 70, y = 190)
btn4Find = ttk.Button(tab2, text="Browse",command=getFolderPath2)
btn4Find.place(x = 300, y = 188)

B2=tkinter.Button(tab2, text = "Run", width = 12, height=1, command = run_water)
B2.place(x=230, y=275, anchor = 'center')



#USGS Logo

PILFile = PIL.Image.open("usgslogo.jpg")
PILFile = PILFile.resize((65, 40), PIL.Image.ANTIALIAS)
Image = ImageTk.PhotoImage(PILFile) # <---
ImageLabel = Label(frame, image=Image)
ImageLabel.image = Image
ImageLabel.place(x=40, y= 25, anchor = 'center')



window.mainloop()
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:19:41 2020

@author: lucab
"""
'''
This program deals with binary two dimensional closed and isolated cellular automata 

'''

#import packages
import numpy as np
from tkinter import Tk, Canvas, Button, LEFT, TOP, Label
import matplotlib.pyplot as plt
import random
import pandas as pd
import ast
from functions import rulegenerator
#---Imput row to be read in Excel file: "Variables and initial configuration"---#                                                                               #
row = 10                                                               #
#imports the programmable variables
Variables = pd.read_excel('Variables and initial configuration.xlsx', sheet_name='Blad1')  

#taking into account the way the file is read
row = row - 2 

#print(Variables.iloc[row])

#Extracting all variables
NbL=400
#NbL = Variables.iloc[row,0]                                   
NbL=int(NbL)                                        
NbC = NbL                                           #sets dimension of map[square]
a=8
#a = Variables.iloc[row,1]                           #pixel-width of squares in grid
Closed_boundaries = True
#Variables.iloc[row,2]           #closes the boundaries if True, else they are periodic 
birth = [3]
death = [0,1,7,8]
#birth,death=rulegenerator()
density = 0.1
#Variables.iloc[row,5]                     #the initial density of alive cells
see_local_s = False
# Variables.iloc[row,6]                #True or False, should activity be counted in a specific area
N=NbL*NbC#total number of cells
E=density*N#initial energy
if see_local_s == True:                             #If True, the coordinates of the rectangular area are needed
    x1, y1, x2, y2 = ast.literal_eval(Variables.iloc[row,7])
    localarea=(x2-x1)*(y2-y1)

if Closed_boundaries == False or Closed_boundaries== True: 
    modify = 1
else:
    raise Exception('"Closed_boundaries" should take argument True or False. Given argument for "Closed_boundaries" was: {}'.format(Closed_boundaries))


#Creation of matrices
StateMatrix = np.zeros((NbL,NbC),dtype=int) #current state
SM1 = np.zeros((NbL,NbC),dtype=int) #update matrix

#state definition
state = {"Dead":0,"Alive":1} #helps with readability


#variable creation
flag=0 #flag defines the status of the automaton: 0=inactive, 1=active
gen=0 #generation counter
s = 0 #current activity
P=0#pressure 
run = 1 #to take resents into account
ActivityList = [] #lists for activity, energy, entropy, temperature and pressure for plots
EnergyList=[]
S=[]
T=[]
PressureList=[]
root = Tk()                 #name change of widget tool
counter = Label(root)       #label for generation counter
activity = Label(root)      #necessary for interface
energy=Label(root)
entropy=Label(root)
localActivity = Label(root)
temperature= Label(root)
pressure=Label(root)

# Calculating next generation

def iterate():      #checks the value of the flag and directs the program to next step, such as rules()
    global flag     #global lets you call a global variable as local
    global gen      #in order to update the generation-counter
    global run
    draw()
    counter.config(text="Gen " + str(gen) + " Run: " + str(run))
                
    if flag==1: #mod2 addition for flag
        gen += 1
        rules()
        root.after(10, iterate)     #waits 10 ms, then runs itearate() again
    elif flag==2:
        gen += 1
        rules()
        flag=0                       #set to inactive 
        root.after(10, iterate)
    elif flag==3:
        run+=1
        initialize_map()
        rules()
        flag=1
        root.after(10, iterate)
    else:
        flag=0

#Initializing the map
def initialize_map():       #run upon starting the program and used to reset the map
    global gen
    global run
    global flag
    global ActivityList
    global localActList
    global E
    global EnergyList
    global S
    global T
    global PressureList
    global h
    T=[]
    PressureList=[]
    flag = 0                #initialize inactive
    
    #counter labels 
    gen = 0
    s = 0
    local_s = 0
    h=(-((E/N)*np.log((E/N))+(1-(E/N))*np.log(1-(E/N))))
    #h=N*np.log(E/N)+N*np.log(N)
    ActivityList=[]
    localActList=[]
    counter.config(text="Gen " + str(gen) + " Run: " + str(run))
    activity.config(text="Average Activity " + str(s))
    localActivity.config(text="Local activity " + str(local_s))
    energy.config(text="Energy "+str(E))
    entropy.config(text="Entropy "+str(h))
    pressure.config(text="Pressure "+str(P))
    #set initial density 
    if density == 0:     
        StateMatrix[0:NbL,0:NbC] = state["Dead"]
        for x in range(NbL):
            for y in range(NbC):
                SM1[x,y] = canvas.create_rectangle((x*a, y*a, (x+1)*a, (y+1)*a), outline="gray", fill="white")            
    elif 0 <= density < 1:
        StateMatrix[0:NbL,0:NbC] = state["Dead"]
        for x in range(NbL):
            for y in range(NbC):
                num=random.random()
                if num<density:
                    StateMatrix[x,y] = state["Alive"]
                SM1[x,y] = canvas.create_rectangle((x*a, y*a, (x+1)*a, (y+1)*a), outline="gray", fill="white")
    elif density == 1:
        StateMatrix[0:NbL,0:NbC] = state["Alive"]
        for x in range(NbL):
            for y in range(NbC):
                SM1[x,y] = canvas.create_rectangle((x*a, y*a, (x+1)*a, (y+1)*a), outline="gray", fill="white")            
    else:
        raise Exception('"density" should be between 0 and 1. The value of "density" was: {}'.format(density))
    
    draw()
    if see_local_s == True:
        canvas.tag_raise(canvas.create_rectangle(x1*a, y1*a, x2*a, y2*a, outline="red", width=2))
    E=density*N #sets initial energy and entropy
    EnergyList=[E]
    S=[h]
# Defining rules of propagation
def rules():  #runs after iterate() if flag=1,2. Creates updated matrix where each cell gets assigned its state for the next gen
    global gen
    global E
    
    E=0
    NEWmatrix = np.zeros((NbL,NbC),dtype=int)
    s = 0
    local_s = 0
    #update ruels
    for i in range(NbL):
            for j in range(NbC):
                global StateMatrix
                n = Moore(i,j,StateMatrix,NbC,NbL,Closed_boundaries); #can call different neighbourhoods (Moore,Moore9 and VN)
                if n in birth:                
                    NEWmatrix[i][j] = 1                    #birth
                elif n in death:         
                    NEWmatrix[i][j] = 0                    #death
                else:                               
                    NEWmatrix[i][j] = StateMatrix[i][j]    #unaffected

    SM1 = StateMatrix.copy()
    StateMatrix = NEWmatrix.copy()          #updates the old matrix to the new
    
    #activity counter
    for i in range(NbL):
            for j in range(NbC):
                if SM1[i][j] != StateMatrix[i][j]:
                    s+=1
    #energy counter  
    for i in range(NbL):
        for j in range(NbC):
            if StateMatrix[i][j]==1:
                E+=1
    #updates for the lists in order to save all values for the plots
    global EnergyList
    global S
    global T
    global h
    #h=N*np.log(E/N)+N*np.log(N)
    h=(-((E/N)*np.log((E/N))+(1-(E/N))*np.log(1-(E/N))))
    EnergyList.append(E)
    energy.config(text="Energy " + str(E))  
    #analytic value for entropy in this specific case, will be generalized later            
    S.append(h) 
    entropy.config(text="Entropy "+str(h))
    #discrete derivative (dE/dS)=T
    t=(EnergyList[-1]-EnergyList[-2])/(S[-1]-S[-2])
    T.append(t)
    temperature.config(text="Temperatue "+str(t))
    
    Z=np.exp(-E/t)#partition function in this specific case, will be used later
    global ActivityList
    
    ActivityList.append(s/N)                         #activity list
    activity.config(text="Average Activity " + str(s/N))      
    
    global P
    global PressureList
    
    P=0
    for i in range (NbL): #counts the number of active cells at the boundaries
    #in analogy to the microscopic definition of pressure
        if StateMatrix[i,0]==1:
            P=P+1
        if StateMatrix[i,NbL-1]==1:
            P=P+1
        if StateMatrix[0,i]==1:
            P=P+1
        if StateMatrix[NbL-1,i]==1:
            P=P+1
    PressureList.append(P)
    pressure.config(text="Pressure "+str(P))
    
    #local activity counter
    if see_local_s == True:
        for i in range(x1,x2):
                for j in range(y1,y2):
                    if SM1[i][j] != StateMatrix[i][j]:
                        local_s += 1
    
        localActList.append(local_s/localarea)               #local activity list
        localActivity.config(text="Local activity " + str(local_s/localarea))



# Drawing the cells (assigning their colors)
def draw():         #called every time the grid should be colored, assigns a color to each state 
    for x in range(NbL):
        for y in range(NbC):
            if StateMatrix[x,y]==state["Dead"]:
                coul = "white"
            elif StateMatrix[x,y]==state["Alive"]:
                coul = "black"
            canvas.itemconfig(SM1[x,y], fill=coul)  
    
def Moore(posX, posY,StateMatrix,NbC,NbL,Closed_boundaries):      #called by rules(), counts the live neighbors to each cell and gives it as output
    N=0    
    temp = StateMatrix[posX][posY]
    if temp == 1:
        N = -1                      #removes central cell
    else:
        N = 0  
    
    for i in range(3):      #range=(0,1,2)
            for j in range(3):
                #considering if periodic boundaries are used or not
                xi = (posX + i - 1) #from (0,1,2) to (-1,0,1)
                yj = (posY + j - 1)
                if Closed_boundaries == True:
                    if xi < 0 or xi >= NbC or yj < 0 or yj >= NbL:
                        continue
                else:
                    xi = xi % NbC
                    yj = yj % NbL
                #count neighbors
                if StateMatrix[xi][yj] == 1:
                    N = N + 1
    return N

def VN(posX, posY,StateMatrix,NbC,NbL,Closed_boundaries):      #called by rules(), counts the live neighbors to each cell and gives it as output
    N=0
 
    for i in range(3):      #modified moore method to von neumann neighbourhood
            for j in range(3):
                #considering if periodic boundaries are used or not
                xi = (posX + i - 1)
                yj = (posY + j - 1)
                if Closed_boundaries == True:
                    if xi < 0 or xi >= NbC or yj < 0 or yj >= NbL:
                        continue
                else:
                    xi = xi % NbC
                    yj = yj % NbL
                #count neighbors
                if xi==yj or abs(xi-yj)==2: #filters out the extra neighbours 
                    continue
                if StateMatrix[xi][yj] == 1:
                    N = N + 1
    return N
def Moore9(posX, posY,StateMatrix,NbC,NbL,Closed_boundaries):      #called by rules(), counts the live neighbors to each cell and gives it as output
    N=0
    
    for i in range(3):      #same as moore without the removal of the central cell
            for j in range(3):
                #considering if periodic boundaries are used or not
                xi = (posX + i - 1)
                yj = (posY + j - 1)
                if Closed_boundaries == True:
                    if xi < 0 or xi >= NbC or yj < 0 or yj >= NbL:
                        continue
                else:
                    xi = xi % NbC
                    yj = yj % NbL
                #count neighbors
                if StateMatrix[xi][yj] == 1:
                    N = N + 1
    return N    

# Animation pause
def pause():        #runs after press of Pause button, sets flag to 0. i.e. stops iterate()
    global flag
    flag=0


# Animation start
def start():        #runs after press of Start button, sets flag to 1. i.e. starts iterate()
    global flag
    if flag==0: 
        flag=1
    iterate()


# Animation step by step
def stepbystep():   #runs after press of Step button, sets flag=2, i.e. makes iterate() put flag=0 after one iteration, thus stop
    global flag
    flag=2
    iterate()
        
    
#Stop animation and do plot
def plotit():       #runs upon press of the Plot button and is last for activities where same. Gives "temp."-plot and activity-gen plot    
    global row
    row=row+2
    g = list(range(1, len(ActivityList)+1) )        #generation list
    plt.plot(g, ActivityList)
    plt.title("temperature  %i" % row )
    plt.xlabel("Generation")
    plt.ylabel("Activity")
    plt.grid()
    plt.show()
    #root.destroy()
    
    g = list(range(1, len(EnergyList)+1) )     
    plt.plot(g, EnergyList)
    plt.title("Total Energy %i" % row)
    plt.xlabel("Generation")
    plt.ylabel("Energy")
    plt.grid()
    plt.show()
    
    g = list(range(1, len(S)+1) )     
    plt.plot(g, S)
    plt.plot([g[0],g[-1]],[np.log(2),np.log(2)]) #log2 limit?
    plt.title("Gibbs Entropy %i" % row)
    plt.xlabel("Generation")
    plt.ylabel("Entropy")
    plt.grid()
    plt.show()
    
    g = list(range(1, len(T)+1) )     
    plt.plot(g, T)
    plt.title("TD Temperature %i" % row)
    plt.xlabel("Generation")
    plt.ylabel("Energy/Entropy")
    plt.grid()
    plt.show()
    
    g = list(range(1, len(PressureList)+1) )     
    plt.plot(g, PressureList)
    plt.title("Pressure - Collisions at the boundaries %i" % row)
    plt.xlabel("Generation")
    plt.ylabel("Pressure")
    plt.grid()
    plt.show()
    
    #local activity over generation
    if see_local_s == True:
        plt.plot(g, localActList)
        plt.title("local temperature %i" % row)
        plt.xlabel("Generation")
        plt.ylabel("Activity in local area")
        plt.grid()
        plt.show()
    
# Mouse click functions
def  CellState(event):              #runs by left-click. Changes current state of a cell and colors it
    x, y = event.x//a, event.y//a
    if StateMatrix[x,y]==state["Dead"]:
        StateMatrix[x,y]=state["Alive"]
        color = "black"
    else:
        StateMatrix[x,y]=state["Dead"]
        color = "white"
    canvas.itemconfig(SM1[x][y], fill=color)
    
def  FillCellState(event):          #Runs by left motion-click. Make cells alive and colors them black
    x, y = event.x//a, event.y//a
    StateMatrix[x,y]=state["Alive"]
    color = "black"
    canvas.itemconfig(SM1[x][y], fill=color)

def  KillCells(event):              #Runs by right motion-click. Make cells dead and colors them white
    x, y = event.x//a, event.y//a
    StateMatrix[x,y]=state["Dead"]
    color = "white"
    canvas.itemconfig(SM1[x][y], fill=color)

# Grpahic Interface
root.title("Cellular automata")
canvas = Canvas(root, width=a*NbC+modify, height=a*NbL+modify, highlightthickness=0)
root.wm_attributes("-topmost", True) 

canvas.bind("<Button-1>", CellState)      #Click on a cell to change its state
canvas.bind("<B1-Motion>", FillCellState) #Left-click and drag to "draw"
canvas.bind("<B3-Motion>", KillCells)     #Right-click and drag to "erase"
canvas.pack()
                        
bou1 = Button(root, text='Start '+ u"\u25B6", fg="green", width=8, command=start) #start button
bou1.pack(side=LEFT)
bou2 = Button(root, text='Pause '+ u"\u23F8", width=8, command=pause)             #Pause button
bou2.pack(side=LEFT)
bou3 = Button(root, text='Step ' + u"\u23ED", width=8, command=stepbystep)        #one step button
bou3.pack(side=LEFT)
bou4 = Button(root, text='Reset', width=8, command=initialize_map)                #Reset/clear screen button
bou4.pack(side=LEFT)
bou5 = Button(root, text='Plot', width=8, command=plotit)                         #Plot button
bou5.pack(side=LEFT)

#packs the labels 
counter.pack(side=TOP)
activity.pack(side=TOP)
energy.pack(side=TOP)
entropy.pack(side=TOP)
temperature.pack(side=TOP)
pressure.pack(side=TOP)
if see_local_s == True:
    localActivity.pack(side=TOP)

# Launching automata
initialize_map()        
iterate()
root.mainloop()


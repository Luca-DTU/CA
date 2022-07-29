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
import matplotlib.style
import matplotlib as mpl
mpl.style.use('default')
from functions import Moore,draw,initialize_map,initialize,rulegenerator,differentiate                 #,plotit
from scipy import signal
import seaborn as sns
kernel = np.ones((3, 3), dtype=np.int8)
kernel[1, 1] = 0
x=1             #update rate of the program in ms
NbL = 400 #size of sides
NbL=int(NbL)                                        
NbC = NbL              #square map
part=5 #partitions of each side
Closed_boundaries = True #boundaries
birth,death=[1,3,6,7],[0,2,4,5] #random rules
print('rule=',birth,death)
density =  0.1                  #initial density 

if Closed_boundaries == False or Closed_boundaries== True:  #input check for boundary
    modify = 1
else:
    raise Exception('"Closed_boundaries" should take argument True or False. Given argument for "Closed_boundaries" was: {}'.format(Closed_boundaries))


#Creation of matrices and initialization of variables
StateMatrix = np.zeros((NbL,NbC),dtype=int) #current state init
SM1 = np.zeros((NbL,NbC),dtype=int) #visualization matrix
flag=0 #flag defines the status of the automaton: 0=inactive, 1=active
gen=0 #generation counter
run = 1 #to take resets into account
N=NbL*NbC#total number of cells
size=int(NbL/part) #lenght of boxes
nn=size*size  #size ofboxes
H=[] #entropy
E=[]#energy
P=[]#pressure
T=[]#temperature

#initialize matrices for each quantities
h=np.zeros([part,part])
e=np.zeros([part,part])
t=np.zeros([part,part])
p=np.zeros([part,part])

root = Tk() #basic widget
counter = Label(root)   
counter.config(text="Gen " + str(gen) + " Run: " + str(run)) #visual counter config  

# Calculating next generation
def iterate():     
    '''
    This is the beating heart of the code, sets up the cellualar automata and
    determines if it should be running or not, visualized or not based on the value of 
    MAP and flag. It calls rules() which implements the generation changes. 
    Cannot be exported into module.
    '''
    global flag     #global lets you call a global variable as local in a function
    global gen      
    global run
    counter.config(text="Gen " + str(gen) + " Run: " + str(run))
                
    if flag==1: #active
        gen += 1
        rules()
        root.after(x, iterate)     #waits x ms, then runs itearate() again
    elif flag==2: #one step
        gen += 1
        rules()
        flag=0                       #set to inactive 
        root.after(x, iterate)
    elif flag==3: #reset
        run+=1
        if MAP ==True:#initialize reads the density and creates a random distribution accordingly
            initialize_map(density,NbL,NbC,StateMatrix,canvas,a,modify)
        else:
            initialize(density,NbL,NbC,StateMatrix)
        rules()
        flag=1
        root.after(x, iterate)
    else:
        flag=0 #if no action, inactive
# implements rules of propagation, calculates the value of thermodynamical variables 
# at the current time-step
def rules(NbL=NbL,NbC=NbC,birth=birth,death=death):    
    '''
    Reads the neighbourhood of every cell and updates the values in the CA matrix
    accordingly. Computes the values of the thermodynamical variables and adds
    them to their relative arrays.
    Cannot be exported into modules 
    '''    
    global StateMatrix,H,E,P,T
    
    NEWmatrix = np.zeros((NbL,NbC),dtype=int)
    #update ruels
    k=signal.convolve(StateMatrix, kernel, mode='same')
    for i in range(NbL):
            for j in range(NbC):
                n = k[i,j] 
                if n in birth:                
                    NEWmatrix[i][j] = 1                    #birth
                elif n in death:         
                    NEWmatrix[i][j] = 0                    #death
                else:                               
                    NEWmatrix[i][j] = StateMatrix[i][j]    #unaffected

    SM1 = StateMatrix.copy() 
    StateMatrix = NEWmatrix.copy()          #updates the old matrix to the new
    
    h=np.zeros([part,part])
    e=np.zeros([part,part])
    t=np.zeros([part,part])
    p=np.zeros([part,part])
    #iterates the variable calculation for each box in the prtition
    for i in range (part): #goes through the partition
        for j in range(part):
            Box=(StateMatrix[size*(i):size*(i+1),  size*(j):size*(j+1)]) #takes the boxes
            e[i][j]=np.sum(Box)/(nn)
            h[i][j]=(-((e[i][j])*np.log((e[i][j]))+(1-(e[i][j]))*np.log(1-(e[i][j]))))
            for k in range(size):
                if Box[k,0]==1:
                    p[i][j]=p[i][j]+1
                if Box[k,size-1]==1:
                    p[i][j]=p[i][j]+1
                if Box[0,k]==1:
                    p[i][j]=p[i][j]+1
                if Box[size-1,k]==1:
                    p[i][j]=p[i][j]+1
    p=p/(size*4)
    H.append(h)
    E.append(e)
    P.append(p)
    #temperature
    if len(E)>3: 
        for i in range (part):
            for j in range(part):
                t[i][j]=(E[-1][i,j]-E[-3][i,j])/(H[-1][i,j]-H[-3][i,j])
            
    T.append(t)
    


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
            
#plot
def plotit():#plot of the lists  

    '''
    Plot of each the values of each box in the partition in both a seaborn heatmap and a
    matplotlib plot.
    The heatmap represents the average values over time of the boxes and fundamentally 
    helps with the identification of nonequilibrium stationary states
    '''
    q=str(birth) +'-'+ str(death)
    mpl.rcParams['font.size'] = 27
    mpl.rcParams['font.family'] ='serif'

    energy=np.zeros([part,part])
    entropy=np.zeros([part,part])
    pressure=np.zeros([part,part])
    temperature=np.zeros([part,part])
    for j in range(part):
        for k in range (part):       
            a1=[]
            for i in range(len(T)):
                a1.append(T[i][j,k])
            plt.plot(a1)
            temperature[j][k]=np.nanmean(a1)
    plt.title("Temperature of boxes " + q)
    plt.grid()
    plt.show()
    
    
    for j in range(part):
        for k in range (part):       
            a1=[]
            for i in range(len(P)):
                a1.append(P[i][j,k])
            plt.plot(a1)
            pressure[j][k]=np.mean(a1)
    plt.title("Pressure of boxes " + q)
    plt.grid()
    plt.show()

    for j in range(part):
        for k in range (part):       
            a1=[]
            for i in range(len(H)):
                a1.append(H[i][j,k])
            plt.plot(a1)
            entropy[j][k]=np.mean(a1)
    plt.title("Entropy density of boxes " + q)
    plt.grid()
    plt.show()

    
    for j in range(part):
        for k in range (part):       
            a1=[]
            for i in range(len(E)):
                a1.append(E[i][j,k])
            plt.plot(a1)
            energy[j][k]=np.mean(a1)
    plt.title("Energy density of bozes " + q)
    plt.grid()
    plt.show()

    plt.subplot() # first heatmap
    plt.title('Energy')
    sns.heatmap(energy,annot=True) 
    plt.show()     
    plt.subplot() # second heatmap
    plt.title('Pressure')
    sns.heatmap(pressure,annot=True) 
    plt.show()     
    plt.subplot() # third heatmap
    plt.title('Entropy')
    sns.heatmap(entropy,annot=True) 
    plt.show()     
    plt.subplot()
    # plt.title('Temperature')
    sns.heatmap(temperature,annot=True) 
    plt.show()   
    
# Mouse click functions
def  CellState(event):              #runs by left-click. Changes current state of a cell and colors it
    x, y = event.x//a, event.y//a
    if StateMatrix[x,y]==0:
        StateMatrix[x,y]=1
        color = "black"
    else:
        StateMatrix[x,y]=0
        color = "white"
    canvas.itemconfig(SM1[x][y], fill=color)
    draw(NbL,NbC,StateMatrix,a,modify,canvas)
    
def  FillCellState(event):          #Runs by left motion-click. Make cells alive and colors them black
    x, y = event.x//a, event.y//a
    StateMatrix[x,y]=1
    color = "black"
    canvas.itemconfig(SM1[x][y], fill=color)
    draw(NbL,NbC,StateMatrix,a,modify,canvas)

def  KillCells(event):              #Runs by right motion-click. Make cells dead and colors them white
    x, y = event.x//a, event.y//a
    StateMatrix[x,y]=0
    color = "white"
    canvas.itemconfig(SM1[x][y], fill=color)
    draw(NbL,NbC,StateMatrix,a,modify,canvas)

# Grpahic Interface
root.title("Cellular automata")
root.wm_attributes("-topmost", True) 

                        
bou1 = Button(root, text='Start '+ u"\u25B6", fg="green", width=8, command=start) #start button
bou1.pack(side=LEFT)
bou2 = Button(root, text='Pause '+ u"\u23F8", width=8, command=pause)             #Pause button
bou2.pack(side=LEFT)
bou3 = Button(root, text='Step ' + u"\u23ED", width=8, command=stepbystep)        #one step button
bou3.pack(side=LEFT)
bou4 = Button(root, text='Plot', width=8, command=plotit)                         #Plot button
bou4.pack(side=LEFT)

counter.pack(side=TOP)

# Launching automata

initialize(density,NbL,NbC,StateMatrix)
iterate()
root.mainloop()


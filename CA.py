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
import pandas as pd
import ast
from functions import Moore,draw,initialize_map,initialize ,rulegenerator,differentiate            #,plotit
from scipy.optimize import curve_fit
import random
'''    #from file
x=1             #update rate of the program in ms
MAP=False       #Boolean for the visualization of the cellular automaton grid
row = 5         #which conditions to read if reading xcel file
#---Imput row to be read in Excel file: "Variables and initial configuration"---#                                                                               #
Variables = pd.read_excel('Variables and initial configuration.xlsx', sheet_name='Blad1')  
row = row - 2           #for one-to-one correspondence with excel file
NbL=Variables.iloc[row,0]                          #length of sides                            
NbL=int(NbL)                                        
NbC = NbL                                           #square map
Closed_boundaries=Variables.iloc[row,2]           #closes the boundaries if True, else they are periodic 
a=Variables.iloc[row,1]  
birth = ast.literal_eval(Variables.iloc[row,3])          #birth rule from file
death = ast.literal_eval(Variables.iloc[row,4])          #death rule from file     
density = Variables.iloc[row,5]                                   #initial density 
'''
x=1             #update rate of the program in ms
MAP=False       #Boolean for the visualization of the cellular automaton grid                                                                         #
NbL = 400 #lenth of sides, note that complexity is O(NbL^2)                     
NbL=int(NbL)                                        
NbC = NbL                                           #square map
if MAP==True: #size of the cells depending on MAP
    a=10
else:
    a=1
Closed_boundaries = True #closes the boundaries if True, else they are periodic 
birth,death=rulegenerator() #random rules
density =  0.1         #initial density 
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
s = 0 #current activity
P=0   #pressure  init
E=density*N #initial energy init

h=N*np.log(E)#ideal gas entropy
hg=(-N*((E/N)*np.log((E/N))+(1-(E/N))*np.log(1-(E/N))))#gibbs entropy
#arrays for activity, energy, entropy, temperature and pressure for plots and calculations
ActivityList = np.array([])
T=np.array([]) #from ideal gas entropy
Tg=np.array([]) #from gibbs entropy
PressureList=np.array([])
EnergyList=np.array([E])    
S=np.array([h]) #ideal
Sg=np.array([hg]) #gibbs
helmholtz=np.array([]) #from Z
Helmholtz=np.array([]) #as defined
IH=np.array([])
C=np.array([]) #ideal gas heat capacity
Cg=np.array([]) #regular (gibbs) heat capacity
loglist=np.array([]) #log(z)
#set up of graphics
root = Tk()                 #name change of widget tool
counter = Label(root)       #label for quantities necessary for interface
activity = Label(root)      
energy=Label(root)
entropy=Label(root)
temperature= Label(root)
pressure=Label(root)
counter.config(text="Gen " + str(gen) + " Run: " + str(run)) #visual counter config
activity.config(text="Average Activity " + str(s))
energy.config(text="Energy "+str(E))
entropy.config(text="Entropy "+str(h))
pressure.config(text="Pressure "+str(P))

# Calculating next generation
def iterate():  
    '''    This is the beating heart of the code, sets up the cellualar automata and
    determines if it should be running or not, visualized or not based on the value of 
    MAP and flag. It calls rules() which implements the generation changes. 
    Cannot be exported into module.
    '''
    global flag     #global lets you call a global variable as local in a function
    global gen      
    global run
    global birth,death
    if MAP ==True:
        draw(NbL,NbC,StateMatrix,a,modify,canvas) #draws the update
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
        birth,death=rulegenerator() #generates new rules
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


def rules(NbL=NbL,NbC=NbC,birth=birth,death=death):  
    '''
    Reads the neighbourhood of every cell and updates the values in the CA matrix
    accordingly. Computes the values of the thermodynamical variables and adds
    them to their relative arrays.
    Cannot be exported into modules 
    '''    
    #initialize counter based variables
    s=0
    E=0
    P=0
    #lists and visual counters need to be global
    global StateMatrix,EnergyList,energy,S,Sg,entropy,T,Tg,temperature,activity,ActivityList,PressureList,pressure,helmholtz,IH,Helmholtz,C,Cg,zE,zvarE,zC,loglist
    
    NEWmatrix = np.zeros((NbL,NbC),dtype=int)
    #update ruels
    for i in range(NbL):
            for j in range(NbC):
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
                if SM1[i][j] != StateMatrix[i][j]: #counts the change of state of each cell
                    s+=1
    #energy counter
    E=np.sum(StateMatrix)
    
    for i in range (NbL): #counts the number of active cells at the boundaries#in analogy to the microscopic definition of pressure
        if StateMatrix[i,0]==1:
            P=P+1
        if StateMatrix[i,NbL-1]==1:
            P=P+1
        if StateMatrix[0,i]==1:
            P=P+1
        if StateMatrix[NbL-1,i]==1:
            P=P+1
    P=P/(4*NbL) #pressure is intensive
    

    h=N*np.log(E/N)+N*np.log(N) #ideal gas entropy

    hg=(-N*((E/N)*np.log((E/N))+(1-(E/N))*np.log(1-(E/N)))) #gibbs entropy
    #addition to arrays and visuals 
    PressureList=np.append(PressureList,P)
    pressure.config(text="Pressure "+str(P))
    EnergyList=np.append(EnergyList,E)
    energy.config(text="Energy " + str(E))             
    S=np.append(S,h)
    Sg=np.append(Sg,hg)
    entropy.config(text="Entropy "+str(h))
    ActivityList=np.append(ActivityList,s/N)                         #average activity list
    activity.config(text="Average Activity " + str(s/N))  
    if len(EnergyList)>3:  #to avoid errors in the first steps in the calculation of temp
        t=differentiate(S,EnergyList,4)
        T=np.append(T,t)
        temperature.config(text="Temperatue "+str(t))  
        tg=differentiate(Sg,EnergyList,4)
        Tg=np.append(Tg,tg)
        Helmholtz=np.append(Helmholtz,E-tg*hg)      
        IH=np.append(IH,E-t*h)
        Es=np.arange(N+1)
        Z=np.sum(np.exp(-Es/(t))) #partition function
        z=np.log(Z)
        loglist=np.append(loglist,z)
        A=-z*t#helmholtz free energy
        helmholtz=np.append(helmholtz,A)
    if len(T)>4:        #to avoid errors
        c=differentiate(T,EnergyList,5)
        C=np.append(C,c)
        cg=differentiate(Tg,EnergyList,5)
        Cg=np.append(Cg,cg)
        
    

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
    
def reset(): 
    #resets lists, updates run, resets flag and gen.
    global ActivityList,T,PressureList,EnergyList,S,birth,death,Sg,Tg,Cg
    ActivityList = np.array([])
    T=np.array([])
    Tg=np.array([])
    Cg=np.array([])
    PressureList=np.array([])  
    EnergyList=np.array([E])   
    S=np.array([h])
    Sg=np.array([hg])
    global gen,run,flag
    gen=0
    flag=3
      
#plot 
def plotit():#plot of the lists  
    
    
    k=str(birth) +'-'+ str(death)
    mpl.rcParams['font.size'] = 18
    #TEMPERATURE
    
    g=np.arange(1, len(ActivityList)+1) #generation list           
    plt.plot(g, ActivityList,label='activity temperature')
    g = np.arange(2, len(T)+2)  
    plt.plot(g, T,label='temperature')
    plt.title("temperature, rule " +k )
    plt.xlabel("Generation")
    plt.ylabel("Temperature")
    plt.legend(loc="best")
    plt.grid()
    #plt.xticks([])
    #plt.yticks([0,max(ActivityList),max(T)])
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Temp.png')
    plt.show()
    
    g = np.arange(2, len(Tg)+2)  
    plt.plot(g, Tg,label='temperature')
    plt.title("temperature, rule " +k )
    plt.xlabel("Generation")
    plt.ylabel("Temperature")
    plt.legend(loc="best")
    plt.grid()
    #plt.xticks([])
    #plt.yticks([0,max(T)])
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Temp.png')
    plt.show()

    #ENERGY
    g = np.arange(1, len(EnergyList)+1)      
    plt.plot(g, EnergyList, label='energy')
    Activity=ActivityList*N
    IE=N*T
    g = np.arange(1, len(Activity)+1) 
    plt.plot(g,Activity,label='Activity') 
    g = np.arange(1, len(IE)+1) 
    plt.plot(g,IE,label='E=NT') 
    plt.title("Total Energy, rule" +k)
    plt.xlabel("Generation")
    plt.ylabel("Energy")
    plt.legend(loc="best")
    plt.grid()
    #plt.xticks([])
    #plt.yticks([0,N/2],('0','N/2'))
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Energy.png')
    plt.savefig('Ideal_gas/'+k+'_Energy.png')
    plt.show()
    
    Activity=np.append(Activity,Activity[-1])
    c=EnergyList/Activity
    g = np.arange(1, len(c)+1)  
    plt.plot(g, c) 
    plt.title("energy metric ratio, rule " +k )
    plt.xlabel("Generation")
    plt.ylabel("E / E_A")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Energy ratio.png')
    plt.show()
    
    g = np.arange(1, len(Sg)+1) 
    plt.plot(g,Sg,label='Gibbs Entropy')
    plt.title("Entropy, rule" +k)
    plt.xlabel("Generation")
    plt.ylabel("Entropy")
    plt.grid()
    plt.legend(loc="best")
    #plt.xticks([])
    #plt.yticks([0,N*np.log(2)],('0','N log(2)'))
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Entropy.png')
    plt.show()
    
    g = np.arange(1, len(S)+1) 
    plt.plot(g,S,label='Ideal Entropy')
    plt.title("Entropy, rule" +k)
    plt.xlabel("Generation")
    plt.ylabel("Entropy")
    plt.grid()
    plt.legend(loc="best")
    #plt.xticks([])
    #plt.yticks([])
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Entropy.png')
    plt.show()
        
    plt.scatter(Tg,Sg[3:])
    plt.title('Gibbs Entropy vs T, RULE '+k)
    plt.xlabel('T')
    plt.ylabel('S')
    plt.grid()
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Entropy_VS_T.png')        
    plt.show()

    plt.scatter(T,S[3:])
    plt.title('Entropy vs T, RULE '+k)
    plt.xlabel('T')
    plt.ylabel('S')
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Entropy_VS_T.png')        
    plt.show()
    
    g = np.arange(1, len(PressureList)+1)    
    plt.plot(g, PressureList)
    plt.title("Pressure, rule" + k)
    plt.xlabel("Generation")
    plt.ylabel("Pressure")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Pressure.png')
    plt.show()
    
    g = np.arange(1, len(helmholtz)+1)        #generation list
    plt.plot(g, helmholtz, label='From partition function') #free energy from partition function
    plt.legend(loc="best")
    plt.title("Helmholtz free energy,rule" +k )
    plt.xlabel("Generation")
    plt.ylabel("Helmholtz")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Partition_function/'+k+'_Helmholtz.png')
    plt.show()
    
    g = np.arange(1, len(Helmholtz)+1)        #generation list
    plt.plot(g, Helmholtz, label='E - TS') #definition of free energy
    plt.legend(loc="best")
    plt.title("Helmholtz free energy,rule" +k )
    plt.xlabel("Generation")
    plt.ylabel("Helmholtz")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Helmholtz.png')
    plt.show()
    
    g = np.arange(1, len(IH)+1)        #generation list
    plt.plot(g, Helmholtz, label='E - TS_A') #definition of free energy
    plt.legend(loc="best")
    plt.title("Helmholtz free energy,rule" +k )
    plt.xlabel("Generation")
    plt.ylabel("Helmholtz")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Helmholtz.png')
    plt.show()
        
    g = np.arange(1, len(C)+1)         #generation list
    plt.plot(g, C) #Heat capacity
    plt.plot([0,len(g)],[N,N],label='C = N')
    plt.legend(loc="best")
    plt.title("Heat capacity, rule" +k )
    plt.xlabel("Generation")
    plt.ylabel("C")
    #plt.xticks([])
    #plt.yticks([0,N],('0','N'))
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Heat capacity.png')
    plt.show()
    
    g = np.arange(1, len(Cg)+1)         #generation list
    plt.plot(g, Cg) #Heat capacity
    plt.plot([0,len(g)],[N,N],label='C = N')
    plt.legend(loc="best")
    plt.title("Heat capacity, rule" +k )
    plt.xlabel("Generation")
    plt.ylabel("C")
   # plt.xticks([])
    plt.grid()
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Heat capacity.png')
    plt.show()
    
    
    #partition function temperature input 
    #Q=ActivityList[2:] 
    Q=T
    zE=-np.gradient(loglist,1/Q)
    zvarE=np.gradient(-zE,1/Q)
    zC=zvarE/Q**2

    
    plt.scatter(Tg,EnergyList[3:])
    plt.title('Energy vs Temperature, rule' +k)
    plt.xlabel('temperature')
    plt.ylabel('Energy')
    plt.grid()
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Energy vs temperature.png')
    plt.show()
    
    plt.scatter(T,EnergyList[3:])
    plt.title('Energy vs Temperature, rule' +k)
    plt.xlabel('temperature')
    plt.ylabel('Energy')
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Energy vs temperature.png')
    plt.show()
    
    plt.scatter(Q, zE)
    plt.title("Expected energy from Z, rule" +k)
    plt.xlabel("Temperature")
    plt.ylabel("E[E]")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Partition_function/'+k+'_Energy.png')
    plt.show()
    
    
    plt.scatter(Q,zC)
    plt.title("Heat Capacity from Z, rule" +k)
    plt.xlabel("Temperature")
    plt.ylabel("Heat capacity")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Partition_function/'+k+'_Heat capacity.png')
    plt.show()
    
    
    zentropy= zE/Q+loglist
    plt.scatter(Q, zentropy, label='partition function entropy')
    plt.title("Entropy, rule " + k)
    plt.xlabel("Temperature")
    plt.ylabel("Entropy")
    plt.grid()
    #plt.xticks([])
    #plt.yticks([])
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig('Partition_function/'+k+'_Entropy.png')
    plt.show()
    
    
#Z Quantities as a function of time
    g = np.arange(1, len(zE)+1)   
    plt.plot(g, zE)
    plt.title("Expected energy from Z, rule" +k)
    plt.xlabel("generation")
    plt.ylabel("E[E]")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Partition_function/'+k+'_Energy-g.png')
    plt.show()
    
    g = np.arange(1, len(zC)+1)   
    plt.plot(g,zC)
    plt.title("Heat Capacity from Z, rule" +k)
    plt.xlabel("generation")
    plt.ylabel("Heat capacity")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Partition_function/'+k+'_Heat capacity-g.png')
    plt.show()
    
    
    zentropy= zE/Q+loglist
    g = np.arange(1, len(zentropy)+1)   
    plt.plot(g, zentropy, label='partition function entropy')
    plt.title("Entropy, rule " + k)
    plt.xlabel("generation")
    plt.ylabel("Entropy")
    plt.grid()
    #plt.xticks([])
    #plt.yticks([])
    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig('Partition_function/'+k+'_Entropy-g.png')
    plt.show()
    
#--    
    

   
    c=T/Activity[3:]
    g = np.arange(1, len(c)+1)         #generation list
    plt.plot(g, c) 
    plt.title("Temperature ratio, rule"+k )
    plt.xlabel("Generation")
    plt.ylabel("T / T_A")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Temp ration.png')
    plt.show()
    
    c=Tg/Activity[3:]
    g = np.arange(1, len(c)+1)         #generation list
    plt.plot(g, c) 
    plt.title("Temperature ratio, rule"+k )
    plt.xlabel("Generation")
    plt.ylabel("T_g / T_A")
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Gibbs_Temp ration.png')
    plt.show()
    
    
    
    zP=1/(np.exp(N/Q)-1)
    plt.scatter(Q,zP)
    plt.title('Partition function pressure'+k)
    plt.xlabel("Temperature")
    plt.ylabel("Pressure")
    plt.grid()
    #plt.xticks([])
    #plt.yticks([])
    plt.tight_layout()
    plt.savefig('Partition_function/'+k+'_pressure.png')
    plt.show()
    
    
    
    plt.scatter(Tg,PressureList[2:])
    plt.title("Pressure, rule"+ k)
    plt.xlabel("Temperature")
    plt.ylabel("Pressure")
    #plt.xticks([])
    #plt.yticks([])
    plt.grid()
    plt.tight_layout()
    plt.savefig('Standard/'+k+'_Pressure vs T.png')
    plt.show()
    
    plt.scatter(T,PressureList[2:])
    plt.title("Pressure, rule"+ k)
    plt.xlabel("Temperature")
    plt.ylabel("Pressure")
    #plt.xticks([])
    #plt.yticks([])
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Pressure vs T.png')
    plt.show()

    plt.scatter(T,T*np.log(N*T))
    plt.title("Pressure, rule"+ k)
    plt.xlabel("Temperature")
    plt.ylabel("Pressure")
    #plt.xticks([])
    #plt.yticks([])
    plt.grid()
    plt.tight_layout()
    plt.savefig('Ideal_gas/'+k+'_Pressure vs T model.png')
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
if MAP==True:    
    canvas = Canvas(Tk(), width=a*NbC+modify, height=a*NbL+modify, highlightthickness=0)
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
bou4 = Button(root, text='Reset', width=8, command=reset)               #Reset/clear screen button
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

# Launching automata
if MAP==True:
    initialize_map(density,NbL,NbC,StateMatrix,canvas,a,modify)   
else:
    initialize(density,NbL,NbC,StateMatrix)
iterate()
root.mainloop()


# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 10:18:52 2020

@author: lucab
"""
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
from functions import draw,initialize_map,initialize ,rulegenerator 
from scipy import signal

#Neighbourhood choice for update rules
kernel = np.ones((3, 3), dtype=np.int8)#Neighourhood
kernel[1, 1] = 0#Reduced moore neighbourhood
x=1             #update rate of the program in ms
MAP=False      #Boolean for the visualization of the cellular automaton grid                                                                         #
NbL = 500     #lenth of sides, note that complexity is O(NbL^2)                     
NbL=int(NbL)                                        
NbC = NbL       #square map
if MAP==True:   #size of the cells depending on MAP
    a=10
else:
    a=1
Closed_boundaries = True #closes the boundaries if True, else they are periodic 
birth,death=[3],[0,1,4,5,6,7,8] #random rules, see functions module
#birth,death = rulegenerator()
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


#arrays for activity, energy, entropy, temperature and pressure for plots and calculations
ActivityList = np.array([])
PressureList=np.array([])
EnergyList=np.array([])    

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
pressure.config(text="Pressure "+str(P))

# Calculating next generation
def iterate():  
    '''Sets up the cellualar automata and
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
        if MAP ==True: #initialize reads the density and creates a random distribution accordingly
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
    NbL,NbC integer like map dimensions
    Birth,death list/array like 
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
    global StateMatrix,EnergyList,energy,activity,ActivityList,PressureList,pressure
    
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
                    NEWmatrix[i][j] = StateMatrix[i][j]  #unaffected

    SM1 = StateMatrix.copy() 
    StateMatrix = NEWmatrix.copy()          #updates the old matrix to the new
    
    #activity counter
    for i in range(NbL):
            for j in range(NbC):
                if SM1[i][j] != StateMatrix[i][j]: #counts the change of state of each cell
                    s+=1
    #energy counter
    E=np.sum(StateMatrix)
    #pressure counter
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
    PressureList=np.append(PressureList,P)
    pressure.config(text="Pressure "+str(P))
    EnergyList=np.append(EnergyList,E)
    energy.config(text="Energy " + str(E))  
    ActivityList=np.append(ActivityList,s/N)                         #average activity list
    activity.config(text="Average Activity " + str(s/N))  

        
    

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
    global ActivityList,EnergyList,S,birth,death,PressureList
    ActivityList = np.array([])
    PressureList=np.array([])  
    EnergyList=np.array([])   
    global gen,run,flag
    gen=0
    flag=3

def normalize(x,levels):
    norm_fem = x*np.mean(levels[:5]/x[:5])
    norm_ett = x*levels[0]/x[0]
    return norm_fem
     
#plot 
def plotit():#plot of the lists  
    k=str(birth) +'-'+ str(death)
    mpl.rcParams['font.size'] = 23
    mpl.rcParams['font.family'] ='serif'
    # mpl.rcParams['font.family'] ='Computer Modern Roman'
    # plt.rc('font', family='serif',size=24)
    S=N*np.log(EnergyList)
    Sg=(-N*((EnergyList/N)*np.log((EnergyList/N))+(1-(EnergyList/N))*np.log(1-(EnergyList/N))))
    T=EnergyList/N
    Tg=1/(np.log(1-EnergyList/N)-np.log(EnergyList/N))
    C=np.full((len(EnergyList)), N)
    Cg=-N*T**2/(EnergyList*(EnergyList-N))
    Activity=ActivityList*N
    Helmholtz=EnergyList-Tg*Sg
    IH=EnergyList-T*S
    Q=Tg#Partition fucntion input temperature
    Z=(np.exp(-N/Q)*(np.exp((N+1)/Q)-1))/(np.exp(1/Q)-1)
    
    Z=np.array([])    
    Es=np.arange(N+1)
    for i in range(len(Q)):
        Z=np.append(Z,np.sum(np.exp(-Es/(Q[i])))) #partition function
    
    
    z=np.log(Z)
    helmholtz=-z*Q#helmholtz free energy
    zE=-np.gradient(z,1/Q)
    zvarE=np.gradient(-zE,1/Q)
    zC=zvarE/Q**2
    zS= zE/Q+z
    
    
    pf=np.std(Helmholtz[:40]/helmholtz[:40])/np.mean(Helmholtz[:40]/helmholtz[:40])
    # print('partition function',pf)
    #TEMPERATURE


    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    g = np.arange(1, len(Tg)+1)  
    plt.plot(g, Tg,label='Benchmark')
    plt.plot(g, normalize(T,Tg),label='Ideal gas', linestyle='dashed')
    plt.plot(g, normalize(ActivityList,Tg),label='Avg Activity', linestyle='dotted')
    # plt.title("Temperature")
    plt.xlabel("Generation")
    plt.yticks([])
    plt.ylabel("T [Arbitrary Units]")
    # plt.legend(loc="best")
    # plt.ylabel("T")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    plt.xticks([0,25,50,75,100])
    #plt.grid()
    #plt.xticks([])
    #plt.yticks([0,max(T)])
    plt.tight_layout()
    plt.savefig('LAST_PLOTS/'+k+'_Temp.png')
    plt.show()

    #ENERGY
    g = np.arange(1, len(EnergyList)+1)      
    plt.plot(g, EnergyList, label='energy')
    #plt.plot( [0, 30.5],[N/2, N/2],':',color="blue")
    #plt.plot([30.5,30.5],[0, N/2],':',color="blue")
    #g = np.arange(1, len(Activity)+1) 
    #plt.plot(g,Activity,label='Activity') 
    # plt.title("Energy")
    plt.xlabel("Generation")
    plt.ylabel("E")
    plt.xticks([0,25,50,75,100])
    # plt.legend(loc="best")
    #plt.grid()
    #plt.xticks([])
    plt.yticks([0,N/4,N/2,3*N/4,N],['0','N/4','N/2','3N/4','N'])

    plt.tight_layout()
    plt.savefig('LAST_PLOTS/'+k+'_Energy.png')
    # plt.savefig('Ideal_gas/'+k+'_Energy.png')
    plt.show()

    #ENTROPY
    g = np.arange(1, len(Sg)+1) 
    plt.plot(g,Sg,label='Benchmark')
    plt.plot(g, normalize(S,Sg), label='Ideal gas', linestyle='dashed')
    plt.plot(g, normalize(zS,Sg), label= 'Partition function', linestyle='dotted')
    # plt.title("Entropy")
    plt.xlabel("Generation")
    plt.ylabel("S [Arbitrary Units]")
    # plt.ylabel("S")
    plt.xticks([0,25,50,75,100])
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #plt.yticks([40000,60000,80000,100000],['4e4','6e4','8e4','1e5'])
    #plt.grid()
    # plt.legend(loc="best")
    plt.yticks([])

    plt.tight_layout()
    plt.savefig('LAST_PLOTS/'+k+'_Entropy.png')
    plt.show()
 
    
    # PRESSURE
    PI=T*np.log(N*T)
    zP=Tg**2/(N-EnergyList)*(np.log(Tg)+np.log(1-(N/EnergyList-1)**-N)+1)+((-N/(N-EnergyList)-1/Tg)*(N/EnergyList-1)**-N)/(1-(N/E-1)**-N)
    zp=np.log((N-2*EnergyList)/EnergyList)/((N-EnergyList)*(np.log((N-EnergyList)/EnergyList)**2))-1/((N-2*EnergyList)*np.log((N-EnergyList)/EnergyList))
    g = np.arange(1, len(PressureList)+1)  
    
    
    plt.plot(g,PressureList,label='Benchmark')
    plt.xticks([0,25,50,75,100])
    plt.plot(g,normalize(PI,PressureList),label='Ideal gas', linestyle='dashed')
    plt.plot(g,normalize(zp,PressureList),label='Partition function', linestyle='dotted')
    #plt.plot(g,normalize(zP,PressureList),label='$P_Z^1$')
    # plt.legend(loc='best')
    # plt.title("Pressure")
    plt.xlabel("Generation")
    plt.ylabel("P [Arbitrary Units]")
    plt.yticks([])
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    # plt.ylabel("P")
    #plt.yticks([0,-0.0001,-0.0002,-0.0003],['0','-1e-4','-2e-4','-3e-4'])
    #plt.grid()
    plt.tight_layout()
    plt.savefig('LAST_PLOTS/'+k+'_Pressure.png')
    plt.show()
    
    # #plt.plot(g,zP,label='$P_Z^1$')
    # plt.plot(g,zp,label='Partition function')
    # #plt.plot(g,zp+zP)
    # plt.legend(loc='best')
    # plt.title("Pressure")
    # plt.xlabel("Generation")
    # plt.ylabel("Partition function")
    # plt.yticks([])
    # #plt.yticks([-0.0002,0,0.0002],['-2e-4','0','2e-4',])
    # #plt.grid()
    # plt.tight_layout()
    # plt.savefig('Standard/'+k+'_Pressure.png')
    # plt.show()
  
    #HELMHOLTZ
    g = np.arange(1, len(helmholtz)+1)        #generation list
    
    plt.plot(g, Helmholtz, label='Benchmark')
    plt.plot(g, normalize(IH,Helmholtz), label='Ideal gas', linestyle='dashed')
    plt.plot(g, normalize(helmholtz,Helmholtz), label='Partition function', linestyle='dotted') 
    # plt.legend(loc="best")
    # plt.title("Helmholtz free energy" )
    plt.xlabel("Generation")
    plt.xticks([0,25,50,75,100])
    plt.ylabel("A [Arbitrary Units]")
    plt.yticks([])
    # plt.ylabel("A")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    #plt.grid()
    #plt.yticks([-200000,-400000,-600000,-800000],['-2e5','-4e5','-6e5','-8e5'])

    plt.tight_layout()
    plt.savefig('LAST_PLOTS/'+k+'_Helmholtz.png')
    plt.show()

    #HEAT CAPACITY
    g = np.arange(1, len(zC)+1)         #generation list
    plt.plot(g, Cg, label='Benchmark')
    plt.plot(g, normalize(C,Cg), label='Ideal gas', linestyle='dashed')
    plt.plot(g, normalize(zC,Cg),label='Partition function', linestyle='dotted') #Heat capacity
    # plt.title('Heat capacity' )
    plt.xlabel("Generation")
    plt.xticks([0,25,50,75,100])
    plt.ylabel("C [Arbitrary Units]")
    plt.yticks([])
    # plt.ylabel("C")
    plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))
    # plt.xlim(0,g[-5])
    # plt.legend(loc='best')

    #plt.xticks([])
    #plt.grid()
    plt.tight_layout()
    plt.savefig('LAST_PLOTS/'+k+'_Heat capacity.png')
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
root.title("Cellular automata,"+str(birth) +'-'+ str(death))
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
pressure.pack(side=TOP)

# Launching automata
if MAP==True:
    initialize_map(density,NbL,NbC,StateMatrix,canvas,a,modify)   
else:
    initialize(density,NbL,NbC,StateMatrix)
iterate()
root.mainloop()


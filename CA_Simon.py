"""
Created on Wed Feb 19 15:33:51 2020

@author: Simon
"""

#import packages
import numpy as np
from tkinter import Tk, Canvas, Button, LEFT, TOP, Label
import matplotlib.pyplot as plt
import random
from scipy.ndimage import measurements
import pylab as py
import pandas as pd
import ast
import PlotFitting as PF



#---Imput row to be read in Excel file: "Variables and initial configuration"---#
                                                                                #
row = 4                                                                     #
                                                                                #
D_S = False #True if you want to plot D(S)                                      #
                                                                                #
D_a = True #Should not be True at the same time as D_S!   Plots D(A,r,t)       #
                                                                                #
entropic = True #If True, don't plot while "+"-sign appears after the Run-label#
                                                                                #
#################################################################################




#imports the programmable variables
Variables = pd.read_excel('Variables and initial configuration.xlsx', sheet_name='Blad1')  

#taking account to the way the file is read
row = row - 2 

print(Variables.iloc[row])

#Extracting all variables
NbL = Variables.iloc[row,0]                         #sets dimension of map
NbC = NbL
a = Variables.iloc[row,1]                           #pixel-width of squares in grid
Closed_boundaries = Variables.iloc[row,2]           #closes the boundaries if True, else they are periodic 
birth = ast.literal_eval(Variables.iloc[row,3])     #list of when cells should born 
death = ast.literal_eval(Variables.iloc[row,4])     #list of when cells hould die
density = Variables.iloc[row,5]                     #the initial density of alive cells
see_local_s =  Variables.iloc[row,6]                #True or False, should activity be counted in a specific area
N=NbL*NbC
if see_local_s == True:                             #If True, the coordinates of the rectangular area are needed
    x1, y1, x2, y2 = ast.literal_eval(Variables.iloc[row,7])
    localarea=(x2-x1)*(y2-y1)

if Closed_boundaries == False or Closed_boundaries== True: #why these lines tho?
    modify = 1
else:
    raise Exception('"Closed_boundaries" should take argument True or False. Given argument for "Closed_boundaries" was: {}'.format(Closed_boundaries))


#Creation of matrices
StateMatrix = np.zeros((NbL,NbC),dtype=int)
SM1 = np.zeros((NbL,NbC),dtype=int)

TempArray = np.zeros((NbL,NbC),dtype=int)


#state definition
state = {"Dead":0,"Alive":1}


#variable creation
flag=0
gen=0
s = 0
ClusterCount=0
run = 1
ActivityList = []

sizeList = []
ClusterFrequency=[]

generationMean=[]
GenFrequency = []

Ent=[]
Ent0=[]
first_run = True

k, l = 0, 0                 #D_a element counters
ATSL = []                   #Avtivity to static list
tau = []                    #list of times to static after perturbation

dist_list = []              #List to store distances between the perturbed site and active sites
D_r = False                 #we don't measure any distances initially

root = Tk()                 #name change of widget tool
counter = Label(root)       #label for generation counter
activity = Label(root)
localActivity = Label(root)
clusters = Label(root)

plt.rcParams['axes.unicode_minus'] = False    #To make the plt font-packege work properly 



# Calculating next generation
def iterate():      #checks the value of the flag and directs the program to next step, such as rules()
    global flag
    global gen      #in order to update the generation-counter
    global run
    draw()
    if entropic == True and first_run == False:
        counter.config(text="Gen " + str(gen) + " Run: " + str(run) + "+") #to display that the perturbation has been made
    else:
        counter.config(text="Gen " + str(gen) + " Run: " + str(run))
                
    if flag==1:
        gen += 1
        rules()
        root.after(100, iterate)     #waits 100 ms, then runs itearate() again
    elif flag==2:
        gen += 1
        rules()
        flag=0                       #set to inactive 
        root.after(100, iterate)
    elif flag==3:
        run+=1
        initialize_map()
        rules()
        flag=1
        root.after(100, iterate)
    else:
        flag=0


#Initializing the map
def initialize_map():       #run upon starting the program and used to reset the map
    global gen
    global run
    global flag
    global ActivityList
    global TempArray
    global localActList
    global ClustersList
    flag = 0                #initialize inactive
    
    #counter labels 
    gen = 0
    s = 0
    local_s = 0
    ClusterCount = 0
    ActivityList=[]
    localActList=[]
    ClustersList=[]
    TempArray = np.zeros((NbL,NbC),dtype=int)
    counter.config(text="Gen " + str(gen) + " Run: " + str(run))
    activity.config(text="Activity " + str(s))
    localActivity.config(text="Local activity " + str(local_s))
    clusters.config(text="Clusters " + str(ClusterCount))
    
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



# Defining rules of propagation
def rules():  #runs after iterate() if flag=1,2. Creates updated matrix where each cell gets assigned its state for the next gen
    global gen
    NEWmatrix = np.zeros((NbL,NbC),dtype=int)
    s = 0
    local_s = 0
    
    #update ruels
    for i in range(NbL):
            for j in range(NbC):
                global StateMatrix
                n = neighborCount(i,j);
                if n in birth:                
                    NEWmatrix[i][j] = 1                    #birth
                elif n in death:         
                    NEWmatrix[i][j] = 0                    #death
                else:                               
                    NEWmatrix[i][j] = StateMatrix[i][j]    #unaffected
    
    global SM1c, SM2, SM3 #Not sure what is happening here
    if gen>4:
        SM4 = SM3.copy()
    if gen>3:
        SM3 = SM2.copy()
    if gen>2:
        SM2 = SM1c.copy()
    SM1c = StateMatrix.copy()
    
    SM1 = StateMatrix.copy()
    StateMatrix = NEWmatrix.copy()          #updates the old matrix to the new

       
    #cluster counter
    global CM
    CM = np.zeros((NbL,NbC),dtype=int)      #empty cluster matrix
    label=1
    for y in range(NbL):                    #consider all sites
            for x in range(NbC):
                if StateMatrix[x][y] == 1:  #if site is alive
                    if CM[x][y] == 0:       #check if it has a label
                        CM[x][y] = label    #if it don't have a label, give it one (OR make sure the label is same size )                     
                    for i in range(5):      #look at neighbors to the cell
                        for j in range(5):
                            d = CM[x][y]
                            xi = (x + i - 2)
                            yj = (y + j - 2)
                            if Closed_boundaries == True:
                                if xi < 0 or xi >= NbC or yj < 0 or yj >= NbL:
                                    continue
                            else:
                                xi = xi % NbC
                                yj = yj % NbL
                            if StateMatrix[xi][yj] == 1:       #if neighbor is alive
                                if 0 < CM[xi][yj] < d:         #give it same label as central one
                                    CM = np.where(CM==d, CM[xi][yj], CM)
                                else:
                                    CM[xi][yj] = d
                    label += 1


    #this ensures that the above cluster counter did not miss-labled any cells
    for y in range(NbL):
        for x in range(NbC):
            if CM[x][y] != 0:
                for i in range(5):
                    for j in range(5):
                        xi = (x + i - 2)
                        yj = (y + j - 2)
                        if Closed_boundaries == True:
                            if xi < 0 or xi >= NbC or yj < 0 or yj >= NbL:
                                continue
                        else:
                            xi = xi % NbC
                            yj = yj % NbL
                        if CM[xi][yj] != 0:
                            if CM[xi][yj] != CM[x][y]:
                                CM = np.where(CM == CM[x][y], CM[xi][yj], CM)

    ClusterCount = len(np.unique(CM))-1     #unique number gives number of clusters, -1 since background (0) is not considered a cluster
    ClustersList.append(ClusterCount)       #save cluster history
    clusters.config(text="Clusters " + str(ClusterCount))

     
    #activity counter
    for i in range(NbL):
            for j in range(NbC):
                if SM1[i][j] != StateMatrix[i][j]:
                    s+=1
                    TempArray[j][i] += 1
                    
                    global D_r
                    if D_r == True:
                        global pert                #retrive the coord. of the perturbed site
                        z = np.array((i, j))       #active site
                        r = np.linalg.norm(pert-z) #distance
                        if r > 0:
                            dist_list.append(r)
                    
    if D_a == True: #if we run D_a, ActivityList is needed to be continuously reset, thus a global update is needed
        global ActivityList
    
    ActivityList.append(s/N)                         #activity list
    activity.config(text="Activity " + str(s/N))      
    
    
    #local activity counter
    if see_local_s == True:
        for i in range(x1,x2):
                for j in range(y1,y2):
                    if SM1[i][j] != StateMatrix[i][j]:
                        local_s += 1
    
        localActList.append(local_s/localarea)               #local activity list
        localActivity.config(text="Local activity " + str(local_s/localarea))

    
    #Identify if periodic
    if D_a == True:
        periodic = False
        if gen>4:
            if np.array_equal(StateMatrix, SM2) == True or np.array_equal(StateMatrix, SM3) or np.array_equal(StateMatrix, SM4):
                periodic = True
                s=0
                   
    #condition to reset once it is 2-, 3- or 4-periodic
    if D_a == False:
        global flag
        if gen>4:
            if np.array_equal(StateMatrix, SM2) == True or np.array_equal(StateMatrix, SM3) or np.array_equal(StateMatrix, SM4):
                generationMean.append(gen)
                s=0   #for a CA with periodic patterns being common (like GoL), else replace with "flag = 3" to restart when peridicity is given
    
    #NOT CLEAR YET
    #Entropy of the first stable state (Ent0), with the entropy of stable state after perturbation
    global entropic
    if s==0 and entropic == True:
        global first_run
        uc, ucSize = np.unique(CM, return_counts = True) #unique clusters
        uc=uc[1:]   #remove the 0 background
        ucSize = ucSize[1:]
        instantEntropy = [ucSize**(-3.484)*np.log(ucSize**(-3.484)) for ucSize in ucSize]
        
        if first_run == True:
            Ent0.append(sum(instantEntropy))
            first_run = False
            
            X=random.randrange(NbL)
            Y=random.randrange(NbC)
            while StateMatrix[X][Y] != 0:
                X=random.randrange(NbL)
                Y=random.randrange(NbC)
            StateMatrix[X][Y] = 1
            flag=1
        
        elif first_run == False and sum(instantEntropy) == Ent0[-1] and not np.any(StateMatrix)!=True:
            X=random.randrange(NbL)
            Y=random.randrange(NbC)
            while StateMatrix[X][Y] != 0:
                X=random.randrange(NbL)
                Y=random.randrange(NbC)
            StateMatrix[X][Y] = 1
            flag=1
         
        else:
            Ent.append(sum(instantEntropy))
            flag = 3
            first_run=True
    
    #Bak re-creation
    if gen > 0 and s == 0:                                  #if nothing changed i.e. static cluster(s) 
        #Cluster size frequency, D(S)
        if D_S == True and entropic != True:
            global sizeList        
            for p in [q for q in np.unique(CM) if q != 0]:  #to run for all different clusters
                size = int(np.count_nonzero(CM==p))
                sizeList.append(size)                       #save the size of the cluster
            
            generationMean.append(gen)
            
            flag=3
        
        
        #Active sites after pertubation, D(a)
        if D_a == True:
            global k, l, originalMatrix, originalSM1
            
            if k == 0 and l == 0 and periodic == True: #if the initial config. lead to periodic, we restart
                flag = 3
                
            else:
                if k == 0 and l == 0:   #before we consider the first element in the matrix we save the current matrix
                    originalMatrix = StateMatrix.copy()
                    originalSM1 = SM1.copy()
        
                StateMatrix = originalMatrix.copy()
                SM1 = originalSM1.copy()
                
                
                if not(k == 0 and l == 0):
                    if periodic == False:                  #if it was not periodic save the avtivity after the perturbation
                        if sum(ActivityList) != 0 and sum(ActivityList) != 1:
                            ATSL.append(sum(ActivityList)) #Activity to static list
                        if gen != 0 and gen != 1:
                            tau.append(gen)
                    ActivityList = []
                    gen=0
                        
                if k < NbL-1 or l < NbC-1:
                    if k < NbL-1:
                        if StateMatrix[k][l] == 0:  #if the cell is dead
                            StateMatrix[k][l] = 1   #make it alive and see what happens
                            pert = np.array((k,l))  #perturbed site
                            D_r = True
                        else:
                            D_r = False
                        k += 1
                        flag=1
                   
                    if k == NbL-1 and l < NbC-1:
                        if StateMatrix[k][l] == 0:  #if the cell is dead
                            StateMatrix[k][l] = 1   #make it alive and see what happens
                            pert= np.array((k,l))   #perturbed site
                            D_r = True
                        else:
                            D_r = False
                        k=0
                        l +=1
                        flag=1
                        
                elif k == NbL-1 and l == NbC-1:
                    if StateMatrix[k][l] == 0:      #if the cell is dead
                        StateMatrix[k][l] = 1       #make it alive and see what happens
                        pert = np.array((k,l))      #perturbed site
                        D_r = True
                    else:
                        D_r = False
                    k += 1
                    l += 1
                    flag = 1
                    
                elif k == NbL and l == NbC:    
                    k=0
                    l=0
                    flag = 3
                    D_r = False
                    
                    #Activity to static
                    UA, ActivityFrequency = np.unique(ATSL, return_counts = True) #unique activities
                    ATSD = [ActivityFrequency / len(ATSL) for ActivityFrequency in ActivityFrequency] #Activity To Static Distribution
                    plt.loglog(UA, ATSD)
                    plt.title("Normal distribution of activity to static")
                    plt.xlabel("Activity, A")
                    plt.ylabel("D(A)")
                    plt.show()
                    
                    #generation mean to reach static  
                    UTL, Times = np.unique(tau, return_counts = True)  #Unique Time List
                    TD = [Times / len(tau) for Times in Times]         #Time Distribution
                    plt.loglog(UTL, TD)
                    plt.title("Normal distribution of timesteps until static")
                    plt.xlabel("Timesteps to static, t")
                    plt.ylabel("D(t)")
                    plt.show()
                    
                    #Distance between the perturbed site and active sites
                    UDL, Distances = np.unique(dist_list, return_counts = True)   #Unique Distance List
                    DD = [Distances / len(dist_list) for Distances in Distances]  #Distance Distribution
                    plt.loglog(UDL, DD)
                    plt.title("Normal distribution of distances between perturbation and active sites")
                    plt.xlabel("Distance, r")
                    plt.ylabel("D(r)")
                    plt.show()                    

           

#Neighbor counting
def neighborCount(posX, posY):      #called by rules(), counts the live neighbors to each cell and gives it as output
    temp = StateMatrix[posX][posY]
    if temp == 1:
        N = -1                      #not counting itself
    else:
        N = 0  
    
    for i in range(3):      #looking at the surrounding 3x3  area
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



# Drawing the cells (assigning their colors)
def draw():         #called every time the grid should be colored, assigns a color to each state 
    for x in range(NbL):
        for y in range(NbC):
            if StateMatrix[x,y]==state["Dead"]:
                coul = "white"
            elif StateMatrix[x,y]==state["Alive"]:
                coul = "black"
            canvas.itemconfig(SM1[x,y], fill=coul)  
    
    

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
    g = list(range(1, len(ActivityList)+1) )        #generation list
    plt.plot(g, ActivityList)
    plt.title("Number of births and deaths between each generation")
    plt.xlabel("Generation")
    plt.ylabel("Activity")
    plt.show()
    #root.destroy()
    
    #local activity over generation
    if see_local_s == True:
        plt.plot(g, localActList)
        plt.title("Number of births and deaths in local area")
        plt.xlabel("Generation")
        plt.ylabel("Activity in local area")
        plt.show()
    
    #plot of global TempArray
    global TempArray
    TA = TempArray.copy()
    plt.figure(figsize=(15,6))
    plt.imshow(TA)
    plt.colorbar()
    plt.title("Activity in each site")
    plt.show()
    
    #Plot of cluster areas
    global CM
    CMr=np.rot90(CM)    #rotating to display same orientation as widget
    plt.figure(figsize=(15,6))
    area = measurements.sum(np.rot90(StateMatrix), CMr, index = py.arange(CMr.max() + 1))
    areaPlot = area[CMr]
    areaPlot[areaPlot == 0] = np.nan
    plt.imshow(areaPlot, origin='lower', interpolation='nearest')
    plt.colorbar()
    plt.title("Clusters by area")
    plt.show()
    
    #plot of cluster count
    plt.plot(g, ClustersList)
    plt.title("Number of clusters in each generation")
    plt.xlabel("Generation")
    plt.ylabel("Clusters")
    plt.show()
    
    if D_S == True:
        #cluster size frequency
        global sizeList
        global ClusterFrequency
        if len(sizeList) != 0:
            USL=np.unique(sizeList)
            for p in USL:
                ClusterFrequency.append(np.count_nonzero(sizeList==p))
            CFD = [ClusterFrequency / len(sizeList) for ClusterFrequency in ClusterFrequency] #Cluster Frequency Distribution
            plt.loglog(USL, CFD)
            plt.title("Normal distribution of cluster size frequency")
            plt.xlabel("Cluster Size (S)")
            plt.ylabel("D(S)")
            plt.show()
            ClusterFrequency = []
    
        #generation mean to reach static  
        global generationMean
        UGML=np.unique(generationMean)  #Unique Generation Mean List
        for p in UGML:
            GenFrequency.append(np.count_nonzero(generationMean==p))
        GFD = [GenFrequency / len(generationMean) for GenFrequency in GenFrequency]  #Gen Frequency Distribution
        plt.plot(UGML, GFD)
        plt.title("Normal distribution of generation frequency")
        plt.xlabel("Number of gerenrations to static")
        plt.ylabel("Generation occurrence")
        plt.show()
        
    if entropic == True:
        aEnt, aEnt0 = np.array(Ent), np.array(Ent0)
        diff = np.subtract(aEnt, aEnt0)
        plt.plot(Ent0, diff, ".")
        plt.title("Entropy change after perturbation")
        plt.xlabel("First stable state entropy, H0")
        plt.ylabel("Entropy change, H1-H0")
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
clusters.pack(side=TOP)
if see_local_s == True:
    localActivity.pack(side=TOP)

# Launching automata
initialize_map()        
iterate()
root.mainloop()


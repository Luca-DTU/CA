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
from functions import initialize
import pandas as pd
from scipy import signal
kernel = np.ones((3, 3), dtype=np.int8)
kernel[1, 1] = 0
filename = 'Classes2_CLO_small.pkl'
#uncomment next two lines to reset/create file
# df = pd.DataFrame(columns=['Rule','Dead','T-Ratio_variance','Variation coefficient of A','(E_f-E_i)/std(E)','TT', 'AA'])
# df.to_pickle(filename)

df=pd.read_pickle(filename)  
NbL = 100 #lenth of sides, note that complexity is O(NbL^2)                     
NbL=int(NbL)                                        
NbC = NbL   
N=NbC*NbL                                        #square map
a=1
Closed_boundaries =     True #closes the boundaries if True, else they are periodic 
density =  0.1         #initial density 
if Closed_boundaries == False or Closed_boundaries== True:  #input check for boundary
    modify = 1
else:
    raise Exception('"Closed_boundaries" should take argument True or False. Given argument for "Closed_boundaries" was: {}'.format(Closed_boundaries))

def run(birth,death):
    '''
    

    Parameters
    ----------
    birth : List
    death : List
    
    Runs a single cellular automata and calls compute to compute the variables

    '''
    global StateMatrix,EnergyList,energy,activity,ActivityList

#Creation of matrices and initialization of variables
    StateMatrix = np.zeros((NbL,NbC),dtype=int) #current state init
    gen=0 #generation counter
    s = 0 #current activity
    E=density*N #initial energy init
#arrays for activity, energy, entropy, temperature and pressure for plots and calculations
    ActivityList = np.array([])
    EnergyList=np.array([])   
    initialize(density,NbL,NbC,StateMatrix)
    for i in range (100):
        rules(birth,death)
        gen += 1
        if ActivityList[-1]==0:
            break
    compute(birth,death)


def rules(birth,death,NbL=NbL,NbC=NbC,Closed_boundaries=Closed_boundaries):  
    '''
    Reads the neighbourhood of every cell and updates the values in the CA matrix
    accordingly. Computes the values of the thermodynamical variables and adds
    them to their relative arrays.
    Cannot be exported into modules 
    '''    
    #initialize counter based variables
    s=0
    E=0
  
    #lists and visual counters need to be global
    global StateMatrix,EnergyList,energy,activity,ActivityList
    
    NEWmatrix = np.zeros((NbL,NbC),dtype=int)
    #update ruels

    k=signal.convolve(StateMatrix, kernel, mode='same')
    for i in range(NbL):
            for j in range(NbC):
                n = k[i,j] #can call different neighbourhoods (Moore,Moore9 and VN)
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
    EnergyList=np.append(EnergyList,E)
    ActivityList=np.append(ActivityList,s/N)                         #average activity list


def compute(birth,death):#plot of the lists  
    k=str(birth) +'-'+ str(death)

    #S=N*np.log(EnergyList)
    Sg=(-N*((EnergyList/N)*np.log((EnergyList/N))+(1-(EnergyList/N))*np.log(1-(EnergyList/N))))
    T=EnergyList/N
    Tg=1/(np.log(1-EnergyList/N)-np.log(EnergyList/N))
    #C=np.full((len(EnergyList)), N)
    #Cg=-N*T**2/(EnergyList*(EnergyList-N))
    #Activity=ActivityList*N
    Helmholtz=EnergyList-Tg*Sg
    #IH=EnergyList-T*S
    Es=np.arange(N+1)
    Z=np.array([])
    Q=Tg#Partition fucntion input temperature
    for i in range(len(Q)):
        Z=np.append(Z,np.sum(np.exp(-Es/(Q[i])))) #partition function
    #Z=(np.exp(-N/Q)*(np.exp((N+1)/Q)-1))/(np.exp(1/Q)-1)
    z=np.log(Z)
    helmholtz=-z*Q#helmholtz free energy
    #zE=-np.gradient(z,1/Q)
    #zvarE=np.gradient(-zE,1/Q)
    #zC=zvarE/Q**2
    
    
    #Classifier
    rule=k
    #ideal gas
    if ActivityList[-1]==0:
        dead=1 
    else:
        dead=0
    Ratio=T[:25]/ActivityList[:25]
    #var= np.var(Ratio)
    var = np.std(Ratio)/np.mean(Ratio)
    #partition function
    pf=np.std(Helmholtz[:40]/helmholtz[:40])/np.mean(Helmholtz[:40]/helmholtz[:40])
    
    #NON-Restricted T/ActivityList
    TT = np.std(T/ActivityList)/np.mean(T/ActivityList)
    AA = np.std(Helmholtz/helmholtz)/np.mean(Helmholtz/helmholtz)
    #conservation
    con=(EnergyList[-1]-EnergyList[0])/np.std(EnergyList)
    
    #write results
    a=df.count()[0]
    df.loc[a]=([rule,dead,var,pf,con,TT,AA])
    df.to_pickle(filename)


    
# Launching automata
cd=pd.read_pickle('rulelist.pkl')
from tqdm import tqdm
from itertools import islice

#for index, row in tqdm(islice(cd.iterrows(), 11219, None)):
for index, row in tqdm(cd.iterrows()):
    birth,death=row[0],row[1]
    run(birth,death)

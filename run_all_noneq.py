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
from tqdm import tqdm
from itertools import islice
import seaborn as sns
import pandas as pd
kernel = np.ones((3, 3), dtype=np.int8)
kernel[1, 1] = 0
x=1             #update rate of the program in ms
NbL = 400 #size of sides
NbL=int(NbL)                                        
NbC = NbL              #square map
part=5 #partitions of each side
Closed_boundaries = True #boundaries
# birth,death=[1,3,6,7],[0,2,4,5] #random rules
# print('rule=',birth,death)
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

# implements rules of propagation, calculates the value of thermodynamical variables 
# at the current time-step
def rules(birth,death,NbL=NbL,NbC=NbC):    
    '''
    Reads the neighbourhood of every cell and updates the values in the CA matrix
    accordingly. Computes the values of the thermodynamical variables and adds
    them to their relative arrays.
    Cannot be exported into modules 
    '''    
    global StateMatrix
    
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
    

    e=np.zeros([part,part])
    #iterates the variable calculation for each box in the prtition
    for i in range (part): #goes through the partition
        for j in range(part):
            Box=(StateMatrix[size*(i):size*(i+1),  size*(j):size*(j+1)]) #takes the boxes
            e[i][j]=np.sum(Box)/(nn)       
    return e

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
    # E=density*N #initial energy init
#arrays for activity, energy, entropy, temperature and pressure for plots and calculations  
    initialize(density,NbL,NbC,StateMatrix)
    E = np.zeros([part,part,100])
    for i in range (100):
        e = rules(birth,death)
        E[:,:,i] = e
        gen += 1
    return E
    
if __name__ == '__main__':
    cd=pd.read_pickle('rulelist.pkl')
    all_values = np.zeros([part,part,100,cd.shape[0]])
    #for index, row in tqdm(islice(cd.iterrows(), 11219, None)):
    for index, row in tqdm(cd.iterrows()):
        birth,death=row[0],row[1]
        E = run(birth,death)
        all_values[:,:,:,index] = E
    #Dimensions of all values are:
        #0,1 are x and y values with the local energy density 
        #2 is the timestep
        #3 is the rule/run
    np.save("all_neq.npy",all_values)
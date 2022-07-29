# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 10:20:39 2020

@author: lucab
"""
from tkinter import Tk, Canvas, Button, LEFT, TOP, Label
import numpy as np
import random
import matplotlib.pyplot as plt

def differentiate(x,y,n=2): #n=2 is onepoint, n=3 is twopoints, etc.
    ''' Rolling derivative of an updating list or array with adjustable step-size
    as opposed to np.gradient which only takes complete lists
    '''
    f=(y[-1]-y[-n])/(x[-1]-x[-n])
    return f  

def rulegenerator():
    '''
    Generates two random lists of random size that do not overlap

    Returns
    -------
    birth : List of size <9 containing random numbers
        
    death : List of size <9 containing random numbers
        

    '''
    n=np.arange(0,9)
    sizeB=np.random.randint(8)
    sizeD=np.random.randint(9-sizeB)
    sizeE=9-sizeB-sizeD
    birth=[]
    death=[]
    for i in range(sizeB+1):
        a=np.random.choice(n)
        birth.append(a)
        n = np.delete(n, np.where(n == a))
    for i in range(sizeD+1):
        b=np.random.choice(n)
        death.append(b)
        n = np.delete(n, np.where(n == b))       
    return birth,death


#Neighbor counting - Moore neighbourhood minus center 
def Moore(posX, posY,StateMatrix,NbC,NbL,Closed_boundaries):  
    #called by rules(), counts the live neighbors to each cell and gives it as output
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


# Drawing the cells (assigning their colors)
def draw(NbL,NbC,StateMatrix,a,modify,canvas):    
    #canvas = Canvas(Tk(), width=a*NbC+modify, height=a*NbL+modify, highlightthickness=0)
    SM1 = np.zeros((NbL,NbC),dtype=int)
    for x in range(NbL):
        for y in range(NbC):
            if StateMatrix[x,y]==0:
                coul = "white"
            elif StateMatrix[x,y]==1:
                coul = "black"
            SM1[x,y] = canvas.create_rectangle((x*a, y*a, (x+1)*a, (y+1)*a), outline="gray", fill="white")
            canvas.itemconfig(SM1[x,y], fill=coul) 
    return canvas

def initialize_map(density,NbL,NbC,StateMatrix,canvas,a,modify):       #creates state matrices based on desnity
    #set initial density 
    SM1 = np.zeros((NbL,NbC),dtype=int)
    if density == 0:     
        StateMatrix[0:NbL,0:NbC] = 0
        for x in range(NbL):
            for y in range(NbC):
                SM1[x,y] = canvas.create_rectangle((x*a, y*a, (x+1)*a, (y+1)*a), outline="gray", fill="white")            
    elif 0 <= density < 1:
        StateMatrix[0:NbL,0:NbC] = 0
        for x in range(NbL):
            for y in range(NbC):
                num=random.random()
                if num<density:
                    StateMatrix[x,y] = 1
                SM1[x,y] = canvas.create_rectangle((x*a, y*a, (x+1)*a, (y+1)*a), outline="gray", fill="white")
    elif density == 1:
        StateMatrix[0:NbL,0:NbC] = 1
        for x in range(NbL):
            for y in range(NbC):
                SM1[x,y] = canvas.create_rectangle((x*a, y*a, (x+1)*a, (y+1)*a), outline="gray", fill="white")            
    else:
        raise Exception('"density" should be between 0 and 1. The value of "density" was: {}'.format(density))   
    
    draw(NbL,NbC,StateMatrix,a,modify,canvas)
    return StateMatrix
    
def initialize(density,NbL,NbC,StateMatrix):
    #set initial density 
    if density == 0:     
        StateMatrix[0:NbL,0:NbC] = 0
    elif 0 <= density < 1:
        StateMatrix[0:NbL,0:NbC] = 0
        for x in range(NbL):
            for y in range(NbC):
                num=random.random()
                if num<density:
                    StateMatrix[x,y] = 1
    elif density == 1:
        StateMatrix[0:NbL,0:NbC] = 1
    else:
        raise Exception('"density" should be between 0 and 1. The value of "density" was: {}'.format(density))
    return StateMatrix


    
        
  
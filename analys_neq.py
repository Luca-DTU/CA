# -*- coding: utf-8 -*-
"""
Created on Wed May 11 15:48:24 2022

@author: luca.bertolani
"""
import numpy as np
val = np.load("all_neq.npy")
#Dimensions of val are:
        #0,1 are x and y values with the local energy density 
        #2 is the timestep
        #3 is the rule/run
val_mean_t = np.mean(val,2)
diff = np.max(val_mean_t,(0,1))-np.min(val_mean_t,(0,1))
plt.hist(diff,bins=20)
diff.max()
len(diff[diff>0.3])

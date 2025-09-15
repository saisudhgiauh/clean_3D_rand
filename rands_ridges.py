#%%
import random
import numpy as np
import sys
import os
sys.path.insert(0,os.getcwd())   #replace this with local directory
import config as c  #no point importing everything here

random.seed()
rands = np.array([random.uniform(-c.max_ridge_angle,c.max_ridge_angle) for _ in range(c.n_ridges)])
np.save("rands",rands)
# %%

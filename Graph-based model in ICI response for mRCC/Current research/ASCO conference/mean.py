import numpy as np
from numpy import mean


DT_s=[0.586 , 0.513 , 0.567 , 0.486]
DT_c=[ 0.513 , 0.621 , 0.540  , 0.513]   
KNN_s=[0.637 ,   0.648  , 0.540  , 0.540]
KNN_c = [ 0.400   , 0.648 , 0.513 , 0.486]    
LR_s=[0.724  , 0.486 ,   0.756  ,   0.756 ]
LR_c=[ 0.729  ,   0.837  ,   0.729  ,   0.864]    
MLP_s=[0.603 , 0.428 , 0.750  , 0.678]
MLP_c=[ 0.486 , 0.459 , 0.459 , 0.729] 
RF_s=[0.693 , 0.567 , 0.594 , 0.594 ]
RF_c=[ 0.567 , 0.648 , 0.540  , 0.540] 

print("DT = " + str(mean(DT_s)-mean(DT_c)))
print("KNN = " + str(mean(KNN_s)-mean(KNN_c)))
print("LR = " + str(mean(LR_s)-mean(LR_c)))
print("MLP = " + str(mean(MLP_s)-mean(MLP_c)))
print("RF = " + str(mean(RF_s)-mean(RF_c)))

print(np.amax(DT_s)-np.amax(DT_c))
print(np.amax(KNN_s)-np.amax(KNN_c))
print(np.amax(LR_s)-np.amax(LR_c))
print(np.amax(MLP_s)-np.amax(MLP_c))
print(np.amax(RF_s)-np.amax(RF_c))
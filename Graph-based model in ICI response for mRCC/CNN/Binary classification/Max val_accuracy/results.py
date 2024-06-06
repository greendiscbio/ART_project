import pandas as pd
import numpy as np
loss = [0.6807, 0.6711, 0.6686, 0.6693, 0.6674]
accuracy = [0.5309, 0.6728, 0.5617, 0.6049, 0.6173]
f1_m = [0.5229, 0.6065, 0.6598, 0.6558, 0.6673]
precision_m = [0.4580, 0.5590, 0.6291, 0.6555, 0.6726]
recall_m = [0.6269, 0.6774, 0.7644, 0.7028, 0.7090]
val_loss = [0.6345, 0.7735,  0.6094, 0.6121, 0.6187]
val_accuracy = [0.6842, 0.3684, 0.7895, 0.6841, 0.6316]
val_f1_m = [0.7273, 0, 0.75, 0.6667, 0.4615]
val_precision_m = [0.6667, 0, 1, 0.75, 1]
val_recall_m = [0.8, 0, 0.6, 0.6, 0.3]

print(np.array(loss).mean())
print(np.array(accuracy).mean())
print(np.array(f1_m).mean())
print(np.array(precision_m).mean())
print(np.array(recall_m).mean())
print(np.array(val_loss).mean())
print(np.array(val_accuracy).mean())
print(np.array(val_f1_m).mean())
print(np.array(val_precision_m).mean())
print(np.array(val_recall_m).mean())

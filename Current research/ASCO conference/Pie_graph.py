from tkinter import font
import matplotlib.pyplot as plt
import numpy as np

y = np.array([4, 8, 4, 67, 2, 1, 2, 4, 8])
mylabels = ["Chromophobe renal cell carcinoma" , "papillary renal cell carcinoma, type 2", "papillary renal cell carcinoma, type 1",  "Clear cell renal cell carcinoma", "Others", "Cystic nephroma", "Multilocular cyst", "Angiomyolipoma", "Oncocytoma"]
myexplode = [0, 0, 0, 0, 0.15, 0.15, 0.15, 0.15, 0.15]
colors = ["lightcoral", "red", "red","firebrick", "darkgrey", "royalblue","royalblue","royalblue", "lightcoral"]
plt.pie(y, labels = mylabels, explode = myexplode, colors=colors, wedgeprops = {"edgecolor" : "black",
                      'linewidth': 0.5,
                      'antialiased': True})
plt.show() 
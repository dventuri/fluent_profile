import numpy as np
import matplotlib.pyplot as plt


values = np.loadtxt('./teste.txt')
n = 64

zz = np.ndarray(shape=(n,n), buffer=np.array(values[:,0]))
yy = np.ndarray(shape=(n,n), buffer=np.array(values[:,1]))
plot_var = -np.ndarray(shape=(n,n), buffer=np.array(values[:,3]))

circle = plt.Circle((0.5,0.5), 0.295275, color='white', fill=False, linewidth=1)
fig, ax = plt.subplots(figsize=(7.5,6))
plot = ax.pcolor(yy, zz, plot_var, edgecolors='k', linewidths=1, cmap='jet')
ax.set_title('Fluid velocity (U) [m/s]')
ax.add_patch(circle)
fig.colorbar(plot)


zz = np.ndarray(shape=(n,n), buffer=np.array(values[:,0]))
yy = np.ndarray(shape=(n,n), buffer=np.array(values[:,1]))
plot_var = np.ndarray(shape=(n,n), buffer=np.array(values[:,2]))

circle = plt.Circle((0.5,0.5), 0.295275, color='white', fill=False, linewidth=1)
fig, ax = plt.subplots(figsize=(7.5,6))
plot = ax.pcolor(yy, zz, plot_var, edgecolors='k', linewidths=1, cmap='jet')
ax.set_title('Fluid velocity (V) [m/s]')
ax.add_patch(circle)
fig.colorbar(plot)


zz = np.ndarray(shape=(n,n), buffer=np.array(values[:,0]))
yy = np.ndarray(shape=(n,n), buffer=np.array(values[:,1]))
plot_var = np.ndarray(shape=(n,n), buffer=np.array(values[:,4]))

circle = plt.Circle((0.5,0.5), 0.295275, color='white', fill=False, linewidth=1)
fig, ax = plt.subplots(figsize=(7.5,6))
plot = ax.pcolor(yy, zz, plot_var, edgecolors='k', linewidths=1, cmap='jet')
ax.set_title('Fluid velocity (W) [m/s]')
ax.add_patch(circle)
fig.colorbar(plot)

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

points = np.loadtxt('profile.axdt',
                        delimiter=',',
                        skiprows=2,
                        max_rows=1080,
                        usecols=(0,1))
points[:,0] -= 3.464  # center x values at 0

vel_z = np.loadtxt('profile.axdt',
                    delimiter=',',
                    skiprows=2,
                    max_rows=1080,
                    usecols=3)

space = np.arange(-0.5,0.5001,0.025)
X,Y = np.meshgrid(space, space)

space_centroid = np.empty(len(space)-1)
for i in range(1,len(space)):
    space_centroid[i-1] = space[i-1] + (space[i]-space[i-1])/2
X_c,Y_c = np.meshgrid(space_centroid, space_centroid)

grid = griddata(points, vel_z,
                (X_c, Y_c),
                method='linear',
                fill_value=0)

fig, ax = plt.subplots()
ax.pcolor(grid.T)
# plt.show()


# circle = plt.Circle((0, 0), 0.33, color='white', fill=False, lw=3)
# fig, ax = plt.subplots(figsize=(7,6))
# c = ax.pcolor(X_c, Y_c, grid,
#             cmap='jet',
#             edgecolors='k',
#             linewidths=1)
# ax.add_patch(circle)
# fig.colorbar(c)
# fig.show()

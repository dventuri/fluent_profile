import sys
import time
import numpy as np
import matplotlib.pyplot as plt

### SIMULATING MFSIM BEHAVIOR AND FUNCTIONS

### read_externalrij(headlev,l1,l2,       -> ins
#                   r11_in,r22_in,r33_in) -> outs

# l1 = 1 // l2 = ltop

r11_in = 0

# my_file=open("profile.in","r",encoding='utf-8')

# nlr11 = my_file.readline()  # number of lines

#allocate
# r11ex = np.empty((nlr11,3))

r11ex = np.loadtxt('profile.in',
                    skiprows=1,
                    delimiter=' ')


def find_closest(data,r_find,theta_find):
    r = data[:,0]
    r_max = np.max(r)

    if(r_find > r_max):
        return 0.0
    
    theta = data[:,1]
    vel_z = data[:,2]
    
    dist = sys.float_info.max

    for i in range(len(r)):
        dr = r[i] - r_find
        dtheta = theta[i] - theta_find
        if(abs(dr) < 0.05 and abs(dtheta) < np.deg2rad(5)):
            distNew = np.sqrt(r_find**2 + r[i]**2 - 2*r_find*r[i]*np.cos(dtheta))
            if(distNew < dist):
                pos = i
                dist = distNew

    # if (dist > 0):
    #     print(dist)

    #     print(r[pos])
    #     print(r_find)
    #     print(theta[pos])
    #     print(theta_find)

    #     print('')

    return vel_z[pos]

def xy2r(x,y):
    r = np.sqrt(np.asarray(x)**2 + np.asarray(y)**2)
    return r

def xy2tetha(x,y):
    theta = np.arctan2(np.asarray(y),np.asarray(x))
    theta = theta % (2*np.pi)
    return theta

r = r11ex[:,0]
theta = r11ex[:,1]
vel_z = r11ex[:,2]
vel_z_from_func = np.empty(len(vel_z))

start_time = time.time()
for i in range(len(r11ex)):
    vel_z_from_func[i] = find_closest(r11ex,r[i],theta[i])
print("--- %s seconds ---" % (time.time() - start_time))

# vel_z_error = vel_z - vel_z_from_func
# fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
# plot = ax.scatter(theta, r, c=vel_z_error,
#                 cmap="seismic")
# ax.set_rmax(0.33)
# fig.colorbar(plot)

# Make 2D grid - simulating MFSim

space = np.arange(-0.5,0.5,0.025)
X,Y = np.meshgrid(space, space)

space_centroid = np.empty(len(space)-1)
for i in range(1,len(space)):
    space_centroid[i-1] = space[i-1] + (space[i]-space[i-1])/2
X_c,Y_c = np.meshgrid(space_centroid, space_centroid)

R = xy2r(X_c, Y_c)
THETA = xy2tetha(X_c, Y_c)

Z = np.empty((len(space_centroid),len(space_centroid)))
for i in range(len(space_centroid)):
    for j in range(len(space_centroid)):
        Z[i,j] = find_closest(r11ex,R[i,j],THETA[i,j])

circle = plt.Circle((0, 0), 0.33, color='k', fill=False, lw=2)
fig, ax = plt.subplots(figsize=(7,6))
c = ax.pcolor(X, Y, Z,
            cmap='jet')
ax.add_patch(circle)
fig.colorbar(c)

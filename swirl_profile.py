import numpy as np
import matplotlib.pyplot as plt

def func_R(y,z):
    return np.sqrt(z**2 + y**2)

def y_vel_MFSim(y, z, aux=None):
    if(not aux):
        aux = func_R(y,z)
    return aux*2*np.pi*150000*(z*np.sin(np.pi/2)+y*np.cos(np.pi/2))

def z_vel_MFSim(y, z, aux=None):
    if(not aux):
        aux = func_R(y,z)
    return aux*2*np.pi*150000*(-y*np.sin(np.pi/2)+z*np.cos(np.pi/2))

def y_vel(y, z, aux=None):
    if(not aux):
        aux = func_R(y,z)
    return -aux*150000*np.sin(np.pi*z)*np.cos(np.pi*y)

def z_vel(y, z, aux=None):
    if(not aux):
        aux = func_R(y,z)
    return aux*150000*np.sin(np.pi*y)*np.cos(np.pi*z)


r_jet = 0.0055

start = 0.5-r_jet
stop = 0.5+r_jet

space = np.linspace(start,stop,50)
space = space-0.5
yy, zz = np.meshgrid(space, space)

R = func_R(yy, zz)

# auxiliary function
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, R)
fig.colorbar(plot)

# MFSim - y and z velocities with aux function
vel_y_MFSim = y_vel_MFSim(yy,zz)
vel_z_MFSim = z_vel_MFSim(yy,zz)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y_MFSim)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z_MFSim)
fig.colorbar(plot)

fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y_MFSim, vel_z_MFSim,
                    density=1,
                    color=np.sqrt(vel_y_MFSim**2 + vel_z_MFSim**2),
                    cmap='viridis'
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)

# MFSim - y and z velocities without aux function
vel_y_MFSim = y_vel_MFSim(yy,zz,1)
vel_z_MFSim = z_vel_MFSim(yy,zz,1)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y_MFSim)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z_MFSim)
fig.colorbar(plot)

fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y_MFSim, vel_z_MFSim,
                    density=1,
                    color=np.sqrt(vel_y_MFSim**2 + vel_z_MFSim**2),
                    cmap='viridis'
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)

# TEST - y and z velocities with aux function
vel_y = y_vel(yy,zz)
vel_z = z_vel(yy,zz)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z)
fig.colorbar(plot)

fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y, vel_z,
                    density=1,
                    color=np.sqrt(vel_y**2 + vel_z**2),
                    cmap='viridis'
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)

# TEST - y and z velocities without aux function
vel_y = y_vel(yy,zz,aux=1)
vel_z = z_vel(yy,zz,aux=1)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z)
fig.colorbar(plot)

fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y, vel_z,
                    density=1,
                    color=np.sqrt(vel_y**2 + vel_z**2),
                    cmap='viridis'
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)
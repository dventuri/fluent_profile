import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors

def func_R(y, z):
    return np.sqrt(y**2 + z**2)

def func_R_norm(y, z, R):
    aux = np.sqrt(y**2 + z**2)
    return aux/R

def y_vel_MFSim(y, z, C=150000, aux=None):
    if(not aux):
        aux = func_R(y,z)
    return aux*2*np.pi*C*z

def z_vel_MFSim(y, z, C=150000, aux=None):
    if(not aux):
        aux = func_R(y,z)
    return aux*2*np.pi*C*(-y)

def y_vel_MFSim_norm(y, z, R, angle, axial_vel, aux=None):
    if(not aux):
        aux = func_R_norm(y,z,R)
    return aux*np.tan(angle/2)*axial_vel*z/R

def z_vel_MFSim_norm(y, z, R, angle, axial_vel, aux=None):
    if(not aux):
        aux = func_R_norm(y,z,R)
    return aux*np.tan(angle/2)*axial_vel*(-y)/R

def y_vel(y, z, C=150000, aux=None):
    if(not aux):
        aux = func_R(y,z)
    return aux*C*2*np.sin(np.pi*z)*np.cos(np.pi*y)
        #  aux*C*2*np.pi*z

def z_vel(y, z, C=150000, aux=None):
    if(not aux):
        aux = func_R(y,z)
    return -aux*C*2*np.sin(np.pi*y)*np.cos(np.pi*z)
        #  -aux*C*2*np.pi*y

# Small angle approximation
# sin(theta) = theta
# cos(theta) = 1 - theta^2/2 = 1

r_jet = 0.0055

start = -r_jet
stop = r_jet

space = np.linspace(start,stop,50)
yy, zz = np.meshgrid(space, space)

R = func_R(yy, zz)

# auxiliary function
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, R)
fig.colorbar(plot)

# MFSim - y and z velocities with aux function
vel_y = y_vel_MFSim(yy,zz)
vel_z = z_vel_MFSim(yy,zz)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z)
fig.colorbar(plot)

max_vel = y_vel_MFSim(0,r_jet)
vel_mag = np.sqrt(vel_y**2 + vel_z**2)
fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y, vel_z,
                    density=1,
                    color=vel_mag,
                    cmap='viridis',
                    norm=colors.Normalize(vmin=0,
                                          vmax=max_vel)
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)

dvy_dz = np.gradient(vel_y, space, axis=0, edge_order=2)
dvz_dy = np.gradient(vel_z, space, axis=1, edge_order=2)
vort = (dvz_dy-dvy_dz)
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vort)
fig.colorbar(plot)

# MFSim - y and z velocities normalized
angle = np.deg2rad(90)
vel_y = y_vel_MFSim_norm(yy,zz,r_jet,angle,12.66)
vel_z = z_vel_MFSim_norm(yy,zz,r_jet,angle,12.66)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z)
fig.colorbar(plot)

max_vel = y_vel_MFSim_norm(0,r_jet,r_jet,angle,12.66)
vel_mag = np.sqrt(vel_y**2 + vel_z**2)
fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y, vel_z,
                    density=1,
                    color=vel_mag,
                    cmap='viridis',
                    norm=colors.Normalize(vmin=0,
                                          vmax=max_vel)
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)

dvy_dz = np.gradient(vel_y, space, axis=0, edge_order=2)
dvz_dy = np.gradient(vel_z, space, axis=1, edge_order=2)
vort = (dvz_dy-dvy_dz)
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vort)
fig.colorbar(plot)

# MFSim - y and z velocities without aux function
vel_y = y_vel_MFSim(yy,zz,150000,1)
vel_z = z_vel_MFSim(yy,zz,150000,1)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z)
fig.colorbar(plot)

max_vel = y_vel_MFSim(0,r_jet,150000,1)
vel_mag = np.sqrt(vel_y**2 + vel_z**2)
fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y, vel_z,
                    density=1,
                    color=vel_mag,
                    cmap='viridis',
                    norm=colors.Normalize(vmin=0,
                                          vmax=max_vel)
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)

dvy_dz = np.gradient(vel_y, space, axis=0, edge_order=2)
dvz_dy = np.gradient(vel_z, space, axis=1, edge_order=2)
vort = (dvz_dy-dvy_dz)
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vort)
fig.colorbar(plot)

# TEST - y and z velocities with aux function
vel_y = y_vel(yy,zz)
vel_z = z_vel(yy,zz)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z)
fig.colorbar(plot)

max_vel = y_vel(0,r_jet)
vel_mag = np.sqrt(vel_y**2 + vel_z**2)
fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y, vel_z,
                    density=1,
                    color=vel_mag,
                    cmap='viridis',
                    norm=colors.Normalize(vmin=0,
                                          vmax=max_vel)
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)

dvy_dz = np.gradient(vel_y, space, axis=0, edge_order=2)
dvz_dy = np.gradient(vel_z, space, axis=1, edge_order=2)
vort = (dvz_dy-dvy_dz)
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vort)
fig.colorbar(plot)

# TEST - y and z velocities without aux function
vel_y = y_vel(yy,zz,150000,1)
vel_z = z_vel(yy,zz,150000,1)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_y)
fig.colorbar(plot)

fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vel_z)
fig.colorbar(plot)

max_vel = y_vel(0,r_jet,150000,1)
vel_mag = np.sqrt(vel_y**2 + vel_z**2)
fig, ax = plt.subplots()
strm = ax.streamplot(yy, zz, vel_y, vel_z,
                    density=1,
                    color=vel_mag,
                    cmap='viridis',
                    norm=colors.Normalize(vmin=0,
                                          vmax=max_vel)
                    )
ax.set_xlim(-r_jet,+r_jet)
ax.set_ylim(-r_jet,+r_jet)
fig.colorbar(strm.lines)

dvy_dz = np.gradient(vel_y, space, axis=0, edge_order=2)
dvz_dy = np.gradient(vel_z, space, axis=1, edge_order=2)
vort = (dvz_dy-dvy_dz)
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vort)
fig.colorbar(plot)

vort = -3*2*np.pi*150000*np.sqrt(yy**2+zz**2)
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vort)
fig.colorbar(plot)

vort = np.ones((len(space),len(space)))*(-2)*2*np.pi*150000
fig, ax = plt.subplots()
plot = ax.contourf(yy, zz, vort)
fig.colorbar(plot)

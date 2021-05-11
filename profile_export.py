import numpy as np


### READ FILE
my_file=open("profile.prof","r",encoding='utf-8')

file_header = my_file.readline()
print(file_header)

print("Reading coordinate: ", my_file.readline()[1])
x = []
for f in my_file:
    if (f == ')\n'):
        break
    x.append(float(f[:-1])-3.464)

print("Reading coordinate: ", my_file.readline()[1])
y = []
for f in my_file:
    if (f == ')\n'):
        break
    y.append(float(f[:-1]))

print("Jumping coordinate: ", my_file.readline()[1])
for f in my_file:
    if (f == ')\n'):
        break

print("Reading: ", my_file.readline()[1:])
vel_z = []
for f in my_file:
    if (f == ')\n'):
        break
    vel_z.append(float(f[:-1])*(-1))

my_file.close()


### CONVERT TO POLAR COORDINATES (AND TO NUMPY ARRAYS)
r = np.sqrt(np.asarray(x)**2 + np.asarray(y)**2)
theta = np.arctan2(np.asarray(y),np.asarray(x))
vel_z = np.asarray(vel_z)

### CONVERT THETA SIGNED ANGLES TO 0-360 INTERVAL
theta = theta % (2*np.pi)


### WRITE FILE
np.savetxt("profile.in",
           np.column_stack([r,theta,vel_z]),
           header=str(len(r)),
           comments='')

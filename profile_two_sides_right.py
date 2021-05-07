import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


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


### CONVERT TO POLAR COORDINATES
r = np.sqrt(np.asarray(x)**2 + np.asarray(y)**2)
theta = np.arctan2(np.asarray(y),np.asarray(x))


### FILTER R > 0.30
filter = r < 0.3
r = r[filter]
theta = theta[filter]
vel_z = np.asarray(vel_z)[filter]


### FILTER QUADRANTS I AND IV
filter = theta <= (np.pi/2+0.0001)    # filter out quadrant II
r = r[filter]
theta = theta[filter]
vel_z = np.asarray(vel_z)[filter]

filter = theta >= -(np.pi/2+0.0001)   #filter out quadrant III
r = r[filter]
theta = theta[filter]
vel_z = np.asarray(vel_z)[filter]


### PLOT FLUENT PROFILE
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
plot = ax.scatter(theta, r, c=vel_z,
            # cmap="YlGnBu",
            s=10)
ax.set_rmax(np.max(r))
fig.colorbar(plot)


### FIT POLYNOMIAL EQUATIO
# make array X as [[r1 theta1], [r2 theta2], [r3 theta3]]
# degree=2 will turn X=[r, theta] into:
# [1, r, theta, r^2, r*theta, theta^2] -> linear!

from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import Pipeline
import numpy as np
model = Pipeline([('poly', PolynomialFeatures(degree=3)),
                  ('linear', LinearRegression(fit_intercept=False))])
# fit to an order-3 polynomial data
x = np.column_stack([r, theta])
y = np.asarray(vel_z)
model = model.fit(x, y)
# model.named_steps['linear'].coef_
print(model.score(x, y))
# model.predict(x)


### PLOT FITTED MODEL (SCATTER)
vel_z_predicted = model.predict(x)
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
plot = ax.scatter(theta, r, c=vel_z_predicted,
            # cmap="YlGnBu",
            s=10)
ax.set_rmax(np.max(r))
fig.colorbar(plot)


### PLOT ERROR (PROFILE - FITTED)
vel_z_error = (vel_z - vel_z_predicted)
divnorm = mcolors.TwoSlopeNorm(vmin=np.min(vel_z_error),
                               vcenter=0,
                               vmax=np.max(vel_z_error))
fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
plot = ax.scatter(theta, r, c=vel_z_error,
            cmap="seismic",
            norm=divnorm,
            s=2)
ax.set_rmax(0.33)
fig.colorbar(plot)


### PLOT FITTED MODEL (CONTOUR MAP)
r_list = np.linspace(np.min(r), np.max(r), 1000)
theta_list = np.linspace(np.min(theta), np.max(theta), 1000)
r_mesh, theta_mesh = np.meshgrid(r_list, theta_list)

def Z(r, theta):
    return r*0.00001 + theta*0.00001

vel_z_predicted = model.predict(
    np.array(
        [r_mesh.ravel(), theta_mesh.ravel()]
    ).T
)
vel_z_predicted = vel_z_predicted.reshape(r_mesh.shape)

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
plot = ax.contourf(theta_mesh, r_mesh,
                    vel_z_predicted, 10,
                    vmin=8, vmax=23
                    )
plot = ax.scatter(theta, r, c=vel_z,
            s=10)
ax.set_rmax(0.33)
fig.colorbar(plot)

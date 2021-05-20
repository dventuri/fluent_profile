# From https://www.particleincell.com/2012/quad-interpolation

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


# First of all, read CFD-Post data as x and y - no need for r and theta
# Also, these data is on vertices and has connectivity info
# Fluent raw data is on centroids but has no connec. info

vertices = np.loadtxt('profile.axdt',
                        delimiter=',',
                        skiprows=2,
                        max_rows=1080,
                        usecols=(0,1))
vertices[:,0] -= 3.464  # center x values at 0

vel_z = np.loadtxt('profile.axdt',
                    delimiter=',',
                    skiprows=2,
                    max_rows=1080,
                    usecols=3)

vertices_connect = np.loadtxt('profile.axdt',
                                dtype=int,
                                delimiter=',',
                                skiprows=1084)



# %create our polygon
# px = [-1, 8, 13, -4];
# py = [-1, 3, 11, 8];

# %compute coefficients
# A=[1 0 0 0;1 1 0 0;1 1 1 1;1 0 1 0];
# AI = inv(A);
# a = AI*px';
# b = AI*py';







# Define function for finding the cell that contains a given point
def find_cell(point,vertices,vertices_connect):

    for n in range(len(vertices_connect)):

        p = []
        for i in range(4):
            vert_index = vertices_connect[n,i]
            x_vert = vertices[vert_index,0]
            y_vert = vertices[vert_index,1]
            p.append((x_vert, y_vert))

        current_cell = Polygon(p)

        if current_cell.contains(point):
            return n, current_cell

    return None, Polygon()


# Define function to find value at (x,y)
def calc_alpha_beta(vertex_1, vertex_2, vertex_3, vertex_4, point):

    a = -vertex_1[0] + vertex_3[0]
    b = -vertex_1[0] + vertex_2[0]
    c = vertex_1[0] - vertex_2[0] - vertex_3[0] + vertex_4[0]
    d = point.x - vertex_1[0]

    e = -vertex_1[1] + vertex_3[1]
    f = -vertex_1[1] + vertex_2[1]
    g = vertex_1[1] - vertex_2[1] - vertex_3[1] + vertex_4[1]
    h = point.y - vertex_1[1]

    # soeq_a =

    def alpha_1(a,b,c,d,e,f,g,h):
        alpha = -(b*e - a*f + d*g -c*h + np.sqrt(
            -4*(c*e - a*g)*(d*f - b*h) + (
                b*e - a*f + d*g - c*h
            )**2
        ))/(2*c*e - 2*a*g)
        return alpha

    def beta_1(a,b,c,d,e,f,g,h):
        beta = (b*e - a*f - d*g + c*h + np.sqrt(
            -4*(c*e - a*g)*(d*f - b*h) + (
                b*e - a*f + d*g - c*h
            )**2
        ))/(2*c*f - 2*b*g)
        return beta

    def alpha_2(a,b,c,d,e,f,g,h):
        alpha = (-b*e + a*f - d*g + c*h + np.sqrt(
            -4*(c*e - a*g)*(d*f - b*h) + (
                b*e - a*f + d*g - c*h
            )**2
        ))/(2*c*e - 2*a*g)
        return alpha

    def beta_2(a,b,c,d,e,f,g,h):
        beta = -((-b*e + a*f + d*g - c*h + np.sqrt(
            -4*(c*e - a*g)*(d*f - b*h) + (
                b*e - a*f + d*g - c*h
            )**2
        ))/(2*c*f - 2*b*g))
        return beta

    def calc_alpha_beta_1(a,b,c,d,e,f,g,h):
        alpha = alpha_1(a,b,c,d,e,f,g,h)
        if(alpha >= 0 and alpha <= 1):
            beta = beta_1(a,b,c,d,e,f,g,h)
            if(beta >= 0 or beta <= 1):
                return (alpha, beta)
        return None

    def calc_alpha_beta_2(a,b,c,d,e,f,g,h):
        alpha = alpha_2(a,b,c,d,e,f,g,h)
        if(alpha >= 0 and alpha <= 1):
            beta = beta_2(a,b,c,d,e,f,g,h)
            if(beta >= 0 or beta <= 1):
                return (alpha, beta)
        return None

    coeff = calc_alpha_beta_1(a,b,c,d,e,f,g,h)
    if(coeff):
        return coeff
    else:
        coeff = calc_alpha_beta_2(a,b,c,d,e,f,g,h)
        if(coeff):
            return coeff
        else:
            return None

def interpolate_value(x, y, vertices, vertices_connect):

    point = Point(x, y)

    n_cell, cell = find_cell(point, vertices, vertices_connect)

    if(n_cell):

        vertices_x = np.empty(4)
        vertices_y = np.empty(4)
        vertices_coord = cell.exterior.coords
        for i in range(4):
            vertices_x[i] = vertices_coord[i][0]
            vertices_y[i] = vertices_coord[i][1]

        A = np.array([[1, 0, 0, 0],
                      [1, 1, 0, 0],
                      [1, 1, 1, 1],
                      [1, 0, 1, 0]])
        AI = np.linalg.inv(A)

        a = AI.dot(vertices_x)
        b = AI.dot(vertices_y)

        # quadratic equation coeffs, aa*mm^2+bb*m+cc=0
        aa = a[3]*b[2] - a[2]*b[3]
        bb = a[3]*b[0] - a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + point.x*b[3] - point.y*a[3]
        cc = a[1]*b[0] - a[0]*b[1] + point.x*b[1] - point.y*a[1]

        # compute m = (-b+sqrt(b^2-4ac))/(2a)
        det = np.sqrt(bb*bb - 4*aa*cc)
        m = (-bb - det)/(2*aa)

        # compute l
        l = (point.x - a[0] - a[2]*m)/(a[1] + a[3]*m)

        p1 = vel_z[vertices_connect[n_cell,0]]
        p2 = vel_z[vertices_connect[n_cell,1]]
        p3 = vel_z[vertices_connect[n_cell,2]]
        p4 = vel_z[vertices_connect[n_cell,3]]

        value = (1 - m)*(
            (1 - l)*p1 + l*p2
        ) + m*(
            l*p3 + (1-l)*p4
        )

        return value

    return 0


# value = interpolate_value(0.5, -0.2, vertices, vertices_connect)
# print(value)

# Define funciton to create MFSim-like grid
# space = np.arange(-0.5,0.5001,0.025)
space = np.arange(-0.5,0.5001,0.003125)
X,Y = np.meshgrid(space, space)

space_centroid = np.empty(len(space)-1)
for i in range(1,len(space)):
    space_centroid[i-1] = space[i-1] + (space[i]-space[i-1])/2
X_c,Y_c = np.meshgrid(space_centroid, space_centroid)

Z = np.empty((len(space_centroid),len(space_centroid)))
for i in range(len(space_centroid)):
    for j in range(len(space_centroid)):
        Z[i,j] = interpolate_value(X_c[i,j], Y_c[i,j], vertices, vertices_connect)

circle = plt.Circle((0, 0), 0.33, color='white', fill=False, lw=3)
fig, ax = plt.subplots(figsize=(7,6))
c = ax.pcolor(X_c, Y_c, -Z,
            cmap='jet',
            edgecolors='k',
            linewidths=1)
ax.add_patch(circle)
fig.colorbar(c)

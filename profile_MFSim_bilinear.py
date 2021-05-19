# From https://stackoverflow.com/a/23921432/3980223
# Bilinear interpolation for any convex tetragon

import numpy as np
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
def calc_alpha_beta(vertex_1, vertex_2, vertex_3, vertex_4, centroid):

    a = -vertex_1[0] + vertex_3[0]
    b = -vertex_1[0] + vertex_2[0]
    c = vertex_1[0] - vertex_2[0] - vertex_3[0] + vertex_4[0]
    d = centroid[0] - vertex_1[0]

    e = -vertex_1[1] + vertex_3[1]
    f = -vertex_1[1] + vertex_2[1]
    g = vertex_1[1] - vertex_2[1] - vertex_3[1] + vertex_4[1]
    h = centroid[1] - vertex_1[1]

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

        vertex = cell.exterior.coords
        v1 = vertex[0]
        v2 = vertex[1]
        v3 = vertex[2]
        v4 = vertex[3]

        centroid = cell.centroid.coords[0]

        alpha, beta = calc_alpha_beta(v1, v2, v3, v4, centroid)

        p1 = vel_z[vertices_connect[n_cell,0]]
        p2 = vel_z[vertices_connect[n_cell,1]]
        p3 = vel_z[vertices_connect[n_cell,2]]
        p4 = vel_z[vertices_connect[n_cell,3]]

        value = (1 - alpha)*(
            (1 - beta)*p1 + beta*p2
        ) + alpha*(
            (1 - beta)*p3 + beta*p4
        )

        return value

    return 0


# value = interpolate_value(0.5, -0.2, vertices, vertices_connect)
# print(value)

# Define funciton to create MFSim-like grid
space = np.arange(-0.5,0.5001,0.025)
X,Y = np.meshgrid(space, space)

space_centroid = np.empty(len(space)-1)
for i in range(1,len(space)):
    space_centroid[i-1] = space[i-1] + (space[i]-space[i-1])/2
X_c,Y_c = np.meshgrid(space_centroid, space_centroid)

Z = np.empty((len(space_centroid),len(space_centroid)))
for i in range(len(space_centroid)):
    for j in range(len(space_centroid)):
        Z[i,j] = interpolate_value(X_c[i,j], Y_c[i,j], vertices, vertices_connect)

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

    return -1


# Define function to find value at (x,y)
def interpolate_value(x, y, vertices, vertices_connect):

    point = Point(x, y)

    n_cell, cell = find_cell(point, vertices, vertices_connect)
    print(n_cell)

    vertex = cell.exterior.coords
    # vertex = vertices
    center = cell.centroid.coords[0]

    a = -vertex[0][0] + vertex[2][0]
    b = -vertex[0][0] + vertex[1][0]
    c = vertex[0][0] - vertex[1][0] - vertex[2][0] + vertex[3][0]
    d = center[0] - vertex[0][0]
    e = -vertex[0][1] + vertex[2][1]
    f = -vertex[0][1] + vertex[1][1]
    g = vertex[0][1] - vertex[1][1] - vertex[2][1] + vertex[3][1]
    h = center[1] - vertex[0][1]

    alpha = -(b*e - a*f + d*g -c*h + np.sqrt(
        -4*(c*e - a*g)*(d*f - b*h) + (
            b*e - a*f + d*g - c*h
        )**2
    ))/(2*c*e - 2*a*g)

    if(alpha >= 0 or alpha <= 1):
        beta = (b*e - a*f - d*g + c*h + np.sqrt(
            -4*(c*e - a*g)*(d*f - b*h) + (
                b*e - a*f + d*g - c*h
            )**2
        ))/(2*c*f - 2*b*g)

        if(beta >= 0 or beta <= 1):
            pass
        else:
            return -1

    else:
        alpha = (-b*e + a*f - d*g + c*h + np.sqrt(
            -4*(c*e - a*g)*(d*f - b*h) + (
                b*e - a*f + d*g - c*h
            )**2
        ))/(2*c*e - 2*a*g)

        if(alpha >= 0 or alpha <= 1):
            beta = -((-b*e + a*f + d*g - c*h + np.sqrt(
                -4*(c*e - a*g)*(d*f - b*h) + (
                    b*e - a*f + d*g - c*h
                )**2
            ))/(2*c*f - 2*b*g))

            if(beta >= 0 or beta <= 1):
                pass
            else:
                return -1
        
        else:
            return -1

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

value = interpolate_value(0.23, -0.2, vertices, vertices_connect)
print(value)

# Define funciton to create MFSim-like grid

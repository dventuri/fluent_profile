# From https://www.particleincell.com/2012/quad-interpolation

from atpbar.main import atpbar
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from atpbar import atpbar


def read_cfdpost_data(filename, number_of_lines, delimiter=',', skiprows=0):
    # Read CFD-Post vertex data as x and y points with connectivity info
    # Fluent raw data is at centroids but has no connectivity info

    vertices = np.loadtxt(filename,
                          delimiter=delimiter,
                          skiprows=skiprows,
                          max_rows=number_of_lines,
                          usecols=(1,2))

    vel_z = np.loadtxt(filename,
                       delimiter=delimiter,
                       skiprows=skiprows,
                       max_rows=number_of_lines,
                       usecols=(4,5,6))

    vertices_connect = np.loadtxt(filename,
                                  dtype=int,
                                  delimiter=delimiter,
                                  skiprows=number_of_lines+skiprows+2)

    return vertices, vertices_connect, vel_z


def find_cell(point,vertices,vertices_connect):
    # Find the cell that contains a given point

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


def calc_normalized_coord(cell, point):
    # calculate normalized coordinates (rectangular coordinates 0-1)

    vertices_x = np.empty(4)
    vertices_y = np.empty(4)
    vertices_coord = cell.exterior.coords
    for i in range(4):
        vertices_x[i] = vertices_coord[i][0]
        vertices_y[i] = vertices_coord[i][1]

    AI = np.array([[ 1, 0, 0, 0],
                    [-1, 1, 0, 0],
                    [-1, 0, 0, 1],
                    [ 1,-1, 1,-1]])

    a = np.dot(AI,vertices_x)
    b = np.dot(AI,vertices_y)

    # quadratic equation coeffs, aa*mm^2+bb*m+cc=0
    aa = a[3]*b[2] - a[2]*b[3]
    bb = a[3]*b[0] - a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + point.x*b[3] - point.y*a[3]
    cc = a[1]*b[0] - a[0]*b[1] + point.x*b[1] - point.y*a[1]

    # compute m = (-b+sqrt(b^2-4ac))/(2a)
    det = np.sqrt(bb*bb - 4*aa*cc)
    m = (-bb - det)/(2*aa)
    if(not 0 <= m <= 1):
        m = (-bb + det)/(2*aa)

    # compute l
    l = (point.x - a[0] - a[2]*m)/(a[1] + a[3]*m)

    return l, m

def calculate_value(n_cell, cell, point, vertices_connect, values):

    if(n_cell):

        l, m = calc_normalized_coord(cell, point)

        p1 = values[vertices_connect[n_cell,0]]
        p2 = values[vertices_connect[n_cell,1]]
        p3 = values[vertices_connect[n_cell,2]]
        p4 = values[vertices_connect[n_cell,3]]

        value = (1 - m)*(
            (1 - l)*p1 + l*p2
        ) + m*(
            l*p3 + (1-l)*p4
        )

        return value

    return 0


def interpolate_value(x, y, vertices, vertices_connect, values):
    # interpolate value from vertices to the inside point

    point = Point(x, y)

    n_cell, cell = find_cell(point, vertices, vertices_connect)

    ndim = values.shape[1]

    if(ndim > 1):
        value = np.empty(ndim)
        for n in range(ndim):
            value[n] = calculate_value(n_cell, cell, point, vertices_connect, values[:,n])
    else:
        value = calculate_value(n_cell, cell, point, vertices_connect, values)

    return value


if __name__ == '__main__':

    # geometric definitions
    r = 0.64140/2
    r2 = r**2
    Ly = Lz = 0.8

    # MFSim mesh definitions
    n_cells = 4096
    dy = dz = Ly/n_cells

    # read cfd-post data
    vertices, vertices_connect, vel = read_cfdpost_data(
        'profile_outlet_fluent_kwSST.csv',
        1080,
        delimiter=',',
        skiprows=15
    )

    # center values at (L/2, L/2)
    vertices[:,0] -= (3.464 - Ly/2) #TODO: Check x center
    vertices[:,1] += Lz/2

    # writing at the centroid of the cell
    values = np.empty((n_cells*n_cells,5))
    cell_range = range(1,n_cells+1)
    i = 0
    for j in atpbar(cell_range):
        for k in cell_range:
            y_c = j*dy
            z_c = k*dz
            line = (j-1)*n_cells + k

            values[i,0] = y_c
            values[i,1] = z_c

            d2 = (y_c-Ly/2)**2 + (z_c-Lz/2)**2
            if(d2 < (r2*1.05)):
                values[i,2:] = interpolate_value(y_c, z_c,
                                                 vertices,
                                                 vertices_connect,
                                                 vel)
            else:
                values[i,2:] = 0
            
            i += 1

    # save data
    np.savetxt("vel.txt", values)

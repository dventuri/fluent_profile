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
                       usecols=(3,4,5))

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


def interpolate_value(x, y, vertices, vertices_connect, values):
    # interpolate value from vertices to the inside point

    point = Point(x, y)

    n_cell, cell = find_cell(point, vertices, vertices_connect)

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


if __name__ == '__main__':

    # geometric definitions
    r = 0.64140/2
    r2 = r**2
    Ly = Lz = 1

    # MFSim mesh definitions
    n_cells_lbot = 40
    mesh_levels = 6

    # read cfd-post data
    vertices, vertices_connect, vel = read_cfdpost_data(
        'profile_outlet_fluent_kwSST.csv',
        1080,
        delimiter=',',
        skiprows=15
    )

    # center values at (L/2, L/2)
    vertices[:,0] -= (3.464 - Ly/2) #TODO: Check x center
    vertices[:,1] += 0.5

    ### Interpolate values from fluent mesh to MFSim mesh ###

    # index 0,1,2 in vel are related to x,y,z directions
    # however, they are translated to z,y,x into MFSim coordinates

    # for each mesh level
    for lvl in range(mesh_levels):

        n_cells = n_cells_lbot*(2**(lvl))
        dy = dz = Ly/n_cells

        cell_range = range(1,n_cells+1)

        # u in MFSim is vel_z from fluent
        # writing at the centroid of the cell
        with open('u_'+str(lvl+1)+'.txt', 'w') as myfile:
            for j in atpbar(cell_range):
                for k in cell_range:
                    y_c = j*dy
                    z_c = k*dz

                    d2 = (y_c-Ly/2)**2 + (z_c-Lz/2)**2
                    if(d2 < (r2*1.05)):
                        value = interpolate_value(y_c, z_c,
                                                  vertices,
                                                  vertices_connect,
                                                  vel[:,2])
                    else:
                        value = 0

                    myfile.write("{}, {}, {}\n".format(
                        j, k, value
                    ))

        # v in MFSim is vel_y from fluent
        # writing at the west face of the cell
        with open('v_'+str(lvl+1)+'.txt', 'w') as myfile:
            for j in atpbar(cell_range):
                for k in cell_range:
                    y_f = (j-0.5)*dy
                    z_c = k*dz

                    d2 = (y_f-Ly/2)**2 + (z_c-Lz/2)**2
                    if(d2 < (r2*1.05)):
                        value = interpolate_value(y_f, z_c,
                                                  vertices,
                                                  vertices_connect,
                                                  vel[:,1])
                    else:
                        value = 0

                    myfile.write("{}, {}, {}\n".format(
                        j, k, value
                    ))

        # w in MFSim is vel_x from fluent
        # writing at the bot face of the cell
        with open('w_'+str(lvl+1)+'.txt', 'w') as myfile:
            for j in atpbar(cell_range):
                for k in cell_range:
                    y_c = j*dy
                    z_f = (k-0.5)*dz

                    d2 = (y_c-Ly/2)**2 + (z_f-Lz/2)**2
                    if(d2 < (r2*1.05)):
                        value = interpolate_value(y_c, z_f,
                                                  vertices,
                                                  vertices_connect,
                                                  vel[:,0])
                    else:
                        value = 0

                    myfile.write("{}, {}, {}\n".format(
                        j, k, value
                    ))

    # # create MFSim-like grid (vertices)
    # space = np.arange(-0.5,0.5001,0.025)
    # X,Y = np.meshgrid(space, space)

    # # define centroids
    # space_centroid = np.empty(len(space)-1)
    # for i in range(1,len(space)):
    #     space_centroid[i-1] = space[i-1] + (space[i]-space[i-1])/2
    # X_c,Y_c = np.meshgrid(space_centroid, space_centroid)

    # # interpolate
    # percent = 0
    # n = 0
    # Z = np.empty((len(space_centroid),len(space_centroid)))
    # total = len(Z)**2
    # for i in range(len(space_centroid)):
    #     for j in range(len(space_centroid)):
    #         Z[i,j] = interpolate_value(X_c[i,j], Y_c[i,j], vertices, vertices_connect, vel_z)
    #         n += 1
    #         percent = n/total*100
    #         print(percent,"%")

    # # plot
    # circle = plt.Circle((0, 0), 0.33, color='white', fill=False, lw=3)
    # fig, ax = plt.subplots(figsize=(7,6))
    # c = ax.pcolor(X, Y, -Z,
    #             cmap='jet',
    #             edgecolors='k',
    #             linewidths=1)
    # ax.add_patch(circle)
    # fig.colorbar(c)

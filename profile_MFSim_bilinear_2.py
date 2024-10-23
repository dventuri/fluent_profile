# From https://www.particleincell.com/2012/quad-interpolation
# Correction for singularity in vertices ordering and time estimation added by Claude.ai

import numpy as np
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
import time


def read_cfdpost_data(filename, number_of_lines, delimiter=',', skiprows=0):
    # Read CFD-Post vertex data as x and y points with connectivity info
    # Fluent raw data is at centroids but has no connectivity info

    vertices = np.loadtxt(filename,
                          delimiter=delimiter,
                          skiprows=skiprows,
                          max_rows=number_of_lines,
                          usecols=(1,3))

    vel_vector = np.loadtxt(filename,
                            delimiter=delimiter,
                            skiprows=skiprows,
                            max_rows=number_of_lines,
                            usecols=(4,5,6))

    vertices_connect = np.loadtxt(filename,
                                  dtype=int,
                                  delimiter=delimiter,
                                  skiprows=number_of_lines+skiprows+2)

    return vertices, vertices_connect, vel_vector


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
    """
    Calculate normalized coordinates (rectangular coordinates 0-1)
    with additional checks for numerical stability
    """
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

    # Try original vertex ordering
    l, m = try_calculate_coord(vertices_x, vertices_y, point, AI)
    if l is not None and 0 <= l <= 1 and 0 <= m <= 1:
        return l, m

    # If original ordering fails, try rotating vertices
    for i in range(1, 4):
        vertices_x_rot = np.roll(vertices_x, i)
        vertices_y_rot = np.roll(vertices_y, i)
        l, m = try_calculate_coord(vertices_x_rot, vertices_y_rot, point, AI)
        if l is not None and 0 <= l <= 1 and 0 <= m <= 1:
            return l, m

    # If all attempts fail, use distance-based interpolation as fallback
    return distance_based_coord(vertices_x, vertices_y, point)

def try_calculate_coord(vertices_x, vertices_y, point, AI, tol=1e-10):
    """
    Try to calculate normalized coordinates with error checking
    """
    try:
        a = np.dot(AI, vertices_x)
        b = np.dot(AI, vertices_y)

        aa = a[3]*b[2] - a[2]*b[3]
        
        # Check if quadratic coefficient is too small
        if abs(aa) < tol:
            return None, None

        bb = a[3]*b[0] - a[0]*b[3] + a[1]*b[2] - a[2]*b[1] + point.x*b[3] - point.y*a[3]
        cc = a[1]*b[0] - a[0]*b[1] + point.x*b[1] - point.y*a[1]

        # Check discriminant
        disc = bb*bb - 4*aa*cc
        if disc < 0:
            return None, None

        # Try both roots
        det = np.sqrt(disc)
        for sign in [-1, 1]:
            m = (-bb + sign*det)/(2*aa)
            
            # Check denominator for l calculation
            denom = (a[1] + a[3]*m)
            if abs(denom) < tol:
                continue
                
            l = (point.x - a[0] - a[2]*m)/denom
            
            if 0 <= l <= 1 and 0 <= m <= 1:
                return l, m

        return None, None

    except (RuntimeWarning, FloatingPointError):
        return None, None

def distance_based_coord(vertices_x, vertices_y, point):
    """
    Fallback method using distance-based interpolation
    """
    # Calculate distances to all vertices
    vertices = np.column_stack((vertices_x, vertices_y))
    point_coord = np.array([point.x, point.y])
    distances = np.linalg.norm(vertices - point_coord, axis=1)
    
    # Avoid division by zero for points exactly at vertices
    eps = 1e-10
    distances = np.maximum(distances, eps)
    
    # Calculate weights based on inverse distance
    weights = 1.0/distances
    weights = weights/np.sum(weights)
    
    # Calculate normalized coordinates using weighted average
    l = weights[1] + weights[2]  # Right side contribution
    m = weights[2] + weights[3]  # Top side contribution
    
    # Ensure coordinates are in [0,1]
    l = np.clip(l, 0, 1)
    m = np.clip(m, 0, 1)
    
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
    start_time = time.time()

    # geometric definitions
    r = 0.295275
    r2 = r**2
    Ly = Lz = 1.0

    # MFSim mesh definitions
    n_cells = 5120
    dy = dz = Ly/n_cells

    print(f"Starting interpolation for {n_cells}x{n_cells} grid...")
    
    # read cfd-post data
    print("Reading input data...")
    vertices, vertices_connect, vel = read_cfdpost_data(
        'export3.csv',
        1276,
        delimiter=',',
        skiprows=6
    )

    # center values at (L/2, L/2)
    vertices[:,0] -= (3.725185 - Ly/2)
    vertices[:,1] += Lz/2

    # writing at the centroid of the cell
    values = np.empty((n_cells*n_cells,5))
    cell_range = range(1,n_cells+1)
    i = 0
    
    total_cells = n_cells * n_cells
    last_percent = -1
    last_time_update = time.time()
    update_interval = 5  # segundos
    
    print(f"\nStarting interpolation of {total_cells} points...")
    print("Progress: [", end='', flush=True)
    
    for j in cell_range:
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
            
            # Update progress
            current_time = time.time()
            if current_time - last_time_update >= update_interval:
                percent_complete = (i / total_cells) * 100
                elapsed_time = current_time - start_time
                estimated_total = elapsed_time / (i / total_cells)
                remaining_time = estimated_total - elapsed_time
                
                print(f"\rProgress: {percent_complete:.1f}% complete. "
                      f"Elapsed: {elapsed_time/60:.1f} min. "
                      f"Remaining: {remaining_time/60:.1f} min. "
                      f"Points processed: {i}/{total_cells}", end='', flush=True)
                
                last_time_update = current_time

    end_time = time.time()
    total_time = end_time - start_time
    print(f"\n\nInterpolation completed in {total_time/60:.1f} minutes")
    
    print("Saving results...")
    np.savetxt("teste.txt", values)
    print("Done!")

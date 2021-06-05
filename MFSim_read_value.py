import numpy as np

n_cells_lbot = 40
lvl_max = 8
n_cells = n_cells_lbot*(2**(lvl_max-1))

r = 0.64140/2
r2 = r**2
Ly = Lz = 1

dy = dz = Ly/n_cells

y = 0.5
z = 0.5

# Do arquivo, lê a primeira linha que diz qual o nível?

# if not dist > r:

j = round(0.5 + y/dy)
k = round(0.5 + z/dy)

linha = (j-1)*n_cells + k
print(linha)

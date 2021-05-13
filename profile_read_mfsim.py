import numpy as np

### SIMULATING MFSIM BEHAVIOR AND FUNCTIONS

### read_externalrij(headlev,l1,l2,       -> ins
#                   r11_in,r22_in,r33_in) -> outs

# l1 = 1 // l2 = ltop

r11_in = 0

# my_file=open("profile.in","r",encoding='utf-8')

# nlr11 = my_file.readline()  # number of lines

#allocate
# r11ex = np.empty((nlr11,3))

r11ex = np.loadtxt('profile.in',
                    skiprows=1,
                    delimiter=' ')


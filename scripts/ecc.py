import numpy as np
import sys

data = [float(x) for x in sys.argv[1:]]
tensor = np.empty((3,3))

# OVITO Format: xx, yy, zz, xy, xz, yz

tensor[0,0] = data[0] # xx
tensor[1,1] = data[1] # yy
tensor[2,2] = data[2] # zz
tensor[0,1] = tensor[1,0] = data[3] # xy
tensor[0,2] = tensor[2,0] = data[4] # xz
tensor[1,2] = tensor[2,1] = data[5] # yz


i = np.sort(np.linalg.eig(tensor)[0])
i_max = max(i)

e = np.empty(3)
for j in range(3):
    e[j] = np.sqrt(1-(i[j]/i_max)**2)
e_avg = np.round(np.sqrt(e[0]**2 + e[1]**2 + e[2]**2), 3)
print(e_avg)

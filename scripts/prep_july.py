import numpy as np
import pandas as pd
from copy import *
import math
import sys
import argparse
import os

#Structure: [["ResNumber","Res" "Atom Name", "N", x, y, z], ["ResNumber","Res" "Atom Name", "N", x, y, z]]


parser = argparse.ArgumentParser(description='Determine the degree of protonation of -sel atoms.')
parser.add_argument('-L', dest='chain_length', type=int, default=68, help='PEO Chain Length')
parser.add_argument('-M', dest='mols', type=int, default=1, help='Number of Molecules')
parser.add_argument('-N', dest='beads', type=int, help='Number of Calculated Beads')
parser.add_argument('-T', dest='tries', type=int, help='Number of tries per bead')
args = parser.parse_args()


x,y,z = [0,0,0]
box = [10,10,10]

def CollisionCheck(x,y,z,M,T):
    for m in M:
        if abs(m[0]-x) <= T and abs(m[1]-y) <= T and abs(m[2]-z) <= T:
            out=True
        else:
            out=False
    return out
def BoxCheck(x,y,z,B):
    if float(x) < 0 or float(y) < 0 or float(z) < 0:
        return True
    if float(x) > B[0] and float(y) > B[1] and float(z) > B[2]:
        return True
    return False

def UpdateCoords():
    global x,y,z
    x = np.round(x+0.2*np.random.uniform(-1, 1), 3)
    y = np.round(y+0.2*np.random.uniform(-1, 1), 3)
    z = np.round(z+0.2*np.random.uniform(-1, 1), 3)
    return x,y,z

positions = []
positions.append([x,y,z])
# V is number of beads in the system
available=0
for i in range(args.beads):
    x_old, y_old, z_old = deepcopy(x), deepcopy(y), deepcopy(z)
    UpdateCoords()
    comp = np.round((available/((args.chain_length+48)*args.mols)*100), 1)
    os.system('cls' if os.name == 'nt' else 'clear')
    print("💥 Collsion Check ✅")
    print("📦 Box Check ✅")
    print(f"{comp} %")
    print(f"{np.round((i/args.beads*100), 1)} %")
    if comp > 100:
        break
    else:
        for j in range(args.tries):
            if CollisionCheck(x,y,z,positions,0.5) == True:
                UpdateCoords()
            elif CollisionCheck(x,y,z,positions,0.5) == False:
                if BoxCheck(x,y,z,box) == True:
                    UpdateCoords()
                elif BoxCheck(x,y,z,box) == False:
                    positions.append([x,y,z])
                    available += 1
                    break            
            elif j == (args.tries -1):
                print("⚠️")
            else:
                UpdateCoords()

gro = pd.DataFrame(columns=["ResNumber","Res", "Atom Name", "N", "x", "y", "z"])
phe = {"resname": "PHE", "V": 4, "Vnames": ["BB", "SC1", "SC2", "SC3"]}
his = {"resname": "HIS", "V": 4, "Vnames": ["BB", "SC1", "SC2", "SC3"]}
lnk = {"resname": "LNK", "V": 4, "Vnames": ["BB", "BB1", "BB2", "BB3"]}
peo = {"resname": "PEO", "V": 1, "Vnames": ["EO"]}

print(f"Available Coordinates: {available}. That corresponds to {available / args.beads * 100} %")
def AddLines(res, resnumber, offset):
    global gro
    for number in range(res["V"]):
        index = offset + number
        gro.loc[index] = [resnumber+1, res["resname"], res["Vnames"][number], index+1, positions[index][0], positions[index][1], positions[index][2]]


#phe, his, phe, his, phe, lnk, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, lnk, phe, his, phe, his, phe,phe, his, phe, his, phe, lnk, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, peo, lnk, phe, his, phe, his, phe

# Setup for Input Structure
general = [phe, his, phe, his, phe]
link = [lnk]
#chain_length = 44
peo_chain = [peo] * args.chain_length
#mols = 5
structure = (general + link + peo_chain + link + general) * args.mols

i,p=0,0
for k in structure:
    try:
        AddLines(k, i, p)
    except:
        print("Too many restraints. Cannot generate coordinate file. ⛔️")
    i+=1
    p+=k["V"]


print(gro)
xmax =float(str(gro.nlargest(1, "x")["x"].tolist())[1:][:-1])
ymax =float(str(gro.nlargest(1, "y")["y"].tolist())[1:][:-1])
zmax =float(str(gro.nlargest(1, "z")["z"].tolist())[1:][:-1])
maxvals = [xmax, ymax, zmax]
temp=[]
for item in maxvals:
    temp.append(math.ceil(item))
boxsize = max(temp)

print(f'Max Values: x {xmax}, y {ymax}, z {zmax}')

column_positions = [5, 5, 5, 5, 8, 8, 8]  # Example: Columns will start at positions 0, 10, and 20
output_file = 'coordinates.gro'

with open(output_file, 'w') as f:
    for col, position in zip(gro.columns, column_positions):
        f.write(f'{col: <{position}}')  # Write column name at specified position
    f.write('\n')  # Move to the next line
    
    for _, row in gro.iterrows():
        for col, position in zip(gro.columns, column_positions):
            value = str(row[col])
            f.write(f'{value: <{position}}')  # Write value at specified position
        f.write('\n')  # Move to the next line

with open(output_file, 'r') as f:
    lines = f.readlines()
header_line="## Generated Coordinates (3D Random Walk without Overlap): projectraccoon.org ##"
lines = lines[1:]
with open(output_file, 'w') as f:
    f.writelines(lines)
with open(output_file, 'r') as f:
    content = f.read()
content = header_line + '\n' + str(p) + '\n' + content
new_content = content + str(boxsize) + " " + str(boxsize) + " " + str(boxsize)
with open(output_file, 'w') as f:
    f.writelines(new_content)


with open('system.top', 'w') as f:
    f.write(f'#include "martini3_MLO.itp" \n')
    f.write(f'#include "PPC.itp" \n')
    f.write(f'#ifdef POSRES \n#include "posre.itp" \n#endif \n')
    f.write(f'[ system ] \n')
    f.write(f'PPC \n')
    f.write(f'[ molecules ] \n')
    f.write(f'PPC {args.mols} \n')

with open('box.dimensions', 'w') as b:
    b.write(f'{str(boxsize)}')

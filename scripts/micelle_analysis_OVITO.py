import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
import argparse


parser = argparse.ArgumentParser(description='OVITO Imort')
parser.add_argument('-f', dest='fileformat', type=str, default='.csv', help='.csv/.txt')
parser.add_argument('-sep', dest='sep', type=str, default=' ', help='Seperator')
parser.add_argument('-dir', dest='directory', type=str, default=os.getcwd(), help='Directory')
parser.add_argument('-bins', dest='bins', type=int, default=25, help='Bins')
parser.add_argument('-out', dest='output', type=str, default='pH_x_', help='output filename')


args = parser.parse_args()


# Import Ovito Cluster Analysis Files

li = []

directory = args.directory
for dirpath, dirnames, filenames in os.walk(directory):
	for filename in filenames:
		if filename.endswith(args.fileformat):
			df = pd.read_csv(filename, sep=args.sep, header=[1], skipinitialspace=True)
			df = df.shift(periods=1, axis="columns")
			df = df.drop('#', axis=1)
			li.append(df)


df = pd.concat(li, axis=0, ignore_index=True)
rg = df.iloc[:,5].tolist()

# Error Correction for RG Calculation

rg_np = np.empty((len(rg), 1))
for i, radius in enumerate(rg):
	if (radius > 10) and (radius < 100):
		rg_np[i] = radius
	else:
		rg_np[i] = np.nan
		
rgmean = np.nanmean(rg_np, axis=0)
rgstd = np.nanstd(rg_np, axis=0)
rg = pd.DataFrame(rg_np)
rg.to_csv(f'{args.output}RG.txt')
plt.hist(rg_np, bins=args.bins)
plt.savefig(f'{args.output}RG.png')
def GetE(i, df):
	tensor = np.empty((3,3))
	tensor[0,0] = df.iloc[i,6] # xx
	tensor[1,1] = df.iloc[i,7] # yy
	tensor[2,2] = df.iloc[i,8] # zz
	tensor[0,1] = tensor[1,0] = df.iloc[i,9] # xy
	tensor[0,2] = tensor[2,0] = df.iloc[i,10] # xz
	tensor[1,2] = tensor[2,1] = df.iloc[i,11] # yz

	i = np.sort(np.linalg.eig(tensor)[0])
	i_max = max(i)

	e = np.empty(3)
	for j in range(3):
	    e[j] = np.sqrt(1-(i[j]/i_max)**2)
	e_avg = np.round(np.sqrt(e[0]**2 + e[1]**2 + e[2]**2), 3)
	return e_avg

e = []
for i, j in enumerate(df.index):
	e.append(GetE(i, df))
emean = np.nanmean(e, axis=0)
estd = np.nanstd(e, axis=0)
output_dic = {'RG': rgmean, 'D RG': rgstd, 'E': emean, 'D E': estd}
output = pd.DataFrame(output_dic)
output.to_csv(f'{args.output}results.txt')
	

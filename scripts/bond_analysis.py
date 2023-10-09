import MDAnalysis
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from tqdm import tqdm
import pandas as pd
import argparse

# Moritz Lennart Obenauer, August 2023, JGU Mainz

u = MDAnalysis.Universe('md_0_1_noPBC.gro', 'md_0_1_noPBC.xtc')

hydrogens1 = u.select_atoms('resname PHE and name H')
oxygens1 = u.select_atoms('resname HISD and name O')

hydrogens2 = u.select_atoms('resname HISD and name H')
oxygens2 = u.select_atoms('resname PHE and name O')

chainA = u.select_atoms('id 1-118')
chainB = u.select_atoms('id 119-236')


distance, deltadistance = [], []
chaindistance, dchaindistance = [], []

def UnitVector(v):
	return v / np.linalg.norm(v)

def GetDistances(a, b):
	temp = []
	for h in a:
		for o in b:
			temp.append(np.linalg.norm(h.position - o.position))
	avg = np.average(temp)
	std = np.std(temp)
	return avg, std

def DistanceChains(a,b):
	temp=[]
	for h in a:
		for o in b:
			temp.append(np.linalg.norm(a.center_of_mass() - b.center_of_mass()))
			k = UnitVector(a.center_of_mass())
			j = UnitVector(b.center_of_mass())
			dot = np.dot(k, j)
	avg = np.average(temp)
	std = np.std(temp)
	return avg, std, dot

def EtED(chain):
	temp=[]
	temp.append(np.linalg.norm(chain[0].position - chain[-1].position))
	avg = np.average(temp)
	std = np.std(temp)
	return avg, std

def DistanceAtoms(h, o):
	temp = []
	avg = np.linalg.norm(h.position - o.position)
	return avg
	


ete1 = []
ete2 = []
x = []
dot = []
time = []


for ts in tqdm(u.trajectory[::1]):
	#x.append(DistanceAtoms(chainA[0], chainB[-1]))
	distance.append(GetDistances(hydrogens1, oxygens1)[0])
	deltadistance.append(GetDistances(hydrogens1, oxygens1)[1])
	ete2.append(GetDistances(hydrogens2, oxygens2)[0])
	#deltadistance.append(GetDistances(hydrogens, oxygens)[1])
	ete1.append(DistanceChains(chainA.select_atoms("resname GLY"), chainB.select_atoms("resname ACE"))[0])
	dot.append(DistanceChains(chainA.select_atoms("resname GLY"), chainB.select_atoms("resname ACE"))[2])
	#ete2.append(DistanceChains(chainA.select_atoms("resname ACE"), chainB.select_atoms("resname GLY"))[0])	
	#dchaindistance.append(DistanceChains(chainA, chainB)[1])
	#ete1.append(EtED(chainA)[0])
	time.append(ts.frame / 100)
	#print(f'Time: {ts.frame/100} ns')

data = {'Time': time, 'HO-Distance 1': distance, 'HO-Distance 2': ete2, 'Dot Product': dot, 'Gly-Ace-Distance': ete1}
df = pd.DataFrame(data)
df.to_csv('results.csv')
print(df)

#plt.plot(x, label="Test")
plt.plot(time, ete1, label="Gly-Ace-Distance")
plt.plot(time, dot, label="Dot Product")
plt.errorbar(time, ete2, yerr=None, label="O-H-Inverse")
plt.errorbar(time, distance, yerr=None, label="O-H-Distance")
#plt.plot(chaindistance, label="Chain-Chain-Distance")
plt.legend()
plt.show()

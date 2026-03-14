from timeit import default_timer
import cantera as ct
import numpy as np

import matplotlib.pyplot as plt

gas = ct.Solution('gri30.xml')
input_file = 'gri30.yaml'
all_species = ct.Species.listFromFile(input_file)
species = []
for S in all_species:
    species.append(S)
species_names = {S.name for S in species}
print('Species: {0}'.format(', '.join(S.name for S in species)))
# standard atmosphere: 78.08 N2, 20.95 oxygen 0.93 argon 0.04 carbon dioxide
gas()
T = 500.0
P = 1.0 * ct.one_atm
#ch4=0.044
#ch4=0.022
ch4=0.011
o2=0.2095
h2o=0.015
ar=0.008056482
n2=0.7808
sum=ch4+o2+h2o+ar+n2
ch4=ch4/sum
o2=o2/sum
h2o=h2o/sum
ar=ar/sum
n2=n2/sum
X=f'CH4:{ch4:.6f}, O2:{o2:.6f}, H2O:{h2o:.5f}, N2:{n2:.5f}, Ar:{ar:.5f}'
print(X)



initial_state = T, P, X

species = gas.species
T=100.0
npoints = 20
#print(gas.n_species)
   
xeq = np.zeros((gas.n_species*npoints))
xH2 = np.zeros(npoints)
xCH4 = np.zeros(npoints)
Teq = np.zeros(npoints)
ii=0
for i in range(npoints):
    gas.TPX = T, P, X
    gas.equilibrate('TP',max_iter=2500)
    print(T,gas.X)
    T = T + 100.0

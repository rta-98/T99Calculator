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
#ch4=0.044
#ch4=0.022
ch4=0.06
#ch4=0.011
#ch4=0.0
o2=0.2095
h2o=0.2
#h2o=0.2095
#h2o=0.0212
#h2o=0.001
#h2o=0.0
ar=0.008056482
n2=0.7808
sum=ch4+o2+h2o+ar+n2
ch4=ch4/sum
o2=o2/sum
h2o=h2o/sum
ar=ar/sum
n2=n2/sum
X=f'CH4:{ch4:.6f}, O2:{o2:.6f}, H2O:{h2o:.5f}, N2:{n2:.5f}, Ar:{ar:.5f}'
#X=f'CH4:0.090743709, O2:0.181487417, H2O:0.043314419, N2:0.676397973, Ar:0.008056482'
print(X)

#T = 1000.0
T = 1000.0+273.15
P = ct.one_atm

#gas.TP=1000, ct.one_atm



initial_state = T, P, X

species = gas.species

npoints = 1
#print(gas.n_species)
   
xeq = np.zeros((gas.n_species*npoints))
xH2 = np.zeros(npoints)
xCH4 = np.zeros(npoints)
Teq = np.zeros(npoints)
ii=0
print("O2 OH H2O CH4")

for i in range(npoints):
    gas.TPX = T, P, X
    #for j in range(gas.n_reactions):
    #print(f"Reaction {j}; {gas.reaction_equation(j)}")
       
    reactor = ct.Reactor(gas)
    sim = ct.ReactorNet([reactor])
    time = 0.0
    for step in range(200):
        time += 1e-2 # increment time
        sim.advance(time)
        #print(f"Time: {time:.3e} s, Temperature: {reactor.T:.3f} K")
        print(time,T,gas.X[3],gas.X[4],gas.X[5],gas.X[13])

    #gas.equilibrate('TP',max_iter=2500)
    #print(T,gas.X[3],gas.X[4],gas.X[5],gas.X[13])
    #xH2[i]=gas.X[0]
    #xCH4[i]=gas.X[13]
    #for j in range(gas.n_species):
     #   xeq[ii+j] = gas.X[j]
    #Teq[i] = gas.T
    T = T + 50.0
    #ii=ii+gas.n_species
    #gas()
#gas()
#print(Teq[10])
#print(Teq)
#print(xH2)
#H2=0.0019872*Teq[10]*np.log(xH2[10])
#CH4=0.0019872*Teq[10]*np.log(xCH4[10])
#print(Teq[10],xH2[10],H2,xCH4[10],CH4)

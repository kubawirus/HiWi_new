"""
Constant-pressure, adiabatic kinetics simulation.
Requires: cantera >= 2.5.0, matplotlib >= 2.0
Keywords: combustion, reactor network, plotting
"""
import os
import datetime
import sys
from tkinter import Y
import matplotlib.pyplot as plt
import cantera as ct
import numpy as np
import pandas as pd

# INPUT
# plot = str(input('Plot biggest components? (y/n)\n'))
# n_plots = int(input('How many components to plot (int) \n'))

plot = 'y'
n_plots = 10

# How many time steps? delta t is 1e-5!
steps = 1e3

# Get the current working directory
cwd = os.getcwd()

# Get daytime
day_id = datetime.datetime.now().strftime("%Y.%m.%d_%H.%M.%S")
output = r"BatchTest_" + day_id + ".csv"
print('Output file: ' + output + ' in: ' + cwd)

# Mechanism?
mech13 = 'mech_13.yaml'
gas = ct.Solution(mech13)

# Other reaction for debug
# gas = ct.Solution('h2o2.yaml')
# gas.TPX = 1001.0, ct.one_atm, 'H2:2,O2:1,N2:4'

# Initial State?
#gas.TPX = 1000, ct.one_atm, 'C3H8:10, H2:1'
gas.TPX = 873.0, ct.one_atm, 'C3H8:5, C3H6:3, H2:2'
r = ct.IdealGasConstPressureReactor(gas)

# Reaktor
sim = ct.ReactorNet([r])
sim.verbose = True

# Time step and number of steps
dt_max = 1.e-5
# steps = 1e3
#steps = 5*1e5
t_end = steps * dt_max
# states = ct.SolutionArray(gas, extra=['t'])

# Create ndarrays for results
X = np.array([r.thermo.X])
Names =np.array(r.thermo.species_names)
M = r.thermo.molecular_weights
# Sort acc to Moll Mass
sorted_to = M.argsort()
M = M[sorted_to[::-1]]
Names = Names[sorted_to[::-1]]

sim_steps = 0
# Simulation loop
while sim.time < t_end:
    sim.advance(sim.time + dt_max)
    x = r.thermo.X
    X = np.vstack([X,x])
    sim_steps += 1

# # Debug prints
# print(Names)
# print(X[0])
# sort X acc to Moll Mass
for i in range(int(sim_steps) + 1):
    X[i] = X[i][sorted_to[::-1]]
# print(X[0])
# Change ndarray type to str
M_str = M.astype('U')
# Create index array for dataframe
columnI = [Names[i] + ' ' + M_str[i] for i in range(len(Names))]
# Names = np.insert(Names,0, 'time')
times = np.arange(0,int(sim_steps) + 1)
df_result = pd.DataFrame(data=X,
             index=times,
             columns=columnI)
df_result.to_csv(cwd + '/' + output)

Id_plot=[]

# Plot biggest compounds that are present in system
for i in range(len(Names)):
    for j in range(int(sim_steps)):
        if X[j][i] > 10e-8:
            if i not in Id_plot:
                Id_plot.append(i)

fig, ax = plt.subplots(figsize=(8,5))
for i in (Id_plot[:n_plots]):
    ax.plot(times[:-1], df_result[(Names[i] + ' ' + M_str[i])][:-1], label=(Names[i]+' '+M_str[i]))
fig.legend()
fig.show()

fig.savefig(cwd + '/figure_BatchTest_' + day_id + '.png', dpi=400)

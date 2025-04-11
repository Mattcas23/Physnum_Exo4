import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = r"/Users/a-x-3/Desktop/Exercice4_2025_student/exe"  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/a-x-3/Desktop/Exercice4_2025_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Name of the input file


N1 = np.array([8,9,10,15,20,40,50,90])
N2 = N1

nsimul = N1.size

energy = np.zeros(nsimul) # added 
paramstr = 'N1'  # Parameter name to scan
param = N1  # Parameter values to scan
paramstr2 = 'N2'
param2 = N2

# Simulations
outputs = []  # List to store output file names
convergence_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}"
    outputs.append(output_file)
    #cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    cmd = f"{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    #cmd = f"{executable}"
    print(cmd)
    print(output_file)
    subprocess.run(cmd, shell=True)
    print('Done.')


phi0 = []

for i in range(nsimul):  # Iterate through the results of all simulations

    data = np.loadtxt(outputs[i]+"_phi.out")  # Load the output file of the i-th simulation
    phi0.append(data[0,1])


lw = 1.5
fs = 16

E = np.loadtxt(outputs[-1]+"_E.out")
D = np.loadtxt(outputs[-1]+"_D.out")
Phi = np.loadtxt(outputs[-1]+"_phi.out")

def Conv ( norder = 2 ) : # étude de convergence

    plt.figure()
    plt.plot(1/pow(N1*2,norder), phi0 , 'k+-')
    plt.xlabel('(1/N)' + f"$^{norder}$", fontsize = fs)
    plt.ylabel('$\\phi_0$', fontsize = fs)
    #plt.legend()

    

def Eplot () : # Plot le champ électrique en fonction de r 

    plt.figure()
    plt.plot(E[:,0],E[:,1], 'k+-')
    plt.xlabel('r', fontsize = fs)
    plt.ylabel('E', fontsize = fs)
    #plt.legend()

def Dplot () : # Plot le champ de déplacement en fonction de r 

    plt.figure()
    plt.plot(D[:,0],D[:,1], 'k+-')
    plt.xlabel('r', fontsize = fs)
    plt.ylabel('D', fontsize = fs)
    #plt.legend()
    
def Phiplot () : # Plot le potentiel en fonction de r 

    plt.figure()
    plt.plot(Phi[:,0],Phi[:,1], 'k+-')
    plt.xlabel('r', fontsize = fs)
    plt.ylabel('$\\phi(r)$', fontsize = fs)
    #plt.legend()

def Distance () : # Contrôle pour les midPoints

    plt.figure()
    plt.scatter(Phi[:,0],np.ones(Phi[:,0].size),label = "$r_i$")
    plt.scatter(D[:,0],np.ones(D[:,0].size),label = "$r_{i+1/2}$" , marker='+')
    plt.legend()

Conv()
Phiplot()
Dplot()
Eplot()
#Distance()


plt.show()

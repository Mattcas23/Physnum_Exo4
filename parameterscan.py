import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = r"/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/exe" #r"/Users/a-x-3/Desktop/Ex3_2024_student/exe"  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4"#r"/Users/a-x-3/Desktop/Ex3_2024_student"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Name of the input file

#----------------------------------------- Valeurs du Configfile --------------------------------- # 

#Values = np.genfromtxt("/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/buildconfiguration.in.example" , comments = '//')

#R = Values[0,-1]
#r1 = Values[1,-1]
#epsilon_a = Values[2,-1]
#epsilon_b= Values[3,-1]
#uniform_rho_case = Values[4,-1]
#VR = Values[5,-1]
#rho0  =  Values[6,-1]
#N1    = Values[7,-1]
#N2  = Values[8,-1]
#maxit  = Values[9,-1]
#output = Values[10,-1]
#sampling = Values[11,-1]
#verbose = Values[12,-1]
# ---------------------------------------------------------------

# ------------------------------------------- Simulations ---------------------------------------------

data = np.loadtxt("/Users/Matte/Desktop/EPFL/BA4/Physnum/Physnum_Exo4/phioutput.out")  # Load the output file of the i-th simulation
t = data[:, 0]
phi = data[:,1] # nombre de pas de temps total pour la simulation
     
print(phi)


#conv()
plt.show()

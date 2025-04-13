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


############################ Valeurs du configfile ############################

Values = np.genfromtxt("configuration.in.example" , comments = '//')

R  = Values[0,-1]
print(R)
r1 = Values[1,-1] 
print(r1)
uniform_rho_case = Values[5,-1]
print(uniform_rho_case)
rho0 = Values[7,-1] # valeur pour la question a) 
print(rho0)
alpha = Values[13,-1] # coeff de proportionnalité
print(alpha)


########################### Simulations ##############################


#N1 = np.array([90,110,120,140,150,170,200])
#N2 = N1 * alpha

N1 = np.array([1500])
N2 = N1*int(alpha)
 
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

################################## Plots #################################

phi0 = []
phir1 = []

for i in range(nsimul):  # Iterate through the results of all simulations

    data = np.loadtxt(outputs[i]+"_phi.out")  # Load the output file of the i-th simulation
    phi0.append(data[0,1])
    phir1.append(data[N1[i],1])

lw = 1.5
fs = 16

E = np.loadtxt(outputs[-1]+"_E.out")
D = np.loadtxt(outputs[-1]+"_D.out")
Phi = np.loadtxt(outputs[-1]+"_phi.out")

def Conv_phi0 ( norder = 3 ) : # étude de convergence ( question a) ordre 2 )

    plt.figure()
    plt.plot(1/pow(N1+N2,norder), phi0 , 'k+-')
    plt.xlabel('(1/N)' + f"$^{norder}$", fontsize = fs)
    plt.ylabel('$\\phi_0$', fontsize = fs)
    #plt.legend()

def Conv_phir1 ( norder = 2 ) : # ( question b) ordre 2 )

    plt.figure()
    plt.plot(1/pow(N1+N2,norder), phir1 , 'k+-')
    plt.xlabel('(1/N)' + f"$^{norder}$", fontsize = fs)
    plt.ylabel('$\\phi(r_1)$', fontsize = fs)  
    

def Eplot () : # Plot le champ électrique en fonction de r

    E0_eq = ''

    plt.figure()
    plt.title(f'N1 = {N1[-1]} et N2 = {N2[-1]}', fontsize = fs - 2)
    plt.plot(E[:,0],E[:,1], color = "black" )

    if uniform_rho_case : # question a) (rho unif) : on affiche la solution analytique 

        print("A faire")#plt.plot(E[:,0], -(Phi[:,0]**2 - R**2)/(2), color = 'orange', linestyle = 'dashed', label = E0_eq)

    else : # sinon on affiche la zone de changement de permitivité 

        plt.vlines(r1, ymin = min(E[:,1]) , ymax = max(E[:,1]) , color = 'red' , linestyle = 'dashed' , label = f"$ r_1 = {r1} $")
        
    plt.xlabel('r [m]', fontsize = fs)
    plt.ylabel('E [V/m]', fontsize = fs)
    plt.legend()

def Dplot () : # Plot le champ de déplacement en fonction de r

    D0_eq = ''

    plt.figure()
    plt.title(f'N1 = {N1[-1]} et N2 = {N2[-1]}', fontsize = fs - 2)
    plt.plot(D[:,0],D[:,1], color = "black")

    if uniform_rho_case : # question a) (rho unif) : on affiche la solution analytique 

        print("A faire") #plt.plot(D[:,0], -(Phi[:,0]**2 - R**2)/(2), color = 'orange', linestyle = 'dashed', label = D0_eq)

    else : # sinon on affiche la zone de changement de permitivité 
        plt.vlines(r1, ymin = min(D[:,1]) , ymax = max(D[:,1]) , color = 'red' , linestyle = 'dashed' , label = f"$ r_1 = {r1} $")

    plt.xlabel('r [m]', fontsize = fs)
    plt.ylabel('D []', fontsize = fs)
    plt.legend()
    
def Phiplot () : # Plot le potentiel en fonction de r

    phi0_eq = '$\\phi(r) = -\\frac{\\rho_0}{2 \\epsilon_0} (r^2-R^2)$'
    
    plt.figure()
    plt.title(f'N1 = {N1[-1]} et N2 = {N2[-1]}', fontsize = fs - 2)
    plt.plot(Phi[:,0],Phi[:,1], color = "black")

    if uniform_rho_case : # question a) (rho unif) : on affiche la solution analytique 

        plt.plot(Phi[:,0], -(Phi[:,0]**2 - R**2)/(2), color = 'orange', linestyle = 'dashed', label = phi0_eq)

    else : # sinon on affiche la zone de changement de permitivité 

        plt.vlines(r1, ymin = min(Phi[:,1]) , ymax = max(Phi[:,1]) , color = 'red' , linestyle = 'dashed' , label = f"$ r_1 = {r1} $")

    plt.xlabel('r [m]', fontsize = fs)
    plt.ylabel('$\\phi(r)$ [V]', fontsize = fs)
    plt.legend()

def Distance () : # Contrôle pour les midPoints

    plt.figure()
    plt.title(f'N1 = {N1[-1]} et N2 = {N2[-1]}', fontsize = fs - 2)
    plt.scatter(Phi[:,0],np.ones(Phi[:,0].size),label = "$r_i$")
    plt.scatter(D[:,0],np.ones(D[:,0].size),label = "$r_{i+1/2}$" , marker='+')
    plt.legend()

#Conv_phi0()
#Conv_phir1()
Phiplot()
Dplot()
Eplot()
#Distance()


plt.show()

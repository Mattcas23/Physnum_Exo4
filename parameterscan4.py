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
print("R : " , R)
r1 = Values[1,-1] 
print("r1 : " ,r1)
uniform_rho_case = Values[5,-1]
print("unirhocase : " , uniform_rho_case)
rho0 = Values[7,-1] # valeur pour la question a) 
print("rho0: " , rho0)



########################### Simulations ##############################


#N1 = np.array([90,110,120,140,150,170,200]) # liste utilisée pour la convergence a)


N1 = np.array([2,5,10,20,50,100])
#N2 = N1*int(alpha)
N2 = np.array([1,2,10,20,50,100]) 

if N1.size != N2.size :
    
    raise ValueError("La taille de N1 doit être la même que celle de N2") 

nsimul = N1.size

energy = np.zeros(nsimul) # added 
paramstr = 'N1'  # Parameter name to scan
param = N1  # Parameter values to scan
paramstr2 = 'N2'
param2 = N2

# Simulations avec N1 et N2
outputs = []  # List to store output file names
convergence_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}_{paramstr2}={param2[i]}"
    outputs.append(output_file)
    #cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    cmd = f"{executable} {input_filename} {paramstr}={param[i]:.15g} {paramstr2}={param2[i]:.15g} output={output_file}"
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

def Conv_phi0 ( norder = 2 ) : # étude de convergence ( question a) ordre 2 )

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
    

def Eplot (multiple) : # Plot le champ électrique en fonction de r

    E0_eq = '$ E(r) = \\rho_0 \\frac{r}{ \\epsilon_0} $'

    plt.figure()
    plt.title(f'N1 = {N1[-1]} et N2 = {N2[-1]}', fontsize = fs - 2)
    plt.plot(E[:,0],E[:,1], color = "black" )

    if uniform_rho_case : # question a) (rho unif) : on affiche la solution analytique 

        plt.plot(E[:,0], - E[:,0] , color = 'orange', linestyle = (0, (5, 7)) , label = E0_eq) # rho0 = eps0 

    else : # sinon on affiche la zone de changement de permitivité 

        plt.vlines(r1, ymin = min(E[:,1]) , ymax = max(E[:,1]) , color = 'red' , linestyle = 'dashed' , label = f"$ r_1 = {r1} $")
        
    plt.xlabel('r [m]', fontsize = fs)
    plt.ylabel('E [V/m]', fontsize = fs)
    plt.legend(fontsize = fs - 2.5)

    if multiple :

        plt.figure()
        for i in range(len(outputs)) :
            Ei = np.loadtxt(outputs[i]+"_E.out")
            plt.plot(Ei[:,0],Ei[:,1], label = f"N1 = {N1[i]} / N2 = {N2[i]}")
            
        plt.xlabel('r [m]', fontsize = fs)
        plt.ylabel('E [V/m]', fontsize = fs)
        plt.legend()    

def Dplot (multiple) : # Plot le champ de déplacement en fonction de r (multiple = true : met toutes les courbes pour diff N)

    D0_eq = '$ D(r) = r \\rho_0  $'

    plt.figure()
    plt.title(f'N1 = {N1[-1]} et N2 = {N2[-1]}', fontsize = fs - 2)
    plt.plot(D[:,0],D[:,1], color = "black")

    if uniform_rho_case : # question a) (rho unif) : on affiche la solution analytique 

        plt.plot(D[:,0], - rho0 * D[:,0] , color = 'orange', linestyle = (0, (5, 7)), label = D0_eq)

    else : # sinon on affiche la zone de changement de permitivité 
        plt.vlines(r1, ymin = min(D[:,1]) , ymax = max(D[:,1]) , color = 'red' , linestyle = 'dashed' , label = f"$ r_1 = {r1} $")

    plt.xlabel('r [m]', fontsize = fs)
    plt.ylabel('D [C$\,$m$^{-2}$]', fontsize = fs)
    plt.legend(fontsize = fs - 2.5)

    if multiple :

        plt.figure()
        for i in range(len(outputs)) :
            Di = np.loadtxt(outputs[i]+"_D.out")
            plt.plot(Di[:,0],Di[:,1], label = f"N1 = {N1[i]} / N2 = {N2[i]}")
            
        plt.xlabel('r [m]', fontsize = fs)
        plt.ylabel('D [C$\,$m$^{-2}$]', fontsize = fs)
        plt.legend()    

        
    
def Phiplot (multiple) : # Plot le potentiel en fonction de r

    phi0_eq = '$\\phi(r) = -\\frac{\\rho_0}{2 \\epsilon_0} (r^2-R^2)$'
    
    plt.figure()
    plt.title(f'N1 = {N1[-1]} et N2 = {N2[-1]}', fontsize = fs - 2)
    plt.plot(Phi[:,0],Phi[:,1], color = "black")

    if uniform_rho_case : # question a) (rho unif) : on affiche la solution analytique 

        plt.plot(Phi[:,0], -(Phi[:,0]**2 - R**2)/(2), color = 'orange', linestyle = (0, (5, 7)), label = phi0_eq)

    else : # sinon on affiche la zone de changement de permitivité 

        plt.vlines(r1, ymin = min(Phi[:,1]) , ymax = max(Phi[:,1]) , color = 'red' , linestyle = 'dashed' , label = f"$ r_1 = {r1} $")

    plt.xlabel('r [m]', fontsize = fs)
    plt.ylabel('$\\phi(r)$ [V]', fontsize = fs)
    plt.legend(fontsize = fs - 2.5)

    if multiple :

        plt.figure()
        for i in range(len(outputs)) :
            Phii = np.loadtxt(outputs[i]+"_Phi.out")
            plt.plot(Phii[:,0],Phii[:,1], label = f"N1 = {N1[i]} / N2 = {N2[i]}")
            
        plt.xlabel('r [m]', fontsize = fs)
        plt.ylabel('$\\phi(r)$ [V]', fontsize = fs)
        plt.legend() 

def Distance () : # Contrôle pour les midPoints

    plt.figure()
    plt.title(f'N1 = {N1[-1]} et N2 = {N2[-1]}', fontsize = fs - 2)
    plt.scatter(Phi[:,0],np.ones(Phi[:,0].size),label = "$r_i$")
    plt.scatter(D[:,0],np.ones(D[:,0].size),label = "$r_{i+1/2}$" , marker='+')
    plt.legend(fontsize = fs - 2.5)

#Conv_phi0()
#Conv_phir1()
Phiplot(True)
Dplot(True)
Eplot(True)
#Distance()


plt.show()

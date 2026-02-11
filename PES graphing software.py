import numpy as np
import re
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
import os

plt.style.use('_mpl-gallery')


def PES_Filereader(Folder): 


    global r_list

    all_files = os.listdir(Folder)

    r_list = []
    theta_list = []

            
    """Now that we have the list of filenames. We can iteratively open each file and read the contents. 
    Using regular expressions, we search for the line with the energyand extract the numerical value of the energy. We then close the file and move on to the next. 
    At the end of the process, a list of energies is output in same order as the files in the directory.  """

    Energies = []
    for name in all_files: 

        if not name.endswith('.out') or name.startswith('.'):
            continue

        r_value = re.search(r'\.r([-+]?\d*\.?\d+)theta', name)
        r_list.append(float(r_value.group(1)))
        

        theta_value = re.search(r'theta([-+]?\d*\.?\d+)\.out', name)
        theta_list.append(float(theta_value.group(1)))
        
        Gaussian_Output = open(f"{Folder}/{name}", 'r')

        #This reads each file but skips up to line 166 
        for line in Gaussian_Output.readlines()[166:]: 
            match = re.search(r"E\(RHF\)\s*=?\s*([-+]?\d+\.\d+)", line)
            if match: 
                Energies.append(float(match.group(1)))
                break

        Gaussian_Output.close()
    
    #The below code makes sure all of the variables exit the function as arrays. 

    Energies = np.array(Energies) 
    theta_array = np.array(theta_list)
    r_array = np.array(r_list)
    return(Energies, r_array, theta_array)

#Comments!
def PES_landscaper(Energy, r, theta):

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'}, figsize = (10,7))
    surf = ax.plot_trisurf(r, theta, Energy, cmap=cm.viridis, edgecolor='none', antialiased=True)

    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

    ax.set_xlabel("Bond Length (r) / Å")
    ax.set_ylabel("Bond Angle (Theta) / °")
    ax.set_zlabel("Energies / Hartrees")

    plt.show()


def PES_Vibrational_Freq(Energy, r, theta): 

    #Extracts the unique values of r/theta from the list of every value of r/theta that filereader provides
    unique_r = np.unique(r)
    unique_theta = np.unique(theta)

    #Reduced masses (hardcoded)
    red_mass_1 = 2*1.66*(10**(-27))
    red_mass_2 = 0.5*1.66*(10**(-27))

    """Energies are grouped by theta and sorted, then sorted by r. This is then put through a standard "C" reshape.
    This ends up generating E_grid as a 2D array of """
    indices = np.lexsort((r, theta)) 
    sorted_Energy = Energy[indices]
    
    E_grid = sorted_Energy.reshape(( len(unique_theta),len(unique_r)))

    
    #This determines the equilibrium geometery by finding the minimum value of the energy
    Equi_row, Equi_coloumn = np.where(E_grid == E_grid.min())
    Equi_theta = Equi_row[0]
    Equi_r = Equi_coloumn[0]

    """This is a set of hardcoded conversion factors to convert the elements of our hessian matrix to SI units"""

    Hartree_to_J = 4.3597447222071e-18
    Angstrom_to_m = 1e-10
    Degree_to_rad = np.pi / 180
    r_equilibrium = unique_r[Equi_r]*Angstrom_to_m

    """Finding the elements of the hessian matrix finding the gradient along each axis (corresponding to theta and r).
    Then taking those two matrixes and taking the derivative again, this time at the equilibrium geometry"""
    
    theta_grad, r_grad  = np.gradient(E_grid, 1, 0.05)

    Second_deriv_of_r = np.gradient(r_grad, 0.05, axis = 1)[Equi_theta, Equi_r]
    Second_deriv_of_theta = np.gradient(theta_grad, 1, axis = 0)[Equi_theta, Equi_r]
    drdtheta = np.gradient(r_grad, 1 , axis = 0)[Equi_theta, Equi_r]


    Second_deriv_of_r = Second_deriv_of_r * Hartree_to_J / ((Angstrom_to_m**2)*red_mass_1)
    Second_deriv_of_theta = Second_deriv_of_theta * Hartree_to_J / ((Degree_to_rad**2)*red_mass_2*r_equilibrium**2)
    drdtheta = drdtheta * Hartree_to_J / (Degree_to_rad * Angstrom_to_m* np.sqrt(red_mass_2*r_equilibrium**2*red_mass_1))

    Hessian_matrix = np.array([[Second_deriv_of_r, drdtheta], 
                               [drdtheta, Second_deriv_of_theta]])
    
    Ei_vals, Eigenvectors = LA.eig(Hessian_matrix)
 
    Ei_vals.sort()
    K_theta = Ei_vals[0]
    K_r = Ei_vals[1]

    

    Frequency_symmetric = ((1/(2*np.pi))*np.sqrt(K_r))/(3*10**10)
    Frequency_bend = ((1/(2*np.pi))*np.sqrt(K_theta))/(2.9979*10**10)

    print(Frequency_symmetric)
    print(Frequency_bend)



if __name__ == "__main__":
    UserFolder = input("What is the directories filepath?")

    Energies_in, r_in, theta_in = PES_Filereader(UserFolder)

    Type_of_query = input("""Which function should I perform?
        1. Calculate and graph the PES
        2. Calculate the vibrational frequencies
        3. Exit""")

    if Type_of_query == "1": 
        PES_landscaper(Energies_in, r_in,  theta_in)
        exit()
    
    elif Type_of_query == "2": 
        PES_Vibrational_Freq(Energies_in, r_in, theta_in)
        exit()

    elif Type_of_query == "3": 
        exit()
    
    else: 
        Type_of_query = input("""Which function should I perform? (Please use the indicies next to each option)
    1. Calculate and graph the PES
    2. Calculate the vibrational frequencies
    3. Exit""")


#Vib frequ  won't match exactly

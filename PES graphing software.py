import numpy as np
import re
import matplotlib.pyplot as plt
from numpy import linalg as LA
from matplotlib import cm
from scipy.optimize import curve_fit
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

    """Energies are grouped by theta and sorted, then sorted by r. This is then put through a standard "C" reshape.
    This ends up generating E_grid as a 2D array of """
    indices = np.lexsort((r, theta)) 
    sorted_Energy = Energy[indices]

    unique_r = np.unique(r)
    unique_theta = np.unique(theta)

    dr = unique_r[1] - unique_r[0]
    dtheta = unique_theta[1] - unique_theta[0]


    red_mass_r = 2*1.66*(10**(-27))
    red_mass_theta = 0.5*1.66*(10**(-27))

    #This creates a grid of the energies, with r on the x axis and theta on the y axis. It is necessary to find the hessian matrix. 
    E_grid = sorted_Energy.reshape((len(unique_theta),len(unique_r)))

    Equi_row, Equi_coloumn = np.where(E_grid == E_grid.min())

    """Finding the elements of the hessian matrix, we do this by taking the derivative with respect to theta and r.
    Then taking those two matrixes and taking the derivative again!
    We also use conversion factors to change the units from hartrees / angstrom^2 to J/m^2 """

    Hartree_to_J = 4.3597447222071e-18
    Angstrom_to_m = 1e-10
    Degree_to_rad = np.pi / 180
    r_equilibrium = unique_r[Equi_coloumn[0]] * Angstrom_to_m
    
    theta_grad, r_grad = np.gradient(E_grid, dtheta, dr)

    Second_deriv_of_r = np.gradient(r_grad, dr, axis = 1)[Equi_row, Equi_coloumn] * Hartree_to_J / (Angstrom_to_m**2)

    Second_deriv_of_theta = np.gradient(theta_grad, dtheta, axis = 0)[Equi_row, Equi_coloumn] * Hartree_to_J / ((Degree_to_rad**2)*red_mass_theta*r_equilibrium**2)

    drdtheta = np.gradient(r_grad, dtheta, axis = 0)[Equi_row, Equi_coloumn] * Hartree_to_J / (Degree_to_rad * Angstrom_to_m) *np.sqrt(red_mass_theta * r_equilibrium ** 2 * red_mass_r)

    Hessian_matrix = np.array([[Second_deriv_of_r, drdtheta],
    [drdtheta, Second_deriv_of_theta]])

    Hessian_matrix = Hessian_matrix.transpose(2, 0, 1)


    eigenvalues, Eigenvectors = LA.eigh(Hessian_matrix)

    kr = eigenvalues[:, 0]
    ktheta = eigenvalues[:, 1]

    Frequency_symmetric = ((1 / ( 2 * np.pi)) * np.sqrt(kr)) / (2.9979*10**10)
    Frequency_bend = ((1 / (2 * np.pi)) * np.sqrt(ktheta)) / (2.9979*10**10)

    print(Frequency_symmetric, Frequency_bend)

    



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


#Vib frequ  won't match exactly

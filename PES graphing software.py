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

    #Extracts the unique values of r/theta from the list of every value of r/theta that filereader provides
    unique_r = np.unique(r)
    unique_theta = np.unique(theta)

    #Extracts the average step size for theta and r
    theta_diffs = np.diff(np.sort(unique_theta))
    r_diffs = np.diff(np.sort(unique_r))

    theta_step = theta_diffs[0]
    r_step = r_diffs[0]

    #Reduced masses (hardcoded)
    red_mass_r = 2*1.66*(10**(-27))
    red_mass_theta = 0.5*1.66*(10**(-27))

    """Energies are grouped by theta and sorted, then sorted by r. This is then put through a standard "C" reshape.
    This ends up generating E_grid as a 2D array of energies with theta on x and r on y. 
    We also generate two 1D arrays to describe the r and theta at each point in the grid for the fit."""

    indices = np.lexsort((r, theta)) 
    sorted_Energy = Energy[indices]
    E_grid = sorted_Energy.reshape((len(unique_theta),len(unique_r)))
    R, THETA = np.meshgrid(unique_r, unique_theta)
   
    r_flat_array = R.flatten()
    theta_flat_array = THETA.flatten()
    coords_flat = np.vstack((r_flat_array, theta_flat_array))
    
    #Now we grab the equilibrium geometries

    Equi_row, Equi_coloumn = np.where(E_grid == E_grid.min())
    Equi_theta_idx = Equi_row[0]
    Equi_r_idx = Equi_coloumn[0]
    

    #Now we find guesses for K theta and K r using a rough mass weighted hessian
    theta_grad, r_grad = np.gradient(E_grid, unique_theta, unique_r)
    Second_deriv_of_r = np.gradient(r_grad, unique_r, axis=1)[Equi_theta_idx, Equi_r_idx]
    Second_deriv_of_theta = np.gradient(theta_grad, unique_theta, axis=0)[Equi_theta_idx, Equi_r_idx]
    drdtheta = np.gradient(r_grad, unique_theta, axis=0)[Equi_theta_idx, Equi_r_idx]

    Guessian_matrix = np.array([[Second_deriv_of_r, drdtheta],
                                  [drdtheta, Second_deriv_of_theta]])
    

    Ei_vals, Eigenvectors = LA.eigh(Guessian_matrix)

    Ei_vals.sort()

    K_theta_approx = Ei_vals[0]
    K_r_approx = Ei_vals[1]

    def PES_fit_func(coords, E_0, K_r, K_theta, r_equiv, theta_equiv):
        r_val, theta_val = coords
        return E_0 + 0.5 * K_r * (r_val - r_equiv)**2 + 0.5 * K_theta * (theta_val - theta_equiv)**2
    
    initial_guesses = [E_grid.min(), K_r_approx, K_theta_approx, unique_r[Equi_r_idx], unique_theta[Equi_theta_idx]]

    r_min = unique_r[Equi_r_idx] - 3 * r_step
    r_max = unique_r[Equi_r_idx] + 3 * r_step
    theta_min = unique_theta[Equi_theta_idx] - 10 * theta_step
    theta_max = unique_theta[Equi_theta_idx] + 10 * theta_step
    
    Masked_coords_flat = (coords_flat[0, :] > r_min) & \
                     (coords_flat[0, :] < r_max) & \
                     (coords_flat[1, :] > theta_min) & \
                     (coords_flat[1, :] < theta_max)


    popt, pcov = curve_fit(PES_fit_func, Masked_coords_flat, sorted_Energy, p0=initial_guesses)

    fitted_K_r, fitted_K_theta = popt[1], popt[2]

    Hartree_to_J = 4.35974e-18
    Angstrom_to_m = 1e-10

    k_r_SI = fitted_K_r * Hartree_to_J / (Angstrom_to_m**2)
    k_theta_SI = fitted_K_theta * Hartree_to_J * ((180 / np.pi)**2)

    r_eq_m = popt[3] * Angstrom_to_m 
    moment_of_inertia = red_mass_theta * (r_eq_m**2)

    omega_r = np.sqrt(k_r_SI / red_mass_r)
    omega_theta = np.sqrt(k_theta_SI / moment_of_inertia)

    c_cm = 2.9979e10

    Frequency_symmetric = (omega_r / (2 * np.pi)) / c_cm
    Frequency_bend = (omega_theta / (2 * np.pi)) / c_cm

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


#Vib frequ  won't match exactly

import numpy as np
import re
import matplotlib.pyplot as plt

from matplotlib import cm

plt.style.use('_mpl-gallery')


def PES_Filereader(Folder): 

    """The first step in reading these directories is generating the list of filenames. I have used loops to generate a list of filenames
    using hardcoded parameters dependening on the input directory."""

    filenames = [] 
    r_list = []
    theta_list = []

    if "H2Ooutfiles" in Folder:

        r_range = np.arange(0.70, 1.95, 0.05)
        molecule = "H2O"

    elif "H2Soutfiles" in Folder: 

        r_range = np.arange(0.70, 1.85, 0.05)
        molecule = "H2S"
    
    else: 
        print("Sorry, this isn't a directory I can parse. Make sure you've pasted the full file path.")
        exit()

    for r in r_range:
        for theta in range(70, 161, 1):
            filenames.append(f"{Folder}/{molecule}.r{r:.2f}theta{theta}.0.out")
            r_list.append(round(r, 5))
            theta_list.append(round(theta,5))
            
    """Now that we have the list of filenames. We can iteratively open each file and read the contents. 
    Using regular expressions, we search for the line with the energyand extract the numerical value of the energy. We then close the file and move on to the next. 
    At the end of the process, a list of energies is output in same order as the files in the directory.  """

    Energies = []
    for name in filenames: 
        Gaussian_Output = open(name, 'r')

        for line in Gaussian_Output: 
            match = re.search(r"E\(RHF\)\s*=?\s*([-+]?\d+\.\d+)", line)
            if match: 
                Energies.append(float(match.group(1)))
                break

        Gaussian_Output.close()
    
    #The below code makes sure all of the variables exit the function as arrays. This is probably not necessary for r and theta. 

    Energies = np.array(Energies) 
    theta_array = np.array(theta_list)
    r_array = np.array(r_list)

    return(Energies, r_array, theta_array, r_range)

#Comments!
def PES_landscaper(Energy, r, theta):

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'}, figsize = (10,7))
    surf = ax.plot_trisurf(r, theta, Energy, cmap=cm.viridis, edgecolor='none', antialiased=True)

    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

    ax.set_xlabel("Bond Length (r) / Å")
    ax.set_ylabel("Bond Angle (Theta) / °")
    ax.set_zlabel("Energies / Hartrees")

    plt.show()

#Comments!
def PES_Vibrational_Freq(Energy, r, theta, r_range): 

    #This creates a grid of the energies, with r on the x axis and theta on the y axis. It is necessary to find the hessian matrix. 
    E_grid = Energy.reshape(len(r_range), 91)
    Equi_row, Equi_coloumn = np.where(E_grid == E_grid.min())
    

    r_grad, theta_grad = np.gradient(E_grid, 0.05, 1)

    Second_deriv_of_r = np.gradient(r_grad, 0.05, axis = 0)
    Second_deriv_of_theta = np.gradient(theta_grad, 1, axis = 0)




if __name__ == "__main__":
    UserFolder = input("What is the folder name?")

    Energies_in, r_in, theta_in, range_in = PES_Filereader(UserFolder)

    Type_of_query = input("""Which function should I perform?
        1. Calculate and graph the PES
        2. Calculate the vibrational frequencies
        3. Exit""")

    if Type_of_query == "1": 
        PES_landscaper(Energies_in, r_in,  theta_in)
        exit()
    
    elif Type_of_query == "2": 
        PES_Vibrational_Freq(Energies_in, r_in, theta_in, range_in)
        exit()

    elif Type_of_query == "3": 
        exit()
    
    else: 
        Type_of_query = input("""Which function should I perform? (Please use the indicies next to each option)
    1. Calculate and graph the PES
    2. Calculate the vibrational frequencies""")




    
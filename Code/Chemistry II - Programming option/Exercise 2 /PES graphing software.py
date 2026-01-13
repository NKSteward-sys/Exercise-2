import numpy as np
import re
import matplotlib.pyplot as plt

#Explain this! 
def PES_Filereader(Folder): 

    """The first step in reading these directories is generating the list of filenames. I have used loops to generate a list of filenames
    using hardcoded parameters dependening on the input directory."""

    filenames = [] 

    if "H2Ooutfiles" in Folder:

        r_range = np.arange(0.70, 1.95, 0.05)
        molecule = "H2O"

    elif "H2Soutfiles" in Folder: 

        r_range = np.arange(0.70, 1.85, 0.05)
        molecule = "H2S"
    
    else: 
        print("Sorry, this isn't a directory I can parse. Make sure you've pasted the full file path.")

    for r in r_range:
        for theta in range(70, 161, 1):
            filenames.append(f"{Folder}/{molecule}.r{r:.2f}theta{theta}.0.out")
            
    """Now that we have the list of filenames. We can iteratively open each file and read the contents. Using regular expressions, we search for the line with the energy
      and extract the numerical value of the energy. We then close the file."""

    Energies = []
    for name in filenames: 
        Gaussian_Output = open(name, 'r')

        for line in Gaussian_Output: 
            match = re.search(r"E\(RHF\)\s*=?\s*([-+]?\d+\.\d+)", line)
            if match: 
                Energies.append(float(match.group(1)))
                break

        Gaussian_Output.close()
    
    print(len(Energies))







if __name__ == "__main__":
    UserFolder = input("What is the folder name?")

    PES_Filereader(UserFolder)

#NB E is not in line 168 with H2O r = 1.9 and theta = 153 
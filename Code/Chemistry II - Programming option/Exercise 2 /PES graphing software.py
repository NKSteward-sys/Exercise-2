import numpy as np

#Explain this! 
def PES_Filereader(Folder): 
    
    filenames = [] 
    idx = 0

    if Folder == "H2Ooutfiles":

        for r in np.arange(0.70, 1.95, 0.05): 
            for theta in range(70, 161, 1):
                filenames.append(f"{UserFolder}/H2O.r{r:.2f}theta{theta}.0.out")
                idx = idx + 1

    elif Folder == "H2Soutfiles": 

        for r in np.arange(0.60, 1.85, 0.05): 
            for theta in range(70, 161, 1):
                filenames.append(f"{UserFolder}/H2S.r{r:.2f}theta{theta}.0.out")
                idx = idx + 1
    


    for name in filenames: 
        ORCA_Output = open(name, r)
        Lines = []
        for line in ORCA_Output: 
            Lines.append(line)
        print(Lines)



if __name__ == "__main__":
    UserFolder = input("What is the folder name?")

    PES_Filereader(UserFolder)


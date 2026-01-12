import numpy as np

#Explain this! 
def PES_Filereader(Folder): 
    
    filenames = [] 
    idx = 0

    if Folder == "H2Ooutfiles":

        for r in np.arange(0.70, 1.95, 0.05): 
            for theta in range(70, 161, 1):
                filenames.append(f"H2O.r{r:.2f}theta{theta}.0.out")
                idx = idx + 1

    elif Folder == "H2Soutfiles": 

        for r in np.arange(0.60, 1.85, 0.05): 
            for theta in range(70, 161, 1):
                filenames.append(f"H2S.r{r:.2f}theta{theta}.0.out")
                idx = idx + 1
    




    print(filenames)



if __name__ == "__main__":
    Folder = input("What is the folder name?")

    PES_Filereader(Folder)


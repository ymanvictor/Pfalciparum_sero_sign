# R Code

The 'src' folder contain R code for identifying serological markers of recent P. falciparum exposure as described in Yman et al. Nat. Comms. 2021

1. Execute 'cloned_repo_init.R'
    
  This will establish a local version of R including the appropriate version of the necessary dependencies for executeing 'Pfalciparum_sero_sign.R'
    
2. Execute 'Pfalciparum_sero_sign.R'

  This script will reproduce results and figures from Yman et al. Nat. Comms. 2021. Figures will be saved to file within the "figures" folder.

  The script includes a demo run of random forest classification of recent malaria exposure based on a combination of antibody responses.
  The demo run includes analysis of 10 combinations of each of two to five responses. Expected runtime of the demo on a standard desktop computer is < 5 minutes.

  Instructions on how to edit the code to perform an exhaustive evaluation of all combinations of 2-5 antibody responses are provided within the script itself. 
  A full exhaustive evaluation will take several days to run on a desktop computer even using parallel backends.




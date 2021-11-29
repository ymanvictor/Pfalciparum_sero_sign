# Repository README

The "Pfalciparum_sero_sign" repository contains code and data for identification of serological markers of recent malaria exposure as described in 
Yman et al. Nature Communications 2021.



**Note: From now on, everything described below will be executed in RStudio**

<a name="init_renv"></a>
## Initiating your R environment (Reproducibility feature)
`src/project_init.R` Run the script to initiate the R environment and connect it to the R Project. Consent with **yes** when asked.
  * From now, all of your used r packages and dependencies will now be recorded in the ```renv.loc``` file when you execute ```renv::snapshot()```.
  * It is recomended to always use ```ren::hydrate()``` when installing new packages to your project enviroment. It will retrive the packages if allready installed on your computer, or install it from the CRAN r repository if missing. The packages that don't exists on CRAN have to be installed as usuall with ```install.packages()```


### 1. Automatic connection to your environment (renv)
Changes in your environment, e.g. by installing new packages, have to be capture by executing `renv::snapshot()` inside your R markdown.


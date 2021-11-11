# Project template to start a new computational project
### Our suggestion for project organisation and better reproducibility

## General info
This repository contains our idea project directory structure to make computational research projects within our group easier to work with and more reproducible. This project structure is currently our best suggestion but we are happy for suggestions and discussions.

## Table of contents
* [General info](#general-info)
* [Repo description](#repo-description)
* [Getting started](#getting-started)
* [Initiating your R environment (Reproducibility feature)](#init_renv)
* [Create R markdowns from template (Reproducibility feature)](#markdown)
* [Thoughts and Comments](#thoughs-comments)


## Repo description
This repository provides a template for new computational projects. It contains a pre-defined structure, as well as some features for reproducible research using an R project and version control via GitHub.
  * Pre-defined directory structure with readme files
  * Description on how to connect a template-based repository with rstudio
  * Several features to help analysis reproducibility

<a name="getting-started"></a>
## :technologist: Getting started

### Creating a your own repository from this repository template
1.	Create a new GitHub repository based on [this template](https://github.com/LautenbachMJ/project_template). Click on the green button (“Use this template”).

2. Go to your new repository and copy the url (green button “Code”). You should have copied something like `https://github.com/yourGitHubName/yourRepoName.git`.
 
### Connecting your GitHub Repository to your RStudio project
3.	You can do this either in RStudio or in the terminal. When using RStudio, click on (```File/New Project..```) in the menu bar. Then select "from VersionControl" and "Git". Paste the copied URL and give the respository a name. This will connect your GitHub repository to your R project and allows version control.


**Note: From now on, everything described below will be executed in RStudio**

<a name="init_renv"></a>
## Initiating your R environment (Reproducibility feature)
`src/project_init.R` Run the script to initiate the R environment and connect it to the R Project. Consent with **yes** when asked.
  * From now, all of your used r packages and dependencies will now be recorded in the ```renv.loc``` file when you execute ```renv::snapshot()```.
  * It is recomended to always use ```ren::hydrate()``` when installing new packages to your project enviroment. It will retrive the packages if allready installed on your computer, or install it from the CRAN r repository if missing. The packages that don't exists on CRAN have to be installed as usuall with ```install.packages()```


### 1. Automatic connection to your environment (renv)
Changes in your environment, e.g. by installing new packages, have to be capture by executing `renv::snapshot()` inside your R markdown.


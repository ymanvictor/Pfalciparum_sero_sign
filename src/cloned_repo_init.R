## project initiation script
## when repo cloned from github with excisting 'renv.lock' file

if (!require("renv")) install.packages("renv")

library("renv")

renv::restore()

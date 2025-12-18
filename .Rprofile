# The GEDI container has its own installation of R and several packages, which will conflict with renv if left alone. Therefore, we need to disable the renv autoloader if we are running in a container.
if (Sys.getenv("SINGULARITY_CONTAINER") != "" || 
    Sys.getenv("APPTAINER_CONTAINER") != "") {
  Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
} else {
  source("renv/activate.R")
}
source("renv/activate.R")

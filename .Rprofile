if (Sys.getenv("SINGULARITY_CONTAINER") != "" || 
    Sys.getenv("APPTAINER_CONTAINER") != "") {
  Sys.setenv(RENV_CONFIG_AUTOLOADER_ENABLED = "FALSE")
} else {
  source("renv/activate.R")
}
source("renv/activate.R")

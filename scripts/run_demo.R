root_dir <- normalizePath(Sys.getenv("AGEING_ROOT", unset = getwd()), mustWork = TRUE)
Sys.setenv(AGEING_ROOT = root_dir, AGEING_DEMO = "TRUE")

source(file.path(root_dir, "scripts", "create_demo_data.R"))
source(file.path(root_dir, "scripts", "01_prepare_gbd_inputs.R"))
source(file.path(root_dir, "scripts", "02_define_age_related_causes.R"))
source(file.path(root_dir, "scripts", "03_equivalent_age.R"))
source(file.path(root_dir, "scripts", "04_projection_nordpred.R"))
source(file.path(root_dir, "scripts", "06_decomposition_shapley.R"))

message("Demo completed. Outputs are in outputs/demo/.")

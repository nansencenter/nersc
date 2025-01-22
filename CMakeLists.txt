add_library(fabm_models_nersc OBJECT
            nersc_model_library.F90
            ecosmo/ecosmo_operational.F90
            ecosmo/ecosmo.F90
            ecosmo/ecosmo.GMD.2023_25.F90
            dvm/conservative_migrator.F90
            dvm/get_dependencies.F90
            dvm/upper_lower_boundaries.F90
            dvm/weight_distribution.F90
            dvm/move.F90
            dvm/apply_move.F90
           )

add_dependencies(fabm_models_nersc fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

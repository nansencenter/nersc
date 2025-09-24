add_library(fabm_models_nersc OBJECT
            nersc_model_library.F90
            ecosmo/ecosmo_operational.F90
            ecosmo/ecosmo.F90
            ecosmo/ecosmo.GMD.2023_25.F90
            dvm/conservative_migrator.F90
            dvm/get_dependencies.F90
            dvm/upper_lower_boundaries_simple.F90
            dvm/weight_distribution.F90
            dvm/move.F90
            # below this is a modular version of ecosmo
            # based on the version that configures
            # mesozooplankton diel vertical migration from 2025 NECCTON hindcast
            modular/shared.F90
            modular/nutrient.F90
            modular/oxygen.F90
            modular/organic_matter.F90
            modular/opal.F90
            modular/phytoplankton.F90
            modular/diatom.F90
            modular/coccolith.F90
            modular/phy.F90
            modular/zooplankton.F90
            modular/caco3.F90
            modular/sediment.F90
           )

add_dependencies(fabm_models_nersc fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

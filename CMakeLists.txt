add_library(fabm_models_nersc OBJECT
            nersc_model_library.F90
            ecosmo/ecosmo_operational.F90
            ecosmo/ecosmo.F90
            ecosmo/ecosmo.GMD.2023_25.F90
            dvm/conservative_migrator.F90
            dvm/get_dependencies.F90
            dvm/upper_lower_boundaries_operational.F90
	        dvm/weight_distribution.F90
            dvm/move.F90
	        ice_algae/icealgae_G.F90
	        light/light.F90
            # MODULAR ECOSMO COMPONENTS
            modular/shared.F90
            modular/nutrient.F90
            modular/oxygen.F90
            modular/organic_matter.F90
            modular/opal.F90
            modular/phy.F90
            modular/zooplankton.F90
            modular/caco3.F90
            modular/sediment.F90
            # ADD OTHER MODELS (e.g. development versions) HERE
           )

add_dependencies(fabm_models_nersc fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------

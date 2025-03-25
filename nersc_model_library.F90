module nersc_model_library

   use fabm_types, only: type_base_model_factory,type_base_model

   implicit none

   private

   type,extends(type_base_model_factory) :: type_factory
      contains
      procedure :: create
   end type

   type (type_factory),save,target,public :: nersc_model_factory

contains

   subroutine create(self,name,model)

      use fabm_nersc_ecosmo_operational
      use fabm_nersc_ecosmo
      use fabm_nersc_ecosmo_gmd_2023_25
      use dvm_conservative_migrator
      use dvm_get_dependencies
      use dvm_upper_lower_boundaries
      use dvm_upper_lower_boundaries_simple
      use dvm_upper_lower_boundaries_complex
      use dvm_weight_distribution
      use dvm_move
      ! Add new models here

      class (type_factory),intent(in) :: self
      character(*),        intent(in) :: name
      class (type_base_model),pointer :: model

      select case (name)
         case ('ecosmo_operational');       allocate(type_nersc_ecosmo_operational::model)
         case ('ecosmo');       allocate(type_nersc_ecosmo::model)
         case ('ecosmo_gmd_2023_25');       allocate(type_nersc_ecosmo_gmd_2023_25::model)
         case ('dvm_conservative_migrator'); allocate(type_conservative_migrator::model)
         case ('dvm_get_dependencies'); allocate(type_get_dependencies::model)
         case ('dvm_upper_lower_boundaries_simple'); allocate(type_upper_lower_boundaries_simple::model)
         case ('dvm_weight_distribution'); allocate(type_weight_distribution::model)
         case ('dvm_move'); allocate(type_move::model)
         ! Add new models here
      end select

   end subroutine



end module

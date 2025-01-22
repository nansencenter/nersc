#include "fabm_driver.h"

module dvm_get_dependencies

use fabm_types
use fabm_expressions
!use ecosmo_shared

implicit none

private

type, extends(type_base_model), public :: type_get_dependencies

    type (type_horizontal_diagnostic_variable_id) :: id_light_present0
!    type (type_horizontal_diagnostic_variable_id) :: id_integrated_food
    type (type_diagnostic_variable_id) :: id_migrator_food
    type (type_horizontal_dependency_id) :: id_par0

    integer  :: nprey
    type (type_state_variable_id),  allocatable,dimension(:) :: id_prey

    contains
        procedure :: initialize
        procedure :: do_surface
        procedure :: do

end type

contains

    subroutine initialize(self, configunit)
        class (type_get_dependencies), intent(inout), target :: self
        integer, intent(in)                                  :: configunit
        integer                                              :: iprey
        character(len=16)                                    :: index

        call self%register_diagnostic_variable(self%id_light_present0,'light_presence','-','light is available at the surface',source=source_do_surface)
        call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)

        call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)
        call self%register_diagnostic_variable(self%id_migrator_food,'migrator_food','mgC/m3','food availability for the migrators', act_as_state_variable=.true., missing_value=0.0_rk, source=source_do)
        ! Get prey-specific coupling links.
        allocate(self%id_prey(self%nprey))
        do iprey=1,self%nprey
            write (index,'(i0)') iprey
            call self%register_state_dependency(self%id_prey(iprey),'prey'//trim(index)//'','mgC/m3', 'prey '//trim(index)//' carbon concentration')
        end do

    end subroutine initialize

    subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

        class (type_get_dependencies),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_SURFACE_

        real(rk) :: par0
        _HORIZONTAL_LOOP_BEGIN_

            _GET_SURFACE_(self%id_par0,par0)
            if (par0 < 5.0_rk) then
                _SET_HORIZONTAL_DIAGNOSTIC_(self%id_light_present0, 0.0_rk)
            else
                _SET_HORIZONTAL_DIAGNOSTIC_(self%id_light_present0, 24.0_rk/86400_rk) 
            end if

        _HORIZONTAL_LOOP_END_

    end subroutine do_surface

    subroutine do(self, _ARGUMENTS_DO_)
        class (type_get_dependencies), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_

        integer  :: iprey
        real(rk),dimension(self%nprey) :: prey

        _LOOP_BEGIN_

            ! Get prey concentrations
            do iprey=1,self%nprey
                _GET_(self%id_prey(iprey), prey(iprey))
            end do
            _SET_DIAGNOSTIC_(self%id_migrator_food, sum(prey))

        _LOOP_END_

    end subroutine do

end module

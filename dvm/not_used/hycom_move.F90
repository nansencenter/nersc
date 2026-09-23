#include "fabm_driver.h"

module dvm_hycom_move

use fabm_types
use fabm_expressions

implicit none 

private 

public type_hycom_move

type, extends(type_base_model) :: type_hycom_move
    type (type_dependency_id)                      :: id_distributed
    type (type_dependency_id)                      :: id_thickness
    type (type_state_variable_id)                  :: id_target
    type (type_bottom_dependency_id)               :: id_integral
    type (type_bottom_dependency_id)               :: id_distributed_integral
    contains
        procedure :: initialize
        procedure :: check_state
    end type

contains

subroutine initialize(self, configunit)

    class (type_hycom_move), intent(inout), target :: self
    integer,         intent(in)            :: configunit

    call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
    call self%register_dependency(self%id_distributed,'migrator_distributed','mgC/m3','migrators final concentration')
    call self%register_dependency(self%id_integral,'integral','mgC/m2','depth-integrated target variable')
    call self%register_dependency(self%id_distributed_integral,'migrator_integrated_mass','mgC/m2','migrators final integrated mass')
    call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)

end subroutine initialize

    subroutine check_state(self,_ARGUMENTS_CHECK_STATE_)
        class (type_hycom_move), intent(in) :: self
        _DECLARE_ARGUMENTS_CHECK_STATE_

        real(rk) :: integral
        real(rk) :: distributed_integral
        real(rk) :: local
        real(rk) :: distributed
        real(rk) :: thickness
        real(rk) :: final_concentration
        real(rk) :: difference

        _LOOP_BEGIN_
           _GET_BOTTOM_(self%id_integral,integral)
           _GET_BOTTOM_(self%id_distributed_integral,distributed_integral)
           _GET_(self%id_target,local)
           _GET_(self%id_distributed,distributed)
           _GET_(self%id_thickness,thickness)

            ! difference = distributed_integral - integral
            ! !write(*,*) difference / distributed_integral
            ! if ( abs(difference) < 0.025_rk * integral) then
                final_concentration = distributed !* (1.0 - difference / distributed_integral )

                _SET_(self%id_target, final_concentration)
            ! else
            !     _SET_(self%id_target,local) 
            ! end if

        _LOOP_END_

    end subroutine check_state

end module
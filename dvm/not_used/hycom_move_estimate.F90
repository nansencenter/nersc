#include "fabm_driver.h"

module dvm_hycom_move_estimate

use fabm_types
use fabm_expressions

implicit none 

private 

public type_hycom_move_estimate

type, extends(type_base_model) :: type_hycom_move_estimate
    type (type_dependency_id)                      :: id_random_weights
    type (type_dependency_id)                      :: id_thickness
    type (type_state_variable_id)                  :: id_target
    type (type_bottom_dependency_id)               :: id_integral
    type (type_bottom_dependency_id)               :: id_integral_random_weights
    type (type_diagnostic_variable_id)             :: id_distributed
    type (type_bottom_diagnostic_variable_id)      :: id_distributed_integral
    real(rk) :: ratioMig
    real(rk) :: allowed_thickness
    contains
        procedure :: initialize
        procedure :: do_column
    end type

contains

subroutine initialize(self, configunit)

    class (type_hycom_move_estimate), intent(inout), target :: self
    integer,         intent(in)            :: configunit

    call self%register_state_dependency(self%id_target, 'target', '', 'variable to apply sources and sinks to')
    call self%register_dependency(self%id_random_weights,'migrator_random_weights','-','migrators distribution random weights')
    call self%register_dependency(self%id_integral,'integral','','depth-integrated target variable')
    call self%register_dependency(self%id_integral_random_weights,'migrator_integral_random_weights','','migrators distribution integral random weights')
    call self%register_diagnostic_variable(self%id_distributed,'migrator_distributed','mgC/m3','migrators final concentration', missing_value=0.0_rk, source=source_do_column)
    call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
    call self%get_parameter( self%ratioMig,  'ratioMig',    '-',      'ratio of moving biomass', default=0.5_rk)
    call self%register_diagnostic_variable(self%id_distributed_integral,'migrator_integrated_mass','-','migrators final integrated mass', missing_value=0.0_rk, source=source_do_column)
    call self%get_parameter(self%allowed_thickness,'allowed_thickness','m','layer thickness where concentration will be calculated, else a valid layer is copied',default=1.0E-20_rk)

end subroutine initialize

    subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)

        class (type_hycom_move_estimate), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_COLUMN_

        real(rk) :: integral, integral_random_weights, integral_random_weights_fabm, random_weights
        real(rk) :: local, distributed, thickness, target0
        real(rk) :: mortality_switch, local_loss
        real(rk) :: valid_concentration
        real(rk) :: calculated_concentration
        real(rk) :: distributed_integral_biomass
        integer  :: first_layer

        distributed_integral_biomass = 0.0_rk
        first_layer = 1
        _VERTICAL_LOOP_BEGIN_
           _GET_BOTTOM_(self%id_integral,integral)
           _GET_BOTTOM_(self%id_integral_random_weights,integral_random_weights)
           _GET_(self%id_target,local)
           _GET_(self%id_random_weights,random_weights)
           _GET_(self%id_thickness,thickness)
    
           distributed = random_weights / integral_random_weights ! this ensures that total weight distribution is 1, so mass conserved
           calculated_concentration = max(0.0_rk,local * (1.0_rk - self%ratioMig) + (distributed * integral/max(thickness,1.0E-20_rk)) * self%ratioMig )
           if (first_layer == 1) then
                first_layer = 0
                valid_concentration = calculated_concentration
                _SET_DIAGNOSTIC_(self%id_distributed,calculated_concentration)
           else
                if (thickness >= self%allowed_thickness .and. local >= 0.0_rk) then
                    valid_concentration = calculated_concentration
                    _SET_DIAGNOSTIC_(self%id_distributed,valid_concentration)
                else
                    _SET_DIAGNOSTIC_(self%id_distributed,valid_concentration) 
                end if
           end if
           distributed_integral_biomass = distributed_integral_biomass + valid_concentration * max(thickness,1.0E-20_rk)

        _VERTICAL_LOOP_END_
        _SET_BOTTOM_DIAGNOSTIC_(self%id_distributed_integral,distributed_integral_biomass)
    end subroutine do_column

end module
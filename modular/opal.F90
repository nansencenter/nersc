#include "fabm_driver.h"

module ecosmo_opal

use fabm_types
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_opal

! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_opal
    type (type_state_variable_id)         :: id_c
    type (type_state_variable_id)         :: id_sil, id_opal
    real(rk) :: regenSi, OMsink
    contains

!     Model procedures
    procedure :: initialize
    procedure :: do

end type type_ecosmo_opal

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_opal), intent(inout),target  :: self
    integer,  intent(in) :: configunit

    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    call self%get_parameter( self%regenSi,  'regenSi',    '1/day',      'Si regeneration rate',            default=0.015_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%OMsink, 'OMsink', 'm/day', 'organic matter sinking rate', default=0.0_rk, scale_factor=1.0_rk/sedy0)
    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'concentration in carbon units', minimum=0.0_rk, vertical_movement=-self%OMsink)
    call self%register_state_dependency(self%id_sil, 'sil', 'mgC/m3', 'silicate')

end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)

    class (type_ecosmo_opal),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    real(rk) :: c

    _LOOP_BEGIN_

    _GET_(self%id_c,c)
    _ADD_SOURCE_(self%id_opal, -self%regenSi * c)
    _ADD_SOURCE_(self%id_sil, self%regenSi * c)
    
    _LOOP_END_
end subroutine do
end module
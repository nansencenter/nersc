#include "fabm_driver.h"

module ecosmo_nutrient

use fabm_types
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_nutrient

! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_nutrient
    type (type_state_variable_id)         :: id_c
    contains

!     Model procedures
    procedure :: initialize

end type type_ecosmo_nutrient

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_nutrient), intent(inout),target  :: self
    integer,  intent(in) :: configunit
    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'concentration in carbon units', minimum=0.0_rk)

end subroutine initialize

end module
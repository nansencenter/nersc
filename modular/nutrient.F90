#include "fabm_driver.h"

! --------- Change Log ---------- !
! Veli Çağlar Yumruktepe (VCY)
!
! VCY - 27/05/2026
! This is the first iteration of the modular version of ECOSMO II(CHL).
! It is based on the ECOSMO II(CHL) code used for Copernicus ARC MFC 2026 operational model code and parameters.
! Missing components: 
!   1) community dependent organic matter sinking speed
!   2) sea-ice algae (and fast sinking detritus implementation) 
! ------------------------------- !

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
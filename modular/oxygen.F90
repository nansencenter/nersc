#include "fabm_driver.h"

module ecosmo_oxygen

use fabm_types
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_oxygen

! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_oxygen
    type (type_state_variable_id)         :: id_c
    type (type_state_variable_id)         :: id_no3, id_nh4, id_alk
    type (type_dependency_id)             :: id_temp

    contains

!     Model procedures
    procedure :: initialize
    procedure :: do
!    procedure :: do_surface

end type type_ecosmo_oxygen

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_oxygen), intent(inout),target  :: self
    integer,  intent(in) :: configunit
    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    call self%register_state_variable(self%id_c, 'c', 'mmol/m3', 'oxygen',minimum=0.0_rk)
    call self%register_state_dependency(self%id_no3, 'no3', 'mgC/m3', 'nitrate')
    call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
    call self%register_dependency(self%id_temp,standard_variables%temperature)
    if (couple_co2) then
        call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
    end if
end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)

    class (type_ecosmo_oxygen),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    real(rk) :: oxy, no3, nh4, temp
    real(rk) :: Onitr, bioom1, bioom6
    real(rk) :: rhs_oxy, rhs_amm, rhs_nit, rhs_alk

    _LOOP_BEGIN_

    _GET_(self%id_c,oxy)
    _GET_(self%id_no3,no3)
    _GET_(self%id_nh4,nh4)
    _GET_(self%id_temp, temp)    

    Onitr = 0.01_rk * O2ml_l_to_O2mmol_m3 !according to Neumann  (Onitr in mlO2/l see also Stigebrand and Wulff)
    bioom1 = 0.0_rk
    bioom6 = 0.0_rk
    if (oxy > 0) then
        bioom1 = 0.1_rk/sedy0 * exp(temp*0.11_rk) * oxy/(Onitr+oxy)
        bioom6 = 1.0_rk
    end if
    rhs_amm = - bioom1 * nh4 
    _ADD_SOURCE_(self%id_nh4, rhs_amm )
    rhs_nit = bioom1 * nh4
    _ADD_SOURCE_(self%id_no3, rhs_nit )

    rhs_oxy = -2.0_rk * bioom1 * nh4 * Cmg_to_Cmmol * Cmmol_to_Nmmol 
    _ADD_SOURCE_( self%id_c, rhs_oxy ) 

    if (couple_co2) then
        rhs_alk = (rhs_amm - rhs_nit) * Cmg_to_Cmmol * Cmmol_to_Nmmol -0.5_rk * rhs_oxy * (1._rk-bioom6)
        _ADD_SOURCE_(self%id_alk, rhs_alk )
    end if

    _LOOP_END_
end subroutine do

end module
#include "fabm_driver.h"

module ecosmo_caco3

use fabm_types
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_caco3

! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_caco3
    type (type_state_variable_id)         :: id_c
    type (type_state_variable_id)         :: id_alk, id_dic
    type (type_dependency_id)         :: id_Om_cal
    real(rk) :: IMsink, calcDis
    contains

!     Model procedures
    procedure :: initialize
    procedure :: do

end type type_ecosmo_caco3

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_caco3), intent(inout),target  :: self
    integer,  intent(in) :: configunit

    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    call self%get_parameter( self%IMsink, 'IMsink', 'm/day', 'inorganic matter sinking rate', default=5.0_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%calcDis , 'calcDis',        '1/day',      'calcite dissolution rate',  default=0.03_rk,  scale_factor=1.0_rk/sedy0)

    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'concentration in carbon units', minimum=0.0_rk, vertical_movement=-self%IMsink)
    call self%register_dependency(self%id_Om_cal,'Om_cal_target','-','calcite saturation')
    call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
    call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')

end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)

    class (type_ecosmo_caco3),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    real(rk) :: c
    real(rk) :: Om_cal, Lstar

    _LOOP_BEGIN_

    _GET_(self%id_c,c)
    _GET_(self%id_Om_cal,Om_cal)


    Lstar = self%calcDis * max(0.0,1.0 - Om_cal)
    _ADD_SOURCE_(self%id_c, -Lstar * c)
    _ADD_SOURCE_(self%id_dic, Lstar * c)
    _ADD_SOURCE_(self%id_alk, Lstar * c * Cmg_to_Cmmol * 2.0_rk)
    
    _LOOP_END_
end subroutine do
end module
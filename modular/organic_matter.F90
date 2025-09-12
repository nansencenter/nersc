#include "fabm_driver.h"

module ecosmo_organic_matter

use fabm_types
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_organic_matter

! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_organic_matter
    type (type_state_variable_id)         :: id_c
    type (type_dependency_id)             :: id_temp
    type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_oxy
    type (type_state_variable_id)         :: id_alk, id_dic

    real(rk) :: remin
    real(rk) :: OMsink
    contains

!     Model procedures
    procedure :: initialize
    procedure :: do
end type type_ecosmo_organic_matter

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_organic_matter), intent(inout),target  :: self
    integer,  intent(in) :: configunit
    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    call self%get_parameter( self%OMsink, 'OMsink', 'm/day', 'organic matter sinking rate', default=0.0_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%remin, 'remin', '1/day', 'organic matter remineralization rate', default=0.003_rk, scale_factor=1.0_rk/sedy0)

    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'concentration in carbon units', minimum=0.0_rk, vertical_movement=-self%OMsink)
    call self%register_dependency(self%id_temp,standard_variables%temperature)
    call self%register_state_dependency(self%id_no3, 'no3', 'mgC/m3', 'nitrate')
    call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
    call self%register_state_dependency(self%id_pho, 'pho', 'mgC/m3', 'phosphate')
    call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol/m3', 'oxygen')
    if (couple_co2) then
        call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
        call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
    end if

end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)

    class (type_ecosmo_organic_matter),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    real(rk) :: frem
    real(rk) :: c, temp
    real(rk) :: no3, pho, nh4, oxy, dic, alk
    real(rk) :: bioom5, bioom6, bioom7
    real(rk) :: rhs_oxy, rhs, rhs_amm, rhs_nit

    _LOOP_BEGIN_

    ! Retrieve current (local) state variable values.
    _GET_(self%id_c,c)
    _GET_(self%id_temp,temp)
    _GET_(self%id_no3,no3)
    _GET_(self%id_nh4,nh4)
    _GET_(self%id_pho,pho)
    _GET_(self%id_oxy,oxy)
    if (couple_co2) then
        _GET_(self%id_dic,dic)
        _GET_(self%id_alk,alk)
    end if

    bioom5 = 0.0_rk
    bioom6 = 0.0_rk
    bioom7 = 0.0_rk
    if (oxy > 0) then
      bioom6 = 1.0_rk
    else
      if (no3>0) then
        bioom5 = 5.0_rk
      else
        bioom7 = 1.0_rk
      end if
    end if

    frem = self%remin * (1._rk+20._rk*(temp**2/(13._rk**2+temp**2)))
    
    _ADD_SOURCE_(self%id_c, -frem * c)
    rhs_amm = frem * c
    rhs_nit = -frem * c * bioom5
    _ADD_SOURCE_(self%id_no3, rhs_nit )
    _ADD_SOURCE_(self%id_nh4, rhs_amm )
    _ADD_SOURCE_(self%id_pho, frem * c )
    rhs_oxy = -frem * c * (bioom6+bioom7) * 6.625 * Cmg_to_Cmmol* Cmmol_to_Nmmol
    _ADD_SOURCE_(self%id_oxy, rhs_oxy)

    if (couple_co2) then
      rhs = frem * c * Cmg_to_Cmmol
      _ADD_SOURCE_(self%id_dic, rhs )
      rhs = (rhs_amm - rhs_nit) * Cmg_to_Cmmol * Cmmol_to_Nmmol - 0.5_rk * rhs_oxy * (1._rk-bioom6)
      _ADD_SOURCE_(self%id_alk, rhs )
    end if

    _LOOP_END_

end subroutine do

end module
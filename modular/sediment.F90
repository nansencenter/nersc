#include "fabm_driver.h"

module ecosmo_sediment

use fabm_types
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_sediment
! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_sediment
    type (type_bottom_state_variable_id)         :: id_sed1, id_sed2, id_sed3, id_sed4
    type (type_state_variable_id)         :: id_no3, id_pho, id_sil, id_oxy, id_det, id_opa, id_alk, id_dic, id_nh4, id_caco3
    type (type_dependency_id)             :: id_temp
    type (type_horizontal_dependency_id)  :: id_tbs
    type (type_dependency_id)             :: id_thickness

    real(rk) :: crBotStr, resuspRt, sedimRt, burialRt
    real(rk) :: reminSED, TctrlDenit, RelSEDp1, RelSEDp2, reminSEDsi
    contains

!     Model procedures
    procedure :: initialize
    procedure :: do_bottom

end type type_ecosmo_sediment

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_sediment), intent(inout),target  :: self
    integer,  intent(in) :: configunit
    ! crBotStr is formulated to be 0.1, see do_bottom subroutine 
    call self%get_parameter( self%crBotStr,  'crBotStr',   'N/m**2',     'critic. bot. stress for resusp.', default=0.1_rk)
    !
    call self%get_parameter( self%resuspRt,  'resuspRt',   '1/day',      'resuspension rate',               default=25.0_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%sedimRt,  'sedimRt',    'm/day',      'sedimentation rate',              default=3.5_rk,   scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%burialRt,  'burialRt',   '1/day',      'burial rate',                     default=1E-5_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%reminSED,  'reminSED',   '1/day',      'sediment remineralization rate',  default=0.001_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%TctrlDenit,  'TctrlDenit', '1/degC',     'temp. control denitrification',   default=0.15_rk)
    call self%get_parameter( self%RelSEDp1,  'RelSEDp1',   'units??',    'P sedim. rel. p1',                default=0.15_rk)
    call self%get_parameter( self%RelSEDp2,  'RelSEDp2',   'units??',    'P sedim. rel. p2',                default=0.10_rk)
    call self%get_parameter( self%reminSEDsi,  'reminSEDsi', '1/day',      'sed. remineralization rate Si',   default=0.0002_rk,scale_factor=1.0_rk/sedy0)

    call self%register_state_variable( self%id_sed1,     'sed1',    'mgC/m2',    'sediment detritus',         minimum=0.0_rk , maximum=1E4_rk, &
        initial_value = 20.0_rk * Nmmol_to_Cmmol * Cmmol_to_Cmg * Nmg_to_Nmmol )
    call self%register_state_variable( self%id_sed2,     'sed2',    'mgC/m2',    'sediment opal',             minimum=0.0_rk , maximum=1E4_rk, &
        initial_value = 20.0_rk * Simmol_to_Cmmol * Cmmol_to_Cmg * Simg_to_Simmol )
    call self%register_state_variable( self%id_sed3,     'sed3',    'mgC/m2',    'sediment adsorbed pho.',    minimum=0.0_rk  , maximum=1E4_rk, &
        initial_value=2.0_rk * Pmmol_to_Cmmol * Cmmol_to_Cmg * Pmg_to_Pmmol )

    call self%register_dependency(self%id_temp,standard_variables%temperature)
    call self%register_state_dependency(self%id_no3, 'no3', 'mgC/m3', 'nitrate')
    call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
    call self%register_state_dependency(self%id_pho, 'pho', 'mgC/m3', 'phosphate')
    call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol/m3', 'oxygen')
    call self%register_state_dependency(self%id_det, 'det', 'mgC/m3', 'detritus')
    call self%register_state_dependency(self%id_opa, 'opal', 'mgC/m3', 'opal')    
    if (couple_co2) then
        call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
        call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
        if (.false.) then
            call self%register_state_variable( self%id_sed4,     'sed4',    'mmol/m2',    'sediment calcite',         minimum=0.0_rk , &
            initial_value=1e-2_rk * Cmmol_to_Cmg )
            call self%register_state_dependency(self%id_caco3, 'caco3', 'mmol/m3', 'calcite') 
        end if
    end if
    call self%register_dependency(self%id_tbs,standard_variables%bottom_stress)
    call self%register_dependency(self%id_thickness, standard_variables%cell_thickness)
end subroutine initialize

subroutine do_bottom(self,_ARGUMENTS_DO_BOTTOM_)
    class (type_ecosmo_sediment),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_BOTTOM_

    real(rk) :: long_time_step_for_assumed_sedimentation_flux = 1200.0_rk
    real(rk) :: temp, tbs, oxy, no3, det, opa, pho

    real(rk) :: caco3
    real(rk) :: caco3_sedimentation
    real(rk) :: sed4, rhs_sed4

    real(rk) :: sed1, sed2, sed3
    real(rk) :: opal_sedimentation
    real(rk) :: det_sedimentation
    real(rk) :: bioom1, bioom2, bioom3, bioom4, bioom5, bioom6, bioom7, bioom8
    real(rk) :: Rsd, Rds, Rsa, Rsdenit, rhs, flux
    real(rk) :: Rsa_p, yt1, yt2
    real(rk) :: rhs_sed2 
    real(rk) :: thickness
    real(rk) :: alk_flux
    _HORIZONTAL_LOOP_BEGIN_
    
    _GET_(self%id_temp,temp)
    _GET_(self%id_oxy,oxy)
    _GET_(self%id_det,det)
    _GET_(self%id_no3,no3)
    _GET_(self%id_opa,opa)
    _GET_(self%id_thickness,thickness)
    _GET_HORIZONTAL_(self%id_sed1,sed1)
    _GET_HORIZONTAL_(self%id_sed2,sed2)
    _GET_HORIZONTAL_(self%id_sed3,sed3)
    _GET_HORIZONTAL_(self%id_tbs,tbs)

    _GET_(self%id_caco3,caco3)  
  
    bioom1 = 0.0_rk
    bioom2 = 0.0_rk
    bioom3 = 0.0_rk
    bioom4 = 0.0_rk
    bioom5 = 0.0_rk
    bioom6 = 0.0_rk
    bioom7 = 0.0_rk
    bioom8 = 0.0_rk
    if (oxy > 0) then
      bioom1 = 0.1/sedy0 * exp(temp * 0.11_rk) * oxy/( (0.1_rk * O2ml_l_to_O2mmol_m3) + oxy )
      bioom2 = bioom1
      bioom6 = 1.0_rk
    else
      if (no3>0) then
        bioom5 = 5.0_rk
        bioom8 = 1.0_rk
      else
        bioom7 = 1.0_rk
      end if
    end if
  
    !----citical bottom shear stress
    if (tbs .ge. self%crBotStr ) then
      Rsd=min(self%resuspRt, self%resuspRt * (tbs-0.1)**2 * 100.) ! sets to max=self%resuspRt when tbs=0.2,
                   ! sets to max=self%resuspRt when tbs=0.2,
                   ! else rapid increase to max when tbs=0.1
                   ! it assumes crBotStr = 0.1 in the fabm.yaml file
      Rds=0.0_rk
    else if (tbs .lt. self%crBotStr) then
        Rsd=0.0_rk
        Rds=self%sedimRt
        !det_loss = max(sign(-1.0_rk,det-0.5_rk),0.0_rk) ! TO BE TESTED, prevents DET going negative
        !Rds=self%sedimRt * det_loss
    end if    

    !---------------------------------------------------------------
    !----denitrification parameter in dependence of available oxygen
    if (oxy .gt. 0.0) then
        Rsa = self%reminSED * exp( self%TctrlDenit * temp ) * 1.0_rk
        Rsdenit = 0.0_rk
    else if (oxy .le. 0.0) then
        Rsdenit = self%reminSED * exp( self%TctrlDenit * temp ) * 2.0_rk
        Rsa=0.0_rk
    end if    

    !--- sediment 1 total sediment biomass and nitrogen pool
    det_sedimentation = Rds*det
    if (det_sedimentation * long_time_step_for_assumed_sedimentation_flux > det * thickness) det_sedimentation = 0.0_rk

    rhs = det_sedimentation - Rsd * sed1 - 2.0_rk * Rsa * sed1 - Rsdenit * sed1 - (2.0E-3 * self%burialRt * sed1 ) * sed1
    _SET_BOTTOM_ODE_(self%id_sed1, rhs)
    
    ! oxygen
    flux = -(&
        bioom6*6.625_rk * 2.0_rk * Rsa * sed1 &
        + bioom7 * 6.625_rk * Rsdenit * sed1 &
        + 2.0_rk * bioom1 * Rsa * sed1 & 
        ) * Cmg_to_Cmmol * Cmmol_to_Nmmol 
    _SET_BOTTOM_EXCHANGE_(self%id_oxy, flux)

    ! nitrate
    _SET_BOTTOM_EXCHANGE_(self%id_no3, -bioom5 * Rsdenit * sed1 )

    ! detritus
    _SET_BOTTOM_EXCHANGE_(self%id_det, Rsd * sed1 - det_sedimentation)

    ! ammonium
    _SET_BOTTOM_EXCHANGE_(self%id_nh4, (Rsdenit + Rsa) * sed1)

    if (couple_co2) then
        _SET_BOTTOM_EXCHANGE_(self%id_dic, Cmg_to_Cmmol * (Rsdenit + 2.0_rk * Rsa) * sed1 )
        flux = Cmg_to_Cmmol * Cmmol_to_Nmmol * ( (Rsdenit + Rsa + bioom5 * Rsdenit) * sed1 ) - 0.5_rk * flux * (1.0_rk-bioom6)
        _SET_BOTTOM_EXCHANGE_(self%id_alk, flux)
    end if

    !---------------------------------------------------------------
    !---- phosphate

    Rsa_p = self%reminSED * exp( self%TctrlDenit * temp ) * 2.0_rk

    if (oxy.gt.0.0) then
        yt2 = oxy/375.0_rk   !normieren des wertes wie in Neumann et al 2002
        yt1 = yt2**2.0_rk / ( self%RelSEDp2**2.0_rk + yt2**2.0_rk )

        _SET_BOTTOM_EXCHANGE_( self%id_pho, Rsa_p * (1.0_rk - self%RelSEDp1 * yt1) * sed3 )

        !--sed 3 phosphate pool sediment+remineralization-P release
        _SET_BOTTOM_ODE_( self%id_sed3, 2.0_rk * Rsa * sed1 - Rsa_p * (1.0_rk - self%RelSEDp1 * yt1 ) * sed3 )

    else if (oxy.le.0.0) then
        _SET_BOTTOM_EXCHANGE_( self%id_pho, Rsa_p * sed3 )
        _SET_BOTTOM_ODE_( self%id_sed3, Rsdenit * sed1 - Rsa_p * sed3 )
    end if    

    !---------------------------------------------------------------
    !---- sediment opal(Si)
    opal_sedimentation = 2.0*Rds*opa
    if (opal_sedimentation * long_time_step_for_assumed_sedimentation_flux > opa * thickness) opal_sedimentation = 0.0_rk

    rhs_sed2 = 2.0 * opal_sedimentation - Rsd * sed2 &
        - self%reminSEDsi * sed2 &
        - (2.0E-3 * self%burialRt * sed2 ) * sed2
    _SET_BOTTOM_ODE_(self%id_sed2, rhs_sed2)

    _SET_BOTTOM_EXCHANGE_(self%id_opa, Rsd * sed2 - opal_sedimentation)
    _SET_BOTTOM_EXCHANGE_(self%id_sil, self%reminSEDsi * sed2)

    !---------------------------------------------------------------
    !---- sediment CaCO3

    if (.false.) then
        caco3_sedimentation = Rds*caco3
        if (caco3_sedimentation * long_time_step_for_assumed_sedimentation_flux > caco3 * thickness) caco3_sedimentation = 0.0_rk

        _GET_HORIZONTAL_(self%id_sed4,sed4)
        _GET_(self%id_caco3, caco3)
        rhs_sed4 = caco3_sedimentation - Rsd*sed4 - self%reminSEDsi*sed4 &
            -( 2.0E-3 * self%burialRt * sed4 ) * sed4
        _SET_BOTTOM_ODE_(self%id_sed4, rhs_sed4)
        _SET_BOTTOM_EXCHANGE_(self%id_caco3, Rsd * sed4 - caco3_sedimentation)
    end if 

    if (couple_co2) then
        _SET_BOTTOM_EXCHANGE_(self%id_dic, Cmg_to_Cmmol * (Rsdenit+2*Rsa)*sed1)
        alk_flux = Cmg_to_Cmmol * Cmmol_to_Nmmol * ((Rsdenit+Rsa+bioom5*Rsdenit)*sed1) - 0.5_rk*flux*(1._rk-bioom6)
        _SET_BOTTOM_EXCHANGE_(self%id_alk, alk_flux)
    end if

    _HORIZONTAL_LOOP_END_
  
end subroutine do_bottom
  
end module
  
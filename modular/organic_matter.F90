#include "fabm_driver.h"

! --------- Change Log ---------- !
! Veli Çağlar Yumruktepe (VCY) github: @caglartac
! Philip Wallhead (PJW) github: @pwallhead
!
! VCY - 27/05/2026
! This is the first iteration of the modular version of ECOSMO II(CHL).
! It is based on the ECOSMO II(CHL) code used for Copernicus ARC MFC 2026 operational model code and parameters.
! Missing components: 
!   1) community dependent organic matter sinking speed
!   2) sea-ice algae (and fast sinking detritus implementation)
!
! VCY 28/05/2026:
!   Following the suggestions from PJW, added optional linear depth dependence parameters for sinking rates (cz_OMsink)
!      to achieve Martin curve flux variations following A15 and Middelburg (2019).
!      A15: Aumont et al. (2015), doi:10.5194/gmd-8-2465-2015
!      Middelburg (2019), doi:10.1007/978-3-030-10822-9
! 
! ------------------------------- !

module ecosmo_organic_matter
  use fabm_types
  use fabm_expressions
  use ecosmo_shared
  implicit none
  private
  public type_ecosmo_organic_matter

  type,extends(type_base_model), public  :: type_ecosmo_organic_matter
    type (type_state_variable_id)         :: id_c
    type (type_dependency_id)             :: id_temp
    type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_oxy
    type (type_state_variable_id)         :: id_alk, id_dic
    type (type_dependency_id)             :: id_depth
    
    real(rk) :: remin
    real(rk) :: OMsink
    real(rk) :: cz_OMsink
    real(rk) :: K_bact_NH4
    real(rk) :: K_bact_PO4

    type (type_state_variable_id)         :: id_dsnk
    type (type_diagnostic_variable_id)    :: id_snkspd
    real(rk)                              :: min_dsnk
    real(rk)                              :: max_dsnk
  contains
      procedure :: initialize
      procedure :: do
      procedure :: get_vertical_movement
      procedure :: check_state
  end type type_ecosmo_organic_matter

contains

  subroutine initialize(self,configunit)
    class (type_ecosmo_organic_matter), intent(inout),target  :: self
    integer,  intent(in) :: configunit

    ! Register parameters
    call self%get_parameter( self%OMsink, 'OMsink', 'm/day', 'organic matter sinking rate', default=0.0_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%remin,  'remin',  '1/day', 'organic matter remineralization rate', default=0.003_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%K_bact_NH4, 'K_bact_NH4', 'mgC/m3', 'Bacterial half-saturation for NH4', default=0.20_rk)
    call self%get_parameter( self%K_bact_PO4, 'K_bact_PO4', 'mgC/m3', 'Bacterial half-saturation for PO4', default=0.05_rk)
    call self%get_parameter( self%min_dsnk, 'min_dsnk', 'm/d', 'minimum community sinking speed', default=0.0_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%max_dsnk, 'max_dsnk', 'm/d', 'maximum community sinking speed', default=50.0_rk, scale_factor=1.0_rk/sedy0)

    ! Register state variables
    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'concentration in carbon units', minimum=0.0_rk)

    if (use_community_sinking) then
        call self%register_state_variable(self%id_dsnk, 'dsnk', 'mgC/m3', 'detritus sinking speed advector', minimum=1.0e-10_rk/sedy0)
        call self%register_diagnostic_variable(self%id_snkspd, 'snkspd', 'm/d', 'community detritus sinking speed', output=output_time_step_averaged)
    end if

    ! Register dependencies
    call self%register_dependency(self%id_temp,standard_variables%temperature)
    call self%register_state_dependency(self%id_no3, 'no3', 'mgC/m3', 'nitrate')
    call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
    call self%register_state_dependency(self%id_pho, 'pho', 'mgC/m3', 'phosphate')
    call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol/m3', 'oxygen')
    if (couple_co2) then
        call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
        call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
    end if

    ! optional routines ---------- !
    ! depth dependent sinking speed
    if (depth_dependent_sinking_speed) then
        call self%get_parameter( self%cz_OMsink,'cz_OMsink',    '1/day',      'detritus sinking rate increase per metre below surface', default=0.0_rk, scale_factor=1.0_rk/sedy0)        
        call self%register_dependency(self%id_depth,standard_variables%depth)
    end if
  end subroutine initialize

  subroutine do(self,_ARGUMENTS_DO_)
    class (type_ecosmo_organic_matter),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    real(rk) :: frem, remineralization
    real(rk) :: bact_lim
    real(rk) :: c, temp
    real(rk) :: no3, pho, nh4, oxy, dic, alk
    real(rk) :: bioom5, bioom6, bioom7
    real(rk) :: rhs_oxy, rhs, rhs_amm, rhs_nit, rhs_dic, rhs_alk, rhs_c, rhs_pho
    real(rk) :: dsnk, rhs_dsnk

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

    bioom5 = merge(5.0_rk, 0.0_rk, (oxy <= 0.0_rk) .and. (no3 > 0.0_rk)) ! anoxic and nitrate available: 5.0, else 0.0
    bioom6 = merge(1.0_rk, 0.0_rk, oxy > 0.0_rk) ! if oxic: 1.0, else 0.0
    bioom7 = merge(1.0_rk, 0.0_rk, (oxy <= 0.0_rk) .and. (no3 <= 0.0_rk)) ! if anoxic and nitrate not available: 1.0, else 0.0

    frem = self%remin * ( 1.0_rk + 20.0_rk * ( (temp*temp) / ( 169.0_rk + (temp*temp) ) ) )

    ! -------------------------------------------------------------------------
    ! Implicit Bacterial Nutrient Limitation
    ! -------------------------------------------------------------------------
    ! Remineralization is mediated by heterotrophic bacteria. Because DOM is 
    ! typically carbon-rich and nutrient-poor, bacteria must consume dissolved 
    ! inorganic nutrients (NH4, PO4) from the water column to satisfy their 
    ! cellular stoichiometry. In oligotrophic conditions, bacterial DOM 
    ! degradation becomes nutrient-limited. 
    !
    ! To capture this without an explicit bacterial state variable, the 
    ! base remineralization rate is modulated by a Michaelis-Menten 
    ! limitation term based on ambient NH4 and PO4. 
    !
    ! References: 
    ! - Aumont et al. (2015), PISCES-v2 (GMD 8, 2465-2513).
    ! - Letscher et al. (2015), (Biogeosciences, 12(1), 209-221).
    !
    ! VCY - 09/10/2026
    ! Added community sinking option for organic matter sinking speed.
    
    ! -------------------------------------------------------------------------
    if (use_bact_nutrient_limitation) then
        bact_lim = min( nh4 / (self%K_bact_NH4 + nh4), pho / (self%K_bact_PO4 + pho) )
        frem = frem * bact_lim
    end if

    remineralization = frem * c
    
    ! Organic matter change
    rhs_c = -remineralization
    _ADD_SOURCE_(self%id_c, rhs_c)

    if (use_community_sinking) then
        _GET_(self%id_dsnk, dsnk)
        rhs_dsnk = -remineralization * (dsnk / max(c, 1e-10_rk))  
        _ADD_SOURCE_(self%id_dsnk, rhs_dsnk)
        _SET_DIAGNOSTIC_(self%id_snkspd, (dsnk / max(c, 1e-10_rk)) * sedy0)
    end if

    ! Ammonium change
    rhs_amm = remineralization
    _ADD_SOURCE_(self%id_nh4, rhs_amm)

    ! Nitrate change
    rhs_nit = -remineralization * bioom5
    _ADD_SOURCE_(self%id_no3, rhs_nit )

    ! Phosphate change
    rhs_pho = remineralization
    _ADD_SOURCE_(self%id_pho, rhs_pho )

    ! Oxygen change
    rhs_oxy = -remineralization * (bioom6+bioom7) * 6.625 * Cmg_to_Cmmol * Cmmol_to_Nmmol
    _ADD_SOURCE_(self%id_oxy, rhs_oxy)

    if (couple_co2) then

      ! CO2 change
      rhs_dic = remineralization * Cmg_to_Cmmol
      _ADD_SOURCE_(self%id_dic, rhs_dic )

      ! Alkalinity change
      rhs = (rhs_amm - rhs_nit) * Cmg_to_Cmmol * Cmmol_to_Nmmol - 0.5_rk * rhs_oxy * (1._rk-bioom6)
      _ADD_SOURCE_(self%id_alk, rhs )

    end if

    _LOOP_END_

  end subroutine do


  subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
    class (type_ecosmo_organic_matter),intent(in) :: self
    _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

    real(rk) :: w_det
    real(rk) :: depth
    real(rk) :: det, dsnk
  
    _LOOP_BEGIN_

    if (depth_dependent_sinking_speed) then
        _GET_(self%id_depth, depth)
        depth = abs(depth) ! some host models may report depth as negative
        w_det = self%OMsink + self%cz_OMsink * depth
    else
        w_det = self%OMsink 
    end if

    if (use_community_sinking) then
        _GET_(self%id_dsnk, dsnk)
        _GET_(self%id_c, det)
        w_det = dsnk / max(det, 1e-10_rk)
        _ADD_VERTICAL_VELOCITY_(self%id_dsnk, -w_det)
    end if

    _ADD_VERTICAL_VELOCITY_(self%id_c, -w_det)

    _LOOP_END_
  end subroutine get_vertical_movement

  subroutine check_state(self,_ARGUMENTS_CHECK_STATE_)
    class (type_ecosmo_organic_matter),intent(in) :: self
    _DECLARE_ARGUMENTS_CHECK_STATE_
    real(rk) :: det, dsnk, meanspd
    
    _LOOP_BEGIN_
    if (use_community_sinking) then
        _GET_(self%id_c, det)
        _GET_(self%id_dsnk, dsnk)
        meanspd = dsnk / max(det, 1e-10_rk)
        if (meanspd < self%min_dsnk) _SET_(self%id_dsnk, max(det, 1e-10_rk) * self%min_dsnk)
        if (meanspd > self%max_dsnk) _SET_(self%id_dsnk, max(det, 1e-10_rk) * self%max_dsnk)
    end if
    _LOOP_END_
  end subroutine check_state

end module
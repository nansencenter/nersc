#include "fabm_driver.h"

module ecosmo_phytoplankton

use fabm_types
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_phytoplankton

! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_phytoplankton
!     Variable identifiers
    type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_sil, id_oxy
    type (type_state_variable_id)         :: id_c, id_chl
    type (type_state_variable_id)         :: id_alk, id_dic
    type (type_state_variable_id)         :: id_det, id_dom, id_opal !, id_caco3

    type (type_dependency_id)             :: id_temp, id_salt, id_par, id_parmean

    type (type_diagnostic_variable_id)    :: id_primprod, id_netpp !, id_pcal
    type (type_horizontal_dependency_id)  :: id_sfpar, id_meansfpar

    real(rk) :: MAXchl2nP, MINchl2nP 
    real(rk) :: alfaP !, aa
    real(rk) :: mu, m, m2, Km2
    real(rk) :: rNH4, rNO3, rPO4, rSi
    real(rk) :: Psink
    real(rk) :: psi
    real(rk) :: exulim, qexcr
    real(rk) :: silicate_flux_multiplier !, caco3_flux_multiplier
    real(rk) :: SiUptLim
    !real(rk) :: calcR

    logical  :: turn_on_additional_diagnostics ! activates additional diagnostics for model debugging
    logical  :: is_diatom ! phytoplankton is a diatom
    logical  :: is_cyano ! phytoplankton is a cyanobacteria
    !logical  :: is_coccolith ! phytoplankton is a coccolith

    real(rk) :: TctrlBG, TrefBG, Bg_fix 
    real(rk) :: nfixation_minimum_daily_par, bg_growth_minimum_daily_rad
    contains

!     Model procedures
    procedure :: initialize
    procedure :: do
!    procedure :: do_surface
!    procedure :: get_light_extinction

end type type_ecosmo_phytoplankton


type (type_bulk_standard_variable), parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)
type (type_bulk_standard_variable), parameter :: total_gpp = type_bulk_standard_variable(name='total_gpp',units='mgC/m3/d',aggregate_variable=.true.)
type (type_bulk_standard_variable), parameter :: total_npp = type_bulk_standard_variable(name='total_npp',units='mgC/m3/d',aggregate_variable=.true.)

!type (type_bulk_standard_variable), parameter :: calcite_production = type_bulk_standard_variable(name='calcite_production', units='mmolN m-3 s-1', aggregate_variable=.true.)
type (type_bulk_standard_variable), parameter :: pbiomass = type_bulk_standard_variable(name='pbiomass',units='mg/m^3',aggregate_variable=.true.)

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_phytoplankton), intent(inout),target  :: self
    integer,  intent(in) :: configunit
 
    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    
    call self%get_parameter( self%is_diatom, 'is_diatom', '', 'use silicate', default=.false. )
!    call self%get_parameter( self%is_cyano, 'is_cyano', '', 'use n-fixation', default=.false. )
!    call self%get_parameter( self%is_coccolith, 'is_coccolith', '', 'calcify', default=.false. )
    call self%get_parameter( self%turn_on_additional_diagnostics, 'turn_on_additional_diagnostics','','activates additional diagnostics for model debugging',default=.false.)

    call self%get_parameter( self%mu , 'mu',  '1/day', 'max growth rate', default=1.30_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%rNH4, 'rNH4', 'mmolN/m**3', 'NH4 half saturation constant', default=0.20_rk, scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg)
    call self%get_parameter( self%rNO3, 'rNO3', 'mmolN/m**3', 'NO3 half saturation constant', default=0.50_rk, scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg)
    call self%get_parameter( self%rPO4, 'rPO4', 'mmolP/m**3', 'PO4 half saturation constant', default=0.05_rk,  scale_factor=Pmmol_to_Cmmol*Cmmol_to_Cmg)
    call self%get_parameter( self%rSi, 'rSi', 'mmolSi/m**3','SiO2 half saturation', default=0.50_rk,  scale_factor=Simmol_to_Cmmol*Cmmol_to_Cmg)
    call self%get_parameter( self%m, 'm', '1/day', 'mortality rate', default=0.04_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%m2 , 'm2',         '1/day',      'P higher mortality rate',               default=0.0_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%Km2, 'Km2',         'mgC/m**3',      'half saturation for P higher mortality rate', default=300.0_rk)
    call self%get_parameter( self%Psink, 'Psink', 'm/day', 'phytoplankton sinking rate', default=0.0_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%MINchl2nP, 'MINchl2nP', 'mgChl/mmolN', 'minimum Chl to N ratio P', default=0.50_rk, scale_factor=Cmmol_to_Nmmol*Cmg_to_Cmmol)
    call self%get_parameter( self%MAXchl2nP, 'MAXchl2nP', 'mgChl/mmolN', 'maximum Chl to N ratio P', default=3.83_rk, scale_factor=Cmmol_to_Nmmol*Cmg_to_Cmmol)    
    call self%get_parameter( self%alfaP, 'alfaP', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve P', default=0.0393_rk, scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg )
    call self%get_parameter( self%psi , 'psi', 'm**3/mmolN', 'NH4 inhibition', default=3.0_rk,   scale_factor=Cmmol_to_Nmmol*Cmg_to_Cmmol )
!    call self%get_parameter( self%aa, 'aa', 'm**2/W', 'photosynthesis ef-cy', default=0.04_rk)
    call self%get_parameter( self%exulim,  'exulim',   '', 'fraction of GPP released as DOC (exudation) due to nutrient limiting conditions',  default=0.0_rk)
    call self%get_parameter( self%qexcr,  'qexcr',   '', 'fraction of GPP released as DOC (excretion) due to activity',  default=0.0_rk) 
    call self%get_parameter( self%SiUptLim,  'SiUptLim',   'mgC/m3', 'Stop Si uptake below this concentration',  default=80.0_rk)

    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'carbon', minimum=1.0e-7_rk, vertical_movement=-self%Psink ,initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg )
!    if (use_chl) then
        call self%register_state_variable(self%id_chl, 'chl', 'mgChl/m3', 'chlorophyll', minimum=1.0e-7_rk/27., vertical_movement=-self%Psink ,initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg/27.)
        call self%add_to_aggregate_variable(total_chlorophyll, self%id_chl)
        call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
        self%id_chl,scale_factor=light_att_chl,include_background=.true.)
    ! else
    !     call self%add_to_aggregate_variable(total_chlorophyll, self%id_c, scale_factor=1.0_rk/60.0_rk)
    !     call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
    !     self%id_c,scale_factor=light_att_phy,include_background=.true.)
    ! end if

    call self%register_diagnostic_variable(self%id_primprod,'primprod','mgC/m3/s', &
         'primary production rate', output=output_time_step_averaged)
    call self%add_to_aggregate_variable(total_gpp, self%id_primprod)

    call self%register_diagnostic_variable(self%id_netpp,'netpp','mgC/m**3/s', &
         'net primary production rate', output=output_time_step_averaged)
    call self%add_to_aggregate_variable(total_npp, self%id_netpp)

    call self%add_to_aggregate_variable(pbiomass, self%id_c)

    call self%register_state_dependency(self%id_no3, 'no3', 'mgC/m3', 'nitrate')
    call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
    call self%register_state_dependency(self%id_pho, 'pho', 'mgC/m3', 'phosphate')
    call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol/m3', 'oxygen')
    call self%register_state_dependency(self%id_det, 'det', 'mgC/m3', 'detritus')
    call self%register_state_dependency(self%id_dom, 'dom', 'mgC/m3', 'dom')
    call self%register_state_dependency(self%id_opal, 'opal', 'mgC/m3', 'opal')

    if (self%is_diatom) then
        call self%register_state_dependency(self%id_sil, 'sil', 'mgC/m3', 'silicate')
        self%silicate_flux_multiplier = 1.0_rk ! do not change
    else
        self%silicate_flux_multiplier = 0.0_rk ! do not change
    end if 
    call self%register_dependency(self%id_temp,standard_variables%temperature)
    call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
    call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)

    if (couple_co2) then
        call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
        call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
    end if

    ! if (self%is_cyano) then
    !     call self%register_dependency(self%id_sfpar,standard_variables%surface_downwelling_photosynthetic_radiative_flux)        
    !     call self%register_dependency(self%id_meansfpar,temporal_mean(self%id_sfpar,period=86400._rk,resolution=3600._rk))
    !     call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk))
    !     call self%get_parameter(self%nfixation_minimum_daily_par, 'nfixation_minimum_daily_par', 'nfixation minimum daily par', default=40.0_rk)
    !     call self%get_parameter(self%bg_growth_minimum_daily_rad, 'bg_growth_minimum_daily_rad', 'bg growth minimum daily rad', default=120.0_rk)     
    !     call self%get_parameter( self%TctrlBG,  'TctrlBG',    '1/degC',     'BG T control beta',               default=1.00_rk)
    !     call self%get_parameter( self%TrefBG,  'TrefBG',     'degC',       'BG reference temperature',        default=0.00_rk)
    ! end if

!     if (self%is_coccolith) then
!         call self%get_parameter( self%calcR, 'calcR', '-', 'Maximum calcification to organic carbon production RCaCO3', default=0.4_rk)
!         call self%register_diagnostic_variable(self%id_pcal, 'pcal', 'mmolN m-3 s-1', 'calcite production')
! !        call self%add_to_aggregate_variable(calcite_production, self%id_pcal)
!         call self%register_state_dependency(self%id_caco3, 'caco3', 'mmol/m3', 'calcite')
!         self%caco3_flux_multiplier = 1.0_rk ! do not change
!     else
!         self%caco3_flux_multiplier = 0.0_rk ! do not change
!     end if
        
end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)

    class (type_ecosmo_phytoplankton),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    real(rk) :: no3,nh4,pho,sil,oxy
    real(rk) :: det, dom
    real(rk) :: temp, salt, par
    real(rk) :: c, chl
    real(rk) :: p_loss, limit, prod, nutlimit
    real(rk) :: up_nh4, up_no3, up_pho, up_sil, t_sil, up_n, psi
    real(rk) :: highMort
    real(rk) :: Tdep ! temperature effect on P growth
    real(rk) :: blight !, aa
    real(rk) :: chl2c
    real(rk) :: rhs, rhs_amm, rhs_nit, rhs_oxy !, rhs_caco3
    real(rk) :: bioom6
    real(rk) :: dic,alk
    real(rk) :: N_or_Si_uptake, P_uptake
    real(rk) :: mean_surface_par, mean_par
    real(rk) :: exu
    !real(rk) :: Bg_fix

!    real(rk) :: Rstar

    _LOOP_BEGIN_

    ! Retrieve current (local) state variable values.
    _GET_(self%id_temp,temp)
    _GET_(self%id_salt,salt)
    _GET_(self%id_par,par)
    _GET_(self%id_no3,no3)
    _GET_(self%id_nh4,nh4)
    _GET_(self%id_pho,pho)
    _GET_(self%id_sil,sil)
    _GET_(self%id_oxy,oxy)
    _GET_(self%id_det,det)
    _GET_(self%id_dom,dom)
    _GET_(self%id_c,c)
!    if (use_chl) then
      _GET_(self%id_chl,chl)
!    end if
    if (couple_co2) then
        _GET_(self%id_dic,dic)
        _GET_(self%id_alk,alk)
    end if
    p_loss = max( sign( -1.0_rk, c - prevent_loss_P ), 0.0_rk ) 

   ! nutrient limitation factors
   ! k denotes half-saturation values
    up_nh4 = nh4/(self%rNH4 + nh4)
    up_no3 = no3/(self%rNO3 + no3) * exp(-self%psi * nh4)
    up_n = up_nh4 + up_no3
    up_pho = pho/(self%rPO4 + pho)
    t_sil = max(sil-self%SiUptLim,0.0_rk)
    up_sil = t_sil/(self%rSi + t_sil)

    ! temperature dependence
    Tdep = 1.0_rk
    ! if (self%is_cyano) then
    !     _GET_(self%id_parmean,mean_par)
    !     _GET_HORIZONTAL_(self%id_meansfpar,mean_surface_par)
    !     Tdep = 0.0_rk
    !     if ((salt<=10.0) .and. (mean_surface_par > self%bg_growth_minimum_daily_rad)) then
    !         Tdep = 1.0_rk/(1.0_rk + exp( self%TctrlBG * (self%TrefBG - temp )))
    !     end if
    ! end if

!    if (use_chl_in_PI_curve) then
        blight  = max( ((chl/c) * self%alfaP * par) / sqrt((self%mu * sedy0)**2 + (chl/c)**2 * self%alfaP**2 * (par**2)) ,0.0_rk)
    ! else
    !     blight = max(tanh(self%aa * par), 0.0_rk)
    ! end if

    if (self%is_diatom) then
        limit = Tdep * min(blight, up_n, up_pho, up_sil) ! limitation in maximum growth (range: 0 - 1)
        nutlimit = min(up_n, up_pho, up_sil)
    else
        limit = Tdep * min(blight, up_n, up_pho)
        nutlimit = min(up_n, up_pho)
    end if

    ! calculate exudation. Default: exulim=0, qexcr=0, thus ignored
    exu = min(1.0_rk,( ( 1.0_rk - nutlimit ) * self%exulim + self%qexcr )) * p_loss

    ! below are all the same primary production in seconds, 
        ! but are also stored as nutrient specific versions if in future some functionality added.
        ! nutrients rhs's later use their respective uptake/prod names 
    N_or_Si_uptake = self%mu * limit * c ! primary production in seconds
    P_uptake = self%mu * limit * c ! primary production in seconds
    prod = self%mu * limit * c ! primary production in seconds

    ! if (self%is_cyano .and. mean_par > self%nfixation_minimum_daily_par) then
    !     Bg_fix = Tdep * min(blight, up_pho) - limit
    !     limit = limit + Bg_fix
    ! else
    !     Bg_fix = 0.0_rk
    ! end if

    !** 
!    prod = self%mu * limit * c ! primary production in seconds
!    P_uptake = prod
    ! this includes the extra phosphate uptake due to n-fixation of cyanobacteria
    ! if phyto is not a cyano, this is essentially equal to N_or_Si_uptake calculated above
    ! but if there is n-fixation, this ensures the extra pho uptake unaccounted for
    ! as limit is modified within n-fixation if/else case above
    !**
    highMort = self%m2 * ( c / ( c + self%Km2 ) )
    _ADD_SOURCE_(self%id_c, (self%mu * limit * (1.0_rk - exu) - (self%m + highMort) * p_loss) * c )

!    if (use_chl) then
        ! chlorophyll-a to C change
        chl2c = self%MAXchl2nP * max(0.01,limit) * self%mu * sedy0 * c / (self%alfaP * par * chl)
        chl2c = max(self%MINchl2nP,chl2c)
        chl2c = min(self%MAXchl2nP,chl2c)

        rhs = self%mu * limit * chl2c * c * (1.0_rk - exu) - ( (self%m + highMort) * p_loss * chl )
        _ADD_SOURCE_(self%id_chl,rhs)
!    end if

    rhs_nit = -(up_no3+0.5d-10)/(up_n+1.0d-10) * N_or_Si_uptake
    _ADD_SOURCE_(self%id_no3, rhs_nit)

    rhs_amm = -(up_nh4+0.5d-10)/(up_n+1.0d-10) * N_or_Si_uptake
    _ADD_SOURCE_(self%id_nh4, rhs_amm)

    _ADD_SOURCE_(self%id_pho, -P_uptake) 

    _ADD_SOURCE_(self%id_sil, -N_or_Si_uptake * self%silicate_flux_multiplier) ! will be effective only for diatoms
    _ADD_SOURCE_(self%id_opal, (self%m + highMort) * p_loss * c * self%silicate_flux_multiplier) ! will be effective only for diatoms

    rhs_oxy = (6.625*up_nh4 + 8.125*up_no3+1.d-10)/(up_n+1.d-10) * prod * Cmg_to_Cmmol * Cmmol_to_Nmmol 
    _ADD_SOURCE_(self%id_oxy, rhs_oxy)

    if (oxy > 0) then
        bioom6 = 1.0_rk
    else
        bioom6 = 0.0_rk
    end if

    ! ! Calcite production. self%caco3_flux_multiplier ensures production only for coccoliths 
    ! Rstar = self%calcR * limit * max(0.0001,temp/(2.0 + temp)) &
    !             * max( 1.0, 0.5 * c * Cmg_to_Cmmol * Cmmol_to_Nmmol ) &
    !             * self%caco3_flux_multiplier

    ! rhs_caco3 = Rstar * ( (self%mu * limit - 0.5 * self%m * p_loss) * c )
    ! _ADD_SOURCE_(self%id_caco3, rhs_caco3)
    ! _SET_DIAGNOSTIC_(self%id_pcal, Rstar)

    if (couple_co2) then
        !_ADD_SOURCE_(self%id_dic, -prod * Cmg_to_Cmmol - rhs_caco3)
        _ADD_SOURCE_(self%id_dic, -prod * Cmg_to_Cmmol)
        !rhs = (rhs_amm - rhs_nit) * Cmg_to_Cmmol * Cmmol_to_Nmmol - 0.5_rk * rhs_oxy * (1._rk-bioom6) - rhs_caco3
        rhs = (rhs_amm - rhs_nit) * Cmg_to_Cmmol * Cmmol_to_Nmmol - 0.5_rk * rhs_oxy * (1._rk-bioom6)
        _ADD_SOURCE_(self%id_alk, rhs )
    end if

    _ADD_SOURCE_(self%id_det, (1._rk - frr) * (self%m+highMort) * c * p_loss)
    _ADD_SOURCE_(self%id_dom, (frr * (self%m+highMort) * c * p_loss) + (self%mu * limit * exu * c) )

    _SET_DIAGNOSTIC_(self%id_primprod, prod * sedy0)

    _SET_DIAGNOSTIC_(self%id_netpp, prod * sedy0 - c * 0.1_rk )   

    _LOOP_END_

end subroutine do

! subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
!     class (type_ecosmo_phytoplankton),intent(in) :: self
!     _DECLARE_ARGUMENTS_DO_SURFACE_
 
!     real(rk) :: o2flux, T, tr, S, o2sat, oxy
!     real(rk) :: no3flux, phoflux

!     _HORIZONTAL_LOOP_BEGIN_
 
!     _GET_(self%id_temp,T)
!     _GET_(self%id_salt,S)
!     _GET_(self%id_oxy,oxy)

!    ! Oxygen saturation micromol/liter__(Benson and Krause, 1984)
!     tr = 1.0_rk/(T + 273.15_rk)
!     o2sat= exp(- 135.90205_rk              &
!         + (1.575701d05 ) * tr               &
!         - (6.642308d07 ) * tr**2            &
!         + (1.243800d10) * tr**3            &
!         - (8.621949d11) * tr**4            &
!         - S*(0.017674_rk-10.754_rk*tr+2140.7_rk*tr**2)  )
 
!  !   o2flux = 5._rk/sedy0 * (o2sat - oxy)
!     o2flux = 1._rk/sedy0 * (o2sat - oxy)
 
!     _ADD_SURFACE_FLUX_(self%id_oxy,o2flux)

!     _HORIZONTAL_LOOP_END_

! end subroutine do_surface

! subroutine get_light_extinction(self,_ARGUMENTS_GET_EXTINCTION_)
!     class (type_ecosmo_phytoplankton), intent(in) :: self
!     _DECLARE_ARGUMENTS_GET_EXTINCTION_
 
!     real(rk)                     :: chl,c
!     real(rk)                     :: my_extinction
 
!     ! Enter spatial loops (if any)
!     _LOOP_BEGIN_
 
!     ! Retrieve current (local) state variable values.
 
! !    my_extinction = 0.0_rk
!     if (use_chl) then
!         _GET_(self%id_chl, chl)
! !        my_extinction = my_extinction + light_att_chl * chl
!         my_extinction = light_att_chl * chl 
!     else
!         _GET_(self%id_c, c)
! !        my_extinction = my_extinction + light_att_phy * c
!     end if
 
!     _SET_EXTINCTION_( my_extinction )

!     _LOOP_END_
 
! end subroutine get_light_extinction

end module
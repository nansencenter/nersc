#include "fabm_driver.h"

module ecosmo_coccolith

use fabm_types
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_coccolith

! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_coccolith
!     Variable identifiers
    type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_oxy
    type (type_state_variable_id)         :: id_c, id_chl
    type (type_state_variable_id)         :: id_alk, id_dic
    type (type_state_variable_id)         :: id_det, id_dom, id_opal , id_caco3

    type (type_dependency_id)             :: id_temp, id_salt, id_par, id_parmean

    type (type_diagnostic_variable_id)    :: id_primprod, id_netpp , id_pcal
    type (type_horizontal_dependency_id)  :: id_sfpar, id_meansfpar

    real(rk) :: MAXchl2nP, MINchl2nP 
    real(rk) :: alfaP !, aa
    real(rk) :: mu, m, m2, Km2
    real(rk) :: rNH4, rNO3, rPO4
    real(rk) :: Psink
    real(rk) :: psi
    real(rk) :: exulim, qexcr
    real(rk) :: calcR

!    logical  :: turn_on_additional_diagnostics ! activates additional diagnostics for model debugging

    contains

!     Model procedures
    procedure :: initialize
    procedure :: do

end type type_ecosmo_coccolith


type (type_bulk_standard_variable), parameter :: total_chlorophyll = type_bulk_standard_variable(name='total_chlorophyll',units='mg/m^3',aggregate_variable=.true.)
type (type_bulk_standard_variable), parameter :: total_gpp = type_bulk_standard_variable(name='total_gpp',units='mgC/m3/d',aggregate_variable=.true.)
type (type_bulk_standard_variable), parameter :: total_npp = type_bulk_standard_variable(name='total_npp',units='mgC/m3/d',aggregate_variable=.true.)
type (type_bulk_standard_variable), parameter :: calcite_production = type_bulk_standard_variable(name='calcite_production', units='mmolN m-3 s-1', aggregate_variable=.true.)
type (type_bulk_standard_variable), parameter :: pbiomass = type_bulk_standard_variable(name='pbiomass',units='mg/m^3',aggregate_variable=.true.)

contains
subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_coccolith), intent(inout),target  :: self
    integer,  intent(in) :: configunit
 
    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    
!    call self%get_parameter( self%turn_on_additional_diagnostics, 'turn_on_additional_diagnostics','','activates additional diagnostics for model debugging',default=.false.)

    call self%get_parameter( self%mu , 'mu',  '1/day', 'max growth rate', default=1.30_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%rNH4, 'rNH4', 'mmolN/m**3', 'NH4 half saturation constant', default=0.20_rk, scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg)
    call self%get_parameter( self%rNO3, 'rNO3', 'mmolN/m**3', 'NO3 half saturation constant', default=0.50_rk, scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg)
    call self%get_parameter( self%rPO4, 'rPO4', 'mmolP/m**3', 'PO4 half saturation constant', default=0.05_rk,  scale_factor=Pmmol_to_Cmmol*Cmmol_to_Cmg)
    call self%get_parameter( self%m, 'm', '1/day', 'mortality rate', default=0.04_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%m2 , 'm2',         '1/day',      'P higher mortality rate',               default=0.0_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%Km2, 'Km2',         'mgC/m**3',      'half saturation for P higher mortality rate', default=300.0_rk)
    call self%get_parameter( self%Psink, 'Psink', 'm/day', 'phytoplankton sinking rate', default=0.0_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%MINchl2nP, 'MINchl2nP', 'mgChl/mmolN', 'minimum Chl to N ratio P', default=0.50_rk, scale_factor=Cmmol_to_Nmmol*Cmg_to_Cmmol)
    call self%get_parameter( self%MAXchl2nP, 'MAXchl2nP', 'mgChl/mmolN', 'maximum Chl to N ratio P', default=3.83_rk, scale_factor=Cmmol_to_Nmmol*Cmg_to_Cmmol)    
    call self%get_parameter( self%alfaP, 'alfaP', 'mmolN m2/(mgChl day W)**-1', 'initial slope P-I curve P', default=0.0393_rk, scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg )
    call self%get_parameter( self%psi , 'psi', 'm**3/mmolN', 'NH4 inhibition', default=3.0_rk,   scale_factor=Cmmol_to_Nmmol*Cmg_to_Cmmol )
    call self%get_parameter( self%exulim,  'exulim',   '', 'fraction of GPP released as DOC (exudation) due to nutrient limiting conditions',  default=0.0_rk)
    call self%get_parameter( self%qexcr,  'qexcr',   '', 'fraction of GPP released as DOC (excretion) due to activity',  default=0.0_rk) 

    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'carbon', minimum=1.0e-7_rk, vertical_movement=-self%Psink ,initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg )
    call self%register_state_variable(self%id_chl, 'chl', 'mgChl/m3', 'chlorophyll', minimum=1.0e-7_rk/27., vertical_movement=-self%Psink ,initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg/27.)
    call self%add_to_aggregate_variable(total_chlorophyll, self%id_chl)
    call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
        self%id_chl,scale_factor=light_att_chl,include_background=.true.)

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

    call self%register_dependency(self%id_temp,standard_variables%temperature)
    call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
    call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)

    if (couple_co2) then
        call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
        call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
    end if

    call self%get_parameter( self%calcR, 'calcR', '-', 'Maximum calcification to organic carbon production RCaCO3', default=0.4_rk)
    call self%register_diagnostic_variable(self%id_pcal, 'pcal', 'mmolN m-3 s-1', 'calcite production')
    call self%add_to_aggregate_variable(calcite_production, self%id_pcal)
    call self%register_state_dependency(self%id_caco3, 'caco3', 'mmol/m3', 'calcite')
        
end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)

    class (type_ecosmo_coccolith),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_
    real(rk) :: no3,nh4,pho,oxy
    real(rk) :: det, dom
    real(rk) :: temp, salt, par
    real(rk) :: c, chl
    real(rk) :: p_loss, limit, prod, nutlimit
    real(rk) :: up_nh4, up_no3, up_pho, up_n, psi
    real(rk) :: highMort
    real(rk) :: Tdep ! temperature effect on P growth
    real(rk) :: blight !, aa
    real(rk) :: chl2c
    real(rk) :: rhs, rhs_amm, rhs_nit, rhs_oxy , rhs_caco3
    real(rk) :: bioom6
    real(rk) :: dic,alk
    real(rk) :: N_or_Si_uptake, P_uptake
    real(rk) :: mean_surface_par, mean_par
    real(rk) :: exu

    real(rk) :: Rstar

    _LOOP_BEGIN_

    ! Retrieve current (local) state variable values.
    _GET_(self%id_temp,temp)
    _GET_(self%id_salt,salt)
    _GET_(self%id_par,par)
    _GET_(self%id_no3,no3)
    _GET_(self%id_nh4,nh4)
    _GET_(self%id_pho,pho)
    _GET_(self%id_oxy,oxy)
    _GET_(self%id_det,det)
    _GET_(self%id_dom,dom)
    _GET_(self%id_c,c)
    _GET_(self%id_chl,chl)
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

    ! temperature dependence
    Tdep = 1.0_rk
    blight  = max( ((chl/c) * self%alfaP * par) / sqrt((self%mu * sedy0)**2 + (chl/c)**2 * self%alfaP**2 * (par**2)) ,0.0_rk)

    limit = Tdep * min(blight, up_n, up_pho)
    nutlimit = min(up_n, up_pho)

    ! calculate exudation. Default: exulim=0, qexcr=0, thus ignored
    exu = min(1.0_rk,( ( 1.0_rk - nutlimit ) * self%exulim + self%qexcr )) * p_loss

    ! below are all the same primary production in seconds, 
        ! but are also stored as nutrient specific versions if in future some functionality added.
        ! nutrients rhs's later use their respective uptake/prod names 
    N_or_Si_uptake = self%mu * limit * c ! primary production in seconds
    P_uptake = self%mu * limit * c ! primary production in seconds
    prod = self%mu * limit * c ! primary production in seconds

    highMort = self%m2 * ( c / ( c + self%Km2 ) )
    _ADD_SOURCE_(self%id_c, (self%mu * limit * (1.0_rk - exu) - (self%m + highMort) * p_loss) * c )

    ! chlorophyll-a to C change
    chl2c = self%MAXchl2nP * max(0.01,limit) * self%mu * sedy0 * c / (self%alfaP * par * chl)
    chl2c = max(self%MINchl2nP,chl2c)
    chl2c = min(self%MAXchl2nP,chl2c)

    rhs = self%mu * limit * chl2c * c * (1.0_rk - exu) - ( (self%m + highMort) * p_loss * chl )
    _ADD_SOURCE_(self%id_chl,rhs)

    rhs_nit = -(up_no3+0.5d-10)/(up_n+1.0d-10) * N_or_Si_uptake
    _ADD_SOURCE_(self%id_no3, rhs_nit)

    rhs_amm = -(up_nh4+0.5d-10)/(up_n+1.0d-10) * N_or_Si_uptake
    _ADD_SOURCE_(self%id_nh4, rhs_amm)

    _ADD_SOURCE_(self%id_pho, -P_uptake) 

    rhs_oxy = (6.625*up_nh4 + 8.125*up_no3+1.d-10)/(up_n+1.d-10) * prod * Cmg_to_Cmmol * Cmmol_to_Nmmol 
    _ADD_SOURCE_(self%id_oxy, rhs_oxy)

    if (oxy > 0) then
        bioom6 = 1.0_rk
    else
        bioom6 = 0.0_rk
    end if

    ! Calcite production. self%caco3_flux_multiplier ensures production only for coccoliths 
    Rstar = self%calcR * limit * max(0.0001,temp/(2.0 + temp)) &
                 * max( 1.0, 0.5 * c * Cmg_to_Cmmol * Cmmol_to_Nmmol ) 

    rhs_caco3 = Rstar * ( (self%mu * limit - 0.5 * (self%m + highMort) * p_loss) * c )
    _ADD_SOURCE_(self%id_caco3, rhs_caco3)
    _SET_DIAGNOSTIC_(self%id_pcal, Rstar)

    if (couple_co2) then
        _ADD_SOURCE_(self%id_dic, -prod * Cmg_to_Cmmol - rhs_caco3)
        rhs = (rhs_amm - rhs_nit) * Cmg_to_Cmmol * Cmmol_to_Nmmol - 0.5_rk * rhs_oxy * (1._rk-bioom6) - rhs_caco3
        _ADD_SOURCE_(self%id_alk, rhs )
    end if

    _ADD_SOURCE_(self%id_det, (1._rk - frr) * (self%m+highMort) * c * p_loss)
    _ADD_SOURCE_(self%id_dom, (frr * (self%m+highMort) * c * p_loss) + (self%mu * limit * exu * c) )

    _SET_DIAGNOSTIC_(self%id_primprod, prod * sedy0)

    _SET_DIAGNOSTIC_(self%id_netpp, prod * sedy0 - c * 0.1_rk )   

    _LOOP_END_

end subroutine do

end module
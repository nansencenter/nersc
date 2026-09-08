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
!
! Note on coding approach: The code is written in a way that the generic structure here is to be 
! utilized by fabm.yaml file.
! Hence, the code here includes if cases such as "if (self%is_diatom) then" to allow for flexibility in the model configuration.
! According to AI suggestion on the convention used here, cases like "if (self%is_diatom) then" are included.
! Here is AI reasoning on optimisation:
! In high-performance computing (HPC), we generally avoid if statements inside massive loops only if the condition changes from grid point to grid point (e.g., if (oxy > 0)),
! because this causes "branch mispredictions" that stall the CPU pipeline.
! However, self%is_diatom and self%is_calcifier are loop invariants. They are defined once during initialization and never change during the spatial _LOOP_BEGIN_ ... _LOOP_END_.
! Modern CPUs will perfectly predict this branch 100% of the time after the first grid point. The cost of evaluating the if statement becomes virtually zero.
!
! VCY - 01/06/2026
! Added temperature dependency on phytoplankton growth rate. 
! By default, temperature dependency is turned off in shared.F90 using use_temp_dependency_phy = .false.
! When set true in fabm.yaml, the model will use q10 parameter (default:1.0, thus inactive for safety).
! The user needs to set an appropriate q10 value in fabm.yaml.
!
! VCY - 22/07/2026
! Chl to biomass ratio was N-based. Modified the parameter names and values to be Chl to C ratio.
!
! VCY - 03/09/2026
! Added a different temperature dependency on calcifier growth rate after Fielding et al. 2013 (https://doi.org/10.4319%2Flo.2013.58.2.0663).
! It assumes much lower growth rates in cold temperatures.
!
! VCY - 08/09/2026 - General structural changes following PJW's suggestions (niva-ecosmo codebase):
!    - Refactored light limitation and Chl:C ratio calculation
!    - Implemented multiplicative light-nutrient limitation (limit = blight * nutlimit)
!    - Removed artificial 0.01 minimum for Chl:C ratio
!    - Refined Geider Chl:C synthesis ratio calculation
!    - Added "use_geider_PI_curve" to control method choice
! ------------------------------- !

module ecosmo_phy
    use fabm_types
    use fabm_expressions
    use ecosmo_shared
    implicit none
    private
    public type_ecosmo_phy


    type,extends(type_base_model), public  :: type_ecosmo_phy

        type (type_state_variable_id)         :: id_no3, id_nh4, id_pho, id_sil, id_oxy
        type (type_state_variable_id)         :: id_c, id_chl
        type (type_state_variable_id)         :: id_alk, id_dic
        type (type_state_variable_id)         :: id_det, id_dom, id_opal , id_caco3
        type (type_dependency_id)             :: id_temp, id_salt, id_par, id_parmean, id_Om_cal
        type (type_diagnostic_variable_id)    :: id_primprod, id_netpp , id_pcal

        real(rk) :: MAXchl2cP, MINchl2cP, alfaP, betaP
        real(rk) :: mu, m, m2, Km2, mu_act   
        real(rk) :: rNH4, rNO3, rPO4, rSi
        real(rk) :: Psink
        real(rk) :: psi
        real(rk) :: exulim, qexcr
        real(rk) :: SiUptLim
        real(rk) :: q10
        real(rk) :: calcR
        real(rk) :: Rain0, Kcalom, dissCmax, ndissC

        logical  :: turn_on_additional_diagnostics ! activates additional diagnostics for model debugging
        logical  :: is_diatom ! phytoplankton is a diatom, thus requires silicate
        logical  :: is_calcifier ! phytoplankton is a calcifier, e.g. coccolithophore

    contains
        procedure :: initialize
        procedure :: do
    end type type_ecosmo_phy

contains
    subroutine initialize(self,configunit)
        class (type_ecosmo_phy), intent(inout),target  :: self
        integer,  intent(in)                           :: configunit

        ! Register parameters
        call self%get_parameter( self%is_diatom,    'is_diatom',    '-',          'Diatom configuration', default=.false. )
        call self%get_parameter( self%is_calcifier, 'is_calcifier', '-',          'Calcifier configuration', default=.false. )
        call self%get_parameter( self%turn_on_additional_diagnostics, 'turn_on_additional_diagnostics','','activates additional diagnostics for model debugging',default=.false.)

        call self%get_parameter( self%mu ,          'mu',           '1/day',      'max growth rate', default=1.30_rk,  scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%rNH4,         'rNH4',         'mmolN/m**3', 'NH4 half saturation constant', default=0.20_rk, scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg)
        call self%get_parameter( self%rNO3,         'rNO3',         'mmolN/m**3', 'NO3 half saturation constant', default=0.50_rk, scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg)
        call self%get_parameter( self%rPO4,         'rPO4',         'mmolP/m**3', 'PO4 half saturation constant', default=0.05_rk,  scale_factor=Pmmol_to_Cmmol*Cmmol_to_Cmg)
        call self%get_parameter( self%rSi,          'rSi',          'mmolSi/m**3','SiO2 half saturation', default=0.50_rk,  scale_factor=Simmol_to_Cmmol*Cmmol_to_Cmg)
        call self%get_parameter( self%m,            'm',            '1/day',      'mortality rate', default=0.04_rk,  scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%m2 ,          'm2',           '1/day',      'quadratic mortality rate',               default=0.0_rk,  scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%Km2,          'Km2',          'mgC/m**3',   'quadratic loss half-sat.', default=300.0_rk)
        call self%get_parameter( self%Psink,        'Psink',        'm/day',      'phytoplankton sinking rate', default=0.0_rk, scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%MINchl2cP,    'MINchl2cP',    'mgChl/mgC',  'minimum Chl to C ratio P', default=0.0063_rk)
        call self%get_parameter( self%MAXchl2cP,    'MAXchl2cP',    'mgChl/mgC',  'maximum Chl to C ratio P', default=0.0370_rk)    
        call self%get_parameter( self%alfaP,        'alfaP',        'mgC m2/(mgChl day W)**-1', 'initial slope P-I curve P', default=4.225_rk, scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%betaP,        'betaP',        'mgC/(mgChl day W/m2)', 'photoinhibition parameter for P', default=0.0_rk, scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%psi ,         'psi',          'm**3/mmolN', 'NH4 inhibition', default=3.0_rk,   scale_factor=Cmmol_to_Nmmol*Cmg_to_Cmmol )
        call self%get_parameter( self%exulim,       'exulim',       '-',          'fraction of GPP released as DOC (exudation) due to nutrient limiting conditions',  default=0.0_rk)
        call self%get_parameter( self%qexcr,        'qexcr',        '-',          'fraction of GPP released as DOC (excretion) due to activity',  default=0.0_rk) 
        call self%get_parameter( self%SiUptLim,     'SiUptLim',     'mgC/m3',     'Stop Si uptake below this concentration',  default=80.0_rk)
        call self%get_parameter( self%q10,          'q10',          '-',          'Q_10 temperature coefficient', default=1.0_rk)
        ! ----------------------------------- !

        ! Register state variables
        call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'carbon', minimum=1.0e-7_rk, vertical_movement=-self%Psink ,initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg )
        call self%register_state_variable(self%id_chl, 'chl', 'mgChl/m3', 'chlorophyll', minimum=1.0e-7_rk/27., vertical_movement=-self%Psink ,initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg/27.)
        ! ----------------------------------- !

        ! Register aggregate variables
        ! light attenuation due to chlorophyll
        call self%add_to_aggregate_variable(standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux, &
                                            self%id_chl,scale_factor=light_att_chl,include_background=.true.)

        ! total chlorophyll for output
        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='Chla',units='mg/m^3',aggregate_variable=.true.), &
                                            self%id_chl,include_background=.true.)

        ! gross primary production
        call self%register_diagnostic_variable(self%id_primprod,'primprod','mgC/m3/s', &
            'primary production rate', output=output_time_step_averaged)
        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_gpp', units='mgC/m3/d',aggregate_variable=.true.), &
            self%id_primprod,include_background=.true.)

        ! estimate of net primary production, subtracting a constant respiration term (0.1/day) from GPP.
        call self%register_diagnostic_variable(self%id_netpp,'netpp','mgC/m**3/s', &
            'net primary production rate', output=output_time_step_averaged)
        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_npp', units='mgC/m3/d',aggregate_variable=.true.), &
            self%id_netpp,include_background=.true.)        

        ! total phytoplankton biomass for output
        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='pbiomass', units='mgC/m3',aggregate_variable=.true.), &
            self%id_c,include_background=.true.)
        ! ----------------------------------- !

        ! Register dependencies
        call self%register_state_dependency(self%id_no3, 'no3', 'mgC/m3', 'nitrate')
        call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
        call self%register_state_dependency(self%id_pho, 'pho', 'mgC/m3', 'phosphate')
        call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol/m3', 'oxygen')
        call self%register_state_dependency(self%id_det, 'det', 'mgC/m3', 'detritus')
        call self%register_state_dependency(self%id_dom, 'dom', 'mgC/m3', 'dom')
        if (self%is_diatom) then
            call self%register_state_dependency(self%id_sil, 'sil', 'mgC/m3', 'silicate')
            call self%register_state_dependency(self%id_opal, 'opal', 'mgC/m3', 'opal')
        end if 
        call self%register_dependency(self%id_temp,standard_variables%temperature)
        call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
        call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)

        if (couple_co2) then
            call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
            call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
            if (self%is_calcifier) then

                call self%get_parameter( self%Kcalom,       'Kcalom',       '-',          'half-saturation constant for calcifier Rain Ratio dependence on calcite saturation state', default=1.0_rk)
                call self%get_parameter( self%dissCmax,     'dissCmax',     '1/day',      'maximum specific dissolution rate', default=0.03_rk, scale_factor=1.0_rk/sedy0)
                call self%get_parameter( self%ndissC,       'ndissC',       '-',          'power of the dissolution law (Keir 1980)', default=2.22_rk)
                call self%get_parameter( self%Rain0,        'Rain0',       '-',          'maximum Rain Ratio (PIC:POC) within calcifiers', default=1.0_rk)
                call self%register_diagnostic_variable(self%id_pcal, 'pcal', 'mmolC m-3 s-1', 'calcite production')
                call self%register_dependency(self%id_Om_cal,    'Om_cal_target', '-','calcite saturation')
                call self%add_to_aggregate_variable(type_bulk_standard_variable(name='calcite_production', units='mmolC m-3 s-1', aggregate_variable=.true.), &
                    self%id_pcal,include_background=.true.)
                call self%register_state_dependency(self%id_caco3, 'caco3', 'mgC/m3', 'calcite')
            end if
        end if

    end subroutine initialize

    subroutine do(self,_ARGUMENTS_DO_)

        class (type_ecosmo_phy),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
        real(rk) :: no3,nh4,pho,sil,oxy
        real(rk) :: det, dom
        real(rk) :: temp, salt, par
        real(rk) :: c, chl
        real(rk) :: p_loss, limit, prod, nutlimit
        real(rk) :: up_nh4, up_no3, up_pho, up_sil, t_sil, up_n, psi
        real(rk) :: Tdep 
        real(rk) :: blight , growth, linear_mort, quad_mort
        real(rk) :: chl2c
        real(rk) :: rhs_c, rhs_amm, rhs_nit, rhs_oxy , rhs_caco3
        real(rk) :: rhs_pho, rhs_sil, rhs_opa, rhs_chl
        real(rk) :: rhs_alk, rhs_dic, rhs_det, rhs_dom
        real(rk) :: bioom6
        real(rk) :: dic,alk
        real(rk) :: exu
        real(rk) :: mu_act
        real(rk) :: exu_loss

        real(rk) :: Om_cal
        real(rk) :: RainR

        _LOOP_BEGIN_
        ! Not all instances will need rhs_caco3 (e.g. CO2 module without coccoliths))
        ! For coding simplicity, the RHSs include rhs_caco3. Here we assign rhs_caco3=0 for instances where it is not needed.
        rhs_caco3 = 0.0_rk

        ! Retrieve current (local) state variable values.
        _GET_(self%id_temp,temp)
        _GET_(self%id_salt,salt)
        _GET_(self%id_par,par)
        _GET_(self%id_no3,no3)
        _GET_(self%id_nh4,nh4)
        _GET_(self%id_pho,pho)
        if (self%is_diatom) then
            _GET_(self%id_sil,sil)
        else
            sil = 0.0_rk
        end if
        _GET_(self%id_oxy,oxy)
        _GET_(self%id_det,det)
        _GET_(self%id_dom,dom)
        _GET_(self%id_c,c)
        _GET_(self%id_chl,chl)
        if (couple_co2) then
            _GET_(self%id_dic,dic)
            _GET_(self%id_alk,alk)
            if (self%is_calcifier) then
                _GET_(self%id_Om_cal,Om_cal)
            else
                Om_cal = 0.0_rk
            end if
        end if

        ! nutrient limitation factors
        ! k denotes half-saturation values
        up_nh4 = nh4/(self%rNH4 + nh4)
        up_no3 = no3/(self%rNO3 + no3) * exp(-self%psi * nh4)
        up_n = up_nh4 + up_no3
        up_pho = pho/(self%rPO4 + pho)
        t_sil = max(sil-self%SiUptLim,0.0_rk)
        up_sil = t_sil/(self%rSi + t_sil)

        ! temperature dependence
        if (use_temp_dependency_phy) then
            if (self%is_calcifier) then
                ! after Fielding et al. 2013 (https://doi.org/10.4319%2Flo.2013.58.2.0663)
                Tdep = 0.1419_rk * (max(0.0_rk,temp))**0.8151_rk ! this assumes very low growth rates for cold temperatures unlike the other ones
            else
                ! the model assumes q10 is for every 10 degrees Celsius increase in temperature
                ! the reference temperature is 0 degree-C. Adjust your self%mu for this reference temperature
                Tdep = self%q10**(temp/10.0_rk)
            end if
        else
            Tdep = 1.0_rk
        end if
        mu_act = self%mu * Tdep ! actual max growth rate after temperature adjustment

        ! 1. Calculate Nutrient Limitation (Liebig Minimum)
        if (self%is_diatom) then
            nutlimit = min(up_n, up_pho, up_sil)
        else
            nutlimit = min(up_n, up_pho)
        end if

        ! 2. Calculate Light Limitation and Chl:C synthesis ratio
        if (use_geider_PI_curve) then
            if (par > 1.0e-8_rk .and. mu_act > 1.0e-8_rk) then
                ! Platt et al. (1980) function with photoinhibition
                blight = (1.0_rk - exp(-self%alfaP * par * (chl/c) / mu_act)) * exp(-self%betaP * par * (chl/c) / mu_act)
                
                ! Geider Chl:C synthesis ratio (uses light-limited growth, no artificial 0.01 cap)
                chl2c = self%MAXchl2cP * (blight * mu_act * c) / (self%alfaP * par * chl)
            else
                blight = 0.0_rk
                chl2c = self%MAXchl2cP ! L'Hopital's rule limit as PAR -> 0
            end if
            
            ! Multiplicative combination of light and nutrients
            limit = blight * nutlimit
            
        else
            ! Legacy Yumruktepe et al. (2022) formulation
            blight  = max( ((chl/c) * self%alfaP * par) / sqrt((mu_act)**2 + (chl/c)**2 * self%alfaP**2 * (par**2)) ,0.0_rk)
            
            ! Liebig minimum of all limitations
            limit = min(blight, nutlimit)
            
            ! Legacy Chl:C synthesis ratio (with 0.01 cap)
            chl2c = self%MAXchl2cP * max(0.01_rk, limit) * mu_act * c / max(self%alfaP * par * chl, 1.0e-10_rk)
        end if

        ! calculate exudation. Default: exulim=0, qexcr=0, thus ignored
        ! if on, still can be turned off for the pelagic community (if (c < prevent_loss_P))
        exu_loss = merge(1.0_rk, 0.0_rk, c > prevent_loss_P)
        exu = min(1.0_rk,( ( 1.0_rk - nutlimit ) * self%exulim + self%qexcr )) * exu_loss

        ! Clamp Chl:C ratio to allowed min/max bounds
        chl2c = max(self%MINchl2cP,chl2c)
        chl2c = min(self%MAXchl2cP,chl2c)

        ! primary production in seconds
        prod = mu_act * limit * c
        ! growth rate after exudation in seconds
        growth = prod * (1.0_rk - exu)

        ! mortality
        linear_mort = merge(self%m, 0.0_rk, c > prevent_loss_P)
        quad_mort   = merge(self%m2 * ( c / ( c + self%Km2 ) ), 0.0_rk, c > prevent_loss_P)
        p_loss = (linear_mort + quad_mort) * c
    
        ! phytoplankton biomass change in seconds
        rhs_c = growth - p_loss
        _ADD_SOURCE_(self%id_c, rhs_c )

        ! chlorophyll change in seconds
        rhs_chl = growth * chl2c - p_loss * (chl/c)
        _ADD_SOURCE_(self%id_chl,rhs_chl)

        ! nitrate change in seconds
        rhs_nit = -(up_no3+0.5e-10_rk)/(up_n+1.0e-10_rk) * prod
        _ADD_SOURCE_(self%id_no3, rhs_nit)

        ! ammonium change in seconds
        rhs_amm = -(up_nh4+0.5e-10_rk)/(up_n+1.0e-10_rk) * prod
        _ADD_SOURCE_(self%id_nh4, rhs_amm)

        ! phosphate change in seconds
        rhs_pho = -prod
        _ADD_SOURCE_(self%id_pho, rhs_pho) 

        if (self%is_diatom) then
            rhs_sil = -prod
            _ADD_SOURCE_(self%id_sil, rhs_sil)

            rhs_opa = p_loss
            _ADD_SOURCE_(self%id_opal, rhs_opa)
        end if

        ! oxygen change in seconds
        rhs_oxy = (6.625_rk*up_nh4 + 8.125_rk*up_no3+1.0e-10_rk)/(up_n+1.0e-10_rk) * prod * Cmg_to_Cmmol * Cmmol_to_Nmmol 
        _ADD_SOURCE_(self%id_oxy, rhs_oxy)

        bioom6 = merge(1.0_rk, 0.0_rk, oxy > 0.0_rk)
        ! Calcite production. self%caco3_flux_multiplier ensures production only for coccoliths 

        if (self%is_calcifier) then
            ! The calcite formulation deviates from the ARC MFC ECOSMO code following: 
            ! https://github.com/NIVANorge/niva-ecosmo/blob/main/src/niva_ecosmo.F90
            ! First we calculate the calcifier 'rain ratio', i.e. the ratio PIC:POC
            ! within coccolithophores as a function of Omega(calcite), based on experimental data
            ! (Gehlen et al., 2007; Zondervan et al., 2002).
            RainR = self%Rain0 * max(0._rk, (Om_cal-1._rk)/(Om_cal-1._rk+self%Kcalom))
            RainR = RainR * (max(temp, 0.0_rk)/(2._rk+max(temp, 0.0_rk)))
            RainR = max( RainR * limit, 0.005_rk) 

            !Next we use the rain ratio to calculate fluxes to the detrital calcite pool
            !arising from particulate fractions of coccolith mortality
            rhs_caco3 = RainR * max(0.0_rk, prod - 0.5 * p_loss ) ! check later, assumption: production shouldn't be negative
            _ADD_SOURCE_(self%id_caco3, rhs_caco3)
            _SET_DIAGNOSTIC_(self%id_pcal, RainR)
        end if

        if (couple_co2) then
            ! Carbon dioxide changes in seconds
            rhs_dic = ( -prod - rhs_caco3 ) * Cmg_to_Cmmol
            _ADD_SOURCE_(self%id_dic, rhs_dic)

            ! Alkalinity changes in seconds
            rhs_alk = (rhs_amm - rhs_nit) * Cmg_to_Cmmol * Cmmol_to_Nmmol - 0.5_rk * rhs_oxy * (1._rk-bioom6) - 2.0_rk * rhs_caco3 * Cmg_to_Cmmol
            _ADD_SOURCE_(self%id_alk, rhs_alk )
        end if

        ! detritus changes in seconds
        rhs_det = p_loss * (1.0_rk - frr)
        _ADD_SOURCE_(self%id_det, rhs_det)

        ! DOM changes in seconds
        rhs_dom = (frr * p_loss) + (exu * prod)
        _ADD_SOURCE_(self%id_dom, rhs_dom) 

        ! Contribution to total primary production diagnostic in mgC/m3/s
        _SET_DIAGNOSTIC_(self%id_primprod, prod)

        ! Contribution to net primary production diagnostic in mgC/m3/s
        _SET_DIAGNOSTIC_(self%id_netpp, prod - (c * 0.1_rk)/sedy0) ! subtracting a constant respiration term (0.1/day) from GPP to get NPP estimate

        _LOOP_END_

    end subroutine do

end module
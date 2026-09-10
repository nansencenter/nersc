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
! Note on coding approach: The code is written in a way that the generic structure here is to be 
! utilized by fabm.yaml file.
! Hence, the code here includes if cases such as "if (self%is_diatom) then" to allow for flexibility in the model configuration.
! According to AI suggestion on the convention used here, cases like "if (self%is_diatom) then" are included.
!
! Here is AI reasoning on optimisation:
! In high-performance computing (HPC), we generally avoid if statements inside massive loops only if the condition changes from grid point to grid point (e.g., if (oxy > 0)),
! because this causes "branch mispredictions" that stall the CPU pipeline.
! However, self%is_diatom and self%is_calcifier are loop invariants. They are defined once during initialization and never change during the spatial _LOOP_BEGIN_ ... _LOOP_END_.
! Modern CPUs will perfectly predict this branch 100% of the time after the first grid point. The cost of evaluating the if statement becomes virtually zero.
!
! VCY - 08/09/2026
! Added optional parameter for zooplankton prey switching when set to True in fabm.yaml 
! (default: false). This allows the model to use the original prey preference-based grazing formulation.
!
! Added optional parameters related to sloppy feeding and waste routing paradigm, 
! This is an alternative to "classical" ECOSMO grazing scheme.
! By default, the model uses the "classical" scheme and this is deactivated in fabm.yaml
! To activate, set use_slp_egest_paradigm = .true. in fabm.yaml
! This approach follows the work of Steinberg & Landry, 2017 (https://doi.org/10.1146/annurev-marine-010814-015924) on sloppy feeding and waste routing
! The key difference is how egestion and exretion are handled and how they contribute to the carbon cycle.
! In the "classical" scheme, egestion and excretion are treated as separate processes that contribute to the carbon cycle.
! In this new scheme, egestion and excretion are treated as separate processes that contribute to the carbon cycle.
! The implementation is taking from PJW approach used in niva-ecosmo
!
! VCY - 09/09/2026
! Added explicit calcifier fraction (fc) to the grazing loop for prey zooplankton. This is only activated when use_virtual_calcite is True.
! 
! VCY - 09/09/2026
! Added community sinking option for zooplankton egestion and mortality
! If use_community_sinking is True, then the sinking speed is determined by the prey 
! (either detritus sinking speed or the prey detritus sinking speed).
! 
! ------------------------------- !

module ecosmo_zooplankton
    use fabm_types
    use fabm_expressions
    use ecosmo_shared
    implicit none
    private
    public type_ecosmo_zooplankton

    type,extends(type_base_model), public  :: type_ecosmo_zooplankton
        type (type_state_variable_id),  allocatable,dimension(:) :: id_preyc,id_preychl
        type (type_state_variable_id)         :: id_c, id_alk, id_dic
        type (type_state_variable_id)         :: id_nh4, id_no3, id_sil, id_pho, id_det, id_dom, id_oxy, id_opal , id_caco3
        type (type_dependency_id)             :: id_temp, id_salt, id_par
        type (type_dependency_id)             :: id_pcal
        type (type_diagnostic_variable_id)    :: id_secprod, id_totalsecprod

        type (type_state_variable_id)         :: id_dsnk
        real(rk)                              :: sinkD
        real(rk), allocatable                 :: prey_sinkD(:)

        type (type_horizontal_dependency_id) :: id_lat
        type (type_global_dependency_id)     :: id_yearday

        real(rk) :: Zsink
        integer  :: nprey
        real(rk) :: Minprey
        real(rk) :: Rg
        real(rk) :: m, m2, Km2
        real(rk) :: exc
        real(rk) :: KsLightDep, scaleRg
!        real(rk) :: zpr
        real(rk),allocatable :: pref(:), grz(:), gamma(:), fGslp(:), fdissC(:)
        real(rk) :: fAexc, fexcdom, freges, frmort
        real(rk),allocatable :: bio_loss(:),bio_loss_limit(:)

        logical, allocatable :: has_chl(:)
        logical, allocatable :: prey_is_detritus(:)
        logical, allocatable :: prey_is_diatom(:)
        logical, allocatable :: prey_is_calcifier(:)
        logical, allocatable :: prey_is_zoo(:)
        logical              :: any_calcifier
        logical  :: turn_on_additional_diagnostics ! activates additional diagnostics for model debugging
    contains
        procedure :: initialize
        procedure :: do
    end type type_ecosmo_zooplankton

contains

    subroutine initialize(self,configunit)
        class (type_ecosmo_zooplankton), intent(inout),target  :: self
        integer,                         intent(in)            :: configunit

        integer           :: iprey
        character(len=16) :: index
        logical           :: prey_is_zoo
        logical           :: prey_is_detritus
        logical           :: prey_is_diatom
        logical           :: prey_is_calcifier
        real(rk)          :: z_loss
        real(rk)          :: scaleRg, KsLightDep

        call self%get_parameter( self%Rg, 'Rg', 'mmolN/m**3', 'Zs, Zl half saturation',  default=0.50_rk,  scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg)
        call self%get_parameter( self%m, 'm', '1/day', 'Z mortality rate', default=0.10_rk,  scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%m2 , 'm2',         '1/day',      'Z higher mortality rate',               default=0.0_rk,  scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%Km2, 'Km2',         'mgC/m**3',      'half saturation for Z higher mortality rate', default=300.0_rk)
        call self%get_parameter( self%exc, 'exc', '1/day', 'Z excretion rate', default=0.06_rk,  scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%Zsink, 'Zsink', 'm/day', 'zooplankton sinking rate', default=0.0_rk, scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%KsLightDep,  'KsLightDep',   'W m-2', 'PAR half saturation for light dependent mortality',  default=1.0e-20_rk) ! do not make default=0.0
        call self%get_parameter( self%scaleRg,  'scaleRg',   'the minimum value for RgZl scaling for faster grazing while DVM active',  default=1.0_rk) ! default is off (full RgZl) ! remember that smaller Rg means faster feeding 
        call self%get_parameter( self%fAexc,   'fAexc',   '-', 'Fraction of absorbed food excreted (Active Respiration)', default=0.0_rk)
        call self%get_parameter( self%fexcdom, 'fexcdom', '-', 'Fraction of excretion routed to DOM', default=0.0_rk)
        call self%get_parameter( self%freges,  'freges',  '-', 'Fraction of egestion routed to DOM', default=0.0_rk)
        call self%get_parameter( self%frmort,  'frmort',  '-', 'Fraction of mortality routed to DOM', default=frr)
!        call self%get_parameter(self%zpr, 'zpr', '1/day', 'zpr_long_name_needed', default=0.001_rk, scale_factor=1.0_rk/sedy0)

        call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'carbon', minimum=1.0e-7_rk, vertical_movement=-self%Zsink ,initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg )

        ! Determine number of prey types.
        call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)

        ! Get prey-specific parameters.
        allocate(self%pref(self%nprey))
        allocate(self%has_chl(self%nprey))
        allocate(self%prey_is_zoo(self%nprey))
        allocate(self%prey_is_detritus(self%nprey))
        allocate(self%prey_is_diatom(self%nprey))
        allocate(self%prey_is_calcifier(self%nprey))
        allocate(self%bio_loss_limit(self%nprey))
        allocate(self%grz(self%nprey))
        allocate(self%gamma(self%nprey))
        allocate(self%fGslp(self%nprey))
        allocate(self%fdissC(self%nprey))
        allocate(self%prey_sinkD(self%nprey))
        do iprey=1,self%nprey
            write (index,'(i0)') iprey
            call self%get_parameter(self%pref(iprey),'pref'//trim(index),'-','relative affinity for prey type '//trim(index))
            call self%get_parameter(self%grz(iprey),'grz'//trim(index),'-','Grazing rate on prey'//trim(index), default=1.0_rk, scale_factor=1.0_rk/sedy0)
            call self%get_parameter(self%gamma(iprey),'gamma'//trim(index),'-','Assim. eff. on plankton '//trim(index), default=0.75_rk)
            call self%get_parameter(self%fGslp(iprey),'fGslp'//trim(index),'-','Fraction of prey '//trim(index)//' lost to sloppy feeding', default=0.0_rk)
            call self%get_parameter(self%fdissC(iprey),'fdissC'//trim(index),'-','Fraction of prey '//trim(index)//' calcite dissolved in gut', default=0.5_rk)
            call self%get_parameter(self%prey_sinkD(iprey), 'prey'//trim(index)//'_sinkD', 'm/d', 'sinking speed of unassimilated waste from prey '//trim(index), default=5.0_rk, scale_factor=1.0_rk/sedy0)
        end do

        call self%get_parameter(self%sinkD, 'sinkD', 'm/d', 'sinking speed of detritus from zooplankton mortality', default=5.0_rk, scale_factor=1.0_rk/sedy0)

        self%any_calcifier = .false.
        do iprey=1,self%nprey
            write (index,'(i0)') iprey
            call self%get_parameter(self%prey_is_zoo(iprey),'prey'//trim(index)//'_is_zoo','','prey type '//trim(index)//' is zooplankton',default=.false.)
            call self%get_parameter(self%prey_is_detritus(iprey),'prey'//trim(index)//'_is_detritus','','prey type '//trim(index)//' is detritus',default=.false.)
            call self%get_parameter(self%prey_is_diatom(iprey),'prey'//trim(index)//'_is_diatom','','prey type '//trim(index)//' is diatom',default=.false.)
            call self%get_parameter(self%prey_is_calcifier(iprey),'prey'//trim(index)//'_is_calcifier','','prey type '//trim(index)//' is calcifier',default=.false.)
            self%has_chl(iprey) = .false. ! below, phyto will get .true. if use_chl 

            self%has_chl(iprey) = .false.

            if (self%prey_is_zoo(iprey)) then
                self%bio_loss_limit(iprey) = prevent_loss_Z
            else if (self%prey_is_detritus(iprey)) then 
                self%bio_loss_limit(iprey) = 0.0_rk 
            else
                self%bio_loss_limit(iprey) = prevent_loss_P 
                self%has_chl(iprey) = .true.
                if (self%prey_is_calcifier(iprey)) then
                    self%any_calcifier = .true.
                end if
            end if   

        end do

        ! Get prey-specific coupling links.
        allocate(self%id_preyc(self%nprey))
        allocate(self%id_preychl(self%nprey))
        do iprey=1,self%nprey
            write (index,'(i0)') iprey
            call self%register_state_dependency(self%id_preyc(iprey),'prey'//trim(index)//'c','mgC/m3', 'prey '//trim(index)//' carbon concentration')
            if (self%has_chl(iprey)) then
                call self%register_state_dependency(self%id_preychl(iprey),'prey'//trim(index)//'chl','mgChl/m3', 'prey '//trim(index)//' chl concentration')
            end if
        end do

        ! Get dependencies 
        call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
        call self%register_state_dependency(self%id_pho, 'pho', 'mgC/m3', 'phosphate')
        call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol/m3', 'oxygen')
        call self%register_state_dependency(self%id_det, 'det', 'mgC/m3', 'detritus')
        call self%register_state_dependency(self%id_dom, 'dom', 'mgC/m3', 'dom')
        call self%register_state_dependency(self%id_opal,'opal', 'mgC/m3', 'opal')

        if (couple_co2) then
            call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
            call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
            if (self%any_calcifier) then ! default is false. Assign true for calcifier in fabm.yaml
                call self%register_dependency(self%id_pcal,'pcal','-','calcite production')
                call self%register_state_dependency(self%id_caco3, 'caco3','mmol m-3','calcite')
            end if
        end if
        call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_lat,standard_variables%latitude)
        call self%register_dependency(self%id_yearday,standard_variables%number_of_days_since_start_of_the_year)

        if (use_community_sinking) then
            call self%register_state_dependency(self%id_dsnk, 'dsnk', 'mgC/m3', 'detritus sinking advector')
        end if

        ! register total zooplankton biomass for output
        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='zbiomass', units='mgC/m3',aggregate_variable=.true.), &
            self%id_c,include_background=.true.)
        ! register secondary production
        call self%register_diagnostic_variable(self%id_secprod,'secprod','mgC/m**3/s', &
            'secondary production rate', output=output_time_step_averaged)
        call self%add_to_aggregate_variable(type_bulk_standard_variable(name='total_secprod', units='mgC/m3/d',aggregate_variable=.true.), &
            self%id_secprod,include_background=.true.)

end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)

    class (type_ecosmo_zooplankton),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    integer  :: iprey,istate
    real(rk) :: c, oxy 
    real(rk),dimension(self%nprey) :: preyc, preychl, bio_loss 
    real(rk),dimension(self%nprey) :: food_each, uptake_each, uptake_rate_each, pref, effective_prey
    real(rk) :: rhs_phy, rhs_chl
    real(rk) :: uptake, uptake_rate, food
    real(rk) :: rhs_z, rhs_oxy, rhs_dic, rhs_opal , rhs_caco3, rhs_det, rhs_dom, rhs_nut, rhs_alk
    real(rk) :: assimilated_uptake_rate, unassimilated_uptake_rate
    real(rk) :: bioom6, rhs_amm
    real(rk) :: highMort
    real(rk) :: pcal

    real(rk) :: latitude, yearday, declination, day_length
    real, parameter :: pi = 3.14159265358979323846_rk
    real(rk) :: inside_acos
    real(rk) :: scale_Rg
    real(rk) :: light_dep_mort, linear_mort, quad_mort, zpr,  mortality, excretion
    real(rk) :: par
    
    
    real(rk) :: denom, eps
    real(rk) :: total_grazing, grazing_on_detritus
    real(rk) :: z_loss
    real(rk) :: dsnk, rhs_dsnk, det_val, spd_det, spd, egestion, waste
    real(rk) :: sloppy, ingested, active_excretion, basal_excretion, total_excretion
    real(rk) :: total_sloppy, total_egestion, total_absorbed, nut_from_excretion
    eps = 1e-12_rk

    _LOOP_BEGIN_

    ! Retrieve current (local) state variable values.
    _GET_(self%id_c,c)
    _GET_(self%id_par,par)
    _GET_(self%id_oxy,oxy)

   ! THIS BLOCK CONCERNS ACTIVATING DIEL VERTICAL MIGRATION
   ! BY DEFAULT, BELOW IS CONFIGURED TO BE INEFFECTIVE
   ! DVM IS ACTIVATED IN YAML FILE, see parameter declarations above for appropriate DVM activation
   ! for scaleRg, KsLightDep 
    _GET_SURFACE_(self%id_lat,latitude) ! degN
    _GET_GLOBAL_(self%id_yearday,yearday) !decimal day of the year

    latitude = latitude * pi / 180.0_rk
    declination = 23.44_rk * pi / 180.0_rk * sin(2.0_rk * pi / 365.0_rk * (yearday - 81.0_rk))
    inside_acos = max( -1.0_rk, min( 1.0_rk,-tan(latitude) * tan(declination) ) )

    day_length = 24.0_rk / pi * acos(inside_acos)
    day_length = max(0.0_rk, min(24.0_rk, day_length))

    scale_Rg = max( self%scaleRg , min(1.0_rk,(24.0_rk - day_length)/24.0_rk) )
    ! light dependent mortality multiplier
    if (self%KsLightDep < 1.0e-18_rk) then 
        light_dep_mort = 1.0_rk ! light-dependent mortality is off by default
    else
        light_dep_mort = par / (par + self%KsLightDep) ! assumes at low light, mortality decreases
    end if
    ! ------

    pcal = 0.0_rk
    if (couple_co2 .and. self%any_calcifier) then
            _GET_(self%id_pcal,pcal)
    end if  
    ! Get prey concentrations (loop)
    ! Notice that prey is not available below certain concentrations defined in shared.F90
    do iprey=1,self%nprey
        _GET_(self%id_preyc(iprey), preyc(iprey))
        if (self%has_chl(iprey)) then
            _GET_(self%id_preychl(iprey), preychl(iprey))
        end if

        bio_loss(iprey) = max(sign(-1.0_rk,preyc(iprey)-self%bio_loss_limit(iprey)),0.0_rk)
    end do

    pref = self%pref


    food_each = self%pref * preyc
    food = sum(food_each)

    if (use_prey_switching) then
        do iprey = 1, self%nprey
            effective_prey(iprey) = food_each(iprey) * preyc(iprey) / (food + eps)
        end do
        food = sum(effective_prey)
        denom = (self%Rg * scale_Rg)**2.0_rk + food**2.0_rk + eps
        
        do iprey = 1, self%nprey
            uptake_rate_each(iprey) = c * self%grz(iprey) * effective_prey(iprey) * preyc(iprey) / denom
            ! threshold applied ONLY here
            uptake_rate_each(iprey) = uptake_rate_each(iprey) * bio_loss(iprey)
        end do
    else
        denom = (self%Rg * scale_Rg)**2.0_rk + food**2.0_rk + eps
        
        do iprey = 1, self%nprey
            uptake_rate_each(iprey) = c * self%grz(iprey) * food_each(iprey) * preyc(iprey) / denom
            ! threshold applied ONLY here
            uptake_rate_each(iprey) = uptake_rate_each(iprey) * bio_loss(iprey)
        end do
    end if

    total_grazing = sum(uptake_rate_each)
    assimilated_uptake_rate = sum( self%gamma * uptake_rate_each)
    unassimilated_uptake_rate = total_grazing - assimilated_uptake_rate

    ! Apply/Calculate prey specific rates.
    rhs_opal = 0.0_rk
    rhs_caco3 = 0.0_rk
    grazing_on_detritus = 0.0_rk
    do iprey=1,self%nprey
        rhs_phy = 0.0_rk
        rhs_chl = 0.0_rk
        ! grazing rate on each prey type (mg C/m3/s)       
        rhs_phy = -uptake_rate_each(iprey)

        if (.not. self%prey_is_detritus(iprey)) then
            _ADD_SOURCE_(self%id_preyc(iprey), rhs_phy)
        else
            grazing_on_detritus = grazing_on_detritus + uptake_rate_each(iprey)
        end if

        ! grazing rate on each prey type chlorophyll-a (mg C/m3/s)
        if (self%has_chl(iprey)) then
            rhs_chl = rhs_phy * preychl(iprey) / max(preyc(iprey), 1e-12_rk) 
            _ADD_SOURCE_(self%id_preychl(iprey), rhs_chl )
        end if

        ! additive input to total rhs for opal
        if (self%prey_is_diatom(iprey)) then 
            rhs_opal = rhs_opal + uptake_rate_each(iprey)
        end if

        ! additive input to total rhs for caco3
        if (couple_co2 .and. self%prey_is_calcifier(iprey)) then
            if (use_virtual_calcite) then
                if (use_slp_egest_paradigm) then
                    rhs_caco3 = rhs_caco3 + (uptake_rate_each(iprey) * pcal) * (1.0_rk - (self%fdissC(iprey) * (1.0_rk - self%fGslp(iprey))))
                else
                    rhs_caco3 = rhs_caco3 + (uptake_rate_each(iprey) * pcal) * (1.0_rk - self%fdissC(iprey))
                end if
            else
                rhs_caco3 = rhs_caco3 + (0.5_rk * uptake_rate_each(iprey))  * pcal 
            end if
        end if
    end do

    ! Below are the rates that do not require a loop like that above
    ! When zooplankton biomass is below a certain threshold defined in shared.F90, zooplankton loss prevented.
    z_loss = merge(1.0_rk, 0.0_rk, c >= prevent_loss_Z)

    if (use_community_sinking) then
        _GET_(self%id_dsnk, dsnk)
        det_val = 0.0_rk
        do iprey = 1, self%nprey
            if (self%prey_is_detritus(iprey)) det_val = preyc(iprey)
        end do
        spd_det = dsnk / max(det_val, 1e-10_rk)
        rhs_dsnk = 0.0_rk
    end if

    if (use_slp_egest_paradigm) then
        ! --- NEW STEINBERG & LANDRY PARADIGM ---
        total_sloppy = 0.0_rk
        total_egestion = 0.0_rk
        total_absorbed = 0.0_rk
        
        do iprey = 1, self%nprey
            uptake = uptake_rate_each(iprey)
            sloppy = uptake * self%fGslp(iprey)
            ingested = uptake - sloppy
            egestion = ingested * (1.0_rk - self%gamma(iprey))
            
            total_sloppy = total_sloppy + sloppy
            total_egestion = total_egestion + egestion
            total_absorbed = total_absorbed + (ingested * self%gamma(iprey))
            
            if (use_community_sinking) then
                ! Explicitly exclude DOM (sloppy feeding and dissolved egestion)
                spd = merge(spd_det, self%prey_sinkD(iprey), self%prey_is_detritus(iprey)) ! if the prey is detritus, speed is the speed of detritus. If the prey is plankton, then the speed is that of prey detritus.
                rhs_dsnk = rhs_dsnk + egestion * (1.0_rk - self%freges) * spd
            end if
        end do
        
        ! Excretion
        basal_excretion  = self%exc * c * z_loss
        active_excretion = total_absorbed * self%fAexc
        total_excretion  = basal_excretion + active_excretion
        
        ! Mortality
        linear_mort = self%m * max(0.5_rk, light_dep_mort)
        quad_mort   = self%m2 * light_dep_mort * (c / (c + self%Km2))
        mortality   = (linear_mort + quad_mort) * c * z_loss
        
        ! Bulk Fluxes
        rhs_z   = total_absorbed - total_excretion - mortality
        rhs_nut = total_excretion * (1.0_rk - self%fexcdom)
        rhs_dom = total_sloppy + (total_excretion * self%fexcdom) + (total_egestion * self%freges) + (mortality * self%frmort)
        rhs_det = (total_egestion * (1.0_rk - self%freges)) + (mortality * (1.0_rk - self%frmort)) - grazing_on_detritus
        
        ! Store mineralized excretion for Oxygen/CO2 coupling
        nut_from_excretion = rhs_nut 
        
    else
        ! --- LEGACY ECOSMO PARADIGM ---
        ! Mortality
        linear_mort = self%m * max(0.5_rk, light_dep_mort)
        quad_mort   = self%m2 * light_dep_mort * (c / (c + self%Km2))
        mortality   = (linear_mort + quad_mort) * c * z_loss
        
        ! Excretion
        excretion = self%exc * c * z_loss
        
        ! Bulk Fluxes
        rhs_z   = assimilated_uptake_rate - mortality - excretion
        rhs_nut = excretion
        rhs_dom = frr * (unassimilated_uptake_rate + mortality)
        rhs_det = (1.0_rk - frr) * (unassimilated_uptake_rate + mortality) - grazing_on_detritus
        
        ! Store mineralized excretion for Oxygen/CO2 coupling
        nut_from_excretion = rhs_nut 
        
        if (use_community_sinking) then
            do iprey = 1, self%nprey
                spd = merge(spd_det, self%prey_sinkD(iprey), self%prey_is_detritus(iprey))  ! if the prey is detritus, speed is the speed of detritus. If the prey is plankton, then the speed is that of prey detritus.
                rhs_dsnk = rhs_dsnk + uptake_rate_each(iprey) * (1.0_rk - self%gamma(iprey)) * (1.0_rk - self%frmort) * spd
            end do
        end if
    end if
    
    if (use_community_sinking) then
        rhs_dsnk = rhs_dsnk + (mortality * (1.0_rk - self%frmort)) * self%sinkD
        rhs_dsnk = rhs_dsnk - grazing_on_detritus * spd_det
        _ADD_SOURCE_(self%id_dsnk, rhs_dsnk)
    end if

    ! Add bulk fluxes to FABM sources (Common to both paradigms)
    _ADD_SOURCE_(self%id_c, rhs_z)
    _ADD_SOURCE_(self%id_nh4, rhs_nut)
    _ADD_SOURCE_(self%id_pho, rhs_nut)
    _ADD_SOURCE_(self%id_dom, rhs_dom)
    _ADD_SOURCE_(self%id_det, rhs_det)

    ! Opal change. Since rhs_opal accounts for whether the prey is silicifier or not (calculated above), we don't need a switch here.
    _ADD_SOURCE_(self%id_opal, rhs_opal )

    ! Oxygen change.
    bioom6 = merge(1.0_rk, 0.0_rk, oxy > 0.0_rk) 
    rhs_oxy = -bioom6 * 6.625_rk * nut_from_excretion * Cmg_to_Cmmol * Cmmol_to_Nmmol 
    _ADD_SOURCE_(self%id_oxy, rhs_oxy)

    ! CaCO3 change. Since rhs_caco3 accounts for whether the prey is calcifier or not (calculated above), we don't need a switch here.
    _ADD_SOURCE_(self%id_caco3, rhs_caco3 )

    if (couple_co2) then
        ! CO2 change
        rhs_dic = nut_from_excretion * Cmg_to_Cmmol  
        ! Alkalinity change
        rhs_alk = nut_from_excretion * Cmg_to_Cmmol * Cmmol_to_Nmmol 
        
        if (use_virtual_calcite) then
            ! Deduct the DIC/ALK now to balance the carbon budget
            ! If coccolithophores were explicitly defined, the feedback to DIC/ALK is already handled within phy.F90
            rhs_dic = rhs_dic - (rhs_caco3 * Cmg_to_Cmmol)
            rhs_alk = rhs_alk - (2.0_rk * rhs_caco3 * Cmg_to_Cmmol)
        end if

        _ADD_SOURCE_(self%id_dic, rhs_dic )
        _ADD_SOURCE_(self%id_alk, rhs_alk)
    end if

    ! Export diagnostic variables
    _SET_DIAGNOSTIC_(self%id_secprod, assimilated_uptake_rate)

    _LOOP_END_

end subroutine do

end module
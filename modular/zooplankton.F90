#include "fabm_driver.h"

module ecosmo_zooplankton

use fabm_types
!use fabm_particle
use fabm_expressions
use ecosmo_shared

implicit none

private

!PUBLIC MEMBER FUNCTIONS:
public type_ecosmo_zooplankton

! !PUBLIC DERIVED TYPES:
type,extends(type_base_model), public  :: type_ecosmo_zooplankton
!type,extends(type_particle_model), public  :: type_ecosmo_zooplankton
!     Variable identifiers
!    type (type_model_id),         allocatable,dimension(:) :: id_prey
!    type (type_dependency_id),    allocatable,dimension(:) :: id_preyc,id_preychl
    type (type_state_variable_id),  allocatable,dimension(:) :: id_preyc,id_preychl

    type (type_state_variable_id)         :: id_c, id_chl
    type (type_state_variable_id)         :: id_alk, id_dic
    type (type_state_variable_id)         :: id_nh4, id_no3, id_sil, id_pho, id_det, id_dom, id_oxy, id_opal !, id_caco3

    type (type_dependency_id)             :: id_temp, id_salt, id_par
!    type (type_dependency_id)         :: id_pcal

    type (type_diagnostic_variable_id)    :: id_secprod

    type (type_horizontal_dependency_id)  :: id_sfpar

    type (type_horizontal_dependency_id) :: id_lat
    type (type_global_dependency_id)     :: id_yearday

    real(rk) :: Zsink
    integer  :: nprey
    real(rk) :: Minprey
    real(rk) :: grzP, grzZ  
    real(rk) :: Rg
    real(rk) :: m, m2, Km2
    real(rk) :: exc
    real(rk) :: gamma
    real(rk) :: KsLightDep, scaleRg
    real(rk),allocatable :: pref(:), opal_multiplier(:), grz(:) !, caco3_multiplier(:)
    real(rk),allocatable :: bio_loss(:),bio_loss_limit(:)
    !real(rk),allocatable :: caco3_loss(:)
    logical, allocatable :: has_chl(:)
    logical  :: turn_on_additional_diagnostics ! activates additional diagnostics for model debugging
    logical :: is_migrator   
    real(rk) :: swimspd
    contains

!     Model procedures
    procedure :: initialize
    procedure :: do
!    procedure :: get_vertical_movement
end type type_ecosmo_zooplankton

type (type_bulk_standard_variable), parameter :: total_secp = type_bulk_standard_variable(name='total_secp',units='mgC/m3/s',aggregate_variable=.true.)
type (type_bulk_standard_variable), parameter :: zbiomass = type_bulk_standard_variable(name='zbiomass',units='mg/m^3',aggregate_variable=.true.)
contains

subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_zooplankton), intent(inout),target  :: self
    integer,  intent(in) :: configunit
    ! !LOCAL VARIABLES:
    integer           :: iprey
    character(len=16) :: index
    logical           :: prey_is_not_phyto
    logical           :: prey_is_diatom
    logical           :: prey_is_coccolith
    real(rk)          :: z_loss
    real(rk)          :: scaleRg, KsLightDep

    !
    ! !REVISION HISTORY
    !
    !  Veli Çağlar Yumruktepe:
    !       XXX
    call self%get_parameter( self%is_migrator, 'is_migrator', '', 'perform_diel_vertical_migration', default=.false. )
    call self%get_parameter( self%swimspd, 'swimspd', 'm/day', 'zooplankton swimming speed', default=0.0_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%grzP, 'grzP', '1/day', 'Grazing rate on P', scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%grzZ, 'grzZ', '1/day', 'Grazing rate on Z', default=0.50_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%Rg, 'Rg', 'mmolN/m**3', 'Zs, Zl half saturation',  default=0.50_rk,  scale_factor=Nmmol_to_Cmmol*Cmmol_to_Cmg)
    call self%get_parameter( self%m, 'm', '1/day', 'Z mortality rate', default=0.10_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%m2 , 'm2',         '1/day',      'Z higher mortality rate',               default=0.0_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%Km2, 'Km2',         'mgC/m**3',      'half saturation for Z higher mortality rate', default=300.0_rk)
    call self%get_parameter( self%exc, 'exc', '1/day', 'Z excretion rate', default=0.06_rk,  scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%gamma, 'gamma', '', 'Z assim. eff. on plankton', default=0.75_rk)
    call self%get_parameter( self%Zsink, 'Zsink', 'm/day', 'zooplankton sinking rate', default=0.0_rk, scale_factor=1.0_rk/sedy0)
    call self%get_parameter( self%KsLightDep,  'KsLightDep',   'W m-2', 'PAR half saturation for light dependent mortality',  default=1.0e-20_rk) ! do not make default=0.0
    call self%get_parameter( self%scaleRg,  'scaleRg',   'the minimum value for RgZl scaling for faster grazing while DVM active',  default=1.0_rk) ! default is off (full RgZl) ! remember that smaller Rg means faster feeding 

!    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'carbon', minimum=1.0e-7_rk, vertical_movement=-self%Zsink ,initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg )
    call self%register_state_variable(self%id_c, 'c', 'mgC/m3', 'carbon', minimum=1.0e-7_rk, initial_value=1e-4_rk*Nmmol_to_Cmmol*Cmmol_to_Cmg )

    ! Determine number of prey types.
    call self%get_parameter(self%nprey,'nprey','','number of prey types',default=0)
    ! Get prey-specific parameters.
    allocate(self%pref(self%nprey))
    allocate(self%has_chl(self%nprey))
    allocate(self%opal_multiplier(self%nprey))
    ! if (use_calcifier) then
    !     allocate(self%caco3_multiplier(self%nprey))
    ! end if
    do iprey=1,self%nprey
        write (index,'(i0)') iprey
        call self%get_parameter(self%pref(iprey),'pref'//trim(index),'-','relative affinity for prey type '//trim(index))
    end do

    allocate(self%bio_loss_limit(self%nprey))
    allocate(self%grz(self%nprey))
    do iprey=1,self%nprey
        write (index,'(i0)') iprey
        call self%get_parameter(prey_is_not_phyto,'prey'//trim(index)//'_is_not_phyto','','prey type '//trim(index)//' is not phytoplankton',default=.false.)
        call self%get_parameter(prey_is_diatom,'prey'//trim(index)//'_is_diatom','','prey type '//trim(index)//' is diatom',default=.false.)
        ! if (use_calcifier) then
        !     call self%get_parameter(prey_is_coccolith,'prey'//trim(index)//'_is_coccolith','','prey type '//trim(index)//' is coccolith',default=.false.)
        ! end if
        self%has_chl(iprey) = .false. ! below, phyto will get .true. if use_chl 

        if (prey_is_not_phyto) then ! true for zoo and det, assign true in fabm.yaml
           self%bio_loss_limit(iprey) = prevent_loss_Z
           self%grz(iprey) = self%grzZ
        else ! false (default) for phyto, no need to assign in fabm.yaml
           self%grz(iprey) = self%grzP
           self%bio_loss_limit(iprey) = prevent_loss_P
        !    if (use_chl) then
                 self%has_chl(iprey) = .true.
        !    end if
        end if

        self%opal_multiplier(iprey) = 0.0_rk
        if (prey_is_diatom) then ! default is false. Assign true for diatoms in fabm.yaml
            self%opal_multiplier(iprey) = 1.0_rk
        end if
        ! if (use_calcifier) then
        !     self%caco3_multiplier(iprey) = 0.0_rk
        !     self%caco3_multiplier(iprey) = 1.0_rk
        !     call self%register_dependency(self%id_pcal,'pcal','-','calcite production')
        ! end if
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

    call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
    call self%register_state_dependency(self%id_pho, 'pho', 'mgC/m3', 'phosphate')
    call self%register_state_dependency(self%id_oxy, 'oxy', 'mmol/m3', 'oxygen')
    call self%register_state_dependency(self%id_det, 'det', 'mgC/m3', 'detritus')
    call self%register_state_dependency(self%id_dom, 'dom', 'mgC/m3', 'dom')
    call self%register_state_dependency(self%id_opal,'opal', 'mgC/m3', 'opal')


    if (couple_co2) then
        call self%register_state_dependency(self%id_dic, 'dic','mmol m-3','dic budget')
        call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
        ! if (use_calcifier) then
        !     call self%register_state_dependency(self%id_caco3, 'caco3','mmol m-3','calcite')
        ! end if
    end if

    call self%register_dependency(self%id_sfpar,standard_variables%surface_downwelling_photosynthetic_radiative_flux)        
    call self%register_dependency(self%id_par,standard_variables%downwelling_photosynthetic_radiative_flux)
    call self%register_dependency(self%id_lat,standard_variables%latitude)
    call self%register_dependency(self%id_yearday,standard_variables%number_of_days_since_start_of_the_year)

    call self%add_to_aggregate_variable(zbiomass, self%id_c)
end subroutine initialize

subroutine do(self,_ARGUMENTS_DO_)

    class (type_ecosmo_zooplankton),intent(in) :: self
    _DECLARE_ARGUMENTS_DO_

    integer  :: iprey,istate
    real(rk) :: c, oxy !, caco3
    real(rk),dimension(self%nprey) :: preyc, preychl, bio_loss !, caco3_loss
    real(rk),dimension(self%nprey) :: food_each, uptake_each, uptake_rate_each, pref
    real(rk) :: uptake, uptake_rate, food, rhs, rhs_oxy, rhs_dic, z_loss, rhs_opal !, rhs_caco3
    real(rk) :: bioom6, rhs_amm
    real(rk) :: highMort
    !real(rk) :: pcal
    ! DVM stuff
    real(rk) :: latitude, yearday, declination, day_length
    real, parameter :: pi = 3.14159265358979323846
    real(rk) :: inside_acos
    real(rk) :: scale_Rg
    real(rk) :: light_dep_mort
    real(rk) :: par
    ! -----


    _LOOP_BEGIN_

    ! Retrieve current (local) state variable values.
    _GET_(self%id_c,c)
    _GET_(self%id_par,par)
    _GET_(self%id_oxy,oxy)

    ! This is for diel vertical migration of mesozooplankton purposes
    ! Default is off
    _GET_SURFACE_(self%id_lat,latitude) ! degN
    _GET_GLOBAL_(self%id_yearday,yearday) !decimal day of the year

    latitude = latitude * pi / 180.0
    declination = 23.44 * pi / 180.0 * sin(2.0 * pi / 365.0 * (yearday - 81.0))
    inside_acos = max( -1.0_rk, min( 1.0_rk,-tan(latitude) * tan(declination) ) )

    day_length = 24.0 / pi * acos(inside_acos)
    day_length = max(0.0_rk, min(24.0_rk, day_length))
    if (self%is_migrator) then
        scale_Rg = max( self%scaleRg , min(1.0_rk,(24.0_rk - day_length)/24.0_rk) )
        ! light dependent mortality multiplier
        if (self%KsLightDep < 1.0e-18_rk) then 
            light_dep_mort = 1.0_rk ! light-dependent mortality is off by default
        else
            light_dep_mort = par / (par + self%KsLightDep) ! assumes at low light, mortality decreases
        end if
    else
        scale_Rg = 1.0_rk
        light_dep_mort = 1.0_rk ! maximum mortality everywhere for non-migrators
    end if
    ! --------------

    !_GET_(self%id_pcal,pcal)
    ! if (use_calcifier) then
    !     _GET_(self%id_caco3, caco3)
    ! end if
    ! Get prey concentrations
    do iprey=1,self%nprey
        _GET_(self%id_preyc(iprey), preyc(iprey))
        !if (use_chl) then
            _GET_(self%id_preychl(iprey), preychl(iprey))
        !end if
        bio_loss(iprey) = max(sign(-1.0_rk,preyc(iprey)-self%bio_loss_limit(iprey)),0.0_rk)
        !caco3_loss(iprey) = max(sign(-1.0_rk,caco3-0.01),0.0_rk)
    end do
    z_loss = max(sign(-1.0_rk,c - prevent_loss_Z),0.0_rk) ! self loss switch 

    pref = self%pref
!    ! Compute total available prey (mg C/m3), weighted according to effective prey preferences.
!    food_each = pref * preyc * bio_loss 
!    food = sum(food_each)
!    ! Compute rates
!    uptake_rate_each = self%grz * food_each**2 / ( (self%Rg * scale_Rg )**2 + food**2)
!    uptake_rate = sum(uptake_rate_each) ! assimilation is included below

    food_each = pref * preyc
    food = sum(food_each)
    
    uptake_rate_each = bio_loss * self%grz * pref * preyc**2/((self%Rg * scale_Rg )**2 + food**2) 
    uptake_rate = sum(uptake_rate_each) ! assimilation is included below

    ! ! Prey uptake based on a Michaelis-Menten/Type II functional response with dynamic preferences "pref".
    ! ! put_u is the relative rate of uptake (1/d), rug the absolute rate of uptake (mg C/m3/d)
    ! uptake_rate = self%grz * c / (self%Rg + food)
    ! uptake = uptake_rate * food
    ! ! Loss rates of individual prey types: sprey is the specific loss rate (1/d), fpreyc the absolute loss rate (mg C/m3/d)
    ! uptake_rate_each = uptake_rate * pref
    ! uptake_each = uptake_rate_each * food_each

    ! Apply/Calculate prey specific rates.
    rhs_opal = 0.0_rk
    !rhs_caco3 = 0.0_rk
    do iprey=1,self%nprey

!        rhs = -uptake_rate_each(iprey) * preyc(iprey) * pref(iprey) * bio_loss(iprey)
        rhs = -uptake_rate_each(iprey) * c
        _ADD_SOURCE_(self%id_preyc(iprey), rhs )

        if (self%has_chl(iprey)) then
            _ADD_SOURCE_(self%id_preychl(iprey), rhs * preychl(iprey) / preyc(iprey) )
        end if 

        rhs_opal = rhs_opal + self%opal_multiplier(iprey) * uptake_rate_each(iprey) * c
        ! if (use_calcifier) then
        !     rhs_caco3 = rhs_caco3 + pcal * (-0.5 * (self%caco3_multiplier(iprey) * caco3_loss(iprey) * uptake_rate_each(iprey) * c) )
        ! end if
    end do

    ! Below are the rates that do not require a loop like that above
    ! zoo


    highMort = self%m2 * ( c/(c + self%Km2) ) * light_dep_mort

    rhs = ( self%gamma * uptake_rate - z_loss * ( self%m * max(0.5_rk,light_dep_mort) + highMort + self%exc ) ) * c 
    _ADD_SOURCE_(self%id_c, rhs)
    ! nutrients
    rhs = z_loss * self%exc * c
    _ADD_SOURCE_(self%id_nh4, rhs)
    _ADD_SOURCE_(self%id_pho, rhs)
    ! det & dom
    rhs = ( (1.0_rk - self%gamma) * uptake_rate + z_loss * ( self%m * max(0.5_rk,light_dep_mort) + highMort ) ) * c
    _ADD_SOURCE_(self%id_det, (1.0_rk - frr) * rhs)
    _ADD_SOURCE_(self%id_dom, frr * rhs)

    _ADD_SOURCE_(self%id_opal, rhs_opal )

    bioom6 = 0.0_rk
    if (oxy > 0) then
        bioom6 = 1.0_rk
    end if 
    rhs_oxy = -bioom6 * 6.625 * ( z_loss * self%exc * c )  * Cmg_to_Cmmol * Cmmol_to_Nmmol 
    _ADD_SOURCE_(self%id_oxy, rhs_oxy)

    if (couple_co2) then
        ! if (use_calcifier) then
        !     _ADD_SOURCE_(self%id_caco3, rhs_caco3 )
        ! end if
        !rhs_dic = z_loss * self%exc * c * Cmg_to_Cmmol  - rhs_caco3
        rhs_dic = z_loss * self%exc * c * Cmg_to_Cmmol  
        
        _ADD_SOURCE_(self%id_dic, rhs_dic )
         rhs_amm = (z_loss * self%exc * c) * Cmg_to_Cmmol * Cmmol_to_Nmmol ! rhs_amm - rhs_nit where rhs_nit = 0 in this part of code
        !_ADD_SOURCE_(self%id_alk, rhs_amm -0.5_rk * rhs_oxy * (1._rk-bioom6) - rhs_caco3)
        _ADD_SOURCE_(self%id_alk, rhs_amm -0.5_rk * rhs_oxy * (1._rk-bioom6) )
    end if

    _LOOP_END_

end subroutine do



! subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)

!     class (type_ecosmo_zooplankton),intent(in) :: self
!     _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
!     real(rk) :: surface_par, par, par_bottom_EZ

!     _LOOP_BEGIN_

!     _GET_HORIZONTAL_(self%id_sfpar,surface_par)
!     _GET_(self%id_par, par)

!     par_bottom_EZ = surface_par / 1000.0_rk

!     if (self%is_migrator) then
!         if (surface_par < 1.0_rk) then
!             _SET_VERTICAL_MOVEMENT_(self%id_c,self%swimspd)
!         else 
!             if (par.gt.par_bottom_EZ) then
!                 _SET_VERTICAL_MOVEMENT_(self%id_c,-self%swimspd)
!             else
!                 _SET_VERTICAL_MOVEMENT_(self%id_c,self%swimspd)
!             end if
!         end if
!     else
!          _SET_VERTICAL_MOVEMENT_(self%id_c,-self%Zsink)
!     endif
!     _LOOP_END_
!  end subroutine get_vertical_movement
! -------------------------------------------------------------------------

end module
#include "fabm_driver.h"

! --------- Change Log ---------- !
! Veli Çağlar Yumruktepe (VCY)
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

module ecosmo_caco3
    use fabm_types
    use fabm_expressions
    use ecosmo_shared
    implicit none
    private
    public type_ecosmo_caco3
    type,extends(type_base_model), public :: type_ecosmo_caco3
        type (type_state_variable_id)     :: id_c
        type (type_state_variable_id)     :: id_alk, id_dic
        type (type_dependency_id)         :: id_Om_cal
        type (type_dependency_id)         :: id_depth
        real(rk) :: IMsink !, calcDis
        real(rk) :: cz_IMsink
        real(rk) :: dissCmax, ndissC
    contains
        procedure :: initialize
        procedure :: do
        procedure :: get_vertical_movement
    end type type_ecosmo_caco3

contains
    subroutine initialize(self,configunit)
        class (type_ecosmo_caco3), intent(inout),target  :: self
        integer,  intent(in) :: configunit
        ! get parameters
        call self%get_parameter( self%IMsink,            'IMsink',   'm/day', 'inorganic matter sinking rate', default=5.0_rk, scale_factor=1.0_rk/sedy0)
!        call self%get_parameter( self%calcDis ,          'calcDis',  '1/day', 'calcite dissolution rate',  default=0.03_rk,  scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%dissCmax,'dissCmax',    '1/day',      'maximum specific dissolution rate', default=0.03_rk, scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%ndissC,  'ndissC',      '-',          'power of the dissolution law (Keir 1980)', default=2.22_rk)
        ! register state variables
        call self%register_state_variable(self%id_c,     'c', 'mgC/m3', 'concentration in carbon units', minimum=0.0_rk, maximum=1000.0_rk)
        ! register dependencies
        call self%register_dependency(self%id_Om_cal,    'Om_cal_target', '-','calcite saturation')
        call self%register_state_dependency(self%id_dic, 'dic','mmol m-3', 'dic budget')
        call self%register_state_dependency(self%id_alk, 'alk','mmol m-3', 'alkalinity budget')

        if (depth_dependent_sinking_speed) then
            call self%get_parameter( self%cz_IMsink,'cz_IMsink',    '1/day',      'CaCO3 sinking rate increase per metre below surface', default=0.0_rk, scale_factor=1.0_rk/sedy0)        
            call self%register_dependency(self%id_depth,standard_variables%depth)
        end if
    end subroutine initialize

    subroutine do(self,_ARGUMENTS_DO_)
        class (type_ecosmo_caco3),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
        real(rk) :: c
        real(rk) :: Om_cal, Lstar, dissolution

        _LOOP_BEGIN_

        ! Retrieve local variables
        _GET_(self%id_c,c)
        _GET_(self%id_Om_cal,Om_cal)

        ! dissolution
        Lstar = self%dissCmax * max(1._rk-Om_cal,0._rk)**self%ndissC
        !Lstar = self%calcDis * max(0.0,1.0 - Om_cal)
        dissolution = Lstar * c
        _ADD_SOURCE_(self%id_c, -dissolution)
        _ADD_SOURCE_(self%id_dic, dissolution * Cmg_to_Cmmol)
        _ADD_SOURCE_(self%id_alk, dissolution * Cmg_to_Cmmol * 2.0_rk)
    
        _LOOP_END_
    end subroutine do

    subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
        class (type_ecosmo_caco3),intent(in) :: self
        _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_

        real(rk) :: w_det
        real(rk) :: depth
  
        _LOOP_BEGIN_

        if (depth_dependent_sinking_speed) then
            _GET_(self%id_depth, depth)
            depth = abs(depth) ! some host models may report depth as negative
            w_det = self%IMsink + self%cz_IMsink * depth
        else
            w_det = self%IMsink 
        end if

        _ADD_VERTICAL_VELOCITY_(self%id_c, -w_det)

        _LOOP_END_
    end subroutine get_vertical_movement

end module
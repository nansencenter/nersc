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

module ecosmo_opal
    use fabm_types
    use fabm_expressions
    use ecosmo_shared
    implicit none
    private
    public type_ecosmo_opal

    type,extends(type_base_model), public  :: type_ecosmo_opal
        type (type_state_variable_id)         :: id_c
        type (type_state_variable_id)         :: id_sil!, id_opal
        type (type_dependency_id)             :: id_depth        
        real(rk) :: regenSi, OMsink, cz_OMsink
    contains
        procedure :: initialize
        procedure :: do
        procedure :: get_vertical_movement
    end type type_ecosmo_opal

contains

    subroutine initialize(self,configunit)
        class (type_ecosmo_opal), intent(inout),target  :: self
        integer,  intent(in) :: configunit

        call self%get_parameter( self%regenSi,       'regenSi',    '1/day',      'Si regeneration rate',            default=0.015_rk, scale_factor=1.0_rk/sedy0)
        call self%get_parameter( self%OMsink,        'OMsink',     'm/day', 'organic matter sinking rate', default=5.0_rk, scale_factor=1.0_rk/sedy0)
        call self%register_state_variable(self%id_c, 'c',          'mgC/m3', 'concentration in carbon units', minimum=0.0_rk, maximum=1000.0_rk)
        call self%register_state_dependency(self%id_sil, 'sil', 'mgC/m3', 'silicate')

        ! optional routines
        if (depth_dependent_sinking_speed) then
            call self%get_parameter( self%cz_OMsink,'cz_OMsink',    '1/day',      'opal sinking rate increase per metre below surface', default=0.0_rk, scale_factor=1.0_rk/sedy0)        
            call self%register_dependency(self%id_depth,standard_variables%depth)
        end if
    end subroutine initialize

    subroutine do(self,_ARGUMENTS_DO_)

        class (type_ecosmo_opal),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
        real(rk) :: c, rhs_opal

        _LOOP_BEGIN_

        ! Retrieve current (local) state variable values.
        _GET_(self%id_c,c)
        ! Opal regeneration
        rhs_opal = self%regenSi * c
        _ADD_SOURCE_(self%id_c, -rhs_opal)
        _ADD_SOURCE_(self%id_sil, rhs_opal)
    
        _LOOP_END_
    end subroutine do

    subroutine get_vertical_movement(self,_ARGUMENTS_GET_VERTICAL_MOVEMENT_)
        class (type_ecosmo_opal),intent(in) :: self
        _DECLARE_ARGUMENTS_GET_VERTICAL_MOVEMENT_
        real(rk) :: w_det
        real(rk) :: depth
    
        _LOOP_BEGIN_

        if (depth_dependent_sinking_speed) then
            _GET_(self%id_depth, depth)
            depth = abs(depth) ! some host models may report depth as negative
            w_det = self%OMsink + self%cz_OMsink * depth
        else
            w_det = self%OMsink 
        end if

        _ADD_VERTICAL_VELOCITY_(self%id_c, -w_det)

        _LOOP_END_
    end subroutine get_vertical_movement

end module
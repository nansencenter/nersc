#include "fabm_driver.h"

module dvm_upper_lower_boundaries_percentage

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_upper_lower_boundaries_percentage

    type (type_dependency_id)                       :: id_par, id_parmean, id_migrator_food, id_depth 
    type (type_horizontal_dependency_id)            :: id_parmean0 
    type (type_horizontal_dependency_id)            :: id_migrator_food0 
!    type (type_horizontal_dependency_id)            :: id_light_present0
!    type (type_horizontal_dependency_id)            :: id_nhours
    type (type_surface_dependency_id)               :: id_par0
!    type (type_horizontal_diagnostic_variable_id)   :: id_nhours_out
    type (type_diagnostic_variable_id)              :: id_present
    type (type_bottom_dependency_id)                :: id_topo
    type (type_horizontal_dependency_id)            :: id_daylength

!type (type_dependency_id)                       :: id_temp
    contains
        procedure :: initialize
!        procedure :: do_surface
        procedure :: do

end type

contains

    subroutine initialize(self, configunit)
        class (type_upper_lower_boundaries_percentage), intent(inout), target :: self
        integer, intent(in)                                  :: configunit
        !real(rk) :: par, par0, parmean, parmean0
!call self%register_dependency(self%id_temp,standard_variables%temperature)
        call self%register_diagnostic_variable(self%id_present,'migrator_presence','-','migrators are present here')

        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_parmean0,temporal_mean(self%id_par0,period=86400._rk,resolution=3600._rk,missing_value=50.0_rk))
        call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk,missing_value=1.0_rk))
!        call self%register_dependency(self%id_light_present0,'light_presence','-','light is available at the surface')
!        call self%register_dependency(self%id_nhours,temporal_mean(self%id_light_present0,period=86400._rk,resolution=3600._rk,missing_value=12.0_rk/86400.0_rk))
        call self%register_dependency(self%id_migrator_food,'migrator_food','mgC/m3','food availability for the migrators')
        call self%register_dependency(self%id_migrator_food0,vertical_integral(self%id_migrator_food))
        call self%register_dependency(self%id_depth,standard_variables%pressure)
        call self%register_dependency(self%id_topo,standard_variables%bottom_depth )
        call self%register_dependency(self%id_daylength,'daylength','hours','number of hours light is available at the surface')

!        call self%register_diagnostic_variable(self%id_nhours_out,'nhours','-','number of daylight hours',source=source_do_surface)
    end subroutine initialize

!     subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

!         class (type_upper_lower_boundaries),intent(in) :: self
!         _DECLARE_ARGUMENTS_DO_SURFACE_

!         real(rk) :: nhours!, lpres,T,par0
!         _HORIZONTAL_LOOP_BEGIN_
!             !_GET_(self%id_temp,T)
!             !_GET_SURFACE_(self%id_par0,par0)
! !            _GET_SURFACE_(self%id_nhours,nhours)
!             !_GET_SURFACE_(self%id_light_present0,lpres)
! !            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_nhours_out, min(24.0_rk, max(0.0_rk,nhours * 86400.0_rk)))
!         _HORIZONTAL_LOOP_END_

!     end subroutine do_surface

    subroutine do(self, _ARGUMENTS_DO_)

        class (type_upper_lower_boundaries_percentage), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
    
        real(rk) :: par, par0, parmean, parmean0, nhours, food
        real(rk) :: parlog, par0log, parmeanlog, parmean0log
        real(rk) :: depth
        real(rk) :: upper_presence, lower_presence
        real(rk) :: topo

        _LOOP_BEGIN_

            _GET_SURFACE_(self%id_parmean0,parmean0)
            _GET_SURFACE_(self%id_par0,par0)
            _GET_SURFACE_(self%id_daylength,nhours)
            _GET_SURFACE_(self%id_migrator_food0,food)
            _GET_BOTTOM_(self%id_topo,topo)

            !nhours = min(24.0_rk, max(0.0_rk,nhours * 86400.0_rk))
            par0log = max(-20.0_rk, log10(par0)) - 4.0_rk
            parmean0log = max(-20.0_rk, log10(parmean0)) - 4.0_rk

            _GET_(self%id_par,par)
            _GET_(self%id_parmean,parmean)
            _GET_(self%id_depth,depth)

            parlog = max(-20.0_rk, log10(par)) - 4.0_rk
            parmeanlog = max(-20.0_rk, log10(parmean)) - 4.0_rk

            ! SPECIFY THE POSSIBLE LOCATIONS OF HIGH MIGRATOR CONCENTRATION !

            ! There are 4 cases
            ! 1. Winter Arctic night (surface parmean < 1E-10)
            ! 2. Summer Arctic night (number of daylight hours > 23.9)
            ! 3. Normal day cycle day time
            ! 4. Normal day cycle night time

            ! CASE 1
            upper_presence = 0.0_rk
            lower_presence = 0.0_rk

            if (nhours == 0.0_rk) then
                upper_presence = 0.0_rk
                lower_presence = 0.0_rk
                
                ! Calculate possibilities above the lower boundary
                ! if (food <= 17.12_rk) then
                !     if (depth < 184.04_rk) then
                !         upper_presence = 1.0_rk
                !     else
                !         upper_presence = 0.0_rk
                !     end if
                ! else
                    if (depth < 240.0_rk) then
                        upper_presence = 1.0_rk
                    else
                        upper_presence = 0.0_rk
                    end if
!                end if
                
                ! Set diagnostic based on presence
                if (upper_presence + lower_presence > 0.9_rk) then
                    _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                else
                    _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                end if
                
            else

                ! CASE 2
                if (nhours == 24.0_rk) then
                    ! there is an upper and a lower light boundary
                    ! first calculate possibilities above the lower boundary

                    ! Initialize presence variables
                    upper_presence = 0.0_rk
                    lower_presence = 0.0_rk
                    
                    ! ! Lowerlight Rules
                    ! if (parmean0log <= -3.49_rk) then
                    !     if (par0log <= -3.09_rk) then
                    !         if (par0log/parmeanlog > 0.18_rk) then
                    !             upper_presence = 1.0_rk
                    !         else
                    !             upper_presence = 0.0_rk
                    !         end if
                    !     else
                    !         if (par0log/parmeanlog > 0.22_rk) then
                    !             upper_presence = 1.0_rk
                    !         else
                    !             upper_presence = 0.0_rk
                    !         end if
                    !     end if
                    ! else
                    !     if (parmean0log <= -2.81_rk) then
                    !         if (par0log/parmeanlog > 0.16_rk) then
                    !             upper_presence = 1.0_rk
                    !         else
                    !             upper_presence = 0.0_rk
                    !         end if
                    !     else
                    !         if (par0log/parmeanlog > 0.20_rk) then
                    !             upper_presence = 1.0_rk
                    !         else
                    !             upper_presence = 0.0_rk
                    !         end if
                    !     end if
                    ! end if

                    ! Lowerlight Rules
                    if (parmean0log <= -3.49_rk) then
                        if (par0log <= -3.09_rk) then
                            if (food <= 23.26_rk) then
                                if (par0log <= -3.82_rk) then
                                    if (par0log/parmeanlog > 0.21_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parmeanlog > 0.20_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 24.12_rk) then
                                    if (par0log/parmeanlog > 0.17_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parmeanlog > 0.19_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= -3.00_rk) then
                                if (par0log <= -3.02_rk) then
                                    if (par0log/parmeanlog > 0.22_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parmeanlog > 0.26_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log/parmeanlog > 0.16_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            end if
                        end if
                    else
                        if (parmean0log <= -2.81_rk) then
                            if (parmean0log <= -2.93_rk) then
                                if (parmean0log <= -2.94_rk) then
                                    if (par0log/parmeanlog > 0.16_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parmeanlog > 0.21_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= -2.38_rk) then
                                    if (par0log/parmeanlog > 0.13_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parmeanlog > 0.18_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= -2.27_rk) then
                                if (par0log <= -2.42_rk) then
                                    if (par0log/parmeanlog > 0.21_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parmeanlog > 0.22_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= -2.25_rk) then
                                    if (par0log/parmeanlog > 0.17_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parmeanlog > 0.19_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    end if

                    ! ! Upperlight Rules
                    ! if (parmean0log <= -2.81_rk) then
                    !     if (parmean0log <= -3.35_rk) then
                    !         if (par0log/parlog < 0.29_rk) then
                    !             lower_presence = 1.0_rk
                    !         else
                    !             lower_presence = 0.0_rk
                    !         end if
                    !     else
                    !         if (par0log/parlog < 0.33_rk) then
                    !             lower_presence = 1.0_rk
                    !         else
                    !             lower_presence = 0.0_rk
                    !         end if
                    !     end if
                    ! else
                    !     if (par0log <= -2.20_rk) then
                    !         if (par0log/parlog < 0.44_rk) then
                    !             lower_presence = 1.0_rk
                    !         else
                    !             lower_presence = 0.0_rk
                    !         end if
                    !     else
                    !         if (par0log/parlog < 0.41_rk) then
                    !                 lower_presence = 1.0_rk
                    !         else
                    !                 lower_presence = 0.0_rk
                    !         end if
                    !     end if
                    ! end if

                    ! Upperlight Rules
                    if (parmean0log <= -2.81_rk) then
                        if (parmean0log <= -3.35_rk) then
                            if (parmean0log <= -4.22_rk) then
                                if (food <= 24.09_rk) then
                                    if (par0log/parlog < 0.31_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parlog < 0.35_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (parmean0log <= -3.45_rk) then
                                    if (par0log/parlog < 0.28_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parlog < 0.24_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= -3.19_rk) then
                                if (par0log <= -3.05_rk) then
                                    if (par0log/parlog < 0.46_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parlog < 0.42_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= -2.61_rk) then
                                    if (par0log/parlog < 0.30_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parlog < 0.33_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    else
                        if (par0log <= -2.20_rk) then
                            if (par0log <= -2.48_rk) then
                                if (par0log <= -2.50_rk) then
                                    if (par0log/parlog < 0.38_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parlog < 0.42_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= -2.21_rk) then
                                    if (par0log/parlog < 0.44_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parlog < 0.50_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= -2.17_rk) then
                                if (par0log <= -2.19_rk) then
                                    if (par0log/parlog < 0.37_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parlog < 0.37_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= -2.16_rk) then
                                    if (par0log/parlog < 0.49_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log/parlog < 0.42_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    end if
                                        
                    ! Set diagnostic based on presence
                    if (upper_presence + lower_presence > 1.0_rk) then
                        _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                    else
                        if (upper_presence > 0.9_rk .and. depth >= max(topo - 20.0_rk, 0.0_rk) ) then 
                            _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                        else 
                            _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                        end if
                    end if
                else

                    ! CASE 3
                    if (par0log > -4.0_rk) then
                        ! there is an upper and a lower light boundary
                        ! first calculate possibilities above the lower boundary
                        
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Lowerlight Rules
                        if (par0log/parmeanlog > 0.17_rk) then
                            upper_presence = 1.0_rk
                        else
                            upper_presence = 0.0_rk
                        end if
                        
                        ! Upperlight Rules
                        if (par0log/parlog < 0.32_rk) then
                            lower_presence = 1.0_rk
                        else
                            lower_presence = 0.0_rk
                        end if
                        
                        
                        ! Set diagnostic based on presence
                        if (upper_presence + lower_presence > 1.0_rk) then
                            _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                        else
                            if (upper_presence > 0.9_rk .and. depth >= max(topo - 20.0_rk, 0.0_rk) ) then 
                                _SET_DIAGNOSTIC_(self%id_present,1.0_rk)
                            else 
                                _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                            end if
                        end if
                        
                    else
                        ! CASE 4
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Calculate possibilities above the lower boundary
                        if (par0log/parmeanlog > 0.17_rk) then
                            upper_presence = 1.0_rk
                        else
                            upper_presence = 0.0_rk
                        end if
                        
                        ! Set diagnostic based on presence
                        if (upper_presence + lower_presence > 0.9_rk) then
                            _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                        else
                            _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                        end if
                        
                    end if
                end if

            end if

            ! This should ensure that each point at least receives the 0.0_rk value
            if (upper_presence + lower_presence < 0.9_rk) then
                _SET_DIAGNOSTIC_(self%id_present,0.0_rk)     
            end if 
            ! 


            ! ------------------------------------------------------------- !


        _LOOP_END_

    end subroutine do

end module

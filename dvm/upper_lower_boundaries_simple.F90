#include "fabm_driver.h"

module dvm_upper_lower_boundaries_simple

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_upper_lower_boundaries_simple

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
        class (type_upper_lower_boundaries_simple), intent(inout), target :: self
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

        class (type_upper_lower_boundaries_simple), intent(in) :: self
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
            par0log = max(-20.0_rk, log10(par0))
            parmean0log = max(-20.0_rk, log10(parmean0))

            _GET_(self%id_par,par)
            _GET_(self%id_parmean,parmean)
            _GET_(self%id_depth,depth)

            parlog = max(-20.0_rk, log10(par))
            parmeanlog = max(-20.0_rk, log10(parmean))

            ! SPECIFY THE POSSIBLE LOCATIONS OF HIGH MIGRATOR CONCENTRATION !

            ! There are 4 cases
            ! 1. Winter Arctic night (surface parmean < 1E-10)
            ! 2. Summer Arctic night (number of daylight hours > 23.9)
            ! 3. Normal day cycle day time
            ! 4. Normal day cycle night time

            ! CASE 1
            upper_presence = 0.0_rk
            lower_presence = 0.0_rk

            if (nhours < 0.05_rk) then
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
                if (nhours > 23.95_rk) then
                    ! there is an upper and a lower light boundary
                    ! first calculate possibilities above the lower boundary

                    ! Initialize presence variables
                    upper_presence = 0.0_rk
                    lower_presence = 0.0_rk
                    
                    ! Lowerlight Rules
                    ! if (food <= 44.51_rk) then
                    !     if (food <= 24.07_rk) then
                            if (parmeanlog > -15.0_rk) then
                                upper_presence = 1.0_rk
                            else
                                upper_presence = 0.0_rk
                            end if
                    !     else
                    !         if (parmeanlog > -17.77_rk) then
                    !             upper_presence = 1.0_rk
                    !         else
                    !             upper_presence = 0.0_rk
                    !         end if
                    !     end if
                    ! else
                    !     if (parmean0log <= 1.19_rk) then
                    !         if (parmeanlog > -13.07_rk) then
                    !             upper_presence = 1.0_rk
                    !         else
                    !             upper_presence = 0.0_rk
                    !         end if
                    !     else
                    !         if (parmeanlog > -7.52_rk) then
                    !             upper_presence = 1.0_rk
                    !         else
                    !             upper_presence = 0.0_rk
                    !         end if
                    !     end if
                    ! end if
                    
                    ! ! Upperlight Rules
                    ! if (parmean0log <= 0.65_rk) then
                    !     !if (parmean0log <= 0.55_rk) then
                    !         if (parlog < -8.0_rk) then
                    !             lower_presence = 1.0_rk
                    !         else
                    !             lower_presence = 0.0_rk
                    !         end if
                    !     !else
                    !     !    if (parlog < -9.08_rk) then
                    !     !        lower_presence = 1.0_rk
                    !     !    else
                    !     !        lower_presence = 0.0_rk
                    !     !    end if
                    !     !end if
                    ! else
                    !     !if (parmean0log <= 1.19_rk) then
                    !         if (parlog < -4.0_rk) then
                    !             lower_presence = 1.0_rk
                    !         else
                    !             lower_presence = 0.0_rk
                    !         end if
                    !     !else
                    !     !    if (parlog < -1.41_rk) then
                    !     !        lower_presence = 1.0_rk
                    !     !    else
                    !     !        lower_presence = 0.0_rk
                    !     !    end if
                    !     !end if
                    ! end if
! ! lower light rules
!                             if (parmean0log <= 0.67_rk) then
!                                 if (parmean0log <= 0.48_rk) then
!                                     if (food <= 22.83_rk) then
!                                         if (food <= 22.04_rk) then
!                                             if (parmeanlog > -13.41_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -14.32_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (parmean0log <= -0.25_rk) then
!                                             if (parmeanlog > -17.29_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -15.23_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 else
!                                     if (parmean0log <= 0.60_rk) then
!                                         if (parmean0log <= 0.57_rk) then
!                                             if (parmeanlog > -17.36_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -15.35_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (food <= 22.06_rk) then
!                                             if (parmeanlog > -18.53_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -17.92_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 end if
!                             else
!                                 if (parmean0log <= 1.11_rk) then
!                                     if (parmean0log <= 0.90_rk) then
!                                         if (parmean0log <= 0.86_rk) then
!                                             if (parmeanlog > -11.75_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -5.86_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (parmean0log <= 1.00_rk) then
!                                             if (parmeanlog > -14.07_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -11.08_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 else
!                                     if (par0log <= 1.32_rk) then
!                                         if (par0log <= 1.18_rk) then
!                                             if (parmeanlog > -8.54_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -14.97_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (parmean0log <= 1.24_rk) then
!                                             if (parmeanlog > -7.83_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -10.14_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 end if
!                             end if


                            ! updated lower light rules
                            if (parmean0log <= 0.67_rk) then
                                if (parmeanlog > -15.84_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            else
                                if (parmeanlog > -11.01_rk) then
                                    upper_presence = 1.0_rk
                                else
                                    upper_presence = 0.0_rk
                                end if
                            end if


                            ! lower light rules
                            ! if (parmean0log <= 0.67_rk) then
                            !     if (parmean0log <= 0.48_rk) then
                            !         if (food <= 22.83_rk) then
                            !             if (parmeanlog > -13.99_rk) then
                            !                 upper_presence = 1.0_rk
                            !             else
                            !                 upper_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parmeanlog > -15.92_rk) then
                            !                 upper_presence = 1.0_rk
                            !             else
                            !                 upper_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     else
                            !         if (parmean0log <= 0.60_rk) then
                            !             if (parmeanlog > -16.69_rk) then
                            !                 upper_presence = 1.0_rk
                            !             else
                            !                 upper_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parmeanlog > -18.33_rk) then
                            !                 upper_presence = 1.0_rk
                            !             else
                            !                 upper_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     end if
                            ! else
                            !     if (parmean0log <= 1.11_rk) then
                            !         if (parmean0log <= 0.90_rk) then
                            !             if (parmeanlog > -10.57_rk) then
                            !                 upper_presence = 1.0_rk
                            !             else
                            !                 upper_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parmeanlog > -12.98_rk) then
                            !                 upper_presence = 1.0_rk
                            !             else
                            !                 upper_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     else
                            !         if (par0log <= 1.32_rk) then
                            !             if (parmeanlog > -14.16_rk) then
                            !                 upper_presence = 1.0_rk
                            !             else
                            !                 upper_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parmeanlog > -8.90_rk) then
                            !                 upper_presence = 1.0_rk
                            !             else
                            !                 upper_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     end if
                            ! end if


                            ! !Upperlight Rules
                            ! if (food <= 49.69_rk) then
                            !     if (parlog < -8.15_rk) then
                            !         lower_presence = 1.0_rk
                            !     else
                            !         lower_presence = 0.0_rk
                            !     end if
                            ! else
                            !     if (parlog < -3.83_rk) then
                            !         lower_presence = 1.0_rk
                            !     else
                            !         lower_presence = 0.0_rk
                            !     end if
                            ! end if

                            ! Upperlight Rules
                            if (parmean0log <= 0.65_rk) then
                                !if (parmean0log <= 0.55_rk) then
                                    if (parlog < -8.0_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                !else
                                !    if (parlog < -9.08_rk) then
                                !        lower_presence = 1.0_rk
                                !    else
                                !        lower_presence = 0.0_rk
                                !    end if
                                !end if
                            else
                                !if (parmean0log <= 1.19_rk) then
                                    if (parlog < -4.0_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                !else
                                !    if (parlog < -1.41_rk) then
                                !        lower_presence = 1.0_rk
                                !    else
                                !        lower_presence = 0.0_rk
                                !    end if
                                !end if
                            end if

                            ! !Upperlight Rules
                            ! if (food <= 49.69_rk) then
                            !     if (parmean0log <= -0.25_rk) then
                            !         if (parlog < -7.29_rk) then
                            !             lower_presence = 1.0_rk
                            !         else
                            !             lower_presence = 0.0_rk
                            !         end if
                            !     else
                            !         if (parlog < -8.41_rk) then
                            !             lower_presence = 1.0_rk
                            !         else
                            !             lower_presence = 0.0_rk
                            !         end if
                            !     end if
                            ! else
                            !     if (parmean0log <= 1.11_rk) then
                            !         if (parlog < -4.47_rk) then
                            !             lower_presence = 1.0_rk
                            !         else
                            !             lower_presence = 0.0_rk
                            !         end if
                            !     else
                            !         if (par0log <= 1.32_rk) then
                            !             if (parlog < -8.16_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parlog < -2.38_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     end if
                            ! end if

                            ! !Upperlight Rules
                            ! if (food <= 49.69_rk) then
                            !     if (parmean0log <= -0.25_rk) then
                            !         if (parmean0log <= -0.34_rk) then
                            !             if (parlog < -7.89_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parlog < -6.09_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     else
                            !         if (food <= 18.18_rk) then
                            !             if (parlog < -0.86_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parlog < -8.47_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     end if
                            ! else
                            !     if (parmean0log <= 1.11_rk) then
                            !         if (parmean0log <= 0.90_rk) then
                            !             if (parlog < -3.91_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parlog < -4.96_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     else
                            !         if (par0log <= 1.32_rk) then
                            !             if (parlog < -8.16_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         else
                            !             if (parlog < -2.15_rk) then
                            !                 lower_presence = 1.0_rk
                            !             else
                            !                 lower_presence = 0.0_rk
                            !             end if
                            !         end if
                            !     end if
                            ! end if

! !Upperlight Rules
!                             if (food <= 49.69_rk) then
!                                 if (parmean0log <= -0.25_rk) then
!                                     if (parmean0log <= -0.34_rk) then
!                                         if (food <= 24.52_rk) then
!                                             if (parlog < -8.30_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parlog < -6.64_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (food <= 23.52_rk) then
!                                             if (parlog < -5.45_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parlog < -6.15_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 else
!                                     if (food <= 18.18_rk) then
!                                         if (par0log <= 0.80_rk) then
!                                             if (parlog > 0.78_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parlog < -2.50_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (food <= 25.48_rk) then
!                                             if (parlog < -8.80_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parlog < -7.57_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 end if
!                             else
!                                 if (parmean0log <= 1.11_rk) then
!                                     if (parmean0log <= 0.90_rk) then
!                                         if (parmean0log <= 0.86_rk) then
!                                             if (parlog < -4.53_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parlog < -1.49_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (par0log <= 1.21_rk) then
!                                             if (parlog < -6.08_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parlog < -4.64_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 else
!                                     if (par0log <= 1.32_rk) then
!                                         if (par0log <= 1.26_rk) then
!                                             if (parlog < -7.94_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parlog < -9.67_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (parmean0log <= 1.24_rk) then
!                                             if (parlog < -1.36_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parlog < -3.08_rk) then
!                                                 lower_presence = 1.0_rk
!                                             else
!                                                 lower_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 end if
!                             end if                            
! ! upper light rules
!                             if (food <= 49.69_rk) then
!                                 if (parmean0log <= -0.25_rk) then
!                                     if (parmean0log <= -0.34_rk) then
!                                         if (food <= 24.52_rk) then
!                                             if (parmeanlog > -8.30_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -6.64_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (food <= 23.52_rk) then
!                                             if (parmeanlog > -5.45_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -6.15_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 else
!                                     if (food <= 18.18_rk) then
!                                         if (food <= 18.17_rk) then
!                                             if (parmeanlog > 0.78_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -2.50_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (food <= 25.48_rk) then
!                                             if (parmeanlog > -8.80_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -7.57_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 end if
!                             else
!                                 if (parmean0log <= 1.11_rk) then
!                                     if (parmean0log <= 0.90_rk) then
!                                         if (parmean0log <= 0.86_rk) then
!                                             if (parmeanlog > -4.53_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -1.49_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (par0log <= 1.21_rk) then
!                                             if (parmeanlog > -6.08_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -4.64_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 else
!                                     if (par0log <= 1.32_rk) then
!                                         if (par0log <= 1.26_rk) then
!                                             if (parmeanlog > -7.94_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -9.67_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     else
!                                         if (parmean0log <= 1.24_rk) then
!                                             if (parmeanlog > -1.36_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         else
!                                             if (parmeanlog > -3.08_rk) then
!                                                 upper_presence = 1.0_rk
!                                             else
!                                                 upper_presence = 0.0_rk
!                                             end if
!                                         end if
!                                     end if
!                                 end if
!                             end if                            

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
                    if (par0log > 0.0_rk) then
                        ! there is an upper and a lower light boundary
                        ! first calculate possibilities above the lower boundary
                        
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Lowerlight Rules
                        if (parmeanlog > -15.0_rk) then
                            upper_presence = 1.0_rk
                        else
                            upper_presence = 0.0_rk
                        end if
                        
                        ! Upperlight Rules
                        if (parlog < -6.5_rk) then
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
                        if (parmeanlog > -15.0_rk) then
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

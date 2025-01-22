#include "fabm_driver.h"

module dvm_upper_lower_boundaries

use fabm_types
use fabm_expressions

implicit none

private 

type, extends(type_base_model), public :: type_upper_lower_boundaries

    type (type_dependency_id)                       :: id_par, id_parmean, id_migrator_food, id_depth 
    type (type_horizontal_dependency_id)            :: id_parmean0 
    type (type_horizontal_dependency_id)            :: id_migrator_food0 
    type (type_horizontal_dependency_id)            :: id_light_present0
    type (type_horizontal_dependency_id)            :: id_nhours
    type (type_surface_dependency_id)               :: id_par0
    type (type_horizontal_diagnostic_variable_id)   :: id_nhours_out
    type (type_diagnostic_variable_id)              :: id_present
    type (type_bottom_dependency_id)                :: id_topo

!type (type_dependency_id)                       :: id_temp
    contains
        procedure :: initialize
        procedure :: do_surface
        procedure :: do

end type

contains

    subroutine initialize(self, configunit)
        class (type_upper_lower_boundaries), intent(inout), target :: self
        integer, intent(in)                                  :: configunit
        !real(rk) :: par, par0, parmean, parmean0
!call self%register_dependency(self%id_temp,standard_variables%temperature)
        call self%register_diagnostic_variable(self%id_present,'migrator_presence','-','migrators are present here')

        call self%register_dependency(self%id_par, standard_variables%downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_par0, standard_variables%surface_downwelling_photosynthetic_radiative_flux)
        call self%register_dependency(self%id_parmean0,temporal_mean(self%id_par0,period=86400._rk,resolution=3600._rk,missing_value=50.0_rk))
        call self%register_dependency(self%id_parmean,temporal_mean(self%id_par,period=86400._rk,resolution=3600._rk,missing_value=1.0_rk))
        call self%register_dependency(self%id_light_present0,'light_presence','-','light is available at the surface')
        call self%register_dependency(self%id_nhours,temporal_mean(self%id_light_present0,period=86400._rk,resolution=3600._rk,missing_value=12.0_rk/86400.0_rk))
        call self%register_dependency(self%id_migrator_food,'migrator_food','mgC/m3','food availability for the migrators')
        call self%register_dependency(self%id_migrator_food0,vertical_integral(self%id_migrator_food))
        call self%register_dependency(self%id_depth,standard_variables%pressure)
        call self%register_dependency(self%id_topo,standard_variables%bottom_depth )

        call self%register_diagnostic_variable(self%id_nhours_out,'nhours','-','number of daylight hours',source=source_do_surface)
    end subroutine initialize

    subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)

        class (type_upper_lower_boundaries),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_SURFACE_

        real(rk) :: nhours!, lpres,T,par0
        _HORIZONTAL_LOOP_BEGIN_
            !_GET_(self%id_temp,T)
            !_GET_SURFACE_(self%id_par0,par0)
            _GET_SURFACE_(self%id_nhours,nhours)
            !_GET_SURFACE_(self%id_light_present0,lpres)
            _SET_HORIZONTAL_DIAGNOSTIC_(self%id_nhours_out, min(24.0_rk, max(0.0_rk,nhours * 86400.0_rk)))
        _HORIZONTAL_LOOP_END_

    end subroutine do_surface

    subroutine do(self, _ARGUMENTS_DO_)

        class (type_upper_lower_boundaries), intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
    
        real(rk) :: par, par0, parmean, parmean0, nhours, food
        real(rk) :: parlog, par0log, parmeanlog, parmean0log
        real(rk) :: depth
        real(rk) :: upper_presence, lower_presence
        real(rk) :: topo

        _LOOP_BEGIN_

            _GET_SURFACE_(self%id_parmean0,parmean0)
            _GET_SURFACE_(self%id_par0,par0)
            _GET_SURFACE_(self%id_nhours,nhours)
            _GET_SURFACE_(self%id_migrator_food0,food)
            _GET_BOTTOM_(self%id_topo,topo)

            nhours = min(24.0_rk, max(0.0_rk,nhours * 86400.0_rk))
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

            if (parmean0 < 1E-3_rk) then
                upper_presence = 0.0_rk
                lower_presence = 0.0_rk
                
                ! Calculate possibilities above the lower boundary
                if (par0log <= -14.04_rk) then
                    if (depth < 214.29_rk) then
                        upper_presence = 1.0_rk
                    else
                        upper_presence = 0.0_rk
                    end if
                else
                    if (depth < 276.80_rk) then
                        upper_presence = 1.0_rk
                    else
                        upper_presence = 0.0_rk
                    end if
                end if
                
                ! Set diagnostic based on presence
                if (upper_presence + lower_presence > 0.9_rk) then
                    _SET_DIAGNOSTIC_(self%id_present, 1.0_rk)
                else
                    _SET_DIAGNOSTIC_(self%id_present, 0.0_rk)
                end if
                
            else

                ! CASE 2
                if (nhours > 23.9_rk) then
                    ! there is an upper and a lower light boundary
                    ! first calculate possibilities above the lower boundary

                    ! Initialize presence variables
                    upper_presence = 0.0_rk
                    lower_presence = 0.0_rk
                    
                    ! Lowerlight Rules
                    if (parmean0log <= 1.22_rk) then
                        if (food <= 790.53_rk) then
                            if (food <= 566.93_rk) then
                                if (par0log <= 1.04_rk) then
                                    if (parmeanlog > -15.17_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -16.39_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= 1.33_rk) then
                                    if (parmeanlog > -18.88_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -18.14_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= 0.92_rk) then
                                if (par0log <= 0.43_rk) then
                                    if (parmeanlog > -9.18_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -13.62_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (parmean0log <= 1.10_rk) then
                                    if (parmeanlog > -17.47_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -14.90_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    else
                        if (parmean0log <= 1.61_rk) then
                            if (food <= 597.60_rk) then
                                if (food <= 568.51_rk) then
                                    if (parmeanlog > -13.88_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -15.96_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (parmean0log <= 1.24_rk) then
                                    if (parmeanlog > -11.85_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -13.78_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= 1.55_rk) then
                                if (par0log <= 1.31_rk) then
                                    if (parmeanlog > -9.02_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -8.16_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= 1.73_rk) then
                                    if (parmeanlog > -9.79_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (parmeanlog > -8.88_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    end if
                    
                    ! Upperlight Rules
                    if (food <= 597.60_rk) then
                        if (par0log <= 0.98_rk) then
                            if (par0log <= 0.45_rk) then
                                if (food <= 566.57_rk) then
                                    if (parlog < -5.26_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -5.63_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 567.18_rk) then
                                    if (parlog < -8.08_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -6.80_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= 1.15_rk) then
                                if (par0log <= 1.04_rk) then
                                    if (parlog < -8.43_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -9.78_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 571.23_rk) then
                                    if (parlog < -8.96_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -8.46_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        end if
                    else
                        if (parmean0log <= 0.92_rk) then
                            if (par0log <= 0.91_rk) then
                                if (parmean0log <= 0.81_rk) then
                                    if (parlog < -9.61_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -7.40_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 780.74_rk) then
                                    if (parlog < -10.54_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -9.88_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (parmean0log <= 1.42_rk) then
                                if (par0log <= 1.04_rk) then
                                    if (parlog < -4.90_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -6.21_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 715.80_rk) then
                                    if (parlog < -3.31_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -2.61_rk) then
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
                    if (par0 > 1E-3_rk) then
                        ! there is an upper and a lower light boundary
                        ! first calculate possibilities above the lower boundary
                        
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Lowerlight Rules
                        if (parmean0log <= 0.60_rk) then
                            if (food <= 660.50_rk) then
                                if (nhours <= 2.81_rk) then
                                    if (parmeanlog > -16.95_rk) then
                                        upper_presence = 1.0_rk
                                    else
                                        upper_presence = 0.0_rk
                                    end if
                                else
                                    if (par0log <= 0.07_rk) then
                                        if (parmeanlog > -14.72_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -15.85_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            else
                                if (par0log <= -1.07_rk) then
                                    if (par0log <= -1.29_rk) then
                                        if (parmeanlog > -18.33_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -16.51_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (par0log <= 0.33_rk) then
                                        if (parmeanlog > -19.77_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -20.36_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= 0.31_rk) then
                                if (parmean0log <= 1.07_rk) then
                                    if (par0log <= -0.95_rk) then
                                        if (parmeanlog > -11.11_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -14.74_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (par0log <= -0.48_rk) then
                                        if (parmeanlog > -9.68_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -10.20_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            else
                                if (parmean0log <= 0.65_rk) then
                                    if (par0log <= 0.92_rk) then
                                        if (parmeanlog > -16.53_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -17.76_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (food <= 563.88_rk) then
                                        if (parmeanlog > -16.00_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -14.53_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            end if
                        end if
                        
                        ! Upperlight Rules
                        if (parmean0log <= 1.07_rk) then
                            if (nhours <= 20.66_rk) then
                                if (parmean0log <= 0.39_rk) then
                                    if (parlog < -9.01_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -7.79_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (food <= 764.27_rk) then
                                    if (parlog < -10.60_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -8.03_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            end if
                        else
                            if (par0log <= 1.33_rk) then
                                if (nhours <= 22.43_rk) then
                                    if (parlog < -5.94_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -4.64_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                end if
                            else
                                if (par0log <= 1.49_rk) then
                                    if (parlog < -6.50_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
                                    end if
                                else
                                    if (parlog < -6.86_rk) then
                                        lower_presence = 1.0_rk
                                    else
                                        lower_presence = 0.0_rk
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
                        ! CASE 4
                        ! Initialize presence variables
                        upper_presence = 0.0_rk
                        lower_presence = 0.0_rk
                        
                        ! Calculate possibilities above the lower boundary
                        if (parmean0log <= -0.03_rk) then
                            if (food <= 660.04_rk) then
                                if (nhours <= 2.80_rk) then
                                    if (food <= 561.85_rk) then
                                        if (parmeanlog > -12.79_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -16.50_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (food <= 561.88_rk) then
                                        if (parmeanlog > -11.04_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -10.82_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            else
                                if (nhours <= 7.15_rk) then
                                    if (food <= 758.27_rk) then
                                        if (parmeanlog > -17.66_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -16.52_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (food <= 758.93_rk) then
                                        if (parmeanlog > -12.92_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -15.81_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            end if
                        else
                            if (food <= 562.58_rk) then
                                if (food <= 562.40_rk) then
                                    if (food <= 562.24_rk) then
                                        if (parmeanlog > -8.22_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -13.08_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (food <= 562.44_rk) then
                                        if (parmeanlog > -8.17_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -9.29_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            else
                                if (food <= 563.75_rk) then
                                    if (food <= 562.71_rk) then
                                        if (parmeanlog > -15.42_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -12.44_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                else
                                    if (nhours <= 17.08_rk) then
                                        if (parmeanlog > -9.68_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    else
                                        if (parmeanlog > -12.40_rk) then
                                            upper_presence = 1.0_rk
                                        else
                                            upper_presence = 0.0_rk
                                        end if
                                    end if
                                end if
                            end if
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

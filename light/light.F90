#include "fabm_driver.h"

! This is a FABM implementation of the two-band light model used in
! the General Ocean Turbulence Model (GOTM)
!
! It distinguishes two wavebands:
! * visible (equivalent to photosynthetically activate radation, PAR, 400 - 700 nm)
! * non-visible, which combines ultraviolet (< 400 nm) and infrared (> 700 nm)
! The non-visible fraction is generally absorbed close to the surface.
! The visible fraction typically penetrates deeper into the water.
! Attenuation of the visible band is influenced by FABM variables that contribute to
! standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux
! (as well as by the background attenuation set by parameter g2)
!
! The model is driven by downwelling shortwave radiation just below the water surface,
! which in FABM is denoted with standard_variables%surface_downwelling_shortwave_flux
! This is the radiation left after reflection by the surface [albedo] is accounted for.

! Added modification of light from surface variables, e.g. ice-algae, and
! moved from gotm to nersc, Jan 2026, Shuang Gao 

module nersc_light

   use fabm_types

   implicit none

   private

   type, extends(type_base_model), public :: type_nersc_light
      ! Identifiers for dependencies [model inputs]
      type (type_surface_dependency_id) :: id_swr0 ! Surface shortwave radiation
      type (type_dependency_id)         :: id_dz   ! Cell thickness
      type (type_dependency_id)         :: id_ext  ! Attentuation coefficient for PAR
      type (type_surface_dependency_id) :: id_ext_sf  ! Surface attentuation coefficient for PAR

      ! Identifiers for diagnostic variables [model outputs]
      type (type_diagnostic_variable_id)         :: id_par  ! Photosynthetically active radiation
      type (type_diagnostic_variable_id)         :: id_swr  ! Shortwave radiation
      type (type_surface_diagnostic_variable_id) :: id_par0 ! Surface photosynthetically active radiation

      ! Parameters
      real(rk) :: a, g1, g2
      logical  :: couple_ice
   contains
      ! Model procedures
      procedure :: initialize
      procedure :: do_column
   end type type_nersc_light

contains

   subroutine initialize(self, configunit)
      class (type_nersc_light), intent(inout), target :: self
      integer,                 intent(in)            :: configunit

      call self%get_parameter(self%a,  'a',  '-','non-visible fraction of shortwave radiation', default=0.58_rk) 
      call self%get_parameter(self%g1, 'g1', 'm','e-folding depth of non-visible fraction',     default=0.35_rk)
      call self%get_parameter(self%g2, 'g2', 'm','e-folding depth of visible fraction',         default=23.0_rk) 
      call self%get_parameter(self%couple_ice, 'couple_ice', 'true/false', 'switch coupling to ice', default=.false.)

      ! Register diagnostic variables
      call self%register_diagnostic_variable(self%id_swr, 'swr', 'W m-2', 'shortwave radiation', &
         standard_variable=standard_variables%downwelling_shortwave_flux, source=source_do_column)
      call self%register_diagnostic_variable(self%id_par, 'par', 'W m-2', 'photosynthetically active radiation', &
         standard_variable=standard_variables%downwelling_photosynthetic_radiative_flux, source=source_do_column)
      call self%register_diagnostic_variable(self%id_par0, 'par0', 'W m-2', 'surface photosynthetically active radiation', &
         standard_variable=standard_variables%surface_downwelling_photosynthetic_radiative_flux, source=source_do_column)

      ! Register environmental dependencies (temperature, shortwave radiation)
      call self%register_dependency(self%id_swr0, standard_variables%surface_downwelling_shortwave_flux)
      call self%register_dependency(self%id_ext,  standard_variables%attenuation_coefficient_of_photosynthetic_radiative_flux)
      if(self%couple_ice) then
        call self%register_dependency(self%id_ext_sf, 'surface_attenuation_coef_IA', '-' ,'surface attenuation by ice algae')
      end if
      call self%register_dependency(self%id_dz,   standard_variables%cell_thickness)
   end subroutine
   
   subroutine do_column(self, _ARGUMENTS_DO_COLUMN_)
      class (type_nersc_light), intent(in) :: self
      _DECLARE_ARGUMENTS_DO_COLUMN_

      real(rk) :: swr0, dz, swr, par, z, ext, bioext, ext_sf

      _GET_SURFACE_(self%id_swr0,swr0)
      if(self%couple_ice) then
      ! Calculate attenuation of surface item (e.g. ice algae)
        _GET_SURFACE_(self%id_ext_sf,ext_sf)   ! attenuation of surface item (unitless)
        swr0 = swr0 * (1.0_rk - self%a) * exp(-ext_sf) + swr0 * self%a
        _SET_SURFACE_DIAGNOSTIC_(self%id_par0,swr0*(1.0_rk - self%a)*exp(-ext_sf)+swr0*self%a)
      else
        _SET_SURFACE_DIAGNOSTIC_(self%id_par0,swr0 * (1.0_rk - self%a))
      end if
      z = 0.0_rk
      bioext = 0.0_rk
      _DOWNWARD_LOOP_BEGIN_
         _GET_(self%id_dz,dz)     ! Layer height (m)
         _GET_(self%id_ext,ext)   ! PAR attenuation (m-1)

         ! Set depth to centre of layer
         z = z + dz * 0.5_rk
         bioext = bioext + ext * dz * 0.5_rk

         ! Calculate photosynthetically active radiation (PAR), shortwave radiation, and PAR attenuation.
         par = swr0 * (1.0_rk - self%a) * exp(-z / self%g2 - bioext)
         swr = par + swr0 * self%a * exp(-z / self%g1)

         ! Move to bottom of layer
         z = z + dz * 0.5_rk
         bioext = bioext + ext * dz * 0.5_rk

         _SET_DIAGNOSTIC_(self%id_swr,swr) ! Shortwave radiation at layer centre
         _SET_DIAGNOSTIC_(self%id_par,par) ! Photosynthetically active radiation at layer centre
      _DOWNWARD_LOOP_END_
   end subroutine do_column

end module nersc_light

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
! VCY - 08/09/2026
! Added optional parameter for zooplankton prey switching when set to True in fabm.yaml 
! (default: false). This allows the model to use the original prey preference-based grazing 
! ------------------------------- !

module ecosmo_shared
   use fabm_types
   implicit none
   public

   real(rk), parameter :: sedy0 = 86400.0_rk
   ! below are conversion rates intended to be multiplied
   ! and converts e.g. N to C in the case of Nmmol_to_Cmmol 
   real(rk), parameter :: Nmmol_to_Cmmol = 6.625_rk                    ! redf(1)
   real(rk), parameter :: Pmmol_to_Cmmol = 106.0_rk                    ! redf(2)
   real(rk), parameter :: Simmol_to_Cmmol = 6.625_rk                   ! redf(3)
   real(rk), parameter :: Pmmol_to_Nmmol = 16.0_rk                     ! redf(4)
   real(rk), parameter :: Simmol_to_Nmmol = 1.0_rk                     ! redf(5)
   real(rk), parameter :: Cmmol_to_Cmg = 12.01_rk                      ! redf(6)
   real(rk), parameter :: O2ml_l_to_O2mmol_m3 = 44.6608009_rk          ! redf(7)
   real(rk), parameter :: Nmmol_to_Nmg = 14.007_rk                     ! redf(8)
   real(rk), parameter :: Pmmol_to_Pmg = 30.97_rk                      ! redf(9)
   real(rk), parameter :: Simmol_to_Simg = 28.09_rk                    ! redf(10)
    
   real(rk), parameter :: Cmmol_to_Nmmol = 1._rk/6.625_rk              ! redf(11)
   real(rk), parameter :: Cmmol_to_Pmmol = 1._rk/106.0_rk              ! redf(12)
   real(rk), parameter :: Cmmol_to_Simmol = 1._rk/6.625_rk             ! redf(13)
   real(rk), parameter :: Nmmol_to_Pmmol = 1._rk/16.0_rk               ! redf(14)
   real(rk), parameter :: Nmmol_to_Simmol = 1.0_rk                     ! redf(15)
   real(rk), parameter :: Cmg_to_Cmmol = 1._rk/12.01_rk                ! redf(16)
   real(rk), parameter :: O2mmol_m3_to_O2ml_l = 1._rk/44.6608009_rk    ! redf(17)
   real(rk), parameter :: Nmg_to_Nmmol = 1._rk/14.007_rk               ! redf(18)
   real(rk), parameter :: Pmg_to_Pmmol = 1._rk/30.97_rk                ! redf(19)
   real(rk), parameter :: Simg_to_Simmol = 1._rk/28.09_rk              ! redf(20)

   real(rk):: frr = 0.4_rk                    ! default fraction of organic carbon from detritus remineralization
   logical  :: couple_co2 = .false.           ! default no CO2 coupling
   real(rk) :: prevent_loss_P = 0.5           ! default minimum biomass where loss terms are turned off for phytoplankton
   real(rk) :: prevent_loss_Z = 0.05          ! default minimum biomass where loss terms are turned off for zooplankton
   logical  :: model_has_silicifier = .false. ! if diatoms are included in the model, set true in fabm.yaml
   logical  :: model_has_calcifier = .false.  ! if coccoliths are included in the model, set true in fabm.yaml
   logical  :: depth_dependent_sinking_speed = .false.  ! optional organic matter sinking speed as a function of depth
                                                        ! see organic_matter.F90, opal.F90, caco3.F90 for implementation
   logical  :: use_niva_ersem_oxygen_exchange = .false. ! optional oxygen exchange following NIVA ERSEM implementation in oxygen.F90
   logical  :: use_temp_dependency_phy = .false.        ! optional temperature dependency for P production, if false, Tdep = 1.0_rk, see phy.F90
   logical  :: use_geider_PI_curve     = .false.        ! optional use of Geider's PI curve for P photoproduction, if false, use Yumruktepe et al., 2023 ECOSMO II(CHL) formulation
   logical  :: use_prey_switching      = .false.        ! optional adaptive prey switching (Murdoch 1969)
   logical  :: use_slp_egest_paradigm  = .false.        ! optional explicit sloppy feeding and waste routing (Steinberg & Landry 2017)
   logical  :: use_virtual_calcite     = .false.        ! optional virtual calcite paradigm where calcite is only formed upon mortality and grazing
   logical  :: use_bact_nutrient_limitation = .false.   ! optional implicit bacterial nutrient limitation on remineralization
   logical  :: use_community_sinking = .false.          ! global switch for community composition dependent sinking rates

   type,extends(type_base_model), public  :: type_ecosmo_shared
   contains
      procedure :: initialize
   end type type_ecosmo_shared

   private initialize
   contains
      subroutine initialize(self,configunit)
         class (type_ecosmo_shared), intent(inout),target  :: self
         integer,  intent(in) :: configunit
    
         call self%get_parameter( frr,                           "frr",                             "-",          "fraction of dissolved from det.", default=frr)
         call self%get_parameter( couple_co2,                    "couple_co2",                      "",           "switch coupling to carbonate module", default=couple_co2)
         call self%get_parameter( prevent_loss_P,                "prevent_loss_P",                  "mgC/m3",     "P biomass low threshold where loss terms are stopped for survival",  default=prevent_loss_P) 
         call self%get_parameter( prevent_loss_Z,                "prevent_loss_Z",                  "mgC/m3",     "Z biomass low threshold where loss terms are stopped for survival",  default=prevent_loss_Z) 
         call self%get_parameter( model_has_silicifier,          "model_has_silicifier",            "",           "global switch for opal", default=model_has_silicifier)
         call self%get_parameter( model_has_calcifier,           "model_has_calcifier",             "",           "global switch for calcite", default=model_has_calcifier)
         call self%get_parameter( depth_dependent_sinking_speed, "depth_dependent_sinking_speed",   "",           "global switch for depth dependent sinking speed", default=depth_dependent_sinking_speed)
         call self%get_parameter( use_niva_ersem_oxygen_exchange, "use_niva_ersem_oxygen_exchange", "",           "global switch for alternative oxygen exchange", default=use_niva_ersem_oxygen_exchange)
         call self%get_parameter( use_temp_dependency_phy,       "use_temp_dependency_phy",         "",           "global switch for phytoplankton temperature dependency", default=use_temp_dependency_phy)
         call self%get_parameter( use_geider_PI_curve,           "use_geider_PI_curve",             "",           "global switch for phytoplankton Geider PI curve", default=use_geider_PI_curve)
         call self%get_parameter( use_prey_switching,            "use_prey_switching",              "",           "global switch for adaptive prey switching", default=use_prey_switching)
         call self%get_parameter( use_slp_egest_paradigm,        "use_slp_egest_paradigm",          "",           "global switch for sloppy feeding paradigm", default=use_slp_egest_paradigm)
         call self%get_parameter( use_virtual_calcite,           "use_virtual_calcite",             "",           "global switch for virtual calcite paradigm", default=use_virtual_calcite)
         call self%get_parameter( use_bact_nutrient_limitation,  "use_bact_nutrient_limitation",    "",           "global switch for implicit bacterial nutrient limitation", default=use_bact_nutrient_limitation)
         call self%get_parameter( use_community_sinking,         "use_community_sinking",           "",           "global switch for community sinking speed",              default=use_community_sinking)
      end subroutine initialize

end module ecosmo_shared
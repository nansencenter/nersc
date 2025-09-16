#include "fabm_driver.h"

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

   real(rk):: frr = 0.4_rk
   logical  :: couple_co2 = .false.
!   logical  :: use_chl = .false.
!   logical  :: use_chl_in_PI_curve = .false.
!   logical  :: use_calcifier = .false.
!   real(rk) :: phyto_turn_off_loss_below_this = 0.5_rk
!   real(rk) :: zoo_turn_off_loss_below_this = 0.05_rk
   real(rk) :: light_att_chl = 0.04_rk
!   real(rk) :: light_att_phy = 0.04_rk / (Nmmol_to_Cmmol * Cmmol_to_Cmg)
   real(rk) :: prevent_loss_P = 0.5 ! mgC/m3
   real(rk) :: prevent_loss_Z = 0.05


   type,extends(type_base_model), public  :: type_ecosmo_shared
    contains

    procedure :: initialize
   end type type_ecosmo_shared

   private initialize
   contains
   subroutine initialize(self,configunit)
    !
    ! !INPUT PARAMETERS:
    class (type_ecosmo_shared), intent(inout),target  :: self
    integer,  intent(in) :: configunit
    !
    call self%get_parameter(frr, "frr", "-", "fraction of dissolved from det.", default=frr)
    call self%get_parameter( couple_co2, "couple_co2", "", "switch coupling to carbonate module", default=couple_co2)
!    call self%get_parameter( use_chl, "use_chl", "", "switch chlorophyll/c dynamics",default=use_chl)
!    call self%get_parameter( use_calcifier, "use_calcifier", "", "include caco3/sediment_caco3 (and coccoliths) in the model",default=use_calcifier)
!    call self%get_parameter( use_chl_in_PI_curve, "use_chl_in_PI_curve","","activated chl dependent light limitation",default=use_chl_in_PI_curve)
!    call self%get_parameter( phyto_turn_off_loss_below_this, "phyto_turn_off_loss_below_this", "-", "turn off loss terms below this concentration", default=phyto_turn_off_loss_below_this)
!    call self%get_parameter( zoo_turn_off_loss_below_this, "zoo_turn_off_loss_below_this", "-", "turn off loss terms below this concentration", default=zoo_turn_off_loss_below_this)
    call self%get_parameter( light_att_chl , 'light_att_chl', 'm**2/mgCHL', 'chl self-shading', default=light_att_chl )
!    call self%get_parameter( light_att_phy , 'light_att_phy', 'm**2/mmolN', 'phyto self-shading', default=light_att_phy, scale_factor=1.0_rk/(Nmmol_to_Cmmol * Cmmol_to_Cmg) )
    call self%get_parameter( prevent_loss_P,  'prevent_loss_P',  'mgC/m3', 'P biomass low threshold where loss terms are stopped for survival',  default=prevent_loss_P) 
    call self%get_parameter( prevent_loss_Z,  'prevent_loss_Z',  'mgC/m3', 'Z biomass low threshold where loss terms are stopped for survival',  default=prevent_loss_Z) 

  end subroutine initialize

  end module ecosmo_shared
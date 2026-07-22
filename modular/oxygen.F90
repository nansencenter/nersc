#include "fabm_driver.h"

! --------- Change Log ---------- !
! Veli Çağlar Yumruktepe (VCY)
!
! VCY - 27/05/2026
! This is the first iteration of the modular version of ECOSMO II(CHL).
! It is based on the ECOSMO II(CHL) code used for Copernicus ARC MFC 2026 operational model code and parameters.
! ------------------------------- !

module ecosmo_oxygen
    use fabm_types
    use fabm_expressions
    use ecosmo_shared
    implicit none
    private
    public type_ecosmo_oxygen
    type,extends(type_base_model), public  :: type_ecosmo_oxygen
        type (type_state_variable_id)         :: id_c
        type (type_state_variable_id)         :: id_no3, id_nh4, id_alk
        type (type_dependency_id)             :: id_temp, id_salt
        type (type_horizontal_dependency_id) :: id_icearea

    contains
        procedure :: initialize
        procedure :: do
        procedure :: do_surface
    end type type_ecosmo_oxygen

contains
    subroutine initialize(self,configunit)
        class (type_ecosmo_oxygen), intent(inout),target  :: self
        integer,  intent(in) :: configunit
        call self%register_state_variable(self%id_c, 'c', 'mmol/m3', 'oxygen',minimum=0.0_rk) ! Why is this not allowed to be negative? 
        call self%register_state_dependency(self%id_no3, 'no3', 'mgC/m3', 'nitrate')
        call self%register_state_dependency(self%id_nh4, 'nh4', 'mgC/m3', 'ammonium')
        call self%register_dependency(self%id_temp,standard_variables%temperature)
        call self%register_dependency(self%id_salt,standard_variables%practical_salinity)
        call self%register_dependency(self%id_icearea,standard_variables%ice_area_fraction)
        if (couple_co2) then
            call self%register_state_dependency(self%id_alk, 'alk','mmol m-3','alkalinity budget')
        end if
    end subroutine initialize

    subroutine do(self,_ARGUMENTS_DO_)
        class (type_ecosmo_oxygen),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_
        real(rk) :: oxy, no3, nh4, temp
        real(rk) :: Onitr, bioom1, bioom6, nitrification
        real(rk) :: rhs_oxy, rhs_amm, rhs_nit, rhs_alk

        _LOOP_BEGIN_

        _GET_(self%id_c,oxy)
        _GET_(self%id_no3,no3)
        _GET_(self%id_nh4,nh4)
        _GET_(self%id_temp, temp)    

        Onitr = 0.01_rk * O2ml_l_to_O2mmol_m3 !according to Neumann  (Onitr in mlO2/l see also Stigebrand and Wulff)
        bioom1 = 0.0_rk
        if (oxy > 0) then
            bioom1 = 0.1_rk/sedy0 * exp(temp*0.11_rk) * oxy/(Onitr+oxy)
        end if

        ! Nitrification
        nitrification = bioom1 * nh4 
        ! Ammonium change
        rhs_amm = -nitrification
        _ADD_SOURCE_(self%id_nh4, rhs_amm )
        ! Nitrate change
        rhs_nit = nitrification
        _ADD_SOURCE_(self%id_no3, rhs_nit )

        ! Oxygen change
        rhs_oxy = -2.0_rk * nitrification * Cmg_to_Cmmol * Cmmol_to_Nmmol 
        _ADD_SOURCE_( self%id_c, rhs_oxy ) 

        ! Alkalinity change
        if (couple_co2) then
            rhs_alk = -2.0_rk * nitrification * Cmg_to_Cmmol * Cmmol_to_Nmmol
            _ADD_SOURCE_(self%id_alk, rhs_alk )
         end if

        _LOOP_END_
    end subroutine do

    subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
        class (type_ecosmo_oxygen),intent(in) :: self
        _DECLARE_ARGUMENTS_DO_SURFACE_
 
        real(rk) :: o2flux, T, tr, S, o2sat, oxy, icearea
        real(rk) :: wnd, OSAT, ko2o, FAIRO2, sc
        _HORIZONTAL_LOOP_BEGIN_
 
        _GET_(self%id_temp,T)
        _GET_(self%id_salt,S)
        _GET_(self%id_c,oxy)
        _GET_SURFACE_(self%id_icearea,icearea)

        ! Optional oxygen exchange formulation (default is Daewel & Shrum 2013)
        if (use_niva_ersem_oxygen_exchange) then
            ! NIVA ERSEM oxygen exchange
            S = max(S, 0.0_rk)
            wnd = max(wnd, 0.0_rk)
            OSAT = oxygen_saturation_concentration(self,T,S)

            ! New formulation for the Schmidt number for O2 following Wanninkhof 2014
            T = max(min(T,40.0_rk), -2.0_rk)
            sc = 1920.4_rk - 135.6_rk*T + 5.2122_rk*T**2 - 0.10939_rk*T**3 + 0.00093777_rk*T**4

            ko2o = 0.251_rk * (wnd**2) * (sc/660._rk)**(-0.5_rk) ! Wanninkhof 2014

            ! units of ko2 converted from cm/hr to m/s
            ko2o = ko2o/360000._rk

            ! ice area must be set to 0.0 for ocean by the physics models
            ko2o = max(0._rk, ( 1.0_rk - icearea )) * ko2o
            o2flux = ko2o * (OSAT - oxy)
        else
            ! Daewel & Shrum 2013 oxygen exchange
            ! Oxygen saturation micromol/liter__(Benson and Krause, 1984)
            tr = 1.0_rk/(T + 273.15_rk)
            o2sat= exp(- 135.90205_rk              &
                + (1.575701e05_rk ) * tr               &
                - (6.642308e07_rk ) * tr**2            &
                + (1.243800e10_rk) * tr**3            &
                - (8.621949e11_rk) * tr**4            &
                - S*(0.017674_rk-10.754_rk*tr+2140.7_rk*tr**2)  )
 
                !   o2flux = 5._rk/sedy0 * (o2sat - oxy)
                o2flux = 1.0_rk/sedy0 * (o2sat - oxy) * ( 1.0_rk - icearea )
        end if 

        _ADD_SURFACE_FLUX_(self%id_c,o2flux)

        _HORIZONTAL_LOOP_END_

    end subroutine do_surface

   ! Optional oxygen saturation concentration function adapted from NIVA ERSEM model
   function oxygen_saturation_concentration(self,ETW,X1X) result(OSAT)
      class (type_ecosmo_oxygen), intent(in) :: self
      real(rk),                 intent(in) :: ETW,X1X
      real(rk)                             :: OSAT

      real(rk),parameter :: A1 = -173.4292_rk
      real(rk),parameter :: A2 = 249.6339_rk
      real(rk),parameter :: A3 = 143.3483_rk
      real(rk),parameter :: A4 = -21.8492_rk
      real(rk),parameter :: B1 = -0.033096_rk
      real(rk),parameter :: B2 = 0.014259_rk
      real(rk),parameter :: B3 = -0.0017_rk
      real(rk),parameter :: R = 8.3145_rk
      real(rk),parameter :: P = 101325_rk
      real(rk),parameter :: T = 273.15_rk

      ! volume of an ideal gas at standard temp (0C) and pressure (1 atm)
      real(rk),parameter :: VIDEAL = (R * 273.15_rk / P) *1000._rk

      real(rk)           :: ABT

      ! calc absolute temperature
      ABT = ETW + T

      ! calc theoretical oxygen saturation for temp + salinity
      ! From WEISS 1970 DEEP SEA RES 17, 721-735.
      ! units of ln(ml(STP)/l)
      OSAT = A1 + A2 * (100._rk/ABT) + A3 * log(ABT/100._rk) &
               + A4 * (ABT/100._rk) &
               + X1X * ( B1 + B2 * (ABT/100._rk) + B3 * ((ABT/100._rk)**2))

      ! convert units to ml(STP)/l then to mMol/m3
      OSAT = exp( OSAT )
      OSAT = OSAT * 1000._rk / VIDEAL
   end function

end module
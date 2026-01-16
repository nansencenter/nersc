! Copyright 2017 Helmholtz-Zentrum Geestacht
! Licensed under the Apache License, Version 2.0
!  (http://www.apache.org/licenses/LICENSE-2.0)
! Author(s): Richard Hofmeister and Deborah Benkort
!
! Modified and moved from nersc to nersc
! Jan 2026, Shuang Gao

#include "fabm_driver.h"

   module shared_variables
      logical :: initialized_icethickness 
   end module shared_variables

!  !MODULE: fabm_nersc_icealgae_G --- new biogeochemical model
!
!  !INTERFACE:
   module fabm_nersc_icealgae_G

! !DESCRIPTION:
!
! This ecosystem model has to be developed
!

! !USES:
   use fabm_types
   use fabm_expressions

   implicit none

!  default: all is private
   private
!
! !PUBLIC MEMBER FUNCTIONS:
   public type_nersc_icealgae_G
!
! !PRIVATE DERIVED TYPES:
   type,extends(type_base_model) :: type_nersc_icealgae_G
!     Variable identifiers
      type (type_surface_state_variable_id)          :: id_Ino3, id_Inh4 , id_Ipho, id_Isil
      type (type_surface_state_variable_id)          :: id_Ialg, id_Idet !,id_Iatt
      type (type_state_variable_id)                  :: id_Wno3, id_Wnh4, id_Wpho, id_Wsil
      type (type_state_variable_id)                  :: id_Wdia, id_Wfla, id_Wdetf, id_Wdet, id_Wmeso, id_Wmicro
      type (type_horizontal_diagnostic_variable_id)  :: id_WZsonIA, id_WZlonIA
      type (type_horizontal_diagnostic_variable_id)  :: id_primprod, id_l_limit, id_n_limit, id_p_limit, id_s_limit, id_denit, id_dhice, id_entIA, id_meltIA, id_entIdet, id_meltIdet, id_sal_limit, id_temp_limit, id_Iatt
      type (type_horizontal_diagnostic_variable_id)  :: id_entno3, id_entnh4, id_entpho, id_entsil, id_meltno3, id_meltnh4, id_meltpho, id_meltsil, id_br_release, id_br_flux, id_br_no3, id_br_nh4, id_br_pho, id_br_sil
      type (type_horizontal_diagnostic_variable_id)  :: id_diffno3, id_diffnh4, id_diffpho, id_diffsil
      type (type_horizontal_diagnostic_variable_id)  :: id_upno3, id_upnh4, id_uppho, id_upsil
      type (type_horizontal_diagnostic_variable_id)  :: id_remineralisation, id_mortality, id_grazing, id_lightIA, id_pno3, id_hice, id_Zice, id_Zbot, id_Zwater, id_iconc, id_muTS !, id_lightBOT
      type (type_horizontal_diagnostic_variable_id)  :: id_no3W, id_nh4W, id_phoW, id_silW
      type (type_dependency_id)                      :: id_temp, id_sal, id_thick_cell
      type (type_horizontal_dependency_id)           :: id_air_temp, id_Fsw
      type (type_horizontal_dependency_id)           :: id_h, id_hs, id_i_con
      type (type_horizontal_dependency_id)           :: id_dh_growth

!     Model parameters
      real(rk) :: muIA
      real(rk) :: TctrlIA
      real(rk) :: aaIA
      real(rk) :: rM
      real(rk) :: rM2
      real(rk) :: rIalg
      real(rk) :: rGr
      real(rk) :: Pm
      real(rk) :: rNO3
      real(rk) :: rPO4
      real(rk) :: rSi
      real(rk) :: rnit
      real(rk) :: rrem
      real(rk) :: vn
      real(rk) :: C0
      real(rk) :: Catt_ice
      real(rk) :: Catt_sno
      real(rk) :: g_ia
      real(rk) :: ice_alb
      real(rk) :: sno_alb
      real(rk) :: Cio
      real(rk) :: diffcoef
      real(rk) :: rgrow
      real(rk) :: rmelt
      real(rk) :: fr
      real(rk) :: Zbot
      real(rk) :: dice
      real(rk) :: dwater 
      real(rk) :: num_melting 
      real(rk) :: gammaZsp
      real(rk) :: gammaZlp
      real(rk) :: prefZsPs
      real(rk) :: prefZlPs
      real(rk) :: prefZsPl
      real(rk) :: prefZlPl
      real(rk) :: prefZsD
      real(rk) :: prefZlD
      real(rk) :: prefZlZs
      real(rk) :: GrZlP
      real(rk) :: GrZsP
      real(rk) :: GrZlZ
      real(rk) :: Rg
      real(rk) :: surv_prop
      logical  :: couple_ice
      logical  :: turn_on_additional_diagnostics

      contains

!     Model procedures
      procedure :: initialize
      procedure :: do_surface

   end type type_nersc_icealgae_G

! !PRIVATE DATA MEMBERS:
   real(rk), parameter :: one_pr_day = 1.0_rk/86400.0_rk
   real(rk), parameter :: one_pr_hour = 1.0_rk/3600.0_rk
   real(rk)            :: redf(20) = 0.0_rk
   
   real(rk), parameter :: sedy0 = 86400.0_rk
   real(rk), parameter :: eps = 1.0e-6_rk

   real(rk) :: baclin
   !namelist /CORE/dt
!-----------------------------------------------------------------------
   !type (type_surface_standard_variable),parameter :: surface_attenuation_coef_of_photosynthetic_radiative_flux = type_surface_standard_variable(name='surface_attenuation_coef_of_photosynthetic_radiative_flux',units='-',aggregate_variable=.true.)

   contains

!-----------------------------------------------------------------------
!
! Initialise the icealgae model
!
   subroutine initialize(self,configunit)
!
! !DESCRIPTION:
!  Here, the icealgae namelist is read and the variables exported
!  by the model are registered with FABM.
!
   use shared_variables
   class (type_nersc_icealgae_G),intent(inout),target  :: self
   integer,                  intent(in)            :: configunit
   real(rk)                                        :: ws

!  Set Redfield ratios:
   redf(1)  = 6.625_rk      !molC/molN
   redf(2)  = 106.0_rk      !C_P
   redf(3)  = 6.625_rk      !C_SiO
   redf(4)  = 16.0_rk       !N_P
   redf(5)  = 1.0_rk        !N_SiO
   redf(6)  = 12.01_rk      !gC/molC
   redf(7)  = 44.6608009_rk !O2mm_ml
   redf(8)  = 14.007_rk     !N_Nmg
   redf(9)  = 30.97_rk      !P_Pmg
   redf(10) = 28.09_rk      !Si_Simg    


      !open(15,file='/work/gg0877/KST/dbenkort/schism-setups/arctic/run_schism/new_param.nml',delim='apostrophe',status='old')
      !read(15,nml=CORE)
      !baclin = 600.0_rk ! TP2
      baclin = 300.0_rk ! TP5
      !close(15)
!------------------------------------------------------------------------------------

   call self%get_parameter(self%muIA,       'muIA',       '1/day',                  'ice algae maximum growth rate',                  default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%TctrlIA,    'TctrlIA',    '1/degC',                 'temperature growth dependency',                  default=0.001_rk)
   call self%get_parameter(self%rM,         'rM',         '1/day',                  'ice algae mortality rate',                       default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rM2,         'rM2',        '1/day',                 'ice algae additional mortality rate',            default=0.3_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rIalg,       'rIalg',       'mgC/m**2',              'ice algae saturation coefficient', default=0.1_rk)
   call self%get_parameter(self%rGr,        'rGr',        '',                  'ice algae grazing rate',                         default=0.001_rk)
   call self%get_parameter(self%aaIA,       'aaIA',       'mgC/mgchla/h/Einst/m/s', 'photosynthesis ef-cy',                           default=0.5_rk,   scale_factor=one_pr_hour)
   call self%get_parameter(self%Pm,         'Pm',         'mgC/mgchla/h',           'maximum photosynhtetic rate',                    default=0.5_rk,   scale_factor=one_pr_hour)
   call self%get_parameter(self%rNO3,       'rNO3',       'mmolN/m**3',             'NO3 half saturation',                            default=0.5_rk)
   call self%get_parameter(self%rPO4,       'rPO4',       'mmolP/m**3',             'PO4 half saturation',                            default=0.5_rk)
   call self%get_parameter(self%rSi,        'rSi',        'mmolSi/m**3',            'SiO2 half saturation',                           default=0.5_rk)
   call self%get_parameter(self%rnit,       'rnit',       '1/day',                  'nitrification rate',                             default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%rrem,       'rrem',       '1/day',                  'remineralisation rate',                          default=0.001_rk, scale_factor=one_pr_day)
   call self%get_parameter(self%vn,         'vn',         'mmolN/m**3',             'preferential uptake of nitrate half saturation', default=0.5_rk)
   call self%get_parameter(self%gammaZlp,   'gammaZlp',   'unitless',               'part of available IA for mesozooplankton grazing' ,  default=0.2_rk)
   call self%get_parameter(self%gammaZsp,   'gammaZsp',   'unitless',               'part of available IA for microzooplankton grazing' ,  default=0.2_rk)
   call self%get_parameter( self%GrZlP, 'GrZlP',       '1/day',      'Grazing rate Zl on Phyto',        default=0.80_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZsP, 'GrZsP',       '1/day',      'Grazing rate Zs on Phyto',        default=1.00_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%GrZlZ, 'GrZlZ',       '1/day',      'Grazing rate Zl on Zs',           default=0.50_rk,  scale_factor=1.0_rk/sedy0)
   call self%get_parameter( self%Rg, 'Rg',          'mmolN/m**3', 'Zs, Zl half saturation',          default=0.50_rk,  scale_factor=redf(1)*redf(6))

   call self%get_parameter( self%prefZsPs,  'prefZsPs',   '-',          'Grazing preference Zs on Ps',     default=0.70_rk)
   call self%get_parameter( self%prefZsPl,  'prefZsPl',   '-',          'Grazing preference Zs on Pl',     default=0.25_rk)
   call self%get_parameter( self%prefZsD,   'prefZsD',    '-',          'Grazing preference Zs on Det.',   default=0.00_rk)
   call self%get_parameter( self%prefZlPs,  'prefZlPs',   '-',          'Grazing preference Zl on Ps',     default=0.10_rk)
   call self%get_parameter( self%prefZlPl,  'prefZlPl',   '-',          'Grazing preference Zl on Pl',     default=0.85_rk)
   call self%get_parameter( self%prefZlZs,  'prefZlZs',   '-',          'Grazing preference Zl on Zs',     default=0.15_rk)
   call self%get_parameter( self%prefZlD,   'prefZlD',    '-',          'Grazing preference Zl on Det.',   default=0.00_rk)

   call self%get_parameter(self%C0,         'C0',         'unitless',               'part of incoming radiation absorbed',            default=0.3_rk)
   call self%get_parameter(self%Catt_ice,   'Catt_ice',   '1/m',                    'ice attenuation coefficent',                     default=0.5_rk)                
   call self%get_parameter(self%Catt_sno,   'Catt_sno',   '1/m',                    'snow attenuation coefficent',                    default=0.5_rk)                
   !call self%get_parameter(self%g_ia,   'g_ia',   '1/(mmolC/m2)/m',             'e-folding depth of radiation by ice algea concentration', default=3.0_rk)              
   call self%get_parameter(self%g_ia,   'g_ia',   'm2/mgC',                      'ice algae shading', default=0.002_rk)
   call self%get_parameter(self%ice_alb,    'ice_alb',    'unitless',               'albedo in the ice',                              default=0.5_rk)                
   call self%get_parameter(self%ice_alb,    'ice_alb',    'unitless',               'albedo in the ice',                              default=0.5_rk)                
   call self%get_parameter(self%sno_alb,    'sno_alb',    'unitless',               'albedo in the snow',                             default=0.5_rk)                
   call self%get_parameter(self%Cio,        'Cio',        'm**2/s',                 'drag coefficient',                               default=0.0_rk)                
   call self%get_parameter(self%diffcoef,   'diffcoef',   'm**2/s',                 'molecular diffusion coefficient',                default=0.0_rk)                
   call self%get_parameter(self%rgrow,      'rgrow',      'm/d',                    'growth rate of ice',                             default=0.0_rk)                
   call self%get_parameter(self%fr,         'fr',         'unitless',               'relative brine fraction',                        default=0.0_rk)                
   call self%get_parameter(self%Zbot,       'Zbot',       'm',                      'ice bottom layer thickness',                     default=0.0_rk)                
   call self%get_parameter(self%dice,       'dice',       'kg/m**3',                'density of sea ice',                             default=0.0_rk)
   call self%get_parameter(self%dwater,     'dwater',     'kg/m**3',                'density of the sea water',                       default=0.0_rk)  
   call self%get_parameter(self%couple_ice, 'couple_ice', 'true/false',             'switch coupling to ice',                         default=.false.)
   call self%get_parameter( self%turn_on_additional_diagnostics, 'turn_on_additional_diagnostics','','activates additional diagnostics for model debugging',default=.false.)
   call self%get_parameter(self%num_melting,'num_melting','1/day',                  'numerical melting rate',                         default=1.0_rk, scale_factor=one_pr_day)
   ! call self%get_parameter(self%prefZs,     'prefZs',     'unitless',               'Small zooplankton preference on ice algae',                       default=0.0_rk)  
   !call self%get_parameter(self%prefZl,     'prefZl',     'unitless',               'Large zooplankton preference on ice algae',                       default=0.0_rk)  
   call self%get_parameter(self%surv_prop,  'surv_prop',  'unitless',               'Proportion of diatom surviving in the sea-ice during entrapment',                       default=0.5_rk)  
   call self%register_surface_state_variable(self%id_Ialg,    'Ialg',   'mgC/m**2',   'ice algae', initial_value=0.1_rk, minimum=0.0_rk)

   call self%register_surface_state_variable(self%id_Idet,    'Idet',   'mgC/m**2',   'detritus',       minimum=0.0_rk)

   call self%register_surface_state_variable(self%id_Ino3,    'Ino3',   'mmolN/m**2',   'nitrate',        minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Inh4,    'Inh4',   'mmolN/m**2',   'nitrate',        minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Ipho,    'Ipho',   'mmolP/m**2',   'phosphate',      minimum=0.0_rk)
   call self%register_surface_state_variable(self%id_Isil,    'Isil',   'mmolSi/m**2',  'silicate',       minimum=0.0_rk)
   call self%register_diagnostic_variable(self%id_Iatt,    'Iatt',   '-',  'attenuation of ice algae', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_primprod,         'priIA',    'mgC/m2/s',      'primary production', source=source_do_surface)
   if (self%turn_on_additional_diagnostics) then
   call self%register_diagnostic_variable(self%id_remineralisation, 'remIA',         'mgC/m2/s',      'remineralisation', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_mortality,        'morIA',   'mgC/m2/s',      'mortality IA', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_grazing,          'grazm',     'mgC/m2/s',      'mortality IA by grazing', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_l_limit,          'l_lim',     'unitless',      'light limitation factor for ia', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_sal_limit,        's_lim',   'unitless',      'salinity limitation factor for ia', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_temp_limit,       't_lim',  'unitless',      'temperature limitation factor for ia', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_n_limit,          'N_lim',     'unitless',      'N limitation factor for ia', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_p_limit,          'P_lim',     'unitless',      'phosphate limitation factor for ia', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_s_limit,          'Silim',     'unitless',      'silicate limitation factor for ia', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_lightIA,          'l_avi',       'unitless',      'light available for IA', source=source_do_surface)
   !call self%register_diagnostic_variable(self%id_lightBOT,         'lBOT',      'unitless',      'light available at the surface under the ice', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_muTS,         'muTS',      'unitless',      'IA growth rate after modified by T S lim', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_iconc,            'iconc',         'unitless',      'ice concentration in the domain cell', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_hice,             'hice',          'm',             'thickness of ice pack', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_dhice,            'dhice',        'm/s',             'rate of change of the ice thickness', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_upno3,            'upno3',    'mmolN/m**2',    'no3 uptake by algae', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_upnh4,            'upnh4',    'mmolN/m**2',    'nh4 uptake by algae', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_uppho,            'uppho',    'mmolN/m**2',    'pho uptake by algae', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_upsil,            'upsil',    'mmolN/m**2',    'sil uptake by algae', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_entIA,            'entIA',     'mgC/m**2',      'algae entrapment', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_meltIA,           'melIA',       'mgC/m**2',      'algae release', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_entIdet,          'endet',   'mgC/m**2',      'detritus entrapment', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_meltIdet,         'medet',     'mgC/m**2',      'detritus release', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_pno3,             'pno3',          'unitless',      '', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_entno3,           'enno3',    'mmolN/m**2',    'no3 entrapment', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_meltno3,          'meno3',      'mmolN/m**2',    'no3 release', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_entnh4,           'ennh4',    'mmolN/m**2',    'nh4 entrapment', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_meltnh4,          'menh4',      'mmolN/m**2',    'nh4 release', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_entpho,           'enpho',    'mmolP/m**2',    'pho entrapment', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_meltpho,          'mepho',      'mmolP/m**2',    'pho release', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_entsil,           'ensil',    'mmolSi/m**2',   'sil entrapment', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_meltsil,          'mesil',      'mmolSi/m**2',   'sil release', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_diffno3,          'dfno3', 'mmolN/m**2',    'no3 diffusive flux', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_diffnh4,          'dfnh4', 'mmolN/m**2',    'nh4 diffusive flux', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_diffpho,          'dfpho', 'mmolN/m**2',    'pho diffusive flux', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_diffsil,          'dfsil', 'mmolN/m**2',    'sil diffusive flux', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_br_no3,           'brno3',        'mmolN/m**3',    'no3 concentration in the brine', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_br_nh4,           'brnh4',        'mmolN/m**3',    'nh4 concentration in the brine', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_br_pho,           'brpho',        'mmolP/m**3',    'pho concentration in the brine', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_br_sil,           'brsil',        'mmolSi/m**3',   'sil concentration in the brine', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_Zice,             'Zice',          'm',             'ice without biological activity thickness', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_Zbot,             'Zbot',          'm',             'active bottom layer thickness', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_Zwater,           'Zwat',        'm',             '1st water layer thickness', source=source_do_surface)

   call self%register_diagnostic_variable(self%id_no3W,             'no3_w',     'mgC/m**3',     'no3 come from water', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_nh4W,             'nh4_w',     'mgC/m**3',     'nh4 come from water', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_phoW,             'pho_w',     'mgC/m**3',     'phosphate come from water', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_silW,             'sil_w',     'mgC/m**3',     'silicate come from water', source=source_do_surface)
   end if

   ! Register dependencies
   call self%register_dependency(self%id_temp, standard_variables%temperature)
   call self%register_dependency(self%id_Fsw, standard_variables%surface_downwelling_shortwave_flux)
   call self%register_dependency(self%id_thick_cell, standard_variables%cell_thickness)
   call self%register_dependency(self%id_sal, standard_variables%practical_salinity)
   call self%register_dependency(self%id_h, standard_variables%ice_thickness)
   call self%register_dependency(self%id_i_con, standard_variables%ice_area_fraction)
   call self%register_dependency(self%id_dh_growth, standard_variables%dh_growth)
   !call self%register_dependency(standard_variable(name='time_step', units='s')

   !call self%register_dependency(self%id_hs, standard_variables%snow_thickness)

   call self%register_state_dependency(self%id_Wno3,   'Wno3',         'mgC/m**3',  'nitrate in surface water column')
   call self%register_state_dependency(self%id_Wnh4,   'Wnh4',         'mgC/m**3',  'ammonium in surface water column')
   call self%register_state_dependency(self%id_Wpho,   'Wpho',         'mgC/m**3',  'phosphate in surface water column')
   call self%register_state_dependency(self%id_Wsil,   'Wsil',         'mgC/m**3',  'silicate in surface water column')
   call self%register_state_dependency(self%id_Wfla,   'Wfla',         'mgC/m**3',    'flagellate in surface water column')
   call self%register_state_dependency(self%id_Wdia,   'Wdia',         'mgC/m**3',    'diatom in surface water column')
   call self%register_state_dependency(self%id_Wdetf,  'Wdetf',        'mgC/m**3',    'ice-algae detritus in surface water column')
   call self%register_state_dependency(self%id_Wdet,   'Wdet',         'mgC/m**3',    'detritus in surface water column')
   call self%register_state_dependency(self%id_Wmicro, 'Wmicro',       'mgC/m**3',    'microzooplankton in surface water column')
   call self%register_state_dependency(self%id_Wmeso,  'Wmeso',        'mgC/m**3',    'mesozooplankton in surface water column')
   if (self%turn_on_additional_diagnostics) then
   call self%register_diagnostic_variable(self%id_WZsonIA, 'WZsIA',  '-', 'grazing term of microzooplankton on ice algae', source=source_do_surface)
   call self%register_diagnostic_variable(self%id_WZlonIA, 'WZlIA',  '-', 'grazing term of mesozooplankton on ice algae', source=source_do_surface)
   end if
   !call self%add_to_aggregate_variable(standard_variables%surface_attenuation_coef_of_photosynthetic_radiative_flux,self%id_Iatt)

   return

   end subroutine initialize

! Surface fluxes for the icealgae model

   subroutine do_surface(self,_ARGUMENTS_DO_SURFACE_)
   
   use shared_variables

   class (type_nersc_icealgae_G),intent(in) :: self
   _DECLARE_ARGUMENTS_DO_SURFACE_

   real(rk) :: Ino3, Inh4, Ipho, Isil
   real(rk) :: BR_no3, BR_nh4, BR_pho, BR_sil
   real(rk) :: temp, sal, cell_thick
   real(rk) :: Wno3, Wnh4, Wpho, Wsil, Wdia,Wfla,  Wdetf, Wdet, Wmicro, Wmeso, Widet, Wialg, Wino3
   real(rk) :: WZsonIA, WZlonIA, Fs, Fl
   real(rk) :: Ialg, Idet
   real(rk) :: up_Ino3, up_Ipho, up_Isil
   real(rk) :: rhs_Init, pno3, rhs_Iamm, rhs_Ipho, rhs_Isil, rhs_Ialg
   real(rk) :: rhs_IA, rhs_det
   real(rk) :: blight, tdep, sdep
   real(rk) :: rem, rem_amn, nit, IA_prod, Prod, IA_loss, IA_gr 
   real(rk) :: icethickness, iceconc, ice_dh, melt_ice
   real(rk) :: snowthickness, hbio
   real(rk) :: alb, Fsw, bio_par,Iatt
   real(rk) :: molN2gC, molP2gC, molS2gC
   real(rk) :: entrap_no3,melt_no3, congFlux_no3
   real(rk) :: entrap_nh4,melt_nh4, congFlux_nh4
   real(rk) :: entrap_pho,melt_pho, congFlux_pho
   real(rk) :: entrap_sil,melt_sil, congFlux_sil
   real(rk) :: diff_flux_no3, diff_flux_nh4, diff_flux_pho, diff_flux_sil
   real(rk) :: br_release, convFlux, Zice, Zbot
   real(rk) :: melt_IA, melt_det
   real(rk) :: entrap_IA, entrap_detf, entrap_dets,entrap_det
   real(rk) :: num_melting, test 

   if(self%couple_ice) then 
   _HORIZONTAL_LOOP_BEGIN_
   
  ! _GET_HORIZONTAL_(self%id_WZsonIA, WZsonIA)
   _GET_(self%id_Wno3, Wno3)
   _GET_(self%id_Wnh4, Wnh4)
   _GET_(self%id_Wpho, Wpho)
   _GET_(self%id_Wsil, Wsil)
   _GET_(self%id_Wdia, Wdia)
   _GET_(self%id_Wfla, Wfla)
   _GET_(self%id_Wdetf, Wdetf)
   _GET_(self%id_Wdet, Wdet)
   _GET_(self%id_Wmicro, Wmicro)
   _GET_(self%id_Wmeso, Wmeso)
   _GET_(self%id_temp, temp)
   _GET_(self%id_sal, sal)
   _GET_(self%id_thick_cell, cell_thick)
   _GET_HORIZONTAL_(self%id_Fsw, Fsw) 
   _GET_HORIZONTAL_(self%id_Ino3, Ino3)
   _GET_HORIZONTAL_(self%id_Inh4, Inh4)
   _GET_HORIZONTAL_(self%id_Ipho, Ipho)
   _GET_HORIZONTAL_(self%id_Isil, Isil)
   _GET_HORIZONTAL_(self%id_Ialg, Ialg)
   _GET_HORIZONTAL_(self%id_Idet, Idet)
   _GET_HORIZONTAL_(self%id_h, icethickness)
   _GET_HORIZONTAL_(self%id_i_con, iceconc)
   !_GET_HORIZONTAL_(self%id_Iatt, Iatt)
   _GET_HORIZONTAL_(self%id_dh_growth, ice_dh)

   !!!-------------------------------!!!
   !!!     ICE AND SNOW THICKNESS    !!!
   !!!-------------------------------!!!
   
   !if(icethickness > 0.01_rk) then
        ! If ice-algae is smaller than 0.01, then add up to 0.01
        ! _SET_SURFACE_ODE_: add surface source
   !     if(Ialg < 0.01_rk .AND. Ialg > 0.0_rk) then
   !          _SET_SURFACE_ODE_(self%id_Ialg, ((0.01-Ialg)/baclin))
   !          write(*,*) 'DEBUG1 Ialg set to 0.01', Ialg
   !     endif
   !else     
   !     _SET_SURFACE_ODE_(self%id_Ialg, (-Ialg)/baclin)
   !     write(*,*) 'DEBUG2 ice <0.01m release Ialg'
   !endif 


   !if(icethickness < 0.01_rk) then 
   !     snowthickness = 0.0_rk
   !endif
   !if(icethickness < 0.001_rk) then 
   !     write(*,*) 'DEBUG0 iceh=', icethickness
   !endif

   if(icethickness == 0.0_rk) then
        melt_ice = 0.0_rk
   else   
        melt_ice = ice_dh/icethickness !????icethickness can be very small!
   end if

   Zbot = min(self%Zbot, icethickness)
   Zice = icethickness - Zbot
   
   !!!-------------------------------!!!
   !!!       CONVERSION FACTOR       !!!
   !!!-------------------------------!!!

   ! molN/gC = 1/(molC/molN * gC/molC) 
   molN2gC=(1/(redf(1)*redf(6)))
   ! molP/gC = 1/(molC/molP * gC/molC)
   molP2gC=(1/(redf(2)*redf(6)))
   ! molSi/gC = 1/(molC/molSi * gC/molC) 
   molS2gC=(1/(redf(3)*redf(6)))

   !!!-------------------------------!!!
   !!!       LIGHT IN THE ICE        !!!
   !!!-------------------------------!!!
   !**! PAR available to ice algae
   !if(icethickness == 0.0_rk) then
   !     bio_par = 0.0_rk  !! if there is no ice there is no need to calculate bio_par
   !else
        !bio_par = Fsw * self%C0 * exp((-self%Catt_ice * icethickness) - (self%Catt_sno * snowthickness))
        !bio_par = bio_par * 0.45 * 4.91
        !! in HYCOM Fsw is already modified by ice and snow
        bio_par = Fsw * 0.45 * 4.91 !! 0.45 and 4.91 are used to convert W m-2 to muEinst m-2 s-1 
   !end if

   !!!-------------------------------!!!
   !!!       NUTRIENTS DYNAMIC       !!!
   !!!-------------------------------!!!

   !**! Diffusion flux
   if(icethickness <= 0.005_rk) then
	BR_no3 = 0.0_rk
	BR_nh4 = 0.0_rk
	BR_pho = 0.0_rk
	BR_sil = 0.0_rk

	diff_flux_no3 = 0.0_rk
	diff_flux_nh4 = 0.0_rk
	diff_flux_pho = 0.0_rk
	diff_flux_sil = 0.0_rk

   else 
	BR_no3 = Ino3 / self%fr / Zbot
	BR_nh4 = Inh4 / icethickness  
	BR_pho = Ipho / self%fr / Zbot
	BR_sil = Isil / self%fr / Zbot 

	diff_flux_no3 = -self%diffcoef * (BR_no3 - Wno3*molN2gC)/(Zbot) * self%fr  
	diff_flux_nh4 = -self%diffcoef * (BR_nh4 - Wnh4*molN2gC)/(Zbot) * self%fr 
	diff_flux_pho = -self%diffcoef * (BR_pho - Wpho*molP2gC)/(Zbot) * self%fr 
	diff_flux_sil = -self%diffcoef * (BR_sil - Wsil*molS2gC)/(Zbot) * self%fr  

   endif


   !**! Entrapment flux
   entrap_no3   =  molN2gC *(max(0.0_rk, ice_dh) * max(0.0_rk, Wno3))  * (self%dice/self%dwater) / self%fr 

   entrap_nh4   =  molN2gC * (max(0.0_rk, ice_dh) * max(0.0_rk, Wnh4))  * (self%dice/self%dwater)

   entrap_sil   =  molS2gC * (max(0.0_rk, ice_dh) * max(0.0_rk, Wsil))  * (self%dice/self%dwater) / self%fr

   entrap_pho   =  molP2gC * (max(0.0_rk, ice_dh) * max(0.0_rk, Wpho))  * (self%dice/self%dwater) / self%fr
   !!!-------------------------------!!!
   !!!         ALGAE DYNAMIC         !!!
   !!!-------------------------------!!!

   !**! Growth / Melt fluxes
   entrap_IA  = (max(0.0_rk, ice_dh) * max(0.0_rk, Wdia))  * (self%dice/self%dwater) 
   _SET_SURFACE_ODE_(self%id_Ialg, entrap_IA*self%surv_prop)
   
   !**! Loss terms
   IA_loss = max(0.0_rk, (self%rM * Ialg)) + self%rM2 * (Ialg/(Ialg + self%rIalg)) * Ialg
   
   !**! Nitrification rate
   nit = self%rnit * Inh4 
   
   !**! Nutrient limitation factors
   up_Ino3 = (BR_no3)/(self%rNO3 + (BR_no3))
   up_Ipho = (BR_pho)/(self%rPO4 + (BR_pho))
   up_Isil = (BR_sil)/(self%rSi + (BR_sil))

   !**! Light limitation
   if(icethickness == 0.0_rk)then 
	blight = 1.0_rk ! 0.0_rk
   else
	blight = max(1 - exp(-(self%aaIA * bio_par) / self%Pm), 0.0_rk)
   endif 

   !**! Temperature and salt dependence
   if(temp < 0.0_rk)then 
        Tdep = exp(self%TctrlIA * temp)
   else
        Tdep = 1.0_rk
   endif
   Sdep = exp(-(2.16 - 8.30 * (10**-5) * (sal ** 2.11) - 0.55 * log(sal)) ** 2)  

   !!!-------------------------------!!!
   !!!       DETRITUS DYNAMIC        !!!
   !!!-------------------------------!!!

   !**! Growth / Melt fluxes
   entrap_detf = (max(0.0_rk, ice_dh) * max(0.0_rk, Wdetf))  * (self%dice/self%dwater) 
   entrap_dets = (max(0.0_rk, ice_dh) * max(0.0_rk, Wdet))  * (self%dice/self%dwater) 
   entrap_det = entrap_detf + entrap_dets
   !**! Remineralisation rate
   rem = self%rrem * Idet

   !!!-------------------------------!!!
   !!!     BIOGEO RATE OF CHANGE     !!!
   !!!-------------------------------!!!
   
   !**! Ice Algae change (production and nutrient uptake)
   IA_prod = Tdep * Sdep * min(blight, up_Ino3, up_Ipho, up_Isil) 
   Prod = self%muIA * Ialg * IA_prod
   !if(Ialg>3000.0_rk.AND.IA_prod/=0.0_rk) then
   !   write(*,*) 'DEBUG6 Ialg=',Ialg
   !   write(*,*) 'DEBUG6 M2=',self%rM2 * (Ialg/(Ialg + self%rIalg)) * Ialg
   !   write(*,*) 'DEBUG6 IA_prod=',IA_prod
   !   write(*,*) 'DEBUG6 Tdep=',Tdep
   !   write(*,*) 'DEBUG6 Sdep=',Sdep
   !   write(*,*) 'DEBUG6 min=',min(blight, up_Ino3, up_Ipho, up_Isil)
   !   write(*,*) 'DEBUG6 blight=',blight
   !   write(*,*) 'DEBUG6 up_Ino3=',up_Ino3
   !   write(*,*) 'DEBUG6 up_Ipho=',up_Ipho
   !   write(*,*) 'DEBUG6 up_Isil=',up_Isil
   !endif

   !**!  Grazing term
   !!IA_gr = max(0.0_rk, (self%rGr * Ialg))
   if (icethickness <= 0.005_rk) then
        Fs = 0.0_rk
        Fl = 0.0_rk
        WZsonIA = 0.0_rk
        WZlonIA = 0.0_rk
   else
        !Fs = self%prefZsPs*Wfla + self%prefZsPl*(Wdia+Ialg) + self%prefZsD*Wdet2 + self%prefZsD*Wdet1
        !Fl = self%prefZlPs*Wfla + self%prefZlPl*(Wdia+Ialg) + self%prefZlZs*Wmicro + self%prefZlD*Wdet2 + self%prefZlD*Wdet1
        Fs = self%prefZsPs*Wfla + self%prefZsPl*(Wdia+Ialg/cell_thick) + self%prefZsD*Wdet + self%prefZsD*Wdetf
        Fl = self%prefZlPs*Wfla + self%prefZlPl*(Wdia+Ialg/cell_thick) + self%prefZlZs*Wmicro + self%prefZlD*Wdet +self%prefZlD*Wdetf

        WZsonIA =             self%GrZsP * self%prefZsPl * Ialg/ cell_thick /(self%Rg  + Fs)
        WZlonIA =             self%GrZlP * self%prefZlPl * Ialg/ cell_thick /(self%Rg  + Fl)
   endif

   !**!  Nutrients change 
   pno3 = self%vn / (self%vn + Inh4)


   rhs_Init = molN2gC * self%fr * (- Prod * pno3) + nit 
   rhs_Iamm = molN2gC * self%fr * (- Prod * (1-pno3) + rem) - nit
   rhs_Ipho = molP2gC * self%fr * (- Prod + rem)
   rhs_Isil = molS2gC * self%fr * (- Prod + rem)

   rhs_IA = Prod - IA_loss !- IA_gr
   rhs_det = IA_loss
 
   if((Ialg + (rhs_IA*baclin)) < 0.0_rk)then
        print*, 'IA < 0', Ialg
        rhs_IA = 0.0_rk
        rhs_det = 0.0_rk
        rhs_Init = 0.0_rk
        rhs_Iamm = 0.0_rk
        rhs_Ipho = 0.0_rk
        rhs_Isil = 0.0_rk
        write(*,*) 'DEBUG3 IA<0'
   endif

   !!!-------------------------------!!!
   !!!    APPLY AND WRITE CHANGE     !!!
   !!!-------------------------------!!!
   !**! ODE function
   !_SET_SURFACE_ODE_(self%id_Ialg, rhs_IA - self%rGr*WZlonIA*Wmeso - self%rGr*WZsonIA*Wmicro)
   !grazing per grid area (mgC/m2 grid/s)
   !_SET_SURFACE_ODE_(self%id_Ialg, rhs_IA - self%rGr*WZlonIA*Wmeso*cell_thick - self%rGr*WZsonIA*Wmicro*cell_thick)
   !convert to grazing per ice area (mgC /m2 ice/s) for removing from Ialg 
   _SET_SURFACE_ODE_(self%id_Ialg, rhs_IA - self%rGr*WZlonIA*Wmeso*cell_thick/max(eps,iceconc) - self%rGr*WZsonIA*Wmicro*cell_thick/max(eps,iceconc))
   _SET_SURFACE_ODE_(self%id_Idet, rhs_det  - rem + entrap_det + entrap_IA*(1.0_rk-self%surv_prop))
   _SET_SURFACE_ODE_(self%id_ino3, rhs_init + diff_flux_no3 + entrap_no3)  
   _SET_SURFACE_ODE_(self%id_inh4, rhs_iamm + diff_flux_nh4 + entrap_nh4)  
   _SET_SURFACE_ODE_(self%id_ipho, rhs_ipho + diff_flux_pho + entrap_pho)  
   _SET_SURFACE_ODE_(self%id_isil, rhs_isil + diff_flux_sil + entrap_sil) 

   if(ice_dh < 0.0_rk .AND. icethickness <= 0.005_rk)then
	melt_no3 = -max(0.0_rk, (Ino3/baclin))
	melt_nh4 = -max(0.0_rk, (Inh4/baclin)) 
	melt_sil = -max(0.0_rk, (Isil/baclin)) 
	melt_pho = -max(0.0_rk, (Ipho/baclin)) 
	melt_IA  = -max(0.0_rk, (Ialg/baclin))
	melt_det = -max(0.0_rk, (Idet/baclin))
   else 
	melt_no3 = min(0.0_rk, melt_ice) * max(0.0_rk, Ino3)
	melt_nh4 = min(0.0_rk, melt_ice) * max(0.0_rk, Inh4)
	melt_sil = min(0.0_rk, melt_ice) * max(0.0_rk, Isil)
	melt_pho = min(0.0_rk, melt_ice) * max(0.0_rk, Ipho)
	melt_IA  = min(0.0_rk, melt_ice) * max(0.0_rk, Ialg)
	melt_det = min(0.0_rk, melt_ice) * max(0.0_rk, Idet)
   end if 
 
   if((Ialg + (melt_IA*baclin)) < 0.0_rk)then
        melt_IA = -Ialg/baclin
        !write(*,*) 'DEBUG4 melt_IA more than Ialg' 
        !write(*,*) 'Ialg=',Ialg 
   endif

   _SET_SURFACE_ODE_(self%id_Ialg, melt_IA)
   _SET_SURFACE_ODE_(self%id_Idet, melt_det)
   _SET_SURFACE_ODE_(self%id_ino3, melt_no3)
   _SET_SURFACE_ODE_(self%id_inh4, melt_nh4)
   _SET_SURFACE_ODE_(self%id_ipho, melt_pho)
   _SET_SURFACE_ODE_(self%id_isil, melt_sil)

   !_GET_HORIZONTAL_(self%id_Ialg, Ialg)
   !**! Calculate attenuation of ice algae
   !Iatt = self%g_ia * Ialg * Zbot
   !if (icethickness <= 0.005_rk) then
   !     Iatt = 0.0_rk
   !else
   Iatt = self%g_ia * Ialg
        !if(Ialg>8000.OR.Iatt>1.e+7) then
        !   write(*,*) 'DEBUG5 Iatt,Ialg=',Iatt,Ialg
        !endif
   !endif
   !_SET_SURFACE_ODE_(self%id_Iatt, Iatt)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Iatt, Iatt)

   !**! Exchange with water column function  
   if(self%couple_ice) then   
   ! flux should be 2D 
   _SET_SURFACE_EXCHANGE_(self%id_Wno3, ((-diff_flux_no3 -entrap_no3 -melt_no3) / molN2gC) * iceconc)
   _SET_SURFACE_EXCHANGE_(self%id_Wnh4, ((-diff_flux_nh4 -entrap_nh4 -melt_nh4) / molN2gC) * iceconc)
   _SET_SURFACE_EXCHANGE_(self%id_Wsil, ((-diff_flux_sil -entrap_sil -melt_sil) / molS2gC) * iceconc)
   _SET_SURFACE_EXCHANGE_(self%id_Wpho, ((-diff_flux_pho -entrap_pho -melt_pho) / molP2gC) * iceconc)
   _SET_SURFACE_EXCHANGE_(self%id_Wdia, (-entrap_IA -melt_IA*self%surv_prop) * iceconc)
   _SET_SURFACE_EXCHANGE_(self%id_Wdetf, (-entrap_detf -melt_det -melt_IA*(1.0_rk-self%surv_prop))* iceconc)
   _SET_SURFACE_EXCHANGE_(self%id_Wdet, (-entrap_det) * iceconc)
  ! _SET_SURFACE_EXCHANGE_(self%id_Wmicro, self%rGr*WZsonIA*self%gammaZsp*iceconc) !IA_gr * self%prefZs * iceconc)
  ! _SET_SURFACE_EXCHANGE_(self%id_Wmeso, self%rGr*WZlonIA*self%gammaZlp*iceconc) ! IA_gr * self%prefZl * iceconc)
  !the following two terms don't need to scale by iceconc because they are already per grid area,not per ice area
   _SET_SURFACE_EXCHANGE_(self%id_Wmicro, self%rGr*WZsonIA*Wmicro*self%gammaZsp*cell_thick) !IA_gr * self%prefZs * iceconc)
   _SET_SURFACE_EXCHANGE_(self%id_Wmeso, self%rGr*WZlonIA*Wmeso*self%gammaZlp*cell_thick) ! IA_gr * self%prefZl * iceconc)
   end if

   !**! Save diagnostic variables
   !__ Algae
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_primprod, Prod)
   if (self%turn_on_additional_diagnostics) then
   !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_mortality, rhs_IA)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_mortality, IA_loss)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_grazing, self%rGr*(WZsonIA*Wmicro+WZlonIA*Wmeso)*cell_thick)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_l_limit, blight)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_sal_limit, Sdep)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_temp_limit, Tdep)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_n_limit, up_Ino3)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_p_limit, up_Ipho)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_s_limit, up_Isil)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_entIA, entrap_IA*self%surv_prop)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meltIA, melt_IA)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_WZlonIA, WZlonIA)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_WZsonIA, WZsonIA)  
   
   !__ Detritus
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_remineralisation, rem)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_entIdet, entrap_det)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meltIdet, melt_det)  
   
   !__ Nutients
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_pno3, pno3)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_upno3, rhs_init)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_upnh4, rhs_iamm)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_uppho, rhs_ipho)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_upsil, rhs_isil)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_entno3, entrap_no3)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_entnh4, entrap_nh4)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_entpho, entrap_pho)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_entsil, entrap_sil)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meltno3, melt_no3)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meltnh4, melt_nh4)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meltpho, melt_pho)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_meltsil, melt_sil)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diffno3, diff_flux_no3)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diffnh4, diff_flux_nh4)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diffpho, diff_flux_pho)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_diffsil, diff_flux_sil)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_br_no3, BR_no3)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_br_nh4, BR_nh4)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_br_pho, BR_pho)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_br_sil, BR_sil)  
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Zice, Zice)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Zbot, Zbot)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_Zwater, cell_thick) 
   
   !__ Environmental
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_lightIA, bio_par) !Fsw in Deborah's version
   !_SET_HORIZONTAL_DIAGNOSTIC_(self%id_lightBOT, bio_par)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_muTS, self%muIA*sedy0*Tdep*Sdep)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_hice, icethickness)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_dhice, ice_dh)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_iconc, iceconc)

   !__ Pelagic 
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_no3W, Wno3)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_nh4W, Wnh4)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_phoW, Wpho)
   _SET_HORIZONTAL_DIAGNOSTIC_(self%id_silW, Wsil)
   end if !additional diagnostic 
   _HORIZONTAL_LOOP_END_
   end if  
   end subroutine do_surface


   end module fabm_nersc_icealgae_G


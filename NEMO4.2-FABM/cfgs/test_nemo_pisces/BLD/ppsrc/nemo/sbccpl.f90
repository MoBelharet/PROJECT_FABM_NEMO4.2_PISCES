










MODULE sbccpl
   !!======================================================================
   !!                       ***  MODULE  sbccpl  ***
   !! Surface Boundary Condition :  momentum, heat and freshwater fluxes in coupled mode
   !!======================================================================
   !! History :  2.0  ! 2007-06  (R. Redler, N. Keenlyside, W. Park) Original code split into flxmod & taumod
   !!            3.0  ! 2008-02  (G. Madec, C Talandier)  surface module
   !!            3.1  ! 2009_02  (G. Madec, S. Masson, E. Maisonave, A. Caubel) generic coupled interface
   !!            3.4  ! 2011_11  (C. Harris) more flexibility + multi-category fields
   !!            4.2  ! 2020-12  (G. Madec, E. Clementi)  wave coupling updates
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   namsbc_cpl      : coupled formulation namlist
   !!   sbc_cpl_init    : initialisation of the coupled exchanges
   !!   sbc_cpl_rcv     : receive fields from the atmosphere over the ocean (ocean only)
   !!                     receive stress from the atmosphere over the ocean (ocean-ice case)
   !!   sbc_cpl_ice_tau : receive stress from the atmosphere over ice
   !!   sbc_cpl_ice_flx : receive fluxes from the atmosphere over ice
   !!   sbc_cpl_snd     : send     fields to the atmosphere
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! Surface boundary condition: ocean fields
   USE trc_oce         ! share SMS/Ocean variables
   USE sbc_ice         ! Surface boundary condition: ice fields
   USE sbcapr          ! Stochastic param. : ???
   USE sbcdcy          ! surface boundary condition: diurnal cycle
   USE sbcwave         ! surface boundary condition: waves
   USE phycst          ! physical constants
   USE isf_oce , ONLY : l_isfoasis, fwfisf_oasis ! ice shelf boundary condition
   USE cpl_oasis3     ! OASIS3 coupling
   USE geo2ocean      !
   USE oce     , ONLY : ts, uu, vv, ssh, fraqsr_1lev
   USE ocealb         !
   USE eosbn2         !
   USE sbcrnf  , ONLY : l_rnfcpl
   !
   USE in_out_manager ! I/O manager
   USE iom            ! NetCDF library
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)


   USE sbc_phy, ONLY : pp_cldf, rpref

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_cpl_init      ! routine called by sbcmod.F90
   PUBLIC   sbc_cpl_rcv       ! routine called by icestp.F90
   PUBLIC   sbc_cpl_snd       ! routine called by step.F90
   PUBLIC   sbc_cpl_ice_tau   ! routine called by icestp.F90
   PUBLIC   sbc_cpl_ice_flx   ! routine called by icestp.F90
   PUBLIC   sbc_cpl_alloc     ! routine called in sbcice_cice.F90

   INTEGER, PARAMETER ::   jpr_otx1   =  1   ! 3 atmosphere-ocean stress components on grid 1
   INTEGER, PARAMETER ::   jpr_oty1   =  2   !
   INTEGER, PARAMETER ::   jpr_otz1   =  3   !
   INTEGER, PARAMETER ::   jpr_otx2   =  4   ! 3 atmosphere-ocean stress components on grid 2
   INTEGER, PARAMETER ::   jpr_oty2   =  5   !
   INTEGER, PARAMETER ::   jpr_otz2   =  6   !
   INTEGER, PARAMETER ::   jpr_itx1   =  7   ! 3 atmosphere-ice   stress components on grid 1
   INTEGER, PARAMETER ::   jpr_ity1   =  8   !
   INTEGER, PARAMETER ::   jpr_itz1   =  9   !
   INTEGER, PARAMETER ::   jpr_itx2   = 10   ! 3 atmosphere-ice   stress components on grid 2
   INTEGER, PARAMETER ::   jpr_ity2   = 11   !
   INTEGER, PARAMETER ::   jpr_itz2   = 12   !
   INTEGER, PARAMETER ::   jpr_qsroce = 13   ! Qsr above the ocean
   INTEGER, PARAMETER ::   jpr_qsrice = 14   ! Qsr above the ice
   INTEGER, PARAMETER ::   jpr_qsrmix = 15
   INTEGER, PARAMETER ::   jpr_qnsoce = 16   ! Qns above the ocean
   INTEGER, PARAMETER ::   jpr_qnsice = 17   ! Qns above the ice
   INTEGER, PARAMETER ::   jpr_qnsmix = 18
   INTEGER, PARAMETER ::   jpr_rain   = 19   ! total liquid precipitation (rain)
   INTEGER, PARAMETER ::   jpr_snow   = 20   ! solid precipitation over the ocean (snow)
   INTEGER, PARAMETER ::   jpr_tevp   = 21   ! total evaporation
   INTEGER, PARAMETER ::   jpr_ievp   = 22   ! solid evaporation (sublimation)
   INTEGER, PARAMETER ::   jpr_sbpr   = 23   ! sublimation - liquid precipitation - solid precipitation
   INTEGER, PARAMETER ::   jpr_semp   = 24   ! solid freshwater budget (sublimation - snow)
   INTEGER, PARAMETER ::   jpr_oemp   = 25   ! ocean freshwater budget (evap - precip)
   INTEGER, PARAMETER ::   jpr_w10m   = 26   ! 10m wind
   INTEGER, PARAMETER ::   jpr_dqnsdt = 27   ! d(Q non solar)/d(temperature)
   INTEGER, PARAMETER ::   jpr_rnf    = 28   ! runoffs
   INTEGER, PARAMETER ::   jpr_cal    = 29   ! calving
   INTEGER, PARAMETER ::   jpr_taum   = 30   ! wind stress module
   INTEGER, PARAMETER ::   jpr_co2    = 31
   INTEGER, PARAMETER ::   jpr_topm   = 32   ! topmeltn
   INTEGER, PARAMETER ::   jpr_botm   = 33   ! botmeltn
   INTEGER, PARAMETER ::   jpr_sflx   = 34   ! salt flux
   INTEGER, PARAMETER ::   jpr_toce   = 35   ! ocean temperature
   INTEGER, PARAMETER ::   jpr_soce   = 36   ! ocean salinity
   INTEGER, PARAMETER ::   jpr_ocx1   = 37   ! ocean current on grid 1
   INTEGER, PARAMETER ::   jpr_ocy1   = 38   !
   INTEGER, PARAMETER ::   jpr_ssh    = 39   ! sea surface height
   INTEGER, PARAMETER ::   jpr_fice   = 40   ! ice fraction
   INTEGER, PARAMETER ::   jpr_e3t1st = 41   ! first T level thickness
   INTEGER, PARAMETER ::   jpr_fraqsr = 42   ! fraction of solar net radiation absorbed in the first ocean level
   INTEGER, PARAMETER ::   jpr_mslp   = 43   ! mean sea level pressure
   !**  surface wave coupling  **
   INTEGER, PARAMETER ::   jpr_hsig   = 44   ! Hsig
   INTEGER, PARAMETER ::   jpr_phioc  = 45   ! Wave=>ocean energy flux
   INTEGER, PARAMETER ::   jpr_sdrftx = 46   ! Stokes drift on grid 1
   INTEGER, PARAMETER ::   jpr_sdrfty = 47   ! Stokes drift on grid 2
   INTEGER, PARAMETER ::   jpr_wper   = 48   ! Mean wave period
   INTEGER, PARAMETER ::   jpr_wnum   = 49   ! Mean wavenumber
   INTEGER, PARAMETER ::   jpr_wstrf  = 50   ! Stress fraction adsorbed by waves
   INTEGER, PARAMETER ::   jpr_wdrag  = 51   ! Neutral surface drag coefficient
   INTEGER, PARAMETER ::   jpr_charn  = 52   ! Chranock coefficient
   INTEGER, PARAMETER ::   jpr_twox   = 53   ! wave to ocean momentum flux
   INTEGER, PARAMETER ::   jpr_twoy   = 54   ! wave to ocean momentum flux
   INTEGER, PARAMETER ::   jpr_tawx   = 55   ! net wave-supported stress
   INTEGER, PARAMETER ::   jpr_tawy   = 56   ! net wave-supported stress
   INTEGER, PARAMETER ::   jpr_bhd    = 57   ! Bernoulli head. waves' induced surface pressure
   INTEGER, PARAMETER ::   jpr_tusd   = 58   ! zonal stokes transport
   INTEGER, PARAMETER ::   jpr_tvsd   = 59   ! meridional stokes tranmport
   INTEGER, PARAMETER ::   jpr_isf    = 60
   INTEGER, PARAMETER ::   jpr_icb    = 61
   INTEGER, PARAMETER ::   jpr_ts_ice = 62   ! Sea ice surface temp
   !!INTEGER, PARAMETER ::   jpr_qtrice = 63   ! Transmitted solar thru sea-ice

   INTEGER, PARAMETER ::   jprcv      = 62   ! total number of fields received

   INTEGER, PARAMETER ::   jps_fice   =  1   ! ice fraction sent to the atmosphere
   INTEGER, PARAMETER ::   jps_toce   =  2   ! ocean temperature
   INTEGER, PARAMETER ::   jps_tice   =  3   ! ice   temperature
   INTEGER, PARAMETER ::   jps_tmix   =  4   ! mixed temperature (ocean+ice)
   INTEGER, PARAMETER ::   jps_albice =  5   ! ice   albedo
   INTEGER, PARAMETER ::   jps_albmix =  6   ! mixed albedo
   INTEGER, PARAMETER ::   jps_hice   =  7   ! ice  thickness
   INTEGER, PARAMETER ::   jps_hsnw   =  8   ! snow thickness
   INTEGER, PARAMETER ::   jps_ocx1   =  9   ! ocean current on grid 1
   INTEGER, PARAMETER ::   jps_ocy1   = 10   !
   INTEGER, PARAMETER ::   jps_ocz1   = 11   !
   INTEGER, PARAMETER ::   jps_ivx1   = 12   ! ice   current on grid 1
   INTEGER, PARAMETER ::   jps_ivy1   = 13   !
   INTEGER, PARAMETER ::   jps_ivz1   = 14   !
   INTEGER, PARAMETER ::   jps_co2    = 15
   INTEGER, PARAMETER ::   jps_soce   = 16   ! ocean salinity
   INTEGER, PARAMETER ::   jps_ssh    = 17   ! sea surface height
   INTEGER, PARAMETER ::   jps_qsroce = 18   ! Qsr above the ocean
   INTEGER, PARAMETER ::   jps_qnsoce = 19   ! Qns above the ocean
   INTEGER, PARAMETER ::   jps_oemp   = 20   ! ocean freshwater budget (evap - precip)
   INTEGER, PARAMETER ::   jps_sflx   = 21   ! salt flux
   INTEGER, PARAMETER ::   jps_otx1   = 22   ! 2 atmosphere-ocean stress components on grid 1
   INTEGER, PARAMETER ::   jps_oty1   = 23   !
   INTEGER, PARAMETER ::   jps_rnf    = 24   ! runoffs
   INTEGER, PARAMETER ::   jps_taum   = 25   ! wind stress module
   INTEGER, PARAMETER ::   jps_fice2  = 26   ! ice fraction sent to OCE (by SAS when doing SAS-OCE coupling)
   INTEGER, PARAMETER ::   jps_e3t1st = 27   ! first level depth (vvl)
   INTEGER, PARAMETER ::   jps_fraqsr = 28   ! fraction of solar net radiation absorbed in the first ocean level
   INTEGER, PARAMETER ::   jps_ficet  = 29   ! total ice fraction
   INTEGER, PARAMETER ::   jps_ocxw   = 30   ! currents on grid 1
   INTEGER, PARAMETER ::   jps_ocyw   = 31   ! currents on grid 2
   INTEGER, PARAMETER ::   jps_wlev   = 32   ! water level
   INTEGER, PARAMETER ::   jps_fice1  = 33   ! first-order ice concentration (for semi-implicit coupling of atmos-ice fluxes)
   INTEGER, PARAMETER ::   jps_a_p    = 34   ! meltpond area fraction
   INTEGER, PARAMETER ::   jps_ht_p   = 35   ! meltpond thickness
   INTEGER, PARAMETER ::   jps_kice   = 36   ! sea ice effective conductivity
   INTEGER, PARAMETER ::   jps_sstfrz = 37   ! sea surface freezing temperature
   INTEGER, PARAMETER ::   jps_ttilyr = 38   ! sea ice top layer temp

   INTEGER, PARAMETER ::   jpsnd      = 38   ! total number of fields sent

   ! Dummy variables to enable compilation when oasis3 is not being used
   INTEGER                    ::   OASIS_Sent        = -1
   INTEGER                    ::   OASIS_SentOut     = -1
   INTEGER                    ::   OASIS_ToRest      = -1
   INTEGER                    ::   OASIS_ToRestOut   = -1

   !                                  !!** namelist namsbc_cpl **
   TYPE ::   FLD_C                     !
      CHARACTER(len = 32) ::   cldes      ! desciption of the coupling strategy
      CHARACTER(len = 32) ::   clcat      ! multiple ice categories strategy
      CHARACTER(len = 32) ::   clvref     ! reference of vector ('spherical' or 'cartesian')
      CHARACTER(len = 32) ::   clvor      ! orientation of vector fields ('eastward-northward' or 'local grid')
      CHARACTER(len = 32) ::   clvgrd     ! grids on which is located the vector fields
   END TYPE FLD_C
   !                                   ! Send to the atmosphere
   TYPE(FLD_C) ::   sn_snd_temp  , sn_snd_alb , sn_snd_thick, sn_snd_crt   , sn_snd_co2,  &
      &             sn_snd_thick1, sn_snd_cond, sn_snd_mpnd , sn_snd_sstfrz, sn_snd_ttilyr
   !                                   ! Received from the atmosphere
   TYPE(FLD_C) ::   sn_rcv_w10m, sn_rcv_taumod, sn_rcv_tau, sn_rcv_dqnsdt, sn_rcv_qsr,  &
      &             sn_rcv_qns , sn_rcv_emp   , sn_rcv_rnf, sn_rcv_ts_ice
   TYPE(FLD_C) ::   sn_rcv_cal, sn_rcv_iceflx, sn_rcv_co2, sn_rcv_mslp, sn_rcv_icb, sn_rcv_isf
   !                                   ! Send to waves
   TYPE(FLD_C) ::   sn_snd_ifrac, sn_snd_crtw, sn_snd_wlev
   !                                   ! Received from waves
   TYPE(FLD_C) ::   sn_rcv_hsig, sn_rcv_phioc, sn_rcv_sdrfx, sn_rcv_sdrfy, sn_rcv_wper, sn_rcv_wnum, &
      &             sn_rcv_wstrf, sn_rcv_wdrag, sn_rcv_charn, sn_rcv_taw, sn_rcv_bhd, sn_rcv_tusd, sn_rcv_tvsd
   !                                   ! Other namelist parameters
!!   TYPE(FLD_C) ::   sn_rcv_qtrice
   INTEGER     ::   nn_cplmodel           ! Maximum number of models to/from which NEMO is potentialy sending/receiving data
   LOGICAL     ::   ln_usecplmask         !  use a coupling mask file to merge data received from several models
                                          !   -> file cplmask.nc with the float variable called cplmask (jpi,jpj,nn_cplmodel)
   LOGICAL     ::   ln_scale_ice_flux     !  use ice fluxes that are already "ice weighted" ( i.e. multiplied ice concentration)

   TYPE ::   DYNARR
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   z3
   END TYPE DYNARR

   TYPE( DYNARR ), SAVE, DIMENSION(jprcv) ::   frcv                ! all fields recieved from the atmosphere

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   alb_oce_mix    ! ocean albedo sent to atmosphere (mix clear/overcast sky)

   INTEGER , ALLOCATABLE, SAVE, DIMENSION(:) ::   nrcvinfo           ! OASIS info argument

   !! Substitution





!!----------------------------------------------------------------------
!!                    ***  domzgr_substitute.h90   ***
!!----------------------------------------------------------------------
!! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
!!      factors depending on the vertical coord. used, using CPP macro.
!!----------------------------------------------------------------------
!! History :  4.2  !  2020-02  (S. Techene, G. Madec)  star coordinate
!!----------------------------------------------------------------------
!! NEMO/OCE 4.2 , NEMO Consortium (2020)
!! $Id$
!! Software governed by the CeCILL license (see ./LICENSE)
!!----------------------------------------------------------------------
!!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbccpl.F90 15551 2021-11-28 20:19:36Z gsamson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_cpl_alloc()
      !!----------------------------------------------------------------------
      !!             ***  FUNCTION sbc_cpl_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(4)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( alb_oce_mix(jpi,jpj), nrcvinfo(jprcv),  STAT=ierr(1) )

      ALLOCATE( a_i(jpi,jpj,1) , STAT=ierr(2) )  ! used in sbcice_if.F90 (done here as there is no sbc_ice_if_init)
      ALLOCATE( xcplmask(jpi,jpj,0:nn_cplmodel) , STAT=ierr(3) )
      !
      IF( .NOT. ln_apr_dyn ) ALLOCATE( ssh_ib(jpi,jpj), ssh_ibb(jpi,jpj), apr(jpi, jpj), STAT=ierr(4) )

      sbc_cpl_alloc = MAXVAL( ierr )
      CALL mpp_sum ( 'sbccpl', sbc_cpl_alloc )
      IF( sbc_cpl_alloc > 0 )   CALL ctl_warn('sbc_cpl_alloc: allocation of arrays failed')
      !
   END FUNCTION sbc_cpl_alloc


   SUBROUTINE sbc_cpl_init( k_ice )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_init  ***
      !!
      !! ** Purpose :   Initialisation of send and received information from
      !!                the atmospheric component
      !!
      !! ** Method  : * Read namsbc_cpl namelist
      !!              * define the receive interface
      !!              * define the send    interface
      !!              * initialise the OASIS coupler
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   k_ice   ! ice management in the sbc (=0/1/2/3)
      !
      INTEGER ::   jn          ! dummy loop index
      INTEGER ::   ios, inum   ! Local integer
      REAL(wp), DIMENSION(jpi,jpj) ::   zacs, zaos
      !!
      NAMELIST/namsbc_cpl/  nn_cplmodel  , ln_usecplmask, nn_cats_cpl , ln_scale_ice_flux,             &
         &                  sn_snd_temp  , sn_snd_alb   , sn_snd_thick, sn_snd_crt   , sn_snd_co2   ,  &
         &                  sn_snd_ttilyr, sn_snd_cond  , sn_snd_mpnd , sn_snd_sstfrz, sn_snd_thick1,  &
         &                  sn_snd_ifrac , sn_snd_crtw  , sn_snd_wlev , sn_rcv_hsig  , sn_rcv_phioc ,  &
         &                  sn_rcv_w10m  , sn_rcv_taumod, sn_rcv_tau  , sn_rcv_dqnsdt, sn_rcv_qsr   ,  &
         &                  sn_rcv_sdrfx , sn_rcv_sdrfy , sn_rcv_wper , sn_rcv_wnum  , sn_rcv_wstrf ,  &
         &                  sn_rcv_charn , sn_rcv_taw   , sn_rcv_bhd  , sn_rcv_tusd  , sn_rcv_tvsd,    &
         &                  sn_rcv_wdrag , sn_rcv_qns   , sn_rcv_emp  , sn_rcv_rnf   , sn_rcv_cal  ,   &
         &                  sn_rcv_iceflx, sn_rcv_co2   , sn_rcv_icb  , sn_rcv_isf   , sn_rcv_ts_ice,  & !!, sn_rcv_qtrice
         &                  sn_rcv_mslp

      !!---------------------------------------------------------------------
      !
      ! ================================ !
      !      Namelist informations       !
      ! ================================ !
      !
      READ  ( numnam_ref, namsbc_cpl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_cpl in reference namelist' )
      !
      READ  ( numnam_cfg, namsbc_cpl, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_cpl in configuration namelist' )
      IF(lwm) WRITE ( numond, namsbc_cpl )
      !
      IF(lwp) THEN                        ! control print
         WRITE(numout,*)
         WRITE(numout,*)'sbc_cpl_init : namsbc_cpl namelist '
         WRITE(numout,*)'~~~~~~~~~~~~'
      ENDIF
      IF( lwp .AND. ln_cpl ) THEN                        ! control print
         WRITE(numout,*)'  nn_cplmodel                         = ', nn_cplmodel
         WRITE(numout,*)'  ln_usecplmask                       = ', ln_usecplmask
         WRITE(numout,*)'  ln_scale_ice_flux                   = ', ln_scale_ice_flux
         WRITE(numout,*)'  nn_cats_cpl                         = ', nn_cats_cpl
         WRITE(numout,*)'  received fields (mutiple ice categogies)'
         WRITE(numout,*)'      10m wind module                 = ', TRIM(sn_rcv_w10m%cldes  ), ' (', TRIM(sn_rcv_w10m%clcat  ), ')'
         WRITE(numout,*)'      stress module                   = ', TRIM(sn_rcv_taumod%cldes), ' (', TRIM(sn_rcv_taumod%clcat), ')'
         WRITE(numout,*)'      surface stress                  = ', TRIM(sn_rcv_tau%cldes   ), ' (', TRIM(sn_rcv_tau%clcat   ), ')'
         WRITE(numout,*)'                     - referential    = ', sn_rcv_tau%clvref
         WRITE(numout,*)'                     - orientation    = ', sn_rcv_tau%clvor
         WRITE(numout,*)'                     - mesh           = ', sn_rcv_tau%clvgrd
         WRITE(numout,*)'      non-solar heat flux sensitivity = ', TRIM(sn_rcv_dqnsdt%cldes), ' (', TRIM(sn_rcv_dqnsdt%clcat), ')'
         WRITE(numout,*)'      solar heat flux                 = ', TRIM(sn_rcv_qsr%cldes   ), ' (', TRIM(sn_rcv_qsr%clcat   ), ')'
         WRITE(numout,*)'      non-solar heat flux             = ', TRIM(sn_rcv_qns%cldes   ), ' (', TRIM(sn_rcv_qns%clcat   ), ')'
         WRITE(numout,*)'      freshwater budget               = ', TRIM(sn_rcv_emp%cldes   ), ' (', TRIM(sn_rcv_emp%clcat   ), ')'
         WRITE(numout,*)'      runoffs                         = ', TRIM(sn_rcv_rnf%cldes   ), ' (', TRIM(sn_rcv_rnf%clcat   ), ')'
         WRITE(numout,*)'      calving                         = ', TRIM(sn_rcv_cal%cldes   ), ' (', TRIM(sn_rcv_cal%clcat   ), ')'
         WRITE(numout,*)'      iceberg                         = ', TRIM(sn_rcv_icb%cldes   ), ' (', TRIM(sn_rcv_icb%clcat   ), ')'
         WRITE(numout,*)'      ice shelf                       = ', TRIM(sn_rcv_isf%cldes   ), ' (', TRIM(sn_rcv_isf%clcat   ), ')'
         WRITE(numout,*)'      sea ice heat fluxes             = ', TRIM(sn_rcv_iceflx%cldes), ' (', TRIM(sn_rcv_iceflx%clcat), ')'
!!       WRITE(numout,*)'      transmitted solar thru sea-ice  = ', TRIM(sn_rcv_qtrice%cldes), ' (', TRIM(sn_rcv_qtrice%clcat), ')'
         WRITE(numout,*)'      atm co2                         = ', TRIM(sn_rcv_co2%cldes   ), ' (', TRIM(sn_rcv_co2%clcat   ), ')'
         WRITE(numout,*)'      Sea ice surface skin temperature= ', TRIM(sn_rcv_ts_ice%cldes), ' (', TRIM(sn_rcv_ts_ice%clcat), ')'
         WRITE(numout,*)'      surface waves:'
         WRITE(numout,*)'      significant wave heigth         = ', TRIM(sn_rcv_hsig%cldes  ), ' (', TRIM(sn_rcv_hsig%clcat  ), ')'
         WRITE(numout,*)'      wave to oce energy flux         = ', TRIM(sn_rcv_phioc%cldes ), ' (', TRIM(sn_rcv_phioc%clcat ), ')'
         WRITE(numout,*)'      Surface Stokes drift grid u     = ', TRIM(sn_rcv_sdrfx%cldes ), ' (', TRIM(sn_rcv_sdrfx%clcat ), ')'
         WRITE(numout,*)'      Surface Stokes drift grid v     = ', TRIM(sn_rcv_sdrfy%cldes ), ' (', TRIM(sn_rcv_sdrfy%clcat ), ')'
         WRITE(numout,*)'      Mean wave period                = ', TRIM(sn_rcv_wper%cldes  ), ' (', TRIM(sn_rcv_wper%clcat  ), ')'
         WRITE(numout,*)'      Mean wave number                = ', TRIM(sn_rcv_wnum%cldes  ), ' (', TRIM(sn_rcv_wnum%clcat  ), ')'
         WRITE(numout,*)'      Stress frac adsorbed by waves   = ', TRIM(sn_rcv_wstrf%cldes ), ' (', TRIM(sn_rcv_wstrf%clcat ), ')'
         WRITE(numout,*)'      Neutral surf drag coefficient   = ', TRIM(sn_rcv_wdrag%cldes ), ' (', TRIM(sn_rcv_wdrag%clcat ), ')'
         WRITE(numout,*)'      Charnock coefficient            = ', TRIM(sn_rcv_charn%cldes ), ' (', TRIM(sn_rcv_charn%clcat ), ')'
         WRITE(numout,*)'  sent fields (multiple ice categories)'
         WRITE(numout,*)'      surface temperature             = ', TRIM(sn_snd_temp%cldes  ), ' (', TRIM(sn_snd_temp%clcat  ), ')'
         WRITE(numout,*)'      top ice layer temperature       = ', TRIM(sn_snd_ttilyr%cldes), ' (', TRIM(sn_snd_ttilyr%clcat), ')'
         WRITE(numout,*)'      albedo                          = ', TRIM(sn_snd_alb%cldes   ), ' (', TRIM(sn_snd_alb%clcat   ), ')'
         WRITE(numout,*)'      ice/snow thickness              = ', TRIM(sn_snd_thick%cldes ), ' (', TRIM(sn_snd_thick%clcat ), ')'
         WRITE(numout,*)'      total ice fraction              = ', TRIM(sn_snd_ifrac%cldes ), ' (', TRIM(sn_snd_ifrac%clcat ), ')'
         WRITE(numout,*)'      surface current                 = ', TRIM(sn_snd_crt%cldes   ), ' (', TRIM(sn_snd_crt%clcat   ), ')'
         WRITE(numout,*)'                      - referential   = ', sn_snd_crt%clvref
         WRITE(numout,*)'                      - orientation   = ', sn_snd_crt%clvor
         WRITE(numout,*)'                      - mesh          = ', sn_snd_crt%clvgrd
         WRITE(numout,*)'      oce co2 flux                    = ', TRIM(sn_snd_co2%cldes   ), ' (', TRIM(sn_snd_co2%clcat   ), ')'
         WRITE(numout,*)'      ice effective conductivity      = ', TRIM(sn_snd_cond%cldes  ), ' (', TRIM(sn_snd_cond%clcat  ), ')'
         WRITE(numout,*)'      meltponds fraction and depth    = ', TRIM(sn_snd_mpnd%cldes  ), ' (', TRIM(sn_snd_mpnd%clcat  ), ')'
         WRITE(numout,*)'      sea surface freezing temp       = ', TRIM(sn_snd_sstfrz%cldes), ' (', TRIM(sn_snd_sstfrz%clcat), ')'
         WRITE(numout,*)'      water level                     = ', TRIM(sn_snd_wlev%cldes  ), ' (', TRIM(sn_snd_wlev%clcat  ), ')'
         WRITE(numout,*)'      mean sea level pressure         = ', TRIM(sn_rcv_mslp%cldes  ), ' (', TRIM(sn_rcv_mslp%clcat  ), ')'
         WRITE(numout,*)'      surface current to waves        = ', TRIM(sn_snd_crtw%cldes  ), ' (', TRIM(sn_snd_crtw%clcat  ), ')'
         WRITE(numout,*)'                      - referential   = ', sn_snd_crtw%clvref
         WRITE(numout,*)'                      - orientation   = ', sn_snd_crtw%clvor
         WRITE(numout,*)'                      - mesh          = ', sn_snd_crtw%clvgrd
      ENDIF
      IF( lwp .AND. ln_wave) THEN                        ! control print
      WRITE(numout,*)'      surface waves:'
         WRITE(numout,*)'      Significant wave heigth         = ', TRIM(sn_rcv_hsig%cldes  ), ' (', TRIM(sn_rcv_hsig%clcat  ), ')'
         WRITE(numout,*)'      Wave to oce energy flux         = ', TRIM(sn_rcv_phioc%cldes ), ' (', TRIM(sn_rcv_phioc%clcat ), ')'
         WRITE(numout,*)'      Surface Stokes drift grid u     = ', TRIM(sn_rcv_sdrfx%cldes ), ' (', TRIM(sn_rcv_sdrfx%clcat ), ')'
         WRITE(numout,*)'      Surface Stokes drift grid v     = ', TRIM(sn_rcv_sdrfy%cldes ), ' (', TRIM(sn_rcv_sdrfy%clcat ), ')'
         WRITE(numout,*)'      Mean wave period                = ', TRIM(sn_rcv_wper%cldes  ), ' (', TRIM(sn_rcv_wper%clcat  ), ')'
         WRITE(numout,*)'      Mean wave number                = ', TRIM(sn_rcv_wnum%cldes  ), ' (', TRIM(sn_rcv_wnum%clcat  ), ')'
         WRITE(numout,*)'      Stress frac adsorbed by waves   = ', TRIM(sn_rcv_wstrf%cldes ), ' (', TRIM(sn_rcv_wstrf%clcat ), ')'
         WRITE(numout,*)'      Neutral surf drag coefficient   = ', TRIM(sn_rcv_wdrag%cldes ), ' (', TRIM(sn_rcv_wdrag%clcat ), ')'
         WRITE(numout,*)'      Charnock coefficient            = ', TRIM(sn_rcv_charn%cldes ), ' (', TRIM(sn_rcv_charn%clcat ), ')'
         WRITE(numout,*)' Transport associated to Stokes drift grid u = ', TRIM(sn_rcv_tusd%cldes ), ' (', TRIM(sn_rcv_tusd%clcat ), ')'
         WRITE(numout,*)' Transport associated to Stokes drift grid v = ', TRIM(sn_rcv_tvsd%cldes ), ' (', TRIM(sn_rcv_tvsd%clcat ), ')'
         WRITE(numout,*)'      Bernouilli pressure head        = ', TRIM(sn_rcv_bhd%cldes   ), ' (', TRIM(sn_rcv_bhd%clcat  ), ')'
         WRITE(numout,*)'Wave to ocean momentum flux and Net wave-supported stress = ', TRIM(sn_rcv_taw%cldes ), ' (', TRIM(sn_rcv_taw%clcat ), ')'
         WRITE(numout,*)'      Surface current to waves        = ', TRIM(sn_snd_crtw%cldes  ), ' (', TRIM(sn_snd_crtw%clcat  ), ')'
         WRITE(numout,*)'                      - referential   = ', sn_snd_crtw%clvref
         WRITE(numout,*)'                      - orientation   = ', sn_snd_crtw%clvor
         WRITE(numout,*)'                      - mesh          = ', sn_snd_crtw%clvgrd
      ENDIF
      !                                   ! allocate sbccpl arrays
      IF( sbc_cpl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_cpl_alloc : unable to allocate arrays' )

      ! ================================ !
      !   Define the receive interface   !
      ! ================================ !
      nrcvinfo(:) = OASIS_idle   ! needed by nrcvinfo(jpr_otx1) if we do not receive ocean stress

      ! for each field: define the OASIS name                              (srcv(:)%clname)
      !                 define receive or not from the namelist parameters (srcv(:)%laction)
      !                 define the north fold type of lbc                  (srcv(:)%nsgn)

      ! default definitions of srcv
      srcv(:)%laction = .FALSE.   ;   srcv(:)%clgrid = 'T'   ;   srcv(:)%nsgn = 1.   ;   srcv(:)%nct = 1

      !                                                      ! ------------------------- !
      !                                                      ! ice and ocean wind stress !
      !                                                      ! ------------------------- !
      !                                                           ! Name
      srcv(jpr_otx1)%clname = 'O_OTaux1'      ! 1st ocean component on grid ONE (T or U)
      srcv(jpr_oty1)%clname = 'O_OTauy1'      ! 2nd   -      -         -     -
      srcv(jpr_otz1)%clname = 'O_OTauz1'      ! 3rd   -      -         -     -
      srcv(jpr_otx2)%clname = 'O_OTaux2'      ! 1st ocean component on grid TWO (V)
      srcv(jpr_oty2)%clname = 'O_OTauy2'      ! 2nd   -      -         -     -
      srcv(jpr_otz2)%clname = 'O_OTauz2'      ! 3rd   -      -         -     -
      !
      srcv(jpr_itx1)%clname = 'O_ITaux1'      ! 1st  ice  component on grid ONE (T, F, I or U)
      srcv(jpr_ity1)%clname = 'O_ITauy1'      ! 2nd   -      -         -     -
      srcv(jpr_itz1)%clname = 'O_ITauz1'      ! 3rd   -      -         -     -
      srcv(jpr_itx2)%clname = 'O_ITaux2'      ! 1st  ice  component on grid TWO (V)
      srcv(jpr_ity2)%clname = 'O_ITauy2'      ! 2nd   -      -         -     -
      srcv(jpr_itz2)%clname = 'O_ITauz2'      ! 3rd   -      -         -     -
      !
      ! Vectors: change of sign at north fold ONLY if on the local grid
      IF(       TRIM( sn_rcv_tau%cldes ) == 'oce only' .OR. TRIM( sn_rcv_tau%cldes ) == 'oce and ice'  &
           .OR. TRIM( sn_rcv_tau%cldes ) == 'mixed oce-ice' ) THEN ! avoid working with the atmospheric fields if they are not coupled
      !
      IF( TRIM( sn_rcv_tau%clvor ) == 'local grid' )   srcv(jpr_otx1:jpr_itz2)%nsgn = -1.

      !                                                           ! Set grid and action
      SELECT CASE( TRIM( sn_rcv_tau%clvgrd ) )      !  'T', 'U,V', 'U,V,I', 'U,V,F', 'T,I', 'T,F', or 'T,U,V'
      CASE( 'T' )
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1
      CASE( 'U,V' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'U'        ! ice components given at U-point
         srcv(jpr_itx2:jpr_itz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_otx1:jpr_itz2)%laction = .TRUE.     ! receive oce and ice components on both grid 1 & 2
      CASE( 'U,V,T' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'T'        ! ice components given at T-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'U,V,I' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'I'        ! ice components given at I-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'U,V,F' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'U'        ! oce components given at U-point
         srcv(jpr_otx2:jpr_otz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'F'        ! ice components given at F-point
         srcv(jpr_otx1:jpr_otz2)%laction = .TRUE.     ! receive oce components on grid 1 & 2
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1 only
      CASE( 'T,I' )
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'I'        ! ice components given at I-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1
      CASE( 'T,F' )
         srcv(jpr_otx1:jpr_itz2)%clgrid  = 'T'        ! oce and ice components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'F'        ! ice components given at F-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1
         srcv(jpr_itx1:jpr_itz1)%laction = .TRUE.     ! receive ice components on grid 1
      CASE( 'T,U,V' )
         srcv(jpr_otx1:jpr_otz1)%clgrid  = 'T'        ! oce components given at T-point
         srcv(jpr_itx1:jpr_itz1)%clgrid  = 'U'        ! ice components given at U-point
         srcv(jpr_itx2:jpr_itz2)%clgrid  = 'V'        !           and           V-point
         srcv(jpr_otx1:jpr_otz1)%laction = .TRUE.     ! receive oce components on grid 1 only
         srcv(jpr_itx1:jpr_itz2)%laction = .TRUE.     ! receive ice components on grid 1 & 2
      CASE default
         CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_tau%clvgrd' )
      END SELECT
      !
      IF( TRIM( sn_rcv_tau%clvref ) == 'spherical' )   &           ! spherical: 3rd component not received
         &     srcv( (/jpr_otz1, jpr_otz2, jpr_itz1, jpr_itz2/) )%laction = .FALSE.
      !
      IF( TRIM( sn_rcv_tau%clvor  ) == 'local grid' ) THEN        ! already on local grid -> no need of the second grid
            srcv(jpr_otx2:jpr_otz2)%laction = .FALSE.
            srcv(jpr_itx2:jpr_itz2)%laction = .FALSE.
            srcv(jpr_oty1)%clgrid = srcv(jpr_oty2)%clgrid   ! not needed but cleaner...
            srcv(jpr_ity1)%clgrid = srcv(jpr_ity2)%clgrid   ! not needed but cleaner...
      ENDIF
      !
      IF( TRIM( sn_rcv_tau%cldes ) /= 'oce and ice' ) THEN        ! 'oce and ice' case ocean stress on ocean mesh used
         srcv(jpr_itx1:jpr_itz2)%laction = .FALSE.    ! ice components not received
         srcv(jpr_itx1)%clgrid = 'U'                  ! ocean stress used after its transformation
         srcv(jpr_ity1)%clgrid = 'V'                  ! i.e. it is always at U- & V-points for i- & j-comp. resp.
      ENDIF
      ENDIF

      !                                                      ! ------------------------- !
      !                                                      !    freshwater budget      !   E-P
      !                                                      ! ------------------------- !
      ! we suppose that atmosphere modele do not make the difference between precipiration (liquide or solid)
      ! over ice of free ocean within the same atmospheric cell.cd
      srcv(jpr_rain)%clname = 'OTotRain'      ! Rain = liquid precipitation
      srcv(jpr_snow)%clname = 'OTotSnow'      ! Snow = solid precipitation
      srcv(jpr_tevp)%clname = 'OTotEvap'      ! total evaporation (over oce + ice sublimation)
      srcv(jpr_ievp)%clname = 'OIceEvap'      ! evaporation over ice = sublimation
      srcv(jpr_sbpr)%clname = 'OSubMPre'      ! sublimation - liquid precipitation - solid precipitation
      srcv(jpr_semp)%clname = 'OISubMSn'      ! ice solid water budget = sublimation - solid precipitation
      srcv(jpr_oemp)%clname = 'OOEvaMPr'      ! ocean water budget = ocean Evap - ocean precip
      SELECT CASE( TRIM( sn_rcv_emp%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(jpr_oemp)%laction = .TRUE.
      CASE( 'conservative'  )
         srcv( (/jpr_rain, jpr_snow, jpr_ievp, jpr_tevp/) )%laction = .TRUE.
         IF( k_ice <= 1 )  srcv(jpr_ievp)%laction = .FALSE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_ievp, jpr_sbpr, jpr_semp, jpr_oemp/) )%laction = .TRUE.
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_emp%cldes' )
      END SELECT
      !
      !                                                      ! ------------------------- !
      !                                                      !     Runoffs & Calving     !
      !                                                      ! ------------------------- !
      srcv(jpr_rnf   )%clname = 'O_Runoff'
      IF( TRIM( sn_rcv_rnf%cldes ) == 'coupled' ) THEN
         srcv(jpr_rnf)%laction = .TRUE.
         l_rnfcpl              = .TRUE.                      ! -> no need to read runoffs in sbcrnf
         ln_rnf                = nn_components /= jp_iam_sas ! -> force to go through sbcrnf if not sas
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   runoffs received from oasis -> force ln_rnf = ', ln_rnf
      ENDIF
      !
      srcv(jpr_cal)%clname = 'OCalving'   ;  IF( TRIM( sn_rcv_cal%cldes) == 'coupled' )   srcv(jpr_cal)%laction = .TRUE.
      srcv(jpr_isf)%clname = 'OIcshelf'   ;  IF( TRIM( sn_rcv_isf%cldes) == 'coupled' )   srcv(jpr_isf)%laction = .TRUE.
      srcv(jpr_icb)%clname = 'OIceberg'   ;  IF( TRIM( sn_rcv_icb%cldes) == 'coupled' )   srcv(jpr_icb)%laction = .TRUE.

      IF( srcv(jpr_isf)%laction ) THEN
         l_isfoasis = .TRUE.  ! -> isf fwf comes from oasis
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   iceshelf received from oasis '
      ENDIF
      !
      !
      !                                                      ! ------------------------- !
      !                                                      !    non solar radiation    !   Qns
      !                                                      ! ------------------------- !
      srcv(jpr_qnsoce)%clname = 'O_QnsOce'
      srcv(jpr_qnsice)%clname = 'O_QnsIce'
      srcv(jpr_qnsmix)%clname = 'O_QnsMix'
      SELECT CASE( TRIM( sn_rcv_qns%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(               jpr_qnsoce   )%laction = .TRUE.
      CASE( 'conservative'  )   ;   srcv( (/jpr_qnsice, jpr_qnsmix/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_qnsice, jpr_qnsoce/) )%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   srcv(               jpr_qnsmix   )%laction = .TRUE.
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_qns%cldes' )
      END SELECT
      IF( TRIM( sn_rcv_qns%cldes ) == 'mixed oce-ice' .AND. nn_cats_cpl > 1 ) &
         CALL ctl_stop( 'sbc_cpl_init: sn_rcv_qns%cldes not currently allowed to be mixed oce-ice for multi-category ice' )
      !
      !                                                      ! ------------------------- !
      !                                                      !    solar radiation        !   Qsr
      !                                                      ! ------------------------- !
      srcv(jpr_qsroce)%clname = 'O_QsrOce'
      srcv(jpr_qsrice)%clname = 'O_QsrIce'
      srcv(jpr_qsrmix)%clname = 'O_QsrMix'
      SELECT CASE( TRIM( sn_rcv_qsr%cldes ) )
      CASE( 'none'          )       ! nothing to do
      CASE( 'oce only'      )   ;   srcv(               jpr_qsroce   )%laction = .TRUE.
      CASE( 'conservative'  )   ;   srcv( (/jpr_qsrice, jpr_qsrmix/) )%laction = .TRUE.
      CASE( 'oce and ice'   )   ;   srcv( (/jpr_qsrice, jpr_qsroce/) )%laction = .TRUE.
      CASE( 'mixed oce-ice' )   ;   srcv(               jpr_qsrmix   )%laction = .TRUE.
      CASE default              ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_rcv_qsr%cldes' )
      END SELECT
      IF( TRIM( sn_rcv_qsr%cldes ) == 'mixed oce-ice' .AND. nn_cats_cpl > 1 ) &
         CALL ctl_stop( 'sbc_cpl_init: sn_rcv_qsr%cldes not currently allowed to be mixed oce-ice for multi-category ice' )
      !
      !                                                      ! ------------------------- !
      !                                                      !   non solar sensitivity   !   d(Qns)/d(T)
      !                                                      ! ------------------------- !
      srcv(jpr_dqnsdt)%clname = 'O_dQnsdT'
      IF( TRIM( sn_rcv_dqnsdt%cldes ) == 'coupled' )   srcv(jpr_dqnsdt)%laction = .TRUE.
      !
      ! non solar sensitivity mandatory for mixed oce-ice solar radiation coupling technique
      IF( TRIM( sn_rcv_dqnsdt%cldes ) == 'none' .AND. TRIM( sn_rcv_qns%cldes ) == 'mixed oce-ice' )  &
         &   CALL ctl_stop( 'sbc_cpl_init: namsbc_cpl namelist mismatch between sn_rcv_qns%cldes and sn_rcv_dqnsdt%cldes' )
      !
      !                                                      ! ------------------------- !
      !                                                      !      10m wind module      !
      !                                                      ! ------------------------- !
      srcv(jpr_w10m)%clname = 'O_Wind10'   ;   IF( TRIM(sn_rcv_w10m%cldes  ) == 'coupled' )   srcv(jpr_w10m)%laction = .TRUE.
      !
      !                                                      ! ------------------------- !
      !                                                      !   wind stress module      !
      !                                                      ! ------------------------- !
      srcv(jpr_taum)%clname = 'O_TauMod'   ;   IF( TRIM(sn_rcv_taumod%cldes) == 'coupled' )   srcv(jpr_taum)%laction = .TRUE.
      !
      !                                                      ! ------------------------- !
      !                                                      !      Atmospheric CO2      !
      !                                                      ! ------------------------- !
      srcv(jpr_co2 )%clname = 'O_AtmCO2'
      IF( TRIM(sn_rcv_co2%cldes   ) == 'coupled' )  THEN
         srcv(jpr_co2 )%laction = .TRUE.
         l_co2cpl = .TRUE.
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   Atmospheric pco2 received from oasis '
         IF(lwp) WRITE(numout,*)
      ENDIF
      !
      !                                                      ! ------------------------- !
      !                                                      ! Mean Sea Level Pressure   !
      !                                                      ! ------------------------- !
      srcv(jpr_mslp)%clname = 'O_MSLP'     ;   IF( TRIM(sn_rcv_mslp%cldes  ) == 'coupled' )    srcv(jpr_mslp)%laction = .TRUE.
      !
      !                                                      ! --------------------------------- !
      !                                                      !  ice topmelt and conduction flux  !   
      !                                                      ! --------------------------------- !
      srcv(jpr_topm )%clname = 'OTopMlt'
      srcv(jpr_botm )%clname = 'OBotMlt'
      IF( TRIM(sn_rcv_iceflx%cldes) == 'coupled' ) THEN
         IF( TRIM( sn_rcv_iceflx%clcat ) == 'yes' ) THEN
            srcv(jpr_topm:jpr_botm)%nct = nn_cats_cpl
         ELSE
            CALL ctl_stop( 'sbc_cpl_init: sn_rcv_iceflx%clcat should always be set to yes currently' )
         ENDIF
         srcv(jpr_topm:jpr_botm)%laction = .TRUE.
      ENDIF
!!      !                                                      ! --------------------------- !
!!      !                                                      ! transmitted solar thru ice  !   
!!      !                                                      ! --------------------------- !
!!      srcv(jpr_qtrice)%clname = 'OQtr'
!!      IF( TRIM(sn_rcv_qtrice%cldes) == 'coupled' ) THEN
!!         IF ( TRIM( sn_rcv_qtrice%clcat ) == 'yes' ) THEN
!!            srcv(jpr_qtrice)%nct = nn_cats_cpl
!!         ELSE
!!           CALL ctl_stop( 'sbc_cpl_init: sn_rcv_qtrice%clcat should always be set to yes currently' )
!!         ENDIF
!!         srcv(jpr_qtrice)%laction = .TRUE.
!!      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !    ice skin temperature   !
      !                                                      ! ------------------------- !
      srcv(jpr_ts_ice)%clname = 'OTsfIce'    ! needed by Met Office
      IF( TRIM( sn_rcv_ts_ice%cldes ) == 'ice' )   srcv(jpr_ts_ice)%laction = .TRUE.
      IF( TRIM( sn_rcv_ts_ice%clcat ) == 'yes' )   srcv(jpr_ts_ice)%nct     = nn_cats_cpl
      IF( TRIM( sn_rcv_emp%clcat    ) == 'yes' )   srcv(jpr_ievp)%nct       = nn_cats_cpl

      !                                                      ! ------------------------- !
      !                                                      !      Wave breaking        !
      !                                                      ! ------------------------- !
      srcv(jpr_hsig)%clname  = 'O_Hsigwa'    ! significant wave height
      IF( TRIM(sn_rcv_hsig%cldes  ) == 'coupled' )  THEN
         srcv(jpr_hsig)%laction = .TRUE.
         cpl_hsig = .TRUE.
      ENDIF
      srcv(jpr_phioc)%clname = 'O_PhiOce'    ! wave to ocean energy
      IF( TRIM(sn_rcv_phioc%cldes ) == 'coupled' )  THEN
         srcv(jpr_phioc)%laction = .TRUE.
         cpl_phioc = .TRUE.
      ENDIF
      srcv(jpr_sdrftx)%clname = 'O_Sdrfx'    ! Stokes drift in the u direction
      IF( TRIM(sn_rcv_sdrfx%cldes ) == 'coupled' )  THEN
         srcv(jpr_sdrftx)%laction = .TRUE.
         cpl_sdrftx = .TRUE.
      ENDIF
      srcv(jpr_sdrfty)%clname = 'O_Sdrfy'    ! Stokes drift in the v direction
      IF( TRIM(sn_rcv_sdrfy%cldes ) == 'coupled' )  THEN
         srcv(jpr_sdrfty)%laction = .TRUE.
         cpl_sdrfty = .TRUE.
      ENDIF
      srcv(jpr_wper)%clname = 'O_WPer'       ! mean wave period
      IF( TRIM(sn_rcv_wper%cldes  ) == 'coupled' )  THEN
         srcv(jpr_wper)%laction = .TRUE.
         cpl_wper = .TRUE.
      ENDIF
      srcv(jpr_wnum)%clname = 'O_WNum'       ! mean wave number
      IF( TRIM(sn_rcv_wnum%cldes ) == 'coupled' )  THEN
         srcv(jpr_wnum)%laction = .TRUE.
         cpl_wnum = .TRUE.
      ENDIF
      srcv(jpr_wstrf)%clname = 'O_WStrf'     ! stress fraction adsorbed by the wave
      IF( TRIM(sn_rcv_wstrf%cldes ) == 'coupled' )  THEN
         srcv(jpr_wstrf)%laction = .TRUE.
         cpl_wstrf = .TRUE.
      ENDIF
      srcv(jpr_wdrag)%clname = 'O_WDrag'     ! neutral surface drag coefficient
      IF( TRIM(sn_rcv_wdrag%cldes ) == 'coupled' )  THEN
         srcv(jpr_wdrag)%laction = .TRUE.
         cpl_wdrag = .TRUE.
      ENDIF
      srcv(jpr_charn)%clname = 'O_Charn'     ! Chranock coefficient
      IF( TRIM(sn_rcv_charn%cldes ) == 'coupled' )  THEN
         srcv(jpr_charn)%laction = .TRUE.
         cpl_charn = .TRUE.
      ENDIF
      srcv(jpr_bhd)%clname = 'O_Bhd'     ! Bernoulli head. waves' induced surface pressure
      IF( TRIM(sn_rcv_bhd%cldes ) == 'coupled' )  THEN
         srcv(jpr_bhd)%laction = .TRUE.
         cpl_bhd = .TRUE.
      ENDIF
      srcv(jpr_tusd)%clname = 'O_Tusd'     ! zonal stokes transport
      IF( TRIM(sn_rcv_tusd%cldes ) == 'coupled' )  THEN
         srcv(jpr_tusd)%laction = .TRUE.
         cpl_tusd = .TRUE.
      ENDIF
      srcv(jpr_tvsd)%clname = 'O_Tvsd'     ! meridional stokes tranmport
      IF( TRIM(sn_rcv_tvsd%cldes ) == 'coupled' )  THEN
         srcv(jpr_tvsd)%laction = .TRUE.
         cpl_tvsd = .TRUE.
      ENDIF

      srcv(jpr_twox)%clname = 'O_Twox'     ! wave to ocean momentum flux in the u direction
      srcv(jpr_twoy)%clname = 'O_Twoy'     ! wave to ocean momentum flux in the v direction
      srcv(jpr_tawx)%clname = 'O_Tawx'     ! Net wave-supported stress in the u direction
      srcv(jpr_tawy)%clname = 'O_Tawy'     ! Net wave-supported stress in the v direction
      IF( TRIM(sn_rcv_taw%cldes ) == 'coupled' )  THEN
         srcv(jpr_twox)%laction = .TRUE.
         srcv(jpr_twoy)%laction = .TRUE.
         srcv(jpr_tawx)%laction = .TRUE.
         srcv(jpr_tawy)%laction = .TRUE.
         cpl_taw = .TRUE.
      ENDIF
      !
      !                                                      ! ------------------------------- !
      !                                                      !   OCE-SAS coupling - rcv by opa !
      !                                                      ! ------------------------------- !
      srcv(jpr_sflx)%clname = 'O_SFLX'
      srcv(jpr_fice)%clname = 'RIceFrc'
      !
      IF( nn_components == jp_iam_oce ) THEN    ! OCE coupled to SAS via OASIS: force received field by OCE (sent by SAS)
         srcv(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         srcv(:)%clgrid  = 'T'       ! force default definition in case of opa <-> sas coupling
         srcv(:)%nsgn    = 1.        ! force default definition in case of opa <-> sas coupling
         srcv( (/jpr_qsroce, jpr_qnsoce, jpr_oemp, jpr_sflx, jpr_fice, jpr_otx1, jpr_oty1, jpr_taum/) )%laction = .TRUE.
         srcv(jpr_otx1)%clgrid = 'U'        ! oce components given at U-point
         srcv(jpr_oty1)%clgrid = 'V'        !           and           V-point
         ! Vectors: change of sign at north fold ONLY if on the local grid
         srcv( (/jpr_otx1,jpr_oty1/) )%nsgn = -1.
         sn_rcv_tau%clvgrd = 'U,V'
         sn_rcv_tau%clvor = 'local grid'
         sn_rcv_tau%clvref = 'spherical'
         sn_rcv_emp%cldes = 'oce only'
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'               Special conditions for SAS-OCE coupling  '
            WRITE(numout,*)'               OCE component  '
            WRITE(numout,*)
            WRITE(numout,*)'  received fields from SAS component '
            WRITE(numout,*)'                  ice cover '
            WRITE(numout,*)'                  oce only EMP  '
            WRITE(numout,*)'                  salt flux  '
            WRITE(numout,*)'                  mixed oce-ice solar flux  '
            WRITE(numout,*)'                  mixed oce-ice non solar flux  '
            WRITE(numout,*)'                  wind stress U,V on local grid and sperical coordinates '
            WRITE(numout,*)'                  wind stress module'
            WRITE(numout,*)
         ENDIF
      ENDIF
      !                                                      ! -------------------------------- !
      !                                                      !   OCE-SAS coupling - rcv by sas  !
      !                                                      ! -------------------------------- !
      srcv(jpr_toce  )%clname = 'I_SSTSST'
      srcv(jpr_soce  )%clname = 'I_SSSal'
      srcv(jpr_ocx1  )%clname = 'I_OCurx1'
      srcv(jpr_ocy1  )%clname = 'I_OCury1'
      srcv(jpr_ssh   )%clname = 'I_SSHght'
      srcv(jpr_e3t1st)%clname = 'I_E3T1st'
      srcv(jpr_fraqsr)%clname = 'I_FraQsr'
      !
      IF( nn_components == jp_iam_sas ) THEN
         IF( .NOT. ln_cpl ) srcv(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         IF( .NOT. ln_cpl ) srcv(:)%clgrid  = 'T'       ! force default definition in case of opa <-> sas coupling
         IF( .NOT. ln_cpl ) srcv(:)%nsgn    = 1.        ! force default definition in case of opa <-> sas coupling
         srcv( (/jpr_toce, jpr_soce, jpr_ssh, jpr_fraqsr, jpr_ocx1, jpr_ocy1/) )%laction = .TRUE.
         srcv( jpr_e3t1st )%laction = .NOT.ln_linssh
         srcv(jpr_ocx1)%clgrid = 'U'        ! oce components given at U-point
         srcv(jpr_ocy1)%clgrid = 'V'        !           and           V-point
         ! Vectors: change of sign at north fold ONLY if on the local grid
         srcv(jpr_ocx1:jpr_ocy1)%nsgn = -1.
         ! Change first letter to couple with atmosphere if already coupled OCE
         ! this is nedeed as each variable name used in the namcouple must be unique:
         ! for example O_Runoff received by OCE from SAS and therefore S_Runoff received by SAS from the Atmosphere
         DO jn = 1, jprcv
            IF( srcv(jn)%clname(1:1) == "O" ) srcv(jn)%clname = "S"//srcv(jn)%clname(2:LEN(srcv(jn)%clname))
         END DO
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'               Special conditions for SAS-OCE coupling  '
            WRITE(numout,*)'               SAS component  '
            WRITE(numout,*)
            IF( .NOT. ln_cpl ) THEN
               WRITE(numout,*)'  received fields from OCE component '
            ELSE
               WRITE(numout,*)'  Additional received fields from OCE component : '
            ENDIF
            WRITE(numout,*)'               sea surface temperature (Celsius) '
            WRITE(numout,*)'               sea surface salinity '
            WRITE(numout,*)'               surface currents '
            WRITE(numout,*)'               sea surface height '
            WRITE(numout,*)'               thickness of first ocean T level '
            WRITE(numout,*)'               fraction of solar net radiation absorbed in the first ocean level'
            WRITE(numout,*)
         ENDIF
      ENDIF

      ! =================================================== !
      ! Allocate all parts of frcv used for received fields !
      ! =================================================== !
      DO jn = 1, jprcv
         IF( srcv(jn)%laction ) ALLOCATE( frcv(jn)%z3(jpi,jpj,srcv(jn)%nct) )
      END DO
      ! Allocate taum part of frcv which is used even when not received as coupling field
      IF( .NOT. srcv(jpr_taum)%laction ) ALLOCATE( frcv(jpr_taum)%z3(jpi,jpj,srcv(jpr_taum)%nct) )
      ! Allocate w10m part of frcv which is used even when not received as coupling field
      IF( .NOT. srcv(jpr_w10m)%laction ) ALLOCATE( frcv(jpr_w10m)%z3(jpi,jpj,srcv(jpr_w10m)%nct) )
      ! Allocate jpr_otx1 part of frcv which is used even when not received as coupling field
      IF( .NOT. srcv(jpr_otx1)%laction ) ALLOCATE( frcv(jpr_otx1)%z3(jpi,jpj,srcv(jpr_otx1)%nct) )
      IF( .NOT. srcv(jpr_oty1)%laction ) ALLOCATE( frcv(jpr_oty1)%z3(jpi,jpj,srcv(jpr_oty1)%nct) )
      ! Allocate itx1 and ity1 as they are used in sbc_cpl_ice_tau even if srcv(jpr_itx1)%laction = .FALSE.
      IF( k_ice /= 0 ) THEN
         IF( .NOT. srcv(jpr_itx1)%laction ) ALLOCATE( frcv(jpr_itx1)%z3(jpi,jpj,srcv(jpr_itx1)%nct) )
         IF( .NOT. srcv(jpr_ity1)%laction ) ALLOCATE( frcv(jpr_ity1)%z3(jpi,jpj,srcv(jpr_ity1)%nct) )
      ENDIF

      ! ================================ !
      !     Define the send interface    !
      ! ================================ !
      ! for each field: define the OASIS name                           (ssnd(:)%clname)
      !                 define send or not from the namelist parameters (ssnd(:)%laction)
      !                 define the north fold type of lbc               (ssnd(:)%nsgn)

      ! default definitions of nsnd
      ssnd(:)%laction = .FALSE.   ;   ssnd(:)%clgrid = 'T'   ;   ssnd(:)%nsgn = 1.  ; ssnd(:)%nct = 1

      !                                                      ! ------------------------- !
      !                                                      !    Surface temperature    !
      !                                                      ! ------------------------- !
      ssnd(jps_toce)%clname   = 'O_SSTSST'
      ssnd(jps_tice)%clname   = 'O_TepIce'
      ssnd(jps_ttilyr)%clname = 'O_TtiLyr'
      ssnd(jps_tmix)%clname   = 'O_TepMix'
      SELECT CASE( TRIM( sn_snd_temp%cldes ) )
      CASE( 'none'                                 )       ! nothing to do
      CASE( 'oce only'                             )   ;   ssnd( jps_toce )%laction = .TRUE.
      CASE( 'oce and ice' , 'weighted oce and ice' , 'oce and weighted ice' )
         ssnd( (/jps_toce, jps_tice/) )%laction = .TRUE.
         IF( TRIM( sn_snd_temp%clcat ) == 'yes' )  ssnd(jps_tice)%nct = nn_cats_cpl
      CASE( 'mixed oce-ice'                        )   ;   ssnd( jps_tmix )%laction = .TRUE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_temp%cldes' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !          Albedo           !
      !                                                      ! ------------------------- !
      ssnd(jps_albice)%clname = 'O_AlbIce'
      ssnd(jps_albmix)%clname = 'O_AlbMix'
      SELECT CASE( TRIM( sn_snd_alb%cldes ) )
      CASE( 'none'                 )     ! nothing to do
      CASE( 'ice' , 'weighted ice' )   ; ssnd(jps_albice)%laction = .TRUE.
      CASE( 'mixed oce-ice'        )   ; ssnd(jps_albmix)%laction = .TRUE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_alb%cldes' )
      END SELECT
      !
      ! Need to calculate oceanic albedo if
      !     1. sending mixed oce-ice albedo or
      !     2. receiving mixed oce-ice solar radiation
      IF( TRIM ( sn_snd_alb%cldes ) == 'mixed oce-ice' .OR. TRIM ( sn_rcv_qsr%cldes ) == 'mixed oce-ice' ) THEN
         CALL oce_alb( zaos, zacs )
         ! Due to lack of information on nebulosity : mean clear/overcast sky
         alb_oce_mix(:,:) = ( zacs(:,:) + zaos(:,:) ) * 0.5
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !  Ice fraction & Thickness !
      !                                                      ! ------------------------- !
      ssnd(jps_fice)%clname  = 'OIceFrc'
      ssnd(jps_ficet)%clname = 'OIceFrcT'
      ssnd(jps_hice)%clname  = 'OIceTck'
      ssnd(jps_a_p)%clname   = 'OPndFrc'
      ssnd(jps_ht_p)%clname  = 'OPndTck'
      ssnd(jps_hsnw)%clname  = 'OSnwTck'
      ssnd(jps_fice1)%clname = 'OIceFrd'
      IF( k_ice /= 0 ) THEN
         ssnd(jps_fice)%laction  = .TRUE.                 ! if ice treated in the ocean (even in climato case)
         ssnd(jps_fice1)%laction = .TRUE.                 ! First-order regridded ice concentration, to be used producing atmos-to-ice fluxes (Met Office requirement)
! Currently no namelist entry to determine sending of multi-category ice fraction so use the thickness entry for now
         IF( TRIM( sn_snd_thick%clcat  ) == 'yes' ) ssnd(jps_fice)%nct  = nn_cats_cpl
         IF( TRIM( sn_snd_thick1%clcat ) == 'yes' ) ssnd(jps_fice1)%nct = nn_cats_cpl
      ENDIF

      IF(TRIM( sn_snd_ifrac%cldes )  == 'coupled') ssnd(jps_ficet)%laction = .TRUE.

      SELECT CASE ( TRIM( sn_snd_thick%cldes ) )
      CASE( 'none'         )       ! nothing to do
      CASE( 'ice and snow' )
         ssnd(jps_hice:jps_hsnw)%laction = .TRUE.
         IF( TRIM( sn_snd_thick%clcat ) == 'yes' ) THEN
            ssnd(jps_hice:jps_hsnw)%nct = nn_cats_cpl
         ENDIF
      CASE ( 'weighted ice and snow' )
         ssnd(jps_hice:jps_hsnw)%laction = .TRUE.
         IF( TRIM( sn_snd_thick%clcat ) == 'yes' ) ssnd(jps_hice:jps_hsnw)%nct = nn_cats_cpl
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_thick%cldes' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !      Ice Meltponds        !
      !                                                      ! ------------------------- !
      ! Needed by Met Office
      ssnd(jps_a_p)%clname  = 'OPndFrc'
      ssnd(jps_ht_p)%clname = 'OPndTck'
      SELECT CASE ( TRIM( sn_snd_mpnd%cldes ) )
      CASE ( 'none' )
         ssnd(jps_a_p)%laction  = .FALSE.
         ssnd(jps_ht_p)%laction = .FALSE.
      CASE ( 'ice only' )
         ssnd(jps_a_p)%laction  = .TRUE.
         ssnd(jps_ht_p)%laction = .TRUE.
         IF( TRIM( sn_snd_mpnd%clcat ) == 'yes' ) THEN
            ssnd(jps_a_p)%nct  = nn_cats_cpl
            ssnd(jps_ht_p)%nct = nn_cats_cpl
         ELSE
            IF( nn_cats_cpl > 1 ) THEN
               CALL ctl_stop( 'sbc_cpl_init: use weighted ice option for sn_snd_mpnd%cldes if not exchanging category fields' )
            ENDIF
         ENDIF
      CASE ( 'weighted ice' )
         ssnd(jps_a_p)%laction  = .TRUE.
         ssnd(jps_ht_p)%laction = .TRUE.
         IF( TRIM( sn_snd_mpnd%clcat ) == 'yes' ) THEN
            ssnd(jps_a_p)%nct  = nn_cats_cpl
            ssnd(jps_ht_p)%nct = nn_cats_cpl
         ENDIF
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_mpnd%cldes; '//sn_snd_mpnd%cldes )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !      Surface current      !
      !                                                      ! ------------------------- !
      !        ocean currents              !            ice velocities
      ssnd(jps_ocx1)%clname = 'O_OCurx1'   ;   ssnd(jps_ivx1)%clname = 'O_IVelx1'
      ssnd(jps_ocy1)%clname = 'O_OCury1'   ;   ssnd(jps_ivy1)%clname = 'O_IVely1'
      ssnd(jps_ocz1)%clname = 'O_OCurz1'   ;   ssnd(jps_ivz1)%clname = 'O_IVelz1'
      ssnd(jps_ocxw)%clname = 'O_OCurxw'
      ssnd(jps_ocyw)%clname = 'O_OCuryw'
      !
      ssnd(jps_ocx1:jps_ivz1)%nsgn = -1.   ! vectors: change of the sign at the north fold

      IF( sn_snd_crt%clvgrd == 'U,V' ) THEN
         ssnd(jps_ocx1)%clgrid = 'U' ; ssnd(jps_ocy1)%clgrid = 'V'
      ELSE IF( sn_snd_crt%clvgrd /= 'T' ) THEN
         CALL ctl_stop( 'sn_snd_crt%clvgrd must be equal to T' )
      ENDIF
      ssnd(jps_ocx1:jps_ivz1)%laction = .TRUE.   ! default: all are send
      IF( TRIM( sn_snd_crt%clvref ) == 'spherical' )   ssnd( (/jps_ocz1, jps_ivz1/) )%laction = .FALSE.
      IF( TRIM( sn_snd_crt%clvor ) == 'eastward-northward' ) ssnd(jps_ocx1:jps_ivz1)%nsgn = 1.
      SELECT CASE( TRIM( sn_snd_crt%cldes ) )
      CASE( 'none'                 )   ;   ssnd(jps_ocx1:jps_ivz1)%laction = .FALSE.
      CASE( 'oce only'             )   ;   ssnd(jps_ivx1:jps_ivz1)%laction = .FALSE.
      CASE( 'weighted oce and ice' )   !   nothing to do
      CASE( 'mixed oce-ice'        )   ;   ssnd(jps_ivx1:jps_ivz1)%laction = .FALSE.
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_crt%cldes' )
      END SELECT

      ssnd(jps_ocxw:jps_ocyw)%nsgn = -1.   ! vectors: change of the sign at the north fold

      IF( sn_snd_crtw%clvgrd == 'U,V' ) THEN
         ssnd(jps_ocxw)%clgrid = 'U' ; ssnd(jps_ocyw)%clgrid = 'V'
      ELSE IF( sn_snd_crtw%clvgrd /= 'T' ) THEN
         CALL ctl_stop( 'sn_snd_crtw%clvgrd must be equal to T' )
      ENDIF
      IF( TRIM( sn_snd_crtw%clvor ) == 'eastward-northward' ) ssnd(jps_ocxw:jps_ocyw)%nsgn = 1.
      SELECT CASE( TRIM( sn_snd_crtw%cldes ) )
         CASE( 'none'                 )   ; ssnd(jps_ocxw:jps_ocyw)%laction = .FALSE.
         CASE( 'oce only'             )   ; ssnd(jps_ocxw:jps_ocyw)%laction = .TRUE.
         CASE( 'weighted oce and ice' )   !   nothing to do
         CASE( 'mixed oce-ice'        )   ; ssnd(jps_ivx1:jps_ivz1)%laction = .FALSE.
         CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_crtw%cldes' )
      END SELECT

      !                                                      ! ------------------------- !
      !                                                      !          CO2 flux         !
      !                                                      ! ------------------------- !
      ssnd(jps_co2)%clname = 'O_CO2FLX' ;  IF( TRIM(sn_snd_co2%cldes) == 'coupled' )    ssnd(jps_co2 )%laction = .TRUE.
      !
      !                                                      ! ------------------------- !
      !                                                      ! Sea surface freezing temp !
      !                                                      ! ------------------------- !
      ! needed by Met Office
      ssnd(jps_sstfrz)%clname = 'O_SSTFrz' ; IF( TRIM(sn_snd_sstfrz%cldes) == 'coupled' )  ssnd(jps_sstfrz)%laction = .TRUE.
      !
      !                                                      ! ------------------------- !
      !                                                      !    Ice conductivity       !
      !                                                      ! ------------------------- !
      ! needed by Met Office
      ! Note that ultimately we will move to passing an ocean effective conductivity as well so there
      ! will be some changes to the parts of the code which currently relate only to ice conductivity
      ssnd(jps_ttilyr )%clname = 'O_TtiLyr'
      SELECT CASE ( TRIM( sn_snd_ttilyr%cldes ) )
      CASE ( 'none' )
         ssnd(jps_ttilyr)%laction = .FALSE.
      CASE ( 'ice only' )
         ssnd(jps_ttilyr)%laction = .TRUE.
         IF( TRIM( sn_snd_ttilyr%clcat ) == 'yes' ) THEN
            ssnd(jps_ttilyr)%nct = nn_cats_cpl
         ELSE
            IF( nn_cats_cpl > 1 ) THEN
               CALL ctl_stop( 'sbc_cpl_init: use weighted ice option for sn_snd_ttilyr%cldes if not exchanging category fields' )
            ENDIF
         ENDIF
      CASE ( 'weighted ice' )
         ssnd(jps_ttilyr)%laction = .TRUE.
         IF( TRIM( sn_snd_ttilyr%clcat ) == 'yes' ) ssnd(jps_ttilyr)%nct = nn_cats_cpl
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_ttilyr%cldes;'//sn_snd_ttilyr%cldes )
      END SELECT

      ssnd(jps_kice )%clname = 'OIceKn'
      SELECT CASE ( TRIM( sn_snd_cond%cldes ) )
      CASE ( 'none' )
         ssnd(jps_kice)%laction = .FALSE.
      CASE ( 'ice only' )
         ssnd(jps_kice)%laction = .TRUE.
         IF( TRIM( sn_snd_cond%clcat ) == 'yes' ) THEN
            ssnd(jps_kice)%nct = nn_cats_cpl
         ELSE
            IF( nn_cats_cpl > 1 ) THEN
               CALL ctl_stop( 'sbc_cpl_init: use weighted ice option for sn_snd_cond%cldes if not exchanging category fields' )
            ENDIF
         ENDIF
      CASE ( 'weighted ice' )
         ssnd(jps_kice)%laction = .TRUE.
         IF( TRIM( sn_snd_cond%clcat ) == 'yes' ) ssnd(jps_kice)%nct = nn_cats_cpl
      CASE default   ;   CALL ctl_stop( 'sbc_cpl_init: wrong definition of sn_snd_cond%cldes;'//sn_snd_cond%cldes )
      END SELECT
      !
      !                                                      ! ------------------------- !
      !                                                      !     Sea surface height    !
      !                                                      ! ------------------------- !
      ssnd(jps_wlev)%clname = 'O_Wlevel' ;  IF( TRIM(sn_snd_wlev%cldes) == 'coupled' )   ssnd(jps_wlev)%laction = .TRUE.

      !                                                      ! ------------------------------- !
      !                                                      !   OCE-SAS coupling - snd by opa !
      !                                                      ! ------------------------------- !
      ssnd(jps_ssh   )%clname = 'O_SSHght'
      ssnd(jps_soce  )%clname = 'O_SSSal'
      ssnd(jps_e3t1st)%clname = 'O_E3T1st'
      ssnd(jps_fraqsr)%clname = 'O_FraQsr'
      !
      IF( nn_components == jp_iam_oce ) THEN
         ssnd(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         ssnd( (/jps_toce, jps_soce, jps_ssh, jps_fraqsr, jps_ocx1, jps_ocy1/) )%laction = .TRUE.
         ssnd( jps_e3t1st )%laction = .NOT.ln_linssh
         ! vector definition: not used but cleaner...
         ssnd(jps_ocx1)%clgrid  = 'U'        ! oce components given at U-point
         ssnd(jps_ocy1)%clgrid  = 'V'        !           and           V-point
         sn_snd_crt%clvgrd = 'U,V'
         sn_snd_crt%clvor = 'local grid'
         sn_snd_crt%clvref = 'spherical'
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            WRITE(numout,*)'  sent fields to SAS component '
            WRITE(numout,*)'               sea surface temperature (T before, Celsius) '
            WRITE(numout,*)'               sea surface salinity '
            WRITE(numout,*)'               surface currents U,V on local grid and spherical coordinates'
            WRITE(numout,*)'               sea surface height '
            WRITE(numout,*)'               thickness of first ocean T level '
            WRITE(numout,*)'               fraction of solar net radiation absorbed in the first ocean level'
            WRITE(numout,*)
         ENDIF
      ENDIF
      !                                                      ! ------------------------------- !
      !                                                      !   OCE-SAS coupling - snd by sas !
      !                                                      ! ------------------------------- !
      ssnd(jps_sflx  )%clname = 'I_SFLX'
      ssnd(jps_fice2 )%clname = 'IIceFrc'
      ssnd(jps_qsroce)%clname = 'I_QsrOce'
      ssnd(jps_qnsoce)%clname = 'I_QnsOce'
      ssnd(jps_oemp  )%clname = 'IOEvaMPr'
      ssnd(jps_otx1  )%clname = 'I_OTaux1'
      ssnd(jps_oty1  )%clname = 'I_OTauy1'
      ssnd(jps_rnf   )%clname = 'I_Runoff'
      ssnd(jps_taum  )%clname = 'I_TauMod'
      !
      IF( nn_components == jp_iam_sas ) THEN
         IF( .NOT. ln_cpl ) ssnd(:)%laction = .FALSE.   ! force default definition in case of opa <-> sas coupling
         ssnd( (/jps_qsroce, jps_qnsoce, jps_oemp, jps_fice2, jps_sflx, jps_otx1, jps_oty1, jps_taum/) )%laction = .TRUE.
         !
         ! Change first letter to couple with atmosphere if already coupled with sea_ice
         ! this is nedeed as each variable name used in the namcouple must be unique:
         ! for example O_SSTSST sent by OCE to SAS and therefore S_SSTSST sent by SAS to the Atmosphere
         DO jn = 1, jpsnd
            IF( ssnd(jn)%clname(1:1) == "O" ) ssnd(jn)%clname = "S"//ssnd(jn)%clname(2:LEN(ssnd(jn)%clname))
         END DO
         !
         IF(lwp) THEN                        ! control print
            WRITE(numout,*)
            IF( .NOT. ln_cpl ) THEN
               WRITE(numout,*)'  sent fields to OCE component '
            ELSE
               WRITE(numout,*)'  Additional sent fields to OCE component : '
            ENDIF
            WRITE(numout,*)'                  ice cover '
            WRITE(numout,*)'                  oce only EMP  '
            WRITE(numout,*)'                  salt flux  '
            WRITE(numout,*)'                  mixed oce-ice solar flux  '
            WRITE(numout,*)'                  mixed oce-ice non solar flux  '
            WRITE(numout,*)'                  wind stress U,V components'
            WRITE(numout,*)'                  wind stress module'
         ENDIF
      ENDIF

      !
      ! ================================ !
      !   initialisation of the coupler  !
      ! ================================ !
      CALL cpl_define(jprcv, jpsnd, nn_cplmodel)

      IF(ln_usecplmask) THEN
         xcplmask(:,:,:) = 0.
         CALL iom_open( 'cplmask', inum )
         CALL iom_get( inum, jpdom_unknown, 'cplmask', xcplmask(1:jpi,1:jpj,1:nn_cplmodel),   &
            &          kstart = (/ mig(1),mjg(1),1 /), kcount = (/ jpi,jpj,nn_cplmodel /) )
         CALL iom_close( inum )
      ELSE
         xcplmask(:,:,:) = 1.
      ENDIF
      xcplmask(:,:,0) = 1. - SUM( xcplmask(:,:,1:nn_cplmodel), dim = 3 )
      !
      !
   END SUBROUTINE sbc_cpl_init


   SUBROUTINE sbc_cpl_rcv( kt, k_fsbc, k_ice, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_rcv  ***
      !!
      !! ** Purpose :   provide the stress over the ocean and, if no sea-ice,
      !!                provide the ocean heat and freshwater fluxes.
      !!
      !! ** Method  : - Receive all the atmospheric fields (stored in frcv array). called at each time step.
      !!                OASIS controls if there is something do receive or not. nrcvinfo contains the info
      !!                to know if the field was really received or not
      !!
      !!              --> If ocean stress was really received:
      !!
      !!                  - transform the received ocean stress vector from the received
      !!                 referential and grid into an atmosphere-ocean stress in
      !!                 the (i,j) ocean referencial and at the ocean velocity point.
      !!                    The received stress are :
      !!                     - defined by 3 components (if cartesian coordinate)
      !!                            or by 2 components (if spherical)
      !!                     - oriented along geographical   coordinate (if eastward-northward)
      !!                            or  along the local grid coordinate (if local grid)
      !!                     - given at U- and V-point, resp.   if received on 2 grids
      !!                            or at T-point               if received on 1 grid
      !!                    Therefore and if necessary, they are successively
      !!                  processed in order to obtain them
      !!                     first  as  2 components on the sphere
      !!                     second as  2 components oriented along the local grid
      !!                     third  as  2 components on the U,V grid
      !!
      !!              -->
      !!
      !!              - In 'ocean only' case, non solar and solar ocean heat fluxes
      !!             and total ocean freshwater fluxes
      !!
      !! ** Method  :   receive all fields from the atmosphere and transform
      !!              them into ocean surface boundary condition fields
      !!
      !! ** Action  :   update  utau, vtau   ocean stress at U,V grid
      !!                        taum         wind stress module at T-point
      !!                        wndm         wind speed  module at T-point over free ocean or leads in presence of sea-ice
      !!                        qns          non solar heat fluxes including emp heat content    (ocean only case)
      !!                                     and the latent heat flux of solid precip. melting
      !!                        qsr          solar ocean heat fluxes   (ocean only case)
      !!                        emp          upward mass flux [evap. - precip. (- runoffs) (- calving)] (ocean only case)
      !!----------------------------------------------------------------------
      USE zdf_oce,  ONLY :   ln_zdfswm
      !
      INTEGER, INTENT(in) ::   kt          ! ocean model time step index
      INTEGER, INTENT(in) ::   k_fsbc      ! frequency of sbc (-> ice model) computation
      INTEGER, INTENT(in) ::   k_ice       ! ice management in the sbc (=0/1/2/3)
      INTEGER, INTENT(in) ::   Kbb, Kmm    ! ocean model time level indices
      !!
      LOGICAL  ::   llnewtx, llnewtau      ! update wind stress components and module??
      INTEGER  ::   ji, jj, jn             ! dummy loop indices
      INTEGER  ::   isec                   ! number of seconds since nit000 (assuming rdt did not change since nit000)
      REAL(wp) ::   zcumulneg, zcumulpos   ! temporary scalars
      REAL(wp) ::   zcoef                  ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22          ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3        ! drag coefficient
      REAL(wp) ::   zzx, zzy               ! temporary variables
      REAL(wp) ::   r1_grau                ! = 1.e0 / (grav * rho0)
      REAL(wp), DIMENSION(jpi,jpj) ::   ztx, zty, zmsk, zemp, zqns, zqsr, zcloud_fra
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
      !   cannot be done in the init phase when we use agrif as cpl_freq requires that oasis_enddef is done
         ncpl_qsr_freq = cpl_freq( 'O_QsrOce' ) + cpl_freq( 'O_QsrMix' ) + cpl_freq( 'I_QsrOce' ) + cpl_freq( 'I_QsrMix' )
         IF( ln_dm2dc .AND. ncpl_qsr_freq /= 86400 )   &
            &   CALL ctl_stop( 'sbc_cpl_rcv: diurnal cycle reconstruction (ln_dm2dc) needs daily couping for solar radiation' )

         IF ( ln_wave .AND. nn_components == 0 ) THEN
            ncpl_qsr_freq = 1;
            WRITE(numout,*) 'ncpl_qsr_freq is set to 1 when coupling NEMO with wave (without SAS) '
         ENDIF
      ENDIF
      !
      IF( ln_mixcpl )   zmsk(:,:) = 1. - xcplmask(:,:,0)
      !
      !                                                      ! ======================================================= !
      !                                                      ! Receive all the atmos. fields (including ice information)
      !                                                      ! ======================================================= !
      isec = ( kt - nit000 ) * NINT( rn_Dt )                      ! date of exchanges
      DO jn = 1, jprcv                                          ! received fields sent by the atmosphere
         IF( srcv(jn)%laction )   CALL cpl_rcv( jn, isec, frcv(jn)%z3, xcplmask(:,:,1:nn_cplmodel), nrcvinfo(jn) )
      END DO

      !                                                      ! ========================= !
      IF( srcv(jpr_otx1)%laction ) THEN                      !  ocean stress components  !
         !                                                   ! ========================= !
         ! define frcv(jpr_otx1)%z3(:,:,1) and frcv(jpr_oty1)%z3(:,:,1): stress at U/V point along model grid
         ! => need to be done only when we receive the field
         IF(  nrcvinfo(jpr_otx1) == OASIS_Rcv ) THEN
            !
            IF( TRIM( sn_rcv_tau%clvref ) == 'cartesian' ) THEN            ! 2 components on the sphere
               !                                                       ! (cartesian to spherical -> 3 to 2 components)
               !
               CALL geo2oce( frcv(jpr_otx1)%z3(:,:,1), frcv(jpr_oty1)%z3(:,:,1), frcv(jpr_otz1)%z3(:,:,1),   &
                  &          srcv(jpr_otx1)%clgrid, ztx, zty )
               frcv(jpr_otx1)%z3(:,:,1) = ztx(:,:)   ! overwrite 1st comp. on the 1st grid
               frcv(jpr_oty1)%z3(:,:,1) = zty(:,:)   ! overwrite 2nd comp. on the 1st grid
               !
               IF( srcv(jpr_otx2)%laction ) THEN
                  CALL geo2oce( frcv(jpr_otx2)%z3(:,:,1), frcv(jpr_oty2)%z3(:,:,1), frcv(jpr_otz2)%z3(:,:,1),   &
                     &          srcv(jpr_otx2)%clgrid, ztx, zty )
                  frcv(jpr_otx2)%z3(:,:,1) = ztx(:,:)   ! overwrite 1st comp. on the 2nd grid
                  frcv(jpr_oty2)%z3(:,:,1) = zty(:,:)   ! overwrite 2nd comp. on the 2nd grid
               ENDIF
               !
            ENDIF
            !
            IF( TRIM( sn_rcv_tau%clvor ) == 'eastward-northward' ) THEN   ! 2 components oriented along the local grid
               !                                                       ! (geographical to local grid -> rotate the components)
               CALL rot_rep( frcv(jpr_otx1)%z3(:,:,1), frcv(jpr_oty1)%z3(:,:,1), srcv(jpr_otx1)%clgrid, 'en->i', ztx )
               IF( srcv(jpr_otx2)%laction ) THEN
                  CALL rot_rep( frcv(jpr_otx2)%z3(:,:,1), frcv(jpr_oty2)%z3(:,:,1), srcv(jpr_otx2)%clgrid, 'en->j', zty )
               ELSE
                  CALL rot_rep( frcv(jpr_otx1)%z3(:,:,1), frcv(jpr_oty1)%z3(:,:,1), srcv(jpr_otx1)%clgrid, 'en->j', zty )
               ENDIF
               frcv(jpr_otx1)%z3(:,:,1) = ztx(:,:)      ! overwrite 1st component on the 1st grid
               frcv(jpr_oty1)%z3(:,:,1) = zty(:,:)      ! overwrite 2nd component on the 2nd grid
            ENDIF
            !
            IF( srcv(jpr_otx1)%clgrid == 'T' ) THEN
               DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)                                        ! T ==> (U,V)
                  frcv(jpr_otx1)%z3(ji,jj,1) = 0.5 * ( frcv(jpr_otx1)%z3(ji+1,jj  ,1) + frcv(jpr_otx1)%z3(ji,jj,1) )
                  frcv(jpr_oty1)%z3(ji,jj,1) = 0.5 * ( frcv(jpr_oty1)%z3(ji  ,jj+1,1) + frcv(jpr_oty1)%z3(ji,jj,1) )
               END DO   ;   END DO
               CALL lbc_lnk( 'sbccpl', frcv(jpr_otx1)%z3(:,:,1), 'U',  -1.0_wp, frcv(jpr_oty1)%z3(:,:,1), 'V',  -1.0_wp )
            ENDIF
            llnewtx = .TRUE.
         ELSE
            llnewtx = .FALSE.
         ENDIF
         !                                                   ! ========================= !
      ELSE                                                   !   No dynamical coupling   !
         !                                                   ! ========================= !
         frcv(jpr_otx1)%z3(:,:,1) = 0.e0                               ! here simply set to zero
         frcv(jpr_oty1)%z3(:,:,1) = 0.e0                               ! an external read in a file can be added instead
         llnewtx = .TRUE.
         !
      ENDIF
      !                                                      ! ========================= !
      !                                                      !    wind stress module     !   (taum)
      !                                                      ! ========================= !
      IF( .NOT. srcv(jpr_taum)%laction ) THEN                    ! compute wind stress module from its components if not received
         ! => need to be done only when otx1 was changed
         IF( llnewtx ) THEN
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               zzx = frcv(jpr_otx1)%z3(ji-1,jj  ,1) + frcv(jpr_otx1)%z3(ji,jj,1)
               zzy = frcv(jpr_oty1)%z3(ji  ,jj-1,1) + frcv(jpr_oty1)%z3(ji,jj,1)
               frcv(jpr_taum)%z3(ji,jj,1) = 0.5 * SQRT( zzx * zzx + zzy * zzy )
            END DO   ;   END DO
            CALL lbc_lnk( 'sbccpl', frcv(jpr_taum)%z3(:,:,1), 'T', 1.0_wp )
            llnewtau = .TRUE.
         ELSE
            llnewtau = .FALSE.
         ENDIF
      ELSE
         llnewtau = nrcvinfo(jpr_taum) == OASIS_Rcv
         ! Stress module can be negative when received (interpolation problem)
         IF( llnewtau ) THEN
            frcv(jpr_taum)%z3(:,:,1) = MAX( 0._wp, frcv(jpr_taum)%z3(:,:,1) )
         ENDIF
      ENDIF
      !
      !                                                      ! ========================= !
      !                                                      !      10 m wind speed      !   (wndm)
      !                                                      ! ========================= !
      IF( .NOT. srcv(jpr_w10m)%laction ) THEN                    ! compute wind spreed from wind stress module if not received
         ! => need to be done only when taumod was changed
         IF( llnewtau ) THEN
            zcoef = 1. / ( zrhoa * zcdrag )
            DO jj = ntsj-( nn_hls), ntej+( nn_hls ) ; DO ji = ntsi-( nn_hls), ntei+( nn_hls)
               frcv(jpr_w10m)%z3(ji,jj,1) = SQRT( frcv(jpr_taum)%z3(ji,jj,1) * zcoef )
            END DO   ;   END DO
         ENDIF
      ENDIF
!!$      !                                                      ! ========================= !
!!$      SELECT CASE( TRIM( sn_rcv_clouds%cldes ) )             !       cloud fraction      !
!!$      !                                                      ! ========================= !
!!$      cloud_fra(:,:) = frcv(jpr_clfra)*z3(:,:,1)
!!$      END SELECT
!!$
      zcloud_fra(:,:) = pp_cldf   ! should be real cloud fraction instead (as in the bulk) but needs to be read from atm.
      IF( ln_mixcpl ) THEN
         cloud_fra(:,:) = cloud_fra(:,:) * xcplmask(:,:,0) + zcloud_fra(:,:)* zmsk(:,:)
      ELSE
         cloud_fra(:,:) = zcloud_fra(:,:)
      ENDIF
      !                                                      ! ========================= !
      ! u(v)tau and taum will be modified by ice model
      ! -> need to be reset before each call of the ice/fsbc
      IF( MOD( kt-1, k_fsbc ) == 0 ) THEN
         !
         IF( ln_mixcpl ) THEN
            utau(:,:) = utau(:,:) * xcplmask(:,:,0) + frcv(jpr_otx1)%z3(:,:,1) * zmsk(:,:)
            vtau(:,:) = vtau(:,:) * xcplmask(:,:,0) + frcv(jpr_oty1)%z3(:,:,1) * zmsk(:,:)
            taum(:,:) = taum(:,:) * xcplmask(:,:,0) + frcv(jpr_taum)%z3(:,:,1) * zmsk(:,:)
            wndm(:,:) = wndm(:,:) * xcplmask(:,:,0) + frcv(jpr_w10m)%z3(:,:,1) * zmsk(:,:)
         ELSE
            utau(:,:) = frcv(jpr_otx1)%z3(:,:,1)
            vtau(:,:) = frcv(jpr_oty1)%z3(:,:,1)
            taum(:,:) = frcv(jpr_taum)%z3(:,:,1)
            wndm(:,:) = frcv(jpr_w10m)%z3(:,:,1)
         ENDIF
         CALL iom_put( "taum_oce", taum )   ! output wind stress module
         !
      ENDIF

      !                                                      ! ================== !
      !                                                      ! atmosph. CO2 (ppm) !
      !                                                      ! ================== !
      IF( srcv(jpr_co2)%laction )   atm_co2(:,:) = frcv(jpr_co2)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      ! Mean Sea Level Pressure   !   (taum)
      !                                                      ! ========================= !
      IF( srcv(jpr_mslp)%laction ) THEN                    ! UKMO SHELF effect of atmospheric pressure on SSH
          IF( kt /= nit000 )   ssh_ibb(:,:) = ssh_ib(:,:)    !* Swap of ssh_ib fields

          r1_grau = 1.e0 / (grav * rho0)               !* constant for optimization
          ssh_ib(:,:) = - ( frcv(jpr_mslp)%z3(:,:,1) - rpref ) * r1_grau    ! equivalent ssh (inverse barometer)
          apr   (:,:) =     frcv(jpr_mslp)%z3(:,:,1)                         !atmospheric pressure

          IF( kt == nit000 ) ssh_ibb(:,:) = ssh_ib(:,:)  ! correct this later (read from restart if possible)
      ENDIF
      !
      IF( ln_sdw ) THEN  ! Stokes Drift correction activated
      !                                                      ! ========================= !
      !                                                      !       Stokes drift u      !
      !                                                      ! ========================= !
         IF( srcv(jpr_sdrftx)%laction ) ut0sd(:,:) = frcv(jpr_sdrftx)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !       Stokes drift v      !
      !                                                      ! ========================= !
         IF( srcv(jpr_sdrfty)%laction ) vt0sd(:,:) = frcv(jpr_sdrfty)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !      Wave mean period     !
      !                                                      ! ========================= !
         IF( srcv(jpr_wper)%laction ) wmp(:,:) = frcv(jpr_wper)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !  Significant wave height  !
      !                                                      ! ========================= !
         IF( srcv(jpr_hsig)%laction ) hsw(:,:) = frcv(jpr_hsig)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !    Vertical mixing Qiao   !
      !                                                      ! ========================= !
         IF( srcv(jpr_wnum)%laction .AND. ln_zdfswm ) wnum(:,:) = frcv(jpr_wnum)%z3(:,:,1)

         ! Calculate the 3D Stokes drift both in coupled and not fully uncoupled mode
         IF( srcv(jpr_sdrftx)%laction .OR. srcv(jpr_sdrfty)%laction .OR. &
             srcv(jpr_wper)%laction .OR. srcv(jpr_hsig)%laction )   THEN
            CALL sbc_stokes( Kmm )
         ENDIF
      ENDIF
      !                                                      ! ========================= !
      !                                                      ! Stress adsorbed by waves  !
      !                                                      ! ========================= !
      IF( srcv(jpr_wstrf)%laction .AND. ln_tauoc )  tauoc_wave(:,:) = frcv(jpr_wstrf)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !   Wave drag coefficient   !
      !                                                      ! ========================= !
      IF( srcv(jpr_wdrag)%laction .AND. ln_cdgw )   cdn_wave(:,:) = frcv(jpr_wdrag)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !   Chranock coefficient    !
      !                                                      ! ========================= !
      IF( srcv(jpr_charn)%laction .AND. ln_charn )  charn(:,:) = frcv(jpr_charn)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      ! net wave-supported stress !
      !                                                      ! ========================= !
      IF( srcv(jpr_tawx)%laction .AND. ln_taw )     tawx(:,:) = frcv(jpr_tawx)%z3(:,:,1)
      IF( srcv(jpr_tawy)%laction .AND. ln_taw )     tawy(:,:) = frcv(jpr_tawy)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !wave to ocean momentum flux!
      !                                                      ! ========================= !
      IF( srcv(jpr_twox)%laction .AND. ln_taw )     twox(:,:) = frcv(jpr_twox)%z3(:,:,1)
      IF( srcv(jpr_twoy)%laction .AND. ln_taw )     twoy(:,:) = frcv(jpr_twoy)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !    wave TKE flux at sfc   !
      !                                                      ! ========================= !
      IF( srcv(jpr_phioc)%laction .AND. ln_phioc )     phioc(:,:) = frcv(jpr_phioc)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !      Bernoulli head       !
      !                                                      ! ========================= !
      IF( srcv(jpr_bhd)%laction .AND. ln_bern_srfc )   bhd_wave(:,:) = frcv(jpr_bhd)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !   Stokes transport u dir  !
      !                                                      ! ========================= !
      IF( srcv(jpr_tusd)%laction .AND. ln_breivikFV_2016 )    tusd(:,:) = frcv(jpr_tusd)%z3(:,:,1)
      !
      !                                                      ! ========================= !
      !                                                      !   Stokes transport v dir  !
      !                                                      ! ========================= !
      IF( srcv(jpr_tvsd)%laction .AND. ln_breivikFV_2016 )     tvsd(:,:) = frcv(jpr_tvsd)%z3(:,:,1)
      !
      !  Fields received by SAS when OASIS coupling
      !  (arrays no more filled at sbcssm stage)
      !                                                      ! ================== !
      !                                                      !        SSS         !
      !                                                      ! ================== !
      IF( srcv(jpr_soce)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         sss_m(:,:) = frcv(jpr_soce)%z3(:,:,1)
         CALL iom_put( 'sss_m', sss_m )
      ENDIF
      !
      !                                                      ! ================== !
      !                                                      !        SST         !
      !                                                      ! ================== !
      IF( srcv(jpr_toce)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         sst_m(:,:) = frcv(jpr_toce)%z3(:,:,1)
         IF( srcv(jpr_soce)%laction .AND. l_useCT ) THEN    ! make sure that sst_m is the potential temperature
            sst_m(:,:) = eos_pt_from_ct( sst_m(:,:), sss_m(:,:) )
         ENDIF
      ENDIF
      !                                                      ! ================== !
      !                                                      !        SSH         !
      !                                                      ! ================== !
      IF( srcv(jpr_ssh )%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         ssh_m(:,:) = frcv(jpr_ssh )%z3(:,:,1)
         CALL iom_put( 'ssh_m', ssh_m )
      ENDIF
      !                                                      ! ================== !
      !                                                      !  surface currents  !
      !                                                      ! ================== !
      IF( srcv(jpr_ocx1)%laction ) THEN                      ! received by sas in case of opa <-> sas coupling
         ssu_m(:,:) = frcv(jpr_ocx1)%z3(:,:,1)
         uu(:,:,1,Kbb) = ssu_m(:,:)                          ! will be used in icestp in the call of ice_forcing_tau
         uu(:,:,1,Kmm) = ssu_m(:,:)                          ! will be used in sbc_cpl_snd if atmosphere coupling
         CALL iom_put( 'ssu_m', ssu_m )
      ENDIF
      IF( srcv(jpr_ocy1)%laction ) THEN
         ssv_m(:,:) = frcv(jpr_ocy1)%z3(:,:,1)
         vv(:,:,1,Kbb) = ssv_m(:,:)                          ! will be used in icestp in the call of ice_forcing_tau
         vv(:,:,1,Kmm) = ssv_m(:,:)                          ! will be used in sbc_cpl_snd if atmosphere coupling
         CALL iom_put( 'ssv_m', ssv_m )
      ENDIF
      !                                                      ! ======================== !
      !                                                      !  first T level thickness !
      !                                                      ! ======================== !
      IF( srcv(jpr_e3t1st )%laction ) THEN                   ! received by sas in case of opa <-> sas coupling
         e3t_m(:,:) = frcv(jpr_e3t1st )%z3(:,:,1)
         CALL iom_put( 'e3t_m', e3t_m(:,:) )
      ENDIF
      !                                                      ! ================================ !
      !                                                      !  fraction of solar net radiation !
      !                                                      ! ================================ !
      IF( srcv(jpr_fraqsr)%laction ) THEN                    ! received by sas in case of opa <-> sas coupling
         frq_m(:,:) = frcv(jpr_fraqsr)%z3(:,:,1)
         CALL iom_put( 'frq_m', frq_m )
      ENDIF

      !                                                      ! ========================= !
      IF( k_ice <= 1 .AND. MOD( kt-1, k_fsbc ) == 0 ) THEN   !  heat & freshwater fluxes ! (Ocean only case)
         !                                                   ! ========================= !
         !
         !                                                       ! total freshwater fluxes over the ocean (emp)
         IF( srcv(jpr_oemp)%laction .OR. srcv(jpr_rain)%laction ) THEN
            SELECT CASE( TRIM( sn_rcv_emp%cldes ) )                                    ! evaporation - precipitation
            CASE( 'conservative' )
               zemp(:,:) = frcv(jpr_tevp)%z3(:,:,1) - ( frcv(jpr_rain)%z3(:,:,1) + frcv(jpr_snow)%z3(:,:,1) )
            CASE( 'oce only', 'oce and ice' )
               zemp(:,:) = frcv(jpr_oemp)%z3(:,:,1)
            CASE default
               CALL ctl_stop( 'sbc_cpl_rcv: wrong definition of sn_rcv_emp%cldes' )
            END SELECT
         ELSE
            zemp(:,:) = 0._wp
         ENDIF
         !
         !                                                        ! runoffs and calving (added in emp)
         IF( srcv(jpr_rnf)%laction )     rnf(:,:) = frcv(jpr_rnf)%z3(:,:,1)
         IF( srcv(jpr_cal)%laction )     zemp(:,:) = zemp(:,:) - frcv(jpr_cal)%z3(:,:,1)

         IF( srcv(jpr_icb)%laction )  THEN
             fwficb(:,:) = frcv(jpr_icb)%z3(:,:,1)
             rnf(:,:)    = rnf(:,:) + fwficb(:,:)   ! iceberg added to runfofs
         ENDIF
         !
         ! ice shelf fwf
         IF( srcv(jpr_isf)%laction )  THEN
            fwfisf_oasis(:,:) = frcv(jpr_isf)%z3(:,:,1)  ! fresh water flux from the isf to the ocean ( > 0 = melting )
         END IF

         IF( ln_mixcpl ) THEN   ;   emp(:,:) = emp(:,:) * xcplmask(:,:,0) + zemp(:,:) * zmsk(:,:)
         ELSE                   ;   emp(:,:) =                              zemp(:,:)
         ENDIF
         !
         !                                                       ! non solar heat flux over the ocean (qns)
         IF(      srcv(jpr_qnsoce)%laction ) THEN   ;   zqns(:,:) = frcv(jpr_qnsoce)%z3(:,:,1)
         ELSE IF( srcv(jpr_qnsmix)%laction ) THEN   ;   zqns(:,:) = frcv(jpr_qnsmix)%z3(:,:,1)
         ELSE                                       ;   zqns(:,:) = 0._wp
         ENDIF
         ! update qns over the free ocean with:
         IF( nn_components /= jp_iam_oce ) THEN
            zqns(:,:) =  zqns(:,:) - zemp(:,:) * sst_m(:,:) * rcp         ! remove heat content due to mass flux (assumed to be at SST)
            IF( srcv(jpr_snow  )%laction ) THEN
               zqns(:,:) = zqns(:,:) - frcv(jpr_snow)%z3(:,:,1) * rLfus   ! energy for melting solid precipitation over the free ocean
            ENDIF
         ENDIF
         !
         IF( srcv(jpr_icb)%laction )  zqns(:,:) = zqns(:,:) - frcv(jpr_icb)%z3(:,:,1) * rLfus ! remove heat content associated to iceberg melting
         !
         IF( ln_mixcpl ) THEN   ;   qns(:,:) = qns(:,:) * xcplmask(:,:,0) + zqns(:,:) * zmsk(:,:)
         ELSE                   ;   qns(:,:) =                              zqns(:,:)
         ENDIF

         !                                                       ! solar flux over the ocean          (qsr)
         IF     ( srcv(jpr_qsroce)%laction ) THEN   ;   zqsr(:,:) = frcv(jpr_qsroce)%z3(:,:,1)
         ELSE IF( srcv(jpr_qsrmix)%laction ) then   ;   zqsr(:,:) = frcv(jpr_qsrmix)%z3(:,:,1)
         ELSE                                       ;   zqsr(:,:) = 0._wp
         ENDIF
         IF( ln_dm2dc .AND. ln_cpl )   zqsr(:,:) = sbc_dcy( zqsr )   ! modify qsr to include the diurnal cycle
         IF( ln_mixcpl ) THEN   ;   qsr(:,:) = qsr(:,:) * xcplmask(:,:,0) + zqsr(:,:) * zmsk(:,:)
         ELSE                   ;   qsr(:,:) =                              zqsr(:,:)
         ENDIF
         !
         ! salt flux over the ocean (received by opa in case of opa <-> sas coupling)
         IF( srcv(jpr_sflx )%laction )   sfx(:,:) = frcv(jpr_sflx  )%z3(:,:,1)
         ! Ice cover  (received by opa in case of opa <-> sas coupling)
         IF( srcv(jpr_fice )%laction )   fr_i(:,:) = frcv(jpr_fice )%z3(:,:,1)
         !
      ENDIF
      !
   END SUBROUTINE sbc_cpl_rcv


   SUBROUTINE sbc_cpl_ice_tau( p_taui, p_tauj )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_ice_tau  ***
      !!
      !! ** Purpose :   provide the stress over sea-ice in coupled mode
      !!
      !! ** Method  :   transform the received stress from the atmosphere into
      !!             an atmosphere-ice stress in the (i,j) ocean referencial
      !!             and at the velocity point of the sea-ice model:
      !!                'C'-grid : i- (j-) components given at U- (V-) point
      !!
      !!                The received stress are :
      !!                 - defined by 3 components (if cartesian coordinate)
      !!                        or by 2 components (if spherical)
      !!                 - oriented along geographical   coordinate (if eastward-northward)
      !!                        or  along the local grid coordinate (if local grid)
      !!                 - given at U- and V-point, resp.   if received on 2 grids
      !!                        or at a same point (T or I) if received on 1 grid
      !!                Therefore and if necessary, they are successively
      !!             processed in order to obtain them
      !!                 first  as  2 components on the sphere
      !!                 second as  2 components oriented along the local grid
      !!                 third  as  2 components on the ice grid point
      !!
      !!                Except in 'oce and ice' case, only one vector stress field
      !!             is received. It has already been processed in sbc_cpl_rcv
      !!             so that it is now defined as (i,j) components given at U-
      !!             and V-points, respectively.
      !!
      !! ** Action  :   return ptau_i, ptau_j, the stress over the ice
      !!----------------------------------------------------------------------
      REAL(wp), INTENT(inout), DIMENSION(:,:) ::   p_taui   ! i- & j-components of atmos-ice stress [N/m2]
      REAL(wp), INTENT(inout), DIMENSION(:,:) ::   p_tauj   ! at I-point (B-grid) or U & V-point (C-grid)
      !!
      INTEGER ::   ji, jj   ! dummy loop indices
      INTEGER ::   itx      ! index of taux over ice
      REAL(wp)                     ::   zztmp1, zztmp2
      REAL(wp), DIMENSION(jpi,jpj) ::   ztx, zty
      !!----------------------------------------------------------------------
      !
      !
   END SUBROUTINE sbc_cpl_ice_tau


   SUBROUTINE sbc_cpl_ice_flx( kt, picefr, palbi, psst, pist, phs, phi )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_ice_flx  ***
      !!
      !! ** Purpose :   provide the heat and freshwater fluxes of the ocean-ice system
      !!
      !! ** Method  :   transform the fields received from the atmosphere into
      !!             surface heat and fresh water boundary condition for the
      !!             ice-ocean system. The following fields are provided:
      !!               * total non solar, solar and freshwater fluxes (qns_tot,
      !!             qsr_tot and emp_tot) (total means weighted ice-ocean flux)
      !!             NB: emp_tot include runoffs and calving.
      !!               * fluxes over ice (qns_ice, qsr_ice, emp_ice) where
      !!             emp_ice = sublimation - solid precipitation as liquid
      !!             precipitation are re-routed directly to the ocean and
      !!             calving directly enter the ocean (runoffs are read but included in trasbc.F90)
      !!               * solid precipitation (sprecip), used to add to qns_tot
      !!             the heat lost associated to melting solid precipitation
      !!             over the ocean fraction.
      !!               * heat content of rain, snow and evap can also be provided,
      !!             otherwise heat flux associated with these mass flux are
      !!             guessed (qemp_oce, qemp_ice)
      !!
      !!             - the fluxes have been separated from the stress as
      !!               (a) they are updated at each ice time step compare to
      !!               an update at each coupled time step for the stress, and
      !!               (b) the conservative computation of the fluxes over the
      !!               sea-ice area requires the knowledge of the ice fraction
      !!               after the ice advection and before the ice thermodynamics,
      !!               so that the stress is updated before the ice dynamics
      !!               while the fluxes are updated after it.
      !!
      !! ** Details
      !!             qns_tot = (1-a) * qns_oce + a * qns_ice               => provided
      !!                     + qemp_oce + qemp_ice                         => recalculated and added up to qns
      !!
      !!             qsr_tot = (1-a) * qsr_oce + a * qsr_ice               => provided
      !!
      !!             emp_tot = emp_oce + emp_ice                           => calving is provided and added to emp_tot (and emp_oce).
      !!                                                                      runoff (which includes rivers+icebergs) and iceshelf
      !!                                                                      are provided but not included in emp here. Only runoff will
      !!                                                                      be included in emp in other parts of NEMO code
      !!
      !! ** Note : In case of the ice-atm coupling with conduction fluxes (such as Jules interface for the Met-Office),
      !!              qsr_ice and qns_ice are not provided and they are not supposed to be used in the ice code.
      !!              However, by precaution we also "fake" qns_ice and qsr_ice this way:
      !!              qns_ice = qml_ice + qcn_ice ??
      !!              qsr_ice = qtr_ice_top ??
      !!
      !! ** Action  :   update at each nf_ice time step:
      !!                   qns_tot, qsr_tot  non-solar and solar total heat fluxes
      !!                   qns_ice, qsr_ice  non-solar and solar heat fluxes over the ice
      !!                   emp_tot           total evaporation - precipitation(liquid and solid) (-calving)
      !!                   emp_ice           ice sublimation - solid precipitation over the ice
      !!                   dqns_ice          d(non-solar heat flux)/d(Temperature) over the ice
      !!                   sprecip           solid precipitation over the ocean
      !!----------------------------------------------------------------------
      INTEGER,  INTENT(in)                                ::   kt         ! ocean model time step index (only for a_i_last_couple)
      REAL(wp), INTENT(in)   , DIMENSION(:,:)             ::   picefr     ! ice fraction                [0 to 1]
      !                                                   !!           ! optional arguments, used only in 'mixed oce-ice' case or for Met-Office coupling
      REAL(wp), INTENT(in)   , DIMENSION(:,:,:), OPTIONAL ::   palbi      ! all skies ice albedo
      REAL(wp), INTENT(in)   , DIMENSION(:,:  ), OPTIONAL ::   psst       ! sea surface temperature     [Celsius]
      REAL(wp), INTENT(inout), DIMENSION(:,:,:), OPTIONAL ::   pist       ! ice surface temperature     [Kelvin] => inout for Met-Office
      REAL(wp), INTENT(in)   , DIMENSION(:,:,:), OPTIONAL ::   phs        ! snow depth                  [m]
      REAL(wp), INTENT(in)   , DIMENSION(:,:,:), OPTIONAL ::   phi        ! ice thickness               [m]
      !
      INTEGER  ::   ji, jj, jl   ! dummy loop index
      REAL(wp), DIMENSION(jpi,jpj)     ::   zcptn, zcptrain, zcptsnw, ziceld, zmsk, zsnw
      REAL(wp), DIMENSION(jpi,jpj)     ::   zemp_tot, zemp_ice, zemp_oce, ztprecip, zsprecip  , zevap_oce, zdevap_ice
      REAL(wp), DIMENSION(jpi,jpj)     ::   zqns_tot, zqns_oce, zqsr_tot, zqsr_oce, zqprec_ice, zqemp_oce, zqemp_ice
      REAL(wp), DIMENSION(jpi,jpj)     ::   zevap_ice_total
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   zqns_ice, zqsr_ice, zdqns_ice, zqevap_ice, zevap_ice, zqtr_ice_top, ztsu
      REAL(wp), DIMENSION(jpi,jpj)     ::   ztri
      !!----------------------------------------------------------------------
      !
      !
   END SUBROUTINE sbc_cpl_ice_flx


   SUBROUTINE sbc_cpl_snd( kt, Kbb, Kmm )
      !!----------------------------------------------------------------------
      !!             ***  ROUTINE sbc_cpl_snd  ***
      !!
      !! ** Purpose :   provide the ocean-ice informations to the atmosphere
      !!
      !! ** Method  :   send to the atmosphere through a call to cpl_snd
      !!              all the needed fields (as defined in sbc_cpl_init)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt
      INTEGER, INTENT(in) ::   Kbb, Kmm    ! ocean model time level index
      !
      INTEGER ::   ji, jj, jl   ! dummy loop indices
      INTEGER ::   isec, info   ! local integer
      REAL(wp) ::   zumax, zvmax
      REAL(wp), DIMENSION(jpi,jpj)     ::   zfr_l, ztmp1, ztmp2, zotx1, zoty1, zotz1, zitx1, zity1, zitz1
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   ztmp3, ztmp4
      !!----------------------------------------------------------------------
      !
      isec = ( kt - nit000 ) * NINT( rn_Dt )        ! date of exchanges
      info = OASIS_idle

      zfr_l(:,:) = 1.- fr_i(:,:)
      !                                                      ! ------------------------- !
      !                                                      !    Surface temperature    !   in Kelvin
      !                                                      ! ------------------------- !
      IF( ssnd(jps_toce)%laction .OR. ssnd(jps_tice)%laction .OR. ssnd(jps_tmix)%laction ) THEN

         IF( nn_components == jp_iam_oce ) THEN
            ztmp1(:,:) = ts(:,:,1,jp_tem,Kmm)   ! send temperature as it is (potential or conservative) -> use of l_useCT on the received part
         ELSE
            ! we must send the surface potential temperature
            IF( l_useCT )  THEN    ;   ztmp1(:,:) =eos_pt_from_ct( ts(:,:,1,jp_tem,Kmm), ts(:,:,1,jp_sal,Kmm) )
            ELSE                   ;   ztmp1(:,:) = ts(:,:,1,jp_tem,Kmm)
            ENDIF
            !
            SELECT CASE( sn_snd_temp%cldes)
            CASE( 'oce only'             )   ;   ztmp1(:,:) =   ztmp1(:,:) + rt0
            CASE( 'oce and ice'          )   ;   ztmp1(:,:) =   ztmp1(:,:) + rt0
               SELECT CASE( sn_snd_temp%clcat )
               CASE( 'yes' )
                  ztmp3(:,:,1:jpl) = tn_ice(:,:,1:jpl)
               CASE( 'no' )
                  WHERE( SUM( a_i, dim=3 ) /= 0. )
                     ztmp3(:,:,1) = SUM( tn_ice * a_i, dim=3 ) / SUM( a_i, dim=3 )
                  ELSEWHERE
                     ztmp3(:,:,1) = rt0
                  END WHERE
               CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%clcat' )
               END SELECT
            CASE( 'weighted oce and ice' )   ;   ztmp1(:,:) = ( ztmp1(:,:) + rt0 ) * zfr_l(:,:)
               SELECT CASE( sn_snd_temp%clcat )
               CASE( 'yes' )
                  ztmp3(:,:,1:jpl) = tn_ice(:,:,1:jpl) * a_i(:,:,1:jpl)
               CASE( 'no' )
                  ztmp3(:,:,:) = 0.0
                  DO jl=1,jpl
                     ztmp3(:,:,1) = ztmp3(:,:,1) + tn_ice(:,:,jl) * a_i(:,:,jl)
                  ENDDO
               CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%clcat' )
               END SELECT
            CASE( 'oce and weighted ice')    ;   ztmp1(:,:) =   ts(:,:,1,jp_tem,Kmm) + rt0
               SELECT CASE( sn_snd_temp%clcat )
               CASE( 'yes' )
                  ztmp3(:,:,1:jpl) = tn_ice(:,:,1:jpl) * a_i(:,:,1:jpl)
               CASE( 'no' )
                  ztmp3(:,:,:) = 0.0
                  DO jl=1,jpl
                     ztmp3(:,:,1) = ztmp3(:,:,1) + tn_ice(:,:,jl) * a_i(:,:,jl)
                  ENDDO
               CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%clcat' )
               END SELECT
            CASE( 'mixed oce-ice'        )
               ztmp1(:,:) = ( ztmp1(:,:) + rt0 ) * zfr_l(:,:)
               DO jl=1,jpl
                  ztmp1(:,:) = ztmp1(:,:) + tn_ice(:,:,jl) * a_i(:,:,jl)
               ENDDO
            CASE default                     ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_temp%cldes' )
            END SELECT
         ENDIF
         IF( ssnd(jps_toce)%laction )   CALL cpl_snd( jps_toce, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info )
         IF( ssnd(jps_tice)%laction )   CALL cpl_snd( jps_tice, isec, ztmp3, info )
         IF( ssnd(jps_tmix)%laction )   CALL cpl_snd( jps_tmix, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info )
      ENDIF
      !
      !                                                      ! ------------------------- !
      !                                                      ! 1st layer ice/snow temp.  !
      !                                                      ! ------------------------- !
      !                                                      ! ------------------------- !
      !                                                      !           Albedo          !
      !                                                      ! ------------------------- !
      IF( ssnd(jps_albice)%laction ) THEN                         ! ice
          SELECT CASE( sn_snd_alb%cldes )
          CASE( 'ice' )
             SELECT CASE( sn_snd_alb%clcat )
             CASE( 'yes' )
                ztmp3(:,:,1:jpl) = alb_ice(:,:,1:jpl)
             CASE( 'no' )
                WHERE( SUM( a_i, dim=3 ) /= 0. )
                   ztmp1(:,:) = SUM( alb_ice (:,:,1:jpl) * a_i(:,:,1:jpl), dim=3 ) / SUM( a_i(:,:,1:jpl), dim=3 )
                ELSEWHERE
                   ztmp1(:,:) = alb_oce_mix(:,:)
                END WHERE
             CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_alb%clcat' )
             END SELECT
          CASE( 'weighted ice' )   ;
             SELECT CASE( sn_snd_alb%clcat )
             CASE( 'yes' )
                ztmp3(:,:,1:jpl) =  alb_ice(:,:,1:jpl) * a_i(:,:,1:jpl)
             CASE( 'no' )
                WHERE( fr_i (:,:) > 0. )
                   ztmp1(:,:) = SUM (  alb_ice(:,:,1:jpl) * a_i(:,:,1:jpl), dim=3 )
                ELSEWHERE
                   ztmp1(:,:) = 0.
                END WHERE
             CASE default   ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_ice%clcat' )
             END SELECT
          CASE default      ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_alb%cldes' )
         END SELECT

         SELECT CASE( sn_snd_alb%clcat )
            CASE( 'yes' )
               CALL cpl_snd( jps_albice, isec, ztmp3, info )      !-> MV this has never been checked in coupled mode
            CASE( 'no'  )
               CALL cpl_snd( jps_albice, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info )
         END SELECT
      ENDIF

      IF( ssnd(jps_albmix)%laction ) THEN                         ! mixed ice-ocean
         ztmp1(:,:) = alb_oce_mix(:,:) * zfr_l(:,:)
         DO jl = 1, jpl
            ztmp1(:,:) = ztmp1(:,:) + alb_ice(:,:,jl) * a_i(:,:,jl)
         END DO
         CALL cpl_snd( jps_albmix, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info )
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !  Ice fraction & Thickness !
      !                                                      ! ------------------------- !
      ! Send ice fraction field to atmosphere
      IF( ssnd(jps_fice)%laction ) THEN
         SELECT CASE( sn_snd_thick%clcat )
         CASE( 'yes' )   ;   ztmp3(:,:,1:jpl) =  a_i(:,:,1:jpl)
         CASE( 'no'  )   ;   ztmp3(:,:,1    ) = fr_i(:,:      )
         CASE default    ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
         END SELECT
         CALL cpl_snd( jps_fice, isec, ztmp3, info )
      ENDIF


      IF( ssnd(jps_fice1)%laction ) THEN
         SELECT CASE( sn_snd_thick1%clcat )
         CASE( 'yes' )   ;   ztmp3(:,:,1:jpl) =  a_i(:,:,1:jpl)
         CASE( 'no'  )   ;   ztmp3(:,:,1    ) = fr_i(:,:      )
         CASE default    ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick1%clcat' )
         END SELECT
         CALL cpl_snd( jps_fice1, isec, ztmp3, info )
      ENDIF

      ! Send ice fraction field to OCE (sent by SAS in SAS-OCE coupling)
      IF( ssnd(jps_fice2)%laction ) THEN
         ztmp3(:,:,1) = fr_i(:,:)
         IF( ssnd(jps_fice2)%laction )   CALL cpl_snd( jps_fice2, isec, ztmp3, info )
      ENDIF

      ! Send ice and snow thickness field
      IF( ssnd(jps_hice)%laction .OR. ssnd(jps_hsnw)%laction ) THEN
         SELECT CASE( sn_snd_thick%cldes)
         CASE( 'none'                  )       ! nothing to do
         CASE( 'weighted ice and snow' )
            SELECT CASE( sn_snd_thick%clcat )
            CASE( 'yes' )
               ztmp3(:,:,1:jpl) =  h_i(:,:,1:jpl) * a_i(:,:,1:jpl)
               ztmp4(:,:,1:jpl) =  h_s(:,:,1:jpl) * a_i(:,:,1:jpl)
            CASE( 'no' )
               ztmp3(:,:,:) = 0.0   ;  ztmp4(:,:,:) = 0.0
               DO jl=1,jpl
                  ztmp3(:,:,1) = ztmp3(:,:,1) + h_i(:,:,jl) * a_i(:,:,jl)
                  ztmp4(:,:,1) = ztmp4(:,:,1) + h_s(:,:,jl) * a_i(:,:,jl)
               ENDDO
            CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
            END SELECT
         CASE( 'ice and snow'         )
            SELECT CASE( sn_snd_thick%clcat )
            CASE( 'yes' )
               ztmp3(:,:,1:jpl) = h_i(:,:,1:jpl)
               ztmp4(:,:,1:jpl) = h_s(:,:,1:jpl)
            CASE( 'no' )
               WHERE( SUM( a_i, dim=3 ) /= 0. )
                  ztmp3(:,:,1) = SUM( h_i * a_i, dim=3 ) / SUM( a_i, dim=3 )
                  ztmp4(:,:,1) = SUM( h_s * a_i, dim=3 ) / SUM( a_i, dim=3 )
               ELSEWHERE
                 ztmp3(:,:,1) = 0.
                 ztmp4(:,:,1) = 0.
               END WHERE
            CASE default                  ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%clcat' )
            END SELECT
         CASE default                     ;   CALL ctl_stop( 'sbc_cpl_snd: wrong definition of sn_snd_thick%cldes' )
         END SELECT
         IF( ssnd(jps_hice)%laction )   CALL cpl_snd( jps_hice, isec, ztmp3, info )
         IF( ssnd(jps_hsnw)%laction )   CALL cpl_snd( jps_hsnw, isec, ztmp4, info )
      ENDIF


      !                                                      ! ------------------------- !
      !                                                      !  CO2 flux from PISCES     !
      !                                                      ! ------------------------- !
      IF( ssnd(jps_co2)%laction .AND. l_co2cpl )   THEN
         ztmp1(:,:) = oce_co2(:,:) * 1000.  ! conversion in molC/m2/s
         CALL cpl_snd( jps_co2, isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ) , info )
      ENDIF
      !
      !                                                      ! ------------------------- !
      IF( ssnd(jps_ocx1)%laction ) THEN                      !      Surface current      !
         !                                                   ! ------------------------- !
         !
         !                                                  j+1   j     -----V---F
         ! surface velocity always sent from T point                     !       |
         !                                                        j      |   T   U
         !                                                               |       |
         !                                                   j    j-1   -I-------|
         !                                               (for I)         |       |
         !                                                              i-1  i   i
         !                                                               i      i+1 (for I)
         IF( nn_components == jp_iam_oce ) THEN
            zotx1(:,:) = uu(:,:,1,Kmm)
            zoty1(:,:) = vv(:,:,1,Kmm)
         ELSE
            SELECT CASE( TRIM( sn_snd_crt%cldes ) )
            CASE( 'oce only'             )      ! C-grid ==> T
               DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
                  zotx1(ji,jj) = 0.5 * ( uu(ji,jj,1,Kmm) + uu(ji-1,jj  ,1,Kmm) )
                  zoty1(ji,jj) = 0.5 * ( vv(ji,jj,1,Kmm) + vv(ji  ,jj-1,1,Kmm) )
               END DO   ;   END DO
            CASE( 'weighted oce and ice' )      ! Ocean and Ice on C-grid ==> T
               DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
                  zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) ) * zfr_l(ji,jj)
                  zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) ) * zfr_l(ji,jj)
                  zitx1(ji,jj) = 0.5 * ( u_ice(ji,jj  )     + u_ice(ji-1,jj    )     ) *  fr_i(ji,jj)
                  zity1(ji,jj) = 0.5 * ( v_ice(ji,jj  )     + v_ice(ji  ,jj-1  )     ) *  fr_i(ji,jj)
               END DO   ;   END DO
               CALL lbc_lnk( 'sbccpl', zitx1, 'T', -1.0_wp, zity1, 'T', -1.0_wp )
            CASE( 'mixed oce-ice'        )      ! Ocean and Ice on C-grid ==> T
               DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
                  zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) ) * zfr_l(ji,jj)   &
                     &         + 0.5 * ( u_ice(ji,jj  )     + u_ice(ji-1,jj    )     ) *  fr_i(ji,jj)
                  zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) ) * zfr_l(ji,jj)   &
                     &         + 0.5 * ( v_ice(ji,jj  )     + v_ice(ji  ,jj-1  )     ) *  fr_i(ji,jj)
               END DO   ;   END DO
            END SELECT
            CALL lbc_lnk( 'sbccpl', zotx1, ssnd(jps_ocx1)%clgrid, -1.0_wp,  zoty1, ssnd(jps_ocy1)%clgrid, -1.0_wp )
            !
         ENDIF
         !
         !
         IF( TRIM( sn_snd_crt%clvor ) == 'eastward-northward' ) THEN             ! Rotation of the components
            !                                                                     ! Ocean component
            CALL rot_rep( zotx1, zoty1, ssnd(jps_ocx1)%clgrid, 'ij->e', ztmp1 )       ! 1st component
            CALL rot_rep( zotx1, zoty1, ssnd(jps_ocx1)%clgrid, 'ij->n', ztmp2 )       ! 2nd component
            zotx1(:,:) = ztmp1(:,:)                                                   ! overwrite the components
            zoty1(:,:) = ztmp2(:,:)
            IF( ssnd(jps_ivx1)%laction ) THEN                                     ! Ice component
               CALL rot_rep( zitx1, zity1, ssnd(jps_ivx1)%clgrid, 'ij->e', ztmp1 )    ! 1st component
               CALL rot_rep( zitx1, zity1, ssnd(jps_ivx1)%clgrid, 'ij->n', ztmp2 )    ! 2nd component
               zitx1(:,:) = ztmp1(:,:)                                                ! overwrite the components
               zity1(:,:) = ztmp2(:,:)
            ENDIF
         ENDIF
         !
         ! spherical coordinates to cartesian -> 2 components to 3 components
         IF( TRIM( sn_snd_crt%clvref ) == 'cartesian' ) THEN
            ztmp1(:,:) = zotx1(:,:)                     ! ocean currents
            ztmp2(:,:) = zoty1(:,:)
            CALL oce2geo ( ztmp1, ztmp2, 'T', zotx1, zoty1, zotz1 )
            !
            IF( ssnd(jps_ivx1)%laction ) THEN           ! ice velocities
               ztmp1(:,:) = zitx1(:,:)
               ztmp1(:,:) = zity1(:,:)
               CALL oce2geo ( ztmp1, ztmp2, 'T', zitx1, zity1, zitz1 )
            ENDIF
         ENDIF
         !
         IF( ssnd(jps_ocx1)%laction )   CALL cpl_snd( jps_ocx1, isec, RESHAPE ( zotx1, (/jpi,jpj,1/) ), info )   ! ocean x current 1st grid
         IF( ssnd(jps_ocy1)%laction )   CALL cpl_snd( jps_ocy1, isec, RESHAPE ( zoty1, (/jpi,jpj,1/) ), info )   ! ocean y current 1st grid
         IF( ssnd(jps_ocz1)%laction )   CALL cpl_snd( jps_ocz1, isec, RESHAPE ( zotz1, (/jpi,jpj,1/) ), info )   ! ocean z current 1st grid
         !
         IF( ssnd(jps_ivx1)%laction )   CALL cpl_snd( jps_ivx1, isec, RESHAPE ( zitx1, (/jpi,jpj,1/) ), info )   ! ice   x current 1st grid
         IF( ssnd(jps_ivy1)%laction )   CALL cpl_snd( jps_ivy1, isec, RESHAPE ( zity1, (/jpi,jpj,1/) ), info )   ! ice   y current 1st grid
         IF( ssnd(jps_ivz1)%laction )   CALL cpl_snd( jps_ivz1, isec, RESHAPE ( zitz1, (/jpi,jpj,1/) ), info )   ! ice   z current 1st grid
         !
      ENDIF
      !
      !                                                      ! ------------------------- !
      !                                                      !  Surface current to waves !
      !                                                      ! ------------------------- !
      IF( ssnd(jps_ocxw)%laction .OR. ssnd(jps_ocyw)%laction ) THEN
          !
          !                                                  j+1  j     -----V---F
          ! surface velocity always sent from T point                    !       |
          !                                                       j      |   T   U
          !                                                              |       |
          !                                                   j   j-1   -I-------|
          !                                               (for I)        |       |
          !                                                             i-1  i   i
          !                                                              i      i+1 (for I)
          SELECT CASE( TRIM( sn_snd_crtw%cldes ) )
          CASE( 'oce only'             )      ! C-grid ==> T
             DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
                zotx1(ji,jj) = 0.5 * ( uu(ji,jj,1,Kmm) + uu(ji-1,jj  ,1,Kmm) )
                zoty1(ji,jj) = 0.5 * ( vv(ji,jj,1,Kmm) + vv(ji , jj-1,1,Kmm) )
             END DO   ;   END DO
          CASE( 'weighted oce and ice' )      ! Ocean and Ice on C-grid ==> T
             DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
                zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) ) * zfr_l(ji,jj)
                zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) ) * zfr_l(ji,jj)
                zitx1(ji,jj) = 0.5 * ( u_ice(ji,jj  ) + u_ice(ji-1,jj    ) ) *  fr_i(ji,jj)
                zity1(ji,jj) = 0.5 * ( v_ice(ji,jj  ) + v_ice(ji  ,jj-1  ) ) *  fr_i(ji,jj)
             END DO   ;   END DO
             CALL lbc_lnk( 'sbccpl', zitx1, 'T', -1.0_wp,  zity1, 'T', -1.0_wp )
          CASE( 'mixed oce-ice'        )      ! Ocean and Ice on C-grid ==> T
             DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
                zotx1(ji,jj) = 0.5 * ( uu   (ji,jj,1,Kmm) + uu   (ji-1,jj  ,1,Kmm) ) * zfr_l(ji,jj)   &
                   &         + 0.5 * ( u_ice(ji,jj  ) + u_ice(ji-1,jj    ) ) *  fr_i(ji,jj)
                zoty1(ji,jj) = 0.5 * ( vv   (ji,jj,1,Kmm) + vv   (ji  ,jj-1,1,Kmm) ) * zfr_l(ji,jj)   &
                   &         + 0.5 * ( v_ice(ji,jj  ) + v_ice(ji  ,jj-1  ) ) *  fr_i(ji,jj)
             END DO   ;   END DO
          END SELECT
         CALL lbc_lnk( 'sbccpl', zotx1, ssnd(jps_ocxw)%clgrid, -1.0_wp, zoty1, ssnd(jps_ocyw)%clgrid, -1.0_wp )
         !
         !
         IF( TRIM( sn_snd_crtw%clvor ) == 'eastward-northward' ) THEN             ! Rotation of the components
         !                                                                        ! Ocean component
            CALL rot_rep( zotx1, zoty1, ssnd(jps_ocxw)%clgrid, 'ij->e', ztmp1 )       ! 1st component
            CALL rot_rep( zotx1, zoty1, ssnd(jps_ocxw)%clgrid, 'ij->n', ztmp2 )       ! 2nd component
            zotx1(:,:) = ztmp1(:,:)                                                   ! overwrite the components
            zoty1(:,:) = ztmp2(:,:)
            IF( ssnd(jps_ivx1)%laction ) THEN                                     ! Ice component
               CALL rot_rep( zitx1, zity1, ssnd(jps_ivx1)%clgrid, 'ij->e', ztmp1 )    ! 1st component
               CALL rot_rep( zitx1, zity1, ssnd(jps_ivx1)%clgrid, 'ij->n', ztmp2 )    ! 2nd component
               zitx1(:,:) = ztmp1(:,:)                                                ! overwrite the components
               zity1(:,:) = ztmp2(:,:)
            ENDIF
         ENDIF
         !
!         ! spherical coordinates to cartesian -> 2 components to 3 components
!         IF( TRIM( sn_snd_crtw%clvref ) == 'cartesian' ) THEN
!            ztmp1(:,:) = zotx1(:,:)                     ! ocean currents
!            ztmp2(:,:) = zoty1(:,:)
!            CALL oce2geo ( ztmp1, ztmp2, 'T', zotx1, zoty1, zotz1 )
!            !
!            IF( ssnd(jps_ivx1)%laction ) THEN           ! ice velocities
!               ztmp1(:,:) = zitx1(:,:)
!               ztmp1(:,:) = zity1(:,:)
!               CALL oce2geo ( ztmp1, ztmp2, 'T', zitx1, zity1, zitz1 )
!            ENDIF
!         ENDIF
         !
         IF( ssnd(jps_ocxw)%laction )   CALL cpl_snd( jps_ocxw, isec, RESHAPE ( zotx1, (/jpi,jpj,1/) ), info )   ! ocean x current 1st grid
         IF( ssnd(jps_ocyw)%laction )   CALL cpl_snd( jps_ocyw, isec, RESHAPE ( zoty1, (/jpi,jpj,1/) ), info )   ! ocean y current 1st grid
         !
      ENDIF
      !
      IF( ssnd(jps_ficet)%laction ) THEN
         CALL cpl_snd( jps_ficet, isec, RESHAPE ( fr_i, (/jpi,jpj,1/) ), info )
      ENDIF
      !                                                      ! ------------------------- !
      !                                                      !   Water levels to waves   !
      !                                                      ! ------------------------- !
      IF( ssnd(jps_wlev)%laction ) THEN
         IF( ln_apr_dyn ) THEN
            IF( kt /= nit000 ) THEN
               ztmp1(:,:) = ssh(:,:,Kbb) - 0.5 * ( ssh_ib(:,:) + ssh_ibb(:,:) )
            ELSE
               ztmp1(:,:) = ssh(:,:,Kbb)
            ENDIF
         ELSE
            ztmp1(:,:) = ssh(:,:,Kmm)
         ENDIF
         CALL cpl_snd( jps_wlev  , isec, RESHAPE ( ztmp1, (/jpi,jpj,1/) ), info )
      ENDIF
      !
      !  Fields sent by OCE to SAS when doing OCE<->SAS coupling
      !                                                        ! SSH
      IF( ssnd(jps_ssh )%laction )  THEN
         !                          ! removed inverse barometer ssh when Patm
         !                          forcing is used (for sea-ice dynamics)
         IF( ln_apr_dyn ) THEN   ;   ztmp1(:,:) = ssh(:,:,Kbb) - 0.5 * ( ssh_ib(:,:) + ssh_ibb(:,:) )
         ELSE                    ;   ztmp1(:,:) = ssh(:,:,Kmm)
         ENDIF
         CALL cpl_snd( jps_ssh   , isec, RESHAPE ( ztmp1            , (/jpi,jpj,1/) ), info )

      ENDIF
      !                                                        ! SSS
      IF( ssnd(jps_soce  )%laction )  THEN
         CALL cpl_snd( jps_soce  , isec, RESHAPE ( ts(:,:,1,jp_sal,Kmm), (/jpi,jpj,1/) ), info )
      ENDIF
      !                                                        ! first T level thickness
      IF( ssnd(jps_e3t1st )%laction )  THEN
 CALL cpl_snd( jps_e3t1st, isec, RESHAPE ( e3t_0(:,:,1) , (/jpi,jpj,1/) ), info )
      ENDIF
      !                                                        ! Qsr fraction
      IF( ssnd(jps_fraqsr)%laction )  THEN
         CALL cpl_snd( jps_fraqsr, isec, RESHAPE ( fraqsr_1lev(:,:) , (/jpi,jpj,1/) ), info )
      ENDIF
      !
      !  Fields sent by SAS to OCE when OASIS coupling
      !                                                        ! Solar heat flux
      IF( ssnd(jps_qsroce)%laction )  CALL cpl_snd( jps_qsroce, isec, RESHAPE ( qsr , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_qnsoce)%laction )  CALL cpl_snd( jps_qnsoce, isec, RESHAPE ( qns , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_oemp  )%laction )  CALL cpl_snd( jps_oemp  , isec, RESHAPE ( emp , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_sflx  )%laction )  CALL cpl_snd( jps_sflx  , isec, RESHAPE ( sfx , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_otx1  )%laction )  CALL cpl_snd( jps_otx1  , isec, RESHAPE ( utau, (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_oty1  )%laction )  CALL cpl_snd( jps_oty1  , isec, RESHAPE ( vtau, (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_rnf   )%laction )  CALL cpl_snd( jps_rnf   , isec, RESHAPE ( rnf , (/jpi,jpj,1/) ), info )
      IF( ssnd(jps_taum  )%laction )  CALL cpl_snd( jps_taum  , isec, RESHAPE ( taum, (/jpi,jpj,1/) ), info )

      !
   END SUBROUTINE sbc_cpl_snd

   !!======================================================================
END MODULE sbccpl

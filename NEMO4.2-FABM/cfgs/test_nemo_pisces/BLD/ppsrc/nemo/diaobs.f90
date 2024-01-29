










MODULE diaobs
   !!======================================================================
   !!                       ***  MODULE diaobs  ***
   !! Observation diagnostics: Computation of the misfit between data and
   !!                          their model equivalent 
   !!======================================================================
   !! History :  1.0  !  2006-03  (K. Mogensen) Original code
   !!             -   !  2006-05  (K. Mogensen, A. Weaver) Reformatted
   !!             -   !  2006-10  (A. Weaver) Cleaning and add controls
   !!             -   !  2007-03  (K. Mogensen) General handling of profiles
   !!             -   !  2007-04  (G. Smith) Generalized surface operators
   !!            2.0  !  2008-10  (M. Valdivieso) obs operator for velocity profiles
   !!            3.4  !  2014-08  (J. While) observation operator for profiles in all vertical coordinates
   !!             -   !                      Incorporated SST bias correction  
   !!            3.6  !  2015-02  (M. Martin) Simplification of namelist and code
   !!             -   !  2015-08  (M. Martin) Combined surface/profile routines.
   !!            4.0  !  2017-11  (G. Madec) style only
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_obs_init  : Reading and prepare observations
   !!   dia_obs       : Compute model equivalent to observations
   !!   dia_obs_wri   : Write observational diagnostics
   !!   calc_date     : Compute the date of timestep in YYYYMMDD.HHMMSS format
   !!   ini_date      : Compute the initial date YYYYMMDD.HHMMSS
   !!   fin_date      : Compute the final date YYYYMMDD.HHMMSS
   !!----------------------------------------------------------------------
   USE par_kind       ! Precision variables
   USE in_out_manager ! I/O manager
   USE par_oce        ! ocean parameter
   USE dom_oce        ! Ocean space and time domain variables
   USE sbc_oce        ! Sea-ice fraction
   !
   USE obs_read_prof  ! Reading and allocation of profile obs
   USE obs_read_surf  ! Reading and allocation of surface obs
   USE obs_sstbias    ! Bias correction routine for SST 
   USE obs_readmdt    ! Reading and allocation of MDT for SLA.
   USE obs_prep       ! Preparation of obs. (grid search etc).
   USE obs_oper       ! Observation operators
   USE obs_write      ! Writing of observation related diagnostics
   USE obs_grid       ! Grid searching
   USE obs_read_altbias ! Bias treatment for altimeter
   USE obs_profiles_def ! Profile data definitions
   USE obs_surf_def   ! Surface data definitions
   USE obs_types      ! Definitions for observation types
   !
   USE mpp_map        ! MPP mapping
   USE lib_mpp        ! For ctl_warn/stop

   IMPLICIT NONE
   PRIVATE

   PUBLIC dia_obs_init     ! Initialize and read observations
   PUBLIC dia_obs          ! Compute model equivalent to observations
   PUBLIC dia_obs_wri      ! Write model equivalent to observations
   PUBLIC dia_obs_dealloc  ! Deallocate dia_obs data
   PUBLIC calc_date        ! Compute the date of a timestep

   LOGICAL, PUBLIC :: ln_diaobs            !: Logical switch for the obs operator
   LOGICAL         :: ln_sstnight          !  Logical switch for night mean SST obs
   LOGICAL         :: ln_default_fp_indegs !  T=> Default obs footprint size specified in degrees, F=> in metres
   LOGICAL         :: ln_sla_fp_indegs     !  T=> SLA obs footprint size specified in degrees, F=> in metres
   LOGICAL         :: ln_sst_fp_indegs     !  T=> SST obs footprint size specified in degrees, F=> in metres
   LOGICAL         :: ln_sss_fp_indegs     !  T=> SSS obs footprint size specified in degrees, F=> in metres
   LOGICAL         :: ln_sic_fp_indegs     !  T=> sea-ice obs footprint size specified in degrees, F=> in metres

   REAL(wp) ::   rn_default_avglamscl      ! E/W diameter of SLA observation footprint (metres)
   REAL(wp) ::   rn_default_avgphiscl      ! N/S diameter of SLA observation footprint (metre
   REAL(wp) ::   rn_sla_avglamscl          ! E/W diameter of SLA observation footprint (metres)
   REAL(wp) ::   rn_sla_avgphiscl          ! N/S diameter of SLA observation footprint (metres)
   REAL(wp) ::   rn_sst_avglamscl          ! E/W diameter of SST observation footprint (metres)
   REAL(wp) ::   rn_sst_avgphiscl          ! N/S diameter of SST observation footprint (metres)
   REAL(wp) ::   rn_sss_avglamscl          ! E/W diameter of SSS observation footprint (metres)
   REAL(wp) ::   rn_sss_avgphiscl          ! N/S diameter of SSS observation footprint (metres)
   REAL(wp) ::   rn_sic_avglamscl          ! E/W diameter of sea-ice observation footprint (metres)
   REAL(wp) ::   rn_sic_avgphiscl          ! N/S diameter of sea-ice observation footprint (metres)

   INTEGER :: nn_1dint                     ! Vertical interpolation method
   INTEGER :: nn_2dint_default             ! Default horizontal interpolation method
   INTEGER :: nn_2dint_sla                 ! SLA horizontal interpolation method 
   INTEGER :: nn_2dint_sst                 ! SST horizontal interpolation method 
   INTEGER :: nn_2dint_sss                 ! SSS horizontal interpolation method 
   INTEGER :: nn_2dint_sic                 ! Seaice horizontal interpolation method 
   INTEGER, DIMENSION(imaxavtypes) ::   nn_profdavtypes   ! Profile data types representing a daily average
   INTEGER :: nproftypes     ! Number of profile obs types
   INTEGER :: nsurftypes     ! Number of surface obs types
   INTEGER , DIMENSION(:), ALLOCATABLE ::   nvarsprof, nvarssurf   ! Number of profile & surface variables
   INTEGER , DIMENSION(:), ALLOCATABLE ::   nextrprof, nextrsurf   ! Number of profile & surface extra variables
   INTEGER , DIMENSION(:), ALLOCATABLE ::   n2dintsurf             ! Interpolation option for surface variables
   REAL(wp), DIMENSION(:), ALLOCATABLE ::   zavglamscl, zavgphiscl ! E/W & N/S diameter of averaging footprint for surface variables
   LOGICAL , DIMENSION(:), ALLOCATABLE ::   lfpindegs              ! T=> surface obs footprint size specified in degrees, F=> in metres
   LOGICAL , DIMENSION(:), ALLOCATABLE ::   llnightav              ! Logical for calculating night-time averages

   TYPE(obs_surf), PUBLIC, POINTER, DIMENSION(:) ::   surfdata     !: Initial surface data
   TYPE(obs_surf), PUBLIC, POINTER, DIMENSION(:) ::   surfdataqc   !: Surface data after quality control
   TYPE(obs_prof), PUBLIC, POINTER, DIMENSION(:) ::   profdata     !: Initial profile data
   TYPE(obs_prof), PUBLIC, POINTER, DIMENSION(:) ::   profdataqc   !: Profile data after quality control

   CHARACTER(len=lca), PUBLIC, DIMENSION(:), ALLOCATABLE ::   cobstypesprof, cobstypessurf   !: Profile & surface obs types

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
   !! $Id: diaobs.F90 15077 2021-07-03 10:16:35Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_obs_init( Kmm )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs_init  ***
      !!          
      !! ** Purpose : Initialize and read observations
      !!
      !! ** Method  : Read the namelist and call reading routines
      !!
      !! ** Action  : Read the namelist and call reading routines
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)                ::   Kmm                      ! ocean time level indices
      INTEGER, PARAMETER                 ::   jpmaxnfiles = 1000       ! Maximum number of files for each obs type
      INTEGER, DIMENSION(:), ALLOCATABLE ::   ifilesprof, ifilessurf   ! Number of profile & surface files
      INTEGER :: ios             ! Local integer output status for namelist read
      INTEGER :: jtype           ! Counter for obs types
      INTEGER :: jvar            ! Counter for variables
      INTEGER :: jfile           ! Counter for files
      INTEGER :: jnumsstbias     ! Number of SST bias files to read and apply
      INTEGER :: n2dint_type     ! Local version of nn_2dint*
      !
      CHARACTER(len=128), DIMENSION(jpmaxnfiles) :: &
         & cn_profbfiles, &      ! T/S profile input filenames
         & cn_sstfbfiles, &      ! Sea surface temperature input filenames
         & cn_sssfbfiles, &      ! Sea surface salinity input filenames
         & cn_slafbfiles, &      ! Sea level anomaly input filenames
         & cn_sicfbfiles, &      ! Seaice concentration input filenames
         & cn_velfbfiles, &      ! Velocity profile input filenames
         & cn_sstbiasfiles       ! SST bias input filenames
      CHARACTER(LEN=128) :: &
         & cn_altbiasfile        ! Altimeter bias input filename
      CHARACTER(len=128), DIMENSION(:,:), ALLOCATABLE :: &
         & clproffiles, &        ! Profile filenames
         & clsurffiles           ! Surface filenames
      CHARACTER(len=8), DIMENSION(:), ALLOCATABLE :: &
         & clvars                ! Expected variable names
         !
      LOGICAL :: ln_t3d          ! Logical switch for temperature profiles
      LOGICAL :: ln_s3d          ! Logical switch for salinity profiles
      LOGICAL :: ln_sla          ! Logical switch for sea level anomalies 
      LOGICAL :: ln_sst          ! Logical switch for sea surface temperature
      LOGICAL :: ln_sss          ! Logical switch for sea surface salinity
      LOGICAL :: ln_sic          ! Logical switch for sea ice concentration
      LOGICAL :: ln_vel3d        ! Logical switch for velocity (u,v) obs
      LOGICAL :: ln_nea          ! Logical switch to remove obs near land
      LOGICAL :: ln_altbias      ! Logical switch for altimeter bias
      LOGICAL :: ln_sstbias      ! Logical switch for bias corection of SST 
      LOGICAL :: ln_ignmis       ! Logical switch for ignoring missing files
      LOGICAL :: ln_s_at_t       ! Logical switch to compute model S at T obs
      LOGICAL :: ln_bound_reject ! Logical to remove obs near boundaries in LAMs.
      LOGICAL :: ltype_fp_indegs ! Local version of ln_*_fp_indegs
      LOGICAL :: ltype_night     ! Local version of ln_sstnight (false for other variables)
      LOGICAL, DIMENSION(:), ALLOCATABLE :: llvar   ! Logical for profile variable read
      LOGICAL, DIMENSION(jpmaxnfiles) :: lmask ! Used for finding number of sstbias files
      !
      REAL(wp) :: rn_dobsini      ! Obs window start date YYYYMMDD.HHMMSS
      REAL(wp) :: rn_dobsend      ! Obs window end date   YYYYMMDD.HHMMSS
      REAL(wp) :: ztype_avglamscl ! Local version of rn_*_avglamscl
      REAL(wp) :: ztype_avgphiscl ! Local version of rn_*_avgphiscl
      REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE :: zglam   ! Model longitudes for profile variables
      REAL(wp), DIMENSION(:,:,:),   ALLOCATABLE :: zgphi   ! Model latitudes  for profile variables
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: zmask   ! Model land/sea mask associated with variables
      !!
      NAMELIST/namobs/ln_diaobs, ln_t3d, ln_s3d, ln_sla,              &
         &            ln_sst, ln_sic, ln_sss, ln_vel3d,               &
         &            ln_altbias, ln_sstbias, ln_nea,                 &
         &            ln_grid_global, ln_grid_search_lookup,          &
         &            ln_ignmis, ln_s_at_t, ln_bound_reject,          &
         &            ln_sstnight, ln_default_fp_indegs,              &
         &            ln_sla_fp_indegs, ln_sst_fp_indegs,             &
         &            ln_sss_fp_indegs, ln_sic_fp_indegs,             &
         &            cn_profbfiles, cn_slafbfiles,                   &
         &            cn_sstfbfiles, cn_sicfbfiles,                   &
         &            cn_velfbfiles, cn_sssfbfiles,                   &
         &            cn_sstbiasfiles, cn_altbiasfile,                &
         &            cn_gridsearchfile, rn_gridsearchres,            &
         &            rn_dobsini, rn_dobsend,                         &
         &            rn_default_avglamscl, rn_default_avgphiscl,     &
         &            rn_sla_avglamscl, rn_sla_avgphiscl,             &
         &            rn_sst_avglamscl, rn_sst_avgphiscl,             &
         &            rn_sss_avglamscl, rn_sss_avgphiscl,             &
         &            rn_sic_avglamscl, rn_sic_avgphiscl,             &
         &            nn_1dint, nn_2dint_default,                     &
         &            nn_2dint_sla, nn_2dint_sst,                     &
         &            nn_2dint_sss, nn_2dint_sic,                     &
         &            nn_msshc, rn_mdtcorr, rn_mdtcutoff,             &
         &            nn_profdavtypes
      !-----------------------------------------------------------------------

      !-----------------------------------------------------------------------
      ! Read namelist parameters
      !-----------------------------------------------------------------------
      ! Some namelist arrays need initialising
      cn_profbfiles  (:) = ''
      cn_slafbfiles  (:) = ''
      cn_sstfbfiles  (:) = ''
      cn_sicfbfiles  (:) = ''
      cn_velfbfiles  (:) = ''
      cn_sssfbfiles  (:) = ''
      cn_sstbiasfiles(:) = ''
      nn_profdavtypes(:) = -1

      CALL ini_date( rn_dobsini )
      CALL fin_date( rn_dobsend )

      ! Read namelist namobs : control observation diagnostics
      READ  ( numnam_ref, namobs, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namobs in reference namelist' )
      READ  ( numnam_cfg, namobs, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namobs in configuration namelist' )
      IF(lwm) WRITE ( numond, namobs )

      IF( .NOT.ln_diaobs ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dia_obs_init : NO Observation diagnostic used'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         RETURN
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_obs_init : Observation diagnostic initialization'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namobs : set observation diagnostic parameters' 
         WRITE(numout,*) '      Logical switch for T profile observations                ln_t3d = ', ln_t3d
         WRITE(numout,*) '      Logical switch for S profile observations                ln_s3d = ', ln_s3d
         WRITE(numout,*) '      Logical switch for SLA observations                      ln_sla = ', ln_sla
         WRITE(numout,*) '      Logical switch for SST observations                      ln_sst = ', ln_sst
         WRITE(numout,*) '      Logical switch for Sea Ice observations                  ln_sic = ', ln_sic
         WRITE(numout,*) '      Logical switch for velocity observations               ln_vel3d = ', ln_vel3d
         WRITE(numout,*) '      Logical switch for SSS observations                      ln_sss = ', ln_sss
         WRITE(numout,*) '      Global distribution of observations              ln_grid_global = ', ln_grid_global
         WRITE(numout,*) '      Logical switch for obs grid search lookup ln_grid_search_lookup = ', ln_grid_search_lookup
         IF (ln_grid_search_lookup) &
            WRITE(numout,*) '      Grid search lookup file header                cn_gridsearchfile = ', cn_gridsearchfile
         WRITE(numout,*) '      Initial date in window YYYYMMDD.HHMMSS               rn_dobsini = ', rn_dobsini
         WRITE(numout,*) '      Final date in window YYYYMMDD.HHMMSS                 rn_dobsend = ', rn_dobsend
         WRITE(numout,*) '      Type of vertical interpolation method                  nn_1dint = ', nn_1dint
         WRITE(numout,*) '      Default horizontal interpolation method        nn_2dint_default = ', nn_2dint_default
         WRITE(numout,*) '      Type of horizontal interpolation method for SLA    nn_2dint_sla = ', nn_2dint_sla
         WRITE(numout,*) '      Type of horizontal interpolation method for SST    nn_2dint_sst = ', nn_2dint_sst
         WRITE(numout,*) '      Type of horizontal interpolation method for SSS    nn_2dint_sss = ', nn_2dint_sss        
         WRITE(numout,*) '      Type of horizontal interpolation method for SIC    nn_2dint_sic = ', nn_2dint_sic
         WRITE(numout,*) '      Default E/W diameter of obs footprint      rn_default_avglamscl = ', rn_default_avglamscl
         WRITE(numout,*) '      Default N/S diameter of obs footprint      rn_default_avgphiscl = ', rn_default_avgphiscl
         WRITE(numout,*) '      Default obs footprint in deg [T] or m [F]  ln_default_fp_indegs = ', ln_default_fp_indegs
         WRITE(numout,*) '      SLA E/W diameter of obs footprint              rn_sla_avglamscl = ', rn_sla_avglamscl
         WRITE(numout,*) '      SLA N/S diameter of obs footprint              rn_sla_avgphiscl = ', rn_sla_avgphiscl
         WRITE(numout,*) '      SLA obs footprint in deg [T] or m [F]          ln_sla_fp_indegs = ', ln_sla_fp_indegs
         WRITE(numout,*) '      SST E/W diameter of obs footprint              rn_sst_avglamscl = ', rn_sst_avglamscl
         WRITE(numout,*) '      SST N/S diameter of obs footprint              rn_sst_avgphiscl = ', rn_sst_avgphiscl
         WRITE(numout,*) '      SST obs footprint in deg [T] or m [F]          ln_sst_fp_indegs = ', ln_sst_fp_indegs
         WRITE(numout,*) '      SIC E/W diameter of obs footprint              rn_sic_avglamscl = ', rn_sic_avglamscl
         WRITE(numout,*) '      SIC N/S diameter of obs footprint              rn_sic_avgphiscl = ', rn_sic_avgphiscl
         WRITE(numout,*) '      SIC obs footprint in deg [T] or m [F]          ln_sic_fp_indegs = ', ln_sic_fp_indegs
         WRITE(numout,*) '      Rejection of observations near land switch               ln_nea = ', ln_nea
         WRITE(numout,*) '      Rejection of obs near open bdys                 ln_bound_reject = ', ln_bound_reject
         WRITE(numout,*) '      MSSH correction scheme                                 nn_msshc = ', nn_msshc
         WRITE(numout,*) '      MDT  correction                                      rn_mdtcorr = ', rn_mdtcorr
         WRITE(numout,*) '      MDT cutoff for computed correction                 rn_mdtcutoff = ', rn_mdtcutoff
         WRITE(numout,*) '      Logical switch for alt bias                          ln_altbias = ', ln_altbias
         WRITE(numout,*) '      Logical switch for sst bias                          ln_sstbias = ', ln_sstbias
         WRITE(numout,*) '      Logical switch for ignoring missing files             ln_ignmis = ', ln_ignmis
         WRITE(numout,*) '      Daily average types                             nn_profdavtypes = ', nn_profdavtypes
         WRITE(numout,*) '      Logical switch for night-time SST obs               ln_sstnight = ', ln_sstnight
      ENDIF
      !-----------------------------------------------------------------------
      ! Set up list of observation types to be used
      ! and the files associated with each type
      !-----------------------------------------------------------------------

      nproftypes = COUNT( (/ln_t3d .OR. ln_s3d, ln_vel3d /) )
      nsurftypes = COUNT( (/ln_sla, ln_sst, ln_sic, ln_sss /) )

      IF( ln_sstbias ) THEN 
         lmask(:) = .FALSE. 
         WHERE( cn_sstbiasfiles(:) /= '' )   lmask(:) = .TRUE. 
         jnumsstbias = COUNT(lmask) 
         lmask(:) = .FALSE. 
      ENDIF      

      IF( nproftypes == 0 .AND. nsurftypes == 0 ) THEN
         CALL ctl_warn( 'dia_obs_init: ln_diaobs is set to true, but all obs operator logical flags',   &
            &           ' (ln_t3d, ln_s3d, ln_sla, ln_sst, ln_sic, ln_vel3d)',                          &
            &           ' are set to .FALSE. so turning off calls to dia_obs'  )
         ln_diaobs = .FALSE.
         RETURN
      ENDIF

      IF( nproftypes > 0 ) THEN
         !
         ALLOCATE( cobstypesprof(nproftypes)             )
         ALLOCATE( ifilesprof   (nproftypes)             )
         ALLOCATE( clproffiles  (nproftypes,jpmaxnfiles) )
         !
         jtype = 0
         IF( ln_t3d .OR. ln_s3d ) THEN
            jtype = jtype + 1
            cobstypesprof(jtype) = 'prof'
            clproffiles(jtype,:) = cn_profbfiles
         ENDIF
         IF( ln_vel3d ) THEN
            jtype = jtype + 1
            cobstypesprof(jtype) = 'vel'
            clproffiles(jtype,:) = cn_velfbfiles
         ENDIF
         !
         CALL obs_settypefiles( nproftypes, jpmaxnfiles, ifilesprof, cobstypesprof, clproffiles )
         !
      ENDIF

      IF( nsurftypes > 0 ) THEN
         !
         ALLOCATE( cobstypessurf(nsurftypes)             )
         ALLOCATE( ifilessurf   (nsurftypes)             )
         ALLOCATE( clsurffiles  (nsurftypes,jpmaxnfiles) )
         ALLOCATE( n2dintsurf   (nsurftypes)             )
         ALLOCATE( zavglamscl   (nsurftypes)             )
         ALLOCATE( zavgphiscl   (nsurftypes)             )
         ALLOCATE( lfpindegs    (nsurftypes)             )
         ALLOCATE( llnightav    (nsurftypes)             )
         !
         jtype = 0
         IF( ln_sla ) THEN
            jtype = jtype + 1
            cobstypessurf(jtype) = 'sla'
            clsurffiles(jtype,:) = cn_slafbfiles
         ENDIF
         IF( ln_sst ) THEN
            jtype = jtype + 1
            cobstypessurf(jtype) = 'sst'
            clsurffiles(jtype,:) = cn_sstfbfiles
         ENDIF
         IF( ln_sss ) THEN
            jtype = jtype + 1
            cobstypessurf(jtype) = 'sss'
            clsurffiles(jtype,:) = cn_sssfbfiles
         ENDIF
         !
         CALL obs_settypefiles( nsurftypes, jpmaxnfiles, ifilessurf, cobstypessurf, clsurffiles )

         DO jtype = 1, nsurftypes

            IF ( TRIM(cobstypessurf(jtype)) == 'sla' ) THEN
               IF ( nn_2dint_sla == -1 ) THEN
                  n2dint_type  = nn_2dint_default
               ELSE
                  n2dint_type  = nn_2dint_sla
               ENDIF
               ztype_avglamscl = rn_sla_avglamscl
               ztype_avgphiscl = rn_sla_avgphiscl
               ltype_fp_indegs = ln_sla_fp_indegs
               ltype_night     = .FALSE.
            ELSE IF ( TRIM(cobstypessurf(jtype)) == 'sst' ) THEN
               IF ( nn_2dint_sst == -1 ) THEN
                  n2dint_type  = nn_2dint_default
               ELSE
                  n2dint_type  = nn_2dint_sst
               ENDIF
               ztype_avglamscl = rn_sst_avglamscl
               ztype_avgphiscl = rn_sst_avgphiscl
               ltype_fp_indegs = ln_sst_fp_indegs
               ltype_night     = ln_sstnight
            ELSE IF ( TRIM(cobstypessurf(jtype)) == 'sic' ) THEN
               IF ( nn_2dint_sic == -1 ) THEN
                  n2dint_type  = nn_2dint_default
               ELSE
                  n2dint_type  = nn_2dint_sic
               ENDIF
               ztype_avglamscl = rn_sic_avglamscl
               ztype_avgphiscl = rn_sic_avgphiscl
               ltype_fp_indegs = ln_sic_fp_indegs
               ltype_night     = .FALSE.
            ELSE IF ( TRIM(cobstypessurf(jtype)) == 'sss' ) THEN
               IF ( nn_2dint_sss == -1 ) THEN
                  n2dint_type  = nn_2dint_default
               ELSE
                  n2dint_type  = nn_2dint_sss
               ENDIF
               ztype_avglamscl = rn_sss_avglamscl
               ztype_avgphiscl = rn_sss_avgphiscl
               ltype_fp_indegs = ln_sss_fp_indegs
               ltype_night     = .FALSE.
            ELSE
               n2dint_type     = nn_2dint_default
               ztype_avglamscl = rn_default_avglamscl
               ztype_avgphiscl = rn_default_avgphiscl
               ltype_fp_indegs = ln_default_fp_indegs
               ltype_night     = .FALSE.
            ENDIF
            
            CALL obs_setinterpopts( nsurftypes, jtype, cobstypessurf(jtype),       &
               &                    nn_2dint_default, n2dint_type,                 &
               &                    ztype_avglamscl, ztype_avgphiscl,              &
               &                    ltype_fp_indegs, ltype_night,                  &
               &                    n2dintsurf, zavglamscl, zavgphiscl,            &
               &                    lfpindegs, llnightav )

         END DO
         !
      ENDIF


      !-----------------------------------------------------------------------
      ! Obs operator parameter checking and initialisations
      !-----------------------------------------------------------------------
      !
      IF( ln_vel3d  .AND.  .NOT.ln_grid_global ) THEN
         CALL ctl_stop( 'Velocity data only works with ln_grid_global=.true.' )
         RETURN
      ENDIF
      !
      IF( ln_grid_global ) THEN
         IF( jpnij < jpni * jpnj ) THEN
            CALL ctl_stop( 'STOP', 'dia_obs_init: ln_grid_global=T is not available when land subdomains are suppressed' )
         END IF
         CALL ctl_warn( 'dia_obs_init: ln_grid_global=T may cause memory issues when used with a large number of processors' )
      ENDIF
      !
      IF( nn_1dint < 0  .OR.  nn_1dint > 1 ) THEN
         CALL ctl_stop('dia_obs_init: Choice of vertical (1D) interpolation method is not available')
      ENDIF
      !
      IF( nn_2dint_default < 0  .OR.  nn_2dint_default > 6  ) THEN
         CALL ctl_stop('dia_obs_init: Choice of default horizontal (2D) interpolation method is not available')
      ENDIF
      !
      CALL obs_typ_init
      IF( ln_grid_global )   CALL mppmap_init
      !
      CALL obs_grid_setup( )

      !-----------------------------------------------------------------------
      ! Depending on switches read the various observation types
      !-----------------------------------------------------------------------
      !
      IF( nproftypes > 0 ) THEN
         !
         ALLOCATE( profdata  (nproftypes) , nvarsprof (nproftypes) )
         ALLOCATE( profdataqc(nproftypes) , nextrprof (nproftypes) )
         !
         DO jtype = 1, nproftypes
            !
            IF ( TRIM(cobstypesprof(jtype)) == 'prof' ) THEN
               nvarsprof(jtype) = 2
               nextrprof(jtype) = 1             
               ALLOCATE( llvar (nvarsprof(jtype)) )
               ALLOCATE( clvars(nvarsprof(jtype)) )
               ALLOCATE( zglam(jpi, jpj,      nvarsprof(jtype)) )
               ALLOCATE( zgphi(jpi, jpj,      nvarsprof(jtype)) )
               ALLOCATE( zmask(jpi, jpj, jpk, nvarsprof(jtype)) )
               llvar(1)       = ln_t3d
               llvar(2)       = ln_s3d
               clvars(1)      = 'POTM'
               clvars(2)      = 'PSAL'
               zglam(:,:,1)   = glamt(:,:)
               zglam(:,:,2)   = glamt(:,:)
               zgphi(:,:,1)   = gphit(:,:)
               zgphi(:,:,2)   = gphit(:,:)
               zmask(:,:,:,1) = tmask(:,:,:)
               zmask(:,:,:,2) = tmask(:,:,:)
            ELSE IF ( TRIM(cobstypesprof(jtype)) == 'vel' )  THEN
               nvarsprof(jtype) = 2
               nextrprof(jtype) = 2
               ALLOCATE( llvar (nvarsprof(jtype)) )
               ALLOCATE( clvars(nvarsprof(jtype)) )
               ALLOCATE( zglam(jpi, jpj,      nvarsprof(jtype)) )
               ALLOCATE( zgphi(jpi, jpj,      nvarsprof(jtype)) )
               ALLOCATE( zmask(jpi, jpj, jpk, nvarsprof(jtype)) )
               llvar(1)       = ln_vel3d
               llvar(2)       = ln_vel3d
               clvars(1)      = 'UVEL'
               clvars(2)      = 'VVEL'
               zglam(:,:,1)   = glamu(:,:)
               zglam(:,:,2)   = glamv(:,:)
               zgphi(:,:,1)   = gphiu(:,:)
               zgphi(:,:,2)   = gphiv(:,:)
               zmask(:,:,:,1) = umask(:,:,:)
               zmask(:,:,:,2) = vmask(:,:,:)
            ELSE
               nvarsprof(jtype) = 1
               nextrprof(jtype) = 0
               ALLOCATE( llvar (nvarsprof(jtype)) )
               ALLOCATE( clvars(nvarsprof(jtype)) )
               ALLOCATE( zglam(jpi, jpj,      nvarsprof(jtype)) )
               ALLOCATE( zgphi(jpi, jpj,      nvarsprof(jtype)) )
               ALLOCATE( zmask(jpi, jpj, jpk, nvarsprof(jtype)) )
               llvar(1)       = .TRUE.
               zglam(:,:,1)   = glamt(:,:)
               zgphi(:,:,1)   = gphit(:,:)
               zmask(:,:,:,1) = tmask(:,:,:)
            ENDIF
            !
            ! Read in profile or profile obs types
            CALL obs_rea_prof( profdata(jtype), ifilesprof(jtype),       &
               &               clproffiles(jtype,1:ifilesprof(jtype)), &
               &               nvarsprof(jtype), nextrprof(jtype), nitend-nit000+2, &
               &               rn_dobsini, rn_dobsend, llvar, &
               &               ln_ignmis, ln_s_at_t, .FALSE., clvars, &
               &               kdailyavtypes = nn_profdavtypes )
               !
            DO jvar = 1, nvarsprof(jtype)
               CALL obs_prof_staend( profdata(jtype), jvar )
            END DO
            !
            CALL obs_pre_prof( profdata(jtype), profdataqc(jtype), &
               &               llvar, &
               &               jpi, jpj, jpk, &
               &               zmask, zglam, zgphi, &
               &               ln_nea, ln_bound_reject, Kmm, &
               &               kdailyavtypes = nn_profdavtypes )
            !
            DEALLOCATE( llvar, clvars, zglam, zgphi, zmask )
            !
         END DO
         !
         DEALLOCATE( ifilesprof, clproffiles )
         !
      ENDIF
      !
      IF( nsurftypes > 0 ) THEN
         !
         ALLOCATE( surfdata  (nsurftypes) , nvarssurf(nsurftypes) )
         ALLOCATE( surfdataqc(nsurftypes) , nextrsurf(nsurftypes) )
         !
         DO jtype = 1, nsurftypes
            !
            nvarssurf(jtype) = 1
            nextrsurf(jtype) = 0
            llnightav(jtype) = .FALSE.
            IF( TRIM(cobstypessurf(jtype)) == 'sla' )   nextrsurf(jtype) = 2
            IF( TRIM(cobstypessurf(jtype)) == 'sst' )   llnightav(jtype) = ln_sstnight
            !
            ALLOCATE( clvars( nvarssurf(jtype) ) )
            IF ( TRIM(cobstypessurf(jtype)) == 'sla' ) THEN
               clvars(1) = 'SLA'
            ELSE IF ( TRIM(cobstypessurf(jtype)) == 'sst' ) THEN
               clvars(1) = 'SST'
            ELSE IF ( TRIM(cobstypessurf(jtype)) == 'sic' ) THEN
               clvars(1) = 'ICECONC'
            ELSE IF ( TRIM(cobstypessurf(jtype)) == 'sss' ) THEN
               clvars(1) = 'SSS'
            ENDIF
            !
            ! Read in surface obs types
            CALL obs_rea_surf( surfdata(jtype), ifilessurf(jtype), &
               &               clsurffiles(jtype,1:ifilessurf(jtype)), &
               &               nvarssurf(jtype), nextrsurf(jtype), nitend-nit000+2, &
               &               rn_dobsini, rn_dobsend, ln_ignmis, .FALSE., llnightav(jtype), &
               &               clvars )
               !
            CALL obs_pre_surf( surfdata(jtype), surfdataqc(jtype), ln_nea, ln_bound_reject )
            !
            IF( TRIM(cobstypessurf(jtype)) == 'sla' ) THEN
               CALL obs_rea_mdt( surfdataqc(jtype), n2dintsurf(jtype), Kmm )
               IF( ln_altbias )   &
                  & CALL obs_rea_altbias ( surfdataqc(jtype), n2dintsurf(jtype), cn_altbiasfile )
            ENDIF
            !
            IF( TRIM(cobstypessurf(jtype)) == 'sst' .AND. ln_sstbias ) THEN
               jnumsstbias = 0
               DO jfile = 1, jpmaxnfiles
                  IF( TRIM(cn_sstbiasfiles(jfile)) /= '' )   jnumsstbias = jnumsstbias + 1
               END DO
               IF( jnumsstbias == 0 )   CALL ctl_stop("ln_sstbias set but no bias files to read in")    
               !
               CALL obs_app_sstbias( surfdataqc(jtype), n2dintsurf(jtype)             ,   & 
                  &                  jnumsstbias      , cn_sstbiasfiles(1:jnumsstbias) ) 
            ENDIF
            !
            DEALLOCATE( clvars )
         END DO
         !
         DEALLOCATE( ifilessurf, clsurffiles )
         !
      ENDIF
      !
   END SUBROUTINE dia_obs_init


   SUBROUTINE dia_obs( kstp, Kmm )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs  ***
      !!          
      !! ** Purpose : Call the observation operators on each time step
      !!
      !! ** Method  : Call the observation operators on each time step to
      !!              compute the model equivalent of the following data:
      !!               - Profile data, currently T/S or U/V
      !!               - Surface data, currently SST, SLA or sea-ice concentration.
      !!
      !! ** Action  :
      !!----------------------------------------------------------------------
      USE phycst , ONLY : rday                ! Physical constants
      USE oce    , ONLY : ts, uu, vv, ssh     ! Ocean dynamics and tracers variables (Kmm time-level only)
      USE phycst , ONLY : rday                ! Physical constants

      IMPLICIT NONE

      !! * Arguments
      INTEGER, INTENT(IN) :: kstp  ! Current timestep
      INTEGER, INTENT(in) :: Kmm   ! ocean time level indices
      !! * Local declarations
      INTEGER :: idaystp           ! Number of timesteps per day
      INTEGER :: jtype             ! Data loop variable
      INTEGER :: jvar              ! Variable number
      INTEGER :: ji, jj, jk        ! Loop counters
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         & zprofvar                ! Model values for variables in a prof ob
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
         & zprofmask               ! Mask associated with zprofvar
      REAL(wp), DIMENSION(jpi,jpj) :: &
         & zsurfvar, &             ! Model values equivalent to surface ob.
         & zsurfmask               ! Mask associated with surface variable
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zglam,    &             ! Model longitudes for prof variables
         & zgphi                   ! Model latitudes for prof variables
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: zdept, zdepw

      !-----------------------------------------------------------------------

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_obs : Call the observation operators', kstp
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      idaystp = NINT( rday / rn_Dt )

      !-----------------------------------------------------------------------
      ! Call the profile and surface observation operators
      !-----------------------------------------------------------------------

      IF ( nproftypes > 0 ) THEN

         ALLOCATE( zdept(jpi,jpj,jpk), zdepw(jpi,jpj,jpk) )
         DO jk = 1, jpk
            zdept(:,:,jk) = gdept_0(:,:,jk)
            zdepw(:,:,jk) = gdepw_0(:,:,jk)
         END DO

         DO jtype = 1, nproftypes

            ! Allocate local work arrays
            ALLOCATE( zprofvar (jpi, jpj, jpk, profdataqc(jtype)%nvar) )
            ALLOCATE( zprofmask(jpi, jpj, jpk, profdataqc(jtype)%nvar) )
            ALLOCATE( zglam    (jpi, jpj,      profdataqc(jtype)%nvar) )
            ALLOCATE( zgphi    (jpi, jpj,      profdataqc(jtype)%nvar) )  
                              
            ! Defaults which might change
            DO jvar = 1, profdataqc(jtype)%nvar
               zprofmask(:,:,:,jvar) = tmask(:,:,:)
               zglam(:,:,jvar)       = glamt(:,:)
               zgphi(:,:,jvar)       = gphit(:,:)
            END DO

            SELECT CASE ( TRIM(cobstypesprof(jtype)) )
            CASE('prof')
               zprofvar(:,:,:,1) = ts(:,:,:,jp_tem,Kmm)
               zprofvar(:,:,:,2) = ts(:,:,:,jp_sal,Kmm)
            CASE('vel')
               zprofvar(:,:,:,1) = uu(:,:,:,Kmm)
               zprofvar(:,:,:,2) = vv(:,:,:,Kmm)
               zprofmask(:,:,:,1) = umask(:,:,:)
               zprofmask(:,:,:,2) = vmask(:,:,:)
               zglam(:,:,1) = glamu(:,:)
               zglam(:,:,2) = glamv(:,:)
               zgphi(:,:,1) = gphiu(:,:)
               zgphi(:,:,2) = gphiv(:,:)
            CASE DEFAULT
               CALL ctl_stop( 'Unknown profile observation type '//TRIM(cobstypesprof(jtype))//' in dia_obs' )
            END SELECT

            DO jvar = 1, profdataqc(jtype)%nvar
               CALL obs_prof_opt( profdataqc(jtype), kstp, jpi, jpj, jpk,  &
                  &               nit000, idaystp, jvar,                   &
                  &               zprofvar(:,:,:,jvar),                    &
                  &               zdept(:,:,:), zdepw(:,:,:),      &
                  &               zprofmask(:,:,:,jvar),                   &
                  &               zglam(:,:,jvar), zgphi(:,:,jvar),        &
                  &               nn_1dint, nn_2dint_default,              &
                  &               kdailyavtypes = nn_profdavtypes )
            END DO
            
            DEALLOCATE( zprofvar, zprofmask, zglam, zgphi )

         END DO

         DEALLOCATE( zdept, zdepw )

      ENDIF

      IF ( nsurftypes > 0 ) THEN

         DO jtype = 1, nsurftypes

            !Defaults which might be changed
            zsurfmask(:,:) = tmask(:,:,1)

            SELECT CASE ( TRIM(cobstypessurf(jtype)) )
            CASE('sst')
               zsurfvar(:,:) = ts(:,:,1,jp_tem,Kmm)
            CASE('sla')
               zsurfvar(:,:) = ssh(:,:,Kmm)
            CASE('sss')
               zsurfvar(:,:) = ts(:,:,1,jp_sal,Kmm)
            CASE('sic')
               IF ( kstp == 0 ) THEN
                  IF ( lwp .AND. surfdataqc(jtype)%nsstpmpp(1) > 0 ) THEN
                     CALL ctl_warn( 'Sea-ice not initialised on zeroth '// &
                        &           'time-step but some obs are valid then.' )
                     WRITE(numout,*)surfdataqc(jtype)%nsstpmpp(1), &
                        &           ' sea-ice obs will be missed'
                  ENDIF
                  surfdataqc(jtype)%nsurfup = surfdataqc(jtype)%nsurfup + &
                     &                        surfdataqc(jtype)%nsstp(1)
                  CYCLE
               ELSE
                  CALL ctl_stop( ' Trying to run sea-ice observation operator', &
                     &           ' but no sea-ice model appears to have been defined' )
               ENDIF

            END SELECT

            CALL obs_surf_opt( surfdataqc(jtype), kstp, jpi, jpj,       &
               &               nit000, idaystp, zsurfvar, zsurfmask,    &
               &               n2dintsurf(jtype), llnightav(jtype),     &
               &               zavglamscl(jtype), zavgphiscl(jtype),     &
               &               lfpindegs(jtype) )

         END DO

      ENDIF

   END SUBROUTINE dia_obs

   SUBROUTINE dia_obs_wri
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE dia_obs_wri  ***
      !!          
      !! ** Purpose : Call observation diagnostic output routines
      !!
      !! ** Method  : Call observation diagnostic output routines
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  06-03  (K. Mogensen) Original code
      !!        !  06-05  (K. Mogensen) Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  07-03  (K. Mogensen) General handling of profiles
      !!        !  08-09  (M. Valdivieso) Velocity component (U,V) profiles
      !!        !  15-08  (M. Martin) Combined writing for prof and surf types
      !!----------------------------------------------------------------------
      !! * Modules used
      USE obs_rot_vel          ! Rotation of velocities

      IMPLICIT NONE

      !! * Local declarations
      INTEGER :: jtype                    ! Data set loop variable
      INTEGER :: jo, jvar, jk
      REAL(wp), DIMENSION(:), ALLOCATABLE :: &
         & zu, &
         & zv

      !-----------------------------------------------------------------------
      ! Depending on switches call various observation output routines
      !-----------------------------------------------------------------------

      IF ( nproftypes > 0 ) THEN

         DO jtype = 1, nproftypes

            IF ( TRIM(cobstypesprof(jtype)) == 'vel' ) THEN

               ! For velocity data, rotate the model velocities to N/S, E/W
               ! using the compressed data structure.
               ALLOCATE( &
                  & zu(profdataqc(jtype)%nvprot(1)), &
                  & zv(profdataqc(jtype)%nvprot(2))  &
                  & )

               CALL obs_rotvel( profdataqc(jtype), nn_2dint_default, zu, zv )

               DO jo = 1, profdataqc(jtype)%nprof
                  DO jvar = 1, 2
                     DO jk = profdataqc(jtype)%npvsta(jo,jvar), profdataqc(jtype)%npvend(jo,jvar)

                        IF ( jvar == 1 ) THEN
                           profdataqc(jtype)%var(jvar)%vmod(jk) = zu(jk)
                        ELSE
                           profdataqc(jtype)%var(jvar)%vmod(jk) = zv(jk)
                        ENDIF

                     END DO
                  END DO
               END DO

               DEALLOCATE( zu )
               DEALLOCATE( zv )

            END IF

            CALL obs_prof_decompress( profdataqc(jtype), &
               &                      profdata(jtype), .TRUE., numout )

            CALL obs_wri_prof( profdata(jtype) )

         END DO

      ENDIF

      IF ( nsurftypes > 0 ) THEN

         DO jtype = 1, nsurftypes

            CALL obs_surf_decompress( surfdataqc(jtype), &
               &                      surfdata(jtype), .TRUE., numout )

            CALL obs_wri_surf( surfdata(jtype) )

         END DO

      ENDIF

   END SUBROUTINE dia_obs_wri

   SUBROUTINE dia_obs_dealloc
      IMPLICIT NONE
      !!----------------------------------------------------------------------
      !!                    *** ROUTINE dia_obs_dealloc ***
      !!
      !!  ** Purpose : To deallocate data to enable the obs_oper online loop.
      !!               Specifically: dia_obs_init --> dia_obs --> dia_obs_wri
      !!
      !!  ** Method : Clean up various arrays left behind by the obs_oper.
      !!
      !!  ** Action :
      !!
      !!----------------------------------------------------------------------
      ! obs_grid deallocation
      CALL obs_grid_deallocate

      ! diaobs deallocation
      IF ( nproftypes > 0 ) &
         &   DEALLOCATE( cobstypesprof, profdata, profdataqc, nvarsprof, nextrprof )

      IF ( nsurftypes > 0 ) &
         &   DEALLOCATE( cobstypessurf, surfdata, surfdataqc, nvarssurf, nextrsurf, &
         &               n2dintsurf, zavglamscl, zavgphiscl, lfpindegs, llnightav )

   END SUBROUTINE dia_obs_dealloc

   SUBROUTINE calc_date( kstp, ddobs )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE calc_date  ***
      !!          
      !! ** Purpose : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Action  : Get date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Action  : Get initial date in double precision YYYYMMDD.HHMMSS format
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  06-10  (G. Smith) Calculates initial date the same as method for final date
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) New generic routine now deals with arbitrary initial time of day
      !!----------------------------------------------------------------------
      USE phycst, ONLY : &            ! Physical constants
         & rday
      USE dom_oce, ONLY : &           ! Ocean space and time domain variables
         & rn_Dt

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=wp), INTENT(OUT) :: ddobs                        ! Date in YYYYMMDD.HHMMSS
      INTEGER, INTENT(IN) :: kstp

      !! * Local declarations
      INTEGER :: iyea        ! date - (year, month, day, hour, minute)
      INTEGER :: imon
      INTEGER :: iday
      INTEGER :: ihou
      INTEGER :: imin
      INTEGER :: imday       ! Number of days in month.
      REAL(wp) :: zdayfrc    ! Fraction of day

      INTEGER, DIMENSION(12) ::   imonth_len    !: length in days of the months of the current year

      !!----------------------------------------------------------------------
      !! Initial date initialization (year, month, day, hour, minute)
      !!----------------------------------------------------------------------
      iyea =   ndate0 / 10000
      imon = ( ndate0 - iyea * 10000 ) / 100
      iday =   ndate0 - iyea * 10000 - imon * 100
      ihou =   nn_time0 / 100
      imin = ( nn_time0 - ihou * 100 ) 

      !!----------------------------------------------------------------------
      !! Compute number of days + number of hours + min since initial time
      !!----------------------------------------------------------------------
      zdayfrc = kstp * rn_Dt / rday
      zdayfrc = zdayfrc - aint(zdayfrc)
      imin = imin + int( zdayfrc * 24 * 60 ) 
      DO WHILE (imin >= 60) 
        imin=imin-60
        ihou=ihou+1
      END DO
      DO WHILE (ihou >= 24)
        ihou=ihou-24
        iday=iday+1
      END DO 
      iday = iday + kstp * rn_Dt / rday 

      !-----------------------------------------------------------------------
      ! Convert number of days (iday) into a real date
      !----------------------------------------------------------------------

      CALL calc_month_len( iyea, imonth_len )

      DO WHILE ( iday > imonth_len(imon) )
         iday = iday - imonth_len(imon)
         imon = imon + 1 
         IF ( imon > 12 ) THEN
            imon = 1
            iyea = iyea + 1
            CALL calc_month_len( iyea, imonth_len )  ! update month lengths
         ENDIF
      END DO

      !----------------------------------------------------------------------
      ! Convert it into YYYYMMDD.HHMMSS format.
      !----------------------------------------------------------------------
      ddobs = iyea * 10000_dp + imon * 100_dp + &
         &    iday + ihou * 0.01_dp + imin * 0.0001_dp

   END SUBROUTINE calc_date

   SUBROUTINE ini_date( ddobsini )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ini_date  ***
      !!          
      !! ** Purpose : Get initial date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) Change to call generic routine calc_date
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(KIND=wp), INTENT(OUT) :: ddobsini                   ! Initial date in YYYYMMDD.HHMMSS

      CALL calc_date( nit000 - 1, ddobsini )

   END SUBROUTINE ini_date

   SUBROUTINE fin_date( ddobsfin )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE fin_date  ***
      !!          
      !! ** Purpose : Get final date in double precision YYYYMMDD.HHMMSS format
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  06-03  (K. Mogensen)  Original code
      !!        !  06-05  (K. Mogensen)  Reformatted
      !!        !  06-10  (A. Weaver) Cleaning
      !!        !  10-05  (D. Lea) Update to month length calculation for NEMO vn3.2
      !!        !  2014-09  (D. Lea) Change to call generic routine calc_date
      !!----------------------------------------------------------------------

      IMPLICIT NONE

      !! * Arguments
      REAL(wp), INTENT(OUT) :: ddobsfin ! Final date in YYYYMMDD.HHMMSS

      CALL calc_date( nitend, ddobsfin )

   END SUBROUTINE fin_date
   
   SUBROUTINE obs_settypefiles( ntypes, jpmaxnfiles, ifiles, cobstypes, cfiles )

      INTEGER, INTENT(IN) :: ntypes      ! Total number of obs types
      INTEGER, INTENT(IN) :: jpmaxnfiles ! Maximum number of files allowed for each type
      INTEGER, DIMENSION(ntypes), INTENT(OUT) :: &
         &                   ifiles      ! Out number of files for each type
      CHARACTER(len=lca), DIMENSION(ntypes), INTENT(IN) :: &
         &                   cobstypes   ! List of obs types
      CHARACTER(len=128), DIMENSION(ntypes, jpmaxnfiles), INTENT(IN) :: &
         &                   cfiles      ! List of files for all types

      !Local variables
      INTEGER :: jfile
      INTEGER :: jtype

      DO jtype = 1, ntypes

         ifiles(jtype) = 0
         DO jfile = 1, jpmaxnfiles
            IF ( trim(cfiles(jtype,jfile)) /= '' ) &
                      ifiles(jtype) = ifiles(jtype) + 1
         END DO

         IF ( ifiles(jtype) == 0 ) THEN
              CALL ctl_stop( 'Logical for observation type '//TRIM(cobstypes(jtype))//   &
                 &           ' set to true but no files available to read' )
         ENDIF

         IF(lwp) THEN    
            WRITE(numout,*) '             '//cobstypes(jtype)//' input observation file names:'
            DO jfile = 1, ifiles(jtype)
               WRITE(numout,*) '                '//TRIM(cfiles(jtype,jfile))
            END DO
         ENDIF

      END DO

   END SUBROUTINE obs_settypefiles

   SUBROUTINE obs_setinterpopts( ntypes, jtype, ctypein,             &
              &                  n2dint_default, n2dint_type,        &
              &                  ravglamscl_type, ravgphiscl_type,   &
              &                  lfp_indegs_type, lavnight_type,     &
              &                  n2dint, ravglamscl, ravgphiscl,     &
              &                  lfpindegs, lavnight )

      INTEGER, INTENT(IN)  :: ntypes             ! Total number of obs types
      INTEGER, INTENT(IN)  :: jtype              ! Index of the current type of obs
      INTEGER, INTENT(IN)  :: n2dint_default     ! Default option for interpolation type
      INTEGER, INTENT(IN)  :: n2dint_type        ! Option for interpolation type
      REAL(wp), INTENT(IN) :: &
         &                    ravglamscl_type, & !E/W diameter of obs footprint for this type
         &                    ravgphiscl_type    !N/S diameter of obs footprint for this type
      LOGICAL, INTENT(IN)  :: lfp_indegs_type    !T=> footprint in degrees, F=> in metres
      LOGICAL, INTENT(IN)  :: lavnight_type      !T=> obs represent night time average
      CHARACTER(len=lca), INTENT(IN) :: ctypein 

      INTEGER, DIMENSION(ntypes), INTENT(INOUT) :: &
         &                    n2dint 
      REAL(wp), DIMENSION(ntypes), INTENT(INOUT) :: &
         &                    ravglamscl, ravgphiscl
      LOGICAL, DIMENSION(ntypes), INTENT(INOUT) :: &
         &                    lfpindegs, lavnight

      lavnight(jtype) = lavnight_type

      IF ( (n2dint_type >= 0) .AND. (n2dint_type <= 6) ) THEN
         n2dint(jtype) = n2dint_type
      ELSE IF ( n2dint_type == -1 ) THEN
         n2dint(jtype) = n2dint_default
      ELSE
         CALL ctl_stop(' Choice of '//TRIM(ctypein)//' horizontal (2D) interpolation method', &
           &                    ' is not available')
      ENDIF

      ! For averaging observation footprints set options for size of footprint 
      IF ( (n2dint(jtype) > 4) .AND. (n2dint(jtype) <= 6) ) THEN
         IF ( ravglamscl_type > 0._wp ) THEN
            ravglamscl(jtype) = ravglamscl_type
         ELSE
            CALL ctl_stop( 'Incorrect value set for averaging footprint '// &
                           'scale (ravglamscl) for observation type '//TRIM(ctypein) )      
         ENDIF

         IF ( ravgphiscl_type > 0._wp ) THEN
            ravgphiscl(jtype) = ravgphiscl_type
         ELSE
            CALL ctl_stop( 'Incorrect value set for averaging footprint '// &
                           'scale (ravgphiscl) for observation type '//TRIM(ctypein) )      
         ENDIF

         lfpindegs(jtype) = lfp_indegs_type 

      ENDIF

      ! Write out info 
      IF(lwp) THEN
         IF ( n2dint(jtype) <= 4 ) THEN
            WRITE(numout,*) '             '//TRIM(ctypein)// &
               &            ' model counterparts will be interpolated horizontally'
         ELSE IF ( n2dint(jtype) <= 6 ) THEN
            WRITE(numout,*) '             '//TRIM(ctypein)// &
               &            ' model counterparts will be averaged horizontally'
            WRITE(numout,*) '             '//'    with E/W scale: ',ravglamscl(jtype)
            WRITE(numout,*) '             '//'    with N/S scale: ',ravgphiscl(jtype)
            IF ( lfpindegs(jtype) ) THEN
                WRITE(numout,*) '             '//'    (in degrees)'
            ELSE
                WRITE(numout,*) '             '//'    (in metres)'
            ENDIF
         ENDIF
      ENDIF

   END SUBROUTINE obs_setinterpopts

END MODULE diaobs
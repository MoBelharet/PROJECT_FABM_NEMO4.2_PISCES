










MODULE dyndmp
   !!======================================================================
   !!                       ***  MODULE  dyndmp  ***
   !! Ocean dynamics: internal restoring trend on momentum (U and V current)
   !!                 This should only be used for C1D case in current form
   !!======================================================================
   !! History :  3.5  ! 2013-08  (D. Calvert)  Original code
   !!            3.6  ! 2014-08  (T. Graham) Modified to use netcdf file of
   !!                             restoration coefficients supplied to tradmp
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_dmp_alloc : allocate dyndmp arrays
   !!   dyn_dmp_init  : namelist read, parameter control and resto coeff.
   !!   dyn_dmp       : update the momentum trend with the internal damping
   !!----------------------------------------------------------------------
   USE oce            ! ocean: variables
   USE dom_oce        ! ocean: domain variables
   USE c1d            ! 1D vertical configuration
   USE tradmp         ! ocean: internal damping
   USE zdf_oce        ! ocean: vertical physics
   USE phycst         ! physical constants
   USE dtauvd         ! data: U & V current
   USE zdfmxl         ! vertical physics: mixed layer depth
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing
   USE iom            ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_dmp_init ! routine called by nemogcm.F90
   PUBLIC   dyn_dmp      ! routine called by step_c1d.F90

   LOGICAL, PUBLIC ::   ln_dyndmp   !: Flag for Newtonian damping

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  utrdmp    !: damping U current trend (m/s2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  vtrdmp    !: damping V current trend (m/s2)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::  resto_uv  !: restoring coeff. on U & V current

   !! * Substitutions




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
   !! $Id: dyndmp.F90 15062 2021-06-28 11:19:48Z jchanut $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dyn_dmp_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION dyn_dmp_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( utrdmp(jpi,jpj,jpk), vtrdmp(jpi,jpj,jpk), resto_uv(jpi,jpj,jpk), STAT= dyn_dmp_alloc )
      !
      CALL mpp_sum ( 'dyndmp', dyn_dmp_alloc )
      IF( dyn_dmp_alloc > 0 )   CALL ctl_warn('dyn_dmp_alloc: allocation of arrays failed')
      !
   END FUNCTION dyn_dmp_alloc


   SUBROUTINE dyn_dmp_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_dmp_init  ***
      !! 
      !! ** Purpose :   Initialization for the Newtonian damping 
      !!
      !! ** Method  : - read the ln_dyndmp parameter from the namc1d_dyndmp namelist
      !!              - allocate damping arrays
      !!              - check the parameters of the namtra_dmp namelist
      !!              - calculate damping coefficient
      !!----------------------------------------------------------------------
      INTEGER ::   ios, imask   ! local integers
      !!
      NAMELIST/namc1d_dyndmp/ ln_dyndmp
      !!----------------------------------------------------------------------
      !
      READ  ( numnam_ref, namc1d_dyndmp, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namc1d_dyndmp in reference namelist' )
      READ  ( numnam_cfg, namc1d_dyndmp, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namc1d_dyndmp in configuration namelist' )
      IF(lwm) WRITE ( numond, namc1d_dyndmp )
      !
      IF(lwp) THEN                           ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_dmp_init : U and V current Newtonian damping'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namc1d_dyndmp : Set damping flag'
         WRITE(numout,*) '      add a damping term or not       ln_dyndmp = ', ln_dyndmp
         WRITE(numout,*) '   Namelist namtra_dmp    : Set damping parameters'
         WRITE(numout,*) '      Apply relaxation   or not       ln_tradmp = ', ln_tradmp
         WRITE(numout,*) '      mixed layer damping option      nn_zdmp   = ', nn_zdmp
         WRITE(numout,*) '      Damping file name               cn_resto  = ', cn_resto
         WRITE(numout,*)
      ENDIF
      !
      IF( ln_dyndmp ) THEN
         !                                   !==   allocate the data arrays   ==!
         IF( dyn_dmp_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'dyn_dmp_init: unable to allocate arrays' )
         !
         SELECT CASE ( nn_zdmp )             !==   control print of vertical option   ==!
         CASE ( 0    )   ;   IF(lwp) WRITE(numout,*) '   momentum damping throughout the water column'
         CASE ( 1    )   ;   IF(lwp) WRITE(numout,*) '   no momentum damping in the turbocline (avt > 5 cm2/s)'
         CASE ( 2    )   ;   IF(lwp) WRITE(numout,*) '   no momentum damping in the mixed layer'
         CASE DEFAULT
            WRITE(ctmp1,*) '          bad flag value for nn_zdmp = ', nn_zdmp
            CALL ctl_stop(ctmp1)
         END SELECT
         !
         IF( .NOT. ln_uvd_dyndmp ) THEN      ! force the initialization of U & V current data for damping
            CALL ctl_warn( 'dyn_dmp_init: U & V current read data not initialized, we force ln_uvd_dyndmp=T' )
            CALL dta_uvd_init( ld_dyndmp=ln_dyndmp )
         ENDIF
         !
         utrdmp(:,:,:) = 0._wp               ! internal damping trends
         vtrdmp(:,:,:) = 0._wp
         !
         !Read in mask from file
         CALL iom_open ( cn_resto, imask)
         CALL iom_get  ( imask, jpdom_auto, 'resto', resto)
         CALL iom_close( imask )
      ENDIF
      !
   END SUBROUTINE dyn_dmp_init


   SUBROUTINE dyn_dmp( kt, Kbb, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_dmp  ***
      !!                  
      !! ** Purpose :   Compute the momentum trends due to a newtonian damping
      !!      of the ocean velocities towards the given data and add it to the 
      !!      general momentum trends.
      !!
      !! ** Method  :   Compute Newtonian damping towards u_dta and v_dta 
      !!      and add to the general momentum trends:
      !!                     puu(Krhs) = puu(Krhs) + resto_uv * (u_dta - puu(Kbb))
      !!                     pvv(Krhs) = pvv(Krhs) + resto_uv * (v_dta - pvv(Kbb))
      !!      The trend is computed either throughout the water column
      !!      (nn_zdmp=0), where the vertical mixing is weak (nn_zdmp=1) or
      !!      below the well mixed layer (nn_zdmp=2)
      !!
      !! ** Action  : - (puu(:,:,:,Krhs),pvv(:,:,:,Krhs))   momentum trends updated with the damping trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::   kt             ! ocean time-step index
      INTEGER                             , INTENT(in   ) ::   Kbb, Kmm, Krhs ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::   puu, pvv       ! ocean velocities and RHS of momentum equation
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(dp), DIMENSION(jpi,jpj,jpk,2) ::   zuv_dta   ! Read in data
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start( 'dyn_dmp' )
      !
      !
      !                           !==   read and interpolate U & V current data at kt   ==!
      CALL dta_uvd( kt, Kmm, zuv_dta(:,:,:,1), zuv_dta(:,:,:,2))
      !
      SELECT CASE ( nn_zdmp )     !==   Calculate/add Newtonian damping to the momentum trend   ==!
      !
      CASE( 0 )                   ! Newtonian damping throughout the water column
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            zua = resto_uv(ji,jj,jk) * ( zuv_dta(ji,jj,jk,1) - puu(ji,jj,jk,Kbb) )
            zva = resto_uv(ji,jj,jk) * ( zuv_dta(ji,jj,jk,2) - pvv(ji,jj,jk,Kbb) )
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zua
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zva
            utrdmp(ji,jj,jk) = zua           ! save the trends
            vtrdmp(ji,jj,jk) = zva      
         END DO   ;   END DO   ;   END DO
         !
      CASE ( 1 )                  ! no damping above the turbocline (avt > 5 cm2/s)
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            IF( avt(ji,jj,jk) <= avt_c ) THEN
               zua = resto_uv(ji,jj,jk) * ( zuv_dta(ji,jj,jk,1) - puu(ji,jj,jk,Kbb) )
               zva = resto_uv(ji,jj,jk) * ( zuv_dta(ji,jj,jk,2) - pvv(ji,jj,jk,Kbb) )
            ELSE
               zua = 0._wp
               zva = 0._wp  
            ENDIF
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zua
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zva
            utrdmp(ji,jj,jk) = zua           ! save the trends
            vtrdmp(ji,jj,jk) = zva
         END DO   ;   END DO   ;   END DO
         !
      CASE ( 2 )                  ! no damping in the mixed layer
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            IF( gdept_0(ji,jj,jk) >= hmlp (ji,jj) ) THEN
               zua = resto_uv(ji,jj,jk) * ( zuv_dta(ji,jj,jk,1) - puu(ji,jj,jk,Kbb) )
               zva = resto_uv(ji,jj,jk) * ( zuv_dta(ji,jj,jk,2) - pvv(ji,jj,jk,Kbb) )
            ELSE
               zua = 0._wp
               zva = 0._wp  
            ENDIF
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) + zua
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) + zva
            utrdmp(ji,jj,jk) = zua           ! save the trends
            vtrdmp(ji,jj,jk) = zva
         END DO   ;   END DO   ;   END DO
         !
      END SELECT
      !
      !                           ! Control print
      IF( sn_cfctl%l_prtctl   )   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' dmp  - Ua: ', mask1=umask,   &
         &                                      tab3d_2=pvv(:,:,:,Krhs), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      !
      IF( ln_timing )   CALL timing_stop( 'dyn_dmp')
      !
   END SUBROUTINE dyn_dmp

   !!======================================================================
END MODULE dyndmp

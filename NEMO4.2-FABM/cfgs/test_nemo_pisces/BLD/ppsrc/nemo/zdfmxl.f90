










MODULE zdfmxl
   !!======================================================================
   !!                       ***  MODULE  zdfmxl  ***
   !! Ocean physics: mixed layer depth 
   !!======================================================================
   !! History :  1.0  ! 2003-08  (G. Madec)  original code
   !!            3.2  ! 2009-07  (S. Masson, G. Madec)  IOM + merge of DO-loop
   !!            3.7  ! 2012-03  (G. Madec)  make public the density criteria for trdmxl 
   !!             -   ! 2014-02  (F. Roquet)  mixed layer depth calculated using N2 instead of rhop 
   !!----------------------------------------------------------------------
   !!   zdf_mxl      : Compute the turbocline and mixed layer depths.
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE isf_oce        ! ice shelf
   USE dom_oce        ! ocean space and time domain variables
   USE trc_oce  , ONLY: l_offline         ! ocean space and time domain variables
   USE zdf_oce        ! ocean vertical physics
   !
   USE in_out_manager ! I/O manager
   USE prtctl         ! Print control
   USE phycst         ! physical constants
   USE iom            ! I/O library
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_mxl, zdf_mxl_turb, zdf_mxl_alloc   ! called by zdfphy.F90

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   nmln    !: number of level in the mixed layer (used by LDF, ZDF, TRD, TOP)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmld    !: mixing layer depth (turbocline)      [m]   (used by TOP)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmlp    !: mixed layer depth  (rho=rho0+zdcrit) [m]   (used by LDF)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmlpt   !: depth of the last T-point inside the mixed layer [m] (used by LDF)

   REAL(wp), PUBLIC ::   rho_c = 0.01_wp    !: density criterion for mixed layer depth
   REAL(wp), PUBLIC ::   avt_c = 5.e-4_wp   ! Kz criterion for the turbocline depth

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
   !! $Id: zdfmxl.F90 15249 2021-09-13 09:59:09Z hadcv $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_mxl_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION zdf_mxl_alloc  ***
      !!----------------------------------------------------------------------
      zdf_mxl_alloc = 0      ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( nmln ) ) THEN
         ALLOCATE( nmln(jpi,jpj), hmld(jpi,jpj), hmlp(jpi,jpj), hmlpt(jpi,jpj), STAT= zdf_mxl_alloc )
         !
         CALL mpp_sum ( 'zdfmxl', zdf_mxl_alloc )
         IF( zdf_mxl_alloc /= 0 )   CALL ctl_stop( 'STOP', 'zdf_mxl_alloc: failed to allocate arrays.' )
         !
      ENDIF
   END FUNCTION zdf_mxl_alloc


   SUBROUTINE zdf_mxl( kt, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdfmxl  ***
      !!                   
      !! ** Purpose :   Compute the mixed layer depth with density criteria.
      !!
      !! ** Method  :   The mixed layer depth is the shallowest W depth with 
      !!      the density of the corresponding T point (just bellow) bellow a
      !!      given value defined locally as rho(10m) + rho_c
      !!
      !! ** Action  :   nmln, hmlp, hmlpt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      INTEGER  ::   iik, ikt        ! local integer
      REAL(wp) ::   zN2_c           ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
         ENDIF
      ENDIF
      !
      ! w-level of the mixing and mixed layers
      DO jj = ntsj-(  nn_hls-( nn_hls+ nn_hls )*nthb), ntej+(  nn_hls -( nn_hls + nn_hls)*ntht) ; DO ji = ntsi-( nn_hls-( nn_hls+ nn_hls)*nthl), ntei+(  nn_hls-( nn_hls+ nn_hls)*nthr)
         nmln(ji,jj)  = nlb10                  ! Initialization to the number of w ocean point
         hmlp(ji,jj)  = 0._wp                  ! here hmlp used as a dummy variable, integrating vertically N^2
      END DO   ;   END DO
      zN2_c = grav * rho_c * r1_rho0      ! convert density criteria into N^2 criteria
      DO jk =  nlb10,  jpkm1  ; DO jj = ntsj-(   nn_hls-(  nn_hls+  nn_hls)*nthb), ntej+(   nn_hls-(  nn_hls+  nn_hls)*ntht) ; DO ji = ntsi-( nn_hls-( nn_hls+  nn_hls)*nthl), ntei+(   nn_hls-(  nn_hls+ nn_hls)*nthr)   ! Mixed layer level: w-level
         ikt = mbkt(ji,jj)
         hmlp(ji,jj) =   &
            & hmlp(ji,jj) + MAX( rn2b(ji,jj,jk) , 0._wp ) * e3w_0(ji,jj,jk)
         IF( hmlp(ji,jj) < zN2_c )   nmln(ji,jj) = MIN( jk , ikt ) + 1   ! Mixed layer level
      END DO   ;   END DO   ;   END DO
      ! depth of the mixed layer
      DO jj = ntsj-(  nn_hls-( nn_hls+ nn_hls )*nthb), ntej+(  nn_hls -( nn_hls + nn_hls)*ntht) ; DO ji = ntsi-( nn_hls-( nn_hls+ nn_hls)*nthl), ntei+(  nn_hls-( nn_hls+ nn_hls)*nthr)
         iik = nmln(ji,jj)
         hmlp (ji,jj) = gdepw_0(ji,jj,iik  ) * ssmask(ji,jj)    ! Mixed layer depth
         hmlpt(ji,jj) = gdept_0(ji,jj,iik-1) * ssmask(ji,jj)    ! depth of the last T-point inside the mixed layer
      END DO   ;   END DO
      !
      IF( .NOT.l_offline .AND. iom_use("mldr10_1") ) THEN
         IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN         ! Do only on the last tile
            IF( ln_isfcav ) THEN  ;  CALL iom_put( "mldr10_1", hmlp - risfdep)   ! mixed layer thickness
            ELSE                  ;  CALL iom_put( "mldr10_1", hmlp )            ! mixed layer depth
            END IF
         ENDIF
      ENDIF
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab2d_1=REAL(nmln,dp), clinfo1=' nmln : ', tab2d_2=hmlp, clinfo2=' hmlp : ' )
      !
   END SUBROUTINE zdf_mxl


   SUBROUTINE zdf_mxl_turb( kt, Kmm )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_mxl_turb  ***
      !!
      !! ** Purpose :   Compute the turbocline depth.
      !!
      !! ** Method  :   The turbocline depth is the depth at which the vertical
      !!      eddy diffusivity coefficient (resulting from the vertical physics
      !!      alone, not the isopycnal part, see trazdf.F) fall below a given
      !!      value defined locally (avt_c here taken equal to 5 cm/s2 by default)
      !!
      !! ** Action  :   hmld
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      INTEGER, INTENT(in) ::   Kmm  ! ocean time level index
      !
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      INTEGER  ::   iik             ! local integer
      INTEGER, DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls)) ::   imld   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      ! w-level of the turbocline and mixing layer (iom_use)
      imld(:,:) = mbkt(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls)) + 1                ! Initialization to the number of w ocean point
      DO jk =  jpkm1,  nlb10,  -1  ; DO jj = ntsj-(  1), ntej+(  1) ; DO ji = ntsi-( 1), ntei+(  1)   ! from the bottom to nlb10
         IF( avt (ji,jj,jk) < avt_c * wmask(ji,jj,jk) )   imld(ji,jj) = jk      ! Turbocline
      END DO   ;   END DO   ;   END DO
      ! depth of the mixing layer
      DO jj = ntsj-(  1-( 1+ 1 )*nthb), ntej+(  1 -( 1 + 1)*ntht) ; DO ji = ntsi-( 1-( 1+ 1)*nthl), ntei+(  1-( 1+ 1)*nthr)
         iik = imld(ji,jj)
         hmld (ji,jj) = gdepw_0(ji,jj,iik  ) * ssmask(ji,jj)    ! Turbocline depth
      END DO   ;   END DO
      !
      IF( .NOT.l_offline .AND. iom_use("mldkz5") ) THEN
         IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN         ! Do only on the last tile
            IF( ln_isfcav ) THEN  ;  CALL iom_put( "mldkz5"  , hmld - risfdep )   ! turbocline thickness
            ELSE                  ;  CALL iom_put( "mldkz5"  , hmld )             ! turbocline depth
            END IF
         ENDIF
      ENDIF
      !
   END SUBROUTINE zdf_mxl_turb
   !!======================================================================
END MODULE zdfmxl












MODULE diawri
   !!======================================================================
   !!                     ***  MODULE  diawri  ***
   !! Ocean diagnostics :  write ocean output files
   !!=====================================================================
   !! History :  OPA  ! 1991-03  (M.-A. Foujols)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!                 ! 1992-06  (M. Imbard)  correction restart file
   !!                 ! 1992-07  (M. Imbard)  split into diawri and rstwri
   !!                 ! 1993-03  (M. Imbard)  suppress writibm
   !!                 ! 1998-01  (C. Levy)  NETCDF format using ioipsl INTERFACE
   !!                 ! 1999-02  (E. Guilyardi)  name of netCDF files + variables
   !!            8.2  ! 2000-06  (M. Imbard)  Original code (diabort.F)
   !!   NEMO     1.0  ! 2002-06  (A.Bozec, E. Durand)  Original code (diainit.F)
   !!             -   ! 2002-09  (G. Madec)  F90: Free form and module
   !!             -   ! 2002-12  (G. Madec)  merge of diabort and diainit, F90
   !!                 ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            3.2  ! 2008-11  (B. Lemaire) creation from old diawri
   !!            3.7  ! 2014-01  (G. Madec) remove eddy induced velocity from no-IOM output
   !!                 !                     change name of output variables in dia_wri_state
   !!            4.0  ! 2020-10  (A. Nasser, S. Techene) add diagnostic for SWE
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dia_wri       : create the standart output files
   !!   dia_wri_state : create an output NetCDF file for a single instantaeous ocean state and forcing fields
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers 
   USE isf_oce
   USE isfcpl
   USE abl            ! abl variables in case ln_abl = .true.
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE dianam         ! build name of file (routine)
   USE diahth         ! thermocline diagnostics
   USE dynadv   , ONLY: ln_dynadv_vec
   USE icb_oce        ! Icebergs
   USE icbdia         ! Iceberg budgets
   USE ldftra         ! lateral physics: eddy diffusivity coef.
   USE ldfdyn         ! lateral physics: eddy viscosity   coef.
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE sbcssr         ! restoring term toward SST/SSS climatology
   USE sbcwave        ! wave parameters
   USE wet_dry        ! wetting and drying
   USE zdf_oce        ! ocean vertical physics
   USE zdfdrg         ! ocean vertical physics: top/bottom friction
   USE zdfmxl         ! mixed layer
   USE zdfosm         ! mixed layer
   !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE in_out_manager ! I/O manager
   USE dia25h         ! 25h Mean output
   USE iom            ! 
   USE ioipsl         ! 

   USE lib_mpp         ! MPP library
   USE timing          ! preformance summary
   USE diu_bulk        ! diurnal warm layer
   USE diu_coolskin    ! Cool skin

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dia_wri                 ! routines called by step.F90
   PUBLIC   dia_wri_state
   PUBLIC   dia_wri_alloc           ! Called by nemogcm module
   INTEGER ::   nid_T, nz_T, nh_T, ndim_T, ndim_hT   ! grid_T file
   INTEGER ::          nb_T              , ndim_bT   ! grid_T file
   INTEGER ::   nid_U, nz_U, nh_U, ndim_U, ndim_hU   ! grid_U file
   INTEGER ::   nid_V, nz_V, nh_V, ndim_V, ndim_hV   ! grid_V file
   INTEGER ::   nid_W, nz_W, nh_W                    ! grid_W file
   INTEGER ::   nid_A, nz_A, nh_A, ndim_A, ndim_hA   ! grid_ABL file   
   INTEGER ::   ndex(1)                              ! ???
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hT, ndex_hU, ndex_hV
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_hA, ndex_A ! ABL
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_T, ndex_U, ndex_V
   INTEGER, SAVE, ALLOCATABLE, DIMENSION(:) :: ndex_bT

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
   !! $Id: diawri.F90 15141 2021-07-23 14:20:12Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   'key_xios'                                        use IOM library
   !!----------------------------------------------------------------------
   INTEGER FUNCTION dia_wri_alloc()
      !
      dia_wri_alloc = 0
      !
   END FUNCTION dia_wri_alloc

   
   SUBROUTINE dia_wri( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_wri  ***
      !!                   
      !! ** Purpose :   Standard output of opa: dynamics and tracer fields 
      !!      NETCDF format is used by default 
      !!
      !! ** Method  :  use iom_put
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER, INTENT( in ) ::   Kmm     ! ocean time level index
      !!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
      INTEGER ::   ikbot            ! local integer
      REAL(wp)::   zztmp , zztmpx   ! local scalar
      REAL(wp)::   zztmp2, zztmpy   !   -      -
      REAL(wp)::   ze3
      REAL(wp), DIMENSION(ntsi-(     0):ntei+(     0),ntsj-(     0):ntej+(     0))     ::   z2d   ! 2D workspace
      REAL(wp), DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls),jpk) ::   z3d   ! 3D workspace
      !!----------------------------------------------------------------------
      ! 
      IF( ln_timing )   CALL timing_start('dia_wri')
      ! 
      ! Output the initial state and forcings
      IF( ninist == 1 ) THEN                       
         CALL dia_wri_state( Kmm, 'output.init' )
         ninist = 0
      ENDIF

      ! initialize arrays
      z2d(:,:)   = 0._wp
      z3d(:,:,:) = 0._wp
      
      ! Output of initial vertical scale factor
      CALL iom_put("e3t_0", e3t_0(:,:,:) )
      CALL iom_put("e3u_0", e3u_0(:,:,:) )
      CALL iom_put("e3v_0", e3v_0(:,:,:) )
      CALL iom_put("e3f_0", e3f_0(:,:,:) )
      !
      IF ( iom_use("tpt_dep") ) THEN
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z3d(ji,jj,jk) = gdept_0(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "tpt_dep", z3d )
      ENDIF

      ! --- vertical scale factors --- !
      IF ( iom_use("e3t") .OR. iom_use("e3tdef") ) THEN  ! time-varying e3t
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z3d(ji,jj,jk) =  e3t_0(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "e3t", z3d )
         IF ( iom_use("e3tdef") ) THEN
            DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z3d(ji,jj,jk) = ( ( z3d(ji,jj,jk) - e3t_0(ji,jj,jk) ) / e3t_0(ji,jj,jk) * 100._wp * tmask(ji,jj,jk) ) ** 2
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "e3tdef", z3d ) 
         ENDIF
      ENDIF 
      IF ( iom_use("e3u") ) THEN                         ! time-varying e3u
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z3d(ji,jj,jk) =  e3u_0(ji,jj,jk)
         END DO   ;   END DO   ;   END DO 
         CALL iom_put( "e3u" , z3d )
      ENDIF
      IF ( iom_use("e3v") ) THEN                         ! time-varying e3v
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z3d(ji,jj,jk) =  e3v_0(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "e3v" , z3d )
      ENDIF
      IF ( iom_use("e3w") ) THEN                         ! time-varying e3w
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z3d(ji,jj,jk) =  e3w_0(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "e3w" , z3d )
      ENDIF
      IF ( iom_use("e3f") ) THEN                         ! time-varying e3f caution here at Kaa
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z3d(ji,jj,jk) =  e3f_0(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "e3f" , z3d )
      ENDIF

      IF ( iom_use("ssh") ) THEN
         IF( ll_wd ) THEN                                ! sea surface height (brought back to the reference used for wetting and drying)
            CALL iom_put( "ssh" , (ssh(:,:,Kmm)+ssh_ref)*ssmask(:,:) )
         ELSE
            CALL iom_put( "ssh" ,  ssh(:,:,Kmm) )        ! sea surface height
         ENDIF
      ENDIF

      IF( iom_use("wetdep") )    CALL iom_put( "wetdep" , ht_0(:,:) + ssh(:,:,Kmm) )   ! wet depth
         

      ! --- tracers T&S --- !      
      CALL iom_put( "toce", ts(:,:,:,jp_tem,Kmm) )    ! 3D temperature
      CALL iom_put(  "sst", ts(:,:,1,jp_tem,Kmm) )    ! surface temperature

      IF ( iom_use("sbt") ) THEN
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            ikbot = mbkt(ji,jj)
            z2d(ji,jj) = ts(ji,jj,ikbot,jp_tem,Kmm)
         END DO   ;   END DO
         CALL iom_put( "sbt", z2d )                ! bottom temperature
      ENDIF
      
      CALL iom_put( "soce", ts(:,:,:,jp_sal,Kmm) )    ! 3D salinity
      CALL iom_put(  "sss", ts(:,:,1,jp_sal,Kmm) )    ! surface salinity
      IF ( iom_use("sbs") ) THEN
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            ikbot = mbkt(ji,jj)
            z2d(ji,jj) = ts(ji,jj,ikbot,jp_sal,Kmm)
         END DO   ;   END DO
         CALL iom_put( "sbs", z2d )                ! bottom salinity
      ENDIF

      IF( .NOT.lk_SWE )   CALL iom_put( "rhop", rhop(:,:,:) )          ! 3D potential density (sigma0)

      ! --- momentum --- !
      IF ( iom_use("taubot") ) THEN                ! bottom stress
         zztmp = rho0 * 0.25_wp
         z2d(:,:) = 0._wp
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zztmp2 = (  ( rCdU_bot(ji+1,jj)+rCdU_bot(ji  ,jj) ) * uu(ji  ,jj,mbku(ji  ,jj),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji  ,jj)+rCdU_bot(ji-1,jj) ) * uu(ji-1,jj,mbku(ji-1,jj),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj  ) ) * vv(ji,jj  ,mbkv(ji,jj  ),Kmm)  )**2   &
               &   + (  ( rCdU_bot(ji,jj  )+rCdU_bot(ji,jj-1) ) * vv(ji,jj-1,mbkv(ji,jj-1),Kmm)  )**2
            z2d(ji,jj) = zztmp * SQRT( zztmp2 ) * tmask(ji,jj,1) 
            !
         END DO   ;   END DO
         CALL iom_put( "taubot", z2d )           
      ENDIF
         
      CALL iom_put( "uoce", uu(:,:,:,Kmm) )            ! 3D i-current
      CALL iom_put(  "ssu", uu(:,:,1,Kmm) )            ! surface i-current
      IF ( iom_use("sbu") ) THEN
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            ikbot = mbku(ji,jj)
            z2d(ji,jj) = uu(ji,jj,ikbot,Kmm)
         END DO   ;   END DO
         CALL iom_put( "sbu", z2d )                ! bottom i-current
      ENDIF
      
      CALL iom_put( "voce", vv(:,:,:,Kmm) )            ! 3D j-current
      CALL iom_put(  "ssv", vv(:,:,1,Kmm) )            ! surface j-current
      IF ( iom_use("sbv") ) THEN
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            ikbot = mbkv(ji,jj)
            z2d(ji,jj) = vv(ji,jj,ikbot,Kmm)
         END DO   ;   END DO
         CALL iom_put( "sbv", z2d )                ! bottom j-current
      ENDIF

      !                                            ! vertical velocity
      IF( ln_zad_Aimp ) THEN
         IF( iom_use('woce') ) THEN
            DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z3d(ji,jj,jk) = ww(ji,jj,jk) + wi(ji,jj,jk)
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "woce", z3d )   ! explicit plus implicit parts
         ENDIF
      ELSE
         CALL iom_put( "woce", ww )
      ENDIF

      IF( iom_use('w_masstr') .OR. iom_use('w_masstr2') ) THEN   ! vertical mass transport & its square value
         !                     ! Caution: in the VVL case, it only correponds to the baroclinic mass transport.
         IF( ln_zad_Aimp ) THEN
            DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z3d(ji,jj,jk) = rho0 * e1e2t(ji,jj) * ( ww(ji,jj,jk) + wi(ji,jj,jk) )
            END DO   ;   END DO   ;   END DO
         ELSE
            DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z3d(ji,jj,jk) = rho0 * e1e2t(ji,jj) * ww(ji,jj,jk)
            END DO   ;   END DO   ;   END DO
         ENDIF
         CALL iom_put( "w_masstr" , z3d )  
         IF( iom_use('w_masstr2') )   CALL iom_put( "w_masstr2", z3d * z3d )
      ENDIF

      CALL iom_put( "avt" , avt )                  ! T vert. eddy diff. coef.
      CALL iom_put( "avs" , avs )                  ! S vert. eddy diff. coef.
      CALL iom_put( "avm" , avm )                  ! T vert. eddy visc. coef.

      IF( iom_use('logavt') )   CALL iom_put( "logavt", LOG( MAX( 1.e-20_wp, avt(:,:,:) ) ) )
      IF( iom_use('logavs') )   CALL iom_put( "logavs", LOG( MAX( 1.e-20_wp, avs(:,:,:) ) ) )

      IF ( iom_use("sssgrad") .OR. iom_use("sssgrad2") ) THEN
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)                       ! sss gradient
            zztmp  = ts(ji,jj,1,jp_sal,Kmm)
            zztmpx = (ts(ji+1,jj,1,jp_sal,Kmm) - zztmp) * r1_e1u(ji,jj) + (zztmp - ts(ji-1,jj  ,1,jp_sal,Kmm)) * r1_e1u(ji-1,jj)
            zztmpy = (ts(ji,jj+1,1,jp_sal,Kmm) - zztmp) * r1_e2v(ji,jj) + (zztmp - ts(ji  ,jj-1,1,jp_sal,Kmm)) * r1_e2v(ji,jj-1)
            z2d(ji,jj) = 0.25_wp * ( zztmpx * zztmpx + zztmpy * zztmpy )   &
               &                 * umask(ji,jj,1) * umask(ji-1,jj,1) * vmask(ji,jj,1) * vmask(ji,jj-1,1)
         END DO   ;   END DO
         CALL iom_put( "sssgrad2",  z2d )          ! square of module of sss gradient
         IF ( iom_use("sssgrad") ) THEN
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               z2d(ji,jj) = SQRT( z2d(ji,jj) )
            END DO   ;   END DO
            CALL iom_put( "sssgrad",  z2d )        ! module of sss gradient
         ENDIF
      ENDIF
         
      IF ( iom_use("sstgrad") .OR. iom_use("sstgrad2") ) THEN
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)                       ! sst gradient
            zztmp  = ts(ji,jj,1,jp_tem,Kmm)
            zztmpx = ( ts(ji+1,jj,1,jp_tem,Kmm) - zztmp ) * r1_e1u(ji,jj) + ( zztmp - ts(ji-1,jj  ,1,jp_tem,Kmm) ) * r1_e1u(ji-1,jj)
            zztmpy = ( ts(ji,jj+1,1,jp_tem,Kmm) - zztmp ) * r1_e2v(ji,jj) + ( zztmp - ts(ji  ,jj-1,1,jp_tem,Kmm) ) * r1_e2v(ji,jj-1)
            z2d(ji,jj) = 0.25_wp * ( zztmpx * zztmpx + zztmpy * zztmpy )   &
               &                 * umask(ji,jj,1) * umask(ji-1,jj,1) * vmask(ji,jj,1) * vmask(ji,jj-1,1)
         END DO   ;   END DO
         CALL iom_put( "sstgrad2",  z2d )          ! square of module of sst gradient
         IF ( iom_use("sstgrad") ) THEN
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               z2d(ji,jj) = SQRT( z2d(ji,jj) )
            END DO   ;   END DO
            CALL iom_put( "sstgrad",  z2d )        ! module of sst gradient
         ENDIF
      ENDIF
         
      ! heat and salt contents
      IF( iom_use("heatc") ) THEN
         z2d(:,:)  = 0._wp 
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z2d(ji,jj) = z2d(ji,jj) + e3t_0(ji,jj,jk) * ts(ji,jj,jk,jp_tem,Kmm) * tmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "heatc", rho0_rcp * z2d )   ! vertically integrated heat content (J/m2)
      ENDIF

      IF( iom_use("saltc") ) THEN
         z2d(:,:)  = 0._wp 
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z2d(ji,jj) = z2d(ji,jj) + e3t_0(ji,jj,jk) * ts(ji,jj,jk,jp_sal,Kmm) * tmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "saltc", rho0 * z2d )       ! vertically integrated salt content (PSU*kg/m2)
      ENDIF
      !
      IF( iom_use("salt2c") ) THEN
         z2d(:,:)  = 0._wp 
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z2d(ji,jj) = z2d(ji,jj) + e3t_0(ji,jj,jk) * ts(ji,jj,jk,jp_sal,Kmm) * ts(ji,jj,jk,jp_sal,Kmm) * tmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "salt2c", rho0 * z2d )      ! vertically integrated square of salt content (PSU2*kg/m2)
      ENDIF
      !
      IF ( iom_use("ke") .OR. iom_use("ke_int") ) THEN
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            zztmpx = uu(ji-1,jj  ,jk,Kmm) + uu(ji,jj,jk,Kmm)
            zztmpy = vv(ji  ,jj-1,jk,Kmm) + vv(ji,jj,jk,Kmm)
            z3d(ji,jj,jk) = 0.25_wp * ( zztmpx*zztmpx + zztmpy*zztmpy )
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "ke", z3d )                 ! kinetic energy

         z2d(:,:)  = 0._wp 
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z2d(ji,jj) = z2d(ji,jj) + e3t_0(ji,jj,jk) * z3d(ji,jj,jk) * e1e2t(ji,jj) * tmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "ke_int", z2d )             ! vertically integrated kinetic energy
      ENDIF
      !
      IF ( iom_use("sKE") ) THEN                   ! surface kinetic energy at T point
         z2d(:,:) = 0._wp
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            z2d(ji,jj) = 0.25_wp * ( uu(ji  ,jj,1,Kmm) * uu(ji  ,jj,1,Kmm) * e1e2u(ji  ,jj) * e3u_0(ji  ,jj,1)  &
               &                   + uu(ji-1,jj,1,Kmm) * uu(ji-1,jj,1,Kmm) * e1e2u(ji-1,jj) * e3u_0(ji-1,jj,1)  &
               &                   + vv(ji,jj  ,1,Kmm) * vv(ji,jj  ,1,Kmm) * e1e2v(ji,jj  ) * e3v_0(ji,jj  ,1)  & 
               &                   + vv(ji,jj-1,1,Kmm) * vv(ji,jj-1,1,Kmm) * e1e2v(ji,jj-1) * e3v_0(ji,jj-1,1)  )  &
               &                 * r1_e1e2t(ji,jj) / e3t_0(ji,jj,1) * ssmask(ji,jj)
         END DO   ;   END DO
         IF ( iom_use("sKE" ) )  CALL iom_put( "sKE" , z2d )   
      ENDIF
      !    
      IF ( iom_use("ssKEf") ) THEN                 ! surface kinetic energy at F point
         z2d(:,:) = 0._wp                          ! CAUTION : only valid in SWE, not with bathymetry
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            z2d(ji,jj) = 0.25_wp * ( uu(ji,jj  ,1,Kmm) * uu(ji,jj  ,1,Kmm) * e1e2u(ji,jj  ) * e3u_0(ji,jj  ,1)  &
               &                   + uu(ji,jj+1,1,Kmm) * uu(ji,jj+1,1,Kmm) * e1e2u(ji,jj+1) * e3u_0(ji,jj+1,1)  &
               &                   + vv(ji  ,jj,1,Kmm) * vv(ji,jj  ,1,Kmm) * e1e2v(ji  ,jj) * e3v_0(ji  ,jj,1)  & 
               &                   + vv(ji+1,jj,1,Kmm) * vv(ji+1,jj,1,Kmm) * e1e2v(ji+1,jj) * e3v_0(ji+1,jj,1)  )  &
               &                 * r1_e1e2f(ji,jj) / e3f_0(ji,jj,1) * ssfmask(ji,jj)
         END DO   ;   END DO
         CALL iom_put( "ssKEf", z2d )                     
      ENDIF
      !
      CALL iom_put( "hdiv", hdiv )                 ! Horizontal divergence
      !
      IF( iom_use("u_masstr") .OR. iom_use("u_masstr_vint") .OR. iom_use("u_heattr") .OR. iom_use("u_salttr") ) THEN
         
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z3d(ji,jj,jk) = rho0 * uu(ji,jj,jk,Kmm) * e2u(ji,jj) * e3u_0(ji,jj,jk) * umask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "u_masstr"     , z3d )      ! mass transport in i-direction
         
         IF( iom_use("u_masstr_vint") ) THEN
            z2d(:,:) = 0._wp 
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z2d(ji,jj) = z2d(ji,jj) + z3d(ji,jj,jk)
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "u_masstr_vint", z2d )   ! mass transport in i-direction vertical sum
         ENDIF
         IF( iom_use("u_heattr") ) THEN
            z2d(:,:) = 0._wp 
            zztmp = 0.5_wp * rcp
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z2d(ji,jj) = z2d(ji,jj) + zztmp * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji+1,jj,jk,jp_tem,Kmm) )
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "u_heattr", z2d )        ! heat transport in i-direction
         ENDIF
         IF( iom_use("u_salttr") ) THEN
            z2d(:,:) = 0._wp 
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z2d(ji,jj) = z2d(ji,jj) +   0.5 * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji+1,jj,jk,jp_sal,Kmm) )
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "u_salttr", z2d )        ! heat transport in i-direction
         ENDIF
         
      ENDIF
      
      IF( iom_use("v_masstr") .OR. iom_use("v_heattr") .OR. iom_use("v_salttr") ) THEN
         
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z3d(ji,jj,jk) = rho0 * vv(ji,jj,jk,Kmm) * e1v(ji,jj) * e3v_0(ji,jj,jk) * vmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "v_masstr", z3d )           ! mass transport in j-direction
         
         IF( iom_use("v_heattr") ) THEN
            z2d(:,:) = 0._wp
            zztmp = 0.5_wp * rcp
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z2d(ji,jj) = z2d(ji,jj) + zztmp * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_tem,Kmm) + ts(ji,jj+1,jk,jp_tem,Kmm) )
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "v_heattr", z2d )        !  heat transport in j-direction
         ENDIF
         IF( iom_use("v_salttr") ) THEN
            z2d(:,:) = 0._wp 
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               z2d(ji,jj) = z2d(ji,jj) +   0.5 * z3d(ji,jj,jk) * ( ts(ji,jj,jk,jp_sal,Kmm) + ts(ji,jj+1,jk,jp_sal,Kmm) )
            END DO   ;   END DO   ;   END DO
            CALL iom_put( "v_salttr", z2d )        !  heat transport in j-direction
         ENDIF

      ENDIF

      IF( iom_use("tosmint") ) THEN
         z2d(:,:) = 0._wp
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z2d(ji,jj) = z2d(ji,jj) + rho0 * e3t_0(ji,jj,jk) * ts(ji,jj,jk,jp_tem,Kmm)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "tosmint", z2d )            ! Vertical integral of temperature
      ENDIF
      IF( iom_use("somint") ) THEN
         z2d(:,:) = 0._wp
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            z2d(ji,jj) = z2d(ji,jj) + rho0 * e3t_0(ji,jj,jk) * ts(ji,jj,jk,jp_sal,Kmm)
         END DO   ;   END DO   ;   END DO
         CALL iom_put( "somint", z2d )             ! Vertical integral of salinity
      ENDIF

      CALL iom_put( "bn2", rn2 )                   ! Brunt-Vaisala buoyancy frequency (N^2)
      
      IF (ln_dia25h)   CALL dia_25h( kt, Kmm )     ! 25h averaging
      
      ! Output of surface vorticity terms
      !
      CALL iom_put( "ssplavor", ff_f )             ! planetary vorticity ( f )
      !
      IF ( iom_use("ssrelvor")    .OR. iom_use("ssEns")    .OR.   &
         & iom_use("ssrelpotvor") .OR. iom_use("ssabspotvor") ) THEN
         !
         z2d(:,:) = 0._wp 
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            z2d(ji,jj) = (   e2v(ji+1,jj  ) * vv(ji+1,jj  ,1,Kmm) - e2v(ji,jj) * vv(ji,jj,1,Kmm)    &
            &              - e1u(ji  ,jj+1) * uu(ji  ,jj+1,1,Kmm) + e1u(ji,jj) * uu(ji,jj,1,Kmm)  ) * r1_e1e2f(ji,jj)
         END DO   ;   END DO
         CALL iom_put( "ssrelvor", z2d )           ! relative vorticity ( zeta ) 
         !
         IF ( iom_use("ssEns") .OR. iom_use("ssrelpotvor") .OR. iom_use("ssabspotvor") ) THEN
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)  
               ze3 = (  e3t_0(ji,jj+1,1) * e1e2t(ji,jj+1) + e3t_0(ji+1,jj+1,1) * e1e2t(ji+1,jj+1)    &
                  &    + e3t_0(ji,jj  ,1) * e1e2t(ji,jj  ) + e3t_0(ji+1,jj  ,1) * e1e2t(ji+1,jj  )  ) * r1_e1e2f(ji,jj)
               IF( ze3 /= 0._wp ) THEN   ;   ze3 = 4._wp / ze3
               ELSE                      ;   ze3 = 0._wp
               ENDIF
               z2d(ji,jj) = ze3 * z2d(ji,jj) 
            END DO   ;   END DO
            CALL iom_put( "ssrelpotvor", z2d )     ! relative potential vorticity (zeta/h)
            !
            IF ( iom_use("ssEns") .OR. iom_use("ssabspotvor") ) THEN
               DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
                  ze3 = (  e3t_0(ji,jj+1,1) * e1e2t(ji,jj+1) + e3t_0(ji+1,jj+1,1) * e1e2t(ji+1,jj+1)    &
                     &    + e3t_0(ji,jj  ,1) * e1e2t(ji,jj  ) + e3t_0(ji+1,jj  ,1) * e1e2t(ji+1,jj  )  ) * r1_e1e2f(ji,jj)
                  IF( ze3 /= 0._wp ) THEN   ;   ze3 = 4._wp / ze3
                  ELSE                      ;   ze3 = 0._wp
                  ENDIF
                  z2d(ji,jj) = ze3 * ff_f(ji,jj) + z2d(ji,jj) 
               END DO   ;   END DO
               CALL iom_put( "ssabspotvor", z2d )  ! absolute potential vorticity ( q )
               !
               IF ( iom_use("ssEns") ) THEN
                  DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)  
                     z2d(ji,jj) = 0.5_wp * z2d(ji,jj) * z2d(ji,jj) 
                  END DO   ;   END DO
                  CALL iom_put( "ssEns", z2d )     ! potential enstrophy ( 1/2*q2 )
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      IF( ln_timing )   CALL timing_stop('dia_wri')
      !
   END SUBROUTINE dia_wri


   SUBROUTINE dia_wri_state( Kmm, cdfile_name )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE dia_wri_state  ***
      !!        
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains 
      !!      the instantaneous ocean state and forcing fields.
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! ** Method  :   NetCDF files using ioipsl
      !!      File 'output.init.nc'  is created if ninist = 1 (namelist)
      !!      File 'output.abort.nc' is created in case of abnormal job end
      !!----------------------------------------------------------------------
      INTEGER           , INTENT( in ) ::   Kmm              ! time level index
      CHARACTER (len=* ), INTENT( in ) ::   cdfile_name      ! name of the file created
      !!
      INTEGER ::   ji, jj, jk       ! dummy loop indices
      INTEGER ::   inum
      REAL(wp), DIMENSION(jpi,jpj)     :: z2d      
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: z3d      
      !!----------------------------------------------------------------------
      ! 
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'dia_wri_state : single instantaneous ocean state'
         WRITE(numout,*) '~~~~~~~~~~~~~   and forcing fields file created '
         WRITE(numout,*) '                and named :', cdfile_name, '...nc'
      ENDIF 
      !
      CALL iom_open( TRIM(cdfile_name), inum, ldwrt = .TRUE. )
      !
      CALL iom_rstput( 0, 0, inum, 'votemper', ts(:,:,:,jp_tem,Kmm) )    ! now temperature
      CALL iom_rstput( 0, 0, inum, 'vosaline', ts(:,:,:,jp_sal,Kmm) )    ! now salinity
      CALL iom_rstput( 0, 0, inum, 'sossheig', ssh(:,:,Kmm)         )    ! sea surface height
      CALL iom_rstput( 0, 0, inum, 'vozocrtx', uu(:,:,:,Kmm)        )    ! now i-velocity
      CALL iom_rstput( 0, 0, inum, 'vomecrty', vv(:,:,:,Kmm)        )    ! now j-velocity
      IF( ln_zad_Aimp ) THEN
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
           z3d(ji,jj,jk) = ww(ji,jj,jk) + wi(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         CALL iom_rstput( 0, 0, inum, 'vovecrtz', z3d            )    ! now k-velocity
      ELSE
         CALL iom_rstput( 0, 0, inum, 'vovecrtz', ww             )    ! now k-velocity
      ENDIF
      CALL iom_rstput( 0, 0, inum, 'risfdep', risfdep            )
      CALL iom_rstput( 0, 0, inum, 'ht'     , ht_0(:,:)            )    ! now water column height
      !
      IF ( ln_isf ) THEN
         IF (ln_isfcav_mlt) THEN
            CALL iom_rstput( 0, 0, inum, 'fwfisf_cav', fwfisf_cav          )
            CALL iom_rstput( 0, 0, inum, 'rhisf_cav_tbl', rhisf_tbl_cav    )
            CALL iom_rstput( 0, 0, inum, 'rfrac_cav_tbl', rfrac_tbl_cav    )
            CALL iom_rstput( 0, 0, inum, 'misfkb_cav', REAL(misfkb_cav,wp) )
            CALL iom_rstput( 0, 0, inum, 'misfkt_cav', REAL(misfkt_cav,wp) )
            CALL iom_rstput( 0, 0, inum, 'mskisf_cav', REAL(mskisf_cav,wp), ktype = jp_i1 )
         END IF
         IF (ln_isfpar_mlt) THEN
            CALL iom_rstput( 0, 0, inum, 'isfmsk_par', REAL(mskisf_par,wp) )
            CALL iom_rstput( 0, 0, inum, 'fwfisf_par', fwfisf_par          )
            CALL iom_rstput( 0, 0, inum, 'rhisf_par_tbl', rhisf_tbl_par    )
            CALL iom_rstput( 0, 0, inum, 'rfrac_par_tbl', rfrac_tbl_par    )
            CALL iom_rstput( 0, 0, inum, 'misfkb_par', REAL(misfkb_par,wp) )
            CALL iom_rstput( 0, 0, inum, 'misfkt_par', REAL(misfkt_par,wp) )
            CALL iom_rstput( 0, 0, inum, 'mskisf_par', REAL(mskisf_par,wp), ktype = jp_i1 )
         END IF
      END IF
      !
      IF( ALLOCATED(ahtu) ) THEN
         CALL iom_rstput( 0, 0, inum,  'ahtu', ahtu              )    ! aht at u-point
         CALL iom_rstput( 0, 0, inum,  'ahtv', ahtv              )    ! aht at v-point
      ENDIF
      IF( ALLOCATED(ahmt) ) THEN
         CALL iom_rstput( 0, 0, inum,  'ahmt', ahmt              )    ! ahmt at u-point
         CALL iom_rstput( 0, 0, inum,  'ahmf', ahmf              )    ! ahmf at v-point
      ENDIF
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
         z2d(ji,jj) = emp(ji,jj) - rnf(ji,jj)
      END DO   ;   END DO
      CALL iom_rstput( 0, 0, inum, 'sowaflup', z2d               )    ! freshwater budget
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
         z2d(ji,jj) = qsr(ji,jj) + qns(ji,jj)
      END DO   ;   END DO
      CALL iom_rstput( 0, 0, inum, 'sohefldo', z2d               )    ! total heat flux
      CALL iom_rstput( 0, 0, inum, 'soshfldo', qsr               )    ! solar heat flux
      CALL iom_rstput( 0, 0, inum, 'soicecov', fr_i              )    ! ice fraction
      CALL iom_rstput( 0, 0, inum, 'sozotaux', utau              )    ! i-wind stress
      CALL iom_rstput( 0, 0, inum, 'sometauy', vtau              )    ! j-wind stress
      IF(  .NOT.ln_linssh  ) THEN
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
           z3d(ji,jj,jk) = gdept_0(ji,jj,jk)   ! 3D workspace for qco substitution
         END DO   ;   END DO   ;   END DO
         CALL iom_rstput( 0, 0, inum, 'vovvldep', z3d            )    !  T-cell depth 
         DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
           z3d(ji,jj,jk) = e3t_0(ji,jj,jk)     ! 3D workspace for qco substitution
         END DO   ;   END DO   ;   END DO
         CALL iom_rstput( 0, 0, inum, 'vovvle3t', z3d            )    !  T-cell thickness  
      END IF
      IF( ln_wave .AND. ln_sdw ) THEN
         CALL iom_rstput( 0, 0, inum, 'sdzocrtx', usd            )    ! now StokesDrift i-velocity
         CALL iom_rstput( 0, 0, inum, 'sdmecrty', vsd            )    ! now StokesDrift j-velocity
         CALL iom_rstput( 0, 0, inum, 'sdvecrtz', wsd            )    ! now StokesDrift k-velocity
      ENDIF
      IF ( ln_abl ) THEN
         CALL iom_rstput ( 0, 0, inum, "uz1_abl",   u_abl(:,:,2,nt_a  ) )   ! now first level i-wind
         CALL iom_rstput ( 0, 0, inum, "vz1_abl",   v_abl(:,:,2,nt_a  ) )   ! now first level j-wind
         CALL iom_rstput ( 0, 0, inum, "tz1_abl",  tq_abl(:,:,2,nt_a,1) )   ! now first level temperature
         CALL iom_rstput ( 0, 0, inum, "qz1_abl",  tq_abl(:,:,2,nt_a,2) )   ! now first level humidity
      ENDIF
      IF( ln_zdfosm ) THEN
         CALL iom_rstput( 0, 0, inum, 'hbl', hbl*tmask(:,:,1)  )      ! now boundary-layer depth
         CALL iom_rstput( 0, 0, inum, 'hml', hml*tmask(:,:,1)  )      ! now mixed-layer depth
         CALL iom_rstput( 0, 0, inum, 'avt_k', avt_k*wmask     )      ! w-level diffusion
         CALL iom_rstput( 0, 0, inum, 'avm_k', avm_k*wmask     )      ! now w-level viscosity
         CALL iom_rstput( 0, 0, inum, 'ghamt', ghamt*wmask     )      ! non-local t forcing
         CALL iom_rstput( 0, 0, inum, 'ghams', ghams*wmask     )      ! non-local s forcing
         CALL iom_rstput( 0, 0, inum, 'ghamu', ghamu*umask     )      ! non-local u forcing
         CALL iom_rstput( 0, 0, inum, 'ghamv', ghamv*vmask     )      ! non-local v forcing
         IF( ln_osm_mle ) THEN
            CALL iom_rstput( 0, 0, inum, 'hmle', hmle*tmask(:,:,1)  ) ! now transition-layer depth
         END IF
      ENDIF
      !
      CALL iom_close( inum )
      ! 
   END SUBROUTINE dia_wri_state

   !!======================================================================
END MODULE diawri

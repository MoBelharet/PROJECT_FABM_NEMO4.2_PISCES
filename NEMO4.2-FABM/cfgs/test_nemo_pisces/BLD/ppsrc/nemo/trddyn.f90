










MODULE trddyn
   !!======================================================================
   !!                       ***  MODULE  trddyn  ***
   !! Ocean diagnostics:  ocean dynamic trends
   !!=====================================================================
   !! History :  3.5  !  2012-02  (G. Madec) creation from trdmod: split DYN and TRA trends
   !!                                        and manage  3D trends output for U, V, and KE
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   trd_dyn       : manage the type of momentum trend diagnostics (3D I/O, domain averaged, KE)
   !!   trd_dyn_iom   : output 3D momentum and/or tracer trends using IOM
   !!   trd_dyn_init  : initialization step
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE dom_oce        ! ocean space and time domain variables
   USE phycst         ! physical constants
   USE sbc_oce        ! surface boundary condition: ocean
   USE zdf_oce        ! ocean vertical physics: variables
!!gm   USE zdfdrg         ! ocean vertical physics: bottom friction
   USE trd_oce        ! trends: ocean variables
   USE trdken         ! trends: Kinetic ENergy 
   USE trdglo         ! trends: global domain averaged
   USE trdvor         ! trends: vertical averaged vorticity 
   USE trdmxl         ! trends: mixed layer averaged 
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! lateral boundary condition 
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC trd_dyn        ! called by all dynXXX modules

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
   !! $Id: trddyn.F90 14433 2021-02-11 08:06:49Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trd_dyn( putrd, pvtrd, ktrd, kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mod  ***
      !! 
      !! ** Purpose :   Dispatch momentum trend computation, e.g. 3D output, 
      !!              integral constraints, barotropic vorticity, kinetic enrgy, 
      !!              and/or mixed layer budget.
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(:,:,:), INTENT(inout) ::   putrd, pvtrd   ! U and V trends 
      INTEGER                   , INTENT(in   ) ::   ktrd           ! trend index
      INTEGER                   , INTENT(in   ) ::   kt             ! time step
      INTEGER                   , INTENT(in   ) ::   Kmm            ! time level index
      !!----------------------------------------------------------------------
      !
      putrd(:,:,:) = putrd(:,:,:) * umask(:,:,:)                       ! mask the trends
      pvtrd(:,:,:) = pvtrd(:,:,:) * vmask(:,:,:)
      !

!!gm NB : here a lbc_lnk should probably be added

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !   3D output of momentum and/or tracers trends using IOM interface
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_dyn_trd )   CALL trd_dyn_iom( putrd, pvtrd, ktrd, kt, Kmm )
         
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  Integral Constraints Properties for momentum and/or tracers trends
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF( ln_glo_trd )   CALL trd_glo( putrd, pvtrd, ktrd, 'DYN', kt, Kmm )

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  Kinetic Energy trends
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF( ln_KE_trd  )   CALL trd_ken( putrd, pvtrd, ktrd, kt, Kmm )

      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      !  Vorticity trends
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      IF( ln_vor_trd )   CALL trd_vor( putrd, pvtrd, ktrd, kt, Kmm )

      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      !  Mixed layer trends for active tracers
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!!gm      IF( ln_dyn_mxl )   CALL trd_mxl_dyn   
      !
   END SUBROUTINE trd_dyn


   SUBROUTINE trd_dyn_iom( putrd, pvtrd, ktrd, kt, Kmm )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE trd_dyn_iom  ***
      !! 
      !! ** Purpose :   output 3D trends using IOM
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(:,:,:), INTENT(inout) ::   putrd, pvtrd   ! U and V trends
      INTEGER                   , INTENT(in   ) ::   ktrd           ! trend index
      INTEGER                   , INTENT(in   ) ::   kt             ! time step
      INTEGER                   , INTENT(in   ) ::   Kmm            ! time level index
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      INTEGER ::   ikbu, ikbv   ! local integers
      REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   z2dx, z2dy   ! 2D workspace 
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   z3dx, z3dy   ! 3D workspace 
      !!----------------------------------------------------------------------
      !
      SELECT CASE( ktrd )
      CASE( jpdyn_hpg )   ;   CALL iom_put( "utrd_hpg", putrd )    ! hydrostatic pressure gradient
                              CALL iom_put( "vtrd_hpg", pvtrd )
      CASE( jpdyn_spg )   ;   CALL iom_put( "utrd_spg", putrd )    ! surface pressure gradient
                              CALL iom_put( "vtrd_spg", pvtrd )
      CASE( jpdyn_pvo )   ;   CALL iom_put( "utrd_pvo", putrd )    ! planetary vorticity
                              CALL iom_put( "vtrd_pvo", pvtrd )
      CASE( jpdyn_rvo )   ;   CALL iom_put( "utrd_rvo", putrd )    ! relative  vorticity     (or metric term)
                              CALL iom_put( "vtrd_rvo", pvtrd )
      CASE( jpdyn_keg )   ;   CALL iom_put( "utrd_keg", putrd )    ! Kinetic Energy gradient (or had)
                              CALL iom_put( "vtrd_keg", pvtrd )
                              ALLOCATE( z3dx(jpi,jpj,jpk) , z3dy(jpi,jpj,jpk) )
                              z3dx(:,:,:) = 0._wp                  ! U.dxU & V.dyV (approximation)
                              z3dy(:,:,:) = 0._wp
                              DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)   ! no mask as un,vn are masked
                                 z3dx(ji,jj,jk) = uu(ji,jj,jk,Kmm) * ( uu(ji+1,jj,jk,Kmm) - uu(ji-1,jj,jk,Kmm) ) / ( 2._wp * e1u(ji,jj) )
                                 z3dy(ji,jj,jk) = vv(ji,jj,jk,Kmm) * ( vv(ji,jj+1,jk,Kmm) - vv(ji,jj-1,jk,Kmm) ) / ( 2._wp * e2v(ji,jj) )
                              END DO   ;   END DO   ;   END DO
                              CALL lbc_lnk( 'trddyn', z3dx, 'U', -1.0_wp, z3dy, 'V', -1.0_wp )
                              CALL iom_put( "utrd_udx", z3dx  )
                              CALL iom_put( "vtrd_vdy", z3dy  )
                              DEALLOCATE( z3dx , z3dy )
      CASE( jpdyn_zad )   ;   CALL iom_put( "utrd_zad", putrd )    ! vertical advection
                              CALL iom_put( "vtrd_zad", pvtrd )
      CASE( jpdyn_ldf )   ;   CALL iom_put( "utrd_ldf", putrd )    ! lateral  diffusion
                              CALL iom_put( "vtrd_ldf", pvtrd )
      CASE( jpdyn_zdf )   ;   CALL iom_put( "utrd_zdf", putrd )    ! vertical diffusion 
                              CALL iom_put( "vtrd_zdf", pvtrd )
                              !
                              !                                    ! wind stress trends
                              ALLOCATE( z2dx(jpi,jpj) , z2dy(jpi,jpj) )
                              z2dx(:,:) = ( utau_b(:,:) + utau(:,:) ) / ( e3u_0(:,:,1) * rho0 )
                              z2dy(:,:) = ( vtau_b(:,:) + vtau(:,:) ) / ( e3v_0(:,:,1) * rho0 )
                              CALL iom_put( "utrd_tau", z2dx )
                              CALL iom_put( "vtrd_tau", z2dy )
                              DEALLOCATE( z2dx , z2dy )
!!gm  to be changed : computation should be done in dynzdf.F90
!!gm                + missing the top friction 
!                              !                                    ! bottom stress tends (implicit case)
!                              IF( ln_drgimp ) THEN
!                                 ALLOCATE( z3dx(jpi,jpj,jpk) , z3dy(jpi,jpj,jpk) )
!                             z3dx(:,:,:) = 0._wp   ;   z3dy(:,:,:) = 0._wp  ! after velocity known (now filed at this stage)
!                            DO jk = 1, jpkm1
!                                    DO jj = 2, jpjm1
!                                       DO ji = 2, jpim1
!                                      ikbu = mbku(ji,jj)          ! deepest ocean u- & v-levels
!                                          ikbv = mbkv(ji,jj)
!                                          z3dx(ji,jj,jk) = 0.5 * ( rCdU_bot(ji+1,jj) + rCdU_bot(ji,jj) ) & 
!                                               &         * uu(ji,jj,ikbu,Kmm) / e3u_0(ji,jj,ikbu)
!                                          z3dy(ji,jj,jk) = 0.5 * ( rCdU_bot(ji,jj+1) + rCdU_bot(ji,jj) ) &
!                                               &         * vv(ji,jj,ikbv,Kmm) / e3v_0(ji,jj,ikbv)
!                                    END DO
!                                 END DO
!                              END DO
!                              CALL lbc_lnk( 'trddyn', z3dx, 'U', -1.0_wp, z3dy, 'V', -1.0_wp )
!                              CALL iom_put( "utrd_bfr", z3dx )
!                              CALL iom_put( "vtrd_bfr", z3dy )
!                                 DEALLOCATE( z3dx , z3dy )
!                              ENDIF
!!gm end
      CASE( jpdyn_bfr )       ! called if ln_drgimp=F
                              CALL iom_put( "utrd_bfr", putrd )    ! bottom friction (explicit case)
                              CALL iom_put( "vtrd_bfr", pvtrd )
      CASE( jpdyn_atf )   ;   CALL iom_put( "utrd_atf", putrd )        ! asselin filter trends 
                              CALL iom_put( "vtrd_atf", pvtrd )
      END SELECT
      !
   END SUBROUTINE trd_dyn_iom

   !!======================================================================
END MODULE trddyn

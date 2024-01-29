










MODULE sshwzv   
   !!==============================================================================
   !!                       ***  MODULE  sshwzv  ***
   !! Ocean dynamics : sea surface height and vertical velocity
   !!==============================================================================
   !! History :  3.1  !  2009-02  (G. Madec, M. Leclair)  Original code
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  modified LF-RA 
   !!             -   !  2010-05  (K. Mogensen, A. Weaver, M. Martin, D. Lea)  Assimilation interface
   !!             -   !  2010-09  (D.Storkey and E.O'Dea)  bug fixes for BDY module
   !!            3.3  !  2011-10  (M. Leclair)  split former ssh_wzv routine and remove all vvl related work
   !!            4.0  !  2018-12  (A. Coward)  add mixed implicit/explicit advection
   !!            4.1  !  2019-08  (A. Coward, D. Storkey)  Rename ssh_nxt -> ssh_atf. Now only does time filtering.
   !!             -   !  2020-08  (S. Techene, G. Madec)  add here ssh initiatlisation
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ssh_nxt       : after ssh
   !!   ssh_atf       : time filter the ssh arrays
   !!   wzv           : compute now vertical velocity
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE isf_oce        ! ice shelf
   USE dom_oce        ! ocean space and time domain variables 
   USE sbc_oce        ! surface boundary condition: ocean
   USE domvvl         ! Variable volume
   USE divhor         ! horizontal divergence
   USE phycst         ! physical constants
   USE bdy_oce , ONLY : ln_bdy, bdytmask   ! Open BounDarY
   USE bdydyn2d       ! bdy_ssh routine
   USE wet_dry        ! Wetting/Drying flux limiting
   !
   USE iom 
   USE in_out_manager ! I/O manager
   USE restart        ! only for lrst_oce
   USE prtctl         ! Print control
   USE lbclnk         ! ocean lateral boundary condition (or mpp link)
   USE lib_mpp        ! MPP library
   USE timing         ! Timing
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC   ssh_nxt        ! called by step.F90
   PUBLIC   wzv            ! called by step.F90
   PUBLIC   wAimp          ! called by step.F90
   PUBLIC   ssh_atf        ! called by step.F90

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
   !! $Id: sshwzv.F90 15150 2021-07-27 10:38:24Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ssh_nxt( kt, Kbb, Kmm, pssh, Kaa )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE ssh_nxt  ***
      !!                   
      !! ** Purpose :   compute the after ssh (ssh(Kaa))
      !!
      !! ** Method  : - Using the incompressibility hypothesis, the ssh increment
      !!      is computed by integrating the horizontal divergence and multiply by
      !!      by the time step.
      !!
      !! ** action  :   ssh(:,:,Kaa), after sea surface height
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt             ! time step
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm, Kaa  ! time level index
      REAL(dp), DIMENSION(jpi,jpj,jpt), INTENT(inout) ::   pssh           ! sea-surface height
      ! 
      INTEGER  ::   ji, jj, jk      ! dummy loop index
      REAL(wp) ::   zcoef   ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zhdiv   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ssh_nxt')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ssh_nxt : after sea surface height'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      zcoef = 0.5_wp * r1_rho0

      !                                           !------------------------------!
      !                                           !   After Sea Surface Height   !
      !                                           !------------------------------!

      CALL div_hor( kt, Kbb, Kmm )                     ! Horizontal divergence
      !
      zhdiv(:,:) = 0._wp
      DO jk =  1,  jpkm1  ; DO jj = ntsj-(  1), ntej+(  nn_hls) ; DO ji = ntsi-( 1), ntei+(  nn_hls)                                 ! Horizontal divergence of barotropic transports
        zhdiv(ji,jj) = zhdiv(ji,jj) + e3t_0(ji,jj,jk) * hdiv(ji,jj,jk)
      END DO   ;   END DO   ;   END DO
      !                                                ! Sea surface elevation time stepping
      ! In time-split case we need a first guess of the ssh after (using the baroclinic timestep) in order to
      ! compute the vertical velocity which can be used to compute the non-linear terms of the momentum equations.
      ! 
      DO jj = ntsj-(  1-( 1+ nn_hls )*nthb), ntej+(  nn_hls -( nn_hls + 1)*ntht) ; DO ji = ntsi-( 1-( 1+ nn_hls)*nthl), ntei+(  nn_hls-( nn_hls+ 1)*nthr)                ! Loop bounds limited by hdiv definition in div_hor
         pssh(ji,jj,Kaa) = (  pssh(ji,jj,Kbb) - rDt * ( zcoef * ( emp_b(ji,jj) + emp(ji,jj) ) + zhdiv(ji,jj) )  ) * ssmask(ji,jj)
      END DO   ;   END DO
      ! pssh must be defined everywhere (true for dyn_spg_ts, not for dyn_spg_exp)
      IF ( .NOT. ln_dynspg_ts .AND. nn_hls == 2 ) CALL lbc_lnk( 'sshwzv', pssh(:,:,Kaa), 'T', 1.0_dp )
      !
      !
      IF ( .NOT.ln_dynspg_ts ) THEN
         IF( ln_bdy ) THEN
            IF (nn_hls==1) CALL lbc_lnk( 'sshwzv', pssh(:,:,Kaa), 'T', 1.0_dp )    ! Not sure that's necessary
            CALL bdy_ssh( pssh(:,:,Kaa) )             ! Duplicate sea level across open boundaries
         ENDIF
      ENDIF
      !                                           !------------------------------!
      !                                           !           outputs            !
      !                                           !------------------------------!
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab2d_1=pssh(:,:,Kaa), clinfo1=' pssh(:,:,Kaa)  - : ', mask1=tmask )
      !
      IF( ln_timing )   CALL timing_stop('ssh_nxt')
      !
   END SUBROUTINE ssh_nxt

   
   SUBROUTINE wzv( kt, Kbb, Kmm, Kaa, pww )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE wzv  ***
      !!                   
      !! ** Purpose :   compute the now vertical velocity
      !!
      !! ** Method  : - Using the incompressibility hypothesis, the vertical 
      !!      velocity is computed by integrating the horizontal divergence  
      !!      from the bottom to the surface minus the scale factor evolution.
      !!        The boundary conditions are w=0 at the bottom (no flux) and.
      !!
      !! ** action  :   pww      : now vertical velocity
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in)    ::   kt             ! time step
      INTEGER                         , INTENT(in)    ::   Kbb, Kmm, Kaa  ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pww            ! vertical velocity at Kmm
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zhdiv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('wzv')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'wzv : now vertical velocity '
         IF(lwp) WRITE(numout,*) '~~~~~ '
         !
         pww(:,:,jpk) = 0._wp                  ! bottom boundary condition: w=0 (set once for all)
      ENDIF
      !                                           !------------------------------!
      !                                           !     Now Vertical Velocity    !
      !                                           !------------------------------!
      !
      !                                               !===============================!
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN      !==  z_tilde and layer cases  ==!
         !                                            !===============================!
         ALLOCATE( zhdiv(jpi,jpj,jpk) ) 
         !
         DO jk = 1, jpkm1
            ! horizontal divergence of thickness diffusion transport ( velocity multiplied by e3t)
            ! - ML - note: computation already done in dom_vvl_sf_nxt. Could be optimized (not critical and clearer this way)
            DO jj = ntsj-( nn_hls-1), ntej+( nn_hls ) ; DO ji = ntsi-( nn_hls-1), ntei+( nn_hls)
               zhdiv(ji,jj,jk) = r1_e1e2t(ji,jj) * ( un_td(ji,jj,jk) - un_td(ji-1,jj,jk) + vn_td(ji,jj,jk) - vn_td(ji,jj-1,jk) )
            END DO   ;   END DO
         END DO
         IF( nn_hls == 1)   CALL lbc_lnk('sshwzv', zhdiv, 'T', 1.0_wp)  ! - ML - Perhaps not necessary: not used for horizontal "connexions"
         !                             ! Is it problematic to have a wrong vertical velocity in boundary cells?
         !                             ! Same question holds for hdiv. Perhaps just for security
         !                             ! clem: yes it is a problem because ww is used in many other places where we need the halos
         !
         DO jk =  jpkm1,  1,  -1  ; DO jj = ntsj-(  nn_hls-1), ntej+(  nn_hls) ; DO ji = ntsi-( nn_hls-1), ntei+(  nn_hls)     ! integrate from the bottom the hor. divergence
            ! computation of w
            pww(ji,jj,jk) = pww(ji,jj,jk+1) - (   e3t_0(ji,jj,jk) * hdiv(ji,jj,jk)   &
               &                                  +                  zhdiv(ji,jj,jk)   &
               &                                  + r1_Dt * (  e3t_0(ji,jj,jk)       &
               &                                             - e3t_0(ji,jj,jk) )   ) * tmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         !          IF( ln_vvl_layer ) pww(:,:,:) = 0.e0
         DEALLOCATE( zhdiv ) 
         !                                            !=================================!
      ELSEIF( ln_linssh )   THEN                      !==  linear free surface cases  ==!
         !                                            !=================================!
         DO jk =  jpkm1,  1,  -1  ; DO jj = ntsj-(  nn_hls-1), ntej+(  nn_hls) ; DO ji = ntsi-( nn_hls-1), ntei+(  nn_hls)     ! integrate from the bottom the hor. divergence
            pww(ji,jj,jk) = pww(ji,jj,jk+1) - (  e3t_0(ji,jj,jk) * hdiv(ji,jj,jk)  ) * tmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         !                                            !==========================================!
      ELSE                                            !==  Quasi-Eulerian vertical coordinate  ==!   ('key_qco')
         !                                            !==========================================!
         DO jk =  jpkm1,  1,  -1  ; DO jj = ntsj-(  nn_hls-1), ntej+(  nn_hls) ; DO ji = ntsi-( nn_hls-1), ntei+(  nn_hls)     ! integrate from the bottom the hor. divergence
            pww(ji,jj,jk) = pww(ji,jj,jk+1) - (  e3t_0(ji,jj,jk) * hdiv(ji,jj,jk)    &
               &                                 + r1_Dt * (  e3t_0(ji,jj,jk)        &
               &                                            - e3t_0(ji,jj,jk)  )   ) * tmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
      ENDIF

      IF( ln_bdy ) THEN
         DO jk = 1, jpkm1
            DO jj = ntsj-( nn_hls-1), ntej+( nn_hls ) ; DO ji = ntsi-( nn_hls-1), ntei+( nn_hls)
               pww(ji,jj,jk) = pww(ji,jj,jk) * bdytmask(ji,jj)
            END DO   ;   END DO
         END DO
      ENDIF
      !
      !
      IF( ln_timing )   CALL timing_stop('wzv')
      !
   END SUBROUTINE wzv


   SUBROUTINE ssh_atf( kt, Kbb, Kmm, Kaa, pssh )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ssh_atf  ***
      !!
      !! ** Purpose :   Apply Asselin time filter to now SSH.
      !!
      !! ** Method  : - apply Asselin time fiter to now ssh (excluding the forcing
      !!              from the filter, see Leclair and Madec 2010) and swap :
      !!                pssh(:,:,Kmm) = pssh(:,:,Kaa) + rn_atfp * ( pssh(:,:,Kbb) -2 pssh(:,:,Kmm) + pssh(:,:,Kaa) )
      !!                            - rn_atfp * rn_Dt * ( emp_b - emp ) / rho0
      !!
      !! ** action  : - pssh(:,:,Kmm) time filtered
      !!
      !! Reference  : Leclair, M., and G. Madec, 2009, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt             ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm, Kaa  ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpt), INTENT(inout) ::   pssh           ! SSH field
      !
      REAL(wp) ::   zcoef   ! local scalar
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ssh_atf')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ssh_atf : Asselin time filter of sea surface height'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !
      IF( .NOT.l_1st_euler ) THEN   ! Apply Asselin time filter on Kmm field (not on euler 1st)
         !
         IF( ln_linssh ) THEN                ! filtered "now" field
            pssh(:,:,Kmm) = pssh(:,:,Kmm) + rn_atfp * ( pssh(:,:,Kbb) - 2 * pssh(:,:,Kmm) + pssh(:,:,Kaa) )
            !
         ELSE                                ! filtered "now" field with forcing removed
            zcoef = rn_atfp * rn_Dt * r1_rho0
            pssh(:,:,Kmm) = pssh(:,:,Kmm) + rn_atfp * ( pssh(:,:,Kbb) - 2 * pssh(:,:,Kmm) + pssh(:,:,Kaa) )   &
               &                          - zcoef   * (         emp_b(:,:) -        emp(:,:)   &
               &                                              - rnf_b(:,:) +        rnf(:,:)   &
               &                                       - fwfisf_cav_b(:,:) + fwfisf_cav(:,:)   &
               &                                       - fwfisf_par_b(:,:) + fwfisf_par(:,:)   ) * ssmask(:,:)

            ! ice sheet coupling
            IF( ln_isf .AND. ln_isfcpl .AND. kt == nit000+1 )   &
               &   pssh(:,:,Kbb) = pssh(:,:,Kbb) - rn_atfp * rn_Dt * ( risfcpl_ssh(:,:) - 0._wp ) * ssmask(:,:)

         ENDIF
      ENDIF
      !
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab2d_1=pssh(:,:,Kmm), clinfo1=' atf  - pssh(:,:,Kmm): ', mask1=tmask )
      !
      IF( ln_timing )   CALL timing_stop('ssh_atf')
      !
   END SUBROUTINE ssh_atf

   
   SUBROUTINE wAimp( kt, Kmm )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE wAimp  ***
      !!                   
      !! ** Purpose :   compute the Courant number and partition vertical velocity
      !!                if a proportion needs to be treated implicitly
      !!
      !! ** Method  : - 
      !!
      !! ** action  :   ww      : now vertical velocity (to be handled explicitly)
      !!            :   wi      : now vertical velocity (for implicit treatment)
      !!
      !! Reference  : Shchepetkin, A. F. (2015): An adaptive, Courant-number-dependent
      !!              implicit scheme for vertical advection in oceanic modeling. 
      !!              Ocean Modelling, 91, 38-69.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time step
      INTEGER, INTENT(in) ::   Kmm  ! time level index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp)             ::   zCu, zcff, z1_e3t, zdt                ! local scalars
      REAL(wp) , PARAMETER ::   Cu_min = 0.15_wp                      ! local parameters
      REAL(wp) , PARAMETER ::   Cu_max = 0.30_wp                      ! local parameters
      REAL(wp) , PARAMETER ::   Cu_cut = 2._wp*Cu_max - Cu_min        ! local parameters
      REAL(wp) , PARAMETER ::   Fcu    = 4._wp*Cu_max*(Cu_max-Cu_min) ! local parameters
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('wAimp')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'wAimp : Courant number-based partitioning of now vertical velocity '
         IF(lwp) WRITE(numout,*) '~~~~~ '
      ENDIF
      !
      ! Calculate Courant numbers
      zdt = 2._wp * rn_Dt                            ! 2*rn_Dt and not rDt (for restartability)
      IF( ln_vvl_ztilde .OR. ln_vvl_layer ) THEN
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  nn_hls-1), ntej+(  nn_hls) ; DO ji = ntsi-( nn_hls-1), ntei+(  nn_hls)
            z1_e3t = 1._wp / e3t_0(ji,jj,jk)
            Cu_adv(ji,jj,jk) =   zdt *                                                         &
               &  ( ( MAX( ww(ji,jj,jk) , 0._wp ) - MIN( ww(ji,jj,jk+1) , 0._wp ) )            &
               &  + ( MAX( e2u(ji  ,jj) * e3u_0(ji  ,jj,jk)                                  &
               &                        * uu (ji  ,jj,jk,Kmm) + un_td(ji  ,jj,jk), 0._wp ) -   &
               &      MIN( e2u(ji-1,jj) * e3u_0(ji-1,jj,jk)                                  &
               &                        * uu (ji-1,jj,jk,Kmm) + un_td(ji-1,jj,jk), 0._wp ) )   &
               &                               * r1_e1e2t(ji,jj)                                                                     &
               &  + ( MAX( e1v(ji,jj  ) * e3v_0(ji,jj  ,jk)                                  &
               &                        * vv (ji,jj  ,jk,Kmm) + vn_td(ji,jj  ,jk), 0._wp ) -   &
               &      MIN( e1v(ji,jj-1) * e3v_0(ji,jj-1,jk)                                  &
               &                        * vv (ji,jj-1,jk,Kmm) + vn_td(ji,jj-1,jk), 0._wp ) )   &
               &                               * r1_e1e2t(ji,jj)                                                                     &
               &                             ) * z1_e3t
         END DO   ;   END DO   ;   END DO
      ELSE
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  nn_hls-1), ntej+(  nn_hls) ; DO ji = ntsi-( nn_hls-1), ntei+(  nn_hls)
            z1_e3t = 1._wp / e3t_0(ji,jj,jk)
            Cu_adv(ji,jj,jk) =   zdt *                                                      &
               &  ( ( MAX( ww(ji,jj,jk) , 0._wp ) - MIN( ww(ji,jj,jk+1) , 0._wp ) )         &
               &                             + ( MAX( e2u(ji  ,jj)*e3u_0(ji  ,jj,jk)*uu(ji  ,jj,jk,Kmm), 0._wp ) -   &
               &                                 MIN( e2u(ji-1,jj)*e3u_0(ji-1,jj,jk)*uu(ji-1,jj,jk,Kmm), 0._wp ) )   &
               &                               * r1_e1e2t(ji,jj)                                                 &
               &                             + ( MAX( e1v(ji,jj  )*e3v_0(ji,jj  ,jk)*vv(ji,jj  ,jk,Kmm), 0._wp ) -   &
               &                                 MIN( e1v(ji,jj-1)*e3v_0(ji,jj-1,jk)*vv(ji,jj-1,jk,Kmm), 0._wp ) )   &
               &                               * r1_e1e2t(ji,jj)                                                 &
               &                             ) * z1_e3t
         END DO   ;   END DO   ;   END DO
      ENDIF
      CALL iom_put("Courant",Cu_adv)
      !
      IF( MAXVAL( Cu_adv(:,:,:) ) > Cu_min ) THEN       ! Quick check if any breaches anywhere
         DO jk =  jpkm1,  2,  -1  ; DO jj = ntsj-(  nn_hls-1), ntej+(  nn_hls) ; DO ji = ntsi-( nn_hls-1), ntei+(  nn_hls)    ! or scan Courant criterion and partition ! w where necessary
            !
            zCu = MAX( Cu_adv(ji,jj,jk) , Cu_adv(ji,jj,jk-1) )
! alt:
!                  IF ( ww(ji,jj,jk) > 0._wp ) THEN 
!                     zCu =  Cu_adv(ji,jj,jk) 
!                  ELSE
!                     zCu =  Cu_adv(ji,jj,jk-1)
!                  ENDIF 
            !
            IF( zCu <= Cu_min ) THEN              !<-- Fully explicit
               zcff = 0._wp
            ELSEIF( zCu < Cu_cut ) THEN           !<-- Mixed explicit
               zcff = ( zCu - Cu_min )**2
               zcff = zcff / ( Fcu + zcff )
            ELSE                                  !<-- Mostly implicit
               zcff = ( zCu - Cu_max )/ zCu
            ENDIF
            zcff = MIN(1._wp, zcff)
            !
            wi(ji,jj,jk) =           zcff   * ww(ji,jj,jk)
            ww(ji,jj,jk) = ( 1._wp - zcff ) * ww(ji,jj,jk)
            !
            Cu_adv(ji,jj,jk) = zcff               ! Reuse array to output coefficient below and in stp_ctl
         END DO   ;   END DO   ;   END DO
         Cu_adv(:,:,1) = 0._wp 
      ELSE
         ! Fully explicit everywhere
         Cu_adv(:,:,:) = 0._wp                    ! Reuse array to output coefficient below and in stp_ctl
         wi    (:,:,:) = 0._wp
      ENDIF
      CALL iom_put("wimp",wi) 
      CALL iom_put("wi_cff",Cu_adv)
      CALL iom_put("wexp",ww)
      !
      IF( ln_timing )   CALL timing_stop('wAimp')
      !
   END SUBROUTINE wAimp
      
   !!======================================================================
END MODULE sshwzv












MODULE dynspg_ts

   !! Includes ROMS wd scheme with diagnostic outputs ; puu(:,:,:,Kmm) and puu(:,:,:,Krhs) updates are commented out ! 

   !!======================================================================
   !!                   ***  MODULE  dynspg_ts  ***
   !! Ocean dynamics:  surface pressure gradient trend, split-explicit scheme
   !!======================================================================
   !! History :   1.0  ! 2004-12  (L. Bessieres, G. Madec)  Original code
   !!              -   ! 2005-11  (V. Garnier, G. Madec)  optimization
   !!              -   ! 2006-08  (S. Masson)  distributed restart using iom
   !!             2.0  ! 2007-07  (D. Storkey) calls to BDY routines
   !!              -   ! 2008-01  (R. Benshila)  change averaging method
   !!             3.2  ! 2009-07  (R. Benshila, G. Madec) Complete revisit associated to vvl reactivation
   !!             3.3  ! 2010-09  (D. Storkey, E. O'Dea) update for BDY for Shelf configurations
   !!             3.3  ! 2011-03  (R. Benshila, R. Hordoir, P. Oddo) update calculation of ub_b
   !!             3.5  ! 2013-07  (J. Chanut) Switch to Forward-backward time stepping
   !!             3.6  ! 2013-11  (A. Coward) Update for z-tilde compatibility
   !!             3.7  ! 2015-11  (J. Chanut) free surface simplification
   !!              -   ! 2016-12  (G. Madec, E. Clementi) update for Stoke-Drift divergence
   !!             4.0  ! 2017-05  (G. Madec)  drag coef. defined at t-point (zdfdrg.F90)
   !!---------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_spg_ts     : compute surface pressure gradient trend using a time-splitting scheme 
   !!   dyn_spg_ts_init: initialisation of the time-splitting scheme
   !!   ts_wgt         : set time-splitting weights for temporal averaging (or not)
   !!   ts_rst         : read/write time-splitting fields in restart file
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE isf_oce         ! ice shelf variable (fwfisf)
   USE zdf_oce         ! vertical physics: variables
   USE zdfdrg          ! vertical physics: top/bottom drag coef.
   USE sbcapr          ! surface boundary condition: atmospheric pressure
   USE dynadv    , ONLY: ln_dynadv_vec
   USE dynvor          ! vortivity scheme indicators
   USE phycst          ! physical constants
   USE dynvor          ! vorticity term
   USE wet_dry         ! wetting/drying flux limter
   USE bdy_oce         ! open boundary
   USE bdyvol          ! open boundary volume conservation
   USE bdytides        ! open boundary condition data
   USE bdydyn2d        ! open boundary conditions on barotropic variables
   USE tide_mod        !
   USE sbcwave         ! surface wave
   !
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE prtctl          ! Print control
   USE iom             ! IOM library
   USE restart         ! only for lrst_oce

   USE iom   ! to remove

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_spg_ts        ! called by dyn_spg 
   PUBLIC dyn_spg_ts_init   !    -    - dyn_spg_init

   !! Time filtered arrays at baroclinic time step:
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   un_adv , vn_adv   !: Advection vel. at "now" barocl. step
   !
   INTEGER, SAVE :: icycle      ! Number of barotropic sub-steps for each internal step nn_e <= 2.5 nn_e
   REAL(dp),SAVE :: rDt_e       ! Barotropic time step
   !
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)    :: wgtbtp2! 1st & 2nd weights used in time filtering of barotropic fields
   REAL(dp), ALLOCATABLE, SAVE, DIMENSION(:)    :: wgtbtp1! 1st & 2nd weights used in time filtering of barotropic fields
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   zwz                ! ff_f/h at F points
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ftnw, ftne         ! triad of coriolis parameter
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ftsw, ftse         ! (only used with een vorticity scheme)

   REAL(wp) ::   r1_12 = 1._wp / 12._wp   ! local ratios
   REAL(wp) ::   r1_8  = 0.125_wp         !
   REAL(wp) ::   r1_4  = 0.25_wp          !
   REAL(wp) ::   r1_2  = 0.5_wp           !

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
   !! $Id: dynspg_ts.F90 15489 2021-11-10 09:18:39Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION dyn_spg_ts_alloc()
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_ts_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(3)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
      ALLOCATE( wgtbtp1(3*nn_e), wgtbtp2(3*nn_e), zwz(jpi,jpj), STAT=ierr(1) )
      IF( ln_dynvor_een .OR. ln_dynvor_eeT )   &
         &     ALLOCATE( ftnw(jpi,jpj) , ftne(jpi,jpj) , ftsw(jpi,jpj) , ftse(jpi,jpj), STAT=ierr(2)   )
         !
      ALLOCATE( un_adv(jpi,jpj), vn_adv(jpi,jpj)                    , STAT=ierr(3) )
      !
      dyn_spg_ts_alloc = MAXVAL( ierr(:) )
      !
      CALL mpp_sum( 'dynspg_ts', dyn_spg_ts_alloc )
      IF( dyn_spg_ts_alloc /= 0 )   CALL ctl_stop( 'STOP', 'dyn_spg_ts_alloc: failed to allocate arrays' )
      !
   END FUNCTION dyn_spg_ts_alloc


   SUBROUTINE dyn_spg_ts( kt, Kbb, Kmm, Krhs, puu, pvv, pssh, puu_b, pvv_b, Kaa, k_only_ADV )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose : - Compute the now trend due to the explicit time stepping
      !!              of the quasi-linear barotropic system, and add it to the
      !!              general momentum trend. 
      !!
      !! ** Method  : - split-explicit schem (time splitting) :
      !!      Barotropic variables are advanced from internal time steps
      !!      "n"   to "n+1" if ln_bt_fw=T
      !!      or from 
      !!      "n-1" to "n+1" if ln_bt_fw=F
      !!      thanks to a generalized forward-backward time stepping (see ref. below).
      !!
      !! ** Action :
      !!      -Update the filtered free surface at step "n+1"      : pssh(:,:,Kaa)
      !!      -Update filtered barotropic velocities at step "n+1" : puu_b(:,:,:,Kaa), vv_b(:,:,:,Kaa)
      !!      -Compute barotropic advective fluxes at step "n"     : un_adv, vn_adv
      !!      These are used to advect tracers and are compliant with discrete
      !!      continuity equation taken at the baroclinic time steps. This 
      !!      ensures tracers conservation.
      !!      - (puu(:,:,:,Krhs), pvv(:,:,:,Krhs)) momentum trend updated with barotropic component.
      !!
      !! References : Shchepetkin and McWilliams, Ocean Modelling, 2005. 
      !!---------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt                  ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Kmm, Krhs, Kaa ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv            ! ocean velocities and RHS of momentum equation
      REAL(wp), DIMENSION(jpi,jpj,jpt)    , INTENT(inout)  :: puu_b, pvv_b! SSH and barotropic velocities at main time levels
      REAL(dp), DIMENSION(jpi,jpj,jpt)    , INTENT(inout)  :: pssh! SSH and barotropic velocities at main time levels
      INTEGER , OPTIONAL                  , INTENT( in )  ::  k_only_ADV          ! only Advection in the RHS
      !
      INTEGER  ::   ji, jj, jk, jn        ! dummy loop indices
      LOGICAL  ::   ll_fw_start           ! =T : forward integration 
      LOGICAL  ::   ll_init               ! =T : special startup of 2d equations
      INTEGER  ::   noffset               ! local integers  : time offset for bdy update
      REAL(wp) ::   r1_Dt_b, z1_hu, z1_hv          ! local scalars
      REAL(wp)  :: za0, za2, za3!   -      -
      REAL(dp)  :: za1!   -      -
      REAL(wp) ::   zztmp, zldg               !   -      -
      REAL(wp)  :: zhu_bck, zhv_bck!   -      -
      REAL(dp)  :: zhdiv!   -      -
      REAL(wp) ::   zun_save, zvn_save              !   -      -
      REAL(wp), DIMENSION(jpi,jpj) :: zu_trd, zu_frc, zu_spg
      REAL(wp), DIMENSION(jpi,jpj) :: zv_trd, zv_frc, zv_spg
      REAL(wp), DIMENSION(jpi,jpj) :: zsshu_a, zhup2_e, zhtp2_e
      REAL(wp), DIMENSION(jpi,jpj) :: zsshv_a, zhvp2_e, zsshp2_e
      REAL(wp), DIMENSION(jpi,jpj) :: zCdU_u, zCdU_v   ! top/bottom stress at u- & v-points
      REAL(wp), DIMENSION(jpi,jpj)  :: zhV! fluxes
      REAL(dp), DIMENSION(jpi,jpj)  :: zhU! fluxes
!!st#if defined key_qco 
!!st      REAL(wp), DIMENSION(jpi, jpj, jpk) :: ze3u, ze3v
!!st#endif
      !
      REAL(wp) ::   zwdramp                     ! local scalar - only used if ln_wd_dl = .True. 

      INTEGER  :: iwdg, jwdg, kwdg   ! short-hand values for the indices of the output point

      REAL(wp) ::   zepsilon, zgamma            !   -      -
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: ztwdmask, zuwdmask, zvwdmask ! ROMS wetting and drying masks at t,u,v points
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zuwdav2, zvwdav2    ! averages over the sub-steps of zuwdmask and zvwdmask
      REAL(wp) ::   zt0substep !   Time of day at the beginning of the time substep
      !!----------------------------------------------------------------------
      !
      !                                         !* Allocate temporary arrays
      IF( ln_wd_dl ) ALLOCATE( ztwdmask(jpi,jpj), zuwdmask(jpi,jpj), zvwdmask(jpi,jpj), zuwdav2(jpi,jpj), zvwdav2(jpi,jpj))
      !
      zwdramp = r_rn_wdmin1               ! simplest ramp 
!     zwdramp = 1._wp / (rn_wdmin2 - rn_wdmin1) ! more general ramp
      !                                         ! inverse of baroclinic time step 
      r1_Dt_b = 1._wp / rDt 
      !
      ll_init     = ln_bt_av                    ! if no time averaging, then no specific restart 
      ll_fw_start = .FALSE.
      !                                         ! time offset in steps for bdy data update
      IF( .NOT.ln_bt_fw ) THEN   ;   noffset = - nn_e
      ELSE                       ;   noffset =   0 
      ENDIF
      !
      IF( kt == nit000 ) THEN                   !* initialisation
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_ts : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~   free surface with time splitting'
         IF(lwp) WRITE(numout,*)
         !
         IF( l_1st_euler )   ll_init=.TRUE.
         !
         IF( ln_bt_fw .OR. l_1st_euler ) THEN
            ll_fw_start =.TRUE.
            noffset     = 0
         ELSE
            ll_fw_start =.FALSE.
         ENDIF
         !                    ! Set averaging weights and cycle length:
         CALL ts_wgt( ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2 )
         !
      ELSEIF( kt == nit000 + 1 ) THEN           !* initialisation 2nd time-step
         !
         IF( .NOT.ln_bt_fw ) THEN
            ! If we did an Euler timestep on the first timestep we need to reset ll_fw_start
            ! and the averaging weights. We don't have an easy way of telling whether we did
            ! an Euler timestep on the first timestep (because l_1st_euler is reset to .false.
            ! at the end of the first timestep) so just do this in all cases. 
            ll_fw_start = .FALSE.
            CALL ts_wgt( ln_bt_av, ll_fw_start, icycle, wgtbtp1, wgtbtp2 )
         ENDIF
         !
      ENDIF
      !
      ! -----------------------------------------------------------------------------
      !  Phase 1 : Coupling between general trend and barotropic estimates (1st step)
      ! -----------------------------------------------------------------------------
      !      
      !
      !                                   !=  zu_frc =  1/H e3*d/dt(Ua)  =!  (Vertical mean of Ua, the 3D trends)
      !                                   !  ---------------------------  !
      zu_frc(:,:) = SUM( e3u_0(:,:,:) * puu(:,:,:,Krhs) * umask(:,:,:), DIM=3 ) * r1_hu_0(:,:)
      zv_frc(:,:) = SUM( e3v_0(:,:,:) * pvv(:,:,:,Krhs) * vmask(:,:,:), DIM=3 ) * r1_hv_0(:,:)
      !
      !
      !                                   !=  U(Krhs) => baroclinic trend  =!   (remove its vertical mean)
      DO jk = 1, jpkm1                    !  -----------------------------  !
         puu(:,:,jk,Krhs) = ( puu(:,:,jk,Krhs) - zu_frc(:,:) ) * umask(:,:,jk)
         pvv(:,:,jk,Krhs) = ( pvv(:,:,jk,Krhs) - zv_frc(:,:) ) * vmask(:,:,jk)
      END DO
      
!!gm  Question here when removing the Vertically integrated trends, we remove the vertically integrated NL trends on momentum....
!!gm  Is it correct to do so ?   I think so...
      
      !                                   !=  remove 2D Coriolis trend  =!
      !                                   !  --------------------------  !
      !
      IF( kt == nit000 .OR. .NOT. ln_linssh )   CALL dyn_cor_2D_init( Kmm )   ! Set zwz, the barotropic Coriolis force coefficient
      !                      ! recompute zwz = f/depth  at every time step for (.NOT.ln_linssh) as the water colomn height changes
      !
      IF( .NOT. PRESENT(k_only_ADV) ) THEN   !* remove the 2D Coriolis trend  
         zhU(:,:) = puu_b(:,:,Kmm) * hu_0(:,:) * e2u(:,:)        ! now fluxes 
         zhV(:,:) = pvv_b(:,:,Kmm) * hv_0(:,:) * e1v(:,:)        ! NB: FULL domain : put a value in last row and column
         !
         CALL dyn_cor_2d( ht_0(:,:), hu_0(:,:), hv_0(:,:), puu_b(:,:,Kmm), pvv_b(:,:,Kmm), zhU, zhV,  &   ! <<== in
            &                                                                          zu_trd, zv_trd   )   ! ==>> out
         !
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)                          ! Remove coriolis term (and possibly spg) from barotropic trend
            zu_frc(ji,jj) = zu_frc(ji,jj) - zu_trd(ji,jj) * ssumask(ji,jj)
            zv_frc(ji,jj) = zv_frc(ji,jj) - zv_trd(ji,jj) * ssvmask(ji,jj)
         END DO   ;   END DO
      ENDIF
      !
      !                                   !=  Add bottom stress contribution from baroclinic velocities  =!
      !                                   !  -----------------------------------------------------------  !
      IF( PRESENT(k_only_ADV) ) THEN         !* only Advection in the RHS : provide the barotropic bottom drag coefficients
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zCdU_u(ji,jj) = r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) )
            zCdU_v(ji,jj) = r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) )
         END DO   ;   END DO
      ELSE				     !* remove baroclinic drag AND provide the barotropic drag coefficients
         CALL dyn_drg_init( Kbb, Kmm, puu, pvv, puu_b, pvv_b, zu_frc, zv_frc, zCdU_u, zCdU_v )
      ENDIF
      !
      !                                   !=  Add atmospheric pressure forcing  =!
      !                                   !  ----------------------------------  !
      IF( ln_apr_dyn ) THEN
         IF( ln_bt_fw ) THEN                          ! FORWARD integration: use kt+1/2 pressure (NOW+1/2)
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               zu_frc(ji,jj) = zu_frc(ji,jj) + grav * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj) ) * r1_e1u(ji,jj)
               zv_frc(ji,jj) = zv_frc(ji,jj) + grav * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj) ) * r1_e2v(ji,jj)
            END DO   ;   END DO
         ELSE                                         ! CENTRED integration: use kt-1/2 + kt+1/2 pressure (NOW)
            zztmp = grav * r1_2
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               zu_frc(ji,jj) = zu_frc(ji,jj) + zztmp * (  ssh_ib (ji+1,jj  ) - ssh_ib (ji,jj)  &
                    &                                   + ssh_ibb(ji+1,jj  ) - ssh_ibb(ji,jj)  ) * r1_e1u(ji,jj)
               zv_frc(ji,jj) = zv_frc(ji,jj) + zztmp * (  ssh_ib (ji  ,jj+1) - ssh_ib (ji,jj)  &
                    &                                   + ssh_ibb(ji  ,jj+1) - ssh_ibb(ji,jj)  ) * r1_e2v(ji,jj)
            END DO   ;   END DO
         ENDIF
      ENDIF
      !
      !                                   !=  Add wind forcing  =!
      !                                   !  ------------------  !
      IF( ln_bt_fw ) THEN
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zu_frc(ji,jj) =  zu_frc(ji,jj) + r1_rho0 * utau(ji,jj) * r1_hu_0(ji,jj)
            zv_frc(ji,jj) =  zv_frc(ji,jj) + r1_rho0 * vtau(ji,jj) * r1_hv_0(ji,jj)
         END DO   ;   END DO
      ELSE
         zztmp = r1_rho0 * r1_2
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zu_frc(ji,jj) =  zu_frc(ji,jj) + zztmp * ( utau_b(ji,jj) + utau(ji,jj) ) * r1_hu_0(ji,jj)
            zv_frc(ji,jj) =  zv_frc(ji,jj) + zztmp * ( vtau_b(ji,jj) + vtau(ji,jj) ) * r1_hv_0(ji,jj)
         END DO   ;   END DO
      ENDIF  
      !
      !              !----------------!
      !              !==  sssh_frc  ==!   Right-Hand-Side of the barotropic ssh equation   (over the FULL domain)
      !              !----------------!
      !                                   !=  Net water flux forcing applied to a water column  =!
      !                                   ! ---------------------------------------------------  !
      IF (ln_bt_fw) THEN                          ! FORWARD integration: use kt+1/2 fluxes (NOW+1/2)
         ssh_frc(:,:) = r1_rho0 * ( emp(:,:) - rnf(:,:) - fwfisf_cav(:,:) - fwfisf_par(:,:) )
      ELSE                                        ! CENTRED integration: use kt-1/2 + kt+1/2 fluxes (NOW)
         zztmp = r1_rho0 * r1_2
         ssh_frc(:,:) = zztmp * (   emp(:,:)        + emp_b(:,:)          &
            &                     - rnf(:,:)        - rnf_b(:,:)          &
            &                     - fwfisf_cav(:,:) - fwfisf_cav_b(:,:)   &
            &                     - fwfisf_par(:,:) - fwfisf_par_b(:,:)   )
      ENDIF
      !                                   !=  Add Stokes drift divergence  =!   (if exist)
      IF( ln_sdw ) THEN                   !  -----------------------------  !
         ssh_frc(:,:) = ssh_frc(:,:) + div_sd(:,:)
      ENDIF
      !
      !                                         ! ice sheet coupling
      IF ( ln_isf .AND. ln_isfcpl ) THEN
         !
         ! ice sheet coupling
         IF( ln_rstart .AND. kt == nit000 ) THEN
            ssh_frc(:,:) = ssh_frc(:,:) + risfcpl_ssh(:,:)
         END IF
         !
         ! conservation option
         IF( ln_isfcpl_cons ) THEN
            ssh_frc(:,:) = ssh_frc(:,:) + risfcpl_cons_ssh(:,:)
         END IF
         !
      END IF
      !
      !                                   != Fill boundary data arrays for AGRIF
      !                                   ! ------------------------------------
      !
      ! -----------------------------------------------------------------------
      !  Phase 2 : Integration of the barotropic equations 
      ! -----------------------------------------------------------------------
      !
      !                                             ! ==================== !
      !                                             !    Initialisations   !
      !                                             ! ==================== !  
      ! Initialize barotropic variables:      
      IF( ll_init )THEN
         sshbb_e(:,:) = 0._wp
         ubb_e  (:,:) = 0._wp
         vbb_e  (:,:) = 0._wp
         sshb_e (:,:) = 0._wp
         ub_e   (:,:) = 0._wp
         vb_e   (:,:) = 0._wp
      ENDIF
      !
      IF( ln_linssh ) THEN    ! mid-step ocean depth is fixed (hup2_e=hu_n=hu_0)
         zhup2_e(:,:) = hu_0(:,:)
         zhvp2_e(:,:) = hv_0(:,:)
         zhtp2_e(:,:) = ht_0(:,:)
      ENDIF
      !
      IF( ln_bt_fw ) THEN                 ! FORWARD integration: start from NOW fields                    
         sshn_e(:,:) =    pssh (:,:,Kmm)            
         un_e  (:,:) =    puu_b(:,:,Kmm)            
         vn_e  (:,:) =    pvv_b(:,:,Kmm)
         !
         hu_e  (:,:) =    hu_0(:,:)       
         hv_e  (:,:) =    hv_0(:,:) 
         hur_e (:,:) = r1_hu_0(:,:)    
         hvr_e (:,:) = r1_hv_0(:,:)
      ELSE                                ! CENTRED integration: start from BEFORE fields
         sshn_e(:,:) =    pssh (:,:,Kbb)
         un_e  (:,:) =    puu_b(:,:,Kbb)         
         vn_e  (:,:) =    pvv_b(:,:,Kbb)
         !
         hu_e  (:,:) =    hu_0(:,:)       
         hv_e  (:,:) =    hv_0(:,:) 
         hur_e (:,:) = r1_hu_0(:,:)    
         hvr_e (:,:) = r1_hv_0(:,:)
      ENDIF
      !
      ! Initialize sums:
      puu_b (:,:,Kaa) = 0._wp       ! After barotropic velocities (or transport if flux form)          
      pvv_b (:,:,Kaa) = 0._wp
      pssh  (:,:,Kaa) = 0._wp       ! Sum for after averaged sea level
      un_adv(:,:)     = 0._wp       ! Sum for now transport issued from ts loop
      vn_adv(:,:)     = 0._wp
      !
      IF( ln_wd_dl ) THEN
         zuwdmask(:,:) = 0._wp  ! set to zero for definiteness (not sure this is necessary) 
         zvwdmask(:,:) = 0._wp  ! 
         zuwdav2 (:,:) = 0._wp 
         zvwdav2 (:,:) = 0._wp   
      END IF 

      !                                             ! ==================== !
      DO jn = 1, icycle                             !  sub-time-step loop  !
         !                                          ! ==================== !
         !
         l_full_nf_update = jn == icycle   ! false: disable full North fold update (performances) for jn = 1 to icycle-1
         !
         !                    !==  Update the forcing ==! (BDY and tides)
         !
         IF( ln_bdy      .AND. ln_tide )   CALL bdy_dta_tides( kt, kit=jn, pt_offset= REAL(noffset+1,wp) )
         ! Update tide potential at the beginning of current time substep
         IF( ln_tide_pot .AND. ln_tide ) THEN
            zt0substep = REAL(nsec_day, wp) - 0.5_wp*rn_Dt + (jn + noffset - 1) * rn_Dt / REAL(nn_e, wp)
            CALL upd_tide(zt0substep, Kmm)
         END IF
         !
         !                    !==  extrapolation at mid-step  ==!   (jn+1/2)
         !
         !                       !* Set extrapolation coefficients for predictor step:
         IF ((jn<3).AND.ll_init) THEN      ! Forward           
           za1 = 1._wp                                          
           za2 = 0._wp                        
           za3 = 0._wp                        
         ELSE                              ! AB3-AM4 Coefficients: bet=0.281105 
           za1 =  1.781105_wp              ! za1 =   3/2 +   bet
           za2 = -1.06221_wp               ! za2 = -(1/2 + 2*bet)
           za3 =  0.281105_wp              ! za3 = bet
         ENDIF
         !
         !                       !* Extrapolate barotropic velocities at mid-step (jn+1/2)
         !--        m+1/2               m                m-1           m-2       --!
         !--       u      = (3/2+beta) u   -(1/2+2beta) u      + beta u          --!
         !-------------------------------------------------------------------------!
         ua_e(:,:) = za1 * un_e(:,:) + za2 * ub_e(:,:) + za3 * ubb_e(:,:)
         va_e(:,:) = za1 * vn_e(:,:) + za2 * vb_e(:,:) + za3 * vbb_e(:,:)

         IF( .NOT.ln_linssh ) THEN                        !* Update ocean depth (variable volume case only)
            !                                             !  ------------------
            ! Extrapolate Sea Level at step jit+0.5:
            !--         m+1/2                 m                  m-1             m-2       --!
            !--      ssh      = (3/2+beta) ssh   -(1/2+2beta) ssh      + beta ssh          --!
            !--------------------------------------------------------------------------------!
            zsshp2_e(:,:) = za1 * sshn_e(:,:)  + za2 * sshb_e(:,:) + za3 * sshbb_e(:,:)
            
            ! set wetting & drying mask at tracer points for this barotropic mid-step
            IF( ln_wd_dl )   CALL wad_tmsk( zsshp2_e, ztwdmask )
            !
            !                          ! ocean t-depth at mid-step
            zhtp2_e(:,:) = ht_0(:,:) + zsshp2_e(:,:)
            !
            !                          ! ocean u- and v-depth at mid-step   (separate DO-loops remove the need of a lbc_lnk)
            !                                ! no 'key_qcoTest_FluxForm' : surface weighted ssh average
            DO jj = ntsj-( 1), ntej+( 1 ) ; DO ji = ntsi-( 1), ntei+( 0)   ! not jpi-column
               zhup2_e(ji,jj) = hu_0(ji,jj) + r1_2 * r1_e1e2u(ji,jj)                        &
                    &                              * (  e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)  &
                    &                                 + e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj)  ) * ssumask(ji,jj)
            END DO   ;   END DO
            DO jj = ntsj-( 1), ntej+( 0 ) ; DO ji = ntsi-( 1), ntei+( 1)   ! not jpj-row
               zhvp2_e(ji,jj) = hv_0(ji,jj) + r1_2 * r1_e1e2v(ji,jj)                        &
                    &                              * (  e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )  &
                    &                                 + e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1)  ) * ssvmask(ji,jj)
            END DO   ;   END DO
            !
         ENDIF
         !
         !                    !==  after SSH  ==!   (jn+1)
         !
         !                             ! update (ua_e,va_e) to enforce volume conservation at open boundaries
         !                             ! values of zhup2_e and zhvp2_e on the halo are not needed in bdy_vol2d
         IF( ln_bdy .AND. ln_vol ) CALL bdy_vol2d( kt, jn, ua_e, va_e, zhup2_e, zhvp2_e )
         !      
         !                             ! resulting flux at mid-step (not over the full domain)
         DO jj = ntsj-( 1), ntej+( 1 ) ; DO ji = ntsi-( 1), ntei+( 0)   ! not jpi-column
            zhU(ji,jj) = e2u(ji,jj) * ua_e(ji,jj) * zhup2_e(ji,jj)
         END DO   ;   END DO
         DO jj = ntsj-( 1), ntej+( 0 ) ; DO ji = ntsi-( 1), ntei+( 1)   ! not jpj-row
            zhV(ji,jj) = e1v(ji,jj) * va_e(ji,jj) * zhvp2_e(ji,jj)
         END DO   ;   END DO
         !

         IF( ln_wd_dl ) THEN           ! un_e and vn_e are set to zero at faces where 
            !                          ! the direction of the flow is from dry cells
            CALL wad_Umsk( ztwdmask, zhU, zhV, un_e, vn_e, zuwdmask, zvwdmask )   ! not jpi colomn for U, not jpj row for V
            !
         ENDIF    
         !
         !
         !     Compute Sea Level at step jit+1
         !--           m+1        m                               m+1/2          --!
         !--        ssh    =  ssh   - delta_t' * [ frc + div( flux      ) ]      --!
         !-------------------------------------------------------------------------!
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zhdiv = (   zhU(ji,jj) - zhU(ji-1,jj) + zhV(ji,jj) - zhV(ji,jj-1)   ) * r1_e1e2t(ji,jj)
            ssha_e(ji,jj) = (  sshn_e(ji,jj) - rDt_e * ( ssh_frc(ji,jj) + zhdiv )  ) * ssmask(ji,jj)
         END DO   ;   END DO
         !
         CALL lbc_lnk( 'dynspg_ts', ssha_e, 'T', 1._dp,  zhU, 'U', -1._dp)
         CALL lbc_lnk( 'dynspg_ts',  zhV, 'V', -1._wp )
         !
         ! Duplicate sea level across open boundaries (this is only cosmetic if linssh=T)
         IF( ln_bdy )   CALL bdy_ssh( ssha_e )
         !
         !                             ! Sum over sub-time-steps to compute advective velocities
         za2 = wgtbtp2(jn)             ! zhU, zhV hold fluxes extrapolated at jn+0.5
         un_adv(:,:) = un_adv(:,:) + za2 * zhU(:,:) * r1_e2u(:,:)
         vn_adv(:,:) = vn_adv(:,:) + za2 * zhV(:,:) * r1_e1v(:,:)
         ! sum over sub-time-steps to decide which baroclinic velocities to set to zero (zuwdav2 is only used when ln_wd_dl_bc=True) 
         IF ( ln_wd_dl_bc ) THEN
            DO jj = ntsj-( 1), ntej+( 1 ) ; DO ji = ntsi-( 1), ntei+( 0)   ! not jpi-column
               zuwdav2(ji,jj) = zuwdav2(ji,jj) + za2 * zuwdmask(ji,jj)
            END DO   ;   END DO
            DO jj = ntsj-( 1), ntej+( 0 ) ; DO ji = ntsi-( 1), ntei+( 1)   ! not jpj-row
               zvwdav2(ji,jj) = zvwdav2(ji,jj) + za2 * zvwdmask(ji,jj)
            END DO   ;   END DO
         END IF
         !
         !  
         ! Sea Surface Height at u-,v-points (vvl case only)
         IF( .NOT.ln_linssh ) THEN
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               zsshu_a(ji,jj) = r1_2 * r1_e1e2u(ji,jj) * ( e1e2t(ji  ,jj  ) * ssha_e(ji  ,jj  )   &
                  &                                      + e1e2t(ji+1,jj  ) * ssha_e(ji+1,jj  ) ) * ssumask(ji,jj)
               zsshv_a(ji,jj) = r1_2 * r1_e1e2v(ji,jj) * ( e1e2t(ji  ,jj  ) * ssha_e(ji  ,jj  )   &
                  &                                      + e1e2t(ji  ,jj+1) * ssha_e(ji  ,jj+1) ) * ssvmask(ji,jj)
            END DO   ;   END DO
         ENDIF
         !         
         ! Half-step back interpolation of SSH for surface pressure computation at step jit+1/2
         !--            m+1/2           m+1              m               m-1              m-2     --!
         !--        ssh'    =  za0 * ssh     +  za1 * ssh   +  za2 * ssh      +  za3 * ssh        --!
         !------------------------------------------------------------------------------------------!
         CALL ts_bck_interp( jn, ll_init, za0, za1, za2, za3 )   ! coeficients of the interpolation
         zsshp2_e(:,:) = za0 *  ssha_e(:,:) + za1 *  sshn_e (:,:)   &
            &          + za2 *  sshb_e(:,:) + za3 *  sshbb_e(:,:)
         !
         !                             ! Surface pressure gradient
         zldg = ( 1._wp - rn_scal_load ) * grav    ! local factor
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zu_spg(ji,jj) = - zldg * ( zsshp2_e(ji+1,jj) - zsshp2_e(ji,jj) ) * r1_e1u(ji,jj)
            zv_spg(ji,jj) = - zldg * ( zsshp2_e(ji,jj+1) - zsshp2_e(ji,jj) ) * r1_e2v(ji,jj)
         END DO   ;   END DO
         !
         ! Add Coriolis trend:
         ! zwz array below or triads normally depend on sea level with ln_linssh=F and should be updated
         ! at each time step. We however keep them constant here for optimization.
         ! Recall that zhU and zhV hold fluxes at jn+0.5 (extrapolated not backward interpolated)
         CALL dyn_cor_2d( zhtp2_e, zhup2_e, zhvp2_e, ua_e, va_e, zhU, zhV,    zu_trd, zv_trd   )
         !
         ! Add tidal astronomical forcing if defined
         IF ( ln_tide .AND. ln_tide_pot ) THEN
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               zu_trd(ji,jj) = zu_trd(ji,jj) + grav * ( pot_astro(ji+1,jj) - pot_astro(ji,jj) ) * r1_e1u(ji,jj)
               zv_trd(ji,jj) = zv_trd(ji,jj) + grav * ( pot_astro(ji,jj+1) - pot_astro(ji,jj) ) * r1_e2v(ji,jj)
            END DO   ;   END DO
         ENDIF
         !
         ! Add bottom stresses:
!jth do implicitly instead
         IF ( .NOT. ll_wd ) THEN ! Revert to explicit for bit comparison tests in non wad runs
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               zu_trd(ji,jj) = zu_trd(ji,jj) + zCdU_u(ji,jj) * un_e(ji,jj) * hur_e(ji,jj)
               zv_trd(ji,jj) = zv_trd(ji,jj) + zCdU_v(ji,jj) * vn_e(ji,jj) * hvr_e(ji,jj)
            END DO   ;   END DO
         ENDIF
         !
         ! Set next velocities:
         !     Compute barotropic speeds at step jit+1    (h : total height of the water colomn)
         !--                              VECTOR FORM
         !--   m+1                 m               /                                                       m+1/2           \    --!
         !--  u     =             u   + delta_t' * \         (1-r)*g * grad_x( ssh') -         f * k vect u      +     frc /    --!
         !--                                                                                                                    --!
         !--                             FLUX FORM                                                                              --!
         !--  m+1   __1__  /  m    m               /  m+1/2                             m+1/2              m+1/2    n      \ \  --!
         !-- u    =   m+1 |  h  * u   + delta_t' * \ h     * (1-r)*g * grad_x( ssh') - h     * f * k vect u      + h * frc /  | --!
         !--         h     \                                                                                                 /  --!
         !------------------------------------------------------------------------------------------------------------------------!
         IF( ln_dynadv_vec .OR. ln_linssh ) THEN      !* Vector form
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               ua_e(ji,jj) = (                                 un_e(ji,jj)   & 
                         &     + rDt_e * (                   zu_spg(ji,jj)   &
                         &                                 + zu_trd(ji,jj)   &
                         &                                 + zu_frc(ji,jj) ) & 
                         &   ) * ssumask(ji,jj)

               va_e(ji,jj) = (                                 vn_e(ji,jj)   &
                         &     + rDt_e * (                   zv_spg(ji,jj)   &
                         &                                 + zv_trd(ji,jj)   &
                         &                                 + zv_frc(ji,jj) ) &
                         &   ) * ssvmask(ji,jj)
            END DO   ;   END DO
            !
         ELSE                           !* Flux form
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               !                    ! hu_e, hv_e hold depth at jn,  zhup2_e, zhvp2_e hold extrapolated depth at jn+1/2
               !                    ! backward interpolated depth used in spg terms at jn+1/2
               zhu_bck = hu_0(ji,jj) + r1_2*r1_e1e2u(ji,jj) * (  e1e2t(ji  ,jj) * zsshp2_e(ji  ,jj)    &
                    &                                          + e1e2t(ji+1,jj) * zsshp2_e(ji+1,jj)  ) * ssumask(ji,jj)
               zhv_bck = hv_0(ji,jj) + r1_2*r1_e1e2v(ji,jj) * (  e1e2t(ji,jj  ) * zsshp2_e(ji,jj  )    &
                    &                                          + e1e2t(ji,jj+1) * zsshp2_e(ji,jj+1)  ) * ssvmask(ji,jj)
               !                    ! inverse depth at jn+1
               z1_hu = ssumask(ji,jj) / ( hu_0(ji,jj) + zsshu_a(ji,jj) + 1._wp - ssumask(ji,jj) )
               z1_hv = ssvmask(ji,jj) / ( hv_0(ji,jj) + zsshv_a(ji,jj) + 1._wp - ssvmask(ji,jj) )
               !
               ua_e(ji,jj) = (               hu_e  (ji,jj) *   un_e (ji,jj)      & 
                    &            + rDt_e * (  zhu_bck        * zu_spg (ji,jj)  &   !
                    &                       + zhup2_e(ji,jj) * zu_trd (ji,jj)  &   !
                    &                       +  hu_0(ji,jj) * zu_frc (ji,jj)  )   ) * z1_hu
               !
               va_e(ji,jj) = (               hv_e  (ji,jj) *   vn_e (ji,jj)      &
                    &            + rDt_e * (  zhv_bck        * zv_spg (ji,jj)  &   !
                    &                       + zhvp2_e(ji,jj) * zv_trd (ji,jj)  &   !
                    &                       +  hv_0(ji,jj) * zv_frc (ji,jj)  )   ) * z1_hv
            END DO   ;   END DO
         ENDIF
!jth implicit bottom friction:
         IF ( ll_wd ) THEN ! revert to explicit for bit comparison tests in non wad runs
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               ua_e(ji,jj) =  ua_e(ji,jj) / ( 1._wp - rDt_e * zCdU_u(ji,jj) * hur_e(ji,jj) )
               va_e(ji,jj) =  va_e(ji,jj) / ( 1._wp - rDt_e * zCdU_v(ji,jj) * hvr_e(ji,jj) )
            END DO   ;   END DO
         ENDIF
       
         IF( .NOT.ln_linssh ) THEN !* Update ocean depth (variable volume case only)
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               hu_e (ji,jj) =    hu_0(ji,jj) + zsshu_a(ji,jj)
               hur_e(ji,jj) = ssumask(ji,jj) / (  hu_e(ji,jj) + 1._wp - ssumask(ji,jj)  )
               hv_e (ji,jj) =    hv_0(ji,jj) + zsshv_a(ji,jj)
               hvr_e(ji,jj) = ssvmask(ji,jj) / (  hv_e(ji,jj) + 1._wp - ssvmask(ji,jj)  )
            END DO   ;   END DO
         ENDIF
         !
         IF( .NOT.ln_linssh ) THEN   !* Update ocean depth (variable volume case only)
            CALL lbc_lnk( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  &
                 &                   , hu_e , 'U',  1._wp, hv_e , 'V',  1._wp  &
                 &                   , hur_e, 'U',  1._wp, hvr_e, 'V',  1._wp  )
         ELSE
            CALL lbc_lnk( 'dynspg_ts', ua_e , 'U', -1._wp, va_e , 'V', -1._wp  )
         ENDIF
         !                                                 ! open boundaries
         IF( ln_bdy )   CALL bdy_dyn2d( jn, ua_e, va_e, un_e, vn_e, hur_e, hvr_e, ssha_e )
         !                                             !* Swap
         !                                             !  ----
         ubb_e  (:,:) = ub_e  (:,:)
         ub_e   (:,:) = un_e  (:,:)
         un_e   (:,:) = ua_e  (:,:)
         !
         vbb_e  (:,:) = vb_e  (:,:)
         vb_e   (:,:) = vn_e  (:,:)
         vn_e   (:,:) = va_e  (:,:)
         !
         sshbb_e(:,:) = sshb_e(:,:)
         sshb_e (:,:) = sshn_e(:,:)
         sshn_e (:,:) = ssha_e(:,:)

         !                                             !* Sum over whole bt loop
         !                                             !  ----------------------
         za1 = wgtbtp1(jn)                                    
         IF( ln_dynadv_vec .OR. ln_linssh ) THEN    ! Sum velocities
            puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) 
            pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) 
         ELSE                                       ! Sum transports
            IF ( .NOT.ln_wd_dl ) THEN  
               puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) * hu_e (:,:)
               pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) * hv_e (:,:)
            ELSE 
               puu_b  (:,:,Kaa) = puu_b  (:,:,Kaa) + za1 * ua_e  (:,:) * hu_e (:,:) * zuwdmask(:,:)
               pvv_b  (:,:,Kaa) = pvv_b  (:,:,Kaa) + za1 * va_e  (:,:) * hv_e (:,:) * zvwdmask(:,:)
            END IF 
         ENDIF
         !                                          ! Sum sea level
         pssh(:,:,Kaa) = pssh(:,:,Kaa) + za1 * ssha_e(:,:)

         !                                                 ! ==================== !
      END DO                                               !        end loop      !
      !                                                    ! ==================== !
      ! -----------------------------------------------------------------------------
      ! Phase 3. update the general trend with the barotropic trend
      ! -----------------------------------------------------------------------------
      !
      ! Set advection velocity correction:
      IF (ln_bt_fw) THEN
         IF( .NOT.( kt == nit000 .AND. l_1st_euler ) ) THEN
            DO jj = ntsj-( nn_hls), ntej+( nn_hls ) ; DO ji = ntsi-( nn_hls), ntei+( nn_hls)
               zun_save = un_adv(ji,jj)
               zvn_save = vn_adv(ji,jj)
               !                          ! apply the previously computed correction 
               un_adv(ji,jj) = r1_2 * ( ub2_b(ji,jj) + zun_save - rn_atfp * un_bf(ji,jj) )
               vn_adv(ji,jj) = r1_2 * ( vb2_b(ji,jj) + zvn_save - rn_atfp * vn_bf(ji,jj) )
               !                          ! Update corrective fluxes for next time step
               un_bf(ji,jj)  = rn_atfp * un_bf(ji,jj) + ( zun_save - ub2_b(ji,jj) )
               vn_bf(ji,jj)  = rn_atfp * vn_bf(ji,jj) + ( zvn_save - vb2_b(ji,jj) )
               !                          ! Save integrated transport for next computation
               ub2_b(ji,jj) = zun_save
               vb2_b(ji,jj) = zvn_save
            END DO   ;   END DO
         ELSE
            un_bf(:,:) = 0._wp            ! corrective fluxes for next time step set to zero
            vn_bf(:,:) = 0._wp
            ub2_b(:,:) = un_adv(:,:)      ! Save integrated transport for next computation
            vb2_b(:,:) = vn_adv(:,:)
         END IF
      ENDIF


      !
      ! Update barotropic trend:
      IF( ln_dynadv_vec .OR. ln_linssh ) THEN
         DO jk=1,jpkm1
            puu(:,:,jk,Krhs) = puu(:,:,jk,Krhs) + ( puu_b(:,:,Kaa) - puu_b(:,:,Kbb) ) * r1_Dt_b
            pvv(:,:,jk,Krhs) = pvv(:,:,jk,Krhs) + ( pvv_b(:,:,Kaa) - pvv_b(:,:,Kbb) ) * r1_Dt_b
         END DO
      ELSE
         ! At this stage, pssh(:,:,:,Krhs) has been corrected: compute new depths at velocity points
         DO jj = ntsj-( 1), ntej+( 0 ) ; DO ji = ntsi-( 1), ntei+( 0)
            zsshu_a(ji,jj) = r1_2 * r1_e1e2u(ji,jj) * ( e1e2t(ji  ,jj) * pssh(ji  ,jj,Kaa)   &
               &                                      + e1e2t(ji+1,jj) * pssh(ji+1,jj,Kaa) ) * ssumask(ji,jj)
            zsshv_a(ji,jj) = r1_2 * r1_e1e2v(ji,jj) * ( e1e2t(ji,jj  ) * pssh(ji,jj  ,Kaa)   &
               &                                      + e1e2t(ji,jj+1) * pssh(ji,jj+1,Kaa) ) * ssvmask(ji,jj)
         END DO   ;   END DO
         CALL lbc_lnk( 'dynspg_ts', zsshu_a, 'U', 1._wp, zsshv_a, 'V', 1._wp ) ! Boundary conditions
         !
         DO jk=1,jpkm1
            puu(:,:,jk,Krhs) = puu(:,:,jk,Krhs) + r1_hu_0(:,:)   &
               &             * ( puu_b(:,:,Kaa) - puu_b(:,:,Kbb) * hu_0(:,:) ) * r1_Dt_b
            pvv(:,:,jk,Krhs) = pvv(:,:,jk,Krhs) + r1_hv_0(:,:)   &
               &             * ( pvv_b(:,:,Kaa) - pvv_b(:,:,Kbb) * hv_0(:,:) ) * r1_Dt_b
         END DO
         ! Save barotropic velocities not transport:
         puu_b(:,:,Kaa) =  puu_b(:,:,Kaa) / ( hu_0(:,:) + zsshu_a(:,:) + 1._wp - ssumask(:,:) )
         pvv_b(:,:,Kaa) =  pvv_b(:,:,Kaa) / ( hv_0(:,:) + zsshv_a(:,:) + 1._wp - ssvmask(:,:) )
      ENDIF


      ! Correct velocities so that the barotropic velocity equals (un_adv, vn_adv) (in all cases)  
      DO jk = 1, jpkm1
         puu(:,:,jk,Kmm) = ( puu(:,:,jk,Kmm) + un_adv(:,:)*r1_hu_0(:,:) - puu_b(:,:,Kmm) ) * umask(:,:,jk)
         pvv(:,:,jk,Kmm) = ( pvv(:,:,jk,Kmm) + vn_adv(:,:)*r1_hv_0(:,:) - pvv_b(:,:,Kmm) ) * vmask(:,:,jk)
      END DO

      IF ( ln_wd_dl .and. ln_wd_dl_bc) THEN 
         DO jk = 1, jpkm1
            puu(:,:,jk,Kmm) = ( un_adv(:,:)*r1_hu_0(:,:) &
                       & + zuwdav2(:,:)*(puu(:,:,jk,Kmm) - un_adv(:,:)*r1_hu_0(:,:)) ) * umask(:,:,jk) 
            pvv(:,:,jk,Kmm) = ( vn_adv(:,:)*r1_hv_0(:,:) & 
                       & + zvwdav2(:,:)*(pvv(:,:,jk,Kmm) - vn_adv(:,:)*r1_hv_0(:,:)) ) * vmask(:,:,jk)  
         END DO
      END IF 

      
      CALL iom_put(  "ubar", un_adv(:,:)*r1_hu_0(:,:) )    ! barotropic i-current
      CALL iom_put(  "vbar", vn_adv(:,:)*r1_hv_0(:,:) )    ! barotropic i-current
      !
      !                                   !* write time-spliting arrays in the restart
      IF( lrst_oce .AND.ln_bt_fw )   CALL ts_rst( kt, 'WRITE' )
      !
      IF( ln_wd_dl )   DEALLOCATE( ztwdmask, zuwdmask, zvwdmask, zuwdav2, zvwdav2 )
      !
      CALL iom_put( "baro_u" , puu_b(:,:,Kmm) )  ! Barotropic  U Velocity
      CALL iom_put( "baro_v" , pvv_b(:,:,Kmm) )  ! Barotropic  V Velocity
      !
   END SUBROUTINE dyn_spg_ts


   SUBROUTINE ts_wgt( ll_av, ll_fw, jpit, zwgt1, zwgt2)
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_wgt  ***
      !!
      !! ** Purpose : Set time-splitting weights for temporal averaging (or not)
      !!----------------------------------------------------------------------
      LOGICAL, INTENT(in) ::   ll_av      ! temporal averaging=.true.
      LOGICAL, INTENT(in) ::   ll_fw      ! forward time splitting =.true.
      INTEGER, INTENT(inout) :: jpit      ! cycle length    
      REAL(wp), DIMENSION(3*nn_e), INTENT(inout)  :: zwgt2
      REAL(dp), DIMENSION(3*nn_e), INTENT(inout)  :: zwgt1

      
      INTEGER ::  jic, jn, ji                      ! temporary integers
      REAL(wp)  :: za2
      REAL(dp)  :: za1
      !!----------------------------------------------------------------------

      zwgt1(:) = 0._wp
      zwgt2(:) = 0._wp

      ! Set time index when averaged value is requested
      IF (ll_fw) THEN 
         jic = nn_e
      ELSE
         jic = 2 * nn_e
      ENDIF

      ! Set primary weights:
      IF (ll_av) THEN
           ! Define simple boxcar window for primary weights 
           ! (width = nn_e, centered around jic)     
         SELECT CASE ( nn_bt_flt )
              CASE( 0 )  ! No averaging
                 zwgt1(jic) = 1._wp
                 jpit = jic

              CASE( 1 )  ! Boxcar, width = nn_e
                 DO jn = 1, 3*nn_e
                    za1 = ABS(float(jn-jic))/float(nn_e) 
                    IF (za1 < 0.5_wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO

              CASE( 2 )  ! Boxcar, width = 2 * nn_e
                 DO jn = 1, 3*nn_e
                    za1 = ABS(float(jn-jic))/float(nn_e) 
                    IF (za1 < 1._wp) THEN
                      zwgt1(jn) = 1._wp
                      jpit = jn
                    ENDIF
                 ENDDO
              CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt' )
         END SELECT

      ELSE ! No time averaging
         zwgt1(jic) = 1._wp
         jpit = jic
      ENDIF
    
      ! Set secondary weights
      DO jn = 1, jpit
        DO ji = jn, jpit
             zwgt2(jn) = zwgt2(jn) + zwgt1(ji)
        END DO
      END DO

      ! Normalize weigths:
      za1 = 1._wp / SUM(zwgt1(1:jpit))
      za2 = 1._wp / SUM(zwgt2(1:jpit))
      DO jn = 1, jpit
        zwgt1(jn) = zwgt1(jn) * za1
        zwgt2(jn) = zwgt2(jn) * za2
      END DO
      !
   END SUBROUTINE ts_wgt


   SUBROUTINE ts_rst( kt, cdrw )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ts_rst  ***
      !!
      !! ** Purpose : Read or write time-splitting arrays in restart file
      !!----------------------------------------------------------------------
      INTEGER         , INTENT(in) ::   kt     ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw   ! "READ"/"WRITE" flag
      !!----------------------------------------------------------------------
      !
      IF( TRIM(cdrw) == 'READ' ) THEN        ! Read/initialise 
         !                                   ! ---------------
         IF( ln_rstart .AND. ln_bt_fw .AND. .NOT.l_1st_euler ) THEN    !* Read the restart file
            CALL iom_get( numror, jpdom_auto, 'ub2_b'  , ub2_b  (:,:), cd_type = 'U', psgn = -1._wp )   
            CALL iom_get( numror, jpdom_auto, 'vb2_b'  , vb2_b  (:,:), cd_type = 'V', psgn = -1._wp ) 
            CALL iom_get( numror, jpdom_auto, 'un_bf'  , un_bf  (:,:), cd_type = 'U', psgn = -1._wp )   
            CALL iom_get( numror, jpdom_auto, 'vn_bf'  , vn_bf  (:,:), cd_type = 'V', psgn = -1._wp ) 
            IF( .NOT.ln_bt_av ) THEN
               CALL iom_get( numror, jpdom_auto, 'sshbb_e'  , sshbb_e(:,:), cd_type = 'T', psgn =  1._wp )   
               CALL iom_get( numror, jpdom_auto, 'ubb_e'    ,   ubb_e(:,:), cd_type = 'U', psgn = -1._wp )   
               CALL iom_get( numror, jpdom_auto, 'vbb_e'    ,   vbb_e(:,:), cd_type = 'V', psgn = -1._wp )
               CALL iom_get( numror, jpdom_auto, 'sshb_e'   ,  sshb_e(:,:), cd_type = 'T', psgn =  1._wp ) 
               CALL iom_get( numror, jpdom_auto, 'ub_e'     ,    ub_e(:,:), cd_type = 'U', psgn = -1._wp )   
               CALL iom_get( numror, jpdom_auto, 'vb_e'     ,    vb_e(:,:), cd_type = 'V', psgn = -1._wp )
            ENDIF
         ELSE                                   !* Start from rest
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '   ==>>>   start from rest: set barotropic values to 0'
            ub2_b  (:,:) = 0._wp   ;   vb2_b  (:,:) = 0._wp   ! used in the 1st interpol of agrif
            un_adv (:,:) = 0._wp   ;   vn_adv (:,:) = 0._wp   ! used in the 1st interpol of agrif
            un_bf  (:,:) = 0._wp   ;   vn_bf  (:,:) = 0._wp   ! used in the 1st update   of agrif
         ENDIF
         !
      ELSEIF( TRIM(cdrw) == 'WRITE' ) THEN   ! Create restart file
         !                                   ! -------------------
         IF(lwp) WRITE(numout,*) '---- ts_rst ----'
         CALL iom_rstput( kt, nitrst, numrow, 'ub2_b'   , ub2_b  (:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'vb2_b'   , vb2_b  (:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'un_bf'   , un_bf  (:,:) )
         CALL iom_rstput( kt, nitrst, numrow, 'vn_bf'   , vn_bf  (:,:) )
         !
         IF (.NOT.ln_bt_av) THEN
            CALL iom_rstput( kt, nitrst, numrow, 'sshbb_e'  , sshbb_e(:,:) ) 
            CALL iom_rstput( kt, nitrst, numrow, 'ubb_e'    ,   ubb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vbb_e'    ,   vbb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'sshb_e'   ,  sshb_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'ub_e'     ,    ub_e(:,:) )
            CALL iom_rstput( kt, nitrst, numrow, 'vb_e'     ,    vb_e(:,:) )
         ENDIF
      ENDIF
      !
   END SUBROUTINE ts_rst


   SUBROUTINE dyn_spg_ts_init
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_spg_ts_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj              ! dummy loop indices
      REAL(wp) ::   zxr2, zyr2, zcmax   ! local scalar
      REAL(wp), DIMENSION(jpi,jpj) ::   zcu
      !!----------------------------------------------------------------------
      !
      ! Max courant number for ext. grav. waves
      !
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
         zxr2 = r1_e1t(ji,jj) * r1_e1t(ji,jj)
         zyr2 = r1_e2t(ji,jj) * r1_e2t(ji,jj)
         zcu(ji,jj) = SQRT( grav * MAX(ht_0(ji,jj),0._wp) * (zxr2 + zyr2) )
      END DO   ;   END DO
      !
      !
      zcmax = MAXVAL( zcu(Nis0:Nie0,Njs0:Nje0) )
      CALL mpp_max( 'dynspg_ts', zcmax )

      ! Estimate number of iterations to satisfy a max courant number= rn_bt_cmax
      IF( ln_bt_auto )   nn_e = CEILING( rn_Dt / rn_bt_cmax * zcmax)
      
      rDt_e = rn_Dt / REAL( nn_e , dp )
      zcmax = zcmax * rDt_e
      ! Print results
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'dyn_spg_ts_init : split-explicit free surface'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      IF( ln_bt_auto ) THEN
         IF(lwp) WRITE(numout,*) '     ln_ts_auto =.true. Automatically set nn_e '
         IF(lwp) WRITE(numout,*) '     Max. courant number allowed: ', rn_bt_cmax
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_ts_auto=.false.: Use nn_e in namelist   nn_e = ', nn_e
      ENDIF

      IF(ln_bt_av) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_av =.true.  ==> Time averaging over nn_e time steps is on '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_av =.false. => No time averaging of barotropic variables '
      ENDIF
      !
      !
      IF(ln_bt_fw) THEN
         IF(lwp) WRITE(numout,*) '     ln_bt_fw=.true.  => Forward integration of barotropic variables '
      ELSE
         IF(lwp) WRITE(numout,*) '     ln_bt_fw =.false.=> Centred integration of barotropic variables '
      ENDIF
      !
      !
      IF(lwp) WRITE(numout,*)    '     Time filter choice, nn_bt_flt: ', nn_bt_flt
      SELECT CASE ( nn_bt_flt )
         CASE( 0 )      ;   IF(lwp) WRITE(numout,*) '           Dirac'
         CASE( 1 )      ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = nn_e'
         CASE( 2 )      ;   IF(lwp) WRITE(numout,*) '           Boxcar: width = 2*nn_e' 
         CASE DEFAULT   ;   CALL ctl_stop( 'unrecognised value for nn_bt_flt: should 0,1, or 2' )
      END SELECT
      !
      IF(lwp) WRITE(numout,*) ' '
      IF(lwp) WRITE(numout,*) '     nn_e = ', nn_e
      IF(lwp) WRITE(numout,*) '     Barotropic time step [s] is :', rDt_e
      IF(lwp) WRITE(numout,*) '     Maximum Courant number is   :', zcmax
      !
      IF(lwp) WRITE(numout,*)    '     Time diffusion parameter rn_bt_alpha: ', rn_bt_alpha
      IF ((ln_bt_av.AND.nn_bt_flt/=0).AND.(rn_bt_alpha>0._wp)) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: if rn_bt_alpha > 0, remove temporal averaging' )
      ENDIF
      !
      IF( .NOT.ln_bt_av .AND. .NOT.ln_bt_fw ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: No time averaging => only forward integration is possible' )
      ENDIF
      IF( zcmax>0.9_wp ) THEN
         CALL ctl_stop( 'dynspg_ts ERROR: Maximum Courant number is greater than 0.9: Inc. nn_e !' )          
      ENDIF
      !
      !                             ! Allocate time-splitting arrays
      IF( dyn_spg_ts_alloc() /= 0    )   CALL ctl_stop('STOP', 'dyn_spg_init: failed to allocate dynspg_ts  arrays' )
      !
      !                             ! read restart when needed
      CALL ts_rst( nit000, 'READ' )
      !
   END SUBROUTINE dyn_spg_ts_init

   
   SUBROUTINE dyn_cor_2D_init( Kmm )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_2D_init  ***
      !!
      !! ** Purpose : Set time splitting options
      !! Set arrays to remove/compute coriolis trend.
      !! Do it once during initialization if volume is fixed, else at each long time step.
      !! Note that these arrays are also used during barotropic loop. These are however frozen
      !! although they should be updated in the variable volume case. Not a big approximation.
      !! To remove this approximation, copy lines below inside barotropic loop
      !! and update depths at T- points (ht) at each barotropic time step
      !!
      !! Compute zwz = f / ( height of the water colomn )
      !!----------------------------------------------------------------------
      INTEGER,  INTENT(in)         ::  Kmm  ! Time index
      INTEGER  ::   ji ,jj, jk              ! dummy loop indices
      REAL(wp) ::   z1_ht
      !!----------------------------------------------------------------------
      !
      SELECT CASE( nvor_scheme )
      CASE( np_EEN, np_ENE, np_ENS , np_MIX )   !=  schemes using the same e3f definition
         SELECT CASE( nn_e3f_typ )                  !* ff_f/e3 at F-point
         CASE ( 0 )                                   ! original formulation  (masked averaging of e3t divided by 4)
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               zwz(ji,jj) = ( ht_0(ji,jj+1) + ht_0(ji+1,jj+1)   &
                    &       + ht_0(ji,jj  ) + ht_0(ji+1,jj  ) ) * 0.25_wp  
               IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = ff_f(ji,jj) / zwz(ji,jj)
            END DO   ;   END DO
         CASE ( 1 )                                   ! new formulation  (masked averaging of e3t divided by the sum of mask)
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               zwz(ji,jj) =     (    ht_0(ji,jj+1) +     ht_0(ji+1,jj+1)      &
                    &            +   ht_0(ji,jj  ) +     ht_0(ji+1,jj  )  )   &
                    &    / ( MAX(ssmask(ji,jj+1) + ssmask(ji+1,jj+1)      &
                    &          + ssmask(ji,jj  ) + ssmask(ji+1,jj  ) , 1._wp  )   )
               IF( zwz(ji,jj) /= 0._wp )   zwz(ji,jj) = ff_f(ji,jj) / zwz(ji,jj)
            END DO   ;   END DO
         END SELECT
         CALL lbc_lnk( 'dynspg_ts', zwz, 'F', 1._wp )
      END SELECT
      !
      SELECT CASE( nvor_scheme )
      CASE( np_EEN )
         !
         ftne(1,:) = 0._wp   ;   ftnw(1,:) = 0._wp   ;   ftse(1,:) = 0._wp   ;   ftsw(1,:) = 0._wp
         DO jj = ntsj-( 0), ntej+( 1 ) ; DO ji = ntsi-( 0), ntei+( 1)
            ftne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
            ftnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
            ftse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
            ftsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
         END DO   ;   END DO
         !
      CASE( np_EET )                            != EEN scheme using e3t energy conserving scheme
         ftne(1,:) = 0._wp   ;   ftnw(1,:) = 0._wp   ;   ftse(1,:) = 0._wp   ;   ftsw(1,:) = 0._wp
         DO jj = ntsj-( 0), ntej+( 1 ) ; DO ji = ntsi-( 0), ntei+( 1)
            z1_ht = ssmask(ji,jj) / ( ht_0(ji,jj) + 1._wp - ssmask(ji,jj) )
            ftne(ji,jj) = ( ff_f(ji-1,jj  ) + ff_f(ji  ,jj  ) + ff_f(ji  ,jj-1) ) * z1_ht
            ftnw(ji,jj) = ( ff_f(ji-1,jj-1) + ff_f(ji-1,jj  ) + ff_f(ji  ,jj  ) ) * z1_ht
            ftse(ji,jj) = ( ff_f(ji  ,jj  ) + ff_f(ji  ,jj-1) + ff_f(ji-1,jj-1) ) * z1_ht
            ftsw(ji,jj) = ( ff_f(ji  ,jj-1) + ff_f(ji-1,jj-1) + ff_f(ji-1,jj  ) ) * z1_ht
         END DO   ;   END DO
         !
      END SELECT
      
   END SUBROUTINE dyn_cor_2d_init


   SUBROUTINE dyn_cor_2d( pht, phu, phv, punb, pvnb, zhU, zhV,    zu_trd, zv_trd   )
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE dyn_cor_2d  ***
      !!
      !! ** Purpose : Compute u and v coriolis trends
      !!----------------------------------------------------------------------
      INTEGER  ::   ji ,jj                             ! dummy loop indices
      REAL(wp) ::   zx1, zx2, zy1, zy2, z1_hu, z1_hv   !   -      -
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   )  :: pht, phu, phv, punb, pvnb, zhV
      REAL(dp), DIMENSION(jpi,jpj), INTENT(in   )  :: zhU
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) :: zu_trd, zv_trd
      !!----------------------------------------------------------------------
      SELECT CASE( nvor_scheme )
      CASE( np_ENT )                ! enstrophy conserving scheme (f-point)
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            z1_hu = ssumask(ji,jj) / ( phu(ji,jj) + 1._wp - ssumask(ji,jj) )
            z1_hv = ssvmask(ji,jj) / ( phv(ji,jj) + 1._wp - ssvmask(ji,jj) )
            zu_trd(ji,jj) = + r1_4 * r1_e1e2u(ji,jj) * z1_hu                    &
               &               * (  e1e2t(ji+1,jj)*pht(ji+1,jj)*ff_t(ji+1,jj) * ( pvnb(ji+1,jj) + pvnb(ji+1,jj-1) )   &
               &                  + e1e2t(ji  ,jj)*pht(ji  ,jj)*ff_t(ji  ,jj) * ( pvnb(ji  ,jj) + pvnb(ji  ,jj-1) )   )
               !
            zv_trd(ji,jj) = - r1_4 * r1_e1e2v(ji,jj) * z1_hv                    &
               &               * (  e1e2t(ji,jj+1)*pht(ji,jj+1)*ff_t(ji,jj+1) * ( punb(ji,jj+1) + punb(ji-1,jj+1) )   & 
               &                  + e1e2t(ji,jj  )*pht(ji,jj  )*ff_t(ji,jj  ) * ( punb(ji,jj  ) + punb(ji-1,jj  ) )   ) 
         END DO   ;   END DO
         !         
      CASE( np_ENE , np_MIX )        ! energy conserving scheme (t-point) ENE or MIX
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zy1 = ( zhV(ji,jj-1) + zhV(ji+1,jj-1) ) * r1_e1u(ji,jj)
            zy2 = ( zhV(ji,jj  ) + zhV(ji+1,jj  ) ) * r1_e1u(ji,jj)
            zx1 = ( zhU(ji-1,jj) + zhU(ji-1,jj+1) ) * r1_e2v(ji,jj)
            zx2 = ( zhU(ji  ,jj) + zhU(ji  ,jj+1) ) * r1_e2v(ji,jj)
            ! energy conserving formulation for planetary vorticity term
            zu_trd(ji,jj) =   r1_4 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
            zv_trd(ji,jj) = - r1_4 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
         END DO   ;   END DO
         !
      CASE( np_ENS )                ! enstrophy conserving scheme (f-point)
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zy1 =   r1_8 * ( zhV(ji  ,jj-1) + zhV(ji+1,jj-1) &
              &            + zhV(ji  ,jj  ) + zhV(ji+1,jj  ) ) * r1_e1u(ji,jj)
            zx1 = - r1_8 * ( zhU(ji-1,jj  ) + zhU(ji-1,jj+1) &
              &            + zhU(ji  ,jj  ) + zhU(ji  ,jj+1) ) * r1_e2v(ji,jj)
            zu_trd(ji,jj)  = zy1 * ( zwz(ji  ,jj-1) + zwz(ji,jj) )
            zv_trd(ji,jj)  = zx1 * ( zwz(ji-1,jj  ) + zwz(ji,jj) )
         END DO   ;   END DO
         !
      CASE( np_EET , np_EEN )      ! energy & enstrophy scheme (using e3t or e3f)         
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zu_trd(ji,jj) = + r1_12 * r1_e1u(ji,jj) * (  ftne(ji,jj  ) * zhV(ji  ,jj  ) &
             &                                         + ftnw(ji+1,jj) * zhV(ji+1,jj  ) &
             &                                         + ftse(ji,jj  ) * zhV(ji  ,jj-1) &
             &                                         + ftsw(ji+1,jj) * zhV(ji+1,jj-1) )
            zv_trd(ji,jj) = - r1_12 * r1_e2v(ji,jj) * (  ftsw(ji,jj+1) * zhU(ji-1,jj+1) &
             &                                         + ftse(ji,jj+1) * zhU(ji  ,jj+1) &
             &                                         + ftnw(ji,jj  ) * zhU(ji-1,jj  ) &
             &                                         + ftne(ji,jj  ) * zhU(ji  ,jj  ) )
         END DO   ;   END DO
         !
      END SELECT
      !
   END SUBROUTINE dyn_cor_2D


   SUBROUTINE wad_tmsk( pssh, ptmsk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_tmsk  ***
      !!                    
      !! ** Purpose :   set wetting & drying mask at tracer points 
      !!              for the current barotropic sub-step 
      !!
      !! ** Method  :   ??? 
      !!
      !! ** Action  :  ptmsk : wetting & drying t-mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pssh    !
      REAL(wp), DIMENSION(jpi,jpj), INTENT(  out) ::   ptmsk   !
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( ln_wd_dl_rmp ) THEN     
         DO jj = ntsj-( 1), ntej+( 1 ) ; DO ji = ntsi-( 1), ntei+( 1)
            IF    ( pssh(ji,jj) + ht_0(ji,jj) >  2._wp * rn_wdmin1 ) THEN 
               !           IF    ( pssh(ji,jj) + ht_0(ji,jj) >          rn_wdmin2 ) THEN 
               ptmsk(ji,jj) = 1._wp
            ELSEIF( pssh(ji,jj) + ht_0(ji,jj) >          rn_wdmin1 ) THEN
               ptmsk(ji,jj) = TANH( 50._wp*( ( pssh(ji,jj) + ht_0(ji,jj) -  rn_wdmin1 )*r_rn_wdmin1) )
            ELSE 
               ptmsk(ji,jj) = 0._wp
            ENDIF
         END DO   ;   END DO
      ELSE  
         DO jj = ntsj-( 1), ntej+( 1 ) ; DO ji = ntsi-( 1), ntei+( 1)
            IF ( pssh(ji,jj) + ht_0(ji,jj) >  rn_wdmin1 ) THEN   ;   ptmsk(ji,jj) = 1._wp
            ELSE                                                 ;   ptmsk(ji,jj) = 0._wp
            ENDIF
         END DO   ;   END DO
      ENDIF
      !
   END SUBROUTINE wad_tmsk


   SUBROUTINE wad_Umsk( pTmsk, phU, phV, pu, pv, pUmsk, pVmsk )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE wad_Umsk  ***
      !!                    
      !! ** Purpose :   set wetting & drying mask at tracer points 
      !!              for the current barotropic sub-step 
      !!
      !! ** Method  :   ??? 
      !!
      !! ** Action  :  ptmsk : wetting & drying t-mask
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in   ) ::   pTmsk              ! W & D t-mask
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout)  :: phV, pu, pv! ocean velocities and transports
      REAL(dp), DIMENSION(jpi,jpj), INTENT(inout)  :: phU! ocean velocities and transports
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pUmsk, pVmsk       ! W & D u- and v-mask
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      DO jj = ntsj-( 1), ntej+( 1 ) ; DO ji = ntsi-( 1), ntei+( 0)   ! not jpi-column
         IF ( phU(ji,jj) > 0._wp ) THEN   ;   pUmsk(ji,jj) = pTmsk(ji  ,jj) 
         ELSE                             ;   pUmsk(ji,jj) = pTmsk(ji+1,jj)  
         ENDIF
         phU(ji,jj) = pUmsk(ji,jj)*phU(ji,jj)
         pu (ji,jj) = pUmsk(ji,jj)*pu (ji,jj)
      END DO   ;   END DO
      !
      DO jj = ntsj-( 1), ntej+( 0 ) ; DO ji = ntsi-( 1), ntei+( 1)   ! not jpj-row
         IF ( phV(ji,jj) > 0._wp ) THEN   ;   pVmsk(ji,jj) = pTmsk(ji,jj  )
         ELSE                             ;   pVmsk(ji,jj) = pTmsk(ji,jj+1)  
         ENDIF
         phV(ji,jj) = pVmsk(ji,jj)*phV(ji,jj) 
         pv (ji,jj) = pVmsk(ji,jj)*pv (ji,jj)
      END DO   ;   END DO
      !
   END SUBROUTINE wad_Umsk

   SUBROUTINE dyn_drg_init( Kbb, Kmm, puu, pvv, puu_b ,pvv_b, pu_RHSi, pv_RHSi, pCdU_u, pCdU_v )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_drg_init  ***
      !!                    
      !! ** Purpose : - add the baroclinic top/bottom drag contribution to 
      !!              the baroclinic part of the barotropic RHS
      !!              - compute the barotropic drag coefficients
      !!
      !! ** Method  :   computation done over the INNER domain only 
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT(in   ) ::  Kbb, Kmm           ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(in   ) ::  puu, pvv           ! ocean velocities and RHS of momentum equation
      REAL(wp), DIMENSION(jpi,jpj,jpt)    , INTENT(in   ) ::  puu_b, pvv_b       ! barotropic velocities at main time levels
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(inout) ::  pu_RHSi, pv_RHSi   ! baroclinic part of the barotropic RHS
      REAL(wp), DIMENSION(jpi,jpj)        , INTENT(  out) ::  pCdU_u , pCdU_v    ! barotropic drag coefficients
      !
      INTEGER  ::   ji, jj   ! dummy loop indices
      INTEGER  ::   ikbu, ikbv, iktu, iktv
      REAL(wp) ::   zztmp
      REAL(wp), DIMENSION(jpi,jpj) ::   zu_i, zv_i
      !!----------------------------------------------------------------------
      !
      !                    !==  Set the barotropic drag coef.  ==!
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN          ! top+bottom friction (ocean cavities)
         
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            pCdU_u(ji,jj) = r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) + rCdU_top(ji+1,jj)+rCdU_top(ji,jj) )
            pCdU_v(ji,jj) = r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) + rCdU_top(ji,jj+1)+rCdU_top(ji,jj) )
         END DO   ;   END DO
      ELSE                          ! bottom friction only
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            pCdU_u(ji,jj) = r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) )
            pCdU_v(ji,jj) = r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) )
         END DO   ;   END DO
      ENDIF
      !
      !                    !==  BOTTOM stress contribution from baroclinic velocities  ==!
      !
      IF( ln_bt_fw ) THEN                 ! FORWARD integration: use NOW bottom baroclinic velocities
         
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            ikbu = mbku(ji,jj)       
            ikbv = mbkv(ji,jj)    
            zu_i(ji,jj) = puu(ji,jj,ikbu,Kmm) - puu_b(ji,jj,Kmm)
            zv_i(ji,jj) = pvv(ji,jj,ikbv,Kmm) - pvv_b(ji,jj,Kmm)
         END DO   ;   END DO
      ELSE                                ! CENTRED integration: use BEFORE bottom baroclinic velocities
         
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            ikbu = mbku(ji,jj)       
            ikbv = mbkv(ji,jj)    
            zu_i(ji,jj) = puu(ji,jj,ikbu,Kbb) - puu_b(ji,jj,Kbb)
            zv_i(ji,jj) = pvv(ji,jj,ikbv,Kbb) - pvv_b(ji,jj,Kbb)
         END DO   ;   END DO
      ENDIF
      !
      ! use "unclipped" drag (even if explicit friction is used in 3D calculation)
         
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
         pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + r1_hu_0(ji,jj) * r1_2*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) * zu_i(ji,jj)
         pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + r1_hv_0(ji,jj) * r1_2*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) * zv_i(ji,jj)
      END DO   ;   END DO
      !
      !                    !==  TOP stress contribution from baroclinic velocities  ==!   (no W/D case)
      !
      IF( ln_isfcav.OR.ln_drgice_imp ) THEN
         !
         IF( ln_bt_fw ) THEN                ! FORWARD integration: use NOW top baroclinic velocity
            
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               iktu = miku(ji,jj)
               iktv = mikv(ji,jj)
               zu_i(ji,jj) = puu(ji,jj,iktu,Kmm) - puu_b(ji,jj,Kmm)
               zv_i(ji,jj) = pvv(ji,jj,iktv,Kmm) - pvv_b(ji,jj,Kmm)
            END DO   ;   END DO
         ELSE                                ! CENTRED integration: use BEFORE top baroclinic velocity
            
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               iktu = miku(ji,jj)
               iktv = mikv(ji,jj)
               zu_i(ji,jj) = puu(ji,jj,iktu,Kbb) - puu_b(ji,jj,Kbb)
               zv_i(ji,jj) = pvv(ji,jj,iktv,Kbb) - pvv_b(ji,jj,Kbb)
            END DO   ;   END DO
         ENDIF
         !
         !                    ! use "unclipped" top drag (even if explicit friction is used in 3D calculation)
         
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            pu_RHSi(ji,jj) = pu_RHSi(ji,jj) + r1_hu_0(ji,jj) * r1_2*( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) * zu_i(ji,jj)
            pv_RHSi(ji,jj) = pv_RHSi(ji,jj) + r1_hv_0(ji,jj) * r1_2*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) ) * zv_i(ji,jj)
         END DO   ;   END DO
         !
      ENDIF
      !
   END SUBROUTINE dyn_drg_init

   SUBROUTINE ts_bck_interp( jn, ll_init,       &   ! <== in
      &                      za0, za1, za2, za3 )   ! ==> out
      !!----------------------------------------------------------------------
      INTEGER ,INTENT(in   ) ::   jn                   ! index of sub time step
      LOGICAL ,INTENT(in   ) ::   ll_init              !
      REAL(wp),INTENT(  out)  :: za0, za2, za3! Half-step back interpolation coefficient
      REAL(dp),INTENT(  out)  :: za1! Half-step back interpolation coefficient
      !
      REAL(wp) ::   zepsilon, zgamma                   !   -      -
      !!----------------------------------------------------------------------
      !                             ! set Half-step back interpolation coefficient
      IF    ( jn==1 .AND. ll_init ) THEN   !* Forward-backward
         za0 = 1._wp                        
         za1 = 0._wp                           
         za2 = 0._wp
         za3 = 0._wp
      ELSEIF( jn==2 .AND. ll_init ) THEN   !* AB2-AM3 Coefficients; bet=0 ; gam=-1/6 ; eps=1/12
         za0 = 1.0833333333333_wp                 ! za0 = 1-gam-eps
         za1 =-0.1666666666666_wp                 ! za1 = gam
         za2 = 0.0833333333333_wp                 ! za2 = eps
         za3 = 0._wp              
      ELSE                                 !* AB3-AM4 Coefficients; bet=0.281105 ; eps=0.013 ; gam=0.0880 
         IF( rn_bt_alpha == 0._wp ) THEN      ! Time diffusion  
            za0 = 0.614_wp                        ! za0 = 1/2 +   gam + 2*eps
            za1 = 0.285_wp                        ! za1 = 1/2 - 2*gam - 3*eps
            za2 = 0.088_wp                        ! za2 = gam
            za3 = 0.013_wp                        ! za3 = eps
         ELSE                                 ! no time diffusion
            zepsilon = 0.00976186_wp - 0.13451357_wp * rn_bt_alpha
            zgamma   = 0.08344500_wp - 0.51358400_wp * rn_bt_alpha
            za0 = 0.5_wp + zgamma + 2._wp * rn_bt_alpha + 2._wp * zepsilon
            za1 = 1._wp - za0 - zgamma - zepsilon
            za2 = zgamma
            za3 = zepsilon
         ENDIF 
      ENDIF
   END SUBROUTINE ts_bck_interp


   !!======================================================================
END MODULE dynspg_ts












MODULE dynzdf
   !!==============================================================================
   !!                 ***  MODULE  dynzdf  ***
   !! Ocean dynamics :  vertical component of the momentum mixing trend
   !!==============================================================================
   !! History :  1.0  !  2005-11  (G. Madec)  Original code
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            4.0  !  2017-06  (G. Madec) remove the explicit time-stepping option + avm at t-point
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf       : compute the after velocity through implicit calculation of vertical mixing
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers variables
   USE phycst         ! physical constants
   USE dom_oce        ! ocean space and time domain variables 
   USE sbc_oce        ! surface boundary condition: ocean
   USE zdf_oce        ! ocean vertical physics variables
   USE zdfdrg         ! vertical physics: top/bottom drag coef.
   USE dynadv    ,ONLY: ln_dynadv_vec    ! dynamics: advection form
   USE dynldf_iso,ONLY: akzu, akzv       ! dynamics: vertical component of rotated lateral mixing 
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef. and type of operator
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf   !  routine called by step.F90

   REAL(wp) ::  r_vvl     ! non-linear free surface indicator: =0 if ln_linssh=T, =1 otherwise 

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
   !! $Id: dynzdf.F90 14834 2021-05-11 09:24:44Z hadcv $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS
   
   SUBROUTINE dyn_zdf( kt, Kbb, Kmm, Krhs, puu, pvv, Kaa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf  ***
      !!
      !! ** Purpose :   compute the trend due to the vert. momentum diffusion
      !!              together with the Leap-Frog time stepping using an 
      !!              implicit scheme.
      !!
      !! ** Method  :  - Leap-Frog time stepping on all trends but the vertical mixing
      !!         u(after) =         u(before) + 2*dt *       u(rhs)                vector form or linear free surf.
      !!         u(after) = ( e3u_b*u(before) + 2*dt * e3u_n*u(rhs) ) / e3u_after   otherwise
      !!               - update the after velocity with the implicit vertical mixing.
      !!      This requires to solver the following system: 
      !!         u(after) = u(after) + 1/e3u_after  dk+1[ mi(avm) / e3uw_after dk[ua] ]
      !!      with the following surface/top/bottom boundary condition:
      !!      surface: wind stress input (averaged over kt-1/2 & kt+1/2)
      !!      top & bottom : top stress (iceshelf-ocean) & bottom stress (cf zdfdrg.F90)
      !!
      !! ** Action :   (puu(:,:,:,Kaa),pvv(:,:,:,Kaa))   after velocity 
      !!---------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt                  ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kbb, Kmm, Krhs, Kaa ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv            ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk         ! dummy loop indices
      INTEGER  ::   iku, ikv           ! local integers
      REAL(wp) ::   zzwi, ze3ua, zdt   ! local scalars
      REAL(wp) ::   zzws, ze3va        !   -      -
      REAL(wp) ::   z1_e3ua, z1_e3va   !   -      -
      REAL(wp) ::   zWu , zWv          !   -      -
      REAL(wp) ::   zWui, zWvi         !   -      -
      REAL(wp) ::   zWus, zWvs         !   -      -
      REAL(wp), DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls),jpk)        ::  zwi, zwd, zws   ! 3D workspace
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE ::   ztrdu, ztrdv   !  -      -
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_zdf')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN       !* initialization
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_zdf_imp : vertical momentum diffusion implicit operator'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
            !
            If( ln_linssh ) THEN   ;    r_vvl = 0._wp    ! non-linear free surface indicator
            ELSE                   ;    r_vvl = 1._wp
            ENDIF
         ENDIF
      ENDIF
      !                             !* explicit top/bottom drag case
      IF( .NOT.ln_drgimp )   CALL zdf_drg_exp( kt, Kmm, puu(:,:,:,Kbb), pvv(:,:,:,Kbb), puu(:,:,:,Krhs), pvv(:,:,:,Krhs) )  ! add top/bottom friction trend to (puu(Kaa),pvv(Kaa))
      !
      !
      IF( l_trddyn )   THEN         !* temporary save of ta and sa trends
         ALLOCATE( ztrdu(jpi,jpj,jpk), ztrdv(jpi,jpj,jpk) ) 
         ztrdu(:,:,:) = puu(:,:,:,Krhs)
         ztrdv(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF
      !
      !              !==  RHS: Leap-Frog time stepping on all trends but the vertical mixing  ==!   (put in puu(:,:,:,Kaa),pvv(:,:,:,Kaa))
      !
      !                    ! time stepping except vertical diffusion
      IF( ln_dynadv_vec .OR. ln_linssh ) THEN   ! applied on velocity
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            puu(ji,jj,jk,Kaa) = ( puu(ji,jj,jk,Kbb) + rDt * puu(ji,jj,jk,Krhs) ) * umask(ji,jj,jk)
            pvv(ji,jj,jk,Kaa) = ( pvv(ji,jj,jk,Kbb) + rDt * pvv(ji,jj,jk,Krhs) ) * vmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
      ELSE                                      ! applied on thickness weighted velocity
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
            puu(ji,jj,jk,Kaa) = (         e3u_0(ji,jj,jk) * puu(ji,jj,jk,Kbb )  &
               &                  + rDt * e3u_0(ji,jj,jk) * puu(ji,jj,jk,Krhs)  ) &
               &                        / e3u_0(ji,jj,jk) * umask(ji,jj,jk)
            pvv(ji,jj,jk,Kaa) = (         e3v_0(ji,jj,jk) * pvv(ji,jj,jk,Kbb )  &
               &                  + rDt * e3v_0(ji,jj,jk) * pvv(ji,jj,jk,Krhs)  ) &
               &                        / e3v_0(ji,jj,jk) * vmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
      ENDIF
      !                    ! add top/bottom friction 
      !     With split-explicit free surface, barotropic stress is treated explicitly Update velocities at the bottom.
      !     J. Chanut: The bottom stress is computed considering after barotropic velocities, which does 
      !                not lead to the effective stress seen over the whole barotropic loop. 
      !     G. Madec : in linear free surface, e3u_0(:,:,:) = e3u_0(:,:,:) = e3u_0, so systematic use of e3u_0(:,:,:)
      IF( ln_drgimp .AND. ln_dynspg_ts ) THEN
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)      ! remove barotropic velocities
            puu(ji,jj,jk,Kaa) = ( puu(ji,jj,jk,Kaa) - uu_b(ji,jj,Kaa) ) * umask(ji,jj,jk)
            pvv(ji,jj,jk,Kaa) = ( pvv(ji,jj,jk,Kaa) - vv_b(ji,jj,Kaa) ) * vmask(ji,jj,jk)
         END DO   ;   END DO   ;   END DO
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)      ! Add bottom/top stress due to barotropic component only
            iku = mbku(ji,jj)         ! ocean bottom level at u- and v-points 
            ikv = mbkv(ji,jj)         ! (deepest ocean u- and v-points)
            ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,iku)    &
               &             + r_vvl   * e3u_0(ji,jj,iku)
            ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,ikv)    &
               &             + r_vvl   * e3v_0(ji,jj,ikv)
            puu(ji,jj,iku,Kaa) = puu(ji,jj,iku,Kaa) + rDt * 0.5*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) * uu_b(ji,jj,Kaa) / ze3ua
            pvv(ji,jj,ikv,Kaa) = pvv(ji,jj,ikv,Kaa) + rDt * 0.5*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) * vv_b(ji,jj,Kaa) / ze3va
         END DO   ;   END DO
         IF( ln_isfcav.OR.ln_drgice_imp ) THEN    ! Ocean cavities (ISF)
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               iku = miku(ji,jj)         ! top ocean level at u- and v-points 
               ikv = mikv(ji,jj)         ! (first wet ocean u- and v-points)
               ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,iku)    &
                  &             + r_vvl   * e3u_0(ji,jj,iku)
               ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,ikv)    &
                  &             + r_vvl   * e3v_0(ji,jj,ikv)
               puu(ji,jj,iku,Kaa) = puu(ji,jj,iku,Kaa) + rDt * 0.5*( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) * uu_b(ji,jj,Kaa) / ze3ua
               pvv(ji,jj,ikv,Kaa) = pvv(ji,jj,ikv,Kaa) + rDt * 0.5*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) ) * vv_b(ji,jj,Kaa) / ze3va
            END DO   ;   END DO
         END IF
      ENDIF
      !
      !              !==  Vertical diffusion on u  ==!
      !
      !                    !* Matrix construction
      zdt = rDt * 0.5
      IF( ln_zad_Aimp ) THEN   !!
         SELECT CASE( nldf_dyn )
         CASE( np_lap_i )           ! rotated lateral mixing: add its vertical mixing (akzu)
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,jk)    &
                  &             + r_vvl   * e3u_0(ji,jj,jk)   ! after scale factor at U-point
               zzwi = - zdt * ( avm(ji+1,jj,jk  ) + avm(ji,jj,jk  ) + akzu(ji,jj,jk  ) )   &
                  &         / ( ze3ua * e3uw_0(ji,jj,jk  ) ) * wumask(ji,jj,jk  )
               zzws = - zdt * ( avm(ji+1,jj,jk+1) + avm(ji,jj,jk+1) + akzu(ji,jj,jk+1) )   &
                  &         / ( ze3ua * e3uw_0(ji,jj,jk+1) ) * wumask(ji,jj,jk+1)
               zWui = ( wi(ji,jj,jk  ) + wi(ji+1,jj,jk  ) ) / ze3ua
               zWus = ( wi(ji,jj,jk+1) + wi(ji+1,jj,jk+1) ) / ze3ua
               zwi(ji,jj,jk) = zzwi + zdt * MIN( zWui, 0._wp ) 
               zws(ji,jj,jk) = zzws - zdt * MAX( zWus, 0._wp )
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws + zdt * ( MAX( zWui, 0._wp ) - MIN( zWus, 0._wp ) )
            END DO   ;   END DO   ;   END DO
         CASE DEFAULT               ! iso-level lateral mixing
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,jk)    &    ! after scale factor at U-point
                  &             + r_vvl   * e3u_0(ji,jj,jk)
               zzwi = - zdt * ( avm(ji+1,jj,jk  ) + avm(ji,jj,jk  ) )   &
                  &         / ( ze3ua * e3uw_0(ji,jj,jk  ) ) * wumask(ji,jj,jk  )
               zzws = - zdt * ( avm(ji+1,jj,jk+1) + avm(ji,jj,jk+1) )   &
                  &         / ( ze3ua * e3uw_0(ji,jj,jk+1) ) * wumask(ji,jj,jk+1)
               zWui = ( wi(ji,jj,jk  ) + wi(ji+1,jj,jk  ) ) / ze3ua
               zWus = ( wi(ji,jj,jk+1) + wi(ji+1,jj,jk+1) ) / ze3ua
               zwi(ji,jj,jk) = zzwi + zdt * MIN( zWui, 0._wp )
               zws(ji,jj,jk) = zzws - zdt * MAX( zWus, 0._wp )
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws + zdt * ( MAX( zWui, 0._wp ) - MIN( zWus, 0._wp ) )
            END DO   ;   END DO   ;   END DO
         END SELECT
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)     !* Surface boundary conditions
            zwi(ji,jj,1) = 0._wp
            ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,1)    &
               &             + r_vvl   * e3u_0(ji,jj,1)
            zzws = - zdt * ( avm(ji+1,jj,2) + avm(ji  ,jj,2) )   &
               &         / ( ze3ua * e3uw_0(ji,jj,2) ) * wumask(ji,jj,2)
            zWus = ( wi(ji  ,jj,2) +  wi(ji+1,jj,2) ) / ze3ua
            zws(ji,jj,1 ) = zzws - zdt * MAX( zWus, 0._wp )
            zwd(ji,jj,1 ) = 1._wp - zzws - zdt * ( MIN( zWus, 0._wp ) )
         END DO   ;   END DO
      ELSE
         SELECT CASE( nldf_dyn )
         CASE( np_lap_i )           ! rotated lateral mixing: add its vertical mixing (akzu)
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,jk)    &
                  &             + r_vvl   * e3u_0(ji,jj,jk)   ! after scale factor at U-point
               zzwi = - zdt * ( avm(ji+1,jj,jk  ) + avm(ji,jj,jk  ) + akzu(ji,jj,jk  ) )   &
                  &         / ( ze3ua * e3uw_0(ji,jj,jk  ) ) * wumask(ji,jj,jk  )
               zzws = - zdt * ( avm(ji+1,jj,jk+1) + avm(ji,jj,jk+1) + akzu(ji,jj,jk+1) )   &
                  &         / ( ze3ua * e3uw_0(ji,jj,jk+1) ) * wumask(ji,jj,jk+1)
               zwi(ji,jj,jk) = zzwi
               zws(ji,jj,jk) = zzws
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO   ;   END DO   ;   END DO
         CASE DEFAULT               ! iso-level lateral mixing
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,jk)    &
                  &             + r_vvl   * e3u_0(ji,jj,jk)   ! after scale factor at U-point
               zzwi = - zdt * ( avm(ji+1,jj,jk  ) + avm(ji,jj,jk  ) )    &
                  &         / ( ze3ua * e3uw_0(ji,jj,jk  ) ) * wumask(ji,jj,jk  )
               zzws = - zdt * ( avm(ji+1,jj,jk+1) + avm(ji,jj,jk+1) )    &
                  &         / ( ze3ua * e3uw_0(ji,jj,jk+1) ) * wumask(ji,jj,jk+1)
               zwi(ji,jj,jk) = zzwi
               zws(ji,jj,jk) = zzws
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO   ;   END DO   ;   END DO
         END SELECT
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)     !* Surface boundary conditions
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO   ;   END DO
      ENDIF
      !
      !
      !              !==  Apply semi-implicit bottom friction  ==!
      !
      !     Only needed for semi-implicit bottom friction setup. The explicit
      !     bottom friction has been included in "u(v)a" which act as the R.H.S
      !     column vector of the tri-diagonal matrix equation
      !
      IF ( ln_drgimp ) THEN      ! implicit bottom friction
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            iku = mbku(ji,jj)       ! ocean bottom level at u- and v-points
            ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,iku)    &
               &             + r_vvl   * e3u_0(ji,jj,iku)   ! after scale factor at T-point
            zwd(ji,jj,iku) = zwd(ji,jj,iku) - rDt * 0.5*( rCdU_bot(ji+1,jj)+rCdU_bot(ji,jj) ) / ze3ua
         END DO   ;   END DO
         IF ( ln_isfcav.OR.ln_drgice_imp ) THEN   ! top friction (always implicit)
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               !!gm   top Cd is masked (=0 outside cavities) no need of test on mik>=2  ==>> it has been suppressed
               iku = miku(ji,jj)       ! ocean top level at u- and v-points 
               ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,iku)    &
                  &             + r_vvl   * e3u_0(ji,jj,iku)   ! after scale factor at T-point
               zwd(ji,jj,iku) = zwd(ji,jj,iku) - rDt * 0.5*( rCdU_top(ji+1,jj)+rCdU_top(ji,jj) ) / ze3ua
            END DO   ;   END DO
         END IF
      ENDIF
      !
      ! Matrix inversion starting from the first level
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and a lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (the after velocity) is in puu(:,:,:,Kaa)
      !-----------------------------------------------------------------------
      !
      DO jk =  2,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)   !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
         zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
      END DO   ;   END DO   ;   END DO
      !
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)             !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==!
         ze3ua =  ( 1._wp - r_vvl ) * e3u_0(ji,jj,1)    &
            &             + r_vvl   * e3u_0(ji,jj,1) 
         puu(ji,jj,1,Kaa) = puu(ji,jj,1,Kaa) + rDt * 0.5_wp * ( utau_b(ji,jj) + utau(ji,jj) )   &
            &                                      / ( ze3ua * rho0 ) * umask(ji,jj,1) 
      END DO   ;   END DO
      DO jk =  2,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
         puu(ji,jj,jk,Kaa) = puu(ji,jj,jk,Kaa) - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * puu(ji,jj,jk-1,Kaa)
      END DO   ;   END DO   ;   END DO
      !
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)             !==  thrid recurrence : SOLk = ( Lk - Uk * Ek+1 ) / Dk  ==!
         puu(ji,jj,jpkm1,Kaa) = puu(ji,jj,jpkm1,Kaa) / zwd(ji,jj,jpkm1)
      END DO   ;   END DO
      DO jk =  jpk-2,  1,  -1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
         puu(ji,jj,jk,Kaa) = ( puu(ji,jj,jk,Kaa) - zws(ji,jj,jk) * puu(ji,jj,jk+1,Kaa) ) / zwd(ji,jj,jk)
      END DO   ;   END DO   ;   END DO
      !
      !              !==  Vertical diffusion on v  ==!
      !
      !                       !* Matrix construction
      zdt = rDt * 0.5
      IF( ln_zad_Aimp ) THEN   !!
         SELECT CASE( nldf_dyn )
         CASE( np_lap_i )           ! rotated lateral mixing: add its vertical mixing (akzv)
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,jk)    &
                  &             + r_vvl   * e3v_0(ji,jj,jk)   ! after scale factor at V-point
               zzwi = - zdt * ( avm(ji,jj+1,jk  ) + avm(ji,jj,jk  ) + akzv(ji,jj,jk  ) )   &
                  &         / ( ze3va * e3vw_0(ji,jj,jk  ) ) * wvmask(ji,jj,jk  )
               zzws = - zdt * ( avm(ji,jj+1,jk+1) + avm(ji,jj,jk+1) + akzv(ji,jj,jk+1) )   &
                  &         / ( ze3va * e3vw_0(ji,jj,jk+1) ) * wvmask(ji,jj,jk+1)
               zWvi = ( wi(ji,jj,jk  ) + wi(ji,jj+1,jk  ) ) / ze3va
               zWvs = ( wi(ji,jj,jk+1) + wi(ji,jj+1,jk+1) ) / ze3va
               zwi(ji,jj,jk) = zzwi + zdt * MIN( zWvi, 0._wp )
               zws(ji,jj,jk) = zzws - zdt * MAX( zWvs, 0._wp )
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws - zdt * ( - MAX( zWvi, 0._wp ) + MIN( zWvs, 0._wp ) )
            END DO   ;   END DO   ;   END DO
         CASE DEFAULT               ! iso-level lateral mixing
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,jk)    &
                  &             + r_vvl   * e3v_0(ji,jj,jk)   ! after scale factor at V-point
               zzwi = - zdt * ( avm(ji,jj+1,jk  ) + avm(ji,jj,jk  ) )    &
                  &         / ( ze3va * e3vw_0(ji,jj,jk  ) ) * wvmask(ji,jj,jk  )
               zzws = - zdt * ( avm(ji,jj+1,jk+1) + avm(ji,jj,jk+1) )    &
                  &         / ( ze3va * e3vw_0(ji,jj,jk+1) ) * wvmask(ji,jj,jk+1)
               zWvi = ( wi(ji,jj,jk  ) + wi(ji,jj+1,jk  ) ) / ze3va
               zWvs = ( wi(ji,jj,jk+1) + wi(ji,jj+1,jk+1) ) / ze3va
               zwi(ji,jj,jk) = zzwi  + zdt * MIN( zWvi, 0._wp )
               zws(ji,jj,jk) = zzws  - zdt * MAX( zWvs, 0._wp )
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws - zdt * ( - MAX( zWvi, 0._wp ) + MIN( zWvs, 0._wp ) )
            END DO   ;   END DO   ;   END DO
         END SELECT
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)   !* Surface boundary conditions
            zwi(ji,jj,1) = 0._wp
            ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,1)    &
               &             + r_vvl   * e3v_0(ji,jj,1)
            zzws = - zdt * ( avm(ji,jj+1,2) + avm(ji,jj,2) )    &
               &         / ( ze3va * e3vw_0(ji,jj,2) ) * wvmask(ji,jj,2)
            zWvs = ( wi(ji,jj  ,2) +  wi(ji,jj+1,2) ) / ze3va
            zws(ji,jj,1 ) = zzws - zdt * MAX( zWvs, 0._wp )
            zwd(ji,jj,1 ) = 1._wp - zzws - zdt * ( MIN( zWvs, 0._wp ) )
         END DO   ;   END DO
      ELSE
         SELECT CASE( nldf_dyn )
         CASE( np_lap_i )           ! rotated lateral mixing: add its vertical mixing (akzu)
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,jk)    &
                  &             + r_vvl   * e3v_0(ji,jj,jk)   ! after scale factor at V-point
               zzwi = - zdt * ( avm(ji,jj+1,jk  ) + avm(ji,jj,jk  ) + akzv(ji,jj,jk  ) )   &
                  &         / ( ze3va * e3vw_0(ji,jj,jk  ) ) * wvmask(ji,jj,jk  )
               zzws = - zdt * ( avm(ji,jj+1,jk+1) + avm(ji,jj,jk+1) + akzv(ji,jj,jk+1) )   &
                  &         / ( ze3va * e3vw_0(ji,jj,jk+1) ) * wvmask(ji,jj,jk+1)
               zwi(ji,jj,jk) = zzwi
               zws(ji,jj,jk) = zzws
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO   ;   END DO   ;   END DO
         CASE DEFAULT               ! iso-level lateral mixing
            DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
               ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,jk)    &
                  &             + r_vvl   * e3v_0(ji,jj,jk)   ! after scale factor at V-point
               zzwi = - zdt * ( avm(ji,jj+1,jk  ) + avm(ji,jj,jk  ) )    &
                  &         / ( ze3va * e3vw_0(ji,jj,jk  ) ) * wvmask(ji,jj,jk  )
               zzws = - zdt * ( avm(ji,jj+1,jk+1) + avm(ji,jj,jk+1) )    &
                  &         / ( ze3va * e3vw_0(ji,jj,jk+1) ) * wvmask(ji,jj,jk+1)
               zwi(ji,jj,jk) = zzwi
               zws(ji,jj,jk) = zzws
               zwd(ji,jj,jk) = 1._wp - zzwi - zzws
            END DO   ;   END DO   ;   END DO
         END SELECT
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)        !* Surface boundary conditions
            zwi(ji,jj,1) = 0._wp
            zwd(ji,jj,1) = 1._wp - zws(ji,jj,1)
         END DO   ;   END DO
      ENDIF
      !
      !              !==  Apply semi-implicit top/bottom friction  ==!
      !
      !     Only needed for semi-implicit bottom friction setup. The explicit
      !     bottom friction has been included in "u(v)a" which act as the R.H.S
      !     column vector of the tri-diagonal matrix equation
      !
      IF( ln_drgimp ) THEN
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            ikv = mbkv(ji,jj)       ! (deepest ocean u- and v-points)
            ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,ikv)    &
               &             + r_vvl   * e3v_0(ji,jj,ikv)   ! after scale factor at T-point
            zwd(ji,jj,ikv) = zwd(ji,jj,ikv) - rDt * 0.5*( rCdU_bot(ji,jj+1)+rCdU_bot(ji,jj) ) / ze3va           
         END DO   ;   END DO
         IF ( ln_isfcav.OR.ln_drgice_imp ) THEN
            DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
               ikv = mikv(ji,jj)       ! (first wet ocean u- and v-points)
               ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,ikv)    &
                  &             + r_vvl   * e3v_0(ji,jj,ikv)   ! after scale factor at T-point
               zwd(ji,jj,ikv) = zwd(ji,jj,ikv) - rDt * 0.5*( rCdU_top(ji,jj+1)+rCdU_top(ji,jj) ) / ze3va
            END DO   ;   END DO
         ENDIF
      ENDIF

      ! Matrix inversion
      !-----------------------------------------------------------------------
      !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
      !
      !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
      !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
      !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
      !        (        ...               )( ...  ) ( ...  )
      !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
      !
      !   m is decomposed in the product of an upper and lower triangular matrix
      !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
      !   The solution (after velocity) is in 2d array va
      !-----------------------------------------------------------------------
      !
      DO jk =  2,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)   !==  First recurrence : Dk = Dk - Lk * Uk-1 / Dk-1   (increasing k)  ==
         zwd(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1) / zwd(ji,jj,jk-1)
      END DO   ;   END DO   ;   END DO
      !
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)             !==  second recurrence:    SOLk = RHSk - Lk / Dk-1  Lk-1  ==!
         ze3va =  ( 1._wp - r_vvl ) * e3v_0(ji,jj,1)    &
            &             + r_vvl   * e3v_0(ji,jj,1) 
         pvv(ji,jj,1,Kaa) = pvv(ji,jj,1,Kaa) + rDt * 0.5_wp * ( vtau_b(ji,jj) + vtau(ji,jj) )   &
            &                                      / ( ze3va * rho0 ) * vmask(ji,jj,1) 
      END DO   ;   END DO
      DO jk =  2,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
         pvv(ji,jj,jk,Kaa) = pvv(ji,jj,jk,Kaa) - zwi(ji,jj,jk) / zwd(ji,jj,jk-1) * pvv(ji,jj,jk-1,Kaa)
      END DO   ;   END DO   ;   END DO
      !
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)             !==  third recurrence : SOLk = ( Lk - Uk * SOLk+1 ) / Dk  ==!
         pvv(ji,jj,jpkm1,Kaa) = pvv(ji,jj,jpkm1,Kaa) / zwd(ji,jj,jpkm1)
      END DO   ;   END DO
      DO jk =  jpk-2,  1,  -1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
         pvv(ji,jj,jk,Kaa) = ( pvv(ji,jj,jk,Kaa) - zws(ji,jj,jk) * pvv(ji,jj,jk+1,Kaa) ) / zwd(ji,jj,jk)
      END DO   ;   END DO   ;   END DO
      !
      IF( l_trddyn )   THEN                      ! save the vertical diffusive trends for further diagnostics
         ztrdu(:,:,:) = ( puu(:,:,:,Kaa) - puu(:,:,:,Kbb) ) / rDt - ztrdu(:,:,:)
         ztrdv(:,:,:) = ( pvv(:,:,:,Kaa) - pvv(:,:,:,Kbb) ) / rDt - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_zdf, kt, Kmm )
         DEALLOCATE( ztrdu, ztrdv ) 
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Kaa), clinfo1=' zdf  - Ua: ', mask1=umask,               &
         &                                  tab3d_2=pvv(:,:,:,Kaa), clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
         !
      IF( ln_timing )   CALL timing_stop('dyn_zdf')
      !
   END SUBROUTINE dyn_zdf

   !!==============================================================================
END MODULE dynzdf

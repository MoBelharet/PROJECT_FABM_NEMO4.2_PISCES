










MODULE dynadv_cen2
   !!======================================================================
   !!                       ***  MODULE  dynadv  ***
   !! Ocean dynamics: Update the momentum trend with the flux form advection
   !!                 using a 2nd order centred scheme
   !!======================================================================
   !! History :  2.0  ! 2006-08  (G. Madec, S. Theetten)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv_cen2  : flux form momentum advection (ln_dynadv_cen2=T) using a 2nd order centred scheme  
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_adv_cen2   ! routine called by step.F90

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
   !! $Id: dynadv_cen2.F90 14834 2021-05-11 09:24:44Z hadcv $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_adv_cen2( kt, Kmm, puu, pvv, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_cen2  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time) 
      !!
      !! ** Action  :   (puu(:,:,:,Krhs),pvv(:,:,:,Krhs)) updated with the now vorticity term trend
      !!----------------------------------------------------------------------
      INTEGER                             , INTENT( in )  ::  kt           ! ocean time-step index
      INTEGER                             , INTENT( in )  ::  Kmm, Krhs    ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) ::  puu, pvv     ! ocean velocities and RHS of momentum equation
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls),jpk)  :: zfu_f, zfu
      REAL(dp), DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls),jpk)  :: zfu_t, zfu_uw
      REAL(wp), DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls),jpk)  :: zfv_f, zfv, zfw
      REAL(dp), DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls),jpk)  :: zfv_t, zfv_vw
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 .AND. lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'dyn_adv_cen2 : 2nd order flux form momentum advection'
            WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      IF( l_trddyn ) THEN           ! trends: store the input trends
         zfu_uw(:,:,:) = puu(:,:,:,Krhs)
         zfv_vw(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF
      !
      !                             !==  Horizontal advection  ==!
      !
      DO jk = 1, jpkm1                    ! horizontal transport
         DO jj = ntsj-( 1), ntej+( 1 ) ; DO ji = ntsi-( 1), ntei+( 1)
            zfu(ji,jj,jk) = 0.25_wp * e2u(ji,jj) * e3u_0(ji,jj,jk) * puu(ji,jj,jk,Kmm)
            zfv(ji,jj,jk) = 0.25_wp * e1v(ji,jj) * e3v_0(ji,jj,jk) * pvv(ji,jj,jk,Kmm)
         END DO   ;   END DO
         DO jj = ntsj-( 1), ntej+( 0 ) ; DO ji = ntsi-( 1), ntei+( 0)              ! horizontal momentum fluxes (at T- and F-point)
            zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji+1,jj  ,jk,Kmm) )
            zfv_f(ji  ,jj  ,jk) = ( zfv(ji,jj,jk) + zfv(ji+1,jj,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji  ,jj+1,jk,Kmm) )
            zfu_f(ji  ,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji+1,jj  ,jk,Kmm) )
            zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji  ,jj+1,jk,Kmm) )
         END DO   ;   END DO
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)              ! divergence of horizontal momentum fluxes
            puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - (  zfu_t(ji+1,jj,jk) - zfu_t(ji,jj  ,jk)    &
               &                           + zfv_f(ji  ,jj,jk) - zfv_f(ji,jj-1,jk)  ) * r1_e1e2u(ji,jj)   &
               &                           / e3u_0(ji,jj,jk)
            pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - (  zfu_f(ji,jj  ,jk) - zfu_f(ji-1,jj,jk)    &
               &                           + zfv_t(ji,jj+1,jk) - zfv_t(ji  ,jj,jk)  ) * r1_e1e2v(ji,jj)   &
               &                           / e3v_0(ji,jj,jk)
         END DO   ;   END DO
      END DO
      !
      IF( l_trddyn ) THEN           ! trends: send trend to trddyn for diagnostic
         zfu_uw(:,:,:) = puu(:,:,:,Krhs) - zfu_uw(:,:,:)
         zfv_vw(:,:,:) = pvv(:,:,:,Krhs) - zfv_vw(:,:,:)
         CALL trd_dyn( zfu_uw, zfv_vw, jpdyn_keg, kt, Kmm )
         zfu_t(:,:,:) = puu(:,:,:,Krhs)
         zfv_t(:,:,:) = pvv(:,:,:,Krhs)
      ENDIF
      !
      !                             !==  Vertical advection  ==!
      !
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)                 ! surface/bottom advective fluxes set to zero
         zfu_uw(ji,jj,jpk) = 0._wp   ;   zfv_vw(ji,jj,jpk) = 0._wp
         zfu_uw(ji,jj, 1 ) = 0._wp   ;   zfv_vw(ji,jj, 1 ) = 0._wp
      END DO   ;   END DO
      IF( ln_linssh ) THEN                ! linear free surface: advection through the surface
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zfu_uw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * ww(ji,jj,1) + e1e2t(ji+1,jj) * ww(ji+1,jj,1) ) * puu(ji,jj,1,Kmm)
            zfv_vw(ji,jj,1) = 0.5_wp * ( e1e2t(ji,jj) * ww(ji,jj,1) + e1e2t(ji,jj+1) * ww(ji,jj+1,1) ) * pvv(ji,jj,1,Kmm)
         END DO   ;   END DO
      ENDIF
      DO jk = 2, jpkm1                    ! interior advective fluxes
         DO jj = ntsj-( 0), ntej+( 1 ) ; DO ji = ntsi-( 0), ntei+( 1)                  ! 1/4 * Vertical transport
            zfw(ji,jj,jk) = 0.25_wp * e1e2t(ji,jj) * ww(ji,jj,jk)
         END DO   ;   END DO
         DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
            zfu_uw(ji,jj,jk) = ( zfw(ji,jj,jk) + zfw(ji+1,jj  ,jk) ) * ( puu(ji,jj,jk,Kmm) + puu(ji,jj,jk-1,Kmm) )
            zfv_vw(ji,jj,jk) = ( zfw(ji,jj,jk) + zfw(ji  ,jj+1,jk) ) * ( pvv(ji,jj,jk,Kmm) + pvv(ji,jj,jk-1,Kmm) )
         END DO   ;   END DO
      END DO
      DO jk =  1,  jpkm1  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)       ! divergence of vertical momentum flux divergence
         puu(ji,jj,jk,Krhs) = puu(ji,jj,jk,Krhs) - ( zfu_uw(ji,jj,jk) - zfu_uw(ji,jj,jk+1) ) * r1_e1e2u(ji,jj)   &
            &                                      / e3u_0(ji,jj,jk)
         pvv(ji,jj,jk,Krhs) = pvv(ji,jj,jk,Krhs) - ( zfv_vw(ji,jj,jk) - zfv_vw(ji,jj,jk+1) ) * r1_e1e2v(ji,jj)   &
            &                                      / e3v_0(ji,jj,jk)
      END DO   ;   END DO   ;   END DO
      !
      IF( l_trddyn ) THEN                 ! trends: send trend to trddyn for diagnostic
         zfu_t(:,:,:) = puu(:,:,:,Krhs) - zfu_t(:,:,:)
         zfv_t(:,:,:) = pvv(:,:,:,Krhs) - zfv_t(:,:,:)
         CALL trd_dyn( zfu_t, zfv_t, jpdyn_zad, kt, Kmm )
      ENDIF
      !                                   ! Control print
      IF(sn_cfctl%l_prtctl)   CALL prt_ctl( tab3d_1=puu(:,:,:,Krhs), clinfo1=' cen2 adv - Ua: ', mask1=umask,   &
         &                                  tab3d_2=pvv(:,:,:,Krhs), clinfo2=           ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
   END SUBROUTINE dyn_adv_cen2

   !!==============================================================================
END MODULE dynadv_cen2

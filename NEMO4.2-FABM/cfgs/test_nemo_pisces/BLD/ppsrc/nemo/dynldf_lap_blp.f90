










MODULE dynldf_lap_blp
   !!======================================================================
   !!                   ***  MODULE  dynldf_lap_blp  ***
   !! Ocean dynamics:  lateral viscosity trend (laplacian and bilaplacian)
   !!======================================================================
   !! History : 3.7  ! 2014-01  (G. Madec, S. Masson)  Original code, re-entrant laplacian
   !!           4.0  ! 2020-04  (A. Nasser, G. Madec)  Add symmetric mixing tensor 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_lap   : update the momentum trend with the lateral viscosity using an iso-level   laplacian operator
   !!   dyn_ldf_blp   : update the momentum trend with the lateral viscosity using an iso-level bilaplacian operator
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE domutl, ONLY : is_tile
   USE ldfdyn         ! lateral diffusion: eddy viscosity coef.
   USE ldfslp         ! iso-neutral slopes 
   USE zdf_oce        ! ocean vertical physics
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_ldf_lap  ! called by dynldf.F90
   PUBLIC dyn_ldf_blp  ! called by dynldf.F90

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
   !! $Id: dynldf_lap_blp.F90 15033 2021-06-21 10:24:45Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_lap( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs, kpass )
      !!
      INTEGER                   , INTENT(in   ) ::   kt               ! ocean time-step index
      INTEGER                   , INTENT(in   ) ::   Kbb, Kmm         ! ocean time level indices
      INTEGER                   , INTENT(in   ) ::   kpass            ! =1/2 first or second passage
      REAL(dp), DIMENSION(:,:,:), INTENT(in   ) ::   pu, pv           ! before velocity  [m/s]
      REAL(dp), DIMENSION(:,:,:), INTENT(inout) ::   pu_rhs, pv_rhs   ! velocity trend   [m/s2]
      !!
      CALL dyn_ldf_lap_t( kt, Kbb, Kmm, pu, pv, is_tile(pu), pu_rhs, pv_rhs, is_tile(pu_rhs), kpass )

   END SUBROUTINE dyn_ldf_lap


   SUBROUTINE dyn_ldf_lap_t( kt, Kbb, Kmm, pu, pv, ktuv, pu_rhs, pv_rhs, ktuv_rhs, kpass )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!                       
      !! ** Purpose :   Compute the before horizontal momentum diffusive 
      !!      trend and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The Laplacian operator apply on horizontal velocity is 
      !!      writen as :   grad_h( ahmt div_h(U )) - curl_h( ahmf curl_z(U) ) 
      !!
      !! ** Action : - pu_rhs, pv_rhs increased by the harmonic operator applied on pu, pv.
      !!
      !! Reference : S.Griffies, R.Hallberg 2000 Mon.Wea.Rev., DOI:/ 
      !!----------------------------------------------------------------------
      INTEGER                                 , INTENT(in   ) ::   kt               ! ocean time-step index
      INTEGER                                 , INTENT(in   ) ::   Kbb, Kmm         ! ocean time level indices
      INTEGER                                 , INTENT(in   ) ::   kpass            ! =1/2 first or second passage
      INTEGER                                 , INTENT(in   ) ::   ktuv, ktuv_rhs
      REAL(dp), DIMENSION((ntsi-nn_hls-1)*ktuv+1:,(ntsj-nn_hls-1)*ktuv+1:    ,:), INTENT(in   ) ::   pu, pv           ! before velocity  [m/s]
      REAL(dp), DIMENSION((ntsi-nn_hls-1)*ktuv_rhs+1:,(ntsj-nn_hls-1)*ktuv_rhs+1:,:), INTENT(inout) ::   pu_rhs, pv_rhs   ! velocity trend   [m/s2]
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   iij
      REAL(wp) ::   zsign        ! local scalars
      REAL(wp) ::   zua, zva     ! local scalars
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zcur, zdiv
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) ::   zten, zshe   ! tension (diagonal) and shearing (anti-diagonal) terms
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 .AND. lwp ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacian) operator, pass=', kpass
            WRITE(numout,*) '~~~~~~~ '
         ENDIF
      ENDIF
      !
      ! Define pu_rhs/pv_rhs halo points for multi-point haloes in bilaplacian case
      IF( nldf_dyn == np_blp .AND. kpass == 1 ) THEN ; iij = nn_hls
      ELSE                                           ; iij = 1
      ENDIF
      !
      IF( kpass == 1 ) THEN   ;   zsign =  1._wp      ! bilaplacian operator require a minus sign
      ELSE                    ;   zsign = -1._wp      !  (eddy viscosity coef. >0)
      ENDIF
      !
      SELECT CASE( nn_dynldf_typ )  
      !              
      CASE ( np_typ_rot )       !==  Vorticity-Divergence operator  ==!
         !
         ALLOCATE( zcur(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls)) , zdiv(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls)) )
         !
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !
            DO jj = ntsj-( iij-1), ntej+( iij ) ; DO ji = ntsi-( iij-1), ntei+( iij)
               !                                      ! ahm * e3 * curl  (warning: computed for ji-1,jj-1)
               zcur(ji-1,jj-1) = ahmf(ji-1,jj-1,jk) * e3f_0(ji-1,jj-1,jk) * r1_e1e2f(ji-1,jj-1)       &   ! ahmf already * by fmask
                  &     * (  e2v(ji  ,jj-1) * pv(ji  ,jj-1,jk) - e2v(ji-1,jj-1) * pv(ji-1,jj-1,jk)  &
                  &        - e1u(ji-1,jj  ) * pu(ji-1,jj  ,jk) + e1u(ji-1,jj-1) * pu(ji-1,jj-1,jk)  )
               !                                      ! ahm * div        (warning: computed for ji,jj)
               zdiv(ji,jj)     = ahmt(ji,jj,jk) * r1_e1e2t(ji,jj) / e3t_0(ji,jj,jk)               &   ! ahmt already * by tmask
                  &     * (  e2u(ji,jj)*e3u_0(ji,jj,jk) * pu(ji,jj,jk) - e2u(ji-1,jj)*e3u_0(ji-1,jj,jk) * pu(ji-1,jj,jk)  &
                  &        + e1v(ji,jj)*e3v_0(ji,jj,jk) * pv(ji,jj,jk) - e1v(ji,jj-1)*e3v_0(ji,jj-1,jk) * pv(ji,jj-1,jk)  )
            END DO   ;   END DO
            !
            DO jj = ntsj-( iij-1), ntej+( iij-1 ) ; DO ji = ntsi-( iij-1), ntei+( iij-1)   ! - curl( curl) + grad( div )
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * umask(ji,jj,jk) * (    &    ! * by umask is mandatory for dyn_ldf_blp use
                  &              - ( zcur(ji  ,jj) - zcur(ji,jj-1) ) * r1_e2u(ji,jj) / e3u_0(ji,jj,jk)   &
                  &              + ( zdiv(ji+1,jj) - zdiv(ji,jj  ) ) * r1_e1u(ji,jj)                      )
               !
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * vmask(ji,jj,jk) * (    &    ! * by vmask is mandatory for dyn_ldf_blp use
                  &                ( zcur(ji,jj  ) - zcur(ji-1,jj) ) * r1_e1v(ji,jj) / e3v_0(ji,jj,jk)   &
                  &              + ( zdiv(ji,jj+1) - zdiv(ji  ,jj) ) * r1_e2v(ji,jj)                      )
            END DO   ;   END DO
            !
         END DO                                           !   End of slab
         !
         DEALLOCATE( zcur , zdiv )
         !
      CASE ( np_typ_sym )       !==  Symmetric operator  ==!
         !
         ALLOCATE( zten(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls)) , zshe(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls)) )
         !
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !
            DO jj = ntsj-( iij-1), ntej+( iij ) ; DO ji = ntsi-( iij-1), ntei+( iij)
               !                                      ! shearing stress component (F-point)   NB : ahmf has already been multiplied by fmask
               zshe(ji-1,jj-1) = ahmf(ji-1,jj-1,jk)                                                              &
                  &     * (    e1f(ji-1,jj-1)    * r1_e2f(ji-1,jj-1)                                             &
                  &         * ( pu(ji-1,jj  ,jk) * r1_e1u(ji-1,jj  )  - pu(ji-1,jj-1,jk) * r1_e1u(ji-1,jj-1) )   &
                  &         +  e2f(ji-1,jj-1)    * r1_e1f(ji-1,jj-1)                                             &
                  &         * ( pv(ji  ,jj-1,jk) * r1_e2v(ji  ,jj-1)  - pv(ji-1,jj-1,jk) * r1_e2v(ji-1,jj-1) )   ) 
               !                                      ! tension stress component (T-point)   NB : ahmt has already been multiplied by tmask
               zten(ji,jj)    = ahmt(ji,jj,jk)                                                       &
                  &     * (    e2t(ji,jj)    * r1_e1t(ji,jj)                                         &
                  &         * ( pu(ji,jj,jk) * r1_e2u(ji,jj)  - pu(ji-1,jj,jk) * r1_e2u(ji-1,jj) )   &
                  &         -  e1t(ji,jj)    * r1_e2t(ji,jj)                                         &
                  &         * ( pv(ji,jj,jk) * r1_e1v(ji,jj)  - pv(ji,jj-1,jk) * r1_e1v(ji,jj-1) )   )   
            END DO   ;   END DO
            !
            DO jj = ntsj-( iij-1), ntej+( iij-1 ) ; DO ji = ntsi-( iij-1), ntei+( iij-1)
               pu_rhs(ji,jj,jk) = pu_rhs(ji,jj,jk) + zsign * r1_e1e2u(ji,jj) / e3u_0(ji,jj,jk)                               &
                  &    * (   (   zten(ji+1,jj  ) * e2t(ji+1,jj  )*e2t(ji+1,jj  ) * e3t_0(ji+1,jj  ,jk)                       &
                  &            - zten(ji  ,jj  ) * e2t(ji  ,jj  )*e2t(ji  ,jj  ) * e3t_0(ji  ,jj  ,jk) ) * r1_e2u(ji,jj)     &                                                    
                  &        + (   zshe(ji  ,jj  ) * e1f(ji  ,jj  )*e1f(ji  ,jj  ) * e3f_0(ji  ,jj  ,jk)                           &
                  &            - zshe(ji  ,jj-1) * e1f(ji  ,jj-1)*e1f(ji  ,jj-1) * e3f_0(ji  ,jj-1,jk)     ) * r1_e1u(ji,jj) )   
               !
               pv_rhs(ji,jj,jk) = pv_rhs(ji,jj,jk) + zsign * r1_e1e2v(ji,jj) / e3v_0(ji,jj,jk)                               &
                  &    * (   (   zshe(ji  ,jj  ) * e2f(ji  ,jj  )*e2f(ji  ,jj  ) * e3f_0(ji  ,jj  ,jk)                           &
                  &            - zshe(ji-1,jj  ) * e2f(ji-1,jj  )*e2f(ji-1,jj  ) * e3f_0(ji-1,jj  ,jk)     ) * r1_e2v(ji,jj)     &
                  &        - (   zten(ji  ,jj+1) * e1t(ji  ,jj+1)*e1t(ji  ,jj+1) * e3t_0(ji  ,jj+1,jk)                       &
                  &            - zten(ji  ,jj  ) * e1t(ji  ,jj  )*e1t(ji  ,jj  ) * e3t_0(ji  ,jj  ,jk) ) * r1_e1v(ji,jj) )
               !
            END DO   ;   END DO
            !
         END DO
         !
         DEALLOCATE( zten , zshe )
         !
      END SELECT
      !
   END SUBROUTINE dyn_ldf_lap_t


   SUBROUTINE dyn_ldf_blp( kt, Kbb, Kmm, pu, pv, pu_rhs, pv_rhs )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dyn_ldf_blp  ***
      !!                    
      !! ** Purpose :   Compute the before lateral momentum viscous trend 
      !!              and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The lateral viscous trends is provided by a bilaplacian
      !!      operator applied to before field (forward in time).
      !!      It is computed by two successive calls to dyn_ldf_lap routine
      !!
      !! ** Action :   pt(:,:,:,:,Krhs)   updated with the before rotated bilaplacian diffusion
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   Kbb, Kmm   ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk), INTENT(in   ) ::   pu, pv     ! before velocity fields
      REAL(dp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pu_rhs, pv_rhs   ! momentum trend
      !
      REAL(dp), DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls),jpk) ::   zulap, zvlap   ! laplacian at u- and v-point
      !!----------------------------------------------------------------------
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 )  THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_ldf_blp : bilaplacian operator momentum '
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~'
         ENDIF
      ENDIF
      !
      zulap(:,:,:) = 0._wp
      zvlap(:,:,:) = 0._wp
      !
      CALL dyn_ldf_lap( kt, Kbb, Kmm, pu, pv, zulap, zvlap, 1 )   ! rotated laplacian applied to pt (output in zlap,Kbb)
      !
      IF (nn_hls==1) CALL lbc_lnk( 'dynldf_lap_blp', zulap, 'U', -1.0_dp, zvlap, 'V', -1.0_dp )             ! Lateral boundary conditions
      !
      CALL dyn_ldf_lap( kt, Kbb, Kmm, zulap, zvlap, pu_rhs, pv_rhs, 2 )   ! rotated laplacian applied to zlap (output in pt(:,:,:,:,Krhs))
      !
   END SUBROUTINE dyn_ldf_blp

   !!======================================================================
END MODULE dynldf_lap_blp

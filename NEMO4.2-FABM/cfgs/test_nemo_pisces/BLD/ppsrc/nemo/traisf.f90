










MODULE traisf
   !!==============================================================================
   !!                       ***  MODULE  traisf  ***
   !! Ocean active tracers:  ice shelf boundary condition
   !!==============================================================================
   !! History :    4.0  !  2019-09  (P. Mathiot) original file
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_isf       : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE isf_oce                                     ! Ice shelf variables
   USE par_oce , ONLY : nijtile, ntile, ntsi, ntei, ntsj, ntej
   USE dom_oce                                     ! ocean space domain variables
   USE isfutils, ONLY : debug                      ! debug option
   USE timing  , ONLY : timing_start, timing_stop  ! Timing
   USE in_out_manager                              ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_isf   ! routine called by step.F90

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
   !! $Id: trasbc.F90 10499 2019-01-10 15:12:24Z deazer $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_isf ( kt, Kmm, pts, Krhs )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_isf  ***
      !!
      !! ** Purpose :  Compute the temperature trend due to the ice shelf melting (qhoce + qhc)
      !!
      !! ** Action  : - update pts(:,:,:,:,Krhs) for cav, par and cpl case
      !!----------------------------------------------------------------------
      INTEGER                                  , INTENT(in   ) :: kt        ! ocean time step
      INTEGER                                  , INTENT(in   ) :: Kmm, Krhs ! ocean time level indices
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpts,jpt), INTENT(inout) :: pts       ! active tracers and RHS of tracer equation
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('tra_isf')
      !
      IF( .NOT. l_istiled .OR. ntile == 1 )  THEN                       ! Do only on the first tile
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'tra_isf : Ice shelf heat fluxes'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
         ENDIF
      ENDIF
      !
      ! cavity case
      IF ( ln_isfcav_mlt ) CALL tra_isf_mlt(misfkt_cav, misfkb_cav, rhisf_tbl_cav, rfrac_tbl_cav, risf_cav_tsc, risf_cav_tsc_b, pts(:,:,:,:,Krhs))
      !
      ! parametrisation case
      IF ( ln_isfpar_mlt ) CALL tra_isf_mlt(misfkt_par, misfkb_par, rhisf_tbl_par, rfrac_tbl_par, risf_par_tsc, risf_par_tsc_b, pts(:,:,:,:,Krhs))
      !
      ! ice sheet coupling case
      IF ( ln_isfcpl ) THEN
         !
         ! Dynamical stability at start up after change in under ice shelf cavity geometry is achieve by correcting the divergence.
         ! This is achieved by applying a volume flux in order to keep the horizontal divergence after remapping
         ! the same as at the end of the latest time step. So correction need to be apply at nit000 (euler time step) and
         ! half of it at nit000+1 (leap frog time step).
         ! in accordance to this, the heat content flux due to injected water need to be added in the temperature and salt trend
         ! at time step nit000 and nit000+1
         IF ( kt == nit000  ) CALL tra_isf_cpl(Kmm, risfcpl_tsc       , pts(:,:,:,:,Krhs))
         IF ( kt == nit000+1) CALL tra_isf_cpl(Kmm, risfcpl_tsc*0.5_wp, pts(:,:,:,:,Krhs))
         !
         ! ensure 0 trend due to unconservation of the ice shelf coupling
         IF ( ln_isfcpl_cons ) CALL tra_isf_cpl(Kmm, risfcpl_cons_tsc, pts(:,:,:,:,Krhs))
         !
      END IF
      !
      IF ( ln_isfdebug ) THEN
         IF( .NOT. l_istiled .OR. ntile == nijtile ) THEN                       ! Do only for the full domain
            CALL debug('tra_isf: pts(:,:,:,:,Krhs) T', pts(:,:,:,1,Krhs))
            CALL debug('tra_isf: pts(:,:,:,:,Krhs) S', pts(:,:,:,2,Krhs))
         ENDIF
      END IF
      !
      IF( ln_timing )   CALL timing_stop('tra_isf')
      !
   END SUBROUTINE tra_isf
   !
   SUBROUTINE tra_isf_mlt(ktop, kbot, phtbl, pfrac, ptsc, ptsc_b, pts)
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_isf_mlt  ***
      !!
      !! *** Purpose :  Compute the temperature trend due to the ice shelf melting (qhoce + qhc) for cav or par case
      !!
      !! *** Action :: Update pts(:,:,:,:,Krhs) with the surface boundary condition trend
      !!
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) :: pts
      !!----------------------------------------------------------------------
      INTEGER , DIMENSION(jpi,jpj)     , INTENT(in   ) :: ktop , kbot
      REAL(wp), DIMENSION(jpi,jpj)     , INTENT(in   ) :: phtbl, pfrac
      REAL(wp), DIMENSION(jpi,jpj,jpts), INTENT(in   ) :: ptsc , ptsc_b
      !!----------------------------------------------------------------------
      INTEGER                      :: ji,jj,jk  ! loop index
      INTEGER                      :: ikt, ikb  ! top and bottom level of the tbl
      REAL(wp), DIMENSION(ntsi-(nn_hls):ntei+(nn_hls),ntsj-(nn_hls):ntej+(nn_hls))     :: ztc       ! total ice shelf tracer trend
      !!----------------------------------------------------------------------
      !
      ! compute 2d total trend due to isf
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
         ztc(ji,jj) = 0.5_wp * ( ptsc(ji,jj,jp_tem) + ptsc_b(ji,jj,jp_tem) ) / phtbl(ji,jj)
      END DO   ;   END DO
      !
      ! update pts(:,:,:,:,Krhs)
      DO jj = ntsj-( 0), ntej+( 0 ) ; DO ji = ntsi-( 0), ntei+( 0)
         !
         ikt = ktop(ji,jj)
         ikb = kbot(ji,jj)
         !
         ! level fully include in the ice shelf boundary layer
         DO jk = ikt, ikb - 1
            pts(ji,jj,jk,jp_tem) = pts(ji,jj,jk,jp_tem) + ztc(ji,jj)
         END DO
         !
         ! level partially include in ice shelf boundary layer
         pts(ji,jj,ikb,jp_tem) = pts(ji,jj,ikb,jp_tem) + ztc(ji,jj) * pfrac(ji,jj)
         !
      END DO   ;   END DO
      !
   END SUBROUTINE tra_isf_mlt
   !
   SUBROUTINE tra_isf_cpl( Kmm, ptsc, ptsa )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_isf_cpl  ***
      !!
      !! *** Action :: Update pts(:,:,:,:,Krhs) with the ice shelf coupling trend
      !!
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(inout) :: ptsa
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) :: Kmm   ! ocean time level index
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpts), INTENT(in   ) :: ptsc
      !!----------------------------------------------------------------------
      INTEGER :: ji, jj, jk
      !!----------------------------------------------------------------------
      !
      DO jk =  1,  jpk  ; DO jj = ntsj-(  0), ntej+(  0) ; DO ji = ntsi-( 0), ntei+(  0)
         ptsa(ji,jj,jk,jp_tem) = ptsa(ji,jj,jk,jp_tem) + ptsc(ji,jj,jk,jp_tem) * r1_e1e2t(ji,jj) / e3t_0(ji,jj,jk)
         ptsa(ji,jj,jk,jp_sal) = ptsa(ji,jj,jk,jp_sal) + ptsc(ji,jj,jk,jp_sal) * r1_e1e2t(ji,jj) / e3t_0(ji,jj,jk)
      END DO   ;   END DO   ;   END DO
      !
   END SUBROUTINE tra_isf_cpl
   !
END MODULE traisf

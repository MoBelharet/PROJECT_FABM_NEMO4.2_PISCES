










MODULE trcatf
   !!======================================================================
   !!                       ***  MODULE  trcatf  ***
   !! Ocean passive tracers:  time stepping on passives tracers
   !!======================================================================
   !! History :  7.0  !  1991-11  (G. Madec)  Original code
   !!                 !  1993-03  (M. Guyon)  symetrical conditions
   !!                 !  1995-02  (M. Levy)   passive tracers
   !!                 !  1996-02  (G. Madec & M. Imbard)  opa release 8.0
   !!            8.0  !  1996-04  (A. Weaver)  Euler forward step
   !!            8.2  !  1999-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!                 !  2002-08  (G. Madec)  F90: Free form and module
   !!                 !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!                 !  2004-03  (C. Ethe) passive tracers
   !!                 !  2007-02  (C. Deltel) Diagnose ML trends for passive tracers
   !!            2.0  !  2006-02  (L. Debreu, C. Mazauric) Agrif implementation
   !!            3.0  !  2008-06  (G. Madec)  time stepping always done in trazdf
   !!            3.1  !  2009-02  (G. Madec, R. Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-06  (C. Ethe, G. Madec) Merge TRA-TRC
   !!            4.1  !  2019-08  (A. Coward, D. Storkey) rename trcnxt.F90 -> trcatf.F90. Now only does time filtering.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   !!   trc_atf       : time stepping on passive tracers
   !!----------------------------------------------------------------------
   USE par_trc        ! need jptra, number of passive tracers
   USE oce_trc        ! ocean dynamics and tracers variables
   USE trc            ! ocean passive tracers variables
   USE trd_oce
   USE trdtra
   USE traatf_qco     ! tracer : Asselin filter (qco)
   USE bdy_oce   , ONLY: ln_bdy
   USE trcbdy         ! BDY open boundaries
   !
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control for debbuging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trc_atf   ! routine called by step.F90

   REAL(wp) ::   rfact1, rfact2

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
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcatf.F90 15090 2021-07-06 14:25:18Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE trc_atf( kt, Kbb, Kmm, Kaa, ptr )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE trcatf  ***
      !!
      !! ** Purpose :   Apply the boundary condition on the after passive tracers fields and
      !!      apply Asselin time filter to the now passive tracer fields if using leapfrog timestep
      !! 
      !! ** Method  :   Apply lateral boundary conditions on (uu(Kaa),vv(Kaa)) through 
      !!      call to lbc_lnk routine
      !!
      !!   For Arakawa or TVD Scheme : 
      !!      A Asselin time filter applied on now tracers tr(Kmm) to avoid
      !!      the divergence of two consecutive time-steps and tr arrays
      !!      to prepare the next time_step:
      !!         (tr(Kmm)) = (tr(Kmm)) + rn_atfp [ (tr(Kbb)) + (tr(Kaa)) - 2 (tr(Kmm)) ]
      !!
      !!
      !! ** Action  : - update tr(Kmm), tr(Kaa)
      !!----------------------------------------------------------------------
      INTEGER                                   , INTENT( in )  :: kt             ! ocean time-step index
      INTEGER                                   , INTENT( in )  :: Kbb, Kmm, Kaa ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra,jpt), INTENT(inout) :: ptr            ! passive tracers
      !
      INTEGER  ::   jk, jn   ! dummy loop indices
      REAL(wp) ::   zfact            ! temporary scalar
      CHARACTER (len=22) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   ztrdt    ! 4D workspace
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_atf')
      !
      IF( kt == nittrc000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'trc_atf : Asselin time filtering on passive tracers'
      ENDIF
      !
      ! Update after tracer on domain lateral boundaries
      CALL lbc_lnk( 'trcatf', ptr(:,:,:,:,Kaa), 'T', 1._wp )   

      IF( ln_bdy )  CALL trc_bdy( kt, Kbb, Kmm, Kaa )

      IF( l_trdtrc )  THEN             ! trends: store now fields before the Asselin filter application
         ALLOCATE( ztrdt(jpi,jpj,jpk,jptra) )
         ztrdt(:,:,:,:)  = 0._wp
         IF( ln_traldf_iso ) THEN                       ! diagnose the "pure" Kz diffusive trend 
            DO jn = 1, jptra
               CALL trd_tra( kt, Kmm, Kaa, 'TRC', jn, jptra_zdfp, ztrdt(:,:,:,jn) )
            ENDDO
         ENDIF

         ! total trend for the non-time-filtered variables. 
         zfact = 1.0 / rn_Dt
         ! G Nurser 23 Mar 2017. Recalculate trend as Delta(e3ta*Ta)/e3tn; e3tn cancel from ts(Kmm) terms
         IF( ln_linssh ) THEN       ! linear sea surface height only
            DO jn = 1, jptra
               DO jk = 1, jpkm1
                  ztrdt(:,:,jk,jn) = ( ptr(:,:,jk,jn,Kaa)*e3t_0(:,:,jk) / e3t_0(:,:,jk) - ptr(:,:,jk,jn,Kmm)) * zfact
               END DO
            END DO
         ELSE
            DO jn = 1, jptra
               DO jk = 1, jpkm1
                  ztrdt(:,:,jk,jn) = ( ptr(:,:,jk,jn,Kaa) - ptr(:,:,jk,jn,Kmm) ) * zfact
               END DO
            END DO
         ENDIF
         !
         DO jn = 1, jptra
            CALL trd_tra( kt, Kmm, Kaa, 'TRC', jn, jptra_tot, ztrdt(:,:,:,jn) )
         ENDDO
         !
         IF( ln_linssh ) THEN       ! linear sea surface height only
            ! Store now fields before applying the Asselin filter 
            ! in order to calculate Asselin filter trend later.
            ztrdt(:,:,:,:) = ptr(:,:,:,:,Kmm) 
         ENDIF

      ENDIF
      !                                ! Leap-Frog + Asselin filter time stepping
      IF( l_1st_euler .OR. ln_top_euler ) THEN    ! Euler time-stepping 
         !
         IF (l_trdtrc .AND. .NOT. ln_linssh ) THEN   ! Zero Asselin filter contribution must be explicitly written out since for vvl
            !                                        ! Asselin filter is output by tra_nxt_vvl that is not called on this time step
            ztrdt(:,:,:,:) = 0._wp            
            DO jn = 1, jptra
               CALL trd_tra( kt, Kmm, Kaa, 'TRC', jn, jptra_atf, ztrdt(:,:,:,jn) )
            ENDDO
         END IF
         !
      ELSE     
         IF( .NOT. l_offline ) THEN ! Leap-Frog + Asselin filter time stepping
            IF( ln_linssh ) THEN   ;   CALL tra_atf_fix_lf( kt, Kbb, Kmm, Kaa, nittrc000,        'TRC', ptr, jptra )                     !     linear ssh
            ELSE                   ;   CALL tra_atf_qco_lf( kt, Kbb, Kmm, Kaa, nittrc000, rn_Dt, 'TRC', ptr, sbc_trc, sbc_trc_b, jptra ) ! non-linear ssh
            ENDIF
         ELSE
                                       CALL trc_atf_off( kt, Kbb, Kmm, Kaa, ptr )       ! offline 
         ENDIF
         !
         CALL lbc_lnk( 'trcatf', ptr(:,:,:,:,Kmm), 'T', 1._wp )
      ENDIF
      !
      IF( l_trdtrc .AND. ln_linssh ) THEN      ! trend of the Asselin filter (tb filtered - tb)/dt )
         DO jn = 1, jptra
            DO jk = 1, jpkm1
               zfact = 1._wp / rDt_trc  
               ztrdt(:,:,jk,jn) = ( ptr(:,:,jk,jn,Kbb) - ztrdt(:,:,jk,jn) ) * zfact 
            END DO
            CALL trd_tra( kt, Kmm, Kaa, 'TRC', jn, jptra_atf, ztrdt(:,:,:,jn) )
         END DO
      END IF
      IF( l_trdtrc ) DEALLOCATE( ztrdt ) 
      !
      IF(sn_cfctl%l_prttrc)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('nxt')")
         CALL prt_ctl_info( charout, cdcomp = 'top' )
         CALL prt_ctl(tab4d_1=ptr(:,:,:,:,Kmm), mask1=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('trc_atf')
      !
   END SUBROUTINE trc_atf

   SUBROUTINE trc_atf_off( kt, Kbb, Kmm, Kaa, ptr )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_atf_off  ***
      !!
      !!          !!!!!!!!!!!!!!!!! REWRITE HEADER COMMENTS !!!!!!!!!!!!!!
      !!
      !! ** Purpose :   Time varying volume: apply the Asselin time filter  
      !! 
      !! ** Method  : - Apply a thickness weighted Asselin time filter on now fields.
      !!              - save in (ta,sa) a thickness weighted average over the three 
      !!             time levels which will be used to compute rdn and thus the semi-
      !!             implicit hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = ( e3t_n*tn + rbcp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )   ln_dynhpg_imp = T
      !!                  /( e3t_0(:,:,jk)    + rbcp*[ e3t_0(:,:,jk)    - 2 e3t_0(:,:,jk)    + e3t_0(:,:,jk)    ] )   
      !!             ztm = 0                                                       otherwise
      !!             tb  = ( e3t_n*tn + rn_atfp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )
      !!                  /( e3t_0(:,:,jk)    + rn_atfp*[ e3t_0(:,:,jk)    - 2 e3t_0(:,:,jk)    + e3t_0(:,:,jk)    ] )
      !!             tn  = ta 
      !!             ta  = zt        (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER                                   , INTENT(in   ) ::  kt            ! ocean time-step index
      INTEGER                                   , INTENT(in   ) ::  Kbb, Kmm, Kaa ! time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jptra,jpt), INTENT(inout) ::  ptr           ! passive tracers
      !!     
      INTEGER  ::   ji, jj, jk, jn              ! dummy loop indices
      REAL(wp) ::   ztc_a , ztc_n , ztc_b , ztc_f , ztc_d    ! local scalar
      REAL(wp) ::   ze3t_b, ze3t_n, ze3t_a, ze3t_f           !   -      -
      !!----------------------------------------------------------------------
      !
      IF( kt == nittrc000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'trc_atf_off : Asselin time filtering'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF( .NOT. ln_linssh ) THEN
            rfact1 = rn_atfp * rn_Dt
            rfact2 = rfact1 / rho0
         ENDIF
        !  
      ENDIF
      !
      DO jn = 1, jptra  
         DO jk =  1,  jpkm1  ; DO jj = ntsj-(  nn_hls), ntej+(  nn_hls) ; DO ji = ntsi-( nn_hls), ntei+(  nn_hls)    
            ze3t_b = 1._wp + r3t(ji,jj,Kbb) * tmask(ji,jj,jk)
            ze3t_n = 1._wp + r3t(ji,jj,Kmm) * tmask(ji,jj,jk)
            ze3t_a = 1._wp + r3t(ji,jj,Kaa) * tmask(ji,jj,jk)
            !                                         ! tracer content at Before, now and after
            ztc_b  = ptr(ji,jj,jk,jn,Kbb) * ze3t_b
            ztc_n  = ptr(ji,jj,jk,jn,Kmm) * ze3t_n
            ztc_a  = ptr(ji,jj,jk,jn,Kaa) * ze3t_a
            !
            ztc_d  = ztc_a  - 2. * ztc_n  + ztc_b
            !
            ze3t_f = 1._wp + r3t_f(ji,jj)*tmask(ji,jj,jk)
            ztc_f  = ztc_n  + rn_atfp * ztc_d
            !
            IF( .NOT. ln_linssh .AND. jk == mikt(ji,jj) ) THEN           ! first level 
               ztc_f  = ztc_f  - rfact1 * ( sbc_trc(ji,jj,jn) - sbc_trc_b(ji,jj,jn) )
            ENDIF

            ze3t_f = 1.e0 / ze3t_f
            ptr(ji,jj,jk,jn,Kmm) = ztc_f * ze3t_f     ! time filtered "now" field
            !
         END DO   ;   END DO   ;   END DO
         ! 
      END DO
      !
   END SUBROUTINE trc_atf_off

   !!======================================================================
END MODULE trcatf
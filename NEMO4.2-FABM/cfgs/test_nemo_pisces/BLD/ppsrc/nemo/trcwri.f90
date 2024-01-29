










MODULE trcwri
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    TOP :   Output of passive tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_top'                                           TOP models
   !!----------------------------------------------------------------------
   !! trc_wri_trc   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE dom_oce     ! ocean space and time domain variables
   USE oce_trc     ! shared variables between ocean and passive tracers
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager
   USE dianam      ! Output file name
   USE trcwri_pisces
   USE trcwri_cfc
   USE trcwri_c14
   USE trcwri_age
   USE trcwri_my_trc

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri      

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

CONTAINS

   SUBROUTINE trc_wri( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri  ***
      !! 
      !! ** Purpose :   output passive tracers fields and dynamical trends
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in )     :: kt
      INTEGER, INTENT( in )     :: Kmm  ! time level indices
      !
      INTEGER                   :: jk, jn
      CHARACTER (len=20)        :: cltra
      CHARACTER (len=40)        :: clhstnam
      INTEGER ::   inum = 11            ! temporary logical unit
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   z3d   ! 3D workspace
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('trc_wri')
      !
      IF( l_offline ) THEN    ! WRITE root name in date.file for use by postpro
         IF(  kt == nittrc000 .AND. lwp ) THEN    ! WRITE root name in date.file for use by postpro
           CALL dia_nam( clhstnam, nn_writetrc,' ' )
           CALL ctl_opn( inum, 'date.file', 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
           WRITE(inum,*) clhstnam
           CLOSE(inum)
         ENDIF

         ! Output of initial vertical scale factor
         CALL iom_put( "e3t_0", e3t_0(:,:,:) )
         CALL iom_put( "e3u_0", e3u_0(:,:,:) )
         CALL iom_put( "e3v_0", e3v_0(:,:,:) )
         !
         IF( .NOT.ln_linssh )  CALL iom_put( "ssh" , ssh(:,:,Kmm) )              ! sea surface height
         !
         IF ( iom_use("e3t") ) THEN  ! time-varying e3t
            DO jk = 1, jpk
               z3d(:,:,jk) =  e3t_0(:,:,jk)
            END DO
            CALL iom_put( "e3t", z3d(:,:,:) )
         ENDIF
         IF ( iom_use("e3u") ) THEN                         ! time-varying e3u
            DO jk = 1, jpk
               z3d(:,:,jk) =  e3u_0(:,:,jk)
            END DO
            CALL iom_put( "e3u", z3d(:,:,:) )
         ENDIF
         IF ( iom_use("e3v") ) THEN                         ! time-varying e3v
            DO jk = 1, jpk
               z3d(:,:,jk) =  e3v_0(:,:,jk)
            END DO
            CALL iom_put( "e3v", z3d(:,:,:) )
         ENDIF
         !
      ENDIF
      !
      ! write the tracer concentrations in the file
      ! ---------------------------------------
      IF( ln_pisces  )   CALL trc_wri_pisces( Kmm )     ! PISCES 
      IF( ll_cfc     )   CALL trc_wri_cfc   ( Kmm )     ! surface fluxes of CFC
      IF( ln_c14     )   CALL trc_wri_c14   ( Kmm )     ! surface fluxes of C14
      IF( ln_age     )   CALL trc_wri_age   ( Kmm )     ! AGE tracer
      IF( ln_my_trc  )   CALL trc_wri_my_trc( Kmm )     ! MY_TRC  tracers
      !
      IF( ln_timing )   CALL timing_stop('trc_wri')
      !
   END SUBROUTINE trc_wri


   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trcwri.F90 14255 2021-01-04 11:35:00Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE trcwri

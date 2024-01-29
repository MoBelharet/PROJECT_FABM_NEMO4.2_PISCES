










MODULE istate
   !!======================================================================
   !!                     ***  MODULE  istate  ***
   !! Ocean state   :  initial state setting
   !!=====================================================================
   !! History :  OPA  !  1989-12  (P. Andrich)  Original code
   !!            5.0  !  1991-11  (G. Madec)  rewritting
   !!            6.0  !  1996-01  (G. Madec)  terrain following coordinates
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_eel
   !!            8.0  !  2001-09  (M. Levy, M. Ben Jelloul)  istate_uvg
   !!   NEMO     1.0  !  2003-08  (G. Madec, C. Talandier)  F90: Free form, modules + EEL R5
   !!             -   !  2004-05  (A. Koch-Larrouy)  istate_gyre 
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!            3.3  !  2010-10  (C. Ethe) merge TRC-TRA
   !!            3.4  !  2011-04  (G. Madec) Merge of dtatem and dtasal & suppression of tb,tn/sb,sn 
   !!            3.7  !  2016-04  (S. Flavoni) introduce user defined initial state 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   istate_init   : initial state setting
   !!   istate_uvg    : initial velocity in geostropic balance
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers 
   USE dom_oce        ! ocean space and time domain 
   USE daymod         ! calendar
   USE dtatsd         ! data temperature and salinity   (dta_tsd routine)
   USE dtauvd         ! data: U & V current             (dta_uvd routine)
   USE domvvl          ! varying vertical mesh
   USE wet_dry         ! wetting and drying (needed for wad_istate)
   USE usrdef_istate   ! User defined initial state
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library
   USE lbclnk         ! lateal boundary condition / mpp exchanges
   USE restart         ! restart


   IMPLICIT NONE
   PRIVATE

   PUBLIC   istate_init   ! routine called by nemogcm.F90

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
   !! $Id: istate.F90 15052 2021-06-24 14:39:14Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE istate_init( Kbb, Kmm, Kaa )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  ***
      !! 
      !! ** Purpose :   Initialization of the dynamics and tracer fields.
      !!
      !! ** Method  :   
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in )  ::  Kbb, Kmm, Kaa   ! ocean time level indices
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   zgdept     ! 3D table for qco substitute
!!gm see comment further down
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) ::   zuvd    ! U & V data workspace
!!gm end
      !!----------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_init : Initialization of the dynamics and tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

      CALL dta_tsd_init                 ! Initialisation of T & S input data
      IF( ln_c1d) CALL dta_uvd_init     ! Initialisation of U & V input data (c1d only)

      ts   (:,:,:,:,Kaa) = 0._wp   ;   rn2  (:,:,:  ) = 0._wp         ! set one for all to 0 at levels 1 and jpk
      IF ( ALLOCATED( rhd ) ) THEN                                    ! SWE, for example, will not have allocated these
         rhd  (:,:,:      ) = 0._wp   ;   rhop (:,:,:  ) = 0._wp      ! set one for all to 0 at level jpk
         rn2b (:,:,:      ) = 0._wp                                   ! set one for all to 0 at level jpk
         rab_b(:,:,:,:    ) = 0._wp   ;   rab_n(:,:,:,:) = 0._wp      ! set one for all to 0 at level jpk
      ENDIF

         IF( ln_rstart ) THEN                    ! Restart from a file
            !                                    ! -------------------
            CALL rst_read( Kbb, Kmm )            ! Read the restart file
            CALL day_init                        ! model calendar (using both namelist and restart infos)
            !
         ELSE                                    ! Start from rest
            !                                    ! ---------------
            numror = 0                           ! define numror = 0 -> no restart file to read
            l_1st_euler = .true.                 ! Set time-step indicator at nit000 (euler forward)
            CALL day_init                        ! model calendar (using both namelist and restart infos)
            !                                    ! Initialization of ocean to zero
            !
            IF( ln_tsd_init ) THEN               
               CALL dta_tsd( nit000, ts(:,:,:,:,Kbb) )                     ! read 3D T and S data at nit000
            ENDIF
            !
            IF( ln_uvd_init .AND. ln_c1d ) THEN               
               CALL dta_uvd( nit000, Kbb, uu(:,:,:,Kbb), vv(:,:,:,Kbb) )   ! read 3D U and V data at nit000
            ELSE
               uu  (:,:,:,Kbb) = 0._wp               ! set the ocean at rest
               vv  (:,:,:,Kbb) = 0._wp  
            ENDIF
               !
               !
            IF( .NOT. ln_tsd_init .AND. .NOT. ln_uvd_init ) THEN
               DO jk = 1, jpk
                  zgdept(:,:,jk) = gdept_0(:,:,jk)
               END DO
               CALL usr_def_istate( zgdept, tmask, ts(:,:,:,:,Kbb), uu(:,:,:,Kbb), vv(:,:,:,Kbb) )
               ! make sure that periodicities are properly applied 
               CALL lbc_lnk( 'istate', ts(:,:,:,jp_tem,Kbb), 'T',  1._dp, ts(:,:,:,jp_sal,Kbb), 'T',  1._dp,   &
                  &                    uu(:,:,:,       Kbb), 'U', -1._dp, vv(:,:,:,       Kbb), 'V', -1._dp )
            ENDIF
            ts  (:,:,:,:,Kmm) = ts (:,:,:,:,Kbb)       ! set now values from to before ones
            uu    (:,:,:,Kmm) = uu   (:,:,:,Kbb)
            vv    (:,:,:,Kmm) = vv   (:,:,:,Kbb)

         ENDIF 
      ! 
      ! Initialize "now" and "before" barotropic velocities:
      ! Do it whatever the free surface method, these arrays being eventually used
      !
      uu_b(:,:,Kmm) = 0._wp   ;   vv_b(:,:,Kmm) = 0._wp
      uu_b(:,:,Kbb) = 0._wp   ;   vv_b(:,:,Kbb) = 0._wp
      !
!!gm  the use of umsak & vmask is not necessary below as uu(:,:,:,Kmm), vv(:,:,:,Kmm), uu(:,:,:,Kbb), vv(:,:,:,Kbb) are always masked
      DO jk =  1,  jpkm1  ; DO jj = ntsj-(  nn_hls), ntej+(  nn_hls) ; DO ji = ntsi-( nn_hls), ntei+(  nn_hls)
         uu_b(ji,jj,Kmm) = uu_b(ji,jj,Kmm) + e3u_0(ji,jj,jk) * uu(ji,jj,jk,Kmm) * umask(ji,jj,jk)
         vv_b(ji,jj,Kmm) = vv_b(ji,jj,Kmm) + e3v_0(ji,jj,jk) * vv(ji,jj,jk,Kmm) * vmask(ji,jj,jk)
         !
         uu_b(ji,jj,Kbb) = uu_b(ji,jj,Kbb) + e3u_0(ji,jj,jk) * uu(ji,jj,jk,Kbb) * umask(ji,jj,jk)
         vv_b(ji,jj,Kbb) = vv_b(ji,jj,Kbb) + e3v_0(ji,jj,jk) * vv(ji,jj,jk,Kbb) * vmask(ji,jj,jk)
      END DO   ;   END DO   ;   END DO
      !
      uu_b(:,:,Kmm) = uu_b(:,:,Kmm) * r1_hu_0(:,:)
      vv_b(:,:,Kmm) = vv_b(:,:,Kmm) * r1_hv_0(:,:)
      !
      uu_b(:,:,Kbb) = uu_b(:,:,Kbb) * r1_hu_0(:,:)
      vv_b(:,:,Kbb) = vv_b(:,:,Kbb) * r1_hv_0(:,:)
      !
   END SUBROUTINE istate_init

   !!======================================================================
END MODULE istate

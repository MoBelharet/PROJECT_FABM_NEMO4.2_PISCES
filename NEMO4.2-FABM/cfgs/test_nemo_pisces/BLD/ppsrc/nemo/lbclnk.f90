










MODULE lbclnk
   !!======================================================================
   !!                       ***  MODULE  lbclnk  ***
   !! NEMO        : lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)  Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment
   !!            3.5  ! 2012     (S.Mocavero, I. Epicoco)  optimization of BDY comm. via lbc_bdy_lnk and lbc_obc_lnk
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie, G. Reffray)  add a C1D case
   !!            3.6  ! 2015-06  (O. Tintó and M. Castrillo)  add lbc_lnk_multi
   !!            4.0  ! 2017-03  (G. Madec) automatique allocation of array size (use with any 3rd dim size)
   !!             -   ! 2017-04  (G. Madec) remove duplicated routines (lbc_lnk_2d_9, lbc_lnk_2d_multiple, lbc_lnk_3d_gather)
   !!             -   ! 2017-05  (G. Madec) create generic.h90 files to generate all lbc and north fold routines
   !!----------------------------------------------------------------------
   !!           define the generic interfaces of lib_mpp routines
   !!----------------------------------------------------------------------
   !!   lbc_lnk       : generic interface for mpp_lnk_3d and mpp_lnk_2d routines defined in lib_mpp
   !!   lbc_bdy_lnk   : generic interface for mpp_lnk_bdy_2d and mpp_lnk_bdy_3d routines defined in lib_mpp
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE lib_mpp        ! distributed memory computing library
   USE lbcnfd         ! north fold
   USE in_out_manager ! I/O manager
   USE MPI

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_lnk
      MODULE PROCEDURE   lbc_lnk_call_2d_sp, lbc_lnk_call_3d_sp, lbc_lnk_call_4d_sp
      MODULE PROCEDURE   lbc_lnk_call_2d_dp, lbc_lnk_call_3d_dp, lbc_lnk_call_4d_dp
   END INTERFACE

   INTERFACE lbc_lnk_pt2pt
      MODULE PROCEDURE   lbc_lnk_pt2pt_sp, lbc_lnk_pt2pt_dp
   END INTERFACE

   INTERFACE lbc_lnk_neicoll
      MODULE PROCEDURE   lbc_lnk_neicoll_sp ,lbc_lnk_neicoll_dp
   END INTERFACE
   !
   INTERFACE lbc_lnk_icb
      MODULE PROCEDURE mpp_lnk_2d_icb_dp, mpp_lnk_2d_icb_sp
   END INTERFACE

   PUBLIC   lbc_lnk            ! ocean/ice lateral boundary conditions
   PUBLIC   lbc_lnk_icb        ! iceberg lateral boundary conditions

   REAL(dp), DIMENSION(:), ALLOCATABLE ::   buffsnd_dp, buffrcv_dp   ! MPI send/recv buffers
   REAL(sp), DIMENSION(:), ALLOCATABLE ::   buffsnd_sp, buffrcv_sp   ! 
   INTEGER,  DIMENSION(8)              ::   nreq_p2p                 ! request id for MPI_Isend in point-2-point communication
   
   !! * Substitutions
   !!#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: lbclnk.F90 14433 2021-02-11 08:06:49Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!                   ***   lbc_lnk_call_[234]d_[sd]p   ***
   !!
   !!   * Dummy Argument :
   !!       in    ==>   cdname     ! name of the calling subroutine (for monitoring)
   !!                   ptab       ! array to be loaded (2D, 3D or 4D)
   !!                   cd_nat     ! nature of pt2d array grid-points
   !!                   psgn       ! sign used across the north fold boundary
   !!       inout <=>   ptab_ptr   ! array of 2D, 3D or 4D pointers
   !!                   cdna_ptr   ! nature of ptab array grid-points
   !!                   psgn_ptr   ! sign used across the north fold boundary
   !!                   kfld       ! number of elements that has been attributed
   !!----------------------------------------------------------------------
   !
   !!----------------------------------------------------------------------
   !!
   !!                  ***   lbc_lnk_call_[234]d_[sd]p   ***
   !!                  ***     load_ptr_[234]d_[sd]p     ***
   !!
   !!----------------------------------------------------------------------
   !!
   !!   ----   SINGLE PRECISION VERSIONS
   !!

   SUBROUTINE lbc_lnk_call_2d_sp(                                                                &
      &                     cdname                                                                                  &
      &                   , pt1 , cdna1 , psgn1 , pt2 , cdna2 , psgn2 , pt3 , cdna3 , psgn3 , pt4 , cdna4 , psgn4   &
      &                   , pt5 , cdna5 , psgn5 , pt6 , cdna6 , psgn6 , pt7 , cdna7 , psgn7 , pt8 , cdna8 , psgn8   &
      &                   , pt9 , cdna9 , psgn9 , pt10, cdna10, psgn10, pt11, cdna11, psgn11, pt12, cdna12, psgn12  &
      &                   , pt13, cdna13, psgn13, pt14, cdna14, psgn14, pt15, cdna15, psgn15, pt16, cdna16, psgn16  &
      &                   , pt17, cdna17, psgn17, pt18, cdna18, psgn18, pt19, cdna19, psgn19, pt20, cdna20, psgn20  &
      &                   , pt21, cdna21, psgn21, pt22, cdna22, psgn22, pt23, cdna23, psgn23, pt24, cdna24, psgn24  &
      &                   , pt25, cdna25, psgn25, pt26, cdna26, psgn26, pt27, cdna27, psgn27, pt28, cdna28, psgn28  &
      &                   , pt29, cdna29, psgn29, pt30, cdna30, psgn30                                              &
      &                   , kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      !!---------------------------------------------------------------------
      CHARACTER(len=*)     ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(sp), DIMENSION(:,:)          , TARGET, CONTIGUOUS, INTENT(inout) ::   pt1        ! arrays on which the lbc is applied
      REAL(sp), DIMENSION(:,:), OPTIONAL, TARGET, CONTIGUOUS, INTENT(inout) ::   pt2 , pt3 , pt4 , pt5 , pt6 , pt7 , pt8 , &
         &                                                                               pt9 , pt10, pt11, pt12, pt13, pt14, pt15, &
         &                                                                               pt16, pt17, pt18, pt19, pt20, pt21, pt22, &
         &                                                                               pt23, pt24, pt25, pt26, pt27, pt28, pt29, &
         &                                                                               pt30
      CHARACTER(len=1)                       , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)     , OPTIONAL        , INTENT(in   ) ::   cdna2 , cdna3 , cdna4 , cdna5 , cdna6 , cdna7 , cdna8 , &
         &                                                        cdna9 , cdna10, cdna11, cdna12, cdna13, cdna14, cdna15, &
         &                                                        cdna16, cdna17, cdna18, cdna19, cdna20, cdna21, cdna22, &
         &                                                        cdna23, cdna24, cdna25, cdna26, cdna27, cdna28, cdna29, &
         &                                                        cdna30
      REAL(sp)                        , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(sp)      , OPTIONAL        , INTENT(in   ) ::   psgn2 , psgn3 , psgn4 , psgn5 , psgn6 , psgn7 , psgn8 , &
         &                                                        psgn9 , psgn10, psgn11, psgn12, psgn13, psgn14, psgn15, &
         &                                                        psgn16, psgn17, psgn18, psgn19, psgn20, psgn21, psgn22, &
         &                                                        psgn23, psgn24, psgn25, psgn26, psgn27, psgn28, psgn29, &
         &                                                        psgn30
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(sp)      , OPTIONAL        , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8), OPTIONAL        , INTENT(in   ) ::   lsend, lrecv   ! indicate how communications are to be carried out
      LOGICAL              , OPTIONAL        , INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !!
      INTEGER                          ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_4d_sp), DIMENSION(30) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(30) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(sp)  , DIMENSION(30) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_2d_sp( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2 ) )   CALL load_ptr_2d_sp( pt2 , cdna2 , psgn2 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3 ) )   CALL load_ptr_2d_sp( pt3 , cdna3 , psgn3 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4 ) )   CALL load_ptr_2d_sp( pt4 , cdna4 , psgn4 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5 ) )   CALL load_ptr_2d_sp( pt5 , cdna5 , psgn5 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6 ) )   CALL load_ptr_2d_sp( pt6 , cdna6 , psgn6 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7 ) )   CALL load_ptr_2d_sp( pt7 , cdna7 , psgn7 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8 ) )   CALL load_ptr_2d_sp( pt8 , cdna8 , psgn8 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9 ) )   CALL load_ptr_2d_sp( pt9 , cdna9 , psgn9 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn10) )   CALL load_ptr_2d_sp( pt10, cdna10, psgn10, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn11) )   CALL load_ptr_2d_sp( pt11, cdna11, psgn11, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn12) )   CALL load_ptr_2d_sp( pt12, cdna12, psgn12, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn13) )   CALL load_ptr_2d_sp( pt13, cdna13, psgn13, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn14) )   CALL load_ptr_2d_sp( pt14, cdna14, psgn14, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn15) )   CALL load_ptr_2d_sp( pt15, cdna15, psgn15, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn16) )   CALL load_ptr_2d_sp( pt16, cdna16, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn17) )   CALL load_ptr_2d_sp( pt17, cdna17, psgn17, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn18) )   CALL load_ptr_2d_sp( pt18, cdna18, psgn18, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn19) )   CALL load_ptr_2d_sp( pt19, cdna19, psgn19, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn20) )   CALL load_ptr_2d_sp( pt20, cdna20, psgn20, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn21) )   CALL load_ptr_2d_sp( pt21, cdna21, psgn21, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn22) )   CALL load_ptr_2d_sp( pt22, cdna22, psgn22, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn23) )   CALL load_ptr_2d_sp( pt23, cdna23, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn24) )   CALL load_ptr_2d_sp( pt24, cdna24, psgn24, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn25) )   CALL load_ptr_2d_sp( pt25, cdna25, psgn25, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn26) )   CALL load_ptr_2d_sp( pt26, cdna26, psgn26, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn27) )   CALL load_ptr_2d_sp( pt27, cdna27, psgn27, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn28) )   CALL load_ptr_2d_sp( pt28, cdna28, psgn28, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn29) )   CALL load_ptr_2d_sp( pt29, cdna29, psgn29, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn30) )   CALL load_ptr_2d_sp( pt30, cdna30, psgn30, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      IF( nn_comm == 1 ) THEN 
         CALL lbc_lnk_pt2pt(   cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ELSE
         CALL lbc_lnk_neicoll( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ENDIF
      !
   END SUBROUTINE lbc_lnk_call_2d_sp


   SUBROUTINE load_ptr_2d_sp( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(sp), DIMENSION(:,:), TARGET, INTENT(inout), CONTIGUOUS ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(sp)               , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_4d_sp), DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(sp) , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt4d(1:SIZE(ptab, dim=1),1:SIZE(ptab, dim=2),1:1,1:1) => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_2d_sp


   SUBROUTINE lbc_lnk_call_3d_sp(                                                                &
      &                     cdname                                                                                  &
      &                   , pt1 , cdna1 , psgn1 , pt2 , cdna2 , psgn2 , pt3 , cdna3 , psgn3 , pt4 , cdna4 , psgn4   &
      &                   , pt5 , cdna5 , psgn5 , pt6 , cdna6 , psgn6 , pt7 , cdna7 , psgn7 , pt8 , cdna8 , psgn8   &
      &                   , pt9 , cdna9 , psgn9 , pt10, cdna10, psgn10, pt11, cdna11, psgn11, pt12, cdna12, psgn12  &
      &                   , pt13, cdna13, psgn13, pt14, cdna14, psgn14, pt15, cdna15, psgn15, pt16, cdna16, psgn16  &
      &                   , pt17, cdna17, psgn17, pt18, cdna18, psgn18, pt19, cdna19, psgn19, pt20, cdna20, psgn20  &
      &                   , pt21, cdna21, psgn21, pt22, cdna22, psgn22, pt23, cdna23, psgn23, pt24, cdna24, psgn24  &
      &                   , pt25, cdna25, psgn25, pt26, cdna26, psgn26, pt27, cdna27, psgn27, pt28, cdna28, psgn28  &
      &                   , pt29, cdna29, psgn29, pt30, cdna30, psgn30                                              &
      &                   , kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      !!---------------------------------------------------------------------
      CHARACTER(len=*)     ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(sp), DIMENSION(:,:,:)          , TARGET, CONTIGUOUS, INTENT(inout) ::   pt1        ! arrays on which the lbc is applied
      REAL(sp), DIMENSION(:,:,:), OPTIONAL, TARGET, CONTIGUOUS, INTENT(inout) ::   pt2 , pt3 , pt4 , pt5 , pt6 , pt7 , pt8 , &
         &                                                                               pt9 , pt10, pt11, pt12, pt13, pt14, pt15, &
         &                                                                               pt16, pt17, pt18, pt19, pt20, pt21, pt22, &
         &                                                                               pt23, pt24, pt25, pt26, pt27, pt28, pt29, &
         &                                                                               pt30
      CHARACTER(len=1)                       , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)     , OPTIONAL        , INTENT(in   ) ::   cdna2 , cdna3 , cdna4 , cdna5 , cdna6 , cdna7 , cdna8 , &
         &                                                        cdna9 , cdna10, cdna11, cdna12, cdna13, cdna14, cdna15, &
         &                                                        cdna16, cdna17, cdna18, cdna19, cdna20, cdna21, cdna22, &
         &                                                        cdna23, cdna24, cdna25, cdna26, cdna27, cdna28, cdna29, &
         &                                                        cdna30
      REAL(sp)                        , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(sp)      , OPTIONAL        , INTENT(in   ) ::   psgn2 , psgn3 , psgn4 , psgn5 , psgn6 , psgn7 , psgn8 , &
         &                                                        psgn9 , psgn10, psgn11, psgn12, psgn13, psgn14, psgn15, &
         &                                                        psgn16, psgn17, psgn18, psgn19, psgn20, psgn21, psgn22, &
         &                                                        psgn23, psgn24, psgn25, psgn26, psgn27, psgn28, psgn29, &
         &                                                        psgn30
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(sp)      , OPTIONAL        , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8), OPTIONAL        , INTENT(in   ) ::   lsend, lrecv   ! indicate how communications are to be carried out
      LOGICAL              , OPTIONAL        , INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !!
      INTEGER                          ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_4d_sp), DIMENSION(30) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(30) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(sp)  , DIMENSION(30) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_3d_sp( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2 ) )   CALL load_ptr_3d_sp( pt2 , cdna2 , psgn2 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3 ) )   CALL load_ptr_3d_sp( pt3 , cdna3 , psgn3 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4 ) )   CALL load_ptr_3d_sp( pt4 , cdna4 , psgn4 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5 ) )   CALL load_ptr_3d_sp( pt5 , cdna5 , psgn5 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6 ) )   CALL load_ptr_3d_sp( pt6 , cdna6 , psgn6 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7 ) )   CALL load_ptr_3d_sp( pt7 , cdna7 , psgn7 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8 ) )   CALL load_ptr_3d_sp( pt8 , cdna8 , psgn8 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9 ) )   CALL load_ptr_3d_sp( pt9 , cdna9 , psgn9 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn10) )   CALL load_ptr_3d_sp( pt10, cdna10, psgn10, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn11) )   CALL load_ptr_3d_sp( pt11, cdna11, psgn11, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn12) )   CALL load_ptr_3d_sp( pt12, cdna12, psgn12, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn13) )   CALL load_ptr_3d_sp( pt13, cdna13, psgn13, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn14) )   CALL load_ptr_3d_sp( pt14, cdna14, psgn14, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn15) )   CALL load_ptr_3d_sp( pt15, cdna15, psgn15, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn16) )   CALL load_ptr_3d_sp( pt16, cdna16, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn17) )   CALL load_ptr_3d_sp( pt17, cdna17, psgn17, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn18) )   CALL load_ptr_3d_sp( pt18, cdna18, psgn18, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn19) )   CALL load_ptr_3d_sp( pt19, cdna19, psgn19, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn20) )   CALL load_ptr_3d_sp( pt20, cdna20, psgn20, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn21) )   CALL load_ptr_3d_sp( pt21, cdna21, psgn21, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn22) )   CALL load_ptr_3d_sp( pt22, cdna22, psgn22, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn23) )   CALL load_ptr_3d_sp( pt23, cdna23, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn24) )   CALL load_ptr_3d_sp( pt24, cdna24, psgn24, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn25) )   CALL load_ptr_3d_sp( pt25, cdna25, psgn25, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn26) )   CALL load_ptr_3d_sp( pt26, cdna26, psgn26, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn27) )   CALL load_ptr_3d_sp( pt27, cdna27, psgn27, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn28) )   CALL load_ptr_3d_sp( pt28, cdna28, psgn28, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn29) )   CALL load_ptr_3d_sp( pt29, cdna29, psgn29, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn30) )   CALL load_ptr_3d_sp( pt30, cdna30, psgn30, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      IF( nn_comm == 1 ) THEN 
         CALL lbc_lnk_pt2pt(   cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ELSE
         CALL lbc_lnk_neicoll( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ENDIF
      !
   END SUBROUTINE lbc_lnk_call_3d_sp


   SUBROUTINE load_ptr_3d_sp( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(sp), DIMENSION(:,:,:), TARGET, INTENT(inout), CONTIGUOUS ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(sp)               , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_4d_sp), DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(sp) , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt4d(1:SIZE(ptab, dim=1),1:SIZE(ptab, dim=2),1:SIZE(ptab, dim=3),1:1) => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_3d_sp


   SUBROUTINE lbc_lnk_call_4d_sp(                                                                &
      &                     cdname                                                                                  &
      &                   , pt1 , cdna1 , psgn1 , pt2 , cdna2 , psgn2 , pt3 , cdna3 , psgn3 , pt4 , cdna4 , psgn4   &
      &                   , pt5 , cdna5 , psgn5 , pt6 , cdna6 , psgn6 , pt7 , cdna7 , psgn7 , pt8 , cdna8 , psgn8   &
      &                   , pt9 , cdna9 , psgn9 , pt10, cdna10, psgn10, pt11, cdna11, psgn11, pt12, cdna12, psgn12  &
      &                   , pt13, cdna13, psgn13, pt14, cdna14, psgn14, pt15, cdna15, psgn15, pt16, cdna16, psgn16  &
      &                   , pt17, cdna17, psgn17, pt18, cdna18, psgn18, pt19, cdna19, psgn19, pt20, cdna20, psgn20  &
      &                   , pt21, cdna21, psgn21, pt22, cdna22, psgn22, pt23, cdna23, psgn23, pt24, cdna24, psgn24  &
      &                   , pt25, cdna25, psgn25, pt26, cdna26, psgn26, pt27, cdna27, psgn27, pt28, cdna28, psgn28  &
      &                   , pt29, cdna29, psgn29, pt30, cdna30, psgn30                                              &
      &                   , kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      !!---------------------------------------------------------------------
      CHARACTER(len=*)     ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(sp), DIMENSION(:,:,:,:)          , TARGET, CONTIGUOUS, INTENT(inout) ::   pt1        ! arrays on which the lbc is applied
      REAL(sp), DIMENSION(:,:,:,:), OPTIONAL, TARGET, CONTIGUOUS, INTENT(inout) ::   pt2 , pt3 , pt4 , pt5 , pt6 , pt7 , pt8 , &
         &                                                                               pt9 , pt10, pt11, pt12, pt13, pt14, pt15, &
         &                                                                               pt16, pt17, pt18, pt19, pt20, pt21, pt22, &
         &                                                                               pt23, pt24, pt25, pt26, pt27, pt28, pt29, &
         &                                                                               pt30
      CHARACTER(len=1)                       , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)     , OPTIONAL        , INTENT(in   ) ::   cdna2 , cdna3 , cdna4 , cdna5 , cdna6 , cdna7 , cdna8 , &
         &                                                        cdna9 , cdna10, cdna11, cdna12, cdna13, cdna14, cdna15, &
         &                                                        cdna16, cdna17, cdna18, cdna19, cdna20, cdna21, cdna22, &
         &                                                        cdna23, cdna24, cdna25, cdna26, cdna27, cdna28, cdna29, &
         &                                                        cdna30
      REAL(sp)                        , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(sp)      , OPTIONAL        , INTENT(in   ) ::   psgn2 , psgn3 , psgn4 , psgn5 , psgn6 , psgn7 , psgn8 , &
         &                                                        psgn9 , psgn10, psgn11, psgn12, psgn13, psgn14, psgn15, &
         &                                                        psgn16, psgn17, psgn18, psgn19, psgn20, psgn21, psgn22, &
         &                                                        psgn23, psgn24, psgn25, psgn26, psgn27, psgn28, psgn29, &
         &                                                        psgn30
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(sp)      , OPTIONAL        , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8), OPTIONAL        , INTENT(in   ) ::   lsend, lrecv   ! indicate how communications are to be carried out
      LOGICAL              , OPTIONAL        , INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !!
      INTEGER                          ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_4d_sp), DIMENSION(30) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(30) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(sp)  , DIMENSION(30) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_4d_sp( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2 ) )   CALL load_ptr_4d_sp( pt2 , cdna2 , psgn2 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3 ) )   CALL load_ptr_4d_sp( pt3 , cdna3 , psgn3 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4 ) )   CALL load_ptr_4d_sp( pt4 , cdna4 , psgn4 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5 ) )   CALL load_ptr_4d_sp( pt5 , cdna5 , psgn5 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6 ) )   CALL load_ptr_4d_sp( pt6 , cdna6 , psgn6 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7 ) )   CALL load_ptr_4d_sp( pt7 , cdna7 , psgn7 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8 ) )   CALL load_ptr_4d_sp( pt8 , cdna8 , psgn8 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9 ) )   CALL load_ptr_4d_sp( pt9 , cdna9 , psgn9 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn10) )   CALL load_ptr_4d_sp( pt10, cdna10, psgn10, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn11) )   CALL load_ptr_4d_sp( pt11, cdna11, psgn11, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn12) )   CALL load_ptr_4d_sp( pt12, cdna12, psgn12, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn13) )   CALL load_ptr_4d_sp( pt13, cdna13, psgn13, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn14) )   CALL load_ptr_4d_sp( pt14, cdna14, psgn14, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn15) )   CALL load_ptr_4d_sp( pt15, cdna15, psgn15, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn16) )   CALL load_ptr_4d_sp( pt16, cdna16, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn17) )   CALL load_ptr_4d_sp( pt17, cdna17, psgn17, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn18) )   CALL load_ptr_4d_sp( pt18, cdna18, psgn18, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn19) )   CALL load_ptr_4d_sp( pt19, cdna19, psgn19, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn20) )   CALL load_ptr_4d_sp( pt20, cdna20, psgn20, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn21) )   CALL load_ptr_4d_sp( pt21, cdna21, psgn21, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn22) )   CALL load_ptr_4d_sp( pt22, cdna22, psgn22, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn23) )   CALL load_ptr_4d_sp( pt23, cdna23, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn24) )   CALL load_ptr_4d_sp( pt24, cdna24, psgn24, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn25) )   CALL load_ptr_4d_sp( pt25, cdna25, psgn25, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn26) )   CALL load_ptr_4d_sp( pt26, cdna26, psgn26, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn27) )   CALL load_ptr_4d_sp( pt27, cdna27, psgn27, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn28) )   CALL load_ptr_4d_sp( pt28, cdna28, psgn28, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn29) )   CALL load_ptr_4d_sp( pt29, cdna29, psgn29, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn30) )   CALL load_ptr_4d_sp( pt30, cdna30, psgn30, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      IF( nn_comm == 1 ) THEN 
         CALL lbc_lnk_pt2pt(   cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ELSE
         CALL lbc_lnk_neicoll( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ENDIF
      !
   END SUBROUTINE lbc_lnk_call_4d_sp


   SUBROUTINE load_ptr_4d_sp( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(sp), DIMENSION(:,:,:,:), TARGET, INTENT(inout), CONTIGUOUS ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(sp)               , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_4d_sp), DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(sp) , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt4d(1:SIZE(ptab, dim=1),1:SIZE(ptab, dim=2),1:SIZE(ptab, dim=3),1:SIZE(ptab, dim=4)) => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_4d_sp

   !!
   !!   ----   DOUBLE PRECISION VERSIONS
   !!

   SUBROUTINE lbc_lnk_call_2d_dp(                                                                &
      &                     cdname                                                                                  &
      &                   , pt1 , cdna1 , psgn1 , pt2 , cdna2 , psgn2 , pt3 , cdna3 , psgn3 , pt4 , cdna4 , psgn4   &
      &                   , pt5 , cdna5 , psgn5 , pt6 , cdna6 , psgn6 , pt7 , cdna7 , psgn7 , pt8 , cdna8 , psgn8   &
      &                   , pt9 , cdna9 , psgn9 , pt10, cdna10, psgn10, pt11, cdna11, psgn11, pt12, cdna12, psgn12  &
      &                   , pt13, cdna13, psgn13, pt14, cdna14, psgn14, pt15, cdna15, psgn15, pt16, cdna16, psgn16  &
      &                   , pt17, cdna17, psgn17, pt18, cdna18, psgn18, pt19, cdna19, psgn19, pt20, cdna20, psgn20  &
      &                   , pt21, cdna21, psgn21, pt22, cdna22, psgn22, pt23, cdna23, psgn23, pt24, cdna24, psgn24  &
      &                   , pt25, cdna25, psgn25, pt26, cdna26, psgn26, pt27, cdna27, psgn27, pt28, cdna28, psgn28  &
      &                   , pt29, cdna29, psgn29, pt30, cdna30, psgn30                                              &
      &                   , kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      !!---------------------------------------------------------------------
      CHARACTER(len=*)     ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(dp), DIMENSION(:,:)          , TARGET, CONTIGUOUS, INTENT(inout) ::   pt1        ! arrays on which the lbc is applied
      REAL(dp), DIMENSION(:,:), OPTIONAL, TARGET, CONTIGUOUS, INTENT(inout) ::   pt2 , pt3 , pt4 , pt5 , pt6 , pt7 , pt8 , &
         &                                                                               pt9 , pt10, pt11, pt12, pt13, pt14, pt15, &
         &                                                                               pt16, pt17, pt18, pt19, pt20, pt21, pt22, &
         &                                                                               pt23, pt24, pt25, pt26, pt27, pt28, pt29, &
         &                                                                               pt30
      CHARACTER(len=1)                       , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)     , OPTIONAL        , INTENT(in   ) ::   cdna2 , cdna3 , cdna4 , cdna5 , cdna6 , cdna7 , cdna8 , &
         &                                                        cdna9 , cdna10, cdna11, cdna12, cdna13, cdna14, cdna15, &
         &                                                        cdna16, cdna17, cdna18, cdna19, cdna20, cdna21, cdna22, &
         &                                                        cdna23, cdna24, cdna25, cdna26, cdna27, cdna28, cdna29, &
         &                                                        cdna30
      REAL(dp)                        , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(dp)      , OPTIONAL        , INTENT(in   ) ::   psgn2 , psgn3 , psgn4 , psgn5 , psgn6 , psgn7 , psgn8 , &
         &                                                        psgn9 , psgn10, psgn11, psgn12, psgn13, psgn14, psgn15, &
         &                                                        psgn16, psgn17, psgn18, psgn19, psgn20, psgn21, psgn22, &
         &                                                        psgn23, psgn24, psgn25, psgn26, psgn27, psgn28, psgn29, &
         &                                                        psgn30
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(dp)      , OPTIONAL        , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8), OPTIONAL        , INTENT(in   ) ::   lsend, lrecv   ! indicate how communications are to be carried out
      LOGICAL              , OPTIONAL        , INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !!
      INTEGER                          ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_4d_dp), DIMENSION(30) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(30) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(dp)  , DIMENSION(30) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_2d_dp( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2 ) )   CALL load_ptr_2d_dp( pt2 , cdna2 , psgn2 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3 ) )   CALL load_ptr_2d_dp( pt3 , cdna3 , psgn3 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4 ) )   CALL load_ptr_2d_dp( pt4 , cdna4 , psgn4 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5 ) )   CALL load_ptr_2d_dp( pt5 , cdna5 , psgn5 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6 ) )   CALL load_ptr_2d_dp( pt6 , cdna6 , psgn6 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7 ) )   CALL load_ptr_2d_dp( pt7 , cdna7 , psgn7 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8 ) )   CALL load_ptr_2d_dp( pt8 , cdna8 , psgn8 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9 ) )   CALL load_ptr_2d_dp( pt9 , cdna9 , psgn9 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn10) )   CALL load_ptr_2d_dp( pt10, cdna10, psgn10, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn11) )   CALL load_ptr_2d_dp( pt11, cdna11, psgn11, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn12) )   CALL load_ptr_2d_dp( pt12, cdna12, psgn12, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn13) )   CALL load_ptr_2d_dp( pt13, cdna13, psgn13, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn14) )   CALL load_ptr_2d_dp( pt14, cdna14, psgn14, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn15) )   CALL load_ptr_2d_dp( pt15, cdna15, psgn15, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn16) )   CALL load_ptr_2d_dp( pt16, cdna16, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn17) )   CALL load_ptr_2d_dp( pt17, cdna17, psgn17, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn18) )   CALL load_ptr_2d_dp( pt18, cdna18, psgn18, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn19) )   CALL load_ptr_2d_dp( pt19, cdna19, psgn19, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn20) )   CALL load_ptr_2d_dp( pt20, cdna20, psgn20, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn21) )   CALL load_ptr_2d_dp( pt21, cdna21, psgn21, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn22) )   CALL load_ptr_2d_dp( pt22, cdna22, psgn22, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn23) )   CALL load_ptr_2d_dp( pt23, cdna23, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn24) )   CALL load_ptr_2d_dp( pt24, cdna24, psgn24, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn25) )   CALL load_ptr_2d_dp( pt25, cdna25, psgn25, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn26) )   CALL load_ptr_2d_dp( pt26, cdna26, psgn26, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn27) )   CALL load_ptr_2d_dp( pt27, cdna27, psgn27, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn28) )   CALL load_ptr_2d_dp( pt28, cdna28, psgn28, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn29) )   CALL load_ptr_2d_dp( pt29, cdna29, psgn29, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn30) )   CALL load_ptr_2d_dp( pt30, cdna30, psgn30, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      IF( nn_comm == 1 ) THEN 
         CALL lbc_lnk_pt2pt(   cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ELSE
         CALL lbc_lnk_neicoll( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ENDIF
      !
   END SUBROUTINE lbc_lnk_call_2d_dp


   SUBROUTINE load_ptr_2d_dp( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(dp), DIMENSION(:,:), TARGET, INTENT(inout), CONTIGUOUS ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(dp)               , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_4d_dp), DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(dp) , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt4d(1:SIZE(ptab, dim=1),1:SIZE(ptab, dim=2),1:1,1:1) => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_2d_dp


   SUBROUTINE lbc_lnk_call_3d_dp(                                                                &
      &                     cdname                                                                                  &
      &                   , pt1 , cdna1 , psgn1 , pt2 , cdna2 , psgn2 , pt3 , cdna3 , psgn3 , pt4 , cdna4 , psgn4   &
      &                   , pt5 , cdna5 , psgn5 , pt6 , cdna6 , psgn6 , pt7 , cdna7 , psgn7 , pt8 , cdna8 , psgn8   &
      &                   , pt9 , cdna9 , psgn9 , pt10, cdna10, psgn10, pt11, cdna11, psgn11, pt12, cdna12, psgn12  &
      &                   , pt13, cdna13, psgn13, pt14, cdna14, psgn14, pt15, cdna15, psgn15, pt16, cdna16, psgn16  &
      &                   , pt17, cdna17, psgn17, pt18, cdna18, psgn18, pt19, cdna19, psgn19, pt20, cdna20, psgn20  &
      &                   , pt21, cdna21, psgn21, pt22, cdna22, psgn22, pt23, cdna23, psgn23, pt24, cdna24, psgn24  &
      &                   , pt25, cdna25, psgn25, pt26, cdna26, psgn26, pt27, cdna27, psgn27, pt28, cdna28, psgn28  &
      &                   , pt29, cdna29, psgn29, pt30, cdna30, psgn30                                              &
      &                   , kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      !!---------------------------------------------------------------------
      CHARACTER(len=*)     ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(dp), DIMENSION(:,:,:)          , TARGET, CONTIGUOUS, INTENT(inout) ::   pt1        ! arrays on which the lbc is applied
      REAL(dp), DIMENSION(:,:,:), OPTIONAL, TARGET, CONTIGUOUS, INTENT(inout) ::   pt2 , pt3 , pt4 , pt5 , pt6 , pt7 , pt8 , &
         &                                                                               pt9 , pt10, pt11, pt12, pt13, pt14, pt15, &
         &                                                                               pt16, pt17, pt18, pt19, pt20, pt21, pt22, &
         &                                                                               pt23, pt24, pt25, pt26, pt27, pt28, pt29, &
         &                                                                               pt30
      CHARACTER(len=1)                       , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)     , OPTIONAL        , INTENT(in   ) ::   cdna2 , cdna3 , cdna4 , cdna5 , cdna6 , cdna7 , cdna8 , &
         &                                                        cdna9 , cdna10, cdna11, cdna12, cdna13, cdna14, cdna15, &
         &                                                        cdna16, cdna17, cdna18, cdna19, cdna20, cdna21, cdna22, &
         &                                                        cdna23, cdna24, cdna25, cdna26, cdna27, cdna28, cdna29, &
         &                                                        cdna30
      REAL(dp)                        , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(dp)      , OPTIONAL        , INTENT(in   ) ::   psgn2 , psgn3 , psgn4 , psgn5 , psgn6 , psgn7 , psgn8 , &
         &                                                        psgn9 , psgn10, psgn11, psgn12, psgn13, psgn14, psgn15, &
         &                                                        psgn16, psgn17, psgn18, psgn19, psgn20, psgn21, psgn22, &
         &                                                        psgn23, psgn24, psgn25, psgn26, psgn27, psgn28, psgn29, &
         &                                                        psgn30
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(dp)      , OPTIONAL        , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8), OPTIONAL        , INTENT(in   ) ::   lsend, lrecv   ! indicate how communications are to be carried out
      LOGICAL              , OPTIONAL        , INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !!
      INTEGER                          ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_4d_dp), DIMENSION(30) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(30) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(dp)  , DIMENSION(30) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_3d_dp( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2 ) )   CALL load_ptr_3d_dp( pt2 , cdna2 , psgn2 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3 ) )   CALL load_ptr_3d_dp( pt3 , cdna3 , psgn3 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4 ) )   CALL load_ptr_3d_dp( pt4 , cdna4 , psgn4 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5 ) )   CALL load_ptr_3d_dp( pt5 , cdna5 , psgn5 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6 ) )   CALL load_ptr_3d_dp( pt6 , cdna6 , psgn6 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7 ) )   CALL load_ptr_3d_dp( pt7 , cdna7 , psgn7 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8 ) )   CALL load_ptr_3d_dp( pt8 , cdna8 , psgn8 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9 ) )   CALL load_ptr_3d_dp( pt9 , cdna9 , psgn9 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn10) )   CALL load_ptr_3d_dp( pt10, cdna10, psgn10, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn11) )   CALL load_ptr_3d_dp( pt11, cdna11, psgn11, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn12) )   CALL load_ptr_3d_dp( pt12, cdna12, psgn12, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn13) )   CALL load_ptr_3d_dp( pt13, cdna13, psgn13, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn14) )   CALL load_ptr_3d_dp( pt14, cdna14, psgn14, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn15) )   CALL load_ptr_3d_dp( pt15, cdna15, psgn15, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn16) )   CALL load_ptr_3d_dp( pt16, cdna16, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn17) )   CALL load_ptr_3d_dp( pt17, cdna17, psgn17, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn18) )   CALL load_ptr_3d_dp( pt18, cdna18, psgn18, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn19) )   CALL load_ptr_3d_dp( pt19, cdna19, psgn19, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn20) )   CALL load_ptr_3d_dp( pt20, cdna20, psgn20, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn21) )   CALL load_ptr_3d_dp( pt21, cdna21, psgn21, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn22) )   CALL load_ptr_3d_dp( pt22, cdna22, psgn22, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn23) )   CALL load_ptr_3d_dp( pt23, cdna23, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn24) )   CALL load_ptr_3d_dp( pt24, cdna24, psgn24, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn25) )   CALL load_ptr_3d_dp( pt25, cdna25, psgn25, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn26) )   CALL load_ptr_3d_dp( pt26, cdna26, psgn26, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn27) )   CALL load_ptr_3d_dp( pt27, cdna27, psgn27, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn28) )   CALL load_ptr_3d_dp( pt28, cdna28, psgn28, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn29) )   CALL load_ptr_3d_dp( pt29, cdna29, psgn29, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn30) )   CALL load_ptr_3d_dp( pt30, cdna30, psgn30, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      IF( nn_comm == 1 ) THEN 
         CALL lbc_lnk_pt2pt(   cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ELSE
         CALL lbc_lnk_neicoll( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ENDIF
      !
   END SUBROUTINE lbc_lnk_call_3d_dp


   SUBROUTINE load_ptr_3d_dp( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(dp), DIMENSION(:,:,:), TARGET, INTENT(inout), CONTIGUOUS ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(dp)               , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_4d_dp), DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(dp) , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt4d(1:SIZE(ptab, dim=1),1:SIZE(ptab, dim=2),1:SIZE(ptab, dim=3),1:1) => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_3d_dp


   SUBROUTINE lbc_lnk_call_4d_dp(                                                                &
      &                     cdname                                                                                  &
      &                   , pt1 , cdna1 , psgn1 , pt2 , cdna2 , psgn2 , pt3 , cdna3 , psgn3 , pt4 , cdna4 , psgn4   &
      &                   , pt5 , cdna5 , psgn5 , pt6 , cdna6 , psgn6 , pt7 , cdna7 , psgn7 , pt8 , cdna8 , psgn8   &
      &                   , pt9 , cdna9 , psgn9 , pt10, cdna10, psgn10, pt11, cdna11, psgn11, pt12, cdna12, psgn12  &
      &                   , pt13, cdna13, psgn13, pt14, cdna14, psgn14, pt15, cdna15, psgn15, pt16, cdna16, psgn16  &
      &                   , pt17, cdna17, psgn17, pt18, cdna18, psgn18, pt19, cdna19, psgn19, pt20, cdna20, psgn20  &
      &                   , pt21, cdna21, psgn21, pt22, cdna22, psgn22, pt23, cdna23, psgn23, pt24, cdna24, psgn24  &
      &                   , pt25, cdna25, psgn25, pt26, cdna26, psgn26, pt27, cdna27, psgn27, pt28, cdna28, psgn28  &
      &                   , pt29, cdna29, psgn29, pt30, cdna30, psgn30                                              &
      &                   , kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      !!---------------------------------------------------------------------
      CHARACTER(len=*)     ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(dp), DIMENSION(:,:,:,:)          , TARGET, CONTIGUOUS, INTENT(inout) ::   pt1        ! arrays on which the lbc is applied
      REAL(dp), DIMENSION(:,:,:,:), OPTIONAL, TARGET, CONTIGUOUS, INTENT(inout) ::   pt2 , pt3 , pt4 , pt5 , pt6 , pt7 , pt8 , &
         &                                                                               pt9 , pt10, pt11, pt12, pt13, pt14, pt15, &
         &                                                                               pt16, pt17, pt18, pt19, pt20, pt21, pt22, &
         &                                                                               pt23, pt24, pt25, pt26, pt27, pt28, pt29, &
         &                                                                               pt30
      CHARACTER(len=1)                       , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)     , OPTIONAL        , INTENT(in   ) ::   cdna2 , cdna3 , cdna4 , cdna5 , cdna6 , cdna7 , cdna8 , &
         &                                                        cdna9 , cdna10, cdna11, cdna12, cdna13, cdna14, cdna15, &
         &                                                        cdna16, cdna17, cdna18, cdna19, cdna20, cdna21, cdna22, &
         &                                                        cdna23, cdna24, cdna25, cdna26, cdna27, cdna28, cdna29, &
         &                                                        cdna30
      REAL(dp)                        , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(dp)      , OPTIONAL        , INTENT(in   ) ::   psgn2 , psgn3 , psgn4 , psgn5 , psgn6 , psgn7 , psgn8 , &
         &                                                        psgn9 , psgn10, psgn11, psgn12, psgn13, psgn14, psgn15, &
         &                                                        psgn16, psgn17, psgn18, psgn19, psgn20, psgn21, psgn22, &
         &                                                        psgn23, psgn24, psgn25, psgn26, psgn27, psgn28, psgn29, &
         &                                                        psgn30
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(dp)      , OPTIONAL        , INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER              , OPTIONAL        , INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8), OPTIONAL        , INTENT(in   ) ::   lsend, lrecv   ! indicate how communications are to be carried out
      LOGICAL              , OPTIONAL        , INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !!
      INTEGER                          ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_4d_dp), DIMENSION(30) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(30) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(dp)  , DIMENSION(30) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_4d_dp( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2 ) )   CALL load_ptr_4d_dp( pt2 , cdna2 , psgn2 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3 ) )   CALL load_ptr_4d_dp( pt3 , cdna3 , psgn3 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4 ) )   CALL load_ptr_4d_dp( pt4 , cdna4 , psgn4 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5 ) )   CALL load_ptr_4d_dp( pt5 , cdna5 , psgn5 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6 ) )   CALL load_ptr_4d_dp( pt6 , cdna6 , psgn6 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7 ) )   CALL load_ptr_4d_dp( pt7 , cdna7 , psgn7 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8 ) )   CALL load_ptr_4d_dp( pt8 , cdna8 , psgn8 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9 ) )   CALL load_ptr_4d_dp( pt9 , cdna9 , psgn9 , ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn10) )   CALL load_ptr_4d_dp( pt10, cdna10, psgn10, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn11) )   CALL load_ptr_4d_dp( pt11, cdna11, psgn11, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn12) )   CALL load_ptr_4d_dp( pt12, cdna12, psgn12, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn13) )   CALL load_ptr_4d_dp( pt13, cdna13, psgn13, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn14) )   CALL load_ptr_4d_dp( pt14, cdna14, psgn14, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn15) )   CALL load_ptr_4d_dp( pt15, cdna15, psgn15, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn16) )   CALL load_ptr_4d_dp( pt16, cdna16, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn17) )   CALL load_ptr_4d_dp( pt17, cdna17, psgn17, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn18) )   CALL load_ptr_4d_dp( pt18, cdna18, psgn18, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn19) )   CALL load_ptr_4d_dp( pt19, cdna19, psgn19, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn20) )   CALL load_ptr_4d_dp( pt20, cdna20, psgn20, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn21) )   CALL load_ptr_4d_dp( pt21, cdna21, psgn21, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn22) )   CALL load_ptr_4d_dp( pt22, cdna22, psgn22, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn23) )   CALL load_ptr_4d_dp( pt23, cdna23, psgn16, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn24) )   CALL load_ptr_4d_dp( pt24, cdna24, psgn24, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn25) )   CALL load_ptr_4d_dp( pt25, cdna25, psgn25, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn26) )   CALL load_ptr_4d_dp( pt26, cdna26, psgn26, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn27) )   CALL load_ptr_4d_dp( pt27, cdna27, psgn27, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn28) )   CALL load_ptr_4d_dp( pt28, cdna28, psgn28, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn29) )   CALL load_ptr_4d_dp( pt29, cdna29, psgn29, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn30) )   CALL load_ptr_4d_dp( pt30, cdna30, psgn30, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      IF( nn_comm == 1 ) THEN 
         CALL lbc_lnk_pt2pt(   cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ELSE
         CALL lbc_lnk_neicoll( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      ENDIF
      !
   END SUBROUTINE lbc_lnk_call_4d_dp


   SUBROUTINE load_ptr_4d_dp( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(dp), DIMENSION(:,:,:,:), TARGET, INTENT(inout), CONTIGUOUS ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(dp)               , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_4d_dp), DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(dp) , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt4d(1:SIZE(ptab, dim=1),1:SIZE(ptab, dim=2),1:SIZE(ptab, dim=3),1:SIZE(ptab, dim=4)) => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_4d_dp

   !
   !!----------------------------------------------------------------------
   !!                   ***  lbc_lnk_pt2pt_[sd]p  ***
   !!                  ***  lbc_lnk_neicoll_[sd]p  ***
   !!
   !!   * Argument : dummy argument use in lbc_lnk_... routines
   !!                cdname    :   name of the calling subroutine (for monitoring)
   !!                ptab      :   pointer of arrays on which the boundary condition is applied
   !!                cd_nat    :   nature of array grid-points
   !!                psgn      :   sign used across the north fold boundary
   !!                kfld      :   number of pt3d arrays
   !!                kfillmode :   optional, method to be use to fill the halos (see jpfill* variables)
   !!                pfillval  :   optional, background value (used with jpfillcopy)
   !!----------------------------------------------------------------------
   !!
   !!   ----   SINGLE PRECISION VERSIONS
   !!
   SUBROUTINE lbc_lnk_pt2pt_sp( cdname, ptab, cd_nat, psgn, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      CHARACTER(len=*)              , INTENT(in   ) ::   cdname      ! name of the calling subroutine
      TYPE(PTR_4d_sp),  DIMENSION(:), INTENT(inout) ::   ptab        ! pointer of arrays on which apply the b.c.
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(sp),  DIMENSION(:), INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER                       , INTENT(in   ) ::   kfld        ! number of pt3d arrays
      INTEGER ,             OPTIONAL, INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(sp),      OPTIONAL, INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER ,             OPTIONAL, INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8),OPTIONAL, INTENT(in   ) ::   lsend, lrecv  ! communication with other 4 proc
      LOGICAL,              OPTIONAL, INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !
      INTEGER  ::     ji,   jj,  jk,  jl,  jf, jn     ! dummy loop indices
      INTEGER  ::    ipi,  ipj, ipk, ipl, ipf         ! dimension of the input array
      INTEGER  ::   ip0i, ip1i, im0i, im1i
      INTEGER  ::   ip0j, ip1j, im0j, im1j
      INTEGER  ::   ishti, ishtj, ishti2, ishtj2
      INTEGER  ::   ifill_nfd, icomm, ierr
      INTEGER  ::   ihls, idxs, idxr, iszS, iszR
      INTEGER, DIMENSION(4)  ::   iwewe, issnn
      INTEGER, DIMENSION(8)  ::   isizei, ishtSi, ishtRi, ishtPi
      INTEGER, DIMENSION(8)  ::   isizej, ishtSj, ishtRj, ishtPj
      INTEGER, DIMENSION(8)  ::   ifill, iszall, ishtS, ishtR
      INTEGER, DIMENSION(8)  ::   ireq             ! mpi_request id
      INTEGER, DIMENSION(8)  ::   iStag, iRtag     ! Send and Recv mpi_tag id
      REAL(sp) ::   zland
      LOGICAL, DIMENSION(8)  ::   llsend, llrecv
      LOGICAL  ::   ll4only                                        ! default: 8 neighbourgs
      !!----------------------------------------------------------------------
      !
      ! ----------------------------------------- !
      !     1. local variables initialization     !
      ! ----------------------------------------- !
      !
      ipi = SIZE(ptab(1)%pt4d,1)
      ipj = SIZE(ptab(1)%pt4d,2)
      ipk = SIZE(ptab(1)%pt4d,3)
      ipl = SIZE(ptab(1)%pt4d,4)
      ipf = kfld
      !
      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, ipk, ipl, ipf, ld_lbc = .TRUE. )
      !
      idxs = 1   ! initalize index for send buffer
      idxr = 1   ! initalize index for recv buffer
      icomm = mpi_comm_oce        ! shorter name
      !
      ! take care of optional parameters
      !
      ihls = nn_hls   ! default definition
      IF( PRESENT( khls ) )   ihls = khls
      IF( ihls > n_hlsmax ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with khls > n_hlsmax : ', khls, '>', n_hlsmax
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      IF( ipi /= Ni_0+2*ihls ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with an input array which does not match ihls along i: ', ipi, ihls, Ni_0
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      IF( ipj /= Nj_0+2*ihls ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with an input array which does not match ihls along j:', ipj, ihls , Nj_0
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      !
      ll4only = .FALSE.    ! default definition
      IF( PRESENT(ld4only) )   ll4only = ld4only
      !
      zland = 0._wp                                     ! land filling value: zero by default
      IF( PRESENT( pfillval ) )   zland = pfillval      ! set land value
      !
      ! define llsend and llrecv: logicals which say if mpi-neibourgs for send or receive exist or not.
      IF     ( PRESENT(lsend) .AND. PRESENT(lrecv) ) THEN   ! localy defined neighbourgs 
         llsend(:) = lsend(:)   ;   llrecv(:) = lrecv(:)
      ELSE IF( PRESENT(lsend) .OR.  PRESENT(lrecv) ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with only one of the two arguments lsend or lrecv'
         CALL ctl_stop( 'STOP', ctmp1 )
      ELSE                                              ! default neighbours
         llsend(:) = mpiSnei(ihls,:) >= 0
         IF( ll4only )   llsend(5:8) = .FALSE.          ! exclude corners
         llrecv(:) = mpiRnei(ihls,:) >= 0
         IF( ll4only )   llrecv(5:8) = .FALSE.          ! exclude corners
      ENDIF
      !
      ! define ifill: which method should be used to fill each parts (sides+corners) of the halos
      ! default definition
      DO jn = 1, 4
         IF(             llrecv(jn) ) THEN   ;   ifill(jn) = jpfillmpi    ! with an mpi communication
         ELSEIF(    l_SelfPerio(jn) ) THEN   ;   ifill(jn) = jpfillperio  ! with self-periodicity
         ELSEIF( PRESENT(kfillmode) ) THEN   ;   ifill(jn) = kfillmode    ! localy defined
         ELSE                                ;   ifill(jn) = jpfillcst    ! constant value (zland)
         ENDIF
      END DO
      DO jn = 5, 8
         IF(             llrecv(jn) ) THEN   ;   ifill(jn) = jpfillmpi    ! with an mpi communication
         ELSE                                ;   ifill(jn) = jpfillnothing! do nothing
         ENDIF
      END DO
         !
      ! north fold treatment
      IF( l_IdoNFold ) THEN
         ifill_nfd = ifill(jpno)             ! if we are here, this means llrecv(jpno) = .false. and l_SelfPerio(jpno) = .false.
         ifill( (/jpno/) ) = jpfillnothing   ! we do north fold -> do nothing for northern halo
      ENDIF
      
      ! We first define the localization and size of the parts of the array that will be sent (s), received (r)
      ! or used for periodocity (p). The localization is defined as "the bottom left corner - 1" in i and j directions.
      ! This is a shift that will be applied later in the do loops to pick-up the appropriate part of the array
      !
      ! all definitions bellow do not refer to N[ij][se]0 so we can use it with any local value of ihls
      !                   !                       ________________________
      ip0i =          0   !          im0j = inner |__|__|__________|__|__|
      ip1i =       ihls   !   im1j = inner - halo |__|__|__________|__|__|
      im1i = ipi-2*ihls   !                       |  |  |          |  |  |
      im0i = ipi - ihls   !                       |  |  |          |  |  |
      ip0j =          0   !                       |  |  |          |  |  |
      ip1j =       ihls   !                       |__|__|__________|__|__|
      im1j = ipj-2*ihls   !           ip1j = halo |__|__|__________|__|__|
      im0j = ipj - ihls   !              ip0j = 0 |__|__|__________|__|__|
      !                   !                    ip0i ip1i        im1i im0i
      !
      iwewe(:) = (/ jpwe,jpea,jpwe,jpea /)   ;   issnn(:) = (/ jpso,jpso,jpno,jpno /)
      !     sides:     west  east south north      ;   corners: so-we, so-ea, no-we, no-ea
      isizei(1:4) = (/ ihls, ihls,  ipi,  ipi /)   ;   isizei(5:8) = ihls              ! i- count
      isizej(1:4) = (/ Nj_0, Nj_0, ihls, ihls /)   ;   isizej(5:8) = ihls              ! j- count
      ishtSi(1:4) = (/ ip1i, im1i, ip0i, ip0i /)   ;   ishtSi(5:8) = ishtSi( iwewe )   ! i- shift send data
      ishtSj(1:4) = (/ ip1j, ip1j, ip1j, im1j /)   ;   ishtSj(5:8) = ishtSj( issnn )   ! j- shift send data
      ishtRi(1:4) = (/ ip0i, im0i, ip0i, ip0i /)   ;   ishtRi(5:8) = ishtRi( iwewe )   ! i- shift received data location
      ishtRj(1:4) = (/ ip1j, ip1j, ip0j, im0j /)   ;   ishtRj(5:8) = ishtRj( issnn )   ! j- shift received data location
      ishtPi(1:4) = (/ im1i, ip1i, ip0i, ip0i /)   ;   ishtPi(5:8) = ishtPi( iwewe )   ! i- shift data used for periodicity
      ishtPj(1:4) = (/ ip1j, ip1j, im1j, ip1j /)   ;   ishtPj(5:8) = ishtPj( issnn )   ! j- shift data used for periodicity
      !
      ! -------------------------------- !
      !     2. Prepare MPI exchanges     !
      ! -------------------------------- !
      !
      iStag = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)   ! any value but each one must be different
      ! define iRtag with the corresponding iStag, e.g. data received at west where sent at east.
      iRtag(jpwe) = iStag(jpea)   ;   iRtag(jpea) = iStag(jpwe)   ;   iRtag(jpso) = iStag(jpno)   ;   iRtag(jpno) = iStag(jpso)
      iRtag(jpsw) = iStag(jpne)   ;   iRtag(jpse) = iStag(jpnw)   ;   iRtag(jpnw) = iStag(jpse)   ;   iRtag(jpne) = iStag(jpsw)
      !
      iszall(:) = isizei(:) * isizej(:) * ipk * ipl * ipf
      ishtS(1) = 0
      DO jn = 2, 8
         ishtS(jn) = ishtS(jn-1) + iszall(jn-1) * COUNT( (/llsend(jn-1)/) )
      END DO
      ishtR(1) = 0
      DO jn = 2, 8
         ishtR(jn) = ishtR(jn-1) + iszall(jn-1) * COUNT( (/llrecv(jn-1)/) )
      END DO

      ! Allocate buffer arrays to be sent/received if needed
      iszS = SUM(iszall, mask = llsend)                             ! send buffer size
      IF( ALLOCATED(buffsnd_sp) ) THEN
         CALL mpi_waitall(8, nreq_p2p, MPI_STATUSES_IGNORE, ierr)   ! wait for Isend from the PREVIOUS call
         IF( SIZE(buffsnd_sp) < iszS )    DEALLOCATE(buffsnd_sp)          ! send buffer is too small
      ENDIF
      IF( .NOT. ALLOCATED(buffsnd_sp) )   ALLOCATE( buffsnd_sp(iszS) )
      iszR = SUM(iszall, mask = llrecv)                             ! recv buffer size
      IF( ALLOCATED(buffrcv_sp) ) THEN
         IF( SIZE(buffrcv_sp) < iszR )    DEALLOCATE(buffrcv_sp)          ! recv buffer is too small
      ENDIF
      IF( .NOT. ALLOCATED(buffrcv_sp) )   ALLOCATE( buffrcv_sp(iszR) )
      !
      ! default definition when no communication is done. understood by mpi_waitall
      nreq_p2p(:) = MPI_REQUEST_NULL   ! WARNING: Must be done after the call to mpi_waitall just above
      !
      ! ----------------------------------------------- !
      !     3. Do east and west MPI_Isend if needed     !
      ! ----------------------------------------------- !
      !
      DO jn = 1, 2

   IF( llsend(jn) ) THEN
      ishti = ishtSi(jn)
      ishtj = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         buffsnd_sp(idxs) = ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl)
         idxs = idxs + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      ! non-blocking send of the west/east side using local buffer
      CALL MPI_ISEND( buffsnd_sp(ishtS(jn)+1), iszall(jn), MPI_REAL, mpiSnei(ihls,jn), iStag(jn), icomm, nreq_p2p(jn), ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
   ENDIF

      END DO
      !
      ! ----------------------------------- !
      !     4. Fill east and west halos     !
      ! ----------------------------------- !
      !
      DO jn = 1, 2


   ishti = ishtRi(jn)
   ishtj = ishtRj(jn)
   SELECT CASE ( ifill(jn) )
   CASE ( jpfillnothing )               ! no filling 
   CASE ( jpfillmpi   )                 ! fill with data received by MPI
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !                                 ! blocking receive of the west/east halo in local temporary arrays
      CALL MPI_RECV( buffrcv_sp(ishtR(jn)+1), iszall(jn), MPI_REAL, mpiRnei(ihls,jn), iRtag(jn), icomm, MPI_STATUS_IGNORE, ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = buffrcv_sp(idxr)
         idxr = idxr + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillperio )                 ! use periodicity
      ishti2 = ishtPi(jn)
      ishtj2 = ishtPj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcopy  )                 ! filling with inner domain values
      ishti2 = ishtSi(jn)
      ishtj2 = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcst   )                 ! filling with constant value
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = zland
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   END SELECT
      END DO
      !
      ! ------------------------------------------------- !
      !     5. Do north and south MPI_Isend if needed     !
      ! ------------------------------------------------- !
      !
      DO jn = 3, 4

   IF( llsend(jn) ) THEN
      ishti = ishtSi(jn)
      ishtj = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         buffsnd_sp(idxs) = ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl)
         idxs = idxs + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      ! non-blocking send of the west/east side using local buffer
      CALL MPI_ISEND( buffsnd_sp(ishtS(jn)+1), iszall(jn), MPI_REAL, mpiSnei(ihls,jn), iStag(jn), icomm, nreq_p2p(jn), ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
   ENDIF

      END DO
      !
      ! ------------------------------- !
      !     6. north fold treatment     !
      ! ------------------------------- !
      !
      ! Must be done after receiving data from East/West neighbourgs (as it is coded in mpp_nfd, to be changed one day...)
      ! Do it after MPI_iSend to south/north neighbourgs so they won't wait (too much) to receive their data
      ! Do if before MPI_Recv from south/north neighbourgs so we have more time to receive data
      !
      IF( l_IdoNFold ) THEN
         IF( jpni == 1 )  THEN   ;   CALL lbc_nfd( ptab, cd_nat, psgn                  , ihls, ipf )   ! self NFold
         ELSE                    ;   CALL mpp_nfd( ptab, cd_nat, psgn, ifill_nfd, zland, ihls, ipf )   ! mpi  NFold
         ENDIF
      ENDIF
      !
      ! ------------------------------------- !
      !     7. Fill south and north halos     !
      ! ------------------------------------- !
      !
      DO jn = 3, 4


   ishti = ishtRi(jn)
   ishtj = ishtRj(jn)
   SELECT CASE ( ifill(jn) )
   CASE ( jpfillnothing )               ! no filling 
   CASE ( jpfillmpi   )                 ! fill with data received by MPI
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !                                 ! blocking receive of the west/east halo in local temporary arrays
      CALL MPI_RECV( buffrcv_sp(ishtR(jn)+1), iszall(jn), MPI_REAL, mpiRnei(ihls,jn), iRtag(jn), icomm, MPI_STATUS_IGNORE, ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = buffrcv_sp(idxr)
         idxr = idxr + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillperio )                 ! use periodicity
      ishti2 = ishtPi(jn)
      ishtj2 = ishtPj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcopy  )                 ! filling with inner domain values
      ishti2 = ishtSi(jn)
      ishtj2 = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcst   )                 ! filling with constant value
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = zland
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   END SELECT
      END DO
      !
      ! ----------------------------------------------- !
      !     8. Specific problem in corner treatment     !
      !              ( very rate case... )              !
      ! ----------------------------------------------- !
      !
      DO jn = 5, 8

   IF( llsend(jn) ) THEN
      ishti = ishtSi(jn)
      ishtj = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         buffsnd_sp(idxs) = ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl)
         idxs = idxs + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      ! non-blocking send of the west/east side using local buffer
      CALL MPI_ISEND( buffsnd_sp(ishtS(jn)+1), iszall(jn), MPI_REAL, mpiSnei(ihls,jn), iStag(jn), icomm, nreq_p2p(jn), ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
   ENDIF

      END DO
      DO jn = 5, 8


   ishti = ishtRi(jn)
   ishtj = ishtRj(jn)
   SELECT CASE ( ifill(jn) )
   CASE ( jpfillnothing )               ! no filling 
   CASE ( jpfillmpi   )                 ! fill with data received by MPI
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !                                 ! blocking receive of the west/east halo in local temporary arrays
      CALL MPI_RECV( buffrcv_sp(ishtR(jn)+1), iszall(jn), MPI_REAL, mpiRnei(ihls,jn), iRtag(jn), icomm, MPI_STATUS_IGNORE, ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = buffrcv_sp(idxr)
         idxr = idxr + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillperio )                 ! use periodicity
      ishti2 = ishtPi(jn)
      ishtj2 = ishtPj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcopy  )                 ! filling with inner domain values
      ishti2 = ishtSi(jn)
      ishtj2 = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcst   )                 ! filling with constant value
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = zland
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   END SELECT
      END DO
      !
      ! -------------------------------------------- !
      !     9. deallocate local temporary arrays     !
      !        if they areg larger than jpi*jpj      !  <- arbitrary max size...
      ! -------------------------------------------- !
      !
      IF( iszR > jpi*jpj )   DEALLOCATE(buffrcv_sp)                    ! blocking receive -> can directly deallocate
      IF( iszS > jpi*jpj ) THEN
         CALL mpi_waitall(8, nreq_p2p, MPI_STATUSES_IGNORE, ierr)   ! must wait before deallocate send buffer
         DEALLOCATE(buffsnd_sp)
      ENDIF
      !
   END SUBROUTINE lbc_lnk_pt2pt_sp


SUBROUTINE lbc_lnk_neicoll_sp( cdname, ptab, cd_nat, psgn, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      CHARACTER(len=*)              , INTENT(in   ) ::   cdname      ! name of the calling subroutine
      TYPE(PTR_4d_sp),  DIMENSION(:), INTENT(inout) ::   ptab        ! pointer of arrays on which apply the b.c.
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(sp),  DIMENSION(:), INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER                       , INTENT(in   ) ::   kfld        ! number of pt3d arrays
      INTEGER ,             OPTIONAL, INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(sp),      OPTIONAL, INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER ,             OPTIONAL, INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8),OPTIONAL, INTENT(in   ) ::   lsend, lrecv  ! communication with other 4 proc
      LOGICAL,              OPTIONAL, INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !
      INTEGER  ::    ji,  jj,  jk , jl,  jf, jn      ! dummy loop indices
      INTEGER  ::   ipi, ipj, ipk, ipl, ipf          ! dimension of the input array
      INTEGER  ::   ip0i, ip1i, im0i, im1i
      INTEGER  ::   ip0j, ip1j, im0j, im1j
      INTEGER  ::   ishti, ishtj, ishti2, ishtj2
      INTEGER  ::   iszS, iszR
      INTEGER  ::   ierr
      INTEGER  ::   ihls, idx
      INTEGER  ::   impi_nc
      INTEGER  ::   ifill_nfd
      INTEGER, DIMENSION(4)  ::   iwewe, issnn
      INTEGER, DIMENSION(8)  ::   isizei, ishtSi, ishtRi, ishtPi
      INTEGER, DIMENSION(8)  ::   isizej, ishtSj, ishtRj, ishtPj
      INTEGER, DIMENSION(8)  ::   ifill, iszall
      INTEGER, DIMENSION(8)  ::   jnf
      INTEGER, DIMENSION(:), ALLOCATABLE  ::   iScnt, iRcnt    ! number of elements to be sent/received
      INTEGER, DIMENSION(:), ALLOCATABLE  ::   iSdpl, iRdpl    ! displacement in halos arrays
      LOGICAL, DIMENSION(8)  ::   llsend, llrecv
      REAL(sp) ::   zland
      LOGICAL  ::   ll4only                                    ! default: 8 neighbourgs
      !!----------------------------------------------------------------------
      !
      ! ----------------------------------------- !
      !     1. local variables initialization     !
      ! ----------------------------------------- !
      !
      ipi = SIZE(ptab(1)%pt4d,1)
      ipj = SIZE(ptab(1)%pt4d,2)
      ipk = SIZE(ptab(1)%pt4d,3)
      ipl = SIZE(ptab(1)%pt4d,4)
      ipf = kfld
      !
      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, ipk, ipl, ipf, ld_lbc = .TRUE. )
      !
      ! take care of optional parameters
      !
      ihls = nn_hls       ! default definition
      IF( PRESENT( khls ) )   ihls = khls
      IF( ihls > n_hlsmax ) THEN
         WRITE(ctmp1,*) TRIM(cdname), '  is calling lbc_lnk with khls > n_hlsmax : ', khls, '>', n_hlsmax
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      IF( ipi /= Ni_0+2*ihls ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with an input array which does not match ihls along i: ', ipi, ihls, Ni_0
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      IF( ipj /= Nj_0+2*ihls ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with an input array which does not match ihls along j:', ipj, ihls , Nj_0
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      !
      ll4only = .FALSE.    ! default definition
      IF( PRESENT(ld4only) )   ll4only = ld4only
      !
      impi_nc = mpi_nc_com8(ihls)   ! default
      IF( ll4only )   impi_nc = mpi_nc_com4(ihls)
      !
      zland = 0._wp                                     ! land filling value: zero by default
      IF( PRESENT( pfillval ) )   zland = pfillval      ! set land value
      !
      ! define llsend and llrecv: logicals which say if mpi-neibourgs for send or receive exist or not.
      IF     ( PRESENT(lsend) .AND. PRESENT(lrecv) ) THEN   ! localy defined neighbourgs 
         CALL ctl_stop( 'STOP', 'mpp_nc_generic+lsend and lrecv not yet implemented')
      ELSE IF( PRESENT(lsend) .OR.  PRESENT(lrecv) ) THEN
         WRITE(ctmp1,*) TRIM(cdname), '  is calling lbc_lnk with only one of the two arguments lsend or lrecv'
         CALL ctl_stop( 'STOP', ctmp1 )
      ELSE                                              ! default neighbours
         llsend(:) = mpiSnei(ihls,:) >= 0
         IF( ll4only )   llsend(5:8) = .FALSE.          ! exclude corners
         llrecv(:) = mpiRnei(ihls,:) >= 0
         IF( ll4only )   llrecv(5:8) = .FALSE.          ! exclude corners
      ENDIF
      !
      ! define ifill: which method should be used to fill each parts (sides+corners) of the halos
      ! default definition
      DO jn = 1, 8
         IF(             llrecv(jn) ) THEN   ;   ifill(jn) = jpfillmpi    ! with an mpi communication
         ELSEIF(    l_SelfPerio(jn) ) THEN   ;   ifill(jn) = jpfillperio  ! with self-periodicity
         ELSEIF( PRESENT(kfillmode) ) THEN   ;   ifill(jn) = kfillmode    ! localy defined
         ELSE                                ;   ifill(jn) = jpfillcst    ! constant value (zland)
         ENDIF
      END DO
      ! take care of "indirect self-periodicity" for the corners
      DO jn = 5, 8
         IF(.NOT.l_SelfPerio(jn) .AND. l_SelfPerio(jpwe))   ifill(jn) = jpfillnothing   ! no bi-perio but ew-perio: do corners later
         IF(.NOT.l_SelfPerio(jn) .AND. l_SelfPerio(jpso))   ifill(jn) = jpfillnothing   ! no bi-perio but ns-perio: do corners later
      END DO
      ! north fold treatment
      IF( l_IdoNFold ) THEN
         ifill_nfd = ifill(jpno)             ! if we are here, this means llrecv(jpno) = .false. and l_SelfPerio(jpno) = .false.
         ifill( (/jpno/) ) = jpfillnothing   ! we do north fold -> do nothing for northern halo
      ENDIF
      
      ! We first define the localization and size of the parts of the array that will be sent (s), received (r)
      ! or used for periodocity (p). The localization is defined as "the bottom left corner - 1" in i and j directions.
      ! This is a shift that will be applied later in the do loops to pick-up the appropriate part of the array
      !
      ! all definitions bellow do not refer to N[ij][se]0 so we can use it with any local value of ihls
      !                   !                       ________________________
      ip0i =          0   !          im0j = inner |__|________________|__|
      ip1i =       ihls   !   im1j = inner - halo |  |__|__________|__|  |
      im1i = ipi-2*ihls   !                       |  |  |          |  |  |
      im0i = ipi - ihls   !                       |  |  |          |  |  |
      ip0j =          0   !                       |  |  |          |  |  |
      ip1j =       ihls   !                       |  |__|__________|__|  |
      im1j = ipj-2*ihls   !           ip1j = halo |__|__|__________|__|__|
      im0j = ipj - ihls   !              ip0j = 0 |__|________________|__|
      !                   !                    ip0i ip1i        im1i im0i
      !
      iwewe(:) = (/ jpwe,jpea,jpwe,jpea /)   ;   issnn(:) = (/ jpso,jpso,jpno,jpno /)
      !     sides:     west  east south north      ;   corners: so-we, so-ea, no-we, no-ea
      isizei(1:4) = (/ ihls, ihls, Ni_0, Ni_0 /)   ;   isizei(5:8) = ihls              ! i- count
      isizej(1:4) = (/ Nj_0, Nj_0, ihls, ihls /)   ;   isizej(5:8) = ihls              ! j- count
      ishtSi(1:4) = (/ ip1i, im1i, ip1i, ip1i /)   ;   ishtSi(5:8) = ishtSi( iwewe )   ! i- shift send data
      ishtSj(1:4) = (/ ip1j, ip1j, ip1j, im1j /)   ;   ishtSj(5:8) = ishtSj( issnn )   ! j- shift send data
      ishtRi(1:4) = (/ ip0i, im0i, ip1i, ip1i /)   ;   ishtRi(5:8) = ishtRi( iwewe )   ! i- shift received data location
      ishtRj(1:4) = (/ ip1j, ip1j, ip0j, im0j /)   ;   ishtRj(5:8) = ishtRj( issnn )   ! j- shift received data location
      ishtPi(1:4) = (/ im1i, ip1i, ip1i, ip1i /)   ;   ishtPi(5:8) = ishtPi( iwewe )   ! i- shift data used for periodicity
      ishtPj(1:4) = (/ ip1j, ip1j, im1j, ip1j /)   ;   ishtPj(5:8) = ishtPj( issnn )   ! j- shift data used for periodicity
      !
      ! -------------------------------- !
      !     2. Prepare MPI exchanges     !
      ! -------------------------------- !
      !
      ! Allocate local temporary arrays to be sent/received.
      iszS = COUNT( llsend )
      iszR = COUNT( llrecv )
      ALLOCATE( iScnt(iszS), iRcnt(iszR), iSdpl(iszS), iRdpl(iszR) )   ! ok if iszS = 0 or iszR = 0
      iszall(:) = isizei(:) * isizej(:) * ipk * ipl * ipf
      iScnt(:) = PACK( iszall, mask = llsend )                                       ! ok if mask = .false.
      iRcnt(:) = PACK( iszall, mask = llrecv )
      IF( iszS > 0 )   iSdpl(1) = 0
      DO jn = 2,iszS
         iSdpl(jn) = iSdpl(jn-1) + iScnt(jn-1)   ! with _alltoallv: in units of sendtype
      END DO
      IF( iszR > 0 )   iRdpl(1) = 0
      DO jn = 2,iszR
         iRdpl(jn) = iRdpl(jn-1) + iRcnt(jn-1)   ! with _alltoallv: in units of sendtype
      END DO
      
      ! Allocate buffer arrays to be sent/received if needed
      iszS = SUM(iszall, mask = llsend)                             ! send buffer size
      IF( ALLOCATED(buffsnd_sp) ) THEN
         IF( SIZE(buffsnd_sp) < iszS )    DEALLOCATE(buffsnd_sp)          ! send buffer is too small
      ENDIF
      IF( .NOT. ALLOCATED(buffsnd_sp) )   ALLOCATE( buffsnd_sp(iszS) )
      iszR = SUM(iszall, mask = llrecv)                             ! recv buffer size
      IF( ALLOCATED(buffrcv_sp) ) THEN
         IF( SIZE(buffrcv_sp) < iszR )    DEALLOCATE(buffrcv_sp)          ! recv buffer is too small
      ENDIF
      IF( .NOT. ALLOCATED(buffrcv_sp) )   ALLOCATE( buffrcv_sp(iszR) )

      ! fill sending buffer with ptab(jf)%pt4d
      idx = 1
      DO jn = 1, 8
         IF( llsend(jn) ) THEN
            ishti = ishtSi(jn)
            ishtj = ishtSj(jn)
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
               buffsnd_sp(idx) = ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl)
               idx = idx + 1
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         ENDIF
      END DO
      !
      ! ------------------------------------------------ !
      !     3. Do all MPI exchanges in 1 unique call     !
      ! ------------------------------------------------ !
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      CALL mpi_neighbor_alltoallv (buffsnd_sp, iScnt, iSdpl, MPI_REAL, buffrcv_sp, iRcnt, iRdpl, MPI_REAL, impi_nc, ierr)
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      ! ------------------------- !
      !     4. Fill all halos     !
      ! ------------------------- !
      !
      idx = 1
      ! MPI3 bug fix when domain decomposition has 2 columns/rows
      IF (jpni .eq. 2) THEN
         IF (jpnj .eq. 2) THEN
            jnf(1:8) = (/ 2, 1, 4, 3, 8, 7, 6, 5 /)
         ELSE
            jnf(1:8) = (/ 2, 1, 3, 4, 6, 5, 8, 7 /)
         ENDIF
      ELSE
         IF (jpnj .eq. 2) THEN
            jnf(1:8) = (/ 1, 2, 4, 3, 7, 8, 5, 6 /)
         ELSE
            jnf(1:8) = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
         ENDIF
      ENDIF

      DO jn = 1, 8
         ishti = ishtRi(jnf(jn))
         ishtj = ishtRj(jnf(jn))
         SELECT CASE ( ifill(jnf(jn)) )
         CASE ( jpfillnothing )               ! no filling 
         CASE ( jpfillmpi   )                 ! fill with data received by MPI
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jnf(jn))  ;  DO ji = 1,isizei(jnf(jn))
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = buffrcv_sp(idx)
               idx = idx + 1
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         CASE ( jpfillperio )                 ! use periodicity
            ishti2 = ishtPi(jnf(jn))
            ishtj2 = ishtPj(jnf(jn))
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jnf(jn))  ;  DO ji = 1,isizei(jnf(jn))
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         CASE ( jpfillcopy  )                 ! filling with inner domain values
            ishti2 = ishtSi(jnf(jn))
            ishtj2 = ishtSj(jnf(jn))
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jnf(jn))  ;  DO ji = 1,isizei(jnf(jn))
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         CASE ( jpfillcst   )                 ! filling with constant value
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jnf(jn))  ;  DO ji = 1,isizei(jnf(jn))
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = zland
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         END SELECT
      END DO

      DEALLOCATE( iScnt, iRcnt, iSdpl, iRdpl )
      IF( iszS > jpi*jpj )   DEALLOCATE(buffsnd_sp)                    ! blocking Send -> can directly deallocate
      IF( iszR > jpi*jpj )   DEALLOCATE(buffrcv_sp)                    ! blocking Recv -> can directly deallocate

      ! potential "indirect self-periodicity" for the corners
      DO jn = 5, 8
         IF( .NOT. l_SelfPerio(jn) .AND. l_SelfPerio(jpwe)  ) THEN   ! no bi-perio but ew-perio: corners indirect definition
            ishti  = ishtRi(jn)
            ishtj  = ishtRj(jn)
            ishti2 = ishtPi(jn)   ! use i- shift periodicity
            ishtj2 = ishtRj(jn)   ! use j- shift recv location: use ew-perio -> ok as filling of the south and north halos now done
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         ENDIF
         IF( .NOT. l_SelfPerio(jn) .AND. l_SelfPerio(jpso)  ) THEN   ! no bi-perio but ns-perio: corners indirect definition
            ishti  = ishtRi(jn)
            ishtj  = ishtRj(jn)
            ishti2 = ishtRi(jn)   ! use i- shift recv location: use ns-perio -> ok as filling of the west and east halos now done
            ishtj2 = ishtPj(jn)   ! use j- shift periodicity
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         ENDIF
      END DO
      !
      ! ------------------------------- !
      !     5. north fold treatment     !
      ! ------------------------------- !
      !
      IF( l_IdoNFold ) THEN
         IF( jpni == 1 )  THEN   ;   CALL lbc_nfd( ptab, cd_nat, psgn                  , ihls, ipf )   ! self NFold
         ELSE                    ;   CALL mpp_nfd( ptab, cd_nat, psgn, ifill_nfd, zland, ihls, ipf )   ! mpi  NFold
         ENDIF
      ENDIF
      !
   END SUBROUTINE lbc_lnk_neicoll_sp
   !!
   !!   ----   DOUBLE PRECISION VERSIONS
   !!
   SUBROUTINE lbc_lnk_pt2pt_dp( cdname, ptab, cd_nat, psgn, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      CHARACTER(len=*)              , INTENT(in   ) ::   cdname      ! name of the calling subroutine
      TYPE(PTR_4d_dp),  DIMENSION(:), INTENT(inout) ::   ptab        ! pointer of arrays on which apply the b.c.
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(dp),  DIMENSION(:), INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER                       , INTENT(in   ) ::   kfld        ! number of pt3d arrays
      INTEGER ,             OPTIONAL, INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(dp),      OPTIONAL, INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER ,             OPTIONAL, INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8),OPTIONAL, INTENT(in   ) ::   lsend, lrecv  ! communication with other 4 proc
      LOGICAL,              OPTIONAL, INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !
      INTEGER  ::     ji,   jj,  jk,  jl,  jf, jn     ! dummy loop indices
      INTEGER  ::    ipi,  ipj, ipk, ipl, ipf         ! dimension of the input array
      INTEGER  ::   ip0i, ip1i, im0i, im1i
      INTEGER  ::   ip0j, ip1j, im0j, im1j
      INTEGER  ::   ishti, ishtj, ishti2, ishtj2
      INTEGER  ::   ifill_nfd, icomm, ierr
      INTEGER  ::   ihls, idxs, idxr, iszS, iszR
      INTEGER, DIMENSION(4)  ::   iwewe, issnn
      INTEGER, DIMENSION(8)  ::   isizei, ishtSi, ishtRi, ishtPi
      INTEGER, DIMENSION(8)  ::   isizej, ishtSj, ishtRj, ishtPj
      INTEGER, DIMENSION(8)  ::   ifill, iszall, ishtS, ishtR
      INTEGER, DIMENSION(8)  ::   ireq             ! mpi_request id
      INTEGER, DIMENSION(8)  ::   iStag, iRtag     ! Send and Recv mpi_tag id
      REAL(dp) ::   zland
      LOGICAL, DIMENSION(8)  ::   llsend, llrecv
      LOGICAL  ::   ll4only                                        ! default: 8 neighbourgs
      !!----------------------------------------------------------------------
      !
      ! ----------------------------------------- !
      !     1. local variables initialization     !
      ! ----------------------------------------- !
      !
      ipi = SIZE(ptab(1)%pt4d,1)
      ipj = SIZE(ptab(1)%pt4d,2)
      ipk = SIZE(ptab(1)%pt4d,3)
      ipl = SIZE(ptab(1)%pt4d,4)
      ipf = kfld
      !
      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, ipk, ipl, ipf, ld_lbc = .TRUE. )
      !
      idxs = 1   ! initalize index for send buffer
      idxr = 1   ! initalize index for recv buffer
      icomm = mpi_comm_oce        ! shorter name
      !
      ! take care of optional parameters
      !
      ihls = nn_hls   ! default definition
      IF( PRESENT( khls ) )   ihls = khls
      IF( ihls > n_hlsmax ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with khls > n_hlsmax : ', khls, '>', n_hlsmax
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      IF( ipi /= Ni_0+2*ihls ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with an input array which does not match ihls along i: ', ipi, ihls, Ni_0
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      IF( ipj /= Nj_0+2*ihls ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with an input array which does not match ihls along j:', ipj, ihls , Nj_0
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      !
      ll4only = .FALSE.    ! default definition
      IF( PRESENT(ld4only) )   ll4only = ld4only
      !
      zland = 0._wp                                     ! land filling value: zero by default
      IF( PRESENT( pfillval ) )   zland = pfillval      ! set land value
      !
      ! define llsend and llrecv: logicals which say if mpi-neibourgs for send or receive exist or not.
      IF     ( PRESENT(lsend) .AND. PRESENT(lrecv) ) THEN   ! localy defined neighbourgs 
         llsend(:) = lsend(:)   ;   llrecv(:) = lrecv(:)
      ELSE IF( PRESENT(lsend) .OR.  PRESENT(lrecv) ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with only one of the two arguments lsend or lrecv'
         CALL ctl_stop( 'STOP', ctmp1 )
      ELSE                                              ! default neighbours
         llsend(:) = mpiSnei(ihls,:) >= 0
         IF( ll4only )   llsend(5:8) = .FALSE.          ! exclude corners
         llrecv(:) = mpiRnei(ihls,:) >= 0
         IF( ll4only )   llrecv(5:8) = .FALSE.          ! exclude corners
      ENDIF
      !
      ! define ifill: which method should be used to fill each parts (sides+corners) of the halos
      ! default definition
      DO jn = 1, 4
         IF(             llrecv(jn) ) THEN   ;   ifill(jn) = jpfillmpi    ! with an mpi communication
         ELSEIF(    l_SelfPerio(jn) ) THEN   ;   ifill(jn) = jpfillperio  ! with self-periodicity
         ELSEIF( PRESENT(kfillmode) ) THEN   ;   ifill(jn) = kfillmode    ! localy defined
         ELSE                                ;   ifill(jn) = jpfillcst    ! constant value (zland)
         ENDIF
      END DO
      DO jn = 5, 8
         IF(             llrecv(jn) ) THEN   ;   ifill(jn) = jpfillmpi    ! with an mpi communication
         ELSE                                ;   ifill(jn) = jpfillnothing! do nothing
         ENDIF
      END DO
         !
      ! north fold treatment
      IF( l_IdoNFold ) THEN
         ifill_nfd = ifill(jpno)             ! if we are here, this means llrecv(jpno) = .false. and l_SelfPerio(jpno) = .false.
         ifill( (/jpno/) ) = jpfillnothing   ! we do north fold -> do nothing for northern halo
      ENDIF
      
      ! We first define the localization and size of the parts of the array that will be sent (s), received (r)
      ! or used for periodocity (p). The localization is defined as "the bottom left corner - 1" in i and j directions.
      ! This is a shift that will be applied later in the do loops to pick-up the appropriate part of the array
      !
      ! all definitions bellow do not refer to N[ij][se]0 so we can use it with any local value of ihls
      !                   !                       ________________________
      ip0i =          0   !          im0j = inner |__|__|__________|__|__|
      ip1i =       ihls   !   im1j = inner - halo |__|__|__________|__|__|
      im1i = ipi-2*ihls   !                       |  |  |          |  |  |
      im0i = ipi - ihls   !                       |  |  |          |  |  |
      ip0j =          0   !                       |  |  |          |  |  |
      ip1j =       ihls   !                       |__|__|__________|__|__|
      im1j = ipj-2*ihls   !           ip1j = halo |__|__|__________|__|__|
      im0j = ipj - ihls   !              ip0j = 0 |__|__|__________|__|__|
      !                   !                    ip0i ip1i        im1i im0i
      !
      iwewe(:) = (/ jpwe,jpea,jpwe,jpea /)   ;   issnn(:) = (/ jpso,jpso,jpno,jpno /)
      !     sides:     west  east south north      ;   corners: so-we, so-ea, no-we, no-ea
      isizei(1:4) = (/ ihls, ihls,  ipi,  ipi /)   ;   isizei(5:8) = ihls              ! i- count
      isizej(1:4) = (/ Nj_0, Nj_0, ihls, ihls /)   ;   isizej(5:8) = ihls              ! j- count
      ishtSi(1:4) = (/ ip1i, im1i, ip0i, ip0i /)   ;   ishtSi(5:8) = ishtSi( iwewe )   ! i- shift send data
      ishtSj(1:4) = (/ ip1j, ip1j, ip1j, im1j /)   ;   ishtSj(5:8) = ishtSj( issnn )   ! j- shift send data
      ishtRi(1:4) = (/ ip0i, im0i, ip0i, ip0i /)   ;   ishtRi(5:8) = ishtRi( iwewe )   ! i- shift received data location
      ishtRj(1:4) = (/ ip1j, ip1j, ip0j, im0j /)   ;   ishtRj(5:8) = ishtRj( issnn )   ! j- shift received data location
      ishtPi(1:4) = (/ im1i, ip1i, ip0i, ip0i /)   ;   ishtPi(5:8) = ishtPi( iwewe )   ! i- shift data used for periodicity
      ishtPj(1:4) = (/ ip1j, ip1j, im1j, ip1j /)   ;   ishtPj(5:8) = ishtPj( issnn )   ! j- shift data used for periodicity
      !
      ! -------------------------------- !
      !     2. Prepare MPI exchanges     !
      ! -------------------------------- !
      !
      iStag = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)   ! any value but each one must be different
      ! define iRtag with the corresponding iStag, e.g. data received at west where sent at east.
      iRtag(jpwe) = iStag(jpea)   ;   iRtag(jpea) = iStag(jpwe)   ;   iRtag(jpso) = iStag(jpno)   ;   iRtag(jpno) = iStag(jpso)
      iRtag(jpsw) = iStag(jpne)   ;   iRtag(jpse) = iStag(jpnw)   ;   iRtag(jpnw) = iStag(jpse)   ;   iRtag(jpne) = iStag(jpsw)
      !
      iszall(:) = isizei(:) * isizej(:) * ipk * ipl * ipf
      ishtS(1) = 0
      DO jn = 2, 8
         ishtS(jn) = ishtS(jn-1) + iszall(jn-1) * COUNT( (/llsend(jn-1)/) )
      END DO
      ishtR(1) = 0
      DO jn = 2, 8
         ishtR(jn) = ishtR(jn-1) + iszall(jn-1) * COUNT( (/llrecv(jn-1)/) )
      END DO

      ! Allocate buffer arrays to be sent/received if needed
      iszS = SUM(iszall, mask = llsend)                             ! send buffer size
      IF( ALLOCATED(buffsnd_dp) ) THEN
         CALL mpi_waitall(8, nreq_p2p, MPI_STATUSES_IGNORE, ierr)   ! wait for Isend from the PREVIOUS call
         IF( SIZE(buffsnd_dp) < iszS )    DEALLOCATE(buffsnd_dp)          ! send buffer is too small
      ENDIF
      IF( .NOT. ALLOCATED(buffsnd_dp) )   ALLOCATE( buffsnd_dp(iszS) )
      iszR = SUM(iszall, mask = llrecv)                             ! recv buffer size
      IF( ALLOCATED(buffrcv_dp) ) THEN
         IF( SIZE(buffrcv_dp) < iszR )    DEALLOCATE(buffrcv_dp)          ! recv buffer is too small
      ENDIF
      IF( .NOT. ALLOCATED(buffrcv_dp) )   ALLOCATE( buffrcv_dp(iszR) )
      !
      ! default definition when no communication is done. understood by mpi_waitall
      nreq_p2p(:) = MPI_REQUEST_NULL   ! WARNING: Must be done after the call to mpi_waitall just above
      !
      ! ----------------------------------------------- !
      !     3. Do east and west MPI_Isend if needed     !
      ! ----------------------------------------------- !
      !
      DO jn = 1, 2

   IF( llsend(jn) ) THEN
      ishti = ishtSi(jn)
      ishtj = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         buffsnd_dp(idxs) = ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl)
         idxs = idxs + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      ! non-blocking send of the west/east side using local buffer
      CALL MPI_ISEND( buffsnd_dp(ishtS(jn)+1), iszall(jn), MPI_DOUBLE_PRECISION, mpiSnei(ihls,jn), iStag(jn), icomm, nreq_p2p(jn), ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
   ENDIF

      END DO
      !
      ! ----------------------------------- !
      !     4. Fill east and west halos     !
      ! ----------------------------------- !
      !
      DO jn = 1, 2


   ishti = ishtRi(jn)
   ishtj = ishtRj(jn)
   SELECT CASE ( ifill(jn) )
   CASE ( jpfillnothing )               ! no filling 
   CASE ( jpfillmpi   )                 ! fill with data received by MPI
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !                                 ! blocking receive of the west/east halo in local temporary arrays
      CALL MPI_RECV( buffrcv_dp(ishtR(jn)+1), iszall(jn), MPI_DOUBLE_PRECISION, mpiRnei(ihls,jn), iRtag(jn), icomm, MPI_STATUS_IGNORE, ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = buffrcv_dp(idxr)
         idxr = idxr + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillperio )                 ! use periodicity
      ishti2 = ishtPi(jn)
      ishtj2 = ishtPj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcopy  )                 ! filling with inner domain values
      ishti2 = ishtSi(jn)
      ishtj2 = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcst   )                 ! filling with constant value
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = zland
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   END SELECT
      END DO
      !
      ! ------------------------------------------------- !
      !     5. Do north and south MPI_Isend if needed     !
      ! ------------------------------------------------- !
      !
      DO jn = 3, 4

   IF( llsend(jn) ) THEN
      ishti = ishtSi(jn)
      ishtj = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         buffsnd_dp(idxs) = ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl)
         idxs = idxs + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      ! non-blocking send of the west/east side using local buffer
      CALL MPI_ISEND( buffsnd_dp(ishtS(jn)+1), iszall(jn), MPI_DOUBLE_PRECISION, mpiSnei(ihls,jn), iStag(jn), icomm, nreq_p2p(jn), ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
   ENDIF

      END DO
      !
      ! ------------------------------- !
      !     6. north fold treatment     !
      ! ------------------------------- !
      !
      ! Must be done after receiving data from East/West neighbourgs (as it is coded in mpp_nfd, to be changed one day...)
      ! Do it after MPI_iSend to south/north neighbourgs so they won't wait (too much) to receive their data
      ! Do if before MPI_Recv from south/north neighbourgs so we have more time to receive data
      !
      IF( l_IdoNFold ) THEN
         IF( jpni == 1 )  THEN   ;   CALL lbc_nfd( ptab, cd_nat, psgn                  , ihls, ipf )   ! self NFold
         ELSE                    ;   CALL mpp_nfd( ptab, cd_nat, psgn, ifill_nfd, zland, ihls, ipf )   ! mpi  NFold
         ENDIF
      ENDIF
      !
      ! ------------------------------------- !
      !     7. Fill south and north halos     !
      ! ------------------------------------- !
      !
      DO jn = 3, 4


   ishti = ishtRi(jn)
   ishtj = ishtRj(jn)
   SELECT CASE ( ifill(jn) )
   CASE ( jpfillnothing )               ! no filling 
   CASE ( jpfillmpi   )                 ! fill with data received by MPI
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !                                 ! blocking receive of the west/east halo in local temporary arrays
      CALL MPI_RECV( buffrcv_dp(ishtR(jn)+1), iszall(jn), MPI_DOUBLE_PRECISION, mpiRnei(ihls,jn), iRtag(jn), icomm, MPI_STATUS_IGNORE, ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = buffrcv_dp(idxr)
         idxr = idxr + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillperio )                 ! use periodicity
      ishti2 = ishtPi(jn)
      ishtj2 = ishtPj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcopy  )                 ! filling with inner domain values
      ishti2 = ishtSi(jn)
      ishtj2 = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcst   )                 ! filling with constant value
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = zland
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   END SELECT
      END DO
      !
      ! ----------------------------------------------- !
      !     8. Specific problem in corner treatment     !
      !              ( very rate case... )              !
      ! ----------------------------------------------- !
      !
      DO jn = 5, 8

   IF( llsend(jn) ) THEN
      ishti = ishtSi(jn)
      ishtj = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         buffsnd_dp(idxs) = ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl)
         idxs = idxs + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      ! non-blocking send of the west/east side using local buffer
      CALL MPI_ISEND( buffsnd_dp(ishtS(jn)+1), iszall(jn), MPI_DOUBLE_PRECISION, mpiSnei(ihls,jn), iStag(jn), icomm, nreq_p2p(jn), ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
   ENDIF

      END DO
      DO jn = 5, 8


   ishti = ishtRi(jn)
   ishtj = ishtRj(jn)
   SELECT CASE ( ifill(jn) )
   CASE ( jpfillnothing )               ! no filling 
   CASE ( jpfillmpi   )                 ! fill with data received by MPI
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !                                 ! blocking receive of the west/east halo in local temporary arrays
      CALL MPI_RECV( buffrcv_dp(ishtR(jn)+1), iszall(jn), MPI_DOUBLE_PRECISION, mpiRnei(ihls,jn), iRtag(jn), icomm, MPI_STATUS_IGNORE, ierr )
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = buffrcv_dp(idxr)
         idxr = idxr + 1
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillperio )                 ! use periodicity
      ishti2 = ishtPi(jn)
      ishtj2 = ishtPj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcopy  )                 ! filling with inner domain values
      ishti2 = ishtSi(jn)
      ishtj2 = ishtSj(jn)
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   CASE ( jpfillcst   )                 ! filling with constant value
      DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
         ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = zland
      END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
   END SELECT
      END DO
      !
      ! -------------------------------------------- !
      !     9. deallocate local temporary arrays     !
      !        if they areg larger than jpi*jpj      !  <- arbitrary max size...
      ! -------------------------------------------- !
      !
      IF( iszR > jpi*jpj )   DEALLOCATE(buffrcv_dp)                    ! blocking receive -> can directly deallocate
      IF( iszS > jpi*jpj ) THEN
         CALL mpi_waitall(8, nreq_p2p, MPI_STATUSES_IGNORE, ierr)   ! must wait before deallocate send buffer
         DEALLOCATE(buffsnd_dp)
      ENDIF
      !
   END SUBROUTINE lbc_lnk_pt2pt_dp


SUBROUTINE lbc_lnk_neicoll_dp( cdname, ptab, cd_nat, psgn, kfld, kfillmode, pfillval, khls, lsend, lrecv, ld4only )
      CHARACTER(len=*)              , INTENT(in   ) ::   cdname      ! name of the calling subroutine
      TYPE(PTR_4d_dp),  DIMENSION(:), INTENT(inout) ::   ptab        ! pointer of arrays on which apply the b.c.
      CHARACTER(len=1), DIMENSION(:), INTENT(in   ) ::   cd_nat      ! nature of array grid-points
      REAL(dp),  DIMENSION(:), INTENT(in   ) ::   psgn        ! sign used across the north fold boundary
      INTEGER                       , INTENT(in   ) ::   kfld        ! number of pt3d arrays
      INTEGER ,             OPTIONAL, INTENT(in   ) ::   kfillmode   ! filling method for halo over land (default = constant)
      REAL(dp),      OPTIONAL, INTENT(in   ) ::   pfillval    ! background value (used at closed boundaries)
      INTEGER ,             OPTIONAL, INTENT(in   ) ::   khls        ! halo size, default = nn_hls
      LOGICAL, DIMENSION(8),OPTIONAL, INTENT(in   ) ::   lsend, lrecv  ! communication with other 4 proc
      LOGICAL,              OPTIONAL, INTENT(in   ) ::   ld4only     ! if .T., do only 4-neighbour comm (ignore corners)
      !
      INTEGER  ::    ji,  jj,  jk , jl,  jf, jn      ! dummy loop indices
      INTEGER  ::   ipi, ipj, ipk, ipl, ipf          ! dimension of the input array
      INTEGER  ::   ip0i, ip1i, im0i, im1i
      INTEGER  ::   ip0j, ip1j, im0j, im1j
      INTEGER  ::   ishti, ishtj, ishti2, ishtj2
      INTEGER  ::   iszS, iszR
      INTEGER  ::   ierr
      INTEGER  ::   ihls, idx
      INTEGER  ::   impi_nc
      INTEGER  ::   ifill_nfd
      INTEGER, DIMENSION(4)  ::   iwewe, issnn
      INTEGER, DIMENSION(8)  ::   isizei, ishtSi, ishtRi, ishtPi
      INTEGER, DIMENSION(8)  ::   isizej, ishtSj, ishtRj, ishtPj
      INTEGER, DIMENSION(8)  ::   ifill, iszall
      INTEGER, DIMENSION(8)  ::   jnf
      INTEGER, DIMENSION(:), ALLOCATABLE  ::   iScnt, iRcnt    ! number of elements to be sent/received
      INTEGER, DIMENSION(:), ALLOCATABLE  ::   iSdpl, iRdpl    ! displacement in halos arrays
      LOGICAL, DIMENSION(8)  ::   llsend, llrecv
      REAL(dp) ::   zland
      LOGICAL  ::   ll4only                                    ! default: 8 neighbourgs
      !!----------------------------------------------------------------------
      !
      ! ----------------------------------------- !
      !     1. local variables initialization     !
      ! ----------------------------------------- !
      !
      ipi = SIZE(ptab(1)%pt4d,1)
      ipj = SIZE(ptab(1)%pt4d,2)
      ipk = SIZE(ptab(1)%pt4d,3)
      ipl = SIZE(ptab(1)%pt4d,4)
      ipf = kfld
      !
      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, ipk, ipl, ipf, ld_lbc = .TRUE. )
      !
      ! take care of optional parameters
      !
      ihls = nn_hls       ! default definition
      IF( PRESENT( khls ) )   ihls = khls
      IF( ihls > n_hlsmax ) THEN
         WRITE(ctmp1,*) TRIM(cdname), '  is calling lbc_lnk with khls > n_hlsmax : ', khls, '>', n_hlsmax
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      IF( ipi /= Ni_0+2*ihls ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with an input array which does not match ihls along i: ', ipi, ihls, Ni_0
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      IF( ipj /= Nj_0+2*ihls ) THEN
         WRITE(ctmp1,*) TRIM(cdname), ' is calling lbc_lnk with an input array which does not match ihls along j:', ipj, ihls , Nj_0
         CALL ctl_stop( 'STOP', ctmp1 )
      ENDIF
      !
      ll4only = .FALSE.    ! default definition
      IF( PRESENT(ld4only) )   ll4only = ld4only
      !
      impi_nc = mpi_nc_com8(ihls)   ! default
      IF( ll4only )   impi_nc = mpi_nc_com4(ihls)
      !
      zland = 0._wp                                     ! land filling value: zero by default
      IF( PRESENT( pfillval ) )   zland = pfillval      ! set land value
      !
      ! define llsend and llrecv: logicals which say if mpi-neibourgs for send or receive exist or not.
      IF     ( PRESENT(lsend) .AND. PRESENT(lrecv) ) THEN   ! localy defined neighbourgs 
         CALL ctl_stop( 'STOP', 'mpp_nc_generic+lsend and lrecv not yet implemented')
      ELSE IF( PRESENT(lsend) .OR.  PRESENT(lrecv) ) THEN
         WRITE(ctmp1,*) TRIM(cdname), '  is calling lbc_lnk with only one of the two arguments lsend or lrecv'
         CALL ctl_stop( 'STOP', ctmp1 )
      ELSE                                              ! default neighbours
         llsend(:) = mpiSnei(ihls,:) >= 0
         IF( ll4only )   llsend(5:8) = .FALSE.          ! exclude corners
         llrecv(:) = mpiRnei(ihls,:) >= 0
         IF( ll4only )   llrecv(5:8) = .FALSE.          ! exclude corners
      ENDIF
      !
      ! define ifill: which method should be used to fill each parts (sides+corners) of the halos
      ! default definition
      DO jn = 1, 8
         IF(             llrecv(jn) ) THEN   ;   ifill(jn) = jpfillmpi    ! with an mpi communication
         ELSEIF(    l_SelfPerio(jn) ) THEN   ;   ifill(jn) = jpfillperio  ! with self-periodicity
         ELSEIF( PRESENT(kfillmode) ) THEN   ;   ifill(jn) = kfillmode    ! localy defined
         ELSE                                ;   ifill(jn) = jpfillcst    ! constant value (zland)
         ENDIF
      END DO
      ! take care of "indirect self-periodicity" for the corners
      DO jn = 5, 8
         IF(.NOT.l_SelfPerio(jn) .AND. l_SelfPerio(jpwe))   ifill(jn) = jpfillnothing   ! no bi-perio but ew-perio: do corners later
         IF(.NOT.l_SelfPerio(jn) .AND. l_SelfPerio(jpso))   ifill(jn) = jpfillnothing   ! no bi-perio but ns-perio: do corners later
      END DO
      ! north fold treatment
      IF( l_IdoNFold ) THEN
         ifill_nfd = ifill(jpno)             ! if we are here, this means llrecv(jpno) = .false. and l_SelfPerio(jpno) = .false.
         ifill( (/jpno/) ) = jpfillnothing   ! we do north fold -> do nothing for northern halo
      ENDIF
      
      ! We first define the localization and size of the parts of the array that will be sent (s), received (r)
      ! or used for periodocity (p). The localization is defined as "the bottom left corner - 1" in i and j directions.
      ! This is a shift that will be applied later in the do loops to pick-up the appropriate part of the array
      !
      ! all definitions bellow do not refer to N[ij][se]0 so we can use it with any local value of ihls
      !                   !                       ________________________
      ip0i =          0   !          im0j = inner |__|________________|__|
      ip1i =       ihls   !   im1j = inner - halo |  |__|__________|__|  |
      im1i = ipi-2*ihls   !                       |  |  |          |  |  |
      im0i = ipi - ihls   !                       |  |  |          |  |  |
      ip0j =          0   !                       |  |  |          |  |  |
      ip1j =       ihls   !                       |  |__|__________|__|  |
      im1j = ipj-2*ihls   !           ip1j = halo |__|__|__________|__|__|
      im0j = ipj - ihls   !              ip0j = 0 |__|________________|__|
      !                   !                    ip0i ip1i        im1i im0i
      !
      iwewe(:) = (/ jpwe,jpea,jpwe,jpea /)   ;   issnn(:) = (/ jpso,jpso,jpno,jpno /)
      !     sides:     west  east south north      ;   corners: so-we, so-ea, no-we, no-ea
      isizei(1:4) = (/ ihls, ihls, Ni_0, Ni_0 /)   ;   isizei(5:8) = ihls              ! i- count
      isizej(1:4) = (/ Nj_0, Nj_0, ihls, ihls /)   ;   isizej(5:8) = ihls              ! j- count
      ishtSi(1:4) = (/ ip1i, im1i, ip1i, ip1i /)   ;   ishtSi(5:8) = ishtSi( iwewe )   ! i- shift send data
      ishtSj(1:4) = (/ ip1j, ip1j, ip1j, im1j /)   ;   ishtSj(5:8) = ishtSj( issnn )   ! j- shift send data
      ishtRi(1:4) = (/ ip0i, im0i, ip1i, ip1i /)   ;   ishtRi(5:8) = ishtRi( iwewe )   ! i- shift received data location
      ishtRj(1:4) = (/ ip1j, ip1j, ip0j, im0j /)   ;   ishtRj(5:8) = ishtRj( issnn )   ! j- shift received data location
      ishtPi(1:4) = (/ im1i, ip1i, ip1i, ip1i /)   ;   ishtPi(5:8) = ishtPi( iwewe )   ! i- shift data used for periodicity
      ishtPj(1:4) = (/ ip1j, ip1j, im1j, ip1j /)   ;   ishtPj(5:8) = ishtPj( issnn )   ! j- shift data used for periodicity
      !
      ! -------------------------------- !
      !     2. Prepare MPI exchanges     !
      ! -------------------------------- !
      !
      ! Allocate local temporary arrays to be sent/received.
      iszS = COUNT( llsend )
      iszR = COUNT( llrecv )
      ALLOCATE( iScnt(iszS), iRcnt(iszR), iSdpl(iszS), iRdpl(iszR) )   ! ok if iszS = 0 or iszR = 0
      iszall(:) = isizei(:) * isizej(:) * ipk * ipl * ipf
      iScnt(:) = PACK( iszall, mask = llsend )                                       ! ok if mask = .false.
      iRcnt(:) = PACK( iszall, mask = llrecv )
      IF( iszS > 0 )   iSdpl(1) = 0
      DO jn = 2,iszS
         iSdpl(jn) = iSdpl(jn-1) + iScnt(jn-1)   ! with _alltoallv: in units of sendtype
      END DO
      IF( iszR > 0 )   iRdpl(1) = 0
      DO jn = 2,iszR
         iRdpl(jn) = iRdpl(jn-1) + iRcnt(jn-1)   ! with _alltoallv: in units of sendtype
      END DO
      
      ! Allocate buffer arrays to be sent/received if needed
      iszS = SUM(iszall, mask = llsend)                             ! send buffer size
      IF( ALLOCATED(buffsnd_dp) ) THEN
         IF( SIZE(buffsnd_dp) < iszS )    DEALLOCATE(buffsnd_dp)          ! send buffer is too small
      ENDIF
      IF( .NOT. ALLOCATED(buffsnd_dp) )   ALLOCATE( buffsnd_dp(iszS) )
      iszR = SUM(iszall, mask = llrecv)                             ! recv buffer size
      IF( ALLOCATED(buffrcv_dp) ) THEN
         IF( SIZE(buffrcv_dp) < iszR )    DEALLOCATE(buffrcv_dp)          ! recv buffer is too small
      ENDIF
      IF( .NOT. ALLOCATED(buffrcv_dp) )   ALLOCATE( buffrcv_dp(iszR) )

      ! fill sending buffer with ptab(jf)%pt4d
      idx = 1
      DO jn = 1, 8
         IF( llsend(jn) ) THEN
            ishti = ishtSi(jn)
            ishtj = ishtSj(jn)
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
               buffsnd_dp(idx) = ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl)
               idx = idx + 1
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         ENDIF
      END DO
      !
      ! ------------------------------------------------ !
      !     3. Do all MPI exchanges in 1 unique call     !
      ! ------------------------------------------------ !
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      CALL mpi_neighbor_alltoallv (buffsnd_dp, iScnt, iSdpl, MPI_DOUBLE_PRECISION, buffrcv_dp, iRcnt, iRdpl, MPI_DOUBLE_PRECISION, impi_nc, ierr)
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      ! ------------------------- !
      !     4. Fill all halos     !
      ! ------------------------- !
      !
      idx = 1
      ! MPI3 bug fix when domain decomposition has 2 columns/rows
      IF (jpni .eq. 2) THEN
         IF (jpnj .eq. 2) THEN
            jnf(1:8) = (/ 2, 1, 4, 3, 8, 7, 6, 5 /)
         ELSE
            jnf(1:8) = (/ 2, 1, 3, 4, 6, 5, 8, 7 /)
         ENDIF
      ELSE
         IF (jpnj .eq. 2) THEN
            jnf(1:8) = (/ 1, 2, 4, 3, 7, 8, 5, 6 /)
         ELSE
            jnf(1:8) = (/ 1, 2, 3, 4, 5, 6, 7, 8 /)
         ENDIF
      ENDIF

      DO jn = 1, 8
         ishti = ishtRi(jnf(jn))
         ishtj = ishtRj(jnf(jn))
         SELECT CASE ( ifill(jnf(jn)) )
         CASE ( jpfillnothing )               ! no filling 
         CASE ( jpfillmpi   )                 ! fill with data received by MPI
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jnf(jn))  ;  DO ji = 1,isizei(jnf(jn))
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = buffrcv_dp(idx)
               idx = idx + 1
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         CASE ( jpfillperio )                 ! use periodicity
            ishti2 = ishtPi(jnf(jn))
            ishtj2 = ishtPj(jnf(jn))
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jnf(jn))  ;  DO ji = 1,isizei(jnf(jn))
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         CASE ( jpfillcopy  )                 ! filling with inner domain values
            ishti2 = ishtSi(jnf(jn))
            ishtj2 = ishtSj(jnf(jn))
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jnf(jn))  ;  DO ji = 1,isizei(jnf(jn))
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         CASE ( jpfillcst   )                 ! filling with constant value
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jnf(jn))  ;  DO ji = 1,isizei(jnf(jn))
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = zland
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         END SELECT
      END DO

      DEALLOCATE( iScnt, iRcnt, iSdpl, iRdpl )
      IF( iszS > jpi*jpj )   DEALLOCATE(buffsnd_dp)                    ! blocking Send -> can directly deallocate
      IF( iszR > jpi*jpj )   DEALLOCATE(buffrcv_dp)                    ! blocking Recv -> can directly deallocate

      ! potential "indirect self-periodicity" for the corners
      DO jn = 5, 8
         IF( .NOT. l_SelfPerio(jn) .AND. l_SelfPerio(jpwe)  ) THEN   ! no bi-perio but ew-perio: corners indirect definition
            ishti  = ishtRi(jn)
            ishtj  = ishtRj(jn)
            ishti2 = ishtPi(jn)   ! use i- shift periodicity
            ishtj2 = ishtRj(jn)   ! use j- shift recv location: use ew-perio -> ok as filling of the south and north halos now done
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         ENDIF
         IF( .NOT. l_SelfPerio(jn) .AND. l_SelfPerio(jpso)  ) THEN   ! no bi-perio but ns-perio: corners indirect definition
            ishti  = ishtRi(jn)
            ishtj  = ishtRj(jn)
            ishti2 = ishtRi(jn)   ! use i- shift recv location: use ns-perio -> ok as filling of the west and east halos now done
            ishtj2 = ishtPj(jn)   ! use j- shift periodicity
            DO jf = 1, ipf  ;  DO jl = 1, ipl  ;  DO jk = 1, ipk  ;  DO jj = 1,isizej(jn)  ;  DO ji = 1,isizei(jn)
               ptab(jf)%pt4d(ishti+ji,ishtj+jj,jk,jl) = ptab(jf)%pt4d(ishti2+ji,ishtj2+jj,jk,jl)
            END DO   ;   END DO   ;   END DO   ;   END DO   ;   END DO
         ENDIF
      END DO
      !
      ! ------------------------------- !
      !     5. north fold treatment     !
      ! ------------------------------- !
      !
      IF( l_IdoNFold ) THEN
         IF( jpni == 1 )  THEN   ;   CALL lbc_nfd( ptab, cd_nat, psgn                  , ihls, ipf )   ! self NFold
         ELSE                    ;   CALL mpp_nfd( ptab, cd_nat, psgn, ifill_nfd, zland, ihls, ipf )   ! mpi  NFold
         ENDIF
      ENDIF
      !
   END SUBROUTINE lbc_lnk_neicoll_dp

   !!======================================================================
     !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_icb  ***
      !!
      !! ** Purpose :   Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1 and for 2d
      !!              array with outer extra halo
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4+kextj northern lines of the global domain on 1
      !!              processor and apply lbc north-fold on this sub array.
      !!              Then we scatter the north fold array back to the processors.
      !!              This routine accounts for an extra halo with icebergs
      !!              and assumes ghost rows and columns have been suppressed.
      !!
      !!----------------------------------------------------------------------

   SUBROUTINE mpp_lbc_north_icb_sp( pt2d, cd_type, psgn, kextj)
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_icb  ***
      !!
      !! ** Purpose :   Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1 and for 2d
      !!              array with outer extra halo
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4+kextj northern lines of the global domain on 1
      !!              processor and apply lbc north-fold on this sub array.
      !!              Then we scatter the north fold array back to the processors.
      !!              This routine accounts for an extra halo with icebergs
      !!              and assumes ghost rows and columns have been suppressed.
      !!
      !!----------------------------------------------------------------------
      REAL(sp), DIMENSION(:,:), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type  ! nature of pt3d grid-points
      !                                                     !   = T ,  U , V , F or W -points
      REAL(sp)         , INTENT(in   ) ::   psgn     ! = -1. the sign change across the
      !!                                                    ! north fold, =  1. otherwise
      INTEGER                 , INTENT(in   ) ::   kextj    ! Extra halo width at north fold
      !
      INTEGER ::   ji, jj, jr
      INTEGER ::   ierr, itaille
      INTEGER ::   ipj, ij, iproc, ijnr, ii1, ipi, impp
      !
      REAL(sp), DIMENSION(:,:)  , ALLOCATABLE  ::  ztab_e, znorthloc_e
      REAL(sp), DIMENSION(:,:,:), ALLOCATABLE  ::  znorthgloio_e
      !!----------------------------------------------------------------------
      !
      ipj=4
      ALLOCATE(        ztab_e(jpiglo, 1-kextj:ipj+kextj)       ,       &
     &            znorthloc_e(jpimax, 1-kextj:ipj+kextj)       ,       &
     &          znorthgloio_e(jpimax, 1-kextj:ipj+kextj,ndim_rank_north)    )
      !
      ztab_e(:,:)      = 0._sp
      znorthloc_e(:,:) = 0._sp
      !
      ij = 1 - kextj
      ! put the last ipj+2*kextj lines of pt2d into znorthloc_e 
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         znorthloc_e(1:jpi,ij)=pt2d(1:jpi,jj)
         ij = ij + 1
      END DO
      !
      itaille = jpimax * ( ipj + 2*kextj )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER( znorthloc_e(1,1-kextj)    , itaille, MPI_REAL,    &
         &                znorthgloio_e(1,1-kextj,1), itaille, MPI_REAL,    &
         &                ncomm_north, ierr )
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      ijnr = 0
      DO jr = 1, ndim_rank_north            ! recover the global north array
         iproc = nfproc(jr)
         IF( iproc /= -1 ) THEN
            impp = nfimpp(jr)
            ipi  = nfjpi(jr)
            ijnr = ijnr + 1
            DO jj = 1-kextj, ipj+kextj
               DO ji = 1, ipi
                  ii1 = impp + ji - 1       ! corresponds to mig(ji) but for subdomain iproc
                  ztab_e(ii1,jj) = znorthgloio_e(ji,jj,ijnr)
               END DO
            END DO
         ENDIF
      END DO

      ! 2. North-Fold boundary conditions
      ! ----------------------------------
      CALL lbc_nfd( ztab_e(:,1-kextj:ipj+kextj), cd_type, psgn, kextj )

      ij = 1 - kextj
      !! Scatter back to pt2d
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         DO ji= 1, jpi
            pt2d(ji,jj) = ztab_e(ji+nimpp-1,ij)
         END DO
         ij  = ij +1
      END DO
      !
      DEALLOCATE( ztab_e, znorthloc_e, znorthgloio_e )
      !
   END SUBROUTINE mpp_lbc_north_icb_sp 


   SUBROUTINE mpp_lbc_north_icb_dp( pt2d, cd_type, psgn, kextj)
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_icb  ***
      !!
      !! ** Purpose :   Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1 and for 2d
      !!              array with outer extra halo
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4+kextj northern lines of the global domain on 1
      !!              processor and apply lbc north-fold on this sub array.
      !!              Then we scatter the north fold array back to the processors.
      !!              This routine accounts for an extra halo with icebergs
      !!              and assumes ghost rows and columns have been suppressed.
      !!
      !!----------------------------------------------------------------------
      REAL(dp), DIMENSION(:,:), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type  ! nature of pt3d grid-points
      !                                                     !   = T ,  U , V , F or W -points
      REAL(dp)         , INTENT(in   ) ::   psgn     ! = -1. the sign change across the
      !!                                                    ! north fold, =  1. otherwise
      INTEGER                 , INTENT(in   ) ::   kextj    ! Extra halo width at north fold
      !
      INTEGER ::   ji, jj, jr
      INTEGER ::   ierr, itaille
      INTEGER ::   ipj, ij, iproc, ijnr, ii1, ipi, impp
      !
      REAL(dp), DIMENSION(:,:)  , ALLOCATABLE  ::  ztab_e, znorthloc_e
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE  ::  znorthgloio_e
      !!----------------------------------------------------------------------
      !
      ipj=4
      ALLOCATE(        ztab_e(jpiglo, 1-kextj:ipj+kextj)       ,       &
     &            znorthloc_e(jpimax, 1-kextj:ipj+kextj)       ,       &
     &          znorthgloio_e(jpimax, 1-kextj:ipj+kextj,ndim_rank_north)    )
      !
      ztab_e(:,:)      = 0._dp
      znorthloc_e(:,:) = 0._dp
      !
      ij = 1 - kextj
      ! put the last ipj+2*kextj lines of pt2d into znorthloc_e 
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         znorthloc_e(1:jpi,ij)=pt2d(1:jpi,jj)
         ij = ij + 1
      END DO
      !
      itaille = jpimax * ( ipj + 2*kextj )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      CALL MPI_ALLGATHER( znorthloc_e(1,1-kextj)    , itaille, MPI_DOUBLE_PRECISION,    &
         &                znorthgloio_e(1,1-kextj,1), itaille, MPI_DOUBLE_PRECISION,    &
         &                ncomm_north, ierr )
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      ijnr = 0
      DO jr = 1, ndim_rank_north            ! recover the global north array
         iproc = nfproc(jr)
         IF( iproc /= -1 ) THEN
            impp = nfimpp(jr)
            ipi  = nfjpi(jr)
            ijnr = ijnr + 1
            DO jj = 1-kextj, ipj+kextj
               DO ji = 1, ipi
                  ii1 = impp + ji - 1       ! corresponds to mig(ji) but for subdomain iproc
                  ztab_e(ii1,jj) = znorthgloio_e(ji,jj,ijnr)
               END DO
            END DO
         ENDIF
      END DO

      ! 2. North-Fold boundary conditions
      ! ----------------------------------
      CALL lbc_nfd( ztab_e(:,1-kextj:ipj+kextj), cd_type, psgn, kextj )

      ij = 1 - kextj
      !! Scatter back to pt2d
      DO jj = jpj - ipj + 1 - kextj , jpj + kextj
         DO ji= 1, jpi
            pt2d(ji,jj) = ztab_e(ji+nimpp-1,ij)
         END DO
         ij  = ij +1
      END DO
      !
      DEALLOCATE( ztab_e, znorthloc_e, znorthgloio_e )
      !
   END SUBROUTINE mpp_lbc_north_icb_dp 



      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_icb  ***
      !!
      !! ** Purpose :   Message passing management for 2d array (with extra halo for icebergs)
      !!                This routine receives a (1-kexti:jpi+kexti,1-kexti:jpj+kextj)
      !!                array (usually (0:jpi+1, 0:jpj+1)) from lbc_lnk_icb calls.
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    jpi    : first dimension of the local subdomain
      !!                    jpj    : second dimension of the local subdomain
      !!                    mpinei : number of neighboring domains (starting at 0, -1 if no neighbourg)
      !!----------------------------------------------------------------------


   SUBROUTINE mpp_lnk_2d_icb_sp( cdname, pt2d, cd_type, psgn, kexti, kextj )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_icb  ***
      !!
      !! ** Purpose :   Message passing management for 2d array (with extra halo for icebergs)
      !!                This routine receives a (1-kexti:jpi+kexti,1-kexti:jpj+kextj)
      !!                array (usually (0:jpi+1, 0:jpj+1)) from lbc_lnk_icb calls.
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    jpi    : first dimension of the local subdomain
      !!                    jpj    : second dimension of the local subdomain
      !!                    mpinei : number of neighboring domains (starting at 0, -1 if no neighbourg)
      !!                    kexti  : number of columns for extra outer halo
      !!                    kextj  : number of rows for extra outer halo
      !!----------------------------------------------------------------------
      CHARACTER(len=*)                                        , INTENT(in   ) ::   cdname      ! name of the calling subroutine
      REAL(sp), DIMENSION(1-kexti:jpi+kexti,1-kextj:jpj+kextj), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)                                        , INTENT(in   ) ::   cd_type  ! nature of ptab array grid-points
      REAL(sp)                                         , INTENT(in   ) ::   psgn     ! sign used across the north fold
      INTEGER                                                 , INTENT(in   ) ::   kexti    ! extra i-halo width
      INTEGER                                                 , INTENT(in   ) ::   kextj    ! extra j-halo width
      !
      INTEGER  ::   jl   ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom        ! local integers
      INTEGER  ::   ipreci, iprecj             !   -       -
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for mpi_isend
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for mpi_isend
      !!
      REAL(sp), DIMENSION(1-kexti:jpi+kexti,nn_hls+kextj,2) ::   r2dns, r2dsn
      REAL(sp), DIMENSION(1-kextj:jpj+kextj,nn_hls+kexti,2) ::   r2dwe, r2dew
      !!----------------------------------------------------------------------
      ipreci = nn_hls + kexti      ! take into account outer extra 2D overlap area
      iprecj = nn_hls + kextj

      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, 1, 1, 1, ld_lbc = .TRUE. )

      ! 1. standard boundary treatment
      ! ------------------------------
      ! Order matters Here !!!!
      !
      !                                      ! East-West boundaries
      !                                           !* Cyclic east-west
      IF( l_Iperio ) THEN
         pt2d(1-kexti:     1   ,:) = pt2d(jpi-1-kexti: jpi-1 ,:)       ! east
         pt2d(  jpi  :jpi+kexti,:) = pt2d(     2     :2+kexti,:)       ! west
         !
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(  1-kexti   :nn_hls   ,:) = 0._sp    ! east except at F-point
                                      pt2d(jpi-nn_hls+1:jpi+kexti,:) = 0._sp    ! west
      ENDIF
      !                                      ! North-South boundaries
      IF( l_Jperio ) THEN                         !* cyclic (only with no mpp j-split)
         pt2d(:,1-kextj:     1   ) = pt2d(:,jpj-1-kextj:  jpj-1)       ! north
         pt2d(:,  jpj  :jpj+kextj) = pt2d(:,     2     :2+kextj)       ! south
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(:,  1-kextj   :nn_hls   ) = 0._sp    ! north except at F-point
                                      pt2d(:,jpj-nn_hls+1:jpj+kextj) = 0._sp    ! south
      ENDIF
      !

      ! north fold treatment
      ! -----------------------
      IF( l_IdoNFold ) THEN
         !
         SELECT CASE ( jpni )
                   CASE ( 1 )     ;   CALL lbc_nfd         ( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
                   CASE DEFAULT   ;   CALL mpp_lbc_north_icb_sp        ( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
         END SELECT
         !
      ENDIF

      ! 2. East and west directions exchange
      ! ------------------------------------
      ! we play with the neigbours AND the row number because of the periodicity
      !
      IF( mpinei(jpwe) >= 0 .OR. mpinei(jpea) >= 0 ) THEN   ! Read Dirichlet lateral conditions: all exept 2 (i.e. close case)
         iihom = jpi - (2 * nn_hls) -kexti
         DO jl = 1, ipreci
            r2dew(:,jl,1) = pt2d(nn_hls+jl,:)
            r2dwe(:,jl,1) = pt2d(iihom +jl,:)
         END DO
      ENDIF
      !
      !                           ! Migrations
      imigr = ipreci * ( jpj + 2*kextj )
      !
      !                           ! Migrations
      imigr = ipreci * ( jpj + 2*kextj )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !
      IF( mpinei(jpwe) >= 0  )   CALL mppsend_sp( 1, r2dew(1-kextj,1,1), imigr, mpinei(jpwe), ml_req1 )
      IF( mpinei(jpea) >= 0  )   CALL mppsend_sp( 2, r2dwe(1-kextj,1,1), imigr, mpinei(jpea), ml_req2 )
      IF( mpinei(jpwe) >= 0  )   CALL mpprecv_sp( 2, r2dwe(1-kextj,1,2), imigr, mpinei(jpwe) )
      IF( mpinei(jpea) >= 0  )   CALL mpprecv_sp( 1, r2dew(1-kextj,1,2), imigr, mpinei(jpea) )
      IF( mpinei(jpwe) >= 0  )   CALL mpi_wait(ml_req1,ml_stat,ml_err)
      IF( mpinei(jpea) >= 0  )   CALL mpi_wait(ml_req2,ml_stat,ml_err)
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      !                           ! Write Dirichlet lateral conditions
      iihom = jpi - nn_hls
      IF( mpinei(jpwe) >= 0  ) THEN
         DO jl = 1, ipreci
            pt2d(jl-kexti,:) = r2dwe(:,jl,2)
         END DO
      ENDIF
      IF( mpinei(jpea) >= 0  ) THEN
         DO jl = 1, ipreci
            pt2d(iihom+jl,:) = r2dew(:,jl,2)
         END DO
      ENDIF

      ! 3. North and south directions
      ! -----------------------------
      ! always closed : we play only with the neigbours
      !
      IF( mpinei(jpso) >= 0 .OR. mpinei(jpno) >= 0 ) THEN   ! Read Dirichlet lateral conditions: all exept 2 (i.e. close case)
         ijhom = jpj - (2 * nn_hls) - kextj
         DO jl = 1, iprecj
            r2dsn(:,jl,1) = pt2d(:,ijhom +jl)
            r2dns(:,jl,1) = pt2d(:,nn_hls+jl)
         END DO
      ENDIF
      !
      !                           ! Migrations
      imigr = iprecj * ( jpi + 2*kexti )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !
      IF( mpinei(jpso) >= 0  )   CALL mppsend_sp( 3, r2dns(1-kexti,1,1), imigr, mpinei(jpso), ml_req1 )
      IF( mpinei(jpno) >= 0  )   CALL mppsend_sp( 4, r2dsn(1-kexti,1,1), imigr, mpinei(jpno), ml_req2 )
      IF( mpinei(jpso) >= 0  )   CALL mpprecv_sp( 4, r2dsn(1-kexti,1,2), imigr, mpinei(jpso) )
      IF( mpinei(jpno) >= 0  )   CALL mpprecv_sp( 3, r2dns(1-kexti,1,2), imigr, mpinei(jpno) )
      IF( mpinei(jpso) >= 0  )   CALL mpi_wait(ml_req1,ml_stat,ml_err)
      IF( mpinei(jpno) >= 0  )   CALL mpi_wait(ml_req2,ml_stat,ml_err)
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      !                           ! Write Dirichlet lateral conditions
      ijhom = jpj - nn_hls
      !
      IF( mpinei(jpso) >= 0  ) THEN
         DO jl = 1, iprecj
            pt2d(:,jl-kextj) = r2dsn(:,jl,2)
         END DO
      ENDIF
       IF( mpinei(jpno) >= 0  ) THEN
        DO jl = 1, iprecj
            pt2d(:,ijhom+jl) = r2dns(:,jl,2)
         END DO
      ENDIF
      !
   END SUBROUTINE mpp_lnk_2d_icb_sp


   SUBROUTINE mpp_lnk_2d_icb_dp( cdname, pt2d, cd_type, psgn, kexti, kextj )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_icb  ***
      !!
      !! ** Purpose :   Message passing management for 2d array (with extra halo for icebergs)
      !!                This routine receives a (1-kexti:jpi+kexti,1-kexti:jpj+kextj)
      !!                array (usually (0:jpi+1, 0:jpj+1)) from lbc_lnk_icb calls.
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    jpi    : first dimension of the local subdomain
      !!                    jpj    : second dimension of the local subdomain
      !!                    mpinei : number of neighboring domains (starting at 0, -1 if no neighbourg)
      !!                    kexti  : number of columns for extra outer halo
      !!                    kextj  : number of rows for extra outer halo
      !!----------------------------------------------------------------------
      CHARACTER(len=*)                                        , INTENT(in   ) ::   cdname      ! name of the calling subroutine
      REAL(dp), DIMENSION(1-kexti:jpi+kexti,1-kextj:jpj+kextj), INTENT(inout) ::   pt2d     ! 2D array with extra halo
      CHARACTER(len=1)                                        , INTENT(in   ) ::   cd_type  ! nature of ptab array grid-points
      REAL(dp)                                         , INTENT(in   ) ::   psgn     ! sign used across the north fold
      INTEGER                                                 , INTENT(in   ) ::   kexti    ! extra i-halo width
      INTEGER                                                 , INTENT(in   ) ::   kextj    ! extra j-halo width
      !
      INTEGER  ::   jl   ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom        ! local integers
      INTEGER  ::   ipreci, iprecj             !   -       -
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for mpi_isend
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for mpi_isend
      !!
      REAL(dp), DIMENSION(1-kexti:jpi+kexti,nn_hls+kextj,2) ::   r2dns, r2dsn
      REAL(dp), DIMENSION(1-kextj:jpj+kextj,nn_hls+kexti,2) ::   r2dwe, r2dew
      !!----------------------------------------------------------------------
      ipreci = nn_hls + kexti      ! take into account outer extra 2D overlap area
      iprecj = nn_hls + kextj

      IF( narea == 1 .AND. numcom == -1 ) CALL mpp_report( cdname, 1, 1, 1, ld_lbc = .TRUE. )

      ! 1. standard boundary treatment
      ! ------------------------------
      ! Order matters Here !!!!
      !
      !                                      ! East-West boundaries
      !                                           !* Cyclic east-west
      IF( l_Iperio ) THEN
         pt2d(1-kexti:     1   ,:) = pt2d(jpi-1-kexti: jpi-1 ,:)       ! east
         pt2d(  jpi  :jpi+kexti,:) = pt2d(     2     :2+kexti,:)       ! west
         !
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(  1-kexti   :nn_hls   ,:) = 0._dp    ! east except at F-point
                                      pt2d(jpi-nn_hls+1:jpi+kexti,:) = 0._dp    ! west
      ENDIF
      !                                      ! North-South boundaries
      IF( l_Jperio ) THEN                         !* cyclic (only with no mpp j-split)
         pt2d(:,1-kextj:     1   ) = pt2d(:,jpj-1-kextj:  jpj-1)       ! north
         pt2d(:,  jpj  :jpj+kextj) = pt2d(:,     2     :2+kextj)       ! south
      ELSE                                        !* closed
         IF( .NOT. cd_type == 'F' )   pt2d(:,  1-kextj   :nn_hls   ) = 0._dp    ! north except at F-point
                                      pt2d(:,jpj-nn_hls+1:jpj+kextj) = 0._dp    ! south
      ENDIF
      !

      ! north fold treatment
      ! -----------------------
      IF( l_IdoNFold ) THEN
         !
         SELECT CASE ( jpni )
                   CASE ( 1 )     ;   CALL lbc_nfd         ( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
                   CASE DEFAULT   ;   CALL mpp_lbc_north_icb_dp        ( pt2d(1:jpi,1:jpj+kextj), cd_type, psgn, kextj )
         END SELECT
         !
      ENDIF

      ! 2. East and west directions exchange
      ! ------------------------------------
      ! we play with the neigbours AND the row number because of the periodicity
      !
      IF( mpinei(jpwe) >= 0 .OR. mpinei(jpea) >= 0 ) THEN   ! Read Dirichlet lateral conditions: all exept 2 (i.e. close case)
         iihom = jpi - (2 * nn_hls) -kexti
         DO jl = 1, ipreci
            r2dew(:,jl,1) = pt2d(nn_hls+jl,:)
            r2dwe(:,jl,1) = pt2d(iihom +jl,:)
         END DO
      ENDIF
      !
      !                           ! Migrations
      imigr = ipreci * ( jpj + 2*kextj )
      !
      !                           ! Migrations
      imigr = ipreci * ( jpj + 2*kextj )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !
      IF( mpinei(jpwe) >= 0  )   CALL mppsend_dp( 1, r2dew(1-kextj,1,1), imigr, mpinei(jpwe), ml_req1 )
      IF( mpinei(jpea) >= 0  )   CALL mppsend_dp( 2, r2dwe(1-kextj,1,1), imigr, mpinei(jpea), ml_req2 )
      IF( mpinei(jpwe) >= 0  )   CALL mpprecv_dp( 2, r2dwe(1-kextj,1,2), imigr, mpinei(jpwe) )
      IF( mpinei(jpea) >= 0  )   CALL mpprecv_dp( 1, r2dew(1-kextj,1,2), imigr, mpinei(jpea) )
      IF( mpinei(jpwe) >= 0  )   CALL mpi_wait(ml_req1,ml_stat,ml_err)
      IF( mpinei(jpea) >= 0  )   CALL mpi_wait(ml_req2,ml_stat,ml_err)
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      !                           ! Write Dirichlet lateral conditions
      iihom = jpi - nn_hls
      IF( mpinei(jpwe) >= 0  ) THEN
         DO jl = 1, ipreci
            pt2d(jl-kexti,:) = r2dwe(:,jl,2)
         END DO
      ENDIF
      IF( mpinei(jpea) >= 0  ) THEN
         DO jl = 1, ipreci
            pt2d(iihom+jl,:) = r2dew(:,jl,2)
         END DO
      ENDIF

      ! 3. North and south directions
      ! -----------------------------
      ! always closed : we play only with the neigbours
      !
      IF( mpinei(jpso) >= 0 .OR. mpinei(jpno) >= 0 ) THEN   ! Read Dirichlet lateral conditions: all exept 2 (i.e. close case)
         ijhom = jpj - (2 * nn_hls) - kextj
         DO jl = 1, iprecj
            r2dsn(:,jl,1) = pt2d(:,ijhom +jl)
            r2dns(:,jl,1) = pt2d(:,nn_hls+jl)
         END DO
      ENDIF
      !
      !                           ! Migrations
      imigr = iprecj * ( jpi + 2*kexti )
      !
      IF( ln_timing ) CALL tic_tac(.TRUE.)
      !
      IF( mpinei(jpso) >= 0  )   CALL mppsend_dp( 3, r2dns(1-kexti,1,1), imigr, mpinei(jpso), ml_req1 )
      IF( mpinei(jpno) >= 0  )   CALL mppsend_dp( 4, r2dsn(1-kexti,1,1), imigr, mpinei(jpno), ml_req2 )
      IF( mpinei(jpso) >= 0  )   CALL mpprecv_dp( 4, r2dsn(1-kexti,1,2), imigr, mpinei(jpso) )
      IF( mpinei(jpno) >= 0  )   CALL mpprecv_dp( 3, r2dns(1-kexti,1,2), imigr, mpinei(jpno) )
      IF( mpinei(jpso) >= 0  )   CALL mpi_wait(ml_req1,ml_stat,ml_err)
      IF( mpinei(jpno) >= 0  )   CALL mpi_wait(ml_req2,ml_stat,ml_err)
      !
      IF( ln_timing ) CALL tic_tac(.FALSE.)
      !
      !                           ! Write Dirichlet lateral conditions
      ijhom = jpj - nn_hls
      !
      IF( mpinei(jpso) >= 0  ) THEN
         DO jl = 1, iprecj
            pt2d(:,jl-kextj) = r2dsn(:,jl,2)
         END DO
      ENDIF
       IF( mpinei(jpno) >= 0  ) THEN
        DO jl = 1, iprecj
            pt2d(:,ijhom+jl) = r2dns(:,jl,2)
         END DO
      ENDIF
      !
   END SUBROUTINE mpp_lnk_2d_icb_dp


END MODULE lbclnk

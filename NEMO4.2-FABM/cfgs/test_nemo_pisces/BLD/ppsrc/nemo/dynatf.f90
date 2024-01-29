










MODULE dynatf
   !!=========================================================================
   !!                       ***  MODULE  dynatf  ***
   !! Ocean dynamics: time filtering
   !!=========================================================================
   !! History :  OPA  !  1987-02  (P. Andrich, D. L Hostis)  Original code
   !!                 !  1990-10  (C. Levy, G. Madec)
   !!            7.0  !  1993-03  (M. Guyon)  symetrical conditions
   !!            8.0  !  1997-02  (G. Madec & M. Imbard)  opa, release 8.0
   !!            8.2  !  1997-04  (A. Weaver)  Euler forward step
   !!             -   !  1997-06  (G. Madec)  lateral boudary cond., lbc routine
   !!    NEMO    1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-10  (C. Talandier, A-M. Treguier) Open boundary cond.
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.3  !  2007-07  (D. Storkey) Calls to BDY routines.
   !!            3.2  !  2009-06  (G. Madec, R.Benshila)  re-introduce the vvl option
   !!            3.3  !  2010-09  (D. Storkey, E.O'Dea) Bug fix for BDY module
   !!            3.3  !  2011-03  (P. Oddo) Bug fix for time-splitting+(BDY-OBC) and not VVL
   !!            3.5  !  2013-07  (J. Chanut) Compliant with time splitting changes
   !!            3.6  !  2014-04  (G. Madec) add the diagnostic of the time filter trends
   !!            3.7  !  2015-11  (J. Chanut) Free surface simplification
   !!            4.1  !  2019-08  (A. Coward, D. Storkey) Rename dynnxt.F90 -> dynatf.F90. Now just does time filtering.
   !!-------------------------------------------------------------------------

   !!----------------------------------------------------------------------------------------------
   !!   dyn_atf       : apply Asselin time filtering to "now" velocities and vertical scale factors
   !!----------------------------------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbcrnf         ! river runoffs
   USE phycst         ! physical constants
   USE dynadv         ! dynamics: vector invariant versus flux form
   USE dynspg_ts      ! surface pressure gradient: split-explicit scheme
   USE domvvl         ! variable volume
   USE bdy_oce , ONLY : ln_bdy
   USE bdydta         ! ocean open boundary conditions
   USE bdydyn         ! ocean open boundary conditions
   USE bdyvol         ! ocean open boundary condition (bdy_vol routines)
   USE trd_oce        ! trends: ocean variables
   USE trddyn         ! trend manager: dynamics
   USE trdken         ! trend manager: kinetic energy
   USE isf_oce   , ONLY: ln_isf     ! ice shelf
   USE isfdynatf , ONLY: isf_dynatf ! ice shelf volume filter correction subroutine
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lbclnk         ! lateral boundary condition (or mpp link)
   USE lib_mpp        ! MPP library
   USE prtctl         ! Print control
   USE timing         ! Timing
   USE zdfdrg ,  ONLY : ln_drgice_imp, rCdU_top

   IMPLICIT NONE
   PRIVATE

   PUBLIC    dyn_atf   ! routine called by step.F90

   !!----------------------------------------------------------------------
   !!   'key_qco'                        Quasi-Eulerian vertical coordinate
   !!       OR         EMPTY MODULE
   !!   'key_linssh'                        Fix in time vertical coordinate
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_atf( kt, Kbb, Kmm, Kaa, puu, pvv, pe3t, pe3u, pe3v )
      INTEGER                             , INTENT(in   ) :: kt               ! ocean time-step index
      INTEGER                             , INTENT(in   ) :: Kbb, Kmm, Kaa    ! before and after time level indices
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) :: puu, pvv         ! velocities to be time filtered
      REAL(wp), DIMENSION(jpi,jpj,jpk,jpt), INTENT(inout) :: pe3t, pe3u, pe3v ! scale factors to be time filtered

      WRITE(*,*) 'dyn_atf: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_atf


   !!=========================================================================
END MODULE dynatf












MODULE zdfswm
   !!======================================================================
   !!                       ***  MODULE  zdfswm  ***
   !! vertical physics :   surface wave-induced mixing 
   !!======================================================================
   !! History :  3.6  !  2014-10  (E. Clementi)  Original code
   !!            4.0  !  2017-04  (G. Madec)  debug + simplifications
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_swm      : update Kz due to surface wave-induced mixing 
   !!   zdf_swm_init : initilisation
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean domain variable
   USE zdf_oce        ! vertical physics: mixing coefficients
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE sbcwave        ! wave module
   !
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)  
   USE lib_mpp        ! distribued memory computing library
   
   IMPLICIT NONE
   PRIVATE

   PUBLIC zdf_swm         ! routine called in zdp_phy
   PUBLIC zdf_swm_init    ! routine called in zdf_phy_init

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
   !! $Id: zdfswm.F90 15062 2021-06-28 11:19:48Z jchanut $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE zdf_swm( kt, Kmm, p_avm, p_avt, p_avs )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE zdf_swm ***
      !!
      !! ** Purpose :Compute the swm term (qbv) to be added to
      !!             vertical viscosity and diffusivity coeffs.  
      !!
      !! ** Method  :   Compute the swm term Bv (zqb) and added it to
      !!               vertical viscosity and diffusivity coefficients
      !!                   zqb = alpha * A * Us(0) * exp (3 * k * z)
      !!               where alpha is set here to 1
      !!             
      !! ** action  :    avt, avs, avm updated by the surface wave-induced mixing
      !!                               (inner domain only)
      !!               
      !! reference : Qiao et al. GRL, 2004
      !!---------------------------------------------------------------------
      INTEGER                    , INTENT(in   ) ::   kt             ! ocean time step
      INTEGER                    , INTENT(in   ) ::   Kmm            ! time level index
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avm          ! momentum Kz (w-points)
      REAL(wp), DIMENSION(:,:,:) , INTENT(inout) ::   p_avt, p_avs   ! tracer   Kz (w-points)
      !
      INTEGER ::   ji, jj, jk   ! dummy loop indices
      REAL(wp)::   zcoef, zqb   ! local scalar
      !!---------------------------------------------------------------------
      !
      zcoef = 1._wp * 0.353553_wp
      DO jk =  2,  jpkm1  ; DO jj = ntsj-(   nn_hls-1-(  nn_hls-1+  nn_hls-1)*nthb), ntej+(   nn_hls-1-(  nn_hls-1+  nn_hls-1)*ntht) ; DO ji = ntsi-( nn_hls-1-( nn_hls-1+  nn_hls-1)*nthl), ntei+(   nn_hls-1-(  nn_hls-1+ nn_hls-1)*nthr)
         zqb = zcoef * hsw(ji,jj) * tsd2d(ji,jj) * EXP( -3. * wnum(ji,jj) * gdepw_0(ji,jj,jk) ) * wmask(ji,jj,jk)
         !
         p_avt(ji,jj,jk) = p_avt(ji,jj,jk) + zqb
         p_avs(ji,jj,jk) = p_avs(ji,jj,jk) + zqb
         p_avm(ji,jj,jk) = p_avm(ji,jj,jk) + zqb
      END DO   ;   END DO   ;   END DO
      !
   END SUBROUTINE zdf_swm
   
   
   SUBROUTINE zdf_swm_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE zdf_swm_init ***
      !!
      !! ** Purpose :   surface wave-induced mixing initialisation  
      !!
      !! ** Method  :   check the availability of surface wave fields
      !!---------------------------------------------------------------------
      !
      IF(lwp) THEN                  ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'zdf_swm_init : surface wave-driven mixing'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      IF(  .NOT.ln_wave .OR.   &
         & .NOT.ln_sdw    )   CALL ctl_stop ( 'zdf_swm_init: ln_zdfswm=T but ln_wave and ln_sdw /= T')
      !
   END SUBROUTINE zdf_swm_init

   !!======================================================================
END MODULE zdfswm

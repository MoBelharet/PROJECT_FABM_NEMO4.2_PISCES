










MODULE domvvl
   !!======================================================================
   !!                       ***  MODULE domvvl   ***
   !! Ocean :
   !!======================================================================
   !! History :  2.0  !  2006-06  (B. Levier, L. Marie)  original code
   !!            3.1  !  2009-02  (G. Madec, M. Leclair, R. Benshila)  pure z* coordinate
   !!            3.3  !  2011-10  (M. Leclair) totally rewrote domvvl: vvl option includes z_star and z_tilde coordinates
   !!            3.6  !  2014-11  (P. Mathiot) add ice shelf capability
   !!            4.1  !  2019-08  (A. Coward, D. Storkey) rename dom_vvl_sf_swp -> dom_vvl_sf_update for new timestepping
   !!             -   !  2020-02  (G. Madec, S. Techene) introduce ssh to h0 ratio
   !!----------------------------------------------------------------------

   USE oce             ! ocean dynamics and tracers
   USE phycst          ! physical constant
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! ocean surface boundary condition
   USE wet_dry         ! wetting and drying
   USE usrdef_istate   ! user defined initial state (wad only)
   USE restart         ! ocean restart
   !
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O manager library
   USE lib_mpp         ! distributed memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   !                                                      !!* Namelist nam_vvl
   LOGICAL , PUBLIC :: ln_vvl_zstar           = .FALSE.    ! zstar  vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_ztilde          = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_layer           = .FALSE.    ! level  vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_ztilde_as_zstar = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_zstar_at_eqtor  = .FALSE.    ! ztilde vertical coordinate
   LOGICAL , PUBLIC :: ln_vvl_kepe            = .FALSE.    ! kinetic/potential energy transfer
   !
   INTEGER          :: nn_vvl_interp = 0                   ! scale factors anomaly interpolation method at U-V-F points
                                                           ! =0 linear with no bottom correction over steps (old)
                                                           ! =1 linear with bottom correction over steps
                                                           ! =2 "qco like", i.e. proportional to thicknesses at rest
   !
   !                                                       ! conservation: not used yet
   REAL(wp)         :: rn_ahe3                             ! thickness diffusion coefficient
   REAL(wp)         :: rn_rst_e3t                          ! ztilde to zstar restoration timescale [days]
   REAL(wp)         :: rn_lf_cutoff                        ! cutoff frequency for low-pass filter  [days]
   REAL(wp)         :: rn_zdef_max                         ! maximum fractional e3t deformation
   LOGICAL , PUBLIC :: ln_vvl_dbg = .FALSE.                ! debug control prints

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: un_td, vn_td                ! thickness diffusion transport
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: hdiv_lf                     ! low frequency part of hz divergence
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tilde_e3t_b, tilde_e3t_n    ! baroclinic scale factors
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: tilde_e3t_a, dtilde_e3t_a   ! baroclinic scale factors
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: frq_rst_e3t                 ! retoring period for scale factors
   REAL(wp)        , ALLOCATABLE, SAVE, DIMENSION(:,:)   :: frq_rst_hdv                 ! retoring period for low freq. divergence

   !!----------------------------------------------------------------------
   !!   'key_qco'                        Quasi-Eulerian vertical coordinate
   !!       OR         EMPTY MODULE
   !!   'key_linssh'                        Fix in time vertical coordinate
   !!----------------------------------------------------------------------

   !!======================================================================
END MODULE domvvl

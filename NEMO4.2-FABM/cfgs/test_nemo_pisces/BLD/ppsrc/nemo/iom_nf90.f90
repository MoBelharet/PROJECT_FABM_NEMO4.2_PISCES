










MODULE iom_nf90
   !!======================================================================
   !!                    ***  MODULE  iom_nf90 ***
   !! Input/Output manager :  Library to read input files with NF90 (only fliocom module)
   !!======================================================================
   !! History :  9.0  ! 05 12  (J. Belier) Original code
   !!            9.0  ! 06 02  (S. Masson) Adaptation to NEMO
   !!             "   ! 07 07  (D. Storkey) Changes to iom_nf90_gettime
   !!            3.6  ! 2015-15  (J. Harle) Added procedure to read REAL attributes
   !!            4.0  ! 2017-11 (M. Andrejczuk) Extend IOM interface to write any 3D fields
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   iom_open       : open a file read only
   !!   iom_close      : close a file or all files opened by iom
   !!   iom_get        : read a field (interfaced to several routines)
   !!   iom_varid      : get the id of a variable in a file
   !!   iom_rstput     : write a field in a restart file (interfaced to several routines)
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce, ONLY: ght_abl ! abl vertical level number and height
   USE lbclnk          ! lateal boundary condition / mpp exchanges
   USE iom_def         ! iom variables definitions
   USE netcdf          ! NetCDF library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC iom_nf90_open  , iom_nf90_close, iom_nf90_varid, iom_nf90_get, iom_nf90_rstput
   PUBLIC iom_nf90_chkatt, iom_nf90_getatt, iom_nf90_putatt
   PUBLIC iom_nf90_check

   INTERFACE iom_nf90_get
      MODULE PROCEDURE iom_nf90_g0d_sp
      MODULE PROCEDURE iom_nf90_g0d_dp, iom_nf90_g123d
   END INTERFACE

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: iom_nf90.F90 14433 2021-02-11 08:06:49Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE iom_nf90_open( cdname, kiomid, ldwrt, ldok, kdlev, cdcomp )
      !!---------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_open  ***
      !!
      !! ** Purpose : open an input file with NF90
      !!---------------------------------------------------------------------
      CHARACTER(len=*)       , INTENT(inout)           ::   cdname      ! File name
      INTEGER                , INTENT(  out)           ::   kiomid      ! nf90 identifier of the opened file
      LOGICAL                , INTENT(in   )           ::   ldwrt       ! read or write the file?
      LOGICAL                , INTENT(in   )           ::   ldok        ! check the existence
      INTEGER                , INTENT(in   ), OPTIONAL ::   kdlev       ! size of the ice/abl third dimension
      CHARACTER(len=3)       , INTENT(in   ), OPTIONAL ::   cdcomp      ! name of component calling iom_nf90_open

      CHARACTER(LEN=256) ::   clinfo           ! info character
      CHARACTER(LEN=256) ::   cltmp            ! temporary character
      CHARACTER(LEN=12 ) ::   clfmt            ! writing format
      CHARACTER(LEN=3  ) ::   clcomp           ! name of component calling iom_nf90_open
      INTEGER            ::   idg              ! number of digits
      INTEGER            ::   iln              ! lengths of character
      INTEGER            ::   istop            ! temporary storage of nstop
      INTEGER            ::   if90id           ! nf90 identifier of the opened file
      INTEGER            ::   idmy             ! dummy variable
      INTEGER            ::   jl               ! loop variable
      INTEGER            ::   ichunk           ! temporary storage of nn_chunksz
      INTEGER            ::   imode            ! creation mode flag: NF90_CLOBBER or NF90_NOCLOBBER or NF90_HDF5
      INTEGER            ::   ihdf5            ! local variable for retrieval of value for NF90_HDF5
      LOGICAL            ::   llclobber        ! local definition of ln_clobber
      !---------------------------------------------------------------------
      !
      clinfo = '                    iom_nf90_open ~~~  '
      istop = nstop     ! store the actual value of nstop
      !
      !                 !number of vertical levels
      IF( PRESENT(cdcomp) )   THEN
         IF( .NOT. PRESENT(kdlev) ) CALL ctl_stop( 'iom_nf90_open: cdcomp and kdlev must both be present' )
         clcomp = cdcomp    ! use input value
      ELSE
         clcomp = 'OCE'     ! by default
      ENDIF
      !
      IF( nn_chunksz > 0 ) THEN   ;   ichunk = nn_chunksz
      ELSE                        ;   ichunk = NF90_SIZEHINT_DEFAULT
      ENDIF
      !
      llclobber = ldwrt .AND. ln_clobber
      IF( ldok .AND. .NOT. llclobber ) THEN      !==  Open existing file ==!
         !                                       !=========================!
         IF( ldwrt ) THEN  ! ... in write mode
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' open existing file: '//TRIM(cdname)//' in WRITE mode'
            IF( snc4set%luse ) THEN
               CALL iom_nf90_check(NF90_OPEN( TRIM(cdname), NF90_WRITE  , if90id ), clinfo)
            ELSE
               CALL iom_nf90_check(NF90_OPEN( TRIM(cdname), NF90_WRITE  , if90id, chunksize = ichunk ), clinfo)
            ENDIF
            CALL iom_nf90_check(NF90_SET_FILL( if90id, NF90_NOFILL, idmy                          ), clinfo)
         ELSE              ! ... in read mode
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' open existing file: '//TRIM(cdname)//' in READ mode'
            CALL iom_nf90_check(NF90_OPEN( TRIM(cdname), NF90_NOWRITE, if90id, chunksize = ichunk ), clinfo)
         ENDIF
      ELSE                                       !== the file doesn't exist ==!   (or we overwrite it)
         !                                       !============================!
         iln = INDEX( cdname, '.nc' )
         IF( ldwrt ) THEN              !* the file should be open in write mode so we create it...
            IF( jpnij > 1 ) THEN
               idg = MAX( INT(LOG10(REAL(MAX(1,jpnij-1),wp))) + 1, 4 )          ! how many digits to we need to write? min=4, max=9
               WRITE(clfmt, "('(a,a,i', i1, '.', i1, ',a)')") idg, idg          ! '(a,a,ix.x,a)'
               WRITE(cltmp,clfmt) cdname(1:iln-1), '_', narea-1, '.nc'
               cdname = TRIM(cltmp)
            ENDIF
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' create new file: '//TRIM(cdname)//' in WRITE mode'

            IF( llclobber ) THEN   ;   imode = IOR( NF90_64BIT_OFFSET, NF90_CLOBBER   )
            ELSE                   ;   imode = IOR( NF90_64BIT_OFFSET, NF90_NOCLOBBER )
            ENDIF
            IF( snc4set%luse ) THEN
               IF(lwp) WRITE(numout,*) TRIM(clinfo)//' creating file: '//TRIM(cdname)//' in hdf5 (netcdf4) mode'
               CALL GET_NF90_SYMBOL("NF90_HDF5", ihdf5)
               IF( llclobber ) THEN   ;   imode = IOR(ihdf5, NF90_CLOBBER)
               ELSE                   ;   imode = IOR(ihdf5, NF90_NOCLOBBER)
               ENDIF
               CALL iom_nf90_check(NF90_CREATE( TRIM(cdname), imode, if90id ), clinfo)
            ELSE
               CALL iom_nf90_check(NF90_CREATE( TRIM(cdname), imode, if90id, chunksize = ichunk ), clinfo)
            ENDIF
            CALL iom_nf90_check(NF90_SET_FILL( if90id, NF90_NOFILL,                   idmy ), clinfo)
            ! define dimensions
                               CALL iom_nf90_check(NF90_DEF_DIM( if90id,            'x',  Ni_0, idmy ), clinfo)
                               CALL iom_nf90_check(NF90_DEF_DIM( if90id,            'y',  Nj_0, idmy ), clinfo)
            SELECT CASE (clcomp)
            CASE ('OCE')   ;   CALL iom_nf90_check(NF90_DEF_DIM( if90id,      'nav_lev',   jpk, idmy ), clinfo)
            CASE ('ICE')   ;   CALL iom_nf90_check(NF90_DEF_DIM( if90id,       'numcat', kdlev, idmy ), clinfo)
            CASE ('ABL')   ;   CALL iom_nf90_check(NF90_DEF_DIM( if90id,      'nav_lev', kdlev, idmy ), clinfo)
            CASE ('SED')   ;   CALL iom_nf90_check(NF90_DEF_DIM( if90id,       'numsed', kdlev, idmy ), clinfo)
            CASE DEFAULT   ;   CALL ctl_stop( 'iom_nf90_open unknown component type' )
            END SELECT
                               CALL iom_nf90_check(NF90_DEF_DIM( if90id, 'time_counter', NF90_UNLIMITED, idmy ), clinfo)
            ! global attributes
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_number_total'   , jpnij                        ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_number'         , narea-1                      ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_dimensions_ids' , (/ 1        , 2           /) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_size_global'    , (/ Ni0glo    , Nj0glo     /) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_size_local'     , (/ Ni_0      , Nj_0       /) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_position_first' , (/ mig0(Nis0), mjg0(Njs0) /) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_position_last'  , (/ mig0(Nie0), mjg0(Nje0) /) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_halo_size_start', (/ 0         , 0          /) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_halo_size_end'  , (/ 0         , 0          /) ), clinfo)
            CALL iom_nf90_check(NF90_PUT_ATT( if90id, NF90_GLOBAL, 'DOMAIN_type'           , 'BOX'                        ), clinfo)
         ELSE                          !* the file should be open for read mode so it must exist...
            CALL ctl_stop( TRIM(clinfo), ' should be impossible case...' )
         ENDIF
      ENDIF
      !
      ! start to fill file informations
      ! =============
      IF( istop == nstop ) THEN   ! no error within this routine
!does not work with some compilers         kiomid = MINLOC(iom_file(:)%nfid, dim = 1)
         kiomid = 0
         DO jl = jpmax_files, 1, -1
            IF( iom_file(jl)%nfid == 0 )   kiomid = jl
         ENDDO
         iom_file(kiomid)%name   = TRIM(cdname)
         iom_file(kiomid)%comp   = clcomp
         iom_file(kiomid)%nfid   = if90id
         iom_file(kiomid)%nvars  = 0
         iom_file(kiomid)%irec   = -1   ! useless for NetCDF files, used to know if the file is in define mode
         CALL iom_nf90_check(NF90_Inquire(if90id, unlimitedDimId = iom_file(kiomid)%iduld), clinfo)
         IF( iom_file(kiomid)%iduld .GE. 0 ) THEN
            CALL iom_nf90_check(NF90_Inquire_Dimension(if90id, iom_file(kiomid)%iduld,    &
               &                                       name = iom_file(kiomid)%uldname,   &
               &                                       len  = iom_file(kiomid)%lenuld ), clinfo )
         ENDIF
         IF(lwp) WRITE(numout,*) '                   ---> '//TRIM(cdname)//' OK'
      ELSE
         kiomid = 0               ! return error flag
      ENDIF
      !
   END SUBROUTINE iom_nf90_open


   SUBROUTINE iom_nf90_close( kiomid )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_nf90_close  ***
      !!
      !! ** Purpose : close an input file with NF90
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kiomid   ! iom identifier of the file to be closed
      CHARACTER(LEN=100)  ::   clinfo   ! info character
      !---------------------------------------------------------------------
      clinfo = '      iom_nf90_close    , file: '//TRIM(iom_file(kiomid)%name)
      CALL iom_nf90_check(NF90_CLOSE(iom_file(kiomid)%nfid), clinfo)
   END SUBROUTINE iom_nf90_close


   FUNCTION iom_nf90_varid ( kiomid, cdvar, kiv, kdimsz, kndims, lduld )
      !!-----------------------------------------------------------------------
      !!                  ***  FUNCTION  iom_varid  ***
      !!
      !! ** Purpose : get the id of a variable in a file with NF90
      !!-----------------------------------------------------------------------
      INTEGER              , INTENT(in   )           ::   kiomid   ! file Identifier
      CHARACTER(len=*)     , INTENT(in   )           ::   cdvar    ! name of the variable
      INTEGER              , INTENT(in   )           ::   kiv   !
      INTEGER, DIMENSION(:), INTENT(  out), OPTIONAL ::   kdimsz   ! size of each dimension
      INTEGER              , INTENT(  out), OPTIONAL ::   kndims   ! number of dimensions
      LOGICAL              , INTENT(  out), OPTIONAL ::   lduld    ! true if the last dimension is unlimited (time)
      !
      INTEGER                        ::   iom_nf90_varid   ! iom variable Id
      INTEGER                        ::   if90id           ! nf90 file identifier
      INTEGER                        ::   ji               ! dummy loop index
      INTEGER                        ::   ivarid           ! NetCDF  variable Id
      INTEGER                        ::   i_nvd            ! number of dimension of the variable
      INTEGER, DIMENSION(jpmax_dims) ::   idimid           ! dimension ids of the variable
      LOGICAL                        ::   llok             ! ok  test
      CHARACTER(LEN=100)             ::   clinfo           ! info character
      !!-----------------------------------------------------------------------
      clinfo = '          iom_nf90_varid, file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(cdvar)
      iom_nf90_varid = 0                    ! default definition
      IF( PRESENT(kdimsz) ) kdimsz(:) = 0   ! default definition
      if90id = iom_file(kiomid)%nfid        ! get back NetCDF file id
      !
      llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr   ! does the variable exist in the file
      IF( llok ) THEN
         iom_nf90_varid = kiv
         iom_file(kiomid)%nvars       = kiv
         iom_file(kiomid)%nvid(kiv)   = ivarid
         iom_file(kiomid)%cn_var(kiv) = TRIM(cdvar)
         CALL iom_nf90_check(NF90_Inquire_Variable(if90id, ivarid, ndims = i_nvd), clinfo)   ! number of dimensions
         iom_file(kiomid)%ndims(kiv)  = i_nvd
         CALL iom_nf90_check(NF90_Inquire_Variable(if90id, ivarid, dimids = idimid(1:i_nvd)), clinfo)   ! dimensions ids
         iom_file(kiomid)%luld(kiv) = .FALSE.   ! default value
         iom_file(kiomid)%dimsz(:,kiv) = 0      ! reset dimsz in case previously used
         DO ji = 1, i_nvd                       ! dimensions size
            CALL iom_nf90_check(NF90_Inquire_Dimension(if90id, idimid(ji), len = iom_file(kiomid)%dimsz(ji,kiv)), clinfo)
            IF( idimid(ji) == iom_file(kiomid)%iduld ) iom_file(kiomid)%luld(kiv) = .TRUE.   ! unlimited dimension?
         END DO
         !---------- Deal with scale_factor and add_offset
         llok = NF90_Inquire_attribute(if90id, ivarid, 'scale_factor') == nf90_noerr
         IF( llok) THEN
            CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, 'scale_factor', iom_file(kiomid)%scf(kiv)), clinfo)
         ELSE
            iom_file(kiomid)%scf(kiv) = 1.
         END IF
         llok = NF90_Inquire_attribute(if90id, ivarid, 'add_offset') == nf90_noerr
         IF( llok ) THEN
            CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, 'add_offset', iom_file(kiomid)%ofs(kiv)), clinfo)
         ELSE
            iom_file(kiomid)%ofs(kiv) = 0.
         END IF
         ! return the simension size
         IF( PRESENT(kdimsz) ) THEN
            IF( i_nvd <= SIZE(kdimsz) ) THEN
               kdimsz(1:i_nvd) = iom_file(kiomid)%dimsz(1:i_nvd,kiv)
            ELSE
               WRITE(ctmp1,*) i_nvd, SIZE(kdimsz)
               CALL ctl_stop( TRIM(clinfo), 'error in kdimsz size'//TRIM(ctmp1) )
            ENDIF
         ENDIF
         IF( PRESENT(kndims) )  kndims = iom_file(kiomid)%ndims(kiv)
         IF( PRESENT( lduld) )  lduld  = iom_file(kiomid)%luld(kiv)
      ELSE
         iom_nf90_varid = -1   !   variable not found, return error code: -1
      ENDIF
      !
   END FUNCTION iom_nf90_varid

   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_nf90_get
   !!----------------------------------------------------------------------

   SUBROUTINE iom_nf90_g0d_sp( kiomid, kvid, pvar, kstart )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_g0d  ***
      !!
      !! ** Purpose : read a scalar with NF90
      !!-----------------------------------------------------------------------
      INTEGER ,               INTENT(in   )            ::   kiomid   ! Identifier of the file
      INTEGER ,               INTENT(in   )            ::   kvid     ! variable id
      REAL(sp),               INTENT(  out)            ::   pvar     ! read field
      INTEGER , DIMENSION(1), INTENT(in   ), OPTIONAL  ::   kstart   ! start position of the reading in each axis
      !
      CHARACTER(LEN=100)      ::   clinfo   ! info character
      !---------------------------------------------------------------------
      clinfo = 'iom_nf90_g0d , file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(iom_file(kiomid)%cn_var(kvid))
      CALL iom_nf90_check(NF90_GET_VAR(iom_file(kiomid)%nfid, iom_file(kiomid)%nvid(kvid), pvar, start = kstart), clinfo )
   END SUBROUTINE iom_nf90_g0d_sp

   SUBROUTINE iom_nf90_g0d_dp( kiomid, kvid, pvar, kstart )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_g0d  ***
      !!
      !! ** Purpose : read a scalar with NF90
      !!-----------------------------------------------------------------------
      INTEGER ,               INTENT(in   )            ::   kiomid   ! Identifier of the file
      INTEGER ,               INTENT(in   )            ::   kvid     ! variable id
      REAL(dp),               INTENT(  out)            ::   pvar     ! read field
      INTEGER , DIMENSION(1), INTENT(in   ), OPTIONAL  ::   kstart   ! start position of the reading in each axis
      !
      CHARACTER(LEN=100)      ::   clinfo   ! info character
      !---------------------------------------------------------------------
      clinfo = 'iom_nf90_g0d , file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(iom_file(kiomid)%cn_var(kvid))
      CALL iom_nf90_check(NF90_GET_VAR(iom_file(kiomid)%nfid, iom_file(kiomid)%nvid(kvid), pvar, start = kstart), clinfo )
   END SUBROUTINE iom_nf90_g0d_dp

   SUBROUTINE iom_nf90_g123d( kiomid, kvid, knbdim, kstart, kcount, kx1, kx2, ky1, ky2,   &
         &                    pvsp1d, pvsp2d, pvsp3d, pvdp1d, pvdp2d, pvdp3d )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_g123d  ***
      !!
      !! ** Purpose : read a 1D/2D/3D variable with NF90
      !!
      !! ** Method : read ONE record at each CALL
      !!-----------------------------------------------------------------------
      INTEGER                    , INTENT(in   )           ::   kiomid    ! iom identifier of the file
      INTEGER                    , INTENT(in   )           ::   kvid      ! Name of the variable
      INTEGER                    , INTENT(in   )           ::   knbdim    ! number of dimensions of the variable
      INTEGER , DIMENSION(:)     , INTENT(in   )           ::   kstart    ! start position of the reading in each axis
      INTEGER , DIMENSION(:)     , INTENT(in   )           ::   kcount    ! number of points to be read in each axis
      INTEGER ,                    INTENT(in   )           ::   kx1, kx2, ky1, ky2   ! subdomain indexes
      REAL(sp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pvsp1d    ! read field (1D case), single precision
      REAL(sp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pvsp2d    ! read field (2D case), single precision
      REAL(sp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pvsp3d    ! read field (3D case), single precision
      REAL(dp), DIMENSION(:)     , INTENT(  out), OPTIONAL ::   pvdp1d    ! read field (1D case), double precision
      REAL(dp), DIMENSION(:,:)   , INTENT(  out), OPTIONAL ::   pvdp2d    ! read field (2D case), double precision
      REAL(dp), DIMENSION(:,:,:) , INTENT(  out), OPTIONAL ::   pvdp3d    ! read field (3D case), double precision
      !
      CHARACTER(LEN=100) ::   clinfo               ! info character
      INTEGER            ::   if90id               ! nf90 identifier of the opened file
      INTEGER            ::   ivid                 ! nf90 variable id
      !---------------------------------------------------------------------
      clinfo = 'iom_nf90_g123d , file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(iom_file(kiomid)%cn_var(kvid))
      if90id = iom_file(kiomid)%nfid         ! get back NetCDF file id
      ivid   = iom_file(kiomid)%nvid(kvid)   ! get back NetCDF var id
      !
      IF(     PRESENT(pvsp1d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pvsp1d(:                ), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ELSEIF( PRESENT(pvdp1d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pvdp1d(:                ), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ELSEIF( PRESENT(pvsp2d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pvsp2d(kx1:kx2,ky1:ky2  ), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ELSEIF( PRESENT(pvdp2d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pvdp2d(kx1:kx2,ky1:ky2  ), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ELSEIF( PRESENT(pvsp3d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pvsp3d(kx1:kx2,ky1:ky2,:), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ELSEIF( PRESENT(pvdp3d) ) THEN
         CALL iom_nf90_check( NF90_GET_VAR(if90id, ivid, pvdp3d(kx1:kx2,ky1:ky2,:), start = kstart(1:knbdim),   &
            &                                                                       count = kcount(1:knbdim)), clinfo )
      ENDIF
      !
   END SUBROUTINE iom_nf90_g123d



   SUBROUTINE iom_nf90_chkatt( kiomid, cdatt, llok, ksize, cdvar )
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_chkatt  ***
      !!
      !! ** Purpose : check existence of attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER         , INTENT(in   ) ::   kiomid   ! Identifier of the file
      CHARACTER(len=*), INTENT(in   ) ::   cdatt    ! attribute name
      LOGICAL         , INTENT(  out) ::   llok     ! error code
      INTEGER         , INTENT(  out), OPTIONAL     &
                      &               ::   ksize    ! attribute size
      CHARACTER(len=*), INTENT(in   ), OPTIONAL     &
                      &               ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   isize    ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      IF( PRESENT(cdvar) ) THEN
         ! check the variable exists in the file
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr
         IF( llok ) &
            ! check the variable has the attribute required
            llok = NF90_Inquire_attribute(if90id, ivarid, cdatt, len=isize ) == nf90_noerr
      ELSE
         llok = NF90_Inquire_attribute(if90id, NF90_GLOBAL, cdatt, len=isize ) == nf90_noerr
      ENDIF
      !
      IF( PRESENT(ksize) ) ksize = isize
      !
      IF( .not. llok) &
         CALL ctl_warn('iom_nf90_chkatt: no attribute '//cdatt//' found')
      !
   END SUBROUTINE iom_nf90_chkatt


   !!----------------------------------------------------------------------
   !!                   INTERFACE iom_nf90_getatt
   !!----------------------------------------------------------------------

   SUBROUTINE iom_nf90_getatt( kiomid, cdatt, katt0d, katt1d, patt0d, patt1d, cdatt0d, cdvar)
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_getatt  ***
      !!
      !! ** Purpose : read an attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER               , INTENT(in   )           ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt    ! attribute name
      INTEGER               , INTENT(  out), OPTIONAL ::   katt0d   ! read scalar integer
      INTEGER, DIMENSION(:) , INTENT(  out), OPTIONAL ::   katt1d   ! read 1d array integer
      REAL(wp)              , INTENT(  out), OPTIONAL ::   patt0d   ! read scalar  real
      REAL(wp), DIMENSION(:), INTENT(  out), OPTIONAL ::   patt1d   ! read 1d array real
      CHARACTER(len=*)      , INTENT(  out), OPTIONAL ::   cdatt0d  ! read character
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      LOGICAL                         ::   llok     ! temporary logical
      CHARACTER(LEN=100)              ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      IF( PRESENT(cdvar) ) THEN
         ! check the variable exists in the file
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr
         IF( llok ) THEN
            ! check the variable has the attribute required
            llok = NF90_Inquire_attribute(if90id, ivarid, cdatt) == nf90_noerr
         ELSE
            CALL ctl_warn('iom_nf90_getatt: no variable '//TRIM(cdvar)//' found')
         ENDIF
      ELSE
         llok = NF90_Inquire_attribute(if90id, NF90_GLOBAL, cdatt) == nf90_noerr
         ivarid = NF90_GLOBAL
      ENDIF
      !
      IF( llok) THEN
         clinfo = 'iom_nf90_getatt, file: '//TRIM(iom_file(kiomid)%name)//', att: '//TRIM(cdatt)
         IF(PRESENT( katt0d))   CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values =  katt0d), clinfo)
         IF(PRESENT( katt1d))   CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values =  katt1d), clinfo)
         IF(PRESENT( patt0d))   CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values =  patt0d), clinfo)
         IF(PRESENT( patt1d))   CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values =  patt1d), clinfo)
         IF(PRESENT(cdatt0d))   CALL iom_nf90_check(NF90_GET_ATT(if90id, ivarid, cdatt, values = cdatt0d), clinfo)
      ELSE
         IF(PRESENT( katt0d))    katt0d    = -999
         IF(PRESENT( katt1d))    katt1d(:) = -999
         IF(PRESENT( patt0d))    patt0d    = -999._wp
         IF(PRESENT( patt1d))    patt1d(:) = -999._wp
         IF(PRESENT(cdatt0d))   cdatt0d    = 'UNKNOWN'
      ENDIF
      !
   END SUBROUTINE iom_nf90_getatt


   SUBROUTINE iom_nf90_putatt( kiomid, cdatt, katt0d, katt1d, patt0d, patt1d, cdatt0d, cdvar)
      !!-----------------------------------------------------------------------
      !!                  ***  ROUTINE  iom_nf90_putatt  ***
      !!
      !! ** Purpose : write an attribute with NF90
      !!              (either a global attribute (default) or a variable
      !!               attribute if optional variable name is supplied (cdvar))
      !!-----------------------------------------------------------------------
      INTEGER               , INTENT(in   )           ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)      , INTENT(in   )           ::   cdatt    ! attribute name
      INTEGER               , INTENT(in   ), OPTIONAL ::   katt0d   ! read scalar integer
      INTEGER, DIMENSION(:) , INTENT(in   ), OPTIONAL ::   katt1d   ! read 1d array integer
      REAL(wp)              , INTENT(in   ), OPTIONAL ::   patt0d   ! read scalar  real
      REAL(wp), DIMENSION(:), INTENT(in   ), OPTIONAL ::   patt1d   ! read 1d array real
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdatt0d  ! read character
      CHARACTER(len=*)      , INTENT(in   ), OPTIONAL ::   cdvar    ! name of the variable
      !
      INTEGER                         ::   if90id   ! temporary integer
      INTEGER                         ::   ivarid   ! NetCDF variable Id
      INTEGER                         ::   isize    ! Attribute size
      INTEGER                         ::   itype    ! Attribute type
      LOGICAL                         ::   llok     ! temporary logical
      LOGICAL                         ::   llatt     ! temporary logical
      LOGICAL                         ::   lldata   ! temporary logical
      CHARACTER(LEN=100)              ::   clinfo   ! info character
      !---------------------------------------------------------------------
      !
      if90id = iom_file(kiomid)%nfid
      IF( PRESENT(cdvar) ) THEN
         llok = NF90_INQ_VARID( if90id, TRIM(cdvar), ivarid ) == nf90_noerr   ! is the variable in the file?
         IF( .NOT. llok ) THEN
            CALL ctl_warn('iom_nf90_putatt: no variable '//TRIM(cdvar)//' found'   &
               &        , '                 no attribute '//cdatt//' written' )
            RETURN
         ENDIF
      ELSE
         ivarid = NF90_GLOBAL
      ENDIF
      llatt = NF90_Inquire_attribute(if90id, ivarid, cdatt, len = isize, xtype = itype ) == nf90_noerr
      !
      ! trick: irec used to know if the file is in define mode or not
      lldata = iom_file(kiomid)%irec /= -1   ! default: go back in define mode if in data mode
      IF( lldata .AND. llatt ) THEN          ! attribute already there. Do we really need to go back in define mode?
         ! do we have the appropriate type?
         IF(PRESENT( katt0d) .OR. PRESENT( katt1d))   llok = itype == NF90_INT
         IF(PRESENT( patt0d) .OR. PRESENT( patt1d))   llok = itype == NF90_DOUBLE
         IF(PRESENT(cdatt0d)                      )   llok = itype == NF90_CHAR
         ! and do we have the appropriate size?
         IF(PRESENT( katt0d))   llok = llok .AND. isize == 1
         IF(PRESENT( katt1d))   llok = llok .AND. isize == SIZE(katt1d)
         IF(PRESENT( patt0d))   llok = llok .AND. isize == 1
         IF(PRESENT( patt1d))   llok = llok .AND. isize == SIZE(patt1d)
         IF(PRESENT(cdatt0d))   llok = llok .AND. isize == LEN_TRIM(cdatt0d)
         !
         lldata = .NOT. llok
      ENDIF
      !
      clinfo = 'iom_nf90_putatt, file: '//TRIM(iom_file(kiomid)%name)//', att: '//TRIM(cdatt)
      IF(lldata)   CALL iom_nf90_check(NF90_REDEF( if90id ), clinfo)   ! leave data mode to define mode
      !
      IF(PRESENT( katt0d))   CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values =       katt0d) , clinfo)
      IF(PRESENT( katt1d))   CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values =       katt1d) , clinfo)
      IF(PRESENT( patt0d))   CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values =       patt0d) , clinfo)
      IF(PRESENT( patt1d))   CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values =       patt1d) , clinfo)
      IF(PRESENT(cdatt0d))   CALL iom_nf90_check(NF90_PUT_ATT(if90id, ivarid, cdatt, values = trim(cdatt0d)), clinfo)
      !
      IF(lldata)   CALL iom_nf90_check(NF90_ENDDEF( if90id ), clinfo)   ! leave define mode to data mode
      !
   END SUBROUTINE iom_nf90_putatt

   SUBROUTINE iom_nf90_rstput( kt, kwrite, kiomid, cdvar , kvid  , ktype ,   &
         &                                 pvsp0d, pvsp1d, pvsp2d, pvsp3d,   &
         &                                 pvdp0d, pvdp1d, pvdp2d, pvdp3d )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE  iom_nf90_rstput  ***
      !!
      !! ** Purpose : read the time axis cdvar in the file
      !!--------------------------------------------------------------------
      INTEGER                     , INTENT(in)           ::   kt       ! ocean time-step
      INTEGER                     , INTENT(in)           ::   kwrite   ! writing time-step
      INTEGER                     , INTENT(in)           ::   kiomid   ! Identifier of the file
      CHARACTER(len=*)            , INTENT(in)           ::   cdvar    ! variable name
      INTEGER                     , INTENT(in)           ::   kvid     ! variable id
      INTEGER                     , INTENT(in), OPTIONAL ::   ktype    ! variable type (default R8)
      REAL(sp)                    , INTENT(in), OPTIONAL ::   pvsp0d   ! written Od field
      REAL(sp), DIMENSION(      :), INTENT(in), OPTIONAL ::   pvsp1d   ! written 1d field
      REAL(sp), DIMENSION(:, :   ), INTENT(in), OPTIONAL ::   pvsp2d   ! written 2d field
      REAL(sp), DIMENSION(:, :, :), INTENT(in), OPTIONAL ::   pvsp3d   ! written 3d field
      REAL(dp)                    , INTENT(in), OPTIONAL ::   pvdp0d   ! written Od field
      REAL(dp), DIMENSION(      :), INTENT(in), OPTIONAL ::   pvdp1d   ! written 1d field
      REAL(dp), DIMENSION(:, :   ), INTENT(in), OPTIONAL ::   pvdp2d   ! written 2d field
      REAL(dp), DIMENSION(:, :, :), INTENT(in), OPTIONAL ::   pvdp3d   ! written 3d field
      !
      INTEGER               :: idims                ! number of dimension
      INTEGER               :: idvar                ! variable id
      INTEGER               :: jd                   ! dimension loop counter
      INTEGER               :: ix1, ix2, iy1, iy2   ! subdomain indexes
      INTEGER, DIMENSION(4) :: idimsz               ! dimensions size
      INTEGER, DIMENSION(4) :: idimid               ! dimensions id
      CHARACTER(LEN=256)    :: clinfo               ! info character
      INTEGER               :: if90id               ! nf90 file identifier
      INTEGER               :: itype                ! variable type
      INTEGER, DIMENSION(4) :: ichunksz             ! NetCDF4 chunk sizes. Will be computed using
      !                                             ! nn_nchunks_[i,j,k,t] namelist parameters
      INTEGER               :: ichunkalg, ishuffle, ideflate, ideflate_level
      !                                             ! NetCDF4 internally fixed parameters
      INTEGER               :: idlv                 ! local variable
      LOGICAL               :: lchunk               ! logical switch to activate chunking and compression
      !                                             ! when appropriate (currently chunking is applied to 4d fields only)
      LOGICAL               :: llis0d, llis1d, llis2d, llis3d
      CHARACTER(LEN=256)    :: ccname               ! local variable
      !---------------------------------------------------------------------
      !
      llis0d = PRESENT(pvsp0d) .OR. PRESENT(pvdp0d)
      llis1d = PRESENT(pvsp1d) .OR. PRESENT(pvdp1d)
      llis2d = PRESENT(pvsp2d) .OR. PRESENT(pvdp2d)
      llis3d = PRESENT(pvsp3d) .OR. PRESENT(pvdp3d)
      !
      clinfo = '          iom_nf90_rp0123d, file: '//TRIM(iom_file(kiomid)%name)//', var: '//TRIM(cdvar)
      if90id = iom_file(kiomid)%nfid
      !
      ! define dimension variables if it is not already done
      ! ==========================
      IF( iom_file(kiomid)%nvars == 0 ) THEN
         ! are we in define mode?
         IF( iom_file(kiomid)%irec /= -1 ) THEN   ! trick: irec used to know if the file is in define mode or not
            CALL iom_nf90_check(NF90_REDEF( if90id ), clinfo)   ;   iom_file(kiomid)%irec = -1
         ENDIF
         ! define the dimension variables if it is not already done
         DO jd = 1, 2
            CALL iom_nf90_check(NF90_INQUIRE_DIMENSION(if90id,jd,iom_file(kiomid)%cn_var(jd),iom_file(kiomid)%dimsz(jd,jd)),clinfo)
            ccname = TRIM(iom_file(kiomid)%cn_var(jd))
            IF ( ccname == 'x') ccname = 'nav_lon'
            IF ( ccname == 'y') ccname = 'nav_lat'
            CALL iom_nf90_check(NF90_DEF_VAR( if90id, ccname, NF90_FLOAT , (/ 1, 2 /),   &
               &                              iom_file(kiomid)%nvid(jd) ), clinfo)
         END DO
         iom_file(kiomid)%dimsz(2,1) = iom_file(kiomid)%dimsz(2,2)   ! second dim of first  variable
         iom_file(kiomid)%dimsz(1,2) = iom_file(kiomid)%dimsz(1,1)   ! first  dim of second variable
         DO jd = 3, 4
            CALL iom_nf90_check(NF90_INQUIRE_DIMENSION(if90id,jd,iom_file(kiomid)%cn_var(jd),iom_file(kiomid)%dimsz(1,jd)), clinfo)
            CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(iom_file(kiomid)%cn_var(jd)), NF90_FLOAT , (/ jd   /),   &
               &                              iom_file(kiomid)%nvid(jd) ), clinfo)
         END DO
         ! update informations structure related the dimension variable we just added...
         iom_file(kiomid)%nvars       = 4
         iom_file(kiomid)%luld(1:4)   = (/ .FALSE., .FALSE., .FALSE., .TRUE. /)
         iom_file(kiomid)%ndims(1:4)  = (/ 2, 2, 1, 1 /)
         IF(lwp) WRITE(numout,*) TRIM(clinfo)//' define dimension variables done'
      ENDIF
      ! define the data if it is not already done
      ! ===============
      IF( kvid <= 0 ) THEN
         !
         ! NetCDF4 chunking and compression fixed settings
         ichunkalg = 0
         ishuffle = 1
         ideflate = 1
         ideflate_level = 1
         !
         idvar = iom_file(kiomid)%nvars + 1
         ! are we in define mode?
         IF( iom_file(kiomid)%irec /= -1 ) THEN   ! trick: irec used to know if the file is in define mode or not
            CALL iom_nf90_check(NF90_REDEF( if90id ), clinfo)   ;   iom_file(kiomid)%irec = -1
         ENDIF
         ! variable definition
         IF(     llis0d ) THEN   ;   idims = 0
         ELSEIF( llis1d ) THEN
                                     idims = 2   ;   idimid(1:idims) = (/3,4/)
         ELSEIF( llis2d ) THEN   ;   idims = 3   ;   idimid(1:idims) = (/1,2,4/)
         ELSEIF( llis3d ) THEN
                                     idims = 4   ;   idimid(1:idims) = (/1,2,3,4/)
         ENDIF
         IF( PRESENT(ktype) ) THEN   ! variable external type
            SELECT CASE (ktype)
            CASE (jp_r8)   ;   itype = NF90_DOUBLE
            CASE (jp_r4)   ;   itype = NF90_FLOAT
            CASE (jp_i4)   ;   itype = NF90_INT
            CASE (jp_i2)   ;   itype = NF90_SHORT
            CASE (jp_i1)   ;   itype = NF90_BYTE
            CASE DEFAULT   ;   CALL ctl_stop( TRIM(clinfo)//' unknown variable type' )
            END SELECT
         ELSE
            itype = NF90_DOUBLE
         ENDIF
         IF( llis0d ) THEN
            CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(cdvar), itype,                    &
               &                              iom_file(kiomid)%nvid(idvar) ), clinfo )
         ELSE
            CALL iom_nf90_check(NF90_DEF_VAR( if90id, TRIM(cdvar), itype, idimid(1:idims),   &
               &                              iom_file(kiomid)%nvid(idvar) ), clinfo )
         ENDIF
         lchunk = .false.
         IF( snc4set%luse .AND. idims == 4 )   lchunk = .true.
         ! update informations structure related the new variable we want to add...
         iom_file(kiomid)%nvars         = idvar
         iom_file(kiomid)%cn_var(idvar) = TRIM(cdvar)
         iom_file(kiomid)%scf(idvar)    = 1.
         iom_file(kiomid)%ofs(idvar)    = 0.
         iom_file(kiomid)%ndims(idvar)  = idims
         iom_file(kiomid)%luld(idvar)   = .NOT. llis0d
         DO jd = 1, idims
            CALL iom_nf90_check(NF90_INQUIRE_DIMENSION( if90id, idimid(jd), len = iom_file(kiomid)%dimsz(jd,idvar) ), clinfo)
            IF ( lchunk ) ichunksz(jd) = iom_file(kiomid)%dimsz(jd,idvar)
         END DO
         IF ( lchunk ) THEN
            ! Calculate chunk sizes by partitioning each dimension as requested in namnc4 namelist
            ! Disallow very small chunk sizes and prevent chunk sizes larger than each individual dimension
            ichunksz(1) = MIN( ichunksz(1),MAX( (ichunksz(1)-1)/snc4set%ni + 1 ,16 ) ) ! Suggested default nc4set%ni=4
            ichunksz(2) = MIN( ichunksz(2),MAX( (ichunksz(2)-1)/snc4set%nj + 1 ,16 ) ) ! Suggested default nc4set%nj=2
            ichunksz(3) = MIN( ichunksz(3),MAX( (ichunksz(3)-1)/snc4set%nk + 1 , 1 ) ) ! Suggested default nc4set%nk=6
            ichunksz(4) = 1                                                            ! Do not allow chunks to span the
            !                                                                          ! unlimited dimension
            CALL iom_nf90_check(SET_NF90_DEF_VAR_CHUNKING(if90id, idvar, ichunkalg, ichunksz), clinfo)
            CALL iom_nf90_check(SET_NF90_DEF_VAR_DEFLATE(if90id, idvar, ishuffle, ideflate, ideflate_level), clinfo)
            IF(lwp) WRITE(numout,*) TRIM(clinfo)//' chunked ok. Chunks sizes: ', ichunksz
         ENDIF
         IF(lwp) WRITE(numout,*) TRIM(clinfo)//' defined ok'
      ELSE
         idvar = kvid
      ENDIF
      !
      ! time step kwrite : write the variable
      IF( kt == kwrite ) THEN
         ! are we in write mode?
         IF( iom_file(kiomid)%irec == -1 ) THEN   ! trick: irec used to know if the file is in define mode or not
            CALL iom_nf90_check(NF90_ENDDEF( if90id ), clinfo)   ;   iom_file(kiomid)%irec = 0
         ENDIF
         ! on what kind of domain must the data be written?
         IF( llis2d .OR. llis3d ) THEN
            idimsz(1:2) = iom_file(kiomid)%dimsz(1:2,idvar)
            IF(     idimsz(1) == Ni_0 .AND. idimsz(2) == Nj_0 ) THEN
               ix1 = Nis0   ;   ix2 = Nie0   ;   iy1 = Njs0   ;   iy2 = Nje0
            ELSEIF( idimsz(1) == jpi  .AND. idimsz(2) == jpj  ) THEN
               ix1 = 1      ;   ix2 = jpi    ;   iy1 = 1      ;   iy2 = jpj
            ELSEIF( idimsz(1) == jpi  .AND. idimsz(2) == jpj  ) THEN
               ix1 = 1      ;   ix2 = jpi    ;   iy1 = 1      ;   iy2 = jpj
            ELSE
               CALL ctl_stop( 'iom_nf90_rp0123d: should have been an impossible case...' )
            ENDIF

            ! write dimension variables if it is not already done
            ! =============
            ! trick: is defined to 0 => dimension variable are defined but not yet written
            IF( iom_file(kiomid)%dimsz(1, 4) == 0 ) THEN   ! time_counter = 0
               CALL iom_nf90_check(    NF90_PUT_VAR( if90id, 1,                            glamt(ix1:ix2, iy1:iy2) ), clinfo )
               CALL iom_nf90_check(    NF90_PUT_VAR( if90id, 2,                            gphit(ix1:ix2, iy1:iy2) ), clinfo )
               SELECT CASE (iom_file(kiomid)%comp)
               CASE ('OCE')
                  CALL iom_nf90_check( NF90_PUT_VAR( if90id, 3,                                           gdept_1d ), clinfo )
               CASE ('ABL')
                  CALL iom_nf90_check( NF90_PUT_VAR( if90id, 3,                                            ght_abl ), clinfo )
               CASE DEFAULT
                  CALL iom_nf90_check( NF90_PUT_VAR( if90id, 3, (/ (idlv, idlv = 1,iom_file(kiomid)%dimsz(1,3)) /) ), clinfo )
               END SELECT
               ! "wrong" value: to be improved but not really useful...
               CALL iom_nf90_check(   NF90_PUT_VAR( if90id, 4,                                                  kt ), clinfo )
               ! update the size of the variable corresponding to the unlimited dimension
               iom_file(kiomid)%dimsz(1, 4) = 1   ! so we don't enter this IF case any more...
               IF(lwp) WRITE(numout,*) TRIM(clinfo)//' write dimension variables done'
            ENDIF
         ENDIF

         ! write the data
         ! =============
         IF(     PRESENT(pvsp0d) ) THEN
            CALL iom_nf90_check( NF90_PUT_VAR( if90id, idvar, pvsp0d                    ), clinfo )
         ELSEIF( PRESENT(pvdp0d) ) THEN
            CALL iom_nf90_check( NF90_PUT_VAR( if90id, idvar, pvdp0d                    ), clinfo )
         ELSEIF( PRESENT(pvsp1d) ) THEN
            CALL iom_nf90_check( NF90_PUT_VAR( if90id, idvar, pvsp1d(:)                 ), clinfo )
         ELSEIF( PRESENT(pvdp1d) ) THEN
            CALL iom_nf90_check( NF90_PUT_VAR( if90id, idvar, pvdp1d(:)                 ), clinfo )
         ELSEIF( PRESENT(pvsp2d) ) THEN
            CALL iom_nf90_check( NF90_PUT_VAR( if90id, idvar, pvsp2d(ix1:ix2,iy1:iy2)   ), clinfo )
         ELSEIF( PRESENT(pvdp2d) ) THEN
            CALL iom_nf90_check( NF90_PUT_VAR( if90id, idvar, pvdp2d(ix1:ix2,iy1:iy2)   ), clinfo )
         ELSEIF( PRESENT(pvsp3d) ) THEN
            CALL iom_nf90_check( NF90_PUT_VAR( if90id, idvar, pvsp3d(ix1:ix2,iy1:iy2,:) ), clinfo )
         ELSEIF( PRESENT(pvdp3d) ) THEN
            CALL iom_nf90_check( NF90_PUT_VAR( if90id, idvar, pvdp3d(ix1:ix2,iy1:iy2,:) ), clinfo )
         ENDIF
         ! add 1 to the size of the temporal dimension (not really useful...)
         IF( iom_file(kiomid)%luld(idvar) )   iom_file(kiomid)%dimsz(iom_file(kiomid)%ndims(idvar), idvar)    &
               &                            = iom_file(kiomid)%dimsz(iom_file(kiomid)%ndims(idvar), idvar) + 1
         IF(lwp) WRITE(numout,*) TRIM(clinfo)//' written ok'
      ENDIF
      !
   END SUBROUTINE iom_nf90_rstput


   SUBROUTINE iom_nf90_check( kstatus, cdinfo )
      !!--------------------------------------------------------------------
      !!                   ***  SUBROUTINE iom_nf90_check  ***
      !!
      !! ** Purpose :   check nf90 errors
      !!--------------------------------------------------------------------
      INTEGER,          INTENT(in) :: kstatus
      CHARACTER(LEN=*), INTENT(in) :: cdinfo
      !---------------------------------------------------------------------
      IF(kstatus /= nf90_noerr)   CALL ctl_stop( 'iom_nf90_check : '//TRIM(nf90_strerror(kstatus)), TRIM(cdinfo) )
   END SUBROUTINE iom_nf90_check

   !!======================================================================
END MODULE iom_nf90
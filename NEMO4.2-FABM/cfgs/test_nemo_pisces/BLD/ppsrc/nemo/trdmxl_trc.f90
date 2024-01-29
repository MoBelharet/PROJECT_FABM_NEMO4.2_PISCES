










MODULE trdmxl_trc
   !!======================================================================
   !!                       ***  MODULE  trdmxl_trc  ***
   !! Ocean diagnostics:  mixed layer passive tracer trends 
   !!======================================================================
   !! History :  9.0  !  06-08  (C. Deltel)  Original code (from trdmxl.F90)
   !!                 !  07-04  (C. Deltel)  Bug fix : add trcrad trends
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trd_mxl_trc( kt, Kmm )                                   ! Empty routine
      INTEGER, INTENT( in) ::   kt
      INTEGER, INTENT( in) ::   Kmm            ! time level index
      WRITE(*,*) 'trd_mxl_trc: You should not have seen this print! error?', kt
   END SUBROUTINE trd_mxl_trc
   SUBROUTINE trd_mxl_trc_zint( ptrc_trdmxl, ktrd, ctype, kjn, Kmm )
      INTEGER               , INTENT( in ) ::  ktrd, kjn              ! ocean trend index and passive tracer rank
      INTEGER               , INTENT( in ) ::  Kmm                    ! time level index
      CHARACTER(len=2)      , INTENT( in ) ::  ctype                  ! surface/bottom (2D) or interior (3D) physics
      REAL, DIMENSION(:,:,:), INTENT( in ) ::  ptrc_trdmxl            ! passive trc trend
      WRITE(*,*) 'trd_mxl_trc_zint: You should not have seen this print! error?', ptrc_trdmxl(1,1,1)
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ctype
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', ktrd
      WRITE(*,*) '  "      "      : You should not have seen this print! error?', kjn
   END SUBROUTINE trd_mxl_trc_zint
   SUBROUTINE trd_mxl_trc_init                                    ! Empty routine
      WRITE(*,*) 'trd_mxl_trc_init: You should not have seen this print! error?'
   END SUBROUTINE trd_mxl_trc_init

   !!======================================================================
END MODULE trdmxl_trc

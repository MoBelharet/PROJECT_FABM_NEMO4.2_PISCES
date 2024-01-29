










MODULE exampl
   !!======================================================================
   !!                       ***  MODULE  exampl  ***
   !! Ocean physics:  brief description of the purpose of the module
   !!                 (please no more than 2 lines)
   !!======================================================================
   !! History : 3.0  !  2008-06  (Author Names)  Original code
   !!            -   !  2008-08  (Author names)  brief description of modifications
   !!           3.3  !  2010-11  (Author names)        -              -
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                                         NO example
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE exa_mpl( kt, pvar1, pvar2, ptab )              ! Empty routine
      INTEGER :: kt
      REAL::   pvar1, pvar2, ptab(:,:)
      WRITE(*,*) 'exa_mpl: You should not have seen this print! error?', kt, pvar1, pvar2, ptab(1,1)
   END SUBROUTINE exa_mpl

   !!======================================================================
END MODULE exampl

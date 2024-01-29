











    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    module ppr_1d

    !    
    ! PPR-1D.F90: 1-d piecewise polynomial reconstructions.
    !
    ! Darren Engwirda 
    ! 31-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    implicit none

    !------------------------------------ method selection !
        
    integer, parameter :: p1e_method = +100
    integer, parameter :: p3e_method = +101
    integer, parameter :: p5e_method = +102

    integer, parameter :: pcm_method = +200
    integer, parameter :: plm_method = +201
    integer, parameter :: ppm_method = +202
    integer, parameter :: pqm_method = +203

    integer, parameter :: null_limit = +300
    integer, parameter :: mono_limit = +301
    integer, parameter :: weno_limit = +302

    integer, parameter :: bcon_loose = +400
    integer, parameter :: bcon_value = +401
    integer, parameter :: bcon_slope = +402
 
    type rmap_tics
    !------------------------------- tCPU timer for RCON1D !
        integer :: rmap_time
        integer :: edge_time
        integer :: cell_time
        integer :: oscl_time
    end type rmap_tics

    type rcon_opts
    !------------------------------- parameters for RCON1D !
        integer :: edge_meth
        integer :: cell_meth
        integer :: cell_lims
        integer :: wall_lims
    end type rcon_opts

    type rcon_ends
    !------------------------------- end-conditions struct !
        integer :: bcopt
        real*8  :: value
        real*8  :: slope
    end type rcon_ends

    type rcon_work
    !------------------------------- work-space for RCON1D !
        real*8, allocatable  :: edge_func(:,:)
        real*8, allocatable  :: edge_dfdx(:,:)
        real*8, allocatable  :: cell_oscl(:,:,:)
    contains
        procedure :: init => init_rcon_work
        procedure :: free => free_rcon_work
    end type rcon_work

    type, extends(rcon_opts) :: rmap_opts
    !------------------------------- parameters for RMAP1D !
    end type rmap_opts

    type, extends(rcon_work) :: rmap_work
    !------------------------------- work-space for RMAP1D !
        real*8, allocatable  :: cell_spac(:)
        real*8, allocatable  :: cell_func(:,:,:)
    contains
        procedure :: init => init_rmap_work
        procedure :: free => free_rmap_work
    end type rmap_work

    contains
 
    !------------------------------------------------------!
    ! INIT-RCON-WORK: init. work-space for RCON1D.         !
    !------------------------------------------------------!
    
    subroutine init_rcon_work(this,npos,nvar,opts)

    !
    ! THIS  work-space structure for RCON1D .
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! OPTS  parameters structure for RCON1D .
    !

        implicit none

    !------------------------------------------- arguments !
        class(rcon_work) , intent(inout) :: this
        integer, intent(in):: npos
        integer, intent(in):: nvar
        class(rcon_opts) , optional      :: opts

    !------------------------------------------- variables !
        integer :: okay

        allocate(this% &
        &   edge_func(  nvar,npos), &
        &        this% &
        &   edge_dfdx(  nvar,npos), &
        &        this% &
        &   cell_oscl(2,nvar,npos), &
        &   stat=okay)
        
    end subroutine

    !------------------------------------------------------!
    ! INIT-RMAP-WORK: init. work-space for RMAP1D.         !
    !------------------------------------------------------!
    
    subroutine init_rmap_work(this,npos,nvar,opts)

    !
    ! THIS  work-space structure for RMAP1D .
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! OPTS  parameters structure for RMAP1D .
    !

        implicit none

    !------------------------------------------- arguments !
        class(rmap_work) , intent(inout) :: this
        integer, intent(in) :: npos
        integer, intent(in) :: nvar
        class(rcon_opts) , optional      :: opts

    !------------------------------------------- variables !
        integer :: okay,ndof

        ndof = ndof1d(opts%cell_meth)

        allocate(this% &
        &   edge_func(  nvar,npos), &
        &        this% &
        &   edge_dfdx(  nvar,npos), &
        &        this% &
        &   cell_oscl(2,nvar,npos), &
        &        this% &
        &   cell_spac(       npos), &
        &        this% &
        &   cell_func(ndof,nvar,npos) , &
        &   stat=okay)
 
    end subroutine

    !------------------------------------------------------!
    ! FREE-RCON-WORK: free work-space for RCON1D .         !
    !------------------------------------------------------!
    
    subroutine free_rcon_work(this)

        implicit none

    !------------------------------------------- arguments !
        class(rcon_work), intent(inout) :: this

        deallocate(this%edge_func, &
        &          this%edge_dfdx, &
        &          this%cell_oscl)
     
    end subroutine

    !------------------------------------------------------!
    ! FREE-RMAP-WORK: free work-space for RMAP1D .         !
    !------------------------------------------------------!
    
    subroutine free_rmap_work(this)

        implicit none

    !------------------------------------------- arguments !
        class(rmap_work), intent(inout) :: this


        deallocate(this%edge_func, &
        &          this%edge_dfdx, &
        &          this%cell_oscl, &
        &          this%cell_func, &
        &          this%cell_spac)

    end subroutine

    !------------------------------------------------------!
    ! NDOF1D : no. degrees-of-freedom per polynomial .     !
    !------------------------------------------------------!
    
    pure function ndof1d(meth) result(rdof)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: meth

    !------------------------------------------- variables !
        integer  :: rdof

        select case(meth)
    !-------------------------------- edge reconstructions !
        case (p1e_method)
            rdof = +2
        case (p3e_method)
            rdof = +4
        case (p5e_method)
            rdof = +6
    !-------------------------------- cell reconstructions !
        case (pcm_method)
            rdof = +1
        case (plm_method)
            rdof = +2
        case (ppm_method)
            rdof = +3
        case (pqm_method)
            rdof = +5
        
        case default 
            rdof = +0

        end select
        
    end function  ndof1d

    !------------------------------------------------------!
    ! BFUN1D : one-dimensional poly. basis-functions .     !
    !------------------------------------------------------!
        

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! BFUN1D.h90: poly. basis-functions for reconstruction.
    !
    ! Darren Engwirda 
    ! 07-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine bfun1d(isel,ndof,sval,bfun)

    !
    ! ISEL  basis-function "order", -1 => integral-basis , 
    !       +0 => function-basis, +1 => 1st deriv.-basis ,
    !       +2 => 2nd deriv.-basis.
    ! NDOF  no. degrees-of-freedom in basis.
    ! SVAL  local coord. at which to evaluate basis-func.,
    !       such that -1.0 <= SVAL <= +1.0 .
    ! BFUN  basis-vector evaluated at SVAL .
    !
    
        implicit none
        
    !------------------------------------------- arguments !
        integer, intent( in) :: isel,ndof
        real*8 , intent( in) :: sval
        real*8 , intent(out) :: bfun(:)
        
        select case (isel)
        case (-1)
    !------------------------------------ -1th-order basis !
            select case (ndof)
            case (+1)
                bfun(1) = sval**1 / 1.e0_8
                
            case (+2)
                bfun(1) = sval**1 / 1.e0_8
                bfun(2) = sval**2 / 2.e0_8
                
            case (+3)
                bfun(1) = sval**1 / 1.e0_8
                bfun(2) = sval**2 / 2.e0_8
                bfun(3) = sval**3 / 3.e0_8
                
            case (+4)
                bfun(1) = sval**1 / 1.e0_8
                bfun(2) = sval**2 / 2.e0_8
                bfun(3) = sval**3 / 3.e0_8
                bfun(4) = sval**4 / 4.e0_8
                
            case (+5)
                bfun(1) = sval**1 / 1.e0_8
                bfun(2) = sval**2 / 2.e0_8
                bfun(3) = sval**3 / 3.e0_8
                bfun(4) = sval**4 / 4.e0_8
                bfun(5) = sval**5 / 5.e0_8
                
            case (+6)
                bfun(1) = sval**1 / 1.e0_8
                bfun(2) = sval**2 / 2.e0_8
                bfun(3) = sval**3 / 3.e0_8
                bfun(4) = sval**4 / 4.e0_8
                bfun(5) = sval**5 / 5.e0_8
                bfun(6) = sval**6 / 6.e0_8
                
            case (+7)
                bfun(1) = sval**1 / 1.e0_8
                bfun(2) = sval**2 / 2.e0_8
                bfun(3) = sval**3 / 3.e0_8
                bfun(4) = sval**4 / 4.e0_8
                bfun(5) = sval**5 / 5.e0_8
                bfun(6) = sval**6 / 6.e0_8
                bfun(7) = sval**7 / 7.e0_8
            
            end select

        case (+0)
    !------------------------------------ +0th-order basis !
            select case (ndof)
            case (+1)
                bfun(1) =           1.e0_8
                
            case (+2)
                bfun(1) =           1.e0_8
                bfun(2) = sval**1 * 1.e0_8
                
            case (+3)
                bfun(1) =           1.e0_8
                bfun(2) = sval**1 * 1.e0_8
                bfun(3) = sval**2 * 1.e0_8
                
            case (+4)
                bfun(1) =           1.e0_8
                bfun(2) = sval**1 * 1.e0_8
                bfun(3) = sval**2 * 1.e0_8
                bfun(4) = sval**3 * 1.e0_8
                
            case (+5)
                bfun(1) =           1.e0_8
                bfun(2) = sval**1 * 1.e0_8
                bfun(3) = sval**2 * 1.e0_8
                bfun(4) = sval**3 * 1.e0_8
                bfun(5) = sval**4 * 1.e0_8
                
            case (+6)
                bfun(1) =           1.e0_8
                bfun(2) = sval**1 * 1.e0_8
                bfun(3) = sval**2 * 1.e0_8
                bfun(4) = sval**3 * 1.e0_8
                bfun(5) = sval**4 * 1.e0_8
                bfun(6) = sval**5 * 1.e0_8
                
            case (+7)
                bfun(1) =           1.e0_8
                bfun(2) = sval**1 * 1.e0_8
                bfun(3) = sval**2 * 1.e0_8
                bfun(4) = sval**3 * 1.e0_8
                bfun(5) = sval**4 * 1.e0_8
                bfun(6) = sval**5 * 1.e0_8
                bfun(7) = sval**6 * 1.e0_8
            
            end select

        case (+1)
    !------------------------------------ +1st-order basis !
            select case (ndof)
            case (+1)
                bfun(1) =           0.e0_8
                
            case (+2)
                bfun(1) =           0.e0_8
                bfun(2) =           1.e0_8
                
            case (+3)
                bfun(1) =           0.e0_8
                bfun(2) =           1.e0_8
                bfun(3) = sval**1 * 2.e0_8
                
            case (+4)
                bfun(1) =           0.e0_8
                bfun(2) =           1.e0_8
                bfun(3) = sval**1 * 2.e0_8
                bfun(4) = sval**2 * 3.e0_8
                
            case (+5)
                bfun(1) =           0.e0_8
                bfun(2) =           1.e0_8
                bfun(3) = sval**1 * 2.e0_8
                bfun(4) = sval**2 * 3.e0_8
                bfun(5) = sval**3 * 4.e0_8
                
            case (+6)
                bfun(1) =           0.e0_8
                bfun(2) =           1.e0_8
                bfun(3) = sval**1 * 2.e0_8
                bfun(4) = sval**2 * 3.e0_8
                bfun(5) = sval**3 * 4.e0_8
                bfun(6) = sval**4 * 5.e0_8
                
            case (+7)
                bfun(1) =           0.e0_8
                bfun(2) =           1.e0_8
                bfun(3) = sval**1 * 2.e0_8
                bfun(4) = sval**2 * 3.e0_8
                bfun(5) = sval**3 * 4.e0_8
                bfun(6) = sval**4 * 5.e0_8
                bfun(7) = sval**5 * 6.e0_8
            
            end select

        case (+2)
    !------------------------------------ +2nd-order basis !
            select case (ndof)
            case (+1)
                bfun(1) =           0.e0_8
                
            case (+2)
                bfun(1) =           0.e0_8
                bfun(2) =           0.e0_8
                
            case (+3)
                bfun(1) =           0.e0_8
                bfun(2) =           0.e0_8
                bfun(3) =           2.e0_8
                
            case (+4)
                bfun(1) =           0.e0_8
                bfun(2) =           0.e0_8
                bfun(3) =           2.e0_8
                bfun(4) = sval**1 * 6.e0_8
                
            case (+5)
                bfun(1) =           0.e0_8
                bfun(2) =           0.e0_8
                bfun(3) =           2.e0_8
                bfun(4) = sval**1 * 6.e0_8
                bfun(5) = sval**2 *12.e0_8
                
            case (+6)
                bfun(1) =           0.e0_8
                bfun(2) =           0.e0_8
                bfun(3) =           2.e0_8
                bfun(4) = sval**1 * 6.e0_8
                bfun(5) = sval**2 *12.e0_8
                bfun(6) = sval**3 *20.e0_8
                
            case (+7)
                bfun(1) =           0.e0_8
                bfun(2) =           0.e0_8
                bfun(3) =           2.e0_8
                bfun(4) = sval**1 * 6.e0_8
                bfun(5) = sval**2 *12.e0_8
                bfun(6) = sval**3 *20.e0_8
                bfun(7) = sval**4 *30.e0_8
            
            end select

        end select
    
    end subroutine
    
    
    

    !------------------------------------------------------!
    ! UTIL1D : one-dimensional grid manip. utilities .     !
    !------------------------------------------------------!
    

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! UTIL1D.h90: util. func. for 1-dim. grid manipulation.
    !
    ! Darren Engwirda 
    ! 31-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine linspace(xxll,xxuu,npos,xpos)

    ! 
    ! XXLL  lower-bound grid position.
    ! NNEW  upper-bound grid position.
    ! NPOS  no. edges in the grid.
    ! XPOS  array of grid edges. XPOS has length NPOS .
    !

        implicit none

        real*8 , intent(in)  :: xxll,xxuu
        integer, intent(in)  :: npos
        real*8 , intent(out) :: xpos(:)

        integer   :: ipos
        real*8    :: xdel

        xpos(   1) = xxll
        xpos(npos) = xxuu

        xdel = (xxuu-xxll) / (npos - 1)

        do  ipos = +2, npos-1

            xpos(ipos) = (ipos-1) * xdel         

        end do

        return
    
    end  subroutine

    subroutine rndspace(xxll,xxuu,npos,xpos, &
        &               frac)

    ! 
    ! XXLL  lower-bound grid position.
    ! NNEW  upper-bound grid position.
    ! NPOS  no. edges in the grid.
    ! XPOS  array of grid edges. XPOS has length NPOS .
    ! FRAC  fractional perturbation of cell, OPTIONAL .
    !

        implicit none

        real*8 , intent(in)  :: xxll,xxuu
        integer, intent(in)  :: npos
        real*8 , intent(out) :: xpos(:)
        real*8 , intent(in), optional :: frac

        integer   :: ipos
        real*8    :: xdel,rand,move

        if (present(frac)) then
            move = +frac
        else
            move = 0.33d0
        end if

        xpos(   1) = xxll
        xpos(npos) = xxuu

        xdel = (xxuu-xxll) / (npos - 1)

        do  ipos = +2, npos-1

            xpos(ipos) = (ipos-1) * xdel         

        end do

        do ipos = +2, npos-1

            call random_number (rand)

            rand = 2.e0_8 * (rand-.5d0)
            
            move = rand *  move

            xpos(ipos) = &
        &       xpos(ipos) + move * xdel         

        end do

        return

    end  subroutine




    !------------------------------------------------------!
    ! WENO1D : "essentially" non-oscillatory limiter .     !
    !------------------------------------------------------!


    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! WENO1D.h90: WENO-style slope-limiting for 1d reconst.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine wenoi (npos,delx,oscl,ipos, &
        &                  ivar,halo,&
        &                  wlim,wval )

    !
    ! NPOS  no. edges over grid.
    ! DELX  grid-cell spacing array. SIZE(DELX) == +1 if
    !       the grid is uniformly spaced .
    ! OSCL  cell-centred oscillation-detectors, where OSCL 
    !       has SIZE = +2-by-NVAR-by-NPOS-1. OSCL is given
    !       by calls to OSCLI().
    ! IPOS  grid-cell index for which to calc. weights .
    ! IVAR  state-var index for which to calc/ weights .
    ! HALO  width of recon. stencil, symmetric about IPOS .
    ! WLIM  limiter treatment at endpoints, monotonic or
    !       otherwise . 
    ! WVAL  WENO weights vector, such that FHAT = WVAL(1) *
    !       UHAT + WVAL(2) * LHAT, where UHAT and LHAT are 
    !       the unlimited and monotonic grid-cell profiles
    !       respectively .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        integer, intent(in)  :: wlim
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: wval(2)
    
    !------------------------------------------- variables !
        real*8 :: omin,omax,wsum
        
        real*8 , parameter :: ZERO = +1.e-16
       
        if (size(delx).gt.+1) then
        
    !------------------- use variable grid spacing variant !
        
        call wenov(npos,delx,oscl, &
        &          ipos,ivar,halo, &
        &          wlim,omin,omax)
        
        else
        
    !------------------- use constant grid spacing variant !
        
        call wenoc(npos,delx,oscl, &
        &          ipos,ivar,halo, &
        &          wlim,omin,omax)
        
        end if

    !------------------ compute WENO-style profile weights !

        omax = omax + ZERO
        omin = omin + ZERO

        if (halo .ge. +3) then

        wval(1) = +1.0d+7 / omax ** 3
        wval(2) = +1.0d+0 / omin ** 3

        else &
    &   if (halo .le. +2) then
    
        wval(1) = +1.0d+5 / omax ** 3
        wval(2) = +1.0d+0 / omin ** 3
    
        end if

        wsum = wval(1) + wval(2) + ZERO
        wval(1) = wval(1) / wsum
    !   wval(2) = wval(2) / wsum
        wval(2) =-wval(1) + 1.e0_8 ! wval(2)/wsum but robust !

        return

    end  subroutine

    pure subroutine wenov (npos,delx,oscl,ipos, &
        &                  ivar,halo,&
        &                  wlim,omin,omax)

    !
    ! *this is the variable grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! DELX  grid-cell spacing array. SIZE(DELX) == +1 if
    !       the grid is uniformly spaced .
    ! OSCL  cell-centred oscillation-detectors, where OSCL 
    !       has SIZE = +2-by-NVAR-by-NPOS-1. OSCL is given
    !       by calls to OSCLI().
    ! IPOS  grid-cell index for which to calc. weights .
    ! IVAR  state-var index for which to calc/ weights .
    ! HALO  width of recon. stencil, symmetric about IPOS .
    ! WLIM  limiter treatment at endpoints, monotonic or
    !       otherwise . 
    ! OMIN  min. and max. oscillation indicators over the
    ! OMAX  local re-con. stencil .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        integer, intent(in)  :: wlim
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: omin,omax

    !------------------------------------------- variables !
        integer :: hpos
        integer :: head,tail
        integer :: imin,imax
        real*8  :: deli,delh
        real*8  :: hh00,hsqr
        real*8  :: dfx1,dfx2
        real*8  :: oval
        
    !------------------- calc. lower//upper stencil bounds !    

        head = 1; tail = npos - 1

        if(wlim.eq.mono_limit) then

    !---------------------- deactivate WENO at boundaries !
        
        if (ipos-halo.lt.head) then
        
            omax = 1.e0_8 
            omin = 0.e0_8 ; return
        
        end if
        
        if (ipos+halo.gt.tail) then
            
            omax = 1.e0_8 
            omin = 0.e0_8 ; return
        
        end if

        end if

    !---------------------- truncate stencil at boundaries !
    
       imin = max(ipos-halo,head)
       imax = min(ipos+halo,tail)
      
    !------------------ find min/max indicators on stencil !
 
        dfx1 = oscl(1,ivar,ipos)
        dfx2 = oscl(2,ivar,ipos)
      
        hh00 = delx(ipos+0)**1
        hsqr = delx(ipos+0)**2
      
        oval =(hh00 * dfx1)**2 &
        &    +(hsqr * dfx2)**2
      
        omin = oval
        omax = oval
      
    !---------------------------------------- "lower" part !
      
        delh = 0.e0_8

        do  hpos = ipos-1, imin, -1
        
    !------------------ calc. derivatives centred on IPOS. !     
   
        deli = delx(hpos+0) &
    &        + delx(hpos+1)
    
        delh = delh + deli*.5d0
        
        dfx1 = oscl(1,ivar,hpos)
        dfx2 = oscl(2,ivar,hpos)

        dfx1 = dfx1 + dfx2*delh

    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

        oval = (hh00 * dfx1)**2 &
    &        + (hsqr * dfx2)**2

        if (oval .lt. omin) then
            omin = oval 
        else &
    &   if (oval .gt. omax) then
            omax = oval
        end if
        
        end do
      
    !---------------------------------------- "upper" part !     
      
        delh = 0.e0_8
        
        do  hpos = ipos+1, imax, +1
        
    !------------------ calc. derivatives centred on IPOS. !     
   
        deli = delx(hpos+0) &
    &        + delx(hpos-1)
    
        delh = delh - deli*.5d0
        
        dfx1 = oscl(1,ivar,hpos)
        dfx2 = oscl(2,ivar,hpos)

        dfx1 = dfx1 + dfx2*delh

    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

        oval = (hh00 * dfx1)**2 &
    &        + (hsqr * dfx2)**2

        if (oval .lt. omin) then
            omin = oval 
        else &
    &   if (oval .gt. omax) then
            omax = oval
        end if
        
        end do

        return

    end  subroutine
    
    pure subroutine wenoc (npos,delx,oscl,ipos, &
        &                  ivar,halo,&
        &                  wlim,omin,omax)

    !
    ! *this is the constant grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! DELX  grid-cell spacing array. SIZE(DELX) == +1 if
    !       the grid is uniformly spaced .
    ! OSCL  cell-centred oscillation-detectors, where OSCL 
    !       has SIZE = +2-by-NVAR-by-NPOS-1. OSCL is given
    !       by calls to OSCLI().
    ! IPOS  grid-cell index for which to calc. weights .
    ! IVAR  state-var index for which to calc/ weights .
    ! HALO  width of recon. stencil, symmetric about IPOS .
    ! WLIM  limiter treatment at endpoints, monotonic or
    !       otherwise . 
    ! OMIN  min. and max. oscillation indicators over the
    ! OMAX  local re-con. stencil .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,halo
        integer, intent(in)  :: ipos,ivar
        integer, intent(in)  :: wlim
        real*8 , intent(in)  :: delx(1)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(out) :: omin,omax

    !------------------------------------------- variables !
        integer :: hpos
        integer :: head,tail
        integer :: imin,imax
        real*8  :: delh
        real*8  :: dfx1,dfx2
        real*8  :: oval
    
    !------------------- calc. lower//upper stencil bounds !    

        head = 1; tail = npos - 1

        if(wlim.eq.mono_limit) then

    !---------------------- deactivate WENO at boundaries !
        
        if (ipos-halo.lt.head) then
        
            omax = 1.e0_8 
            omin = 0.e0_8 ; return
        
        end if
        
        if (ipos+halo.gt.tail) then
            
            omax = 1.e0_8 
            omin = 0.e0_8 ; return
        
        end if

        end if

    !---------------------- truncate stencil at boundaries !
        
        imin = max(ipos-halo,head)
        imax = min(ipos+halo,tail)
      
    !------------------ find min/max indicators on stencil !
 
        dfx1 = oscl(1,ivar,ipos)
        dfx2 = oscl(2,ivar,ipos)
      
        oval = (2.e0_8**1*dfx1)**2 &
        &    + (2.e0_8**2*dfx2)**2
      
        omin = oval
        omax = oval
      
    !---------------------------------------- "lower" part !
      
        delh = 0.e0_8

        do  hpos = ipos-1, imin, -1
        
    !------------------ calc. derivatives centred on IPOS. !     

        delh = delh + 2.e0_8
        
        dfx1 = oscl(1,ivar,hpos)
        dfx2 = oscl(2,ivar,hpos)

        dfx1 = dfx1 + dfx2*delh

    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

        oval = (2.e0_8**1*dfx1)**2 &
    &        + (2.e0_8**2*dfx2)**2

        if (oval .lt. omin) then
            omin = oval 
        else &
    &   if (oval .gt. omax) then
            omax = oval
        end if
        
        end do
      
    !---------------------------------------- "upper" part !     
      
        delh = 0.e0_8
        
        do  hpos = ipos+1, imax, +1
        
    !------------------ calc. derivatives centred on IPOS. !     

        delh = delh - 2.e0_8
        
        dfx1 = oscl(1,ivar,hpos)
        dfx2 = oscl(2,ivar,hpos)

        dfx1 = dfx1 + dfx2*delh
    
    !------------------ indicator: NORM(H^N * D^N/DX^N(F)) !

        oval = (2.e0_8**1*dfx1)**2 &
    &        + (2.e0_8**2*dfx2)**2
  
        if (oval .lt. omin) then
            omin = oval 
        else &
    &   if (oval .gt. omax) then
            omax = oval
        end if

        end do

        return

    end  subroutine
    
    
    

    
    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! OSCL1D.h90: "oscillation-indicators" for WENO interp.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !
    
    pure subroutine oscli (npos,nvar,ndof,delx,&
        &                  fdat,oscl,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  (constant) grid-cell spacing. LENGTH(DELX)==+1 .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: ivar,ipos
        
        if (npos.lt.3) then
    !------------------------------- at least 3 grid-cells !
        do  ipos = +1, npos-1
        do  ivar = +1, nvar-0
            oscl(1,ivar,ipos) = +0.e0_8
            oscl(2,ivar,ipos) = +0.e0_8
        end do
        end do
        end if
        
        if (npos.lt.3) return
        if (nvar.lt.1) return
        if (ndof.lt.1) return

        if (size(delx).gt.+1) then

    !------------------------------- variable grid-spacing !

            call osclv(npos,nvar,ndof,delx, &
        &              fdat,oscl,dmin)
        
        else

    !------------------------------- constant grid-spacing !
        
            call osclc(npos,nvar,ndof,delx, &
        &              fdat,oscl,dmin)
                
        end if

        return
        
    end  subroutine
    
    pure subroutine osclv (npos,nvar,ndof,delx,&
        &                  fdat,oscl,dmin)

    !
    ! *this is the variable grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  (variable) grid-cell spacing. LENGTH(DELX)!=+1 .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: head,tail
        integer :: ipos,ivar
        real*8  :: hhll,hhcc,hhrr
        real*8  :: hhmm,hhrc,hhlc
        real*8  :: cmat(2,3)

        head = +1 ; tail = npos-1

    !--------------------------------------- centred point !

        do  ipos = head+1, tail-1

        hhll = max(delx(ipos-1),dmin)
        hhcc = max(delx(ipos+0),dmin)
        hhrr = max(delx(ipos+1),dmin)

        hhrc = hhrr + hhcc
        hhlc = hhll + hhcc
        hhmm = hhll + hhcc + hhrr

        cmat(1,1) = -(hhcc+2.e0_8*hhrr)/(hhlc*hhmm)
        cmat(1,2) = -(hhll-hhrr)* &
        &          (3.e0_8*hhcc+2.e0_8*(hhll+hhrr))/&
        &            (hhlc*hhrc*hhmm)
        cmat(1,3) = +(hhcc+2.e0_8*hhll)/(hhrc*hhmm)

        cmat(2,1) = +3.e0_8/(hhlc*hhmm)
        cmat(2,2) = -3.e0_8*(2.e0_8*hhcc+hhll+hhrr)/&
        &            (hhlc*hhrc*hhmm)
        cmat(2,3) = +3.e0_8/(hhrc*hhmm)

        do  ivar = 1, nvar

            oscl(1,ivar,ipos) = +1.e0_8 * ( &
        & + cmat(1,1)*fdat(1,ivar,ipos-1) &
        & + cmat(1,2)*fdat(1,ivar,ipos+0) &
        & + cmat(1,3)*fdat(1,ivar,ipos+1) )

            oscl(2,ivar,ipos) = +2.e0_8 * ( &
        & + cmat(2,1)*fdat(1,ivar,ipos-1) &
        & + cmat(2,2)*fdat(1,ivar,ipos+0) &
        & + cmat(2,3)*fdat(1,ivar,ipos+1) )

        end do
            
        end do

    !-------------------------------------- lower endpoint !

        hhll = max(delx(head+0),dmin)
        hhcc = max(delx(head+1),dmin)
        hhrr = max(delx(head+2),dmin)

        cmat(1,1) = -2.e0_8 / (hhll+hhcc)
        cmat(1,2) = +2.e0_8 / (hhll+hhcc)

        do  ivar = 1, nvar

            oscl(1,ivar,head) = &
        & + cmat(1,1)*fdat(1,ivar,head+0) &
        & + cmat(1,2)*fdat(1,ivar,head+1)
             
            oscl(2,ivar,head) = +0.e0_8

        end do

    !-------------------------------------- upper endpoint !

        hhll = max(delx(tail-2),dmin)
        hhcc = max(delx(tail-1),dmin)
        hhrr = max(delx(tail-0),dmin)

        cmat(1,2) = -2.e0_8 / (hhrr+hhcc)
        cmat(1,3) = +2.e0_8 / (hhrr+hhcc)

        do  ivar = 1, nvar

            oscl(1,ivar,tail) = &
        & + cmat(1,2)*fdat(1,ivar,tail-1) &
        & + cmat(1,3)*fdat(1,ivar,tail+0)

            oscl(2,ivar,tail) = +0.e0_8

        end do
        
        return

    end  subroutine
    
    pure subroutine osclc (npos,nvar,ndof,delx,&
        &                  fdat,oscl,dmin)

    !
    ! *this is the constant grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  (constant) grid-cell spacing. LENGTH(DELX)==+1 .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(1)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: oscl(:,:,:)

    !------------------------------------------- variables !
        integer :: head,tail,ipos,ivar
        
        head = +1; tail = npos - 1

    !-------------------------------------- centred points !

        do  ipos = 2, npos-2
        do  ivar = 1, nvar-0

            oscl(1,ivar,ipos) = &
        & + .25d+0 * fdat(1,ivar,ipos+1) &
        & - .25d+0 * fdat(1,ivar,ipos-1)
        
            oscl(2,ivar,ipos) = &
        & + .25d+0 * fdat(1,ivar,ipos+1) &
        & - .50d+0 * fdat(1,ivar,ipos+0) &
        & + .25d+0 * fdat(1,ivar,ipos-1)

        end do
        end do

    !-------------------------------------- lower endpoint !

        do  ivar = 1, nvar

            oscl(1,ivar,head) = &
        & + .50d+0 * fdat(1,ivar,head+1) &
        & - .50d+0 * fdat(1,ivar,head+0)
        
            oscl(2,ivar,head) = +0.e0_8

        end do

    !-------------------------------------- upper endpoint !

        do  ivar = 1, nvar

            oscl(1,ivar,tail) = &
        & + .50d+0 * fdat(1,ivar,tail+0) &
        & - .50d+0 * fdat(1,ivar,tail-1)
            
            oscl(2,ivar,tail) = +0.e0_8

        end do
        
        return

    end  subroutine
    
    


    !------------------------------------------------------!
    ! RCON1D : one-dimensional poly. reconstructions .     !
    !------------------------------------------------------!


    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! RCON1D.h90: conservative, polynomial reconstructions.
    !
    ! Darren Engwirda 
    ! 07-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !
    
    subroutine rcon1d(npos,nvar,ndof,delx,fdat, &
        &             bclo,bchi,fhat,work,opts, &
        &             tCPU)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 . 
    ! BCLO  boundary condition at lower endpoint.
    ! BCHI  boundary condition at upper endpoint.
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! WORK  method work-space. See RCON-WORK for details .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! TCPU  method tcpu-timer.
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        class(rcon_work), intent(inout):: work
        class(rcon_opts), intent(in)   :: opts
        real*8 , intent( in) :: delx(:)
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bchi(:)
        type (rmap_tics), &
        &   intent(inout) , optional :: tCPU

    !------------------------------------------- variables !
        integer :: halo,ipos
        real*8  :: dmin,dmid


        if (ndof.lt.1) return
        if (npos.lt.2) return
        if (nvar.lt.1) return
     
    !-------------------------- compute min grid-tolerance !
     
        dmid = delx(1)
     
        if (size(delx).gt.+1) then
        
            do  ipos = 2, npos-1
                dmid = &
            &   dmid + delx (ipos)
            end do
        
            dmid = dmid /(npos-1)
        
        end if

        dmin = +1.0d-14 * dmid
        
    !-------------------------- compute edge values/slopes !


        if ( (opts%cell_meth.eq.ppm_method) &
    &  .or.  (opts%cell_meth.eq.pqm_method) ) then

        select case (opts%edge_meth)
            case(p1e_method)
    !------------------------------------ 2nd-order method !
            halo = +1
            call p1e(npos,nvar,ndof, &
            &        delx,fdat,      &
            &        bclo,bchi,      &
            &        work%edge_func, &
            &        work%edge_dfdx, &
            &        opts,dmin)

            case(p3e_method)
    !------------------------------------ 4th-order method !           
            halo = +2
            call p3e(npos,nvar,ndof, &
            &        delx,fdat,      &
            &        bclo,bchi,      &
            &        work%edge_func, &
            &        work%edge_dfdx, &
            &        opts,dmin)

            case(p5e_method)
    !------------------------------------ 6th-order method !           
            halo = +3
            call p5e(npos,nvar,ndof, &
            &        delx,fdat,      &
            &        bclo,bchi,      &
            &        work%edge_func, &
            &        work%edge_dfdx, &
            &        opts,dmin)

        end select

        end if
        

    !-------------------------- compute oscil. derivatives !


        if (opts%cell_lims.eq.weno_limit) then
    
            call oscli(npos,nvar,ndof, &
            &          delx,fdat, &
            &          work%cell_oscl, &
            &          dmin)
     
        end if
    

    !-------------------------- compute grid-cell profiles !


        select case (opts%cell_meth)
            case(pcm_method)
    !------------------------------------ 1st-order method !
            call pcm(npos,nvar,ndof, &
            &        fdat,fhat)

            case(plm_method)
    !------------------------------------ 2nd-order method !
            call plm(npos,nvar,ndof, &
            &        delx,fdat,fhat, &
            &        dmin,&
            &        opts%cell_lims)

            case(ppm_method)
    !------------------------------------ 3rd-order method !
            call ppm(npos,nvar,ndof, &
            &        delx,fdat,fhat, &
            &        work%edge_func, &
            &        work%cell_oscl, &
            &        dmin,&
            &        opts%cell_lims, &
            &        opts%wall_lims, &
            &        halo )

            case(pqm_method)
    !------------------------------------ 5th-order method !
            call pqm(npos,nvar,ndof, &
            &        delx,fdat,fhat, &
            &        work%edge_func, &
            &        work%edge_dfdx, &
            &        work%cell_oscl, &
            &        dmin,&
            &        opts%cell_lims, &
            &        opts%wall_lims, &
            &        halo )

        end select


    end subroutine
    
    
    


    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! INV.h90: block-wise solution of small linear systems.
    !
    ! Darren Engwirda 
    ! 25-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine inv_2x2(amat,adim,ainv,vdim, &
        &                   adet)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: adim
        real*8 , intent( in) :: amat(adim,*)
        integer, intent( in) :: vdim
        real*8 , intent(out) :: ainv(vdim,*)
        real*8 , intent(out) :: adet
        
    !------------------------------------------- form A^-1 !

        adet   = amat(1,1) * amat(2,2) &
               - amat(1,2) * amat(2,1)

        ainv(1,1) =          amat(2,2)
        ainv(1,2) =        - amat(1,2)
        ainv(2,1) =        - amat(2,1)
        ainv(2,2) =          amat(1,1)

        return

    end  subroutine

    pure subroutine inv_3x3(amat,adim,ainv,vdim, &
        &                   adet)

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: adim
        real*8 , intent( in) :: amat(adim,*)
        integer, intent( in) :: vdim
        real*8 , intent(out) :: ainv(vdim,*)
        real*8 , intent(out) :: adet
        
    !------------------------------------------- variables !
        real*8 :: &
        aa2233,aa2332,aa2133,aa2331,aa2132,&
        aa2231,aa1233,aa1332,aa1223,aa1322,&
        aa1133,aa1331,aa1123,aa1321,aa1132,&
        aa1231,aa1122,aa1221

    !------------------------------------------- form A^-1 !

        aa2233 = amat(2,2) * amat(3,3)
        aa2332 = amat(2,3) * amat(3,2)
        aa2133 = amat(2,1) * amat(3,3)
        aa2331 = amat(2,3) * amat(3,1)
        aa2132 = amat(2,1) * amat(3,2)
        aa2231 = amat(2,2) * amat(3,1)

        adet =  &
        amat(1,1) *  (aa2233 - aa2332) - &
	    amat(1,2) *  (aa2133 - aa2331) + &
	    amat(1,3) *  (aa2132 - aa2231)

        aa1233 = amat(1,2) * amat(3,3)
        aa1332 = amat(1,3) * amat(3,2)
        aa1223 = amat(1,2) * amat(2,3)
        aa1322 = amat(1,3) * amat(2,2)
        aa1133 = amat(1,1) * amat(3,3)
        aa1331 = amat(1,3) * amat(3,1)
        aa1123 = amat(1,1) * amat(2,3)
        aa1321 = amat(1,3) * amat(2,1)
        aa1132 = amat(1,1) * amat(3,2)
        aa1231 = amat(1,2) * amat(3,1)
        aa1122 = amat(1,1) * amat(2,2)
        aa1221 = amat(1,2) * amat(2,1)

        ainv(1,1) =  (aa2233 - aa2332)
        ainv(1,2) = -(aa1233 - aa1332)
        ainv(1,3) =  (aa1223 - aa1322)
        
        ainv(2,1) = -(aa2133 - aa2331)
        ainv(2,2) =  (aa1133 - aa1331)
        ainv(2,3) = -(aa1123 - aa1321)
        
        ainv(3,1) =  (aa2132 - aa2231)
        ainv(3,2) = -(aa1132 - aa1231)
        ainv(3,3) =  (aa1122 - aa1221)

        return

    end  subroutine

    pure subroutine mul_2x2(amat,adim,bmat,bdim, &
        &                   scal,cmat,cdim)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: bdim
        real*8 , intent(in)    :: bmat(bdim,*)
        real*8 , intent(in)    :: scal
        integer, intent(in)    :: cdim
        real*8 , intent(inout) :: cmat(cdim,*)

    !-------------------------------- C = C + scal * A * B !

        if (scal .eq. +1.e0_8) then

        cmat(1,1) = cmat(1,1) & 
             + ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) )
        cmat(2,1) = cmat(2,1) & 
             + ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) )

        cmat(1,2) = cmat(1,2) & 
             + ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) )
        cmat(2,2) = cmat(2,2) & 
             + ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) )

        else &
        if (scal .eq. -1.e0_8) then

        cmat(1,1) = cmat(1,1) & 
             - ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) )
        cmat(2,1) = cmat(2,1) & 
             - ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) )

        cmat(1,2) = cmat(1,2) & 
             - ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) )
        cmat(2,2) = cmat(2,2) & 
             - ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) )

        else 

        cmat(1,1) = cmat(1,1) + & 
        scal * ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) )
        cmat(2,1) = cmat(2,1) + & 
        scal * ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) )

        cmat(1,2) = cmat(1,2) + & 
        scal * ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) )
        cmat(2,2) = cmat(2,2) + & 
        scal * ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) )

        end if

        return

    end  subroutine

    pure subroutine mul_3x3(amat,adim,bmat,bdim, &
        &                   scal,cmat,cdim)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: bdim
        real*8 , intent(in)    :: bmat(bdim,*)
        real*8 , intent(in)    :: scal
        integer, intent(in)    :: cdim
        real*8 , intent(inout) :: cmat(cdim,*)

    !-------------------------------- C = C + scal * A * B !

        if (scal .eq. +1.e0_8) then

        cmat(1,1) = cmat(1,1) &
             + ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) &
               + amat(1,3) * bmat(3,1) )
        cmat(2,1) = cmat(2,1) &
             + ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) &
               + amat(2,3) * bmat(3,1) )
        cmat(3,1) = cmat(3,1) &
             + ( amat(3,1) * bmat(1,1) &
               + amat(3,2) * bmat(2,1) &
               + amat(3,3) * bmat(3,1) )

        cmat(1,2) = cmat(1,2) &
             + ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) &
               + amat(1,3) * bmat(3,2) )
        cmat(2,2) = cmat(2,2) &
             + ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) &
               + amat(2,3) * bmat(3,2) )
        cmat(3,2) = cmat(3,2) &
             + ( amat(3,1) * bmat(1,2) &
               + amat(3,2) * bmat(2,2) &
               + amat(3,3) * bmat(3,2) )

        cmat(1,3) = cmat(1,3) &
             + ( amat(1,1) * bmat(1,3) &
               + amat(1,2) * bmat(2,3) &
               + amat(1,3) * bmat(3,3) )
        cmat(2,3) = cmat(2,3) &
             + ( amat(2,1) * bmat(1,3) &
               + amat(2,2) * bmat(2,3) &
               + amat(2,3) * bmat(3,3) )
        cmat(3,3) = cmat(3,3) &
             + ( amat(3,1) * bmat(1,3) &
               + amat(3,2) * bmat(2,3) &
               + amat(3,3) * bmat(3,3) )

        else &
        if (scal .eq. -1.e0_8) then

        cmat(1,1) = cmat(1,1) &
             - ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) &
               + amat(1,3) * bmat(3,1) )
        cmat(2,1) = cmat(2,1) &
             - ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) &
               + amat(2,3) * bmat(3,1) )
        cmat(3,1) = cmat(3,1) &
             - ( amat(3,1) * bmat(1,1) &
               + amat(3,2) * bmat(2,1) &
               + amat(3,3) * bmat(3,1) )

        cmat(1,2) = cmat(1,2) &
             - ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) &
               + amat(1,3) * bmat(3,2) )
        cmat(2,2) = cmat(2,2) &
             - ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) &
               + amat(2,3) * bmat(3,2) )
        cmat(3,2) = cmat(3,2) &
             - ( amat(3,1) * bmat(1,2) &
               + amat(3,2) * bmat(2,2) &
               + amat(3,3) * bmat(3,2) )

        cmat(1,3) = cmat(1,3) &
             - ( amat(1,1) * bmat(1,3) &
               + amat(1,2) * bmat(2,3) &
               + amat(1,3) * bmat(3,3) )
        cmat(2,3) = cmat(2,3) &
             - ( amat(2,1) * bmat(1,3) &
               + amat(2,2) * bmat(2,3) &
               + amat(2,3) * bmat(3,3) )
        cmat(3,3) = cmat(3,3) &
             - ( amat(3,1) * bmat(1,3) &
               + amat(3,2) * bmat(2,3) &
               + amat(3,3) * bmat(3,3) )

        else 

        cmat(1,1) = cmat(1,1) + & 
        scal * ( amat(1,1) * bmat(1,1) &
               + amat(1,2) * bmat(2,1) &
               + amat(1,3) * bmat(3,1) )
        cmat(2,1) = cmat(2,1) + & 
        scal * ( amat(2,1) * bmat(1,1) &
               + amat(2,2) * bmat(2,1) &
               + amat(2,3) * bmat(3,1) )
        cmat(3,1) = cmat(3,1) + & 
        scal * ( amat(3,1) * bmat(1,1) &
               + amat(3,2) * bmat(2,1) &
               + amat(3,3) * bmat(3,1) )

        cmat(1,2) = cmat(1,2) + & 
        scal * ( amat(1,1) * bmat(1,2) &
               + amat(1,2) * bmat(2,2) &
               + amat(1,3) * bmat(3,2) )
        cmat(2,2) = cmat(2,2) + & 
        scal * ( amat(2,1) * bmat(1,2) &
               + amat(2,2) * bmat(2,2) &
               + amat(2,3) * bmat(3,2) )
        cmat(3,2) = cmat(3,2) + & 
        scal * ( amat(3,1) * bmat(1,2) &
               + amat(3,2) * bmat(2,2) &
               + amat(3,3) * bmat(3,2) )

        cmat(1,3) = cmat(1,3) + & 
        scal * ( amat(1,1) * bmat(1,3) &
               + amat(1,2) * bmat(2,3) &
               + amat(1,3) * bmat(3,3) )
        cmat(2,3) = cmat(2,3) + & 
        scal * ( amat(2,1) * bmat(1,3) &
               + amat(2,2) * bmat(2,3) &
               + amat(2,3) * bmat(3,3) )
        cmat(3,3) = cmat(3,3) + & 
        scal * ( amat(3,1) * bmat(1,3) &
               + amat(3,2) * bmat(2,3) &
               + amat(3,3) * bmat(3,3) )

        end if

        return

    end  subroutine

    pure subroutine slv_2x2(amat,adim,vrhs,vdim, & 
        &                   nrhs,fEPS,okay)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: vdim
        real*8 , intent(inout) :: vrhs(vdim,*)
        integer, intent(in)    :: nrhs
        real*8 , intent(in)    :: fEPS
        logical, intent(inout) :: okay

    !------------------------------------------- variables !
        real*8                 :: ainv(2,2)
        real*8                 :: adet
        real*8                 :: vtmp(  2)
        integer                :: irhs

        integer, parameter     :: LDIM = 2

    !---------------------------------------- calc. inv(A) !

        call inv_2x2(amat,adim,ainv,LDIM,&
                     adet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return
    
    !---------------------------------------- v = A^-1 * v !
    
        do irhs = 1, nrhs

        vtmp(1) =  &
          + ( & 
        ainv(1, 1) *   vrhs(1,irhs) &
      + ainv(1, 2) *   vrhs(2,irhs) &
            ) / adet

        vtmp(2) =  &
          + ( & 
        ainv(2, 1) *   vrhs(1,irhs) &
      + ainv(2, 2) *   vrhs(2,irhs) &
            ) / adet

        vrhs(1,irhs) = vtmp(1)
        vrhs(2,irhs) = vtmp(2)

        end do

        return

    end  subroutine

    pure subroutine slv_3x3(amat,adim,vrhs,vdim, & 
        &                   nrhs,fEPS,okay)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: vdim
        real*8 , intent(inout) :: vrhs(vdim,*)
        integer, intent(in)    :: nrhs
        real*8 , intent(in)    :: fEPS
        logical, intent(inout) :: okay

    !------------------------------------------- variables !
        real*8                 :: ainv(3,3)
        real*8                 :: adet
        real*8                 :: vtmp(  3)
        integer                :: irhs

        integer, parameter     :: LDIM = 3

    !---------------------------------------- calc. inv(A) !

        call inv_3x3(amat,adim,ainv,LDIM,&
                     adet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return
    
    !---------------------------------------- v = A^-1 * v !
    
        do irhs = 1, nrhs

        vtmp(1) =  &
          + ( & 
        ainv(1, 1) *   vrhs(1,irhs) &
      + ainv(1, 2) *   vrhs(2,irhs) &
      + ainv(1, 3) *   vrhs(3,irhs) &
            ) / adet

        vtmp(2) =  &
          + ( & 
        ainv(2, 1) *   vrhs(1,irhs) &
      + ainv(2, 2) *   vrhs(2,irhs) &
      + ainv(2, 3) *   vrhs(3,irhs) &
            ) / adet

        vtmp(3) =  &
          + ( & 
        ainv(3, 1) *   vrhs(1,irhs) &
      + ainv(3, 2) *   vrhs(2,irhs) &
      + ainv(3, 3) *   vrhs(3,irhs) &
            ) / adet

        vrhs(1,irhs) = vtmp(1)
        vrhs(2,irhs) = vtmp(2)
        vrhs(3,irhs) = vtmp(3)

        end do

        return

    end  subroutine

    pure subroutine slv_4x4(amat,adim,vrhs,vdim, & 
        &                   nrhs,fEPS,okay)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: vdim
        real*8 , intent(inout) :: vrhs(vdim,*)
        integer, intent(in)    :: nrhs
        real*8 , intent(in)    :: fEPS
        logical, intent(inout) :: okay

    !------------------------------------------- variables !
        real*8                 :: ainv(2,2)
        real*8                 :: lmat(2,2)
        real*8                 :: umat(2,2)
        real*8                 :: smat(2,2)
        real*8                 :: sinv(2,2)
        real*8                 :: adet,sdet        
        real*8                 :: vtmp(  2)
        integer                :: irhs

        integer, parameter     :: LDIM = 2

    !---------------------- form a block LDU factorisation !

        call inv_2x2(amat(1,1),adim,ainv,LDIM, &
                     adet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return

    !---------------------------------------- L = C * A^-1 !

        lmat(1,1) = +0.e0_8
        lmat(1,2) = +0.e0_8
        lmat(2,1) = +0.e0_8
        lmat(2,2) = +0.e0_8

        call mul_2x2(amat(3,1),adim,ainv,LDIM, &
                    +1.e0_8,lmat,LDIM)
        
    !---------------------------------------- U = A^-1 * B !       

        umat(1,1) = +0.e0_8
        umat(1,2) = +0.e0_8
        umat(2,1) = +0.e0_8
        umat(2,2) = +0.e0_8
 
        call mul_2x2(ainv,LDIM,amat(1,3),adim, &
                    +1.e0_8,umat,LDIM)
        
    !-------------------------------- S = D - C * A^-1 * B !

        smat(1,1) = amat(3,3)
        smat(1,2) = amat(3,4)
        smat(2,1) = amat(4,3)
        smat(2,2) = amat(4,4)

        call mul_2x2(lmat,LDIM,amat(1,3),adim, &
                    -1.e0_8/adet,smat,LDIM)

        call inv_2x2(smat,LDIM,sinv,LDIM,sdet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return

    !-------------------------------- back-solve LDU = rhs !

        do irhs = 1, nrhs

    !---------------------------------------- solve L part !

        vrhs(3,irhs) = vrhs(3,irhs) &
          - ( &
        lmat(1, 1) *   vrhs(1,irhs) &
      + lmat(1, 2) *   vrhs(2,irhs) &
            ) / adet

        vrhs(4,irhs) = vrhs(4,irhs) &
          - ( &
        lmat(2, 1) *   vrhs(1,irhs) &
      + lmat(2, 2) *   vrhs(2,irhs) &
            ) / adet
        
    !---------------------------------------- solve D part !

        vtmp(1) =  &
          + ( & 
        ainv(1, 1) *   vrhs(1,irhs) &
      + ainv(1, 2) *   vrhs(2,irhs) &
            ) / adet

        vtmp(2) =  &
          + ( & 
        ainv(2, 1) *   vrhs(1,irhs) &
      + ainv(2, 2) *   vrhs(2,irhs) &
            ) / adet

        vrhs(1,irhs) = vtmp(1)
        vrhs(2,irhs) = vtmp(2)

        vtmp(1) =  &
          + ( &
        sinv(1, 1) *   vrhs(3,irhs) &
      + sinv(1, 2) *   vrhs(4,irhs) & 
            ) / sdet

        vtmp(2) =  &
          + ( &
        sinv(2, 1) *   vrhs(3,irhs) &
      + sinv(2, 2) *   vrhs(4,irhs) & 
            ) / sdet

        vrhs(3,irhs) = vtmp(1)
        vrhs(4,irhs) = vtmp(2)

    !---------------------------------------- solve U part !

        vrhs(1,irhs) = vrhs(1,irhs) &
          - ( &
        umat(1, 1) *   vrhs(3,irhs) &
      + umat(1, 2) *   vrhs(4,irhs) & 
            ) / adet

        vrhs(2,irhs) = vrhs(2,irhs) &
          - ( &
        umat(2, 1) *   vrhs(3,irhs) &
      + umat(2, 2) *   vrhs(4,irhs) & 
            ) / adet

        end do

        return

    end  subroutine

    pure subroutine slv_6x6(amat,adim,vrhs,vdim, & 
        &                   nrhs,fEPS,okay)

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)    :: adim
        real*8 , intent(in)    :: amat(adim,*)
        integer, intent(in)    :: vdim
        real*8 , intent(inout) :: vrhs(vdim,*)
        integer, intent(in)    :: nrhs
        real*8 , intent(in)    :: fEPS
        logical, intent(inout) :: okay

    !------------------------------------------- variables !
        real*8                 :: ainv(3,3)
        real*8                 :: lmat(3,3)
        real*8                 :: umat(3,3)
        real*8                 :: smat(3,3)
        real*8                 :: sinv(3,3)
        real*8                 :: adet,sdet        
        real*8                 :: vtmp(  3)
        integer                :: irhs

        integer, parameter     :: LDIM = 3

    !---------------------- form a block LDU factorisation !

        call inv_3x3(amat(1,1),adim,ainv,LDIM, &
                     adet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return

    !---------------------------------------- L = C * A^-1 !

        lmat(1,1) = +0.e0_8
        lmat(1,2) = +0.e0_8
        lmat(1,3) = +0.e0_8
        lmat(2,1) = +0.e0_8
        lmat(2,2) = +0.e0_8
        lmat(2,3) = +0.e0_8
        lmat(3,1) = +0.e0_8
        lmat(3,2) = +0.e0_8
        lmat(3,3) = +0.e0_8

        call mul_3x3(amat(4,1),adim,ainv,LDIM, &
                    +1.e0_8,lmat,LDIM)
        
    !---------------------------------------- U = A^-1 * B !       

        umat(1,1) = +0.e0_8
        umat(1,2) = +0.e0_8
        umat(1,3) = +0.e0_8
        umat(2,1) = +0.e0_8
        umat(2,2) = +0.e0_8
        umat(2,3) = +0.e0_8
        umat(3,1) = +0.e0_8
        umat(3,2) = +0.e0_8
        umat(3,3) = +0.e0_8
 
        call mul_3x3(ainv,LDIM,amat(1,4),adim, &
                    +1.e0_8,umat,LDIM)
        
    !-------------------------------- S = D - C * A^-1 * B !

        smat(1,1) = amat(4,4)
        smat(1,2) = amat(4,5)
        smat(1,3) = amat(4,6)
        smat(2,1) = amat(5,4)
        smat(2,2) = amat(5,5)
        smat(2,3) = amat(5,6)
        smat(3,1) = amat(6,4)
        smat(3,2) = amat(6,5)
        smat(3,3) = amat(6,6)

        call mul_3x3(lmat,LDIM,amat(1,4),adim, &
                    -1.e0_8/adet,smat,LDIM)

        call inv_3x3(smat,LDIM,sinv,LDIM,sdet)

        okay = (abs(adet) .gt. fEPS)
        
        if (okay.eqv..false.) return

    !-------------------------------- back-solve LDU = rhs !

        do irhs = 1, nrhs

    !---------------------------------------- solve L part !

        vrhs(4,irhs) = vrhs(4,irhs) &
          - ( &
        lmat(1, 1) *   vrhs(1,irhs) &
      + lmat(1, 2) *   vrhs(2,irhs) &
      + lmat(1, 3) *   vrhs(3,irhs) &
            ) / adet

        vrhs(5,irhs) = vrhs(5,irhs) &
          - ( &
        lmat(2, 1) *   vrhs(1,irhs) &
      + lmat(2, 2) *   vrhs(2,irhs) &
      + lmat(2, 3) *   vrhs(3,irhs) &
            ) / adet

        vrhs(6,irhs) = vrhs(6,irhs) &
          - ( &
        lmat(3, 1) *   vrhs(1,irhs) &
      + lmat(3, 2) *   vrhs(2,irhs) &
      + lmat(3, 3) *   vrhs(3,irhs) &
            ) / adet
        
    !---------------------------------------- solve D part !

        vtmp(1) =  &
          + ( & 
        ainv(1, 1) *   vrhs(1,irhs) &
      + ainv(1, 2) *   vrhs(2,irhs) &
      + ainv(1, 3) *   vrhs(3,irhs) &
            ) / adet

        vtmp(2) =  &
          + ( & 
        ainv(2, 1) *   vrhs(1,irhs) &
      + ainv(2, 2) *   vrhs(2,irhs) &
      + ainv(2, 3) *   vrhs(3,irhs) &
            ) / adet

        vtmp(3) =  &
          + ( & 
        ainv(3, 1) *   vrhs(1,irhs) &
      + ainv(3, 2) *   vrhs(2,irhs) &
      + ainv(3, 3) *   vrhs(3,irhs) &
            ) / adet

        vrhs(1,irhs) = vtmp(1)
        vrhs(2,irhs) = vtmp(2)
        vrhs(3,irhs) = vtmp(3)

        vtmp(1) =  &
          + ( &
        sinv(1, 1) *   vrhs(4,irhs) &
      + sinv(1, 2) *   vrhs(5,irhs) &
      + sinv(1, 3) *   vrhs(6,irhs) & 
            ) / sdet

        vtmp(2) =  &
          + ( &
        sinv(2, 1) *   vrhs(4,irhs) &
      + sinv(2, 2) *   vrhs(5,irhs) &
      + sinv(2, 3) *   vrhs(6,irhs) & 
            ) / sdet

        vtmp(3) =  &
          + ( &
        sinv(3, 1) *   vrhs(4,irhs) &
      + sinv(3, 2) *   vrhs(5,irhs) &
      + sinv(3, 3) *   vrhs(6,irhs) & 
            ) / sdet

        vrhs(4,irhs) = vtmp(1)
        vrhs(5,irhs) = vtmp(2)
        vrhs(6,irhs) = vtmp(3)

    !---------------------------------------- solve U part !

        vrhs(1,irhs) = vrhs(1,irhs) &
          - ( &
        umat(1, 1) *   vrhs(4,irhs) &
      + umat(1, 2) *   vrhs(5,irhs) &
      + umat(1, 3) *   vrhs(6,irhs) & 
            ) / adet

        vrhs(2,irhs) = vrhs(2,irhs) &
          - ( &
        umat(2, 1) *   vrhs(4,irhs) &
      + umat(2, 2) *   vrhs(5,irhs) &
      + umat(2, 3) *   vrhs(6,irhs) & 
            ) / adet

        vrhs(3,irhs) = vrhs(3,irhs) &
          - ( &
        umat(3, 1) *   vrhs(4,irhs) &
      + umat(3, 2) *   vrhs(5,irhs) &
      + umat(3, 3) *   vrhs(6,irhs) & 
            ) / adet

        end do

        return

    end  subroutine





    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! PBC.h90: setup polynomial B.C.'s at domain endpoints.
    !
    ! Darren Engwirda 
    ! 09-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine pbc(npos,nvar,ndof,delx, &
        &          fdat,bcon,edge,dfdx, &
        &          iend,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCON  boundary condition data for endpoint .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! IEND  domain endpoint, IEND < +0 for lower end-point 
    !       and IEND > +0 for upper endpoint .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        integer, intent( in) :: iend
        real*8 , intent( in) :: dmin
        type(rcon_ends), intent(in) :: bcon(:)

    !------------------------------------------- variables !
        integer :: ivar,nlse,nval,nslp

        nlse = 0 ; nval = 0 ; nslp = 0

        do  ivar = +1, nvar

            select case (bcon(ivar)%bcopt)
    !------------------------------------------- find BC's !
            case(bcon_loose)
                nlse = nlse + 1
            
            case(bcon_value)
                nval = nval + 1
            
            case(bcon_slope)
                nslp = nslp + 1

            end select

        end do

    !---------------------------- setup "lower" conditions !

        if (iend.lt.+0) then

        if (nlse.gt.+0) then
    !---------------------------- setup "unset" conditions !
        call lbc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_loose    , &
    &            edge,dfdx,dmin)

        end if

        if (nval.gt.+0) then
    !---------------------------- setup "value" conditions !
        call lbc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_value    , &
    &            edge,dfdx,dmin)

        end if

        if (nslp.gt.+0) then
    !---------------------------- setup "slope" conditions !
        call lbc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_slope    , &
    &            edge,dfdx,dmin)

        end if
    
        end if
        
    !---------------------------- setup "upper" conditions !

        if (iend.gt.+0) then

        if (nlse.gt.+0) then
    !---------------------------- setup "unset" conditions !
        call ubc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_loose    , &
    &            edge,dfdx,dmin)

        end if

        if (nval.gt.+0) then
    !---------------------------- setup "value" conditions !
        call ubc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_value    , &
    &            edge,dfdx,dmin)

        end if

        if (nslp.gt.+0) then
    !---------------------------- setup "slope" conditions !
        call ubc(npos,nvar,ndof, &
    &            delx,fdat,bcon, &
    &            bcon_slope    , &
    &            edge,dfdx,dmin)

        end if
    
        end if
        
        return

    end  subroutine
       
    ! LBC: impose a single B.C.-type at the lower endpoint !
    
    subroutine lbc(npos,nvar,ndof,delx, &
        &          fdat,bcon,bopt,edge, &
        &          dfdx,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCON  boundary condition data for endpoint .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: bopt
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        real*8 , intent( in) :: dmin
        type(rcon_ends), intent(in) :: bcon(:)

    !------------------------------------------- variables !
        integer :: ivar,idof,isel, &
        &          head,tail,nsel
        logical :: okay        
        real*8  :: xhat
        real*8  :: delh(-1:+1)
        real*8  :: xmap(-1:+2)
        real*8  :: bvec(+3,-1:+2)
        real*8  :: gvec(+3,-1:+2)
        real*8  :: cmat(+3,+3)
        real*8  :: fhat(+3, nvar)
        real*8  :: eval(-1:+2)
        real*8  :: gval(-1:+2)
        
        integer, parameter :: NSIZ = +3
        real*8 , parameter :: ZERO = +1.e-14  

        head = +2; tail = npos - 2

        if (size(delx).gt.+1) then

    !------------------ mean grid spacing about ii-th cell !

        xhat = max(delx(head),dmin) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(head-1)
        delh(+0) = delx(head+0)
        delh(+1) = delx(head+1)

        else
        
    !------------------ mean grid spacing about ii-th cell !

        xhat = max(delx(  +1),dmin) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(    +1)
        delh(+0) = delx(    +1)
        delh(+1) = delx(    +1)
        
        end if

    !---------- local coordinate mapping for stencil edges !

        xmap(-1) =-(delh(-1) + &
        &           delh(+0)*0.5d0)/xhat
        xmap(+0) = -1.e0_8
        xmap(+1) = +1.e0_8
        xmap(+2) = (delh(+1) + &
        &           delh(+0)*0.5d0)/xhat

    !------------ linear system: lhs reconstruction matrix !

        select case(bopt )
        case( bcon_loose )

        call bfun1d(-1,+3,xmap(-1),bvec(:,-1))
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))
        call bfun1d(-1,+3,xmap(+2),bvec(:,+2))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+0) &
        &                - bvec(idof,-1)
            cmat(2,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)
            cmat(3,idof) = bvec(idof,+2) &
        &                - bvec(idof,+1)

        end do

        case( bcon_value )

        call bfun1d(+0,+3,xmap(-1),gvec(:,-1))
        
        call bfun1d(-1,+3,xmap(-1),bvec(:,-1))
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+0) &
        &                - bvec(idof,-1)
            cmat(2,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)

            cmat(3,idof) = gvec(idof,-1)

        end do

        case( bcon_slope )

        call bfun1d(+1,+3,xmap(-1),gvec(:,-1))
        
        call bfun1d(-1,+3,xmap(-1),bvec(:,-1))
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+0) &
        &                - bvec(idof,-1)
            cmat(2,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)

            cmat(3,idof) = gvec(idof,-1)

        end do

        end select

    !------------ linear system: rhs reconstruction vector !

        isel = 0 ; nsel = 0

        select case( bopt )
        case ( bcon_loose )

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_loose)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(-1) * &
        &       fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = delh(+0) * &
        &       fdat(1,ivar,head+0) / xhat
            fhat(3,isel) = delh(+1) * &
        &       fdat(1,ivar,head+1) / xhat
 
        end if
        
        end do       
 
        case ( bcon_value )
        
        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_value)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(-1) * &
        &       fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = delh(+0) * &
        &       fdat(1,ivar,head+0) / xhat

            fhat(3,isel) = bcon(ivar)%value
 
        end if
        
        end do
        
        case ( bcon_slope )

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_slope)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(-1) * &
        &       fdat(1,ivar,head-1) / xhat
            fhat(2,isel) = delh(+0) * &
        &       fdat(1,ivar,head+0) / xhat

            fhat(3,isel) = &
        &       bcon(ivar)%slope * xhat
 
        end if
        
        end do
        
        end select

    !------------------------- factor/solve linear systems !

        call slv_3x3(cmat,NSIZ,fhat , &
        &            NSIZ,nvar,       &
        &            ZERO*dmin,okay)

        if (okay .eqv..false.) then


        end if

        if (okay .eqv. .true.) then

    !------------- extrapolate values/slopes at lower edge !

        isel  = +0

        call bfun1d(+0,+3,xmap(-1),bvec(:,-1))
        call bfun1d(+0,+3,xmap(+0),bvec(:,+0))
        call bfun1d(+0,+3,xmap(+1),bvec(:,+1))
        
        call bfun1d(+1,+3,xmap(-1),gvec(:,-1))
        call bfun1d(+1,+3,xmap(+0),gvec(:,+0))
        call bfun1d(+1,+3,xmap(+1),gvec(:,+1))

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            isel = isel  + 1

            eval(-1) = dot_product( &
        &       bvec(:,-1),fhat(:,isel))
            eval(+0) = dot_product( &
        &       bvec(:,+0),fhat(:,isel))        
            eval(+1) = dot_product( &
        &       bvec(:,+1),fhat(:,isel))
        
            gval(-1) = dot_product( &
        &       gvec(:,-1),fhat(:,isel))
            gval(+0) = dot_product( &
        &       gvec(:,+0),fhat(:,isel))        
            gval(+1) = dot_product( &
        &       gvec(:,+1),fhat(:,isel))

            edge(ivar,head-1) = eval(-1)
            edge(ivar,head+0) = eval(+0)
            edge(ivar,head+1) = eval(+1)

            dfdx(ivar,head-1) = gval(-1) &
        &                     / xhat
            dfdx(ivar,head+0) = gval(+0) &
        &                     / xhat
            dfdx(ivar,head+1) = gval(+1) &
        &                     / xhat

        end if

        end do
   
        else

    !------------- low-order if re-con. matrix is singular !

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            eval(-1) = &
        &   fdat(1,ivar,head-1) * 1.e0_8
            eval(+0) = &
        &   fdat(1,ivar,head-1) * .5d0 + &
        &   fdat(1,ivar,head+0) * .5d0
            eval(+1) = &
        &   fdat(1,ivar,head+0) * .5d0 + &
        &   fdat(1,ivar,head+1) * .5d0
   
            gval(-1) = &
        &   fdat(1,ivar,head+0) * .5d0 - &
        &   fdat(1,ivar,head-1) * .5d0
            gval(+0) = &
        &   fdat(1,ivar,head+0) * .5d0 - &
        &   fdat(1,ivar,head-1) * .5d0
            gval(+1) = &
        &   fdat(1,ivar,head+1) * .5d0 - &
        &   fdat(1,ivar,head+0) * .5d0
   
            edge(ivar,head-1) = eval(-1)
            edge(ivar,head+0) = eval(+0)
            edge(ivar,head+1) = eval(+1)

            dfdx(ivar,head-1) = gval(-1) &
        &                     / xhat
            dfdx(ivar,head+0) = gval(+0) &
        &                     / xhat
            dfdx(ivar,head+1) = gval(+1) &
        &                     / xhat

        end if
        
        end do
    
        end if
        
        return

    end  subroutine
    
    ! UBC: impose a single B.C.-type at the upper endpoint !
    
    subroutine ubc(npos,nvar,ndof,delx, &
        &          fdat,bcon,bopt,edge, &
        &          dfdx,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCON  boundary condition data for endpoint .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        integer, intent( in) :: bopt
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        real*8 , intent( in) :: dmin
        type(rcon_ends), intent(in) :: bcon(:)

    !------------------------------------------- variables !
        integer :: ivar,idof,isel, &
        &          head,tail,nsel
        logical :: okay
        real*8  :: xhat
        real*8  :: delh(-1:+1)
        real*8  :: xmap(-1:+2)
        real*8  :: bvec(+3,-1:+2)
        real*8  :: gvec(+3,-1:+2)
        real*8  :: cmat(+3,+3)
        real*8  :: fhat(+3, nvar)
        real*8  :: eval(-1:+2)
        real*8  :: gval(-1:+2)

        integer, parameter :: NSIZ = +3
        real*8 , parameter :: ZERO = +1.e-14

        head = +2; tail = npos - 2

        if (size(delx).gt.+1) then

    !------------------ mean grid spacing about ii-th cell !

        xhat = max(delx(tail),dmin) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(tail-1)
        delh(+0) = delx(tail+0)
        delh(+1) = delx(tail+1)

        else
        
    !------------------ mean grid spacing about ii-th cell !

        xhat = max(delx(  +1),dmin) * 0.5d+0

    !------------------ grid spacing for all stencil cells !

        delh(-1) = delx(    +1)
        delh(+0) = delx(    +1)
        delh(+1) = delx(    +1)
        
        end if

    !---------- local coordinate mapping for stencil edges !

        xmap(-1) =-(delh(-1) + &
        &           delh(+0)*0.5d0)/xhat
        xmap(+0) = -1.e0_8
        xmap(+1) = +1.e0_8
        xmap(+2) = (delh(+1) + &
        &           delh(+0)*0.5d0)/xhat

    !------------ linear system: lhs reconstruction matrix !

        select case(bopt )
        case( bcon_loose )

        call bfun1d(-1,+3,xmap(-1),bvec(:,-1))
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))
        call bfun1d(-1,+3,xmap(+2),bvec(:,+2))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+0) &
        &                - bvec(idof,-1)
            cmat(2,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)
            cmat(3,idof) = bvec(idof,+2) &
        &                - bvec(idof,+1)

        end do

        case( bcon_value )

        call bfun1d(+0,+3,xmap(+2),gvec(:,+2))
        
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))
        call bfun1d(-1,+3,xmap(+2),bvec(:,+2))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)
            cmat(2,idof) = bvec(idof,+2) &
        &                - bvec(idof,+1)

            cmat(3,idof) = gvec(idof,+2)

        end do

        case( bcon_slope )

        call bfun1d(+1,+3,xmap(+2),gvec(:,+2))
        
        call bfun1d(-1,+3,xmap(+0),bvec(:,+0))
        call bfun1d(-1,+3,xmap(+1),bvec(:,+1))
        call bfun1d(-1,+3,xmap(+2),bvec(:,+2))

        do  idof = +1 , +3

            cmat(1,idof) = bvec(idof,+1) &
        &                - bvec(idof,+0)
            cmat(2,idof) = bvec(idof,+2) &
        &                - bvec(idof,+1)

            cmat(3,idof) = gvec(idof,+2)

        end do

        end select

    !------------ linear system: rhs reconstruction vector !

        isel = 0 ; nsel = 0

        select case( bopt )
        case ( bcon_loose )

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_loose)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(-1) * &
        &       fdat(1,ivar,tail-1) / xhat
            fhat(2,isel) = delh(+0) * &
        &       fdat(1,ivar,tail+0) / xhat
            fhat(3,isel) = delh(+1) * &
        &       fdat(1,ivar,tail+1) / xhat
 
        end if
        
        end do       
 
        case ( bcon_value )
        
        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_value)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(+0) * &
        &       fdat(1,ivar,tail+0) / xhat
            fhat(2,isel) = delh(+1) * &
        &       fdat(1,ivar,tail+1) / xhat

            fhat(3,isel) = bcon(ivar)%value
 
        end if
        
        end do
        
        case ( bcon_slope )

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bcon_slope)  then

            isel = isel + 1
            nsel = nsel + 1

            fhat(1,isel) = delh(+0) * &
        &       fdat(1,ivar,tail+0) / xhat
            fhat(2,isel) = delh(+1) * &
        &       fdat(1,ivar,tail+1) / xhat

            fhat(3,isel) = &
        &       bcon(ivar)%slope * xhat
 
        end if
        
        end do
        
        end select

    !------------------------- factor/solve linear systems !

        call slv_3x3(cmat,NSIZ,fhat , &
        &            NSIZ,nvar,       &
        &            ZERO*dmin,okay)

        if (okay .eqv..false.) then


        end if

        if (okay .eqv. .true.) then

    !------------- extrapolate values/slopes at lower edge !

        isel  = +0

        call bfun1d(+0,+3,xmap(+0),bvec(:,+0))
        call bfun1d(+0,+3,xmap(+1),bvec(:,+1))
        call bfun1d(+0,+3,xmap(+2),bvec(:,+2))
        
        call bfun1d(+1,+3,xmap(+0),gvec(:,+0))
        call bfun1d(+1,+3,xmap(+1),gvec(:,+1))
        call bfun1d(+1,+3,xmap(+2),gvec(:,+2))

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            isel = isel  + 1

            eval(+0) = dot_product( &
        &       bvec(:,+0),fhat(:,isel))
            eval(+1) = dot_product( &
        &       bvec(:,+1),fhat(:,isel))        
            eval(+2) = dot_product( &
        &       bvec(:,+2),fhat(:,isel))
        
            gval(+0) = dot_product( &
        &       gvec(:,+0),fhat(:,isel))
            gval(+1) = dot_product( &
        &       gvec(:,+1),fhat(:,isel))        
            gval(+2) = dot_product( &
        &       gvec(:,+2),fhat(:,isel))

            edge(ivar,tail+0) = eval(+0)
            edge(ivar,tail+1) = eval(+1)
            edge(ivar,tail+2) = eval(+2)

            dfdx(ivar,tail+0) = gval(+0) &
        &                     / xhat
            dfdx(ivar,tail+1) = gval(+1) &
        &                     / xhat
            dfdx(ivar,tail+2) = gval(+2) &
        &                     / xhat

        end if

        end do
   
        else

    !------------- low-order if re-con. matrix is singular !

        do  ivar = +1, nvar

        if (bcon(ivar)%bcopt.eq.bopt)  then

            eval(+0) = &
        &   fdat(1,ivar,tail-1) * .5d0 + &
        &   fdat(1,ivar,tail+0) * .5d0
            eval(+1) = &
        &   fdat(1,ivar,tail+0) * .5d0 + &
        &   fdat(1,ivar,tail+1) * .5d0
            eval(+2) = &
        &   fdat(1,ivar,tail+1) * 1.e0_8
   
            gval(+0) = &
        &   fdat(1,ivar,tail+0) * .5d0 - &
        &   fdat(1,ivar,tail-1) * .5d0
            gval(+1) = &
        &   fdat(1,ivar,tail+1) * .5d0 - &
        &   fdat(1,ivar,tail+0) * .5d0
            gval(+2) = &
        &   fdat(1,ivar,tail+1) * .5d0 - &
        &   fdat(1,ivar,tail+0) * .5d0
   
            edge(ivar,tail+0) = eval(+0)
            edge(ivar,tail+1) = eval(+1)
            edge(ivar,tail+2) = eval(+2)

            dfdx(ivar,tail+0) = gval(+0) &
        &                     / xhat
            dfdx(ivar,tail+1) = gval(+1) &
        &                     / xhat
            dfdx(ivar,tail+2) = gval(+2) &
        &                     / xhat

        end if
        
        end do
    
        end if
        
        return

    end  subroutine
    
    
    

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! P1E.h90: set edge estimates via degree-1 polynomials.
    !
    ! Darren Engwirda 
    ! 09-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine p1e(npos,nvar,ndof,delx, &
        &          fdat,bclo,bchi,edge, &
        &          dfdx,opts,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCLO  boundary condition at lower endpoint.
    ! BCHI  boundary condition at upper endpoint.
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bchi(:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        real*8 , intent( in) :: dmin
        class(rcon_opts), intent(in) :: opts

    !------------------------------------------- variables !
        integer :: ipos,ivar,head,tail
        real*8  :: dd10
        real*8  :: delh(-1:+0)

        head = +2; tail = npos-1

        if (npos.lt.2) return
        if (npos.eq.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = 1,nvar
        
            edge(ivar,1) = fdat(1,ivar,1)
            dfdx(ivar,1) = 0.e0_8
            
            edge(ivar,2) = fdat(1,ivar,1)
            dfdx(ivar,2) = 0.e0_8
            
        end do
        end if

        if (npos.le.2) return
   
    ! Reconstruct edge-centred 2nd-order polynomials. Com- !
    ! pute values/slopes at edges directly. Full-order ex- !
    ! trapolation at endpoints.
    
        if (size(delx).eq.+1) then
        
            do  ipos = head , tail
    
    !--------------- reconstruction: constant grid-spacing !
            
            dd10 = delx(+1) * 2.e0_8
            
            do  ivar = +1, nvar

                edge(ivar,ipos) = &
        &         + delx(+1) * &
        &       fdat(1,ivar,ipos-1) &
        &         + delx(+1) * &
        &       fdat(1,ivar,ipos+0)

                dfdx(ivar,ipos) = &
        &         - 2.0d+0 *  &
        &       fdat(1,ivar,ipos-1) &
        &         + 2.0d+0 *  &
        &       fdat(1,ivar,ipos+0)

                edge(ivar,ipos) = &
        &       edge(ivar,ipos) / dd10
                dfdx(ivar,ipos) = &
        &       dfdx(ivar,ipos) / dd10

            end do
            
            end do
            
        else
        
            do  ipos = head , tail
    
    !--------------- reconstruction: variable grid-spacing !
            
            delh(-1) = &
        &       max(delx(ipos-1),dmin)
            delh(+0) = &
        &       max(delx(ipos+0),dmin)

            dd10 = delh(-1)+delh(+0)
            
            do  ivar = +1, nvar

                edge(ivar,ipos) = &
        &         + delh(+0) * &
        &       fdat(1,ivar,ipos-1) &
        &         + delh(-1) * &
        &       fdat(1,ivar,ipos+0)

                dfdx(ivar,ipos) = &
        &         - 2.0d+0 *  &
        &       fdat(1,ivar,ipos-1) &
        &         + 2.0d+0 *  &
        &       fdat(1,ivar,ipos+0)

                edge(ivar,ipos) = &
        &       edge(ivar,ipos) / dd10
                dfdx(ivar,ipos) = &
        &       dfdx(ivar,ipos) / dd10

            end do
            
            end do  
            
        end if
    
    !------------- 1st-order value/slope BC's at endpoints !

        do  ivar = +1, nvar

            edge(ivar,head-1) = &
        &       fdat(+1,ivar,head-1)           
            edge(ivar,tail+1) = &
        &       fdat(+1,ivar,tail+0)
        
            dfdx(ivar,head-1) = 0.e0_8
            dfdx(ivar,tail+1) = 0.e0_8

        end do
    
        return
        
    end subroutine
    
    
    

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! P3E.h90: set edge estimates via degree-3 polynomials.
    !
    ! Darren Engwirda 
    ! 09-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !
    
    subroutine p3e(npos,nvar,ndof,delx, &
        &          fdat,bclo,bchi,edge, &
        &          dfdx,opts,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCLO  boundary condition at lower endpoint.
    ! BCHI  boundary condition at upper endpoint.
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bchi(:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        real*8 , intent( in) :: dmin
        class(rcon_opts), intent(in) :: opts

    !------------------------------------------- variables !
        integer :: ipos,ivar,idof,head,tail
        logical :: okay
        real*8  :: xhat,fEPS
        real*8  :: delh(-2:+1)
        real*8  :: xmap(-2:+2)
        real*8  :: fhat(+4, nvar)
        real*8  :: ivec(+4,-2:+2)
        real*8  :: cmat(+4,+4)
        
        integer, parameter :: NSIZ = +4
        real*8 , parameter :: ZERO = 1.e-14  
   
        head = +3 ; tail = npos - 2

        if (npos.le.4) then
    !----- default to reduced order if insufficient points !
            call p1e (npos,nvar,ndof, &
        &             delx,fdat,bclo, &
        &             bchi,edge,dfdx, &
        &             opts,dmin)
        end if

        if (npos.le.4) return

    !------ impose value/slope B.C.'s about lower endpoint !

        call pbc(npos,nvar,ndof,delx, &
        &        fdat,bclo,edge,dfdx, &
        &        -1  ,dmin)

    !------ impose value/slope B.C.'s about upper endpoint !

        call pbc(npos,nvar,ndof,delx, &
        &        fdat,bchi,edge,dfdx, &
        &        +1  ,dmin)

    ! Reconstruct edge-centred 4th-order polynomials. Com- !
    ! pute values/slopes at edges directly. Mid.-order ex- !
    ! trapolation at endpoints.                            !
    
        if (size(delx).eq.+1) then
        
            do  ipos = head , tail
    
    !--------------- reconstruction: constant grid-spacing !
            
            do  ivar = 1, nvar

                edge(ivar,ipos) = ( &
        &     -      1.e0_8 * &
        &       fdat(1,ivar,ipos-2) &
        &     +      7.e0_8 * &
        &       fdat(1,ivar,ipos-1) &
        &     +      7.e0_8 * &
        &       fdat(1,ivar,ipos+0) &
        &     -      1.e0_8 * &
        &       fdat(1,ivar,ipos+1) ) / 12.e0_8
            
                dfdx(ivar,ipos) = ( &
        &     +      1.e0_8 * &
        &       fdat(1,ivar,ipos-2) &
        &     -     15.e0_8 * &
        &       fdat(1,ivar,ipos-1) &
        &     +     15.e0_8 * &
        &       fdat(1,ivar,ipos+0) &
        &     -      1.e0_8 * &
        &       fdat(1,ivar,ipos+1) ) / 12.e0_8

                dfdx(ivar,ipos) = &
        &       dfdx(ivar,ipos) / delx(+1)

            end do
                      
            end do
        
        else
        
            fEPS     = ZERO * dmin

            do  ipos = head , tail
    
    !--------------- reconstruction: variable grid-spacing !
            
            delh(-2) = delx(ipos-2)
            delh(-1) = delx(ipos-1)
            delh(+0) = delx(ipos+0)
            delh(+1) = delx(ipos+1)

            xhat = .5d0 * max(delh(-1),dmin) + &
        &          .5d0 * max(delh(+0),dmin)

            xmap(-2) = -( delh(-2) &
        &              +  delh(-1) ) / xhat
            xmap(-1) = -  delh(-1)   / xhat
            xmap(+0) = +  0.e0_8
            xmap(+1) = +  delh(+0)   / xhat
            xmap(+2) = +( delh(+0) &
        &              +  delh(+1) ) / xhat
            
    !--------------------------- calc. integral basis vec. !

            do  idof = -2, +2

            ivec(1,idof) = &
        &       xmap(idof) ** 1 / 1.0d+0
            ivec(2,idof) = &
        &       xmap(idof) ** 2 / 2.0d+0
            ivec(3,idof) = &
        &       xmap(idof) ** 3 / 3.0d+0
            ivec(4,idof) = &
        &       xmap(idof) ** 4 / 4.0d+0
        
            end do

    !--------------------------- linear system: lhs matrix !

            do  idof = +1, +4

            cmat(1,idof) = ivec(idof,-1) &
        &                - ivec(idof,-2)
            cmat(2,idof) = ivec(idof,+0) &
        &                - ivec(idof,-1)
            cmat(3,idof) = ivec(idof,+1) &
        &                - ivec(idof,+0)
            cmat(4,idof) = ivec(idof,+2) &
        &                - ivec(idof,+1)

            end do

    !--------------------------- linear system: rhs vector !

            do  ivar = +1, nvar

            fhat(+1,ivar) = &
        &       delx(ipos-2) * &
        &   fdat(+1,ivar,ipos-2) / xhat
            fhat(+2,ivar) = &
        &       delx(ipos-1) * &
        &   fdat(+1,ivar,ipos-1) / xhat
            fhat(+3,ivar) = &
        &       delx(ipos+0) * &
        &   fdat(+1,ivar,ipos+0) / xhat
            fhat(+4,ivar) = &
        &       delx(ipos+1) * &
        &   fdat(+1,ivar,ipos+1) / xhat
        
            end do
            
    !------------------------- factor/solve linear systems !

            call slv_4x4(cmat,NSIZ,fhat, &
       &                 NSIZ,nvar,fEPS, &
       &                 okay)

            if (okay .eqv. .true.) then

            do  ivar =  +1, nvar

            edge(ivar,ipos) = fhat(1,ivar) 
                  
            dfdx(ivar,ipos) = fhat(2,ivar) &
        &                   / xhat
        
            end do

            else

    !------------------------- fallback if system singular !


            do  ivar =  +1, nvar

            edge(ivar,ipos) =   &
        &   fdat(1,ivar,ipos-1) * 0.5d+0 + &
        &   fdat(1,ivar,ipos-0) * 0.5d+0
        
            dfdx(ivar,ipos) =   &
        &   fdat(1,ivar,ipos-0) * 1.0d+0 - &
        &   fdat(1,ivar,ipos-1) * 1.0d+0
                  
            dfdx(ivar,ipos) =   &
        &       dfdx(ivar,ipos) / xhat
        
            end do

            end if
        
            end do
        
        end if

        return
        
    end subroutine




    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! P5E.h90: set edge estimates via degree-5 polynomials.
    !
    ! Darren Engwirda 
    ! 25-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine p5e(npos,nvar,ndof,delx, &
        &          fdat,bclo,bchi,edge, &
        &          dfdx,opts,dmin)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! BCLO  boundary condition at lower endpoint.
    ! BCHI  boundary condition at upper endpoint.
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! DMIN  min. grid-cell spacing thresh . 
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent( in) :: delx(:)
        real*8 , intent( in) :: fdat(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bchi(:)
        real*8 , intent(out) :: edge(:,:)
        real*8 , intent(out) :: dfdx(:,:)
        real*8 , intent( in) :: dmin
        class(rcon_opts), intent(in) :: opts

    !------------------------------------------- variables !
        integer :: ipos,ivar,idof,head,tail
        logical :: okay        
        real*8  :: xhat,fEPS
        real*8  :: delh(-3:+2)
        real*8  :: xmap(-3:+3)
        real*8  :: fhat(+6, nvar)
        real*8  :: ivec(+6,-3:+3)
        real*8  :: cmat(+6,+6)
        
        integer, parameter :: NSIZ = +6
        real*8 , parameter :: ZERO = 1.e-14   

        head = +4 ; tail = npos - 3

        if (npos.le.6) then
    !----- default to reduced order if insufficient points !
            call p3e (npos,nvar,ndof, &
        &             delx,fdat,bclo, &
        &             bchi,edge,dfdx, &
        &             opts,dmin)
        end if

        if (npos.le.6) return

    !------ impose value/slope B.C.'s about lower endpoint !

        call pbc(npos,nvar,ndof,delx, &
        &        fdat,bclo,edge,dfdx, &
        &        -1  ,dmin)

    !------ impose value/slope B.C.'s about upper endpoint !

        call pbc(npos,nvar,ndof,delx, &
        &        fdat,bchi,edge,dfdx, &
        &        +1  ,dmin)

    ! Reconstruct edge-centred 6th-order polynomials. Com- !
    ! pute values/slopes at edges directly. Mid.-order ex- !
    ! trapolation at endpoints.                            !
    
        if (size(delx).eq.+1) then
        
            do  ipos = head , tail
    
    !--------------- reconstruction: constant grid-spacing !
            
            do  ivar = 1, nvar

                edge(ivar,ipos) = &
        &     + ( 1.e0_8 / 60.e0_8) * & 
        &       fdat(1,ivar,ipos-3) &
        &     - ( 8.e0_8 / 60.e0_8) * &
        &       fdat(1,ivar,ipos-2) &
        &     + (37.e0_8 / 60.e0_8) * &
        &       fdat(1,ivar,ipos-1) &
        &     + (37.e0_8 / 60.e0_8) * &
        &       fdat(1,ivar,ipos+0) &
        &     - ( 8.e0_8 / 60.e0_8) * &
        &       fdat(1,ivar,ipos+1) &
        &     + ( 1.e0_8 / 60.e0_8) * &
        &       fdat(1,ivar,ipos+2)
            
                dfdx(ivar,ipos) = &
        &     - ( 1.e0_8 / 90.e0_8) * & 
        &       fdat(1,ivar,ipos-3) &
        &     + ( 5.e0_8 / 36.e0_8) * &
        &       fdat(1,ivar,ipos-2) &
        &     - (49.e0_8 / 36.e0_8) * &
        &       fdat(1,ivar,ipos-1) &
        &     + (49.e0_8 / 36.e0_8) * &
        &       fdat(1,ivar,ipos+0) &
        &     - ( 5.e0_8 / 36.e0_8) * &
        &       fdat(1,ivar,ipos+1) &
        &     + ( 1.e0_8 / 90.e0_8) * &
        &       fdat(1,ivar,ipos+2)

                dfdx(ivar,ipos) = &
                dfdx(ivar,ipos) / delx(+1)

            end do
                      
            end do
        
        else
        
            fEPS     = ZERO * dmin

            do  ipos = head , tail
    
    !--------------- reconstruction: variable grid-spacing !
            
            delh(-3) = &
        &       max(delx(ipos-3),dmin)
            delh(-2) = &
        &       max(delx(ipos-2),dmin)
            delh(-1) = &
        &       max(delx(ipos-1),dmin)
            delh(+0) = &
        &       max(delx(ipos+0),dmin)
            delh(+1) = &
        &       max(delx(ipos+1),dmin)
            delh(+2) = &
        &       max(delx(ipos+2),dmin)

            xhat = .5d0 * delh(-1) + &
        &          .5d0 * delh(+0)

            xmap(-3) = -( delh(-3) &
        &              +  delh(-2) &
        &              +  delh(-1) ) / xhat
            xmap(-2) = -( delh(-2) &
        &              +  delh(-1) ) / xhat
            xmap(-1) = -  delh(-1)   / xhat
            xmap(+0) = +  0.e0_8
            xmap(+1) = +  delh(+0)   / xhat
            xmap(+2) = +( delh(+0) &
        &              +  delh(+1) ) / xhat
            xmap(+3) = +( delh(+0) &
        &              +  delh(+1) &
        &              +  delh(+2) ) / xhat
            
    !--------------------------- calc. integral basis vec. !

            do  idof = -3, +3

            ivec(1,idof) = &
        &       xmap(idof) ** 1 / 1.0d+0
            ivec(2,idof) = &
        &       xmap(idof) ** 2 / 2.0d+0
            ivec(3,idof) = &
        &       xmap(idof) ** 3 / 3.0d+0
            ivec(4,idof) = &
        &       xmap(idof) ** 4 / 4.0d+0
            ivec(5,idof) = &
        &       xmap(idof) ** 5 / 5.0d+0
            ivec(6,idof) = &
        &       xmap(idof) ** 6 / 6.0d+0
        
            end do

    !--------------------------- linear system: lhs matrix !

            do  idof = +1, +6

            cmat(1,idof) = ivec(idof,-2) &
        &                - ivec(idof,-3)
            cmat(2,idof) = ivec(idof,-1) &
        &                - ivec(idof,-2)
            cmat(3,idof) = ivec(idof,+0) &
        &                - ivec(idof,-1)
            cmat(4,idof) = ivec(idof,+1) &
        &                - ivec(idof,+0)
            cmat(5,idof) = ivec(idof,+2) &
        &                - ivec(idof,+1)
            cmat(6,idof) = ivec(idof,+3) &
        &                - ivec(idof,+2)

            end do

    !--------------------------- linear system: rhs vector !

            do  ivar = +1, nvar

            fhat(+1,ivar) = &
        &       delx(ipos-3) * &
        &   fdat(+1,ivar,ipos-3) / xhat
            fhat(+2,ivar) = &
        &       delx(ipos-2) * &
        &   fdat(+1,ivar,ipos-2) / xhat
            fhat(+3,ivar) = &
        &       delx(ipos-1) * &
        &   fdat(+1,ivar,ipos-1) / xhat
            fhat(+4,ivar) = &
        &       delx(ipos+0) * &
        &   fdat(+1,ivar,ipos+0) / xhat
            fhat(+5,ivar) = &
        &       delx(ipos+1) * &
        &   fdat(+1,ivar,ipos+1) / xhat
            fhat(+6,ivar) = &
        &       delx(ipos+2) * &
        &   fdat(+1,ivar,ipos+2) / xhat
        
            end do
            
    !------------------------- factor/solve linear systems !

            call slv_6x6(cmat,NSIZ,fhat, &
       &                 NSIZ,nvar,fEPS, &
       &                 okay)

            if (okay .eqv. .true.) then

            do  ivar =  +1, nvar

            edge(ivar,ipos) = fhat(1,ivar) 
                  
            dfdx(ivar,ipos) = fhat(2,ivar) &
        &                   / xhat
        
            end do

            else

    !------------------------- fallback if system singular !


            do  ivar =  +1, nvar

            edge(ivar,ipos) =   &
        &   fdat(1,ivar,ipos-1) * 0.5d+0 + &
        &   fdat(1,ivar,ipos-0) * 0.5d+0
        
            dfdx(ivar,ipos) =   &
        &   fdat(1,ivar,ipos-0) * 0.5d+0 - &
        &   fdat(1,ivar,ipos-1) * 0.5d+0
                  
            dfdx(ivar,ipos) =   &
        &       dfdx(ivar,ipos) / xhat
        
            end do

            end if
        
            end do
        
        end if

        return
        
    end subroutine





    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! ROOT1D.h90: find the "roots" of degree-k polynomials.
    !
    ! Darren Engwirda 
    ! 25-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine roots_2(aa,bb,cc,xx,haveroot)

    !
    ! solve:: aa * xx**2 + bb * xx**1 + cc = +0.0 .
    !

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent( in) :: aa,bb,cc
        real*8 , intent(out) :: xx(1:2)
        logical, intent(out) :: haveroot

    !------------------------------------------- variables !
        real*8 :: sq,ia,a0,b0,c0,x0

        real*8, parameter :: rt = +1.e-14

        a0 = abs(aa)
        b0 = abs(bb)
        c0 = abs(cc)

        sq = bb * bb - 4.0d+0 * aa * cc

        if (sq .ge. 0.0d+0) then

            sq = sqrt (sq)

            xx(1) =  - bb + sq
            xx(2) =  - bb - sq

            x0 = max(abs(xx(1)), &
        &            abs(xx(2)))

            if (a0 .gt. (rt*x0)) then

    !-------------------------------------- degree-2 roots !
    
            haveroot =  .true.

            ia = 0.5d+0   / aa

            xx(1) = xx(1) * ia
            xx(2) = xx(2) * ia
            
            else &
        &   if (b0 .gt. (rt*c0)) then

    !-------------------------------------- degree-1 roots !

            haveroot =  .true.
            
            xx(1) =  - cc / bb
            xx(2) =  - cc / bb
            
            else
            
            haveroot = .false.
            
            end if

        else

            haveroot = .false.

        end if

        return

    end subroutine



    

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! PCM.h90: 1d piecewise constant reconstruction .
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine pcm(npos,nvar,ndof,fdat, &
        &               fhat)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar,ndof
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)

    !------------------------------------------- variables !
        integer:: ipos,ivar,idof

        do  ipos = +1, npos - 1
        do  ivar = +1, nvar + 0
        do  idof = +1, ndof + 0

            fhat(idof,ivar,ipos) = fdat(idof,ivar,ipos)

        end do
        end do
        end do
        
        return

    end subroutine
    
    
    

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! PLM.h90: a 1d, slope-limited piecewise linear method.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    pure subroutine plm(npos,nvar,ndof,delx, &
        &               fdat,fhat,dmin,ilim)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    ! ILIM  cell slope-limiting selection .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar
        integer, intent( in) :: ndof,ilim
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(:)
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)

        if (size(delx).gt.+1) then
        
    !------------------------------- variable grid-spacing !
        
            call plmv(npos,nvar,ndof,delx,&
        &             fdat,fhat,&
        &             dmin,ilim )
        
        else
        
    !------------------------------- constant grid-spacing !
        
            call plmc(npos,nvar,ndof,delx,&
        &             fdat,fhat,&
        &             dmin,ilim )
        
        end if

        return

    end  subroutine
    
    !------------------------- assemble PLM reconstruction !

    pure subroutine plmv(npos,nvar,ndof,delx, &
        &                fdat,fhat,dmin,ilim)

    !
    ! *this is the variable grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    ! ILIM  cell slope-limiting selection .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar
        integer, intent( in) :: ndof,ilim
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(:)
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)
        
    !------------------------------------------- variables !
        integer :: ipos,ivar,head,tail
        real*8  :: dfds(-1:+1)

        head = +1; tail = npos - 1

        if (npos.eq.2) then
    !----------------------- reduce order if small stencil !
        do  ivar = +1, nvar
            fhat(1,ivar,1) = &
        &   fdat(1,ivar,1)
            fhat(2,ivar,1) = 0.e+0
        end do
        end if

        if (npos.le.2) return
  
    !-------------------------------------- lower-endpoint !
        
        do  ivar = +1 , nvar-0

            call plsv( &
        &   fdat(1,ivar,head+0) , &
        &   delx(head+0), &
        &   fdat(1,ivar,head+0) , &
        &   delx(head+0), &
        &   fdat(1,ivar,head+1) , &
        &   delx(head+1), dfds)

            fhat(1,ivar,head) = &
        &   fdat(1,ivar,head)
            fhat(2,ivar,head) = dfds(0)

        end do
        
    !-------------------------------------- upper-endpoint !
        
        do  ivar = +1 , nvar-0

            call plsv( &
        &   fdat(1,ivar,tail-1) , &
        &   delx(tail-1), &
        &   fdat(1,ivar,tail+0) , &
        &   delx(tail+0), &
        &   fdat(1,ivar,tail+0) , &
        &   delx(tail+0), dfds)

            fhat(1,ivar,tail) = &
        &   fdat(1,ivar,tail)
            fhat(2,ivar,tail) = dfds(0)

        end do

    !-------------------------------------- interior cells !

        do  ipos = +2 , npos-2
        do  ivar = +1 , nvar-0

            call plsv( &
        &   fdat(1,ivar,ipos-1) , &
        &   delx(ipos-1), &
        &   fdat(1,ivar,ipos+0) , &
        &   delx(ipos+0), &
        &   fdat(1,ivar,ipos+1) , &
        &   delx(ipos+1), dfds)

            fhat(1,ivar,ipos) = &
        &   fdat(1,ivar,ipos)
            fhat(2,ivar,ipos) = dfds(0)

        end do
        end do
        
        return
        
    end  subroutine

    !------------------------- assemble PLM reconstruction !

    pure subroutine plmc(npos,nvar,ndof,delx, &
        &                fdat,fhat,dmin,ilim)

    !
    ! *this is the constant grid-spacing variant .
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell .
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh .
    ! ILIM  cell slope-limiting selection .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nvar
        integer, intent( in) :: ndof,ilim
        real*8 , intent( in) :: dmin
        real*8 , intent( in) :: delx(1)
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)
        
    !------------------------------------------- variables !
        integer :: ipos,ivar,head,tail
        real*8  :: dfds(-1:+1)

        head = +1; tail = npos - 1

        if (npos.eq.2) then
    !----------------------- reduce order if small stencil !
        do  ivar = +1, nvar
            fhat(1,ivar,1) = &
        &   fdat(1,ivar,1)
            fhat(2,ivar,1) = 0.e+0
        end do
        end if

        if (npos.le.2) return
  
    !-------------------------------------- lower-endpoint !
        
        do  ivar = +1 , nvar-0

            call plsc( &
        &   fdat(1,ivar,head+0) , &
        &   fdat(1,ivar,head+0) , &
        &   fdat(1,ivar,head+1) , &
        &        dfds)

            fhat(1,ivar,head) = &
        &   fdat(1,ivar,head)
            fhat(2,ivar,head) = dfds(0)

        end do
        
    !-------------------------------------- upper-endpoint !
        
        do  ivar = +1 , nvar-0

            call plsc( &
        &   fdat(1,ivar,tail-1) , &
        &   fdat(1,ivar,tail+0) , &
        &   fdat(1,ivar,tail+0) , &
        &        dfds)

            fhat(1,ivar,tail) = &
        &   fdat(1,ivar,tail)
            fhat(2,ivar,tail) = dfds(0)

        end do

    !-------------------------------------- interior cells !

        do  ipos = +2 , npos-2
        do  ivar = +1 , nvar-0

            call plsc( &
        &   fdat(1,ivar,ipos-1) , &
        &   fdat(1,ivar,ipos+0) , &
        &   fdat(1,ivar,ipos+1) , &
        &        dfds)

            fhat(1,ivar,ipos) = &
        &   fdat(1,ivar,ipos)
            fhat(2,ivar,ipos) = dfds(0)

        end do
        end do
        
        return
        
    end  subroutine 
    
    !------------------------------- assemble PLM "slopes" !
    
    pure subroutine plsv(ffll,hhll,ff00,hh00,&
        &                ffrr,hhrr,dfds)

    !
    ! *this is the variable grid-spacing variant .
    !
    ! FFLL  left -biased grid-cell mean.
    ! HHLL  left -biased grid-cell spac.
    ! FF00  centred grid-cell mean.
    ! HH00  centred grid-cell spac.
    ! FFRR  right-biased grid-cell mean.
    ! HHRR  right-biased grid-cell spac.
    ! DFDS  piecewise linear gradients in local co-ord.'s.
    !       DFDS(+0) is a centred, slope-limited estimate,
    !       DFDS(-1), DFDS(+1) are left- and right-biased
    !       estimates (unlimited).
    !

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent( in) :: ffll,ff00,ffrr
        real*8 , intent( in) :: hhll,hh00,hhrr
        real*8 , intent(out) :: dfds(-1:+1)

    !------------------------------------------- variables !
        real*8  :: fell,ferr,scal

        real*8 , parameter :: ZERO = 1.e-14

    !---------------------------- 2nd-order approximations !

        dfds(-1) = ff00-ffll
        dfds(+1) = ffrr-ff00

        if (dfds(-1) * &
        &   dfds(+1) .gt. 0.0d+0) then

    !---------------------------- calc. ll//rr edge values !

            fell = (hh00*ffll+hhll*ff00) & 
        &        / (hhll+hh00)        
            ferr = (hhrr*ff00+hh00*ffrr) &
        &        / (hh00+hhrr)

    !---------------------------- calc. centred derivative !
            
            dfds(+0) = &
        &       0.5d+0 * (ferr - fell)

    !---------------------------- monotonic slope-limiting !
            
            scal = min(abs(dfds(-1)), &
        &              abs(dfds(+1))) &
        &        / max(abs(dfds(+0)), &
                       ZERO)
            scal = min(scal,+1.0d+0)

            dfds(+0) = scal * dfds(+0)

        else

    !---------------------------- flatten if local extrema ! 
      
            dfds(+0) =      +0.0d+0
        
        end if
        
    !---------------------------- scale onto local co-ord. !
        
        dfds(-1) = dfds(-1) &
        &       / (hhll + hh00) * hh00
        dfds(+1) = dfds(+1) &
        &       / (hh00 + hhrr) * hh00

        return

    end  subroutine
    
    !------------------------------- assemble PLM "slopes" !
    
    pure subroutine plsc(ffll,ff00,ffrr,dfds)

    !
    ! *this is the constant grid-spacing variant .
    !
    ! FFLL  left -biased grid-cell mean.
    ! FF00  centred grid-cell mean.
    ! FFRR  right-biased grid-cell mean.
    ! DFDS  piecewise linear gradients in local co-ord.'s.
    !       DFDS(+0) is a centred, slope-limited estimate,
    !       DFDS(-1), DFDS(+1) are left- and right-biased
    !       estimates (unlimited).
    !

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent( in) :: ffll,ff00,ffrr
        real*8 , intent(out) :: dfds(-1:+1)

    !------------------------------------------- variables !
        real*8  :: fell,ferr,scal

        real*8 , parameter :: ZERO = 1.e-14

    !---------------------------- 2nd-order approximations !

        dfds(-1) = ff00-ffll
        dfds(+1) = ffrr-ff00

        if (dfds(-1) * &
        &   dfds(+1) .gt. 0.0d+0) then

    !---------------------------- calc. ll//rr edge values !

            fell = (ffll+ff00) * .5d+0        
            ferr = (ff00+ffrr) * .5d+0

    !---------------------------- calc. centred derivative !
            
            dfds(+0) = &
        &       0.5d+0 * (ferr - fell)

    !---------------------------- monotonic slope-limiting !
            
            scal = min(abs(dfds(-1)), &
        &              abs(dfds(+1))) &
        &        / max(abs(dfds(+0)), &
                       ZERO)
            scal = min(scal,+1.0d+0)

            dfds(+0) = scal * dfds(+0)

        else

    !---------------------------- flatten if local extrema ! 
      
            dfds(+0) =      +0.0d+0
        
        end if
        
    !---------------------------- scale onto local co-ord. !
        
        dfds(-1) = + 0.5d+0 * dfds(-1) 
        dfds(+1) = + 0.5d+0 * dfds(+1)

        return

    end  subroutine   
    
    
    

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! PPM.h90: 1d slope-limited, piecewise parabolic recon.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !
  
    ! P. Colella and PR. Woodward, The Piecewise Parabolic 
    ! Method (PPM) for gas-dynamical simulations., J. Comp. 
    ! Phys., 54 (1), 1984, 174-201,
    ! https://doi.org/10.1016/0021-9991(84)90143-8
    !

    pure subroutine ppm(npos,nvar,ndof,delx, &
        &               fdat,fhat,edge,oscl, &
        &               dmin,ilim,wlim,halo)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh . 
    ! ILIM  cell slope-limiting selection .
    ! WLIM  wall slope-limiting selection .
    ! HALO  width of re-con. stencil, symmetric about mid. .
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        real*8 , intent(in)  :: dmin
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: fdat(:,:,:)
        real*8 , intent(in)  :: edge(:,:)
        integer, intent(in)  :: ilim,wlim,halo

    !------------------------------------------- variables !
        integer :: ipos,ivar,iill,iirr,head,tail
        real*8  :: ff00,ffll,ffrr,hh00,hhll,hhrr
        integer :: mono
        real*8  :: fell,ferr
        real*8  :: dfds(-1:+1)
        real*8  :: wval(+1:+2)
        real*8  :: uhat(+1:+3)
        real*8  :: lhat(+1:+3)
        
        head = +1; tail = npos - 1

        if (npos.eq.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = +1, nvar
            fhat(1,ivar,+1) = &
        &   fdat(1,ivar,+1)
            fhat(2,ivar,+1) = 0.e0_8
            fhat(3,ivar,+1) = 0.e0_8
        end do
        end if

        if (npos.le.2) return

    !------------------- reconstruct function on each cell !

        uhat = +0.e+0
        lhat = +0.e+0

        do  ipos = +1 , npos-1

        iill = max(head,ipos-1)
        iirr = min(tail,ipos+1)

        do  ivar = +1 , nvar-0

    !----------------------------- cell mean + edge values !
    
            ff00 = fdat(1,ivar,ipos)
            ffll = fdat(1,ivar,iill)
            ffrr = fdat(1,ivar,iirr)

            fell = edge(ivar,ipos+0)
            ferr = edge(ivar,ipos+1)

    !----------------------------- calc. LL/00/RR gradient !
    
            if (size(delx).gt.+1) then

            hh00 = delx(ipos)
            hhll = delx(iill)
            hhrr = delx(iirr)

            call plsv (ffll,hhll,ff00, &
    &                  hh00,ffrr,hhrr, &
    &                  dfds)
            else
            
            call plsc (ffll,ff00,ffrr, &
    &                  dfds)
    
            end if

    !----------------------------- calc. cell-wise profile !
    
            select case(ilim)
            case (null_limit)

    !----------------------------- calc. unlimited profile !   
               
            call ppmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dfds, &
    &                  uhat,lhat,mono)
            
    !----------------------------- pref. unlimited profile !

            wval(1) = +1.e+0
            wval(2) = +0.e+0
            
            case (mono_limit)
            
    !----------------------------- calc. monotonic profile ! 
             
            call ppmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dfds, &
    &                  uhat,lhat,mono)
            
    !----------------------------- pref. monotonic profile !

            wval(1) = +0.e+0
            wval(2) = +1.e+0
            
            case (weno_limit)

    !----------------------------- calc. unlimited profile !  
            
            call ppmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dfds, &
    &                  uhat,lhat,mono)

            if (mono.gt.+0) then
  
    !----------------------------- calc. WENO-type weights ! 
     
            call wenoi(npos,delx,oscl, &
    &                  ipos,ivar,halo, &
    &                  wlim,wval)
            
            else
                
    !----------------------------- pref. unlimited profile !

            wval(1) = +1.e+0
            wval(2) = +0.e+0
            
            end if
            
            end select
 
    !----------------------------- blend "null" and "mono" !
               
            fhat(1,ivar,ipos) = &
    &       wval(1) * uhat(1) + &
    &       wval(2) * lhat(1)
            fhat(2,ivar,ipos) = &
    &       wval(1) * uhat(2) + &
    &       wval(2) * lhat(2)
            fhat(3,ivar,ipos) = &
    &       wval(1) * uhat(3) + &
    &       wval(2) * lhat(3)

        end do

        end do
        
        return

    end  subroutine

    !--------- assemble piecewise parabolic reconstruction !
    
    pure subroutine ppmfn(ff00,ffll,ffrr,fell,&
        &                 ferr,dfds,uhat,lhat,&
        &                 mono)

    !
    ! FF00  centred grid-cell mean.
    ! FFLL  left -biased grid-cell mean.
    ! FFRR  right-biased grid-cell mean.
    ! FELL  left -biased edge interp.
    ! FERR  right-biased edge interp.
    ! DFDS  piecewise linear gradients in local co-ord.'s.
    !       DFDS(+0) is a centred, slope-limited estimate,
    !       DFDS(-1), DFDS(+1) are left- and right-biased
    !       estimates (unlimited). 
    ! UHAT  unlimited PPM reconstruction coefficients .
    ! LHAT  monotonic PPM reconstruction coefficients .
    ! MONO  slope-limiting indicator, MONO > +0 if some 
    !       limiting has occured .
    !

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent(in)    :: ff00
        real*8 , intent(in)    :: ffll,ffrr
        real*8 , intent(inout) :: fell,ferr
        real*8 , intent(in)    :: dfds(-1:+1)
        real*8 , intent(out)   :: uhat(+1:+3)
        real*8 , intent(out)   :: lhat(+1:+3)
        integer, intent(out)   :: mono
          
    !------------------------------------------- variables !
        real*8  :: turn
        
        mono  = 0
        
    !-------------------------------- "null" slope-limiter !

        uhat( 1 ) = &
    & + (3.0d+0 / 2.0d+0) * ff00 &
    & - (1.0d+0 / 4.0d+0) *(ferr+fell)
        uhat( 2 ) = &
    & + (1.0d+0 / 2.0d+0) *(ferr-fell)
        uhat( 3 ) = &
    & - (3.0d+0 / 2.0d+0) * ff00 &
    & + (3.0d+0 / 4.0d+0) *(ferr+fell)
    
    !-------------------------------- "mono" slope-limiter !
        
        if((ffrr - ff00) * & 
    &      (ff00 - ffll) .lt. 0.e+0) then

    !----------------------------------- "flatten" extrema !
    
            mono = +1

            lhat(1) = ff00
            lhat(2) = 0.e0_8
            lhat(3) = 0.e0_8
             
            return
              
        end if

    !----------------------------------- limit edge values !
    
        if((ffll - fell) * &
    &      (fell - ff00) .le. 0.e+0) then

            mono = +1
       
            fell = ff00 - dfds(0)

        end if

        if((ffrr - ferr) * &
    &      (ferr - ff00) .le. 0.e+0) then

            mono = +1

            ferr = ff00 + dfds(0)
         
        end if
   
    !----------------------------------- update ppm coeff. !

        lhat( 1 ) = &
    & + (3.0d+0 / 2.0d+0) * ff00 &
    & - (1.0d+0 / 4.0d+0) *(ferr+fell)
        lhat( 2 ) = &
    & + (1.0d+0 / 2.0d+0) *(ferr-fell)
        lhat( 3 ) = &
    & - (3.0d+0 / 2.0d+0) * ff00 &
    & + (3.0d+0 / 4.0d+0) *(ferr+fell)

    !----------------------------------- limit cell values !
    
        if (abs(lhat(3)) .gt. & 
    &       abs(lhat(2))*.5d+0) then

        turn = -0.5d+0 * lhat(2) &
    &                  / lhat(3)

        if ((turn .ge. -1.e+0)&
    &  .and.(turn .le. +0.e+0)) then

        mono =   +2

    !--------------------------- push TURN onto lower edge !

        ferr =   +3.0d+0  * ff00 &
    &            -2.0d+0  * fell

        lhat( 1 ) = &
    & + (3.0d+0 / 2.0d+0) * ff00 &
    & - (1.0d+0 / 4.0d+0) *(ferr+fell)
        lhat( 2 ) = &
    & + (1.0d+0 / 2.0d+0) *(ferr-fell)
        lhat( 3 ) = &
    & - (3.0d+0 / 2.0d+0) * ff00 &
    & + (3.0d+0 / 4.0d+0) *(ferr+fell)

        else &
    &   if ((turn .gt. +0.e+0)&
    &  .and.(turn .le. +1.e+0)) then

        mono =   +2

    !--------------------------- push TURN onto upper edge !
    
        fell =   +3.0d+0  * ff00 &
    &            -2.0d+0  * ferr

        lhat( 1 ) = &
    & + (3.0d+0 / 2.0d+0) * ff00 &
    & - (1.0d+0 / 4.0d+0) *(ferr+fell)
        lhat( 2 ) = &
    & + (1.0d+0 / 2.0d+0) *(ferr-fell)
        lhat( 3 ) = &
    & - (3.0d+0 / 2.0d+0) * ff00 &
    & + (3.0d+0 / 4.0d+0) *(ferr+fell)
   
        end if
      
        end if
   
        return
    
    end subroutine




    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! PQM.h90: a 1d slope-limited, piecewise quartic recon.
    !
    ! Darren Engwirda 
    ! 08-Sep-2016
    ! de2363 [at] columbia [dot] edu
    !
    !

    ! White, L. and Adcroft, A., A high-order finite volume
    ! remapping scheme for nonuniform grids: The piecewise 
    ! quartic method (PQM), J. Comp. Phys., 227 (15), 2008, 
    ! 7394-7422, https://doi.org/10.1016/j.jcp.2008.04.026.
    !

    pure subroutine pqm(npos,nvar,ndof,delx, &
        &               fdat,fhat,edge,dfdx, &
        &               oscl,dmin,ilim,wlim, &
        &               halo)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! DELX  grid-cell spacing array. LENGTH(DELX) == +1 if 
    !       spacing is uniform .
    ! FDAT  grid-cell moments array. FDAT is an array with
    !       SIZE = NDOF-by-NVAR-by-NPOS-1 .
    ! FHAT  grid-cell re-con. array. FHAT is an array with
    !       SIZE = MDOF-by-NVAR-by-NPOS-1 .
    ! EDGE  edge-centred interp. for function-value. EDGE
    !       is an array with SIZE = NVAR-by-NPOS .
    ! DFDX  edge-centred interp. for 1st-derivative. DFDX
    !       is an array with SIZE = NVAR-by-NPOS .
    ! OSCL  grid-cell oscil. dof.'s. OSCL is an array with
    !       SIZE = +2  -by-NVAR-by-NPOS-1 .
    ! DMIN  min. grid-cell spacing thresh . 
    ! ILIM  cell slope-limiting selection .
    ! WLIM  wall slope-limiting selection .
    ! HALO  width of re-con. stencil, symmetric about mid. .
    !

        implicit none
        
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        integer, intent(in)  :: ilim,wlim,halo
        real*8 , intent(in)  :: dmin
        real*8 , intent(out) :: fhat(:,:,:)
        real*8 , intent(in)  :: oscl(:,:,:)
        real*8 , intent(in)  :: delx(:)
        real*8 , intent(in)  :: fdat(:,:,:)
        real*8 , intent(in)  :: edge(:,:)
        real*8 , intent(in)  :: dfdx(:,:)

    !------------------------------------------- variables !
        integer :: ipos,ivar,iill,iirr,head,tail
        real*8  :: ff00,ffll,ffrr,hh00,hhll,hhrr
        real*8  :: xhat
        integer :: mono
        real*8  :: fell,ferr
        real*8  :: dell,derr
        real*8  :: dfds(-1:+1)
        real*8  :: uhat(+1:+5)
        real*8  :: lhat(+1:+5)
        real*8  :: wval(+1:+2)
        
        head = +1; tail = npos - 1

        if (npos.le.2) then
    !----- default to reduced order if insufficient points !
        do  ivar = +1, nvar
            fhat(1,ivar,+1) = fdat(1,ivar,+1)
            fhat(2,ivar,+1) = 0.e0_8
            fhat(3,ivar,+1) = 0.e0_8
            fhat(4,ivar,+1) = 0.e0_8
            fhat(5,ivar,+1) = 0.e0_8
        end do
        end if

        if (npos.le.2) return

    !------------------- reconstruct function on each cell !

        do  ipos = +1 , npos-1

        iill = max(head,ipos-1)
        iirr = min(tail,ipos+1)

        do  ivar = +1 , nvar-0

    !----------------------------- cell mean + edge values !
            
            ff00 = fdat(1,ivar,ipos)
            ffll = fdat(1,ivar,iill)
            ffrr = fdat(1,ivar,iirr)
           
            fell = edge(ivar,ipos+0)
            ferr = edge(ivar,ipos+1)
        
    !----------------------------- calc. LL/00/RR gradient !
            
            if (size(delx).gt.+1) then

            hh00 = delx(ipos)
            hhll = delx(iill)
            hhrr = delx(iirr)

            xhat = delx(ipos+0)*.5d+0

            call plsv (ffll,hhll,ff00, &
    &                  hh00,ffrr,hhrr, &
    &                  dfds)
            else
            
            xhat = delx(    +1)*.5d+0
            
            call plsc (ffll,ff00,ffrr, &
    &                  dfds)

            end if              

            dell = dfdx (ivar,ipos+0)
            dell = dell * xhat
            
            derr = dfdx (ivar,ipos+1)
            derr = derr * xhat

    !----------------------------- calc. cell-wise profile !
            
            select case(ilim)
            case (null_limit)

    !----------------------------- calc. unlimited profile !              
            
            call pqmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dell, &
    &                  derr,dfds,uhat, &
    &                  lhat,mono)
        
    !----------------------------- pref. unlimited profile !
            
            wval(1) = +1.e+0
            wval(2) = +0.e+0
            
            case (mono_limit)
            
    !----------------------------- calc. monotonic profile !              
            
            call pqmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dell, &
    &                  derr,dfds,uhat, &
    &                  lhat,mono)
            
    !----------------------------- pref. monotonic profile !
            
            wval(1) = +0.e+0
            wval(2) = +1.e+0
            
            case (weno_limit)

    !----------------------------- calc. monotonic profile !              
            
            call pqmfn(ff00,ffll,ffrr, &
    &                  fell,ferr,dell, &
    &                  derr,dfds,uhat, &
    &                  lhat,mono)
            
            if (mono.gt.+0) then
            
    !----------------------------- calc. WENO-type weights !      
            
            call wenoi(npos,delx,oscl, &
    &                  ipos,ivar,halo, &
    &                  wlim,wval)
            
            else
            
    !----------------------------- pref. unlimited profile !
            
            wval(1) = +1.e+0
            wval(2) = +0.e+0
            
            end if
            
            end select

    !----------------------------- blend "null" and "mono" !
                   
            fhat(1,ivar,ipos) = &
    &       wval(1) * uhat(1) + &
    &       wval(2) * lhat(1)
            fhat(2,ivar,ipos) = &
    &       wval(1) * uhat(2) + &
    &       wval(2) * lhat(2)
            fhat(3,ivar,ipos) = &
    &       wval(1) * uhat(3) + &
    &       wval(2) * lhat(3)
            fhat(4,ivar,ipos) = &
    &       wval(1) * uhat(4) + &
    &       wval(2) * lhat(4)
            fhat(5,ivar,ipos) = &
    &       wval(1) * uhat(5) + &
    &       wval(2) * lhat(5)

        end do
        
        end do
        
        return

    end subroutine
    
    !----------- assemble piecewise quartic reconstruction !
    
    pure subroutine pqmfn(ff00,ffll,ffrr,fell, &
        &                 ferr,dell,derr,dfds, &
        &                 uhat,lhat,mono)

    !
    ! FF00  centred grid-cell mean.
    ! FFLL  left -biased grid-cell mean.
    ! FFRR  right-biased grid-cell mean.
    ! FELL  left -biased edge interp.
    ! FERR  right-biased edge interp.
    ! DELL  left -biased edge df//dx.
    ! DERR  right-biased edge df//dx.
    ! DFDS  piecewise linear gradients in local co-ord.'s.
    !       DFDS(+0) is a centred, slope-limited estimate,
    !       DFDS(-1), DFDS(+1) are left- and right-biased
    !       estimates (unlimited). 
    ! UHAT  unlimited PPM reconstruction coefficients .
    ! LHAT  monotonic PPM reconstruction coefficients .
    ! MONO  slope-limiting indicator, MONO > +0 if some 
    !       limiting has occured .
    !

        implicit none

    !------------------------------------------- arguments !
        real*8 , intent(in)    :: ff00
        real*8 , intent(in)    :: ffll,ffrr
        real*8 , intent(inout) :: fell,ferr
        real*8 , intent(inout) :: dell,derr
        real*8 , intent(in)    :: dfds(-1:+1)
        real*8 , intent(out)   :: uhat(+1:+5)
        real*8 , intent(out)   :: lhat(+1:+5)
        integer, intent(out)   :: mono
          
    !------------------------------------------- variables !
        integer :: turn
        real*8  :: grad, iflx(+1:+2)
        logical :: haveroot
        
    !-------------------------------- "null" slope-limiter !
        
        mono    = 0
        
        uhat(1) = &
    & + (30.e+0 / 16.e+0) * ff00 &
    & - ( 7.e+0 / 16.e+0) *(ferr+fell) &
    & + ( 1.e+0 / 16.e+0) *(derr-dell)
        uhat(2) = &
    & + ( 3.e+0 /  4.e+0) *(ferr-fell) &
    & - ( 1.e+0 /  4.e+0) *(derr+dell)
        uhat(3) = &
    & - (30.e+0 /  8.e+0) * ff00 &
    & + (15.e+0 /  8.e+0) *(ferr+fell) &
    & - ( 3.e+0 /  8.e+0) *(derr-dell)
        uhat(4) = &
    & - ( 1.e+0 /  4.e+0) *(ferr-fell  &
    &                      -derr-dell)
        uhat(5) = &
    & + (30.e+0 / 16.e+0) * ff00 &
    & - (15.e+0 / 16.e+0) *(ferr+fell) &
    & + ( 5.e+0 / 16.e+0) *(derr-dell)
  
    !-------------------------------- "mono" slope-limiter !
    
        if((ffrr - ff00) * & 
    &      (ff00 - ffll) .le. 0.e+0) then

    !----------------------------------- "flatten" extrema !
    
            mono = +1

            lhat(1) = ff00
            lhat(2) = 0.e0_8
            lhat(3) = 0.e0_8
            lhat(4) = 0.e0_8
            lhat(5) = 0.e0_8
             
            return
              
        end if

    !----------------------------------- limit edge values !

        if((ffll - fell) * &
    &      (fell - ff00) .le. 0.e+0) then

            mono = +1
       
            fell = ff00 - dfds(0)

        end if

        if (dell * dfds(0) .lt. 0.e+0) then

            mono = +1

            dell = dfds(0)
 
        end if

        if((ffrr - ferr) * &
    &      (ferr - ff00) .le. 0.e+0) then

            mono = +1

            ferr = ff00 + dfds(0)
         
        end if

        if (derr * dfds(0) .lt. 0.e+0) then

            mono = +1

            derr = dfds(0)

        end if
    
    !----------------------------------- limit cell values !
    
        lhat(1) = &
    & + (30.e+0 / 16.e+0) * ff00 &
    & - ( 7.e+0 / 16.e+0) *(ferr+fell) &
    & + ( 1.e+0 / 16.e+0) *(derr-dell)
        lhat(2) = &
    & + ( 3.e+0 /  4.e+0) *(ferr-fell) &
    & - ( 1.e+0 /  4.e+0) *(derr+dell)
        lhat(3) = &
    & - (30.e+0 /  8.e+0) * ff00 &
    & + (15.e+0 /  8.e+0) *(ferr+fell) &
    & - ( 3.e+0 /  8.e+0) *(derr-dell)
        lhat(4) = &
    & - ( 1.e+0 /  4.e+0) *(ferr-fell  &
    &                      -derr-dell)
        lhat(5) = &
    & + (30.e+0 / 16.e+0) * ff00 &
    & - (15.e+0 / 16.e+0) *(ferr+fell) &
    & + ( 5.e+0 / 16.e+0) *(derr-dell)

    !------------------ calc. inflexion via 2nd-derivative !
        
        call roots_2(12.e+0 * lhat(5), &
    &                 6.e+0 * lhat(4), &
    &                 2.e+0 * lhat(3), &
    &                 iflx  , haveroot )

        if (haveroot) then

        turn = +0

        if ( ( iflx(1) .gt. -1.e+0 ) &
    &  .and. ( iflx(1) .lt. +1.e+0 ) ) then

    !------------------ check for non-monotonic inflection !

        grad = lhat(2) &
    &+ (iflx(1)**1) * 2.e+0* lhat(3) &
    &+ (iflx(1)**2) * 3.e+0* lhat(4) &
    &+ (iflx(1)**3) * 4.e+0* lhat(5)

        if (grad * dfds(0) .lt. 0.e+0) then

            if (abs(dfds(-1)) &
    &      .lt. abs(dfds(+1)) ) then

                turn = -1   ! modify L

            else

                turn = +1   ! modify R

            end if

        end if

        end if
            
        if ( ( iflx(2) .gt. -1.e+0 ) &
    &  .and. ( iflx(2) .lt. +1.e+0 ) ) then

    !------------------ check for non-monotonic inflection !
                
        grad = lhat(2) &
    &+ (iflx(2)**1) * 2.e+0* lhat(3) &
    &+ (iflx(2)**2) * 3.e+0* lhat(4) &
    &+ (iflx(2)**3) * 4.e+0* lhat(5)

        if (grad * dfds(0) .lt. 0.e+0) then

            if (abs(dfds(-1)) &
    &      .lt. abs(dfds(+1)) ) then

                turn = -1   ! modify L

            else

                turn = +1   ! modify R

            end if

        end if

        end if
            
    !------------------ pop non-monotone inflexion to edge !
              
        if (turn .eq. -1) then

    !------------------ pop inflection points onto -1 edge !

        mono = +2

        derr = &
    &- ( 5.e+0 /  1.e+0) * ff00 &
    &+ ( 3.e+0 /  1.e+0) * ferr &
    &+ ( 2.e+0 /  1.e+0) * fell
        dell = &
    &+ ( 5.e+0 /  3.e+0) * ff00 &
    &- ( 1.e+0 /  3.e+0) * ferr &
    &- ( 4.e+0 /  3.e+0) * fell

        if (dell*dfds(+0) .lt. 0.e+0) then

        dell =   0.e+0

        ferr = &
    &+ ( 5.e+0 /  1.e+0) * ff00 &
    &- ( 4.e+0 /  1.e+0) * fell
        derr = &
    &+ (10.e+0 /  1.e+0) * ff00 &
    &- (10.e+0 /  1.e+0) * fell

        else &
    &   if (derr*dfds(+0) .lt. 0.e+0) then

        derr =   0.e+0

        fell = &
    &+ ( 5.e+0 /  2.e+0) * ff00 &
    &- ( 3.e+0 /  2.e+0) * ferr
        dell = &
    &- ( 5.e+0 /  3.e+0) * ff00 &
    &+ ( 5.e+0 /  3.e+0) * ferr

        end if

        lhat(1) = &
    &+ (30.e+0 / 16.e+0) * ff00 &
    &- ( 7.e+0 / 16.e+0) *(ferr+fell) &
    &+ ( 1.e+0 / 16.e+0) *(derr-dell)
        lhat(2) = &
    &+ ( 3.e+0 /  4.e+0) *(ferr-fell) &
    &- ( 1.e+0 /  4.e+0) *(derr+dell)
        lhat(3) = &
    &- (30.e+0 /  8.e+0) * ff00 &
    &+ (15.e+0 /  8.e+0) *(ferr+fell) &
    &- ( 3.e+0 /  8.e+0) *(derr-dell)
        lhat(4) = &
    &- ( 1.e+0 /  4.e+0) *(ferr-fell  &
    &                     -derr-dell)
        lhat(5) = &
    &+ (30.e+0 / 16.e+0) * ff00 &
    &- (15.e+0 / 16.e+0) *(ferr+fell) &
    &+ ( 5.e+0 / 16.e+0) *(derr-dell)

        end if

        if (turn .eq. +1) then

    !------------------ pop inflection points onto -1 edge !
    
        mono = +2

        derr = &
    &- ( 5.e+0 /  3.e+0) * ff00 &
    &+ ( 4.e+0 /  3.e+0) * ferr &
    &+ ( 1.e+0 /  3.e+0) * fell
        dell = &
    &+ ( 5.e+0 /  1.e+0) * ff00 &
    &- ( 2.e+0 /  1.e+0) * ferr &
    &- ( 3.e+0 /  1.e+0) * fell

        if (dell*dfds(+0) .lt. 0.e+0) then

        dell =   0.e+0

        ferr = &
    &+ ( 5.e+0 /  2.e+0) * ff00 &
    &- ( 3.e+0 /  2.e+0) * fell
        derr = &
    &+ ( 5.e+0 /  3.e+0) * ff00 &
    &- ( 5.e+0 /  3.e+0) * fell

        else & 
    &   if (derr*dfds(+0) .lt. 0.e+0) then

        derr =   0.e+0

        fell = &
    &+ ( 5.e+0 /  1.e+0) * ff00 &
    &- ( 4.e+0 /  1.e+0) * ferr
        dell = &
    &- (10.e+0 /  1.e+0) * ff00 &
    &+ (10.e+0 /  1.e+0) * ferr

        end if

        lhat(1) = &
    &+ (30.e+0 / 16.e+0) * ff00 &
    &- ( 7.e+0 / 16.e+0) *(ferr+fell) &
    &+ ( 1.e+0 / 16.e+0) *(derr-dell)
        lhat(2) = &
    &+ ( 3.e+0 /  4.e+0) *(ferr-fell) &
    &- ( 1.e+0 /  4.e+0) *(derr+dell)
        lhat(3) = &
    &- (30.e+0 /  8.e+0) * ff00 &
    &+ (15.e+0 /  8.e+0) *(ferr+fell) &
    &- ( 3.e+0 /  8.e+0) *(derr-dell)
        lhat(4) = &
    &- ( 1.e+0 /  4.e+0) *(ferr-fell  &
    &                     -derr-dell)
        lhat(5) = &
    &+ (30.e+0 / 16.e+0) * ff00 &
    &- (15.e+0 / 16.e+0) *(ferr+fell) &
    &+ ( 5.e+0 / 16.e+0) *(derr-dell)

        end if
    
        end if        ! haveroot

        return

    end subroutine
    

    
    
    !------------------------------------------------------!
    ! RMAP1D : one-dimensional conservative "re-map" .     !
    !------------------------------------------------------!
       

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! RMAP1D.h90: high-order integral re-mapping operators.
    !
    ! Darren Engwirda 
    ! 31-Mar-2019
    ! ​de2363 [at] columbia [dot] edu
    !
    !

    subroutine rmap1d(npos,nnew,nvar,ndof,xpos, &
        &             xnew,fdat,fnew,bclo,bcup, &
        &             work,opts,tCPU)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FDAT  grid-cell moments on old grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! FNEW  grid-cell moments on new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! BCLO  boundary condition at lower endpoint .
    ! BCHI  boundary condition at upper endpoint . 
    ! WORK  method work-space. See RCON-WORK for details .
    ! OPTS  method parameters. See RCON-OPTS for details .
    ! TCPU  method tcpu-timer.
    !

        implicit none

    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar,ndof
        class(rmap_work), intent(inout):: work
        class(rmap_opts), intent(inout):: opts
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fdat(:,:,:)
        real*8 , intent(out) :: fnew(:,:,:)
        type (rcon_ends), intent(in) :: bclo(:)
        type (rcon_ends), intent(in) :: bcup(:)
        type (rmap_tics), &
        &   intent(inout) , optional :: tCPU

        real*8 , parameter :: RTOL = +1.e-14

    !------------------------------------------- variables !
        integer :: ipos
        real*8  :: diff,spac,same,xtol
        real*8  :: delx(1)
        logical :: uniform
        

        if (ndof.lt.1) return
        if (npos.lt.2) return
        if (nnew.lt.2) return
        if (nvar.lt.1) return

    !------------- calc. grid-spacing and check uniformity !

        same = (xpos(npos)& 
             -  xpos(  +1)) / (npos-1)

        uniform = .true.
             
        xtol = same * RTOL

        do  ipos = +1 , npos-1, +1
        
            spac = xpos(ipos+1) &
        &        - xpos(ipos+0)
  
            diff = abs(spac - same)
        
            if (diff.gt.xtol) then
            
                uniform = .false.
            
            end if
  
            work% &
        &   cell_spac(ipos) = spac
        
        end do

       !uniform = .false.
    
    !------------- reconstruct FHAT over all cells in XPOS !

        if (.not. uniform) then

    !------------------------------------ variable spacing !
        call rcon1d(npos,nvar,ndof, &
        &           work%cell_spac, &
        &           fdat,bclo,bcup, &
        &           work%cell_func, &
        &           work,opts,tCPU)

        else
        
    !------------------------------------ constant spacing !
        delx(1) = work%cell_spac(1)
        
        call rcon1d(npos,nvar,ndof, &
        &           delx,    &
        &           fdat,bclo,bcup, &
        &           work%cell_func, &
        &           work,opts,tCPU)
        
        end if

    !------------- remap FDAT from XPOS to XNEW using FHAT !


        select case(opts%cell_meth)
        case(pcm_method)
    !------------------------------------ 1st-order method !
        call imap1d(npos,nnew,nvar, &
        &           ndof,  +1,      &
        &           xpos,xnew,      &
        &           work%cell_func, &
        &           fdat,fnew,xtol)
    
        case(plm_method)
    !------------------------------------ 2nd-order method !
        call imap1d(npos,nnew,nvar, &
        &           ndof,  +2,      &
        &           xpos,xnew,      &
        &           work%cell_func, &
        &           fdat,fnew,xtol)

        case(ppm_method)
    !------------------------------------ 3rd-order method !
        call imap1d(npos,nnew,nvar, &
        &           ndof,  +3,      &
        &           xpos,xnew,      &
        &           work%cell_func, &
        &           fdat,fnew,xtol)

        case(pqm_method)
    !------------------------------------ 5th-order method !
        call imap1d(npos,nnew,nvar, &
        &           ndof,  +5,      &
        &           xpos,xnew,      &
        &           work%cell_func, &
        &           fdat,fnew,xtol)
        
        end select
        
        
        return
    
    end  subroutine
    
    !------------ IMAP1D: 1-dimensional degree-k remapping !
    
    pure subroutine imap1d(npos,nnew,nvar,ndof, &
        &                  mdof,xpos,xnew,fhat, &
        &                  fdat,fnew,XTOL)

    !
    ! NPOS  no. edges in old grid.
    ! NNEW  no. edges in new grid.
    ! NVAR  no. discrete variables to remap.
    ! NDOF  no. degrees-of-freedom per cell.
    ! MDOF  no. degrees-of-freedom per FHAT.    
    ! XPOS  old grid edge positions. XPOS is a length NPOS
    !       array .
    ! XNEW  new grid edge positions. XNEW is a length NNEW
    !       array .
    ! FHAT  reconstruction over old grid. FHAT has SIZE =
    !       MDOF-by-NVAR-by-NPOS-1 .
    ! FDAT  grid-cell moments on old grid. FDAT has SIZE = 
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! FNEW  grid-cell moments on new grid. FNEW has SIZE = 
    !       NDOF-by-NVAR-by-NNEW-1 .
    ! XTOL  min. grid-cell thickness . 
    !

        implicit none    
    
    !------------------------------------------- arguments !
        integer, intent( in) :: npos,nnew
        integer, intent( in) :: nvar
        integer, intent( in) :: ndof,mdof
        real*8 , intent( in) :: xpos(:)
        real*8 , intent( in) :: xnew(:)
        real*8 , intent( in) :: fhat(:,:,:)
        real*8 , intent( in) :: fdat(:,:,:)        
        real*8 , intent(out) :: fnew(:,:,:)
        real*8 , intent( in) :: XTOL
        
    !------------------------------------------- variables !
        integer :: kpos,ipos,ivar,idof
        integer :: pos0,pos1,vmin,vmax
        real*8  :: xmid,xhat,khat,stmp
        real*8  :: xxlo,xxhi,sslo,sshi,intf
        real*8  :: vvlo(  +1:+5)
        real*8  :: vvhi(  +1:+5)        
        real*8  :: ivec(  +1:+5)
        real*8  :: sdat(  +1:nvar)
        real*8  :: snew(  +1:nvar)
        real*8  :: serr(  +1:nvar)
        integer :: kmin(  +1:nvar)
        integer :: kmax(  +1:nvar)
        
        integer, parameter :: INTB = -1   ! integral basis  

    !------------------------------------- initializations !

        vvlo(+1:+5) = 0.e0_8
        vvhi(+1:+5) = 0.e0_8

    !------------- remap FDAT from XPOS to XNEW using FHAT !

        kmin = +1 ; kmax = +1
        pos0 = +1 ; pos1 = +1

        do  kpos = +1, nnew-1

    !------ first cell in XPOS overlapping with XNEW(KPOS) !

            pos1 = max(pos1,1)

            do  pos0 = pos1, npos-1

                if (xpos(pos0+1)&
            &       .gt. xnew(kpos+0)) exit

            end do

    !------ final cell in XPOS overlapping with XNEW(KPOS) !
    
            do  pos1 = pos0, npos-1

                if (xpos(pos1+0)&
            &       .ge. xnew(kpos+1)) exit
    
            end do
            
            pos1 = pos1 - 1

    !------------- integrate FHAT across overlapping cells !

            khat = xnew(kpos+1) &
            &    - xnew(kpos+0)
            khat = max (khat , XTOL)

            do  idof = +1,ndof
            do  ivar = +1,nvar
                
                fnew(idof,ivar,kpos) = 0.e0_8
            
            end do
            end do

            do  ipos = pos0, pos1

    !------------------------------- integration endpoints !
    
                xxlo = max (xpos(ipos+0) , &
            &               xnew(kpos+0))
                xxhi = min (xpos(ipos+1) , &
            &               xnew(kpos+1))

    !------------------------------- local endpoint coords !
    
                xmid = xpos(ipos+1) * .5d0 &
            &        + xpos(ipos+0) * .5d0    
                xhat = xpos(ipos+1) * .5d0 &
            &        - xpos(ipos+0) * .5d0
     
                sslo = &
            &  (xxlo-xmid) / max(xhat,XTOL)
                sshi = &
            &  (xxhi-xmid) / max(xhat,XTOL)

    !------------------------------- integral basis vector !
    
                call bfun1d(INTB,mdof, &
                            sslo,vvlo)
                call bfun1d(INTB,mdof, &
                            sshi,vvhi)
                
                ivec =  vvhi - vvlo

    !--------- integrate FHAT across the overlap XXLO:XXHI !
    
                do  ivar = +1, nvar
    
                intf =  dot_product (  &
            &       ivec(+1:mdof),  &
            &   fhat(+1:mdof,ivar,ipos-0) )

                intf =  intf * xhat
        
    !--------- accumulate integral contributions from IPOS !
    
                fnew(  +1,ivar,kpos) = &
            &   fnew(  +1,ivar,kpos) + intf

                end do

            end do

    !------------------------------- finalise KPOS profile !

            do  ivar = +1, nvar

                fnew(  +1,ivar,kpos) = &
            &   fnew(  +1,ivar,kpos) / khat

    !--------- keep track of MIN/MAX for defect correction !

                vmax =    kmax(ivar)
                vmin =    kmin(ivar)

                if(fnew(1,ivar,kpos) &
            &  .gt.fnew(1,ivar,vmax) ) then
                
                kmax(ivar) =   kpos
            
                else &
            &   if(fnew(1,ivar,kpos) &
            &  .lt.fnew(1,ivar,vmin) ) then
                
                kmin(ivar) =   kpos
            
                end if

            end do

        end do

    !--------- defect corrections: Kahan/Babuska/Neumaier. !

    !   Carefully compute column sums, leading to a defect
    !   wrt. column-wise conservation. Use KBN approach to
    !   account for FP roundoff.

        sdat = 0.e0_8; serr = 0.e0_8
        do  ipos = +1, npos-1
        do  ivar = +1, nvar-0
        
    !------------------------------- integrate old profile !

            xhat = xpos(ipos+1) &
        &        - xpos(ipos+0)

            intf = xhat*fdat(1,ivar,ipos)
            
            stmp = sdat(ivar) + intf

            if (abs(sdat(ivar)) &
        &           .ge. abs(intf)) then

            serr(ivar) = &
        &   serr(ivar) + ((sdat(ivar)-stmp)+intf)

            else

            serr(ivar) = &
        &   serr(ivar) + ((intf-stmp)+sdat(ivar))

            end if

            sdat(ivar) = stmp
        
        end do
        end do

        sdat =  sdat + serr

        snew = 0.e0_8; serr = 0.e0_8
        do  ipos = +1, nnew-1
        do  ivar = +1, nvar-0
        
    !------------------------------- integrate new profile !

            khat = xnew(ipos+1) &
        &        - xnew(ipos+0)

            intf = khat*fnew(1,ivar,ipos)
            
            stmp = snew(ivar) + intf

            if (abs(snew(ivar)) &
        &           .ge. abs(intf)) then

            serr(ivar) = &
        &   serr(ivar) + ((snew(ivar)-stmp)+intf)

            else

            serr(ivar) = &
        &   serr(ivar) + ((intf-stmp)+snew(ivar))

            end if

            snew(ivar) = stmp
        
        end do
        end do

        snew =  snew + serr
        serr =  sdat - snew

    !--------- defect corrections: nudge away from extrema !

    !   Add a correction to remapped state to impose exact
    !   conservation. Via sign(correction), nudge min/max.
    !   cell means, such that monotonicity is not violated 
    !   near extrema...
    
        do  ivar = +1, nvar-0

            if (serr(ivar) .gt. 0.e0_8) then

                vmin = kmin(ivar)

                fnew(1,ivar,vmin) = &
        &       fnew(1,ivar,vmin) + &
        &  serr(ivar)/(xnew(vmin+1)-xnew(vmin+0))

            else &
        &   if (serr(ivar) .lt. 0.e0_8) then

                vmax = kmax(ivar)

                fnew(1,ivar,vmax) = &
        &       fnew(1,ivar,vmax) + &
        &  serr(ivar)/(xnew(vmax+1)-xnew(vmax+0))

            end if

        end do
       
    !------------------------------- new profile now final !

        return
    
    end  subroutine




    !------------------------------------------------------!
    ! FFSL1D : one-dimensional FFSL scalar transport .     !
    !------------------------------------------------------!
    

    !
    ! This program may be freely redistributed under the 
    ! condition that the copyright notices (including this 
    ! entire header) are not removed, and no compensation 
    ! is received through use of the software.  Private, 
    ! research, and institutional use is free.  You may 
    ! distribute modified versions of this code UNDER THE 
    ! CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
    ! TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
    ! ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
    ! MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
    ! NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
    ! of this code as part of a commercial system is 
    ! permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
    ! AUTHOR.  (If you are not directly supplying this 
    ! code to a customer, and you are instead telling them 
    ! how they can obtain it for free, then you are not 
    ! required to make any arrangement with me.) 
    !
    ! Disclaimer:  Neither I nor: Columbia University, the 
    ! National Aeronautics and Space Administration, nor 
    ! the Massachusetts Institute of Technology warrant 
    ! or certify this code in any way whatsoever.  This 
    ! code is provided "as-is" to be used at your own risk.
    !
    !

    !    
    ! FFSL1D.h90: upwind-biased flux-reconstruction scheme.
    !
    ! Darren Engwirda 
    ! 31-Mar-2019
    ! de2363 [at] columbia [dot] edu
    !
    !

    subroutine ffsl1d(npos,nvar,ndof,spac,tDEL, &
        &             mask,uvel,qbar,qedg,bclo, &
        &             bchi,work,opts)

    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! NDOF  no. degrees-of-freedom per grid-cell.
    ! SPAC  grid-cell spacing array. LENGTH(SPAC) == +1 if 
    !       spacing is uniform .
    ! TDEL  time-step .
    ! MASK  logical grid-cell masking array.
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QBAR  cell-centred integral moments. QBAR has SIZE =
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind flux eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    ! BCLO  boundary condition at lower endpoint .
    ! BCHI  boundary condition at upper endpoint . 
    ! WORK  method work-space. See RCON-WORK for details .
    ! OPTS  method parameters. See RCON-OPTS for details .
    !
    
        implicit none
    
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,ndof
        class(rmap_work), intent(inout):: work
        class(rmap_opts), intent(inout):: opts
        real*8 , intent(in)  :: spac(:)
        real*8 , intent(in)  :: tDEL
        logical, intent(in)  :: mask(:)
        real*8 , intent(in)  :: qbar(:,:,:)
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(out) :: qedg(:,:)
        class(rcon_ends), intent(in) :: bclo(:)
        class(rcon_ends), intent(in) :: bchi(:)
        
    !------------------------------------------- variables !      
        integer :: head,tail,nprt
    
        head = +0 ; tail = +0 ; qedg = 0.e+0
    
        do while (.true.)

    !--------------------------------- 1. find active part !
           
            do  head = tail+1, npos-1
            if (mask(head) .eqv..true.) exit
            end do
            
            do  tail = head+1, npos-1
            if (mask(tail).neqv..true.) exit
            end do
            tail = tail - 1

            if (head.ge.npos) exit
           
    !--------------------------------- 2. rcon active part !
            
            nprt = tail - head + 1
            
            if (size(spac).ne.+1) then
            
            call rcon1d(nprt+1,nvar,ndof , &
            &    spac(    head:tail), &
            &    qbar(:,:,head:tail), &
            &    bclo,bchi,work%cell_func, &
            &    work,opts )
            
            else
            
            call rcon1d(nprt+1,nvar,ndof , &
            &    spac,qbar(:,:,head:tail), &
            &    bclo,bchi,work%cell_func, &
            &    work,opts )
            
            end if
            
    !--------------------------------- 3. int. active part !

            select case(opts%cell_meth)
                case(pcm_method)       !! 1st-order scheme
   
                if (size(spac).ne.+1) then
                
                call flux1d(nprt+1,nvar,1, &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call flux1d(nprt+1,nvar,1, &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if

                case(plm_method)       !! 2nd-order scheme
   
                if (size(spac).ne.+1) then
                
                call flux1d(nprt+1,nvar,2, &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call flux1d(nprt+1,nvar,2, &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
       
                case(ppm_method)       !! 3rd-order scheme
   
                if (size(spac).ne.+1) then
                
                call flux1d(nprt+1,nvar,3, &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call flux1d(nprt+1,nvar,3, &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
        
                case(pqm_method)       !! 5th-order scheme
   
                if (size(spac).ne.+1) then
                
                call flux1d(nprt+1,nvar,5, &
                &    spac(  head:tail+0) , &
                &    tDEL, &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                else
                
                call flux1d(nprt+1,nvar,5, &
                &    spac,tDEL , &
                &    uvel(  head:tail+1) , &
                &    work%cell_func, &
                &    qedg(:,head:tail+1) )
                
                end if
            
            end select

        end do    
  
        return
        
    end  subroutine
    
    ! FLUX1D: a degree-k, upwind-type flux reconstruction. !

    pure subroutine flux1d(npos,nvar,mdof,SPAC, &
        &                  tDEL,uvel,QHAT,qedg)
    
    !
    ! NPOS  no. edges over grid.
    ! NVAR  no. state variables.
    ! MDOF  no. degrees-of-freedom per QHAT.
    ! SPAC  grid spacing vector. SIZE(SPAC)==+1 if uniform .
    ! TDEL  time-step .
    ! UVEL  edge-centred velocity vectors. UVEL has SIZE = 
    !       NPOS-by-1 .
    ! QHAT  cell-centred polynomial recon. QHAT has SIZE =  
    !       NDOF-by-NVAR-by-NPOS-1 .
    ! QEDG  edge-centred upwind flux eval. QEDG has SIZE = 
    !       NVAR-by-NPOS .
    !

        implicit none
        
    !------------------------------------------- arguments !
        integer, intent(in)  :: npos,nvar,mdof
        real*8 , intent(in)  :: SPAC(:)
        real*8 , intent(in)  :: tDEL
        real*8 , intent(in)  :: uvel(:)
        real*8 , intent(in)  :: QHAT(:,:,:)
        real*8 , intent(out) :: qedg(:,:)  
  
    !------------------------------------------- variables !      
        integer :: ipos,ivar
        real*8  :: uCFL,xhat,ss11,ss22,flux
        real*8  :: vv11(1:5)
        real*8  :: vv22(1:5)
        real*8  :: ivec(1:5)

    !----------- single-cell, lagrangian-type upwind rcon. !
        
        do  ipos = +2 , npos - 1
        
            if (uvel(ipos) .gt. +0.e0_8) then
            
    !----------- integrate profile over upwind cell IPOS-1 !
 
            if (size(SPAC).ne.+1) then
            xhat = .5d0 * SPAC(ipos-1)          
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(ipos-1)
            else
            xhat = .5d0 * SPAC(    +1)
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(    +1)
            end if
 
            ss11 = +1.e0_8 - 2.e0_8 * uCFL
            ss22 = +1.e0_8

            call bfun1d(-1,mdof,ss11,vv11)
            call bfun1d(-1,mdof,ss22,vv22)
            
            ivec =  vv22 - vv11
    
            do  ivar = +1, nvar

                flux =  dot_product (  &
        &           ivec(1:mdof), & 
        &       QHAT(1:mdof,ivar,ipos-1) )

                flux = flux * xhat

                qedg(ivar,ipos) = flux
                
            end do
                      
            else &
        &   if (uvel(ipos) .lt. -0.e0_8) then
            
    !----------- integrate profile over upwind cell IPOS+0 !
            
            if (size(SPAC).ne.+1) then
            xhat = .5d0 * SPAC(ipos-0)          
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(ipos-0)
            else
            xhat = .5d0 * SPAC(    +1)
            uCFL = uvel(ipos) &
        &        * tDEL / SPAC(    +1)
            end if

            ss11 = -1.e0_8 - 2.e0_8 * uCFL
            ss22 = -1.e0_8
        
            call bfun1d(-1,mdof,ss11,vv11)
            call bfun1d(-1,mdof,ss22,vv22)
            
            ivec =  vv22 - vv11
    
            do  ivar = +1, nvar

                flux =  dot_product (  &
        &           ivec(1:mdof), & 
        &       QHAT(1:mdof,ivar,ipos-0) )

                flux = flux * xhat

                qedg(ivar,ipos) = flux
                
            end do
            
            end if
        
        end do
               
        return

    end  subroutine





    !------------------------------------------ end ppr_1d !
      
        
    end module
    
    
    

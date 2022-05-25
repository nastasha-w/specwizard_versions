! this is the spline kernel used in gadget.
! it is modified from the W_4 spline kernel in 
! http://adsabs.harvard.edu/abs/1985A%26A...149..135M  eq. 21
! so that the kernel goes to zero at r=h instead of r=2h and
! normalized so that a [volume,surface,line] integral over the whole
! kernel in [3,2,1] dimensions gives 1.
!---------------------------------------------------------------------
!
! 
module w4_gadget_spline_kernel_class
implicit none
private

public :: kernel
public :: gaussian_kernel
public :: test_kernels
public :: gaussian_variance
public :: integrate_G1
public :: integrate_G2
public :: integrate_G3
public :: integrate_G3_line

! selected_int_kind(r)     -10^r < n < 10^r
!------------------------------------------------------------------------

integer, parameter :: i1b = selected_int_kind(2)  !< 1 byte integer 
integer, parameter :: i2b = selected_int_kind(4)  !< 2 byte integer 
integer, parameter :: i4b = selected_int_kind(9)  !< 4 byte integer 
integer, parameter :: i8b = selected_int_kind(18) !< 8 byte integer 

! selected_real_kind(p,r)  p is decimal precision, r is exponent range
!------------------------------------------------------------------------
integer, parameter :: r4b  = selected_real_kind(p=6,r=37)    !< 4 byte real 
integer, parameter :: r8b  = selected_real_kind(p=15,r=307)  !< 8 byte real 
integer, parameter :: r16b = selected_real_kind(p=33,r=4931) !< 16 byte real

! physical constants
!------------------------------------------------------------------------
real(r8b), parameter :: PI = 3.141592653589793238462643383279502884197d0
real(r8b), parameter :: SQRT_PI = 1.7724538509055158d0
real(r8b), parameter :: PI_one_third = 1.4645918875615231d0 ! = PI^(1/3)
real(r8b), parameter :: PI_one_sixth = 1.2102032422537642d0 ! = PI^(1/6)
real(r8b), parameter :: zero = 0.0d0
real(r8b), parameter :: half = 0.5d0
real(r8b), parameter :: one  = 1.0d0
real(r8b), parameter :: two = 2.0d0
real(r8b), parameter :: three = 3.0d0
real(r8b), parameter :: four = 4.0d0
real(r8b), parameter :: seven = 7.0d0
real(r8b), parameter :: eight = 8.0d0
real(r8b), parameter :: nine = 9.0d0
real(r8b), parameter :: thirty_two = 32.d0
real(r8b), parameter :: eighty = 80.d0

real(r8b), parameter :: four_thirds = 4.0d0 / 3.0d0
real(r8b), parameter :: sixteen_pi_over_nine = 16.0d0 * PI / 9.0d0

real(r8b), parameter :: forty_sevenths_over_pi = 40.0d0 / (7.0d0 * PI)
real(r8b), parameter :: forty_sevenths = 40.0d0 / 7.0d0 

real(r8b), parameter :: eight_over_pi = 8.0d0 / PI
real(r8b), parameter :: four_PI_one_third = 4.0d0 * PI_one_third


real(r8b), parameter :: one_third = 1.0d0 / 3.0d0
real(r8b), parameter :: eight_PI_one_third = eight * PI_one_third
real(r8b), parameter :: seven_over_eighty = seven / eighty
real(r8b), parameter :: half_SQRT_PI = half * SQRT_PI

real(r8b), parameter :: t1 = four * SQRT_PI / three
real(r8b), parameter :: C1 = 0.9991687417652023d0  ! = erf(t1)

real(r8b), parameter :: C2 = 0.996701494244061d0   ! = 1.0d0 - exp(-40.d0/7.d0)

real(r8b), parameter :: t3 = two * PI_one_sixth
real(r8b), parameter :: C3 = 0.9915807372143807d0  ! = erf(t3) - 2 / sqrt(PI) * exp(-t3^2) * t3


interface integrate_G1
   module procedure integrate_G1_s
end interface

interface integrate_G2
   module procedure integrate_G2_s, &
                    integrate_G2_v, &
                    integrate_G2_m
end interface

interface integrate_G3
   module procedure integrate_G3_s
end interface

interface kernel
   module procedure kernel_s, &
                    kernel_v
end interface


interface gaussian_kernel
   module procedure gaussian_kernel_s, &
                    gaussian_kernel_v
end interface

  
contains


  ! returns the 1-D integral of the equivalent 1-D gaussian from
  ! x1 to x2  serial version
  !--------------------------------------------------------------------
  function integrate_G1_s(hsml,x1,x2) result(S)
    real(8) :: hsml
    real(8) :: x1,x2
    real(8) :: N, A2, A
    real(8) :: Sx,S

    if (x1 > x2) stop "x1 gt x2 in integrate_G1_s"
    N = four_thirds / hsml
    A2 = sixteen_pi_over_nine / (hsml * hsml)
    A = sqrt(A2)
    Sx = ( erf(A*x2) - erf(A*x1) ) * half_SQRT_PI
    S = N * Sx / A

  end function integrate_G1_s


  ! returns the 2-D integral of the equivalent 2-D gaussian over 
  ! the area element whose lower corner is (x1,y1) and whose upper 
  ! corner is (x2,y2) serial version
  !--------------------------------------------------------------------
  function integrate_G2_s(hsml,x1,y1,x2,y2) result(S)
    real(8) :: hsml
    real(8) :: x1,y1,x2,y2
    real(8) :: N,A2,A
    real(8) :: Sx,Sy,S

    if (x1 > x2) stop "x1 gt x2 in integrate_G2_s"
    if (y1 > y2) stop "y1 gt y2 in integrate_G2_s" 

    N  = forty_sevenths_over_pi / (hsml * hsml)
    A2 = forty_sevenths / (hsml * hsml)
    A  = sqrt(A2)
    Sx = ( erf(A*x2) - erf(A*x1) ) * half_SQRT_PI
    Sy = ( erf(A*y2) - erf(A*y1) ) * half_SQRT_PI
    S = N * Sx * Sy / A2

  end function integrate_G2_s

  ! returns the 2-D integral of the equivalent 2-D gaussian over 
  ! the area element whose lower corner is (x1,y1) and whose upper 
  ! corner is (x2,y2) vector version
  !--------------------------------------------------------------------
  function integrate_G2_v(hsml,x1,y1,x2,y2) result(S)
    real :: hsml
    real, dimension(:) :: x1,y1,x2,y2
    real, dimension(size(x1)) :: Sx,Sy,S
    real :: N,A2,A

    if (any(x1 > x2)) stop "x1 gt x2 in integrate_G2_v"
    if (any(y1 > y2)) stop "y1 gt y2 in integrate_G2_v"

    N  = forty_sevenths_over_pi / (hsml * hsml)
    A2 = forty_sevenths / (hsml * hsml)
    A  = sqrt(A2)
    Sx = ( erf(A*x2) - erf(A*x1) ) * half_SQRT_PI
    Sy = ( erf(A*y2) - erf(A*y1) ) * half_SQRT_PI
    S = N * Sx * Sy / A2

  end function integrate_G2_v


  ! returns the 2-D integral of the equivalent 2-D gaussian over 
  ! the area element whose lower corner is (x1,y1) and whose upper 
  ! corner is (x2,y2) matrix version
  !--------------------------------------------------------------------
  function integrate_G2_m(hsml,x1,y1,x2,y2) result(S)
    real :: hsml
    real, dimension(:,:) :: x1,y1,x2,y2
    real, dimension(size(x1,1),size(x1,2)) :: Sx,Sy,S
    real :: N,A2,A

    if (any(x1 > x2)) stop "x1 gt x2 in integrate_G2_m"
    if (any(y1 > y2)) stop "y1 gt y2 in integrate_G2_m"

    N  = forty_sevenths_over_pi / (hsml * hsml)
    A2 = forty_sevenths / (hsml * hsml)
    A  = sqrt(A2)
    Sx = ( erf(A*x2) - erf(A*x1) ) * half_SQRT_PI
    Sy = ( erf(A*y2) - erf(A*y1) ) * half_SQRT_PI
    S = N * Sx * Sy / A2

  end function integrate_G2_m


  ! returns the 3-D integral of the equivalent 3-D gaussian over 
  ! the volume element whose lower corner is (x1,y1,z1) and whose upper 
  ! corner is (x2,y2,z1) serial version
  !--------------------------------------------------------------------
  function integrate_G3_s(hsml,x1,y1,z1,x2,y2,z2) result(S)
    real(8) :: hsml
    real(8) :: x1,y1,z1,x2,y2,z2
    real(8) :: N,A2,A3,A
    real(8) :: Sx,Sy,Sz,S

    if (x1 > x2) stop "x1 gt x2 in integrate_G3_s"
    if (y1 > y2) stop "y1 gt y2 in integrate_G3_s" 
    if (z1 > z2) stop "z1 gt z2 in integrate_G3_s" 

    N  = eight_over_pi / (hsml * hsml * hsml)
    A2 = four_pi_one_third / (hsml * hsml)    
    A  = sqrt(A2)
    A3 = A2 * A
    Sx = ( erf(A*x2) - erf(A*x1) ) * half_SQRT_PI
    Sy = ( erf(A*y2) - erf(A*y1) ) * half_SQRT_PI
    Sz = ( erf(A*z2) - erf(A*z1) ) * half_SQRT_PI
    S = N * Sx * Sy * Sz / A3

  end function integrate_G3_s





  ! returns the line integral thru a 3-D gaussian at impact parameter b 
  ! from zi to zf.  this gaussian is normalized so the 3-D volume 
  ! integral over the whole smoothing volume = 1
  !--------------------------------------------------------------------
  function integrate_G3_line(h,b,zi,zf) result(S)
    real(8), intent(in) :: h
    real(8), intent(in) :: b
    real(8), intent(in) :: zi
    real(8), intent(in) :: zf
    real(8) :: S
    real(8) :: var,fac, N, A2, A
    real(8) :: xi,xf,zedge,derf

    if (zi > zf) stop "zi > zf in integrate_G3_line"
    if (b > h) stop "b > h in integrate_G3_line"

    N = eight_over_PI / (h*h*h)
    var = h * h / (eight_PI_one_third)
    A2 = 1.0d0 / (2.0d0*var)
    A = sqrt(A2)
    fac = exp(-A2*b*b)
    
    zedge = sqrt(h*h-b*b)
    
    if ( zf .ge. 0.d0 ) then 
        xf = min(zf,zedge)
    else 
        xf = max(zf,-zedge)
    endif
    if ( zi .ge. 0.d0 ) then 
        xi = min(zi,zedge)
    else 
        xi = max(zi,-zedge)
    endif
    
    xi = A * xi
    xf = A * xf
    derf = erf(xf) - erf(xi)
    S = half * SQRT_PI * N * fac / (C3 * A) * (erf(xf) - erf(xi))

!    write(*,*) 'xi,xf: ', xi,xf
!    write(*,*) 'S = ', S

  end function integrate_G3_line







  ! performs tests on the kernels
  !--------------------------------------------------------------
  function test_kernels() result(err)
    integer(8), parameter :: N = 1000*4 + 1
    real(8), parameter :: hsml = 4.0d0

    real(8), dimension(0:N-1) :: r, f1, f2, f3
    real(8) :: a, b, dr
    real(8) :: sum1_t, sum2_t, sum3_t
    real(8) :: sum1_b, sum2_b, sum3_b
    real(8) :: zi, zf
    real(8) :: x1, x2, y1, y2, z1, z2
    real(8) :: line_int, surf_int, vol_int

    integer :: i, err

    a = 0.0d0 
    b = hsml
    dr = (b-a)/(N-1)

    write(*,*) 
    write(*,*) 'Running kernel tests ...'
    write(*,*) 
    write(*,*) 'Constants'
    write(*,*) '---------'
    write(*,*) "PI       = ", PI
    write(*,*) "Sqrt[PI] = ", sqrt(PI)
    write(*,*) "PI^(1/3) = ", PI**(1.0d0/3.0d0)
    write(*,*) "PI^(1/6) = ", PI**(1.0d0/6.0d0)
    write(*,*) "erf(2.0) = ", erf(2.0d0)
    write(*,*) 
    write(*,*) "t1 = 4 sqrt(PI) / 3 = ", t1
    write(*,*) "t3 = 2 PI^(1/6)     = ", t3
    write(*,*) 
    write(*,*) "C1 = erf(t1)                                  = ", erf(t1)
    write(*,*) "C2 = 1.0d0 - exp(-40.d0/7.d0)                 = ", 1.0d0 - exp(-40.d0/7.d0)
    write(*,*) "C3 = erf(t3) - 2 / sqrt(PI) * exp(-t3^2) * t3 = ", erf(t3) - 2 / sqrt(PI) * exp(-t3*t3) * t3
    write(*,*) 
    write(*,*) 'Volume Integrals'
    write(*,*) '----------------'
    write(*,*) 
    write(*,*) 'N = ', N
    write(*,*) 'hsml = ', hsml
    write(*,*) 
    write(*,*) 'W4(r=0,hsml)  1D: ', kernel(0.0d0, hsml, ndim=1)
    write(*,*) 'W4(r=0,hsml)  2D: ', kernel(0.0d0, hsml, ndim=2)
    write(*,*) 'W4(r=0,hsml)  3D: ', kernel(0.0d0, hsml, ndim=3)
    write(*,*) 
    write(*,*) 'W4(r=hsml/2,hsml)  1D: ', kernel(hsml/2.0d0, hsml, ndim=1)
    write(*,*) 'W4(r=hsml/2,hsml)  2D: ', kernel(hsml/2.0d0, hsml, ndim=2)
    write(*,*) 'W4(r=hsml/2,hsml)  3D: ', kernel(hsml/2.0d0, hsml, ndim=3)
    write(*,*) 

    do i = 0, N-1
       r(i) = i*dr
    end do

    f1 = kernel(r, hsml, ndim=1)
    f2 = 2.0d0 * PI * r * kernel(r,hsml,ndim=2)
    f3 = 4.0d0 * PI * r * r * kernel(r,hsml,ndim=3)       

    sum1_b = int_tabulated_boole( r, f1 )
    sum2_b = int_tabulated_boole( r, f2 )
    sum3_b = int_tabulated_boole( r, f3 )

    write(*,*) "W4 Volume Integral (1D) = ", 2*sum1_b
    write(*,*) "W4 Volume Integral (2D) = ", sum2_b
    write(*,*) "W4 Volume Integral (3D) = ", sum3_b
    write(*,*)

    f1 = gaussian_kernel_v(r, hsml, ndim=1)
    f2 = 2.0d0 * PI * r * gaussian_kernel_v(r,hsml,ndim=2)
    f3 = 4.0d0 * PI * r * r * gaussian_kernel_v(r,hsml,ndim=3)       

    sum1_b = int_tabulated_boole( r, f1 )
    sum2_b = int_tabulated_boole( r, f2 )
    sum3_b = int_tabulated_boole( r, f3 )

    write(*,*) "Wg Volume Integral (1D) = ", 2*sum1_b
    write(*,*) "Wg Volume Integral (2D) = ", sum2_b
    write(*,*) "Wg Volume Integral (3D) = ", sum3_b
    write(*,*) 
    
    b = 0.0d0
    zi = -hsml
    zf = hsml
    line_int = integrate_G3_line(hsml,b,zi,zf) 
    write(*,*) "Wg I(h,b=0,zi=-h,zf=h): ", line_int

    f3 = kernel_v(r,hsml,ndim=3)       
    line_int = 2.0d0 * int_tabulated_boole( r, f3 )
    write(*,*) "W4 I(h,b=0,zi=-h,zf=h): ", line_int
    write(*,*) 

    x1=-hsml
    x2=hsml
    line_int = integrate_G1_s(hsml,x1,x2)
    write(*,*) 'line int: ', line_int
    write(*,*) 

    x1=-hsml
    x2=hsml
    y1=-hsml
    y2=hsml
    surf_int = integrate_G2_s(hsml,x1,y1,x2,y2)
    write(*,*) 'surf int: ', surf_int
    write(*,*) 

    x1=-hsml
    x2=hsml
    y1=-hsml
    y2=hsml
    z1=-hsml
    z2=hsml
    vol_int = integrate_G3_s(hsml,x1,y1,z1,x2,y2,z2)
    write(*,*) 'vol int: ', vol_int
    write(*,*) 

    err = 0


  end function test_kernels


  ! integrates tabulated data using the trap. rule ( 2 pt N-C )
  !--------------------------------------------------------------
  function int_tabulated_trap(x,y) result(sum)
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: y(size(x))
    real(8)  :: sum
     
    real(8) :: dx
    integer :: N, i

    N = size(x)
    sum = 0.0d0
    do i = 1,N-1
       dx = x(i+1) - x(i)
       sum = sum + 0.5d0 * (y(i) + y(i+1)) * dx 
    end do

  end function int_tabulated_trap


  ! integrates tabulated data using Booles rule ( 5 pt N-C )
  !--------------------------------------------------------------
  function int_tabulated_boole(x,y) result(sum)
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: y(size(x))
    real(8)  :: fac
    real(8)  :: sum
     
    real(8) :: dx
    integer :: N, i, j, Ngroup
  
    N = size(x)

    if ( mod(N,4) /= 1 ) then
       write(*,*) 'For NC 5 pt, array must satisfy mod(N,4) = 1'
       stop
    endif

    Ngroup = (N-1)/4

    sum = 0.0d0
    do i = 1,Ngroup
       j = (i-1) * 4 + 1
       dx = x(j+1) - x(j)
       fac = 7.0d0 * (y(j)   + y(j+4)) + &
            32.0d0 * (y(j+1) + y(j+3)) + 12.0d0 * y(j+2)
       sum = sum + 2 * dx / 45.0d0 * fac
    end do

  end function int_tabulated_boole




  ! returns variance = sigma^2 of the equivalent gaussian
  ! (i.e. the normalized gaussian that has the same y-intercept)
  !--------------------------------------------------------------
  function gaussian_variance(hsml,ndim) result(var)
    real :: hsml
    integer, optional :: ndim
    real :: var
    integer :: nd

    if (present(ndim)) then
       nd = ndim
    else
       nd = 3
    endif
    
    select case (nd)
    case(3)
       var = hsml * hsml / (eight * PI**(one_third))
    case(2)
       var = hsml * hsml * seven / eighty
    case(1)
       var = hsml * hsml * nine / (thirty_two * PI)
    case default
       stop "ndim out of range in gaussian variance"
    end select


  end function gaussian_variance

  ! serial version suitable for binding with a particle class
  !-----------------------------------------------------------
  function kernel_s(r, h, ndim) result(W)
    real(8), intent(in) :: r
    real(8) :: h
    integer, intent(in), optional :: ndim
    real(8) :: W

    integer :: nd
    real(8) :: norm
    real(8) :: x, x2, x3

    x  = r / h
    x2 = x * x
    x3 = x2 * x
    W = zero    

    if (present(ndim)) then
       nd=ndim
    else
       nd=3
    endif

    select case (nd)

    case(1)
       norm = four_thirds / h
       if (abs(x) <= half ) then
          W = 1 - 6 * x2 + 6 * x3
       else if (abs(x) > half .and. abs(x) < one ) then
          W = 2 * (1-x)*(1-x)*(1-x)
       endif
       W = W * norm

    case(2)
       norm = forty_sevenths_over_pi / (h*h)
       if ( r < zero ) then
          write(*,*) "ndim > 1 and r < 0.0"
          stop
       endif

       if ( x >= zero .and. x <= half ) then
          W = 1 - 6 * x2 + 6 * x3
       else if ( x > half .and. x < one ) then
          W = 2 * (1-x)*(1-x)*(1-x)
       endif
       W = W * norm

    case(3)
       norm = eight_over_pi / (h * h * h)
       if ( r < zero ) then
          write(*,*) "ndim > 1 and r < 0.0"
          stop
       endif
       
       if ( x >= zero .and. x <= half ) then
          W = 1 - 6 * x2 + 6 * x3
       else if ( x > half .and. x < one ) then
          W = 2 * (1-x)*(1-x)*(1-x)
       endif
       W = W * norm

    case default
       write(*,*) "ndim out of range"
       stop

    end select
    
    
    return


  end function kernel_s

   
  ! accepts array r, single h, and optional ndim
  !-----------------------------------------------
  function kernel_v(r, h, ndim) result(W)
    real(8), intent(in), dimension(:) :: r
    real(8), intent(in) :: h
    integer, intent(in), optional :: ndim
    real(8), dimension(size(r)) :: W

    integer :: nd
    real(8) :: norm
    real(8), dimension(size(r)) :: x, x2, x3

    x  = r / h
    x2 = x * x
    x3 = x2 * x
    W = zero    

    if (present(ndim)) then
       nd=ndim
    else
       nd=3
    endif

    select case (nd)

    case(1)
       norm = four_thirds / h
       where( abs(x) <= half )
          W = 1 - 6 * x2 + 6 * x3
       elsewhere( abs(x) > half .and. abs(x) < one )
          W = 2 * (1-x)*(1-x)*(1-x)
       endwhere
       W = W * norm

    case(2)
       norm = forty_sevenths_over_pi / (h*h)
       if ( count( r < zero ) > 0 ) then
          write(*,*) "ndim > 1 and r < 0.0"
          stop
       endif

       where( x >= zero .and. x <= half )
          W = 1 - 6 * x2 + 6 * x3
       elsewhere( x > half .and. x < one )
          W = 2 * (1-x)*(1-x)*(1-x)
       endwhere
       W = W * norm

    case(3)
       norm = eight_over_pi / (h * h * h)
       if ( count( r < zero ) > 0 ) then
          write(*,*) "ndim > 1 and r < 0.0"
          stop
       endif

       where( x >= zero .and. x <= half )
          W = 1 - 6 * x2 + 6 * x3
       elsewhere( x > half .and. x < one )
          W = 2 * (1-x)*(1-x)*(1-x)
       endwhere
       W = W * norm

    case default
       write(*,*) "ndim out of range"
       stop

    end select
    
    
    return


  end function kernel_v








  ! accepts scalar r, single h, and optional ndim
  !-----------------------------------------------
  function gaussian_kernel_s(r, h, ndim) result(Wg)
    real(8), intent(in) :: r
    real(8), intent(in) :: h
    integer, intent(in), optional :: ndim
    real(8) :: Wg

    integer :: nd
    real(8) :: norm
    real(8) :: var
    real(8) :: r2

    r2 = r * r
    Wg = zero    

    if (present(ndim)) then
       nd=ndim
    else
       nd=3
    endif

    select case (nd)

    case(1)
       var = h * h * nine / (thirty_two * PI)
       norm = four_thirds / h / C1
       if ( abs(r) <= h ) then
          Wg = norm * exp( - r2 / (2.0d0*var) )
       else
          Wg = zero
       endif

    case(2)
       var = h * h * seven / eighty
       norm = forty_sevenths_over_pi / (h*h) / C2
       if ( r < zero ) then
          write(*,*) "ndim > 1 and r < 0.0"
          stop
       endif

       if( r <= h ) then
          Wg = norm * exp( - r2 / (2.0d0*var) )
       else
          Wg = zero
       endif

    case(3)
       var = h * h / (eight * PI_one_third)
       norm = eight_over_pi / (h * h * h) / C3
       if ( r < zero ) then
          write(*,*) "ndim > 1 and r < 0.0"
          stop
       endif

       if( r <= h ) then
          Wg = norm * exp( - r2 / (2.0d0*var) )
       else
          Wg = zero
       endif

    case default
       write(*,*) "ndim out of range"
       stop

    end select
    
    
    return


  end function gaussian_kernel_s







  ! accepts array r, single h, and optional ndim
  !-----------------------------------------------
  function gaussian_kernel_v(r, h, ndim) result(Wg)
    real(8), intent(in), dimension(:) :: r
    real(8), intent(in) :: h
    integer, intent(in), optional :: ndim
    real(8), dimension(size(r)) :: Wg

    integer :: nd
    real(8) :: norm
    real(8) :: var
    real(8), dimension(size(r)) :: r2

    r2 = r * r
    Wg = zero    

    if (present(ndim)) then
       nd=ndim
    else
       nd=3
    endif

    select case (nd)

    case(1)
       var = h * h * nine / (thirty_two * PI)
       norm = four_thirds / h / C1
       where( abs(r) <= h )
          Wg = norm * exp( - r2 / (2.0d0*var) )
       elsewhere
          Wg = zero
       endwhere

    case(2)
       var = h * h * seven / eighty
       norm = forty_sevenths_over_pi / (h*h) / C2
       if ( count( r < zero ) > 0 ) then
          write(*,*) "ndim > 1 and r < 0.0"
          stop
       endif

       where( r <= h )
          Wg = norm * exp( - r2 / (2.0d0*var) )
       elsewhere
          Wg = zero
       endwhere

    case(3)
       var = h * h / (eight * PI_one_third)
       norm = eight_over_pi / (h * h * h) / C3
       if ( count( r < zero ) > 0 ) then
          write(*,*) "ndim > 1 and r < 0.0"
          stop
       endif

       where( r <= h )
          Wg = norm * exp( - r2 / (2.0d0*var) )
       elsewhere
          Wg = zero
       endwhere

    case default
       write(*,*) "ndim out of range"
       stop

    end select
    
    
    return


  end function gaussian_kernel_v



end module w4_gadget_spline_kernel_class
 



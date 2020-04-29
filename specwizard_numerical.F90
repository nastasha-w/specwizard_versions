!
! Source file contains 'numerical' modules, that is:  spline interpolation, rebinning, random numbers, convolution
!
module my_spline_interpolation
  use numbers
  interface spline_interpolation
     module procedure spline_interpolate
  end interface
contains
  subroutine spline_interpolate(nin,inx,iny,minvoc,maxvoc,ion,minbother,nout,outx,outy,loginterpolate,is_positive)
    use numbers
    use spectra, only : work, work2, ions
    use cpu_timers
    implicit none
    integer, intent(in)               :: nin,nout
    real(kind=doubleR), intent(in), dimension(:)    :: inx, iny, outx
    real(kind=doubleR), intent(inout), dimension(:) :: outy

    real(kind=doubleR), intent(in) :: minvoc, maxvoc, minbother

    integer, intent(in)               :: ion
    logical, optional, intent(in)     :: loginterpolate, is_positive
    ! local
    real(kind=doubleR), parameter     :: big=2.d30
    real(kind=doubleR)                :: result
    integer i
    logical :: loginterpolate_use, is_positive_use
    !    
    ! set defaults
    if(present(is_positive))then
       is_positive_use = is_positive
    else
       is_positive_use = .false.
    endif
    if(present(loginterpolate))then
       loginterpolate_use = loginterpolate
    else
       loginterpolate_use = .false.
    endif
    !
    call cpu_timer_start(dointerpolate)
    call spline(inx,iny,nin,big,big,work,work2)
    do i=1, nout
       if(outx(i) .ge. minvoc .and. outx(i) .lt. maxvoc)then
          call splint(inx,iny,work,nin,outx(i),result)
          if(loginterpolate_use)then
             if(result .lt. -minbother) then
                result = linint(nin,inx,iny,outx(i),ions(ion))
             endif
          endif
          if(is_positive_use)then
             if(result .lt. 0) result = 0.
          endif
       else
          result = 0.0
       endif
       outy(i) = outy(i) + result
    enddo
    call cpu_timer_stop(dointerpolate)

  end subroutine spline_interpolate
  SUBROUTINE spline(x,y,n,yp1,ypn,y2,u)
    ! set-up spline interpolation coefficients
    use numbers
    implicit none 

    INTEGER(kind=singleI) ::  n
    REAL(kind=doubleR), intent(in)  :: x(:), y(:), yp1, ypn
    real(kind=doubleR), intent(out) :: y2(:), u(:)
    !
    INTEGER i,k
    REAL(kind=doubleR) ::  p,qn,sig,un
    !
    if (yp1.gt..99e30) then
       y2(1)=0.
       u(1)=0.
    else
       y2(1)=-0.5
       u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
    endif
    do i=2,n-1
       sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
       p=sig*y2(i-1)+2.
       y2(i)=(sig-1.)/p
       u(i)=(6.*((y(i+1)-y(i))/(x(i+ &
            1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig* &
            u(i-1))/p
    enddo
    if (ypn.gt..99e30) then
       qn=0.
       un=0.
    else
       qn=0.5
       un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
    endif
    y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
    do k=n-1,1,-1
       y2(k)=y2(k)*y2(k+1)+u(k)
    enddo
    return

  END subroutine spline

  SUBROUTINE splint(xa,ya,y2a,n,x,y)
    ! do spline interpolation
    use numbers
    implicit none 

    INTEGER(kind=singleI), intent(in) ::  n

    REAL(kind=doubleR), intent(in)    ::  x,xa(:),y2a(:),ya(:)	
    real(kind=doubleR), intent(out)   ::  y
    ! local variables
    INTEGER k
    ! CMB bugfix, declare type only once
    ! ,khi,klo
    REAL*8 a,b,h
    !     Added by JS:
    integer(kind=singleI) ::  klo, khi
    common /splint_ind/ klo, khi ! These indices bracket the input value of x

    klo=1
    khi=n
1   if (khi-klo.gt.1) then
       k=(khi+klo)/2
       if(xa(k).gt.x)then
          khi=k
       else
          klo=k
       endif
       goto 1
    endif
    h=xa(khi)-xa(klo)
    if (h.eq.0.) stop 'bad xa input in splint'
    a=(xa(khi)-x)/h
    b=(x-xa(klo))/h
    y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))* &
         (h**2)/6.

    return
  END SUBROUTINE splint

  function linint(n,x,y,xx,ionname)
    use numbers
    use runtime
    use my_mpi
    implicit none 
    real(kind=doubleR) :: linint
    integer(kind=singleI), intent(in) :: n
    character(*), intent(in)         :: ionname
    real(kind=doubleR), intent(in)   :: x(n),y(n)
    real(kind=doubleR), intent(in)   :: xx

    !  
    integer(kind=singleI) ::  klo, khi
    common /splint_ind/ klo, khi

    linint = ((xx-x(klo)) * y(khi) + (x(khi)-xx) * y(klo)) / (x(khi) - x(klo))

    return
  end function linint
end module my_spline_interpolation

module my_rebin
   contains

   subroutine rebin(input,output)
     use numbers
     use spectra, only: lambda, binned_lambda,nvpix,n_binned_flux,pixsize
     implicit none
     real(kind=doubleR), intent(in)  :: input(:)
     real(kind=doubleR), intent(out) :: output(:)
   
     ! local
     integer :: i,j,k,ninbin
     !write(*,*) 'Rebin (wavelength space)'
     !write(*,'("lambda: ",f7.4,", ", f7.4, ", ", f7.4, " ... ", f7.4)') & 
     !     lambda(1), lambda(2), lambda(3), lambda(nvpix)
     !write(*,'("binned_lambda: ",f7.4,", ", f7.4, ", ",f7.4," ... ",f7.4)') & 
     !     binned_lambda(1), binned_lambda(2), binned_lambda(3), binned_lambda(n_binned_flux)

     j = 1
     do i = 1, n_binned_flux
        ninbin = 0
        do while (lambda(j) .lt. binned_lambda(i)-0.5*pixsize .and. j.lt. nvpix)
           j = j + 1
        enddo
        k = j
        do while (lambda(k) .le. binned_lambda(i)+0.5*pixsize .and. k .le. nvpix)
           ninbin = ninbin + 1
           output(i) = output(i) + input(k)
           k = k + 1
        enddo
        if (ninbin .eq. 0) then
           if (i .lt. nvpix) then 
              write(*,'("ERROR: Grid too coarse! Binned_lambda = ", f10.3)') binned_lambda(i)
              stop
           endif
        else 
           output(i) = output(i) / dble(ninbin)
        endif
     enddo
   end subroutine rebin
   
   subroutine rebin_realspace(input,output)
     use numbers
     use spectra, only: redshift_realspace, binned_redshift_realspace, nppix,&
                        n_binned_realspace, pixsize, maxlambda
     implicit none
     real(kind=doubleR), intent(in)  :: input(:)
     real(kind=doubleR), intent(out) :: output(:)
   
     ! local
     integer :: i,j,k,ninbin
     real(kind=doubleR) :: pixsize_real
     !
     pixsize_real = pixsize / maxlambda
     !
     !write(*, *) 'Rebin (z-space, no pec. vel.)'
     !write(*,'("redshift_realspace: ",f10.6,", ", f10.6, ", ", f10.6, " ... ", f10.6)') & 
     !     redshift_realspace(1), redshift_realspace(2), redshift_realspace(3), &
     !     redshift_realspace(nppix)
     !write(*,'("binned_redshift_realspace: ",f10.6,", ", f10.6, ", ",f10.6," ... ",f10.6)') & 
     !     binned_redshift_realspace(1), binned_redshift_realspace(2), &
     !     binned_redshift_realspace(3), binned_redshift_realspace(n_binned_realspace)
     !
     j = 1
     do i = 1, n_binned_realspace
        ninbin = 0
        do while (redshift_realspace(j) .lt. binned_redshift_realspace(i)-0.5*pixsize_real &
                  .and. j .lt. nppix)
           j = j + 1
        enddo
        k = j
        do while (redshift_realspace(k) .le. binned_redshift_realspace(i)+0.5*pixsize_real &
                 .and. k .le. nppix)
           ninbin = ninbin + 1
           output(i) = output(i) + input(k)
           k = k + 1
        enddo
        if (ninbin .eq. 0) then
           if (i .lt. nppix) then 
              write(*,'("ERROR: Grid too coarse! Binned_lambda = ", f10.3)') binned_redshift_realspace(i)
              stop
           endif
        else 
           output(i) = output(i) / dble(ninbin)
        endif
     enddo
   end subroutine rebin_realspace

end module my_rebin


module uniform_deviate_module

  contains


   DOUBLE PRECISION FUNCTION RAN3(IDUM)
     SAVE
     PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.E-9)
     DIMENSION MA(55)
     DATA IFF /0/
     
     IF(IDUM.LT.0.OR.IFF.EQ.0)THEN
        IFF=1
        MJ=MSEED-IABS(IDUM)
        MJ=MOD(MJ,MBIG)
        MA(55)=MJ
        MK=1
        DO 11 I=1,54
           II=MOD(21*I,55)
           MA(II)=MK
           MK=MJ-MK
           IF(MK.LT.MZ)MK=MK+MBIG
           MJ=MA(II)
   11      CONTINUE
           DO 13 K=1,4
              DO 12 I=1,55
                 MA(I)=MA(I)-MA(1+MOD(I+30,55))
                 IF(MA(I).LT.MZ)MA(I)=MA(I)+MBIG
   12            CONTINUE
   13            CONTINUE
                 INEXT=0
                 INEXTP=31
                 IDUM=1
      ENDIF
      INEXT=INEXT+1
      IF(INEXT.EQ.56)INEXT=1
      INEXTP=INEXTP+1
      IF(INEXTP.EQ.56)INEXTP=1
      MJ=MA(INEXT)-MA(INEXTP)
      IF(MJ.LT.MZ)MJ=MJ+MBIG
      MA(INEXT)=MJ
      RAN3=MJ*FAC
      !     Added by JS to get around a bug (ran3 rarely returns a neg. number)
      if (ran3 .lt. 0.) ran3=abs(ran3)
      if (ran3 .gt. 1.) ran3=mod(ran3,1.)
      RETURN
   END function ran3

end module uniform_deviate_module

module gaussian_deviate_module

   contains

   FUNCTION gasdev(idum)
     use numbers
     use uniform_deviate_module
     implicit none 
     INTEGER idum
     real(kind=doubleR) :: gasdev
     !U    USES ran3
     INTEGER iset
     REAL(kind=doubleR)  fac,gset,rsq,v1,v2
     SAVE iset,gset
     DATA iset/0/
     if (iset.eq.0) then
   1    v1=2.*ran3(idum)-1.
        v2=2.*ran3(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
     else
        gasdev=gset
        iset=0
     endif
     return
   END FUNCTION gasdev
   

   end module gaussian_deviate_module

module my_random_numbers

   contains


   subroutine init_random_numbers(ispec)
     use numbers
     use random_numbers
     implicit none
     integer, intent(in) :: ispec
     ! local variables
     integer :: i
     !
     cycles = 0
     ! 
     seed   = seed0 + ispec - 1
     
   end subroutine init_random_numbers
     
   function random(itype)
     use numbers
     use random_numbers
     use gaussian_deviate_module, only: gasdev
     use uniform_deviate_module, only: ran3
     implicit none
     integer(kind=singleI), intent(in) :: itype
     real(kind=doubleR)                :: random
     integer :: i, iseed, icyc, type
     real(kind=doubleR)    :: dummy

     type = itype
     if(abs(type) .gt. ntypes) &
          call abortrun(' wrong random function')
   
     if(type .lt. 0)then
        ! re-initialize this counter
        type          = - type
        cycles(type)  = 0
     endif
   
     if(cycles(type) .eq. 0)then
        ! initialize
        iseed = seed(type)
        if(random_type(type) == 0)then
           do i=1, nran_max
              randoms(type,i) = ran3(iseed)
           enddo
        else
           do i=1, nran_max
              randoms(type,i) = gasdev(iseed)
           enddo
        endif
        cycles(type)  = 1
        current(type) = 0
     endif
   
     ! 
     if(current(type) .ge. nran_max)then
        ! need to re-initialize
        iseed = seed(type)
        do icyc=1, cycles(type)
           if(random_type(type) == 0)then
              do i=1, nran_max
                 dummy = ran3(iseed)
              enddo
           else
              do i=1, nran_max
                 dummy = gasdev(iseed)
              enddo
           endif
        enddo
        cycles(type)  = cycles(type)+1
        current(type) = 0
        if(random_type(type) == 0)then
           do i=1, nran_max
              randoms(type,i) = ran3(iseed)
           enddo
        else
           do i=1, nran_max
              randoms(type,i) = gasdev(iseed)
           enddo
        endif
     endif
     !
     current(type)  = current(type) + 1
     random         = randoms(type,current(type))
   end function random


end module my_random_numbers


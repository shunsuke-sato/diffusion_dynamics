module global_variables
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)


  integer :: nx, nt
  real(8) :: length_x, dx
  real(8) :: Tprop, dt
  real(8),allocatable :: rho(:), Dx(:)
  real(8) :: D_coeff, tau
  
  
end module global_variables
!------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input_variables
  call initialize

  call time_propagation
  
end program main
!------------------------------------------------------------------------------
subroutine input_variables
  use global_variables
  implicit none

  length_x = 100d0
  nx = 64
  dx = length_x/nx

  Tprop = 100d0
  dt = 0.1d0
  nt = aint(Tprop/dt)+1
  
end subroutine input_variables
!------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  integer :: ix
  real(8) :: xx

  
  D_coeff = 1d0
  tau = 1d40

  allocate(rho(0:nx-1), Dx(0:nx-1))

  do ix = 0, nx-1
     xx = dx*ix
     rho(ix) = 2d0*sin(pi*xx/length_x)**2/length_x
  end do
  
end subroutine initialize
!------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  character(256) :: cit, cfilename
  real(8) :: n0, n1, n2

  it = 0
  write(cit, "(I7.7)")it
  cfilename = trim(cit)//"_rho.out"
  open(20,file=cfilename)
  do ix = 0, nx-1
     write(20,"(99e26.16e3)")dx*ix,rho(ix)
  end do
  close(20)

  open(30,file="rho_dft.out")
  call calc_dft(n0,n1,n2,n3)
  write(30,"(999e26.16e3)")dt*it,n0,n1,n2,n3
  
  do it = 0, nt
     call dt_evolve
     call calc_dft(n0,n1,n2,n3)
     write(30,"(999e26.16e3)")dt*(it+1),n0,n1,n2,n3     

     if(mod(it+1,100)==0)then
        write(cit, "(I7.7)")it
        cfilename = trim(cit)//"_rho.out"
        open(20,file=cfilename)
        do ix = 0, nx-1
           write(20,"(99e26.16e3)")dx*ix,rho(ix)
        end do
        close(20)
     end if
     
  end do
  close(30)
  
  
end subroutine time_propagation
!------------------------------------------------------------------------------
subroutine dt_evolve
  use global_variables
  implicit none
  real(8) :: vec1_t(0:nx-1),vec2_t(0:nx-1)
  integer :: ix

  Dx = D_coeff


  vec1_t(0) = 0.5d0*(rho(1)-rho(nx-1))/dx
  do ix = 1, nx-2
     vec1_t(ix) = 0.5d0*(rho(ix+1)-rho(ix-1))/dx
  end do
  vec1_t(nx-1) = 0.5d0*(rho(0)-rho(nx-2))/dx

  vec1_t = Dx*vec1_t

  vec2_t(0) = 0.5d0*(vec1_t(1)-vec1_t(nx-1))/dx
  do ix = 1, nx-2
     vec2_t(ix) = 0.5d0*(vec1_t(ix+1)-vec1_t(ix-1))/dx
  end do
  vec2_t(nx-1) = 0.5d0*(vec1_t(0)-vec1_t(nx-2))/dx

  rho = vec2_t
  
  
end subroutine dt_evolve
!------------------------------------------------------------------------------
subroutine calc_dft(n0,n1,n2)
  use global_variables
  implicit none
  real(8),intent(out) :: n0,n1,n2
  integer,parameter :: np = 3
  complex(8) :: zn(0:np)
  integer :: ix, n
  real(8) :: xx

  do n = 0, np
     zn(n) = 0d0
     do ix = 0, nx-1
        xx = dx*ix
        zn(n) = zn(n) + rho(ix)*exp(zi*n*2d0*pi*xx/length_x)
     end do
  end do
  zn = zn*dx
  
  n0 = abs(zn(0))**2
  n1 = abs(zn(1))**2
  n2 = abs(zn(2))**2
  n3 = abs(zn(3))**2

  
end subroutine calc_dft
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

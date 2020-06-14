module global
implicit none
real*8 :: j
end module global

!--------------------------------------------------------------------

program Isotropic_Scatt_V3
  use global
  implicit none

  integer :: i, npackets,k,ok
  real*8 :: x,y,z,nx,ny,nz,rmax,L,albedo, ran, r
  integer, Dimension (1) :: Seed
  real*8 :: taumax,t,t1
  character (len=80) :: name
  
  name = "IsotropicCheck.dat"
  open(1,file=name,status='replace',iostat=ok)
  close(1)
  
  albedo=1.0 !always a scattering event here

  
  Seed(1) = 12345
  Call Random_Seed

  npackets = 1000
  rmax = 5.
  

do k = 1, 500
  taumax = 0.1*Real(k)
  j=0.
  do i=1, npackets
     
    call emit_packet(x,y,z,nx,ny,nz)
    r = SQRT(x**2+y**2+z**2)
    do while(r<rmax)!packet is in slab
      !update packet position, x,y,z,r
      call Random_Number(ran)
      L = -log(ran)*(rmax/taumax)
      x = x + L*nx
      y = y + L*ny
      z = z + L*nz
      r = SQRT(x**2 + y**2 + z**2)
      
      if (r > rmax) exit
      
      call random_Number(ran)
      if(ran < albedo) then
        call scatter(nx,ny,nz)
      else
        exit
      end if
    enddo
  end do

   !print*, taumax, (j/(npackets)) 
   open(1,file=name,position='append',iostat=ok)
   write(1,*) taumax,(j/(npackets))
   close(1)
  
end do
call CPU_TIME(t)
  print*, t
end program Isotropic_Scatt_V3

!-----------------------------------------------------------------


subroutine emit_packet(x,y,z,nx,ny,nz)
  implicit none

  real*8, intent(out) :: x,y,z,nx,ny,nz
    real*8 :: Pi, ran, phi, theta


  Pi = 3.1415926535897932385


  call Random_Number(ran)
  theta = acos(2.*ran - 1.)
  
  call Random_Number(ran)
  phi = 2.*Pi*ran

  nx = Sin(theta)*Cos(phi)
  ny = Sin(theta)*Sin(phi)
  nz = Cos(theta)

  x = 0
  y = 0
  z = 0

  !return

end subroutine emit_packet

!-------------------------------------------------------------------------

subroutine scatter(nx,ny,nz)
  use global
  implicit none
  real*8, intent(inout) :: nx,ny,nz
  real*8 :: Pi, ran, phi, theta


  Pi = 3.1415926535897932385


  call Random_Number(ran)
  theta = acos(2.*ran - 1.)
  
  call Random_Number(ran)
  phi = 2.*Pi*ran

  nx = Sin(theta)*Cos(phi)
  ny = Sin(theta)*Sin(phi)
  nz = Cos(theta)
 
  j=j+1.
  
end subroutine scatter

     

      

  

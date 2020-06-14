module bins
  implicit none
  
  real*8 :: thetaarray(50)
  real*8 :: error(50)

end module bins

!------------------------------------------------------

program Monte_Carlo_Check
  use bins
  
  implicit none

  integer :: i, npackets, iseed,ok,k
  real*8 :: x,y,z,nx,ny,nz,L,taumax,albedo, ran, thetamax,thetamin
  real*8 :: xmax, ymax, zmax, t1, pi
  real*8 :: theta
  integer, Dimension (1) :: Seed
  logical :: inx, iny, inz
  character (len=80) :: name
  integer :: thetabin

  pi=4.*Atan(1.)
  
  name = "Angles.dat"
  open(1,file=name,status='replace',iostat=ok)
  close(1)
  

  taumax=100.0 !n*cross section* zmax
  albedo=1.0 !always a scattering event here
  
  
  Seed(1) = 12345
  Call Random_Seed
  
  npackets=10**10
  
  
  thetamax=0.00075
  
  
  do i=1, npackets
  
    
      call scatter(theta)
      thetabin = Int((theta/thetamax)*50) + 1
      If (thetabin < 50) then
        thetaarray(thetabin) = thetaarray(thetabin) + 1.
      endif

    
  end do
  

  
  do k = 1, 50
     error(k) = sqrt(thetaarray(k)/(2.d0*Pi*(k-0.5d0)* &
      & (thetamax/50.d0)))
        thetaarray(k) = thetaarray(k)/(2.*Pi*(k-0.5d0)* &
     &  (thetamax/50.d0))
  enddo
  open(1,file=name,position='append',iostat=ok)
  do k=1,50
  write(1,*) thetaarray(k), error(k)
  enddo
  close(1)
  
  call cpu_time(t1)
  write (*,fmt="(F6.2)",advance="no") t1
    
  
end program Monte_Carlo_Check

!----------------------------------------------------

subroutine emit_packet(x,y,z,nx,ny,nz)
  implicit none


  real*8, intent(out) :: x,y,z,nx,ny,nz

  x = 5.
  y = 5.
  z = 10. !10cm cube

  nx = 0.
  ny = 0.
  nz = -1. !Set for theta is Pi and phi is 0

return

end subroutine emit_packet

!--------------------------------------------------

subroutine scatter(theta)
  
  implicit none
  
  real*8 :: Pi, ran, phi,a, thetam
  real*8 :: T
  real*8, intent(out) :: theta
  a = 2.33333333*10d-9
  Pi = 3.14159
  thetam = 7.536d-5
  

  call Random_Number(ran)
  theta = 2.*asin(SQRT(2./((2./Sin(thetam/2.)**2)-ran*(2./(Sin(thetam/2.)**2))+2*ran)))

    

end subroutine scatter

     



  

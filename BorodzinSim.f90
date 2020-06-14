!------------------------------------------------------
!*********    Borodzin Simulation Program     ********!
!*********      for muons through iron        ********!
!*********         10 March 2020              ********!
!*********        Sebastian Taylor            ********!
!*********Advisors Antje Kohnle and Kenny Wood********!
!------------------------------------------------------
module bins
  implicit none
  real*8, Allocatable :: xval(:)
  real*8, allocatable :: yval(:)
  real*8,allocatable :: img(:,:)
  real*8, Allocatable :: thetaarray(:)
  real*8, Allocatable :: error(:)
  real*8, Allocatable :: Counts(:)
end module bins
module global
  implicit none
  integer :: j
end module global
!------------------------------------------------------
program main
  use bins
  use global
  implicit none
  integer :: i, npackets, iseed,ok,thetabin,k
  real*8 :: x,y,z,L,taumax,albedo, ran, d,Pi
  real*8 :: xmax, ymax, zmax, t1,xmin,ymin
  real*16 :: theta,nx,ny,nz,thetamax,thetamin,phi
  real*16 :: xim,yim,r,nxim,thetafin, phifin
  integer, Dimension (1) :: Seed
  logical :: inx, iny, inz
  character (len=80) :: name, name2, name3, name4
  !Material Dependent Paramters
!***********************************************!
  taumax=26389.4d0 !n*cross section* zmax       
  thetamin = 2.49932d-5                                                                         
  thetamax = 5.d3*thetamin                      
  npackets = 10000                              
  d = 5.d1 !Binning Parameter 
  zmax = 1.d-1   
!***********************************************!
  albedo=1.0 !always a scattering event here
  !Inititialize Random number generator
  Seed(1) = 12345
  Call Random_Seed
  !Dimensions of Material
  xmax = 3.
  ymax = 3.
  xmin = -3.
  ymin= -3.
  j=0
  Pi = 4*ATAN(1.d0)
  call initit(npackets, d)
  name = "RutherfordImg.dat"
  open(1,file=name,status='replace',iostat=ok)
  name2 = "Borodzin.dat"
  name3 = "ErrorHist.dat"
  name4 = "Counts.dat"
  !fire muons
  do i=1, npackets
    print *,  i
    call emit_packet(x,y,z,nx,ny,nz,zmax)
    !Slab dimensions
    inx = ((x >= xmin) .and. (x <= xmax))
    iny = ((y >= xmin) .and. (y <= ymax))
    inz = ((z >= 0) .and. (z <= zmax))
    do while(inx .and. iny .and. inz)!packet is in slab
      call Random_Number(ran)
      L = -log(ran)*zmax/taumax
      x = x + L*nx
      y = y + L*ny
      z = z + L*nz !update packet position, x,y,z
      r = z - L*nz !Projection vector
      if ((z < 0) .or. (z > zmax) .or. (x < xmin) .or. (x > xmax) .or. (y < ymin) .or. (y > ymax)) exit
      call Random_Number(ran)
      if(ran<albedo) then
        call scatter(nx,ny,nz, theta,phi)
      else
        exit
      end if
    enddo
    if (z <= 0) then
     xim = (x - L*nx) + r*nx
     yim = (y - L*ny) + r*ny
     thetafin = Pi-ACos(nz)
     phifin = Atan(ny/nx)
     !Binning Function
     thetabin = Int((thetafin/thetamax)*d) + 1
     If (thetabin < d) then
       thetaarray(thetabin) = thetaarray(thetabin) + 1.
     else
       thetaarray(70) = thetaarray(70) + 1.
     endif
     xval(i) = xim
     yval(i) = yim
     img = Reshape((/ (/xval/),(/yval/)/)  ,(/npackets,2/) ) 
     !Save img coordinates to file
     write(1,*) img(i,1),img(i,2)
    endif
  end do
  !Close image file
  close(1)
  do k=1, Int(d)
    error(k) = Sqrt(thetaarray(k))/(2.d0*Pi*(k-0.5d0)*(thetamax/d))
    thetaarray(k) = thetaarray(k)/(2.d0*Pi*(k-0.5d0)*(thetamax/d))
    Counts(k) = thetaarray(k)
  enddo 
  !Write angle data
  open(2,file=name2,status='replace',iostat=ok)
  write(2,*) thetaarray
  close(2)
  open(3,file=name3,status='replace',iostat=ok)
  write(3,*) error
  close(3)
  open(4,file=name4,status='replace',iostat=ok)
  write(4,*) Counts
  close(4)
  print*, "the average number of scatterings is,", Real(j)/npackets
  !print the time to run the progra,
  call cpu_time(t1)
  write (*,fmt="(F6.2)",advance="no") t1
end program main
!----------------------------------------------------
subroutine emit_packet(x,y,z,nx,ny,nz,zmax)
  implicit none
  real*8, intent(in) :: zmax
  real*8, intent(out) :: x,y,z
  real*16, intent(out) :: nx,ny,nz
  !starting positions
  x = 0.
  y = 0.
  z = zmax 
  !starting angles
  nx = 0.
  ny = 0.
  nz = -1. !Set for theta is Pi and phi is 0
end subroutine emit_packet
!--------------------------------------------------
subroutine scatter(nx,ny,nz,theta,phi)
  use global
  implicit none
  real*16, intent(inout) :: nx,ny,nz
  real*8 :: Pi, ran,a, thetam
  real*16 :: T,nx1,ny1,nz1
  real*16, intent(out) :: theta,phi
  real*16 :: l
  l=9.9999999d-1
  Pi = 4*ATAN(1.d0)
  thetam = 7.536d-5
  a = (2./Sin(thetam/2.)**2)
  call Random_Number(ran)
  theta = 2.d0*asin(SQRT(2.d0/(a-ran*a + 2.d0*ran)))
  call Random_Number(ran)
  phi = 2.d0*Pi*ran
  !If the event is close enough hto the normal use this
  if (Abs(nz) > l) then
    nx = Sin(theta)*Cos(phi)
    ny = Sin(theta)*Sin(phi)
    nz = Cos(theta)*nz/Abs(nz)
    !print*, "inner", l, nx,ny,nz
  !Otherwise Change of Frame of Reference
  else
    T = SQRT(1.d0 - nz**2)
    nx1 = (Sin(theta)*(nx*nz*Cos(phi) - ny*Sin(phi))/T) + nx*Cos(theta)
    ny1 = (Sin(theta)*(ny*nz*Cos(phi) + nx*Sin(phi))/T) + ny*Cos(theta)
    nz1 = -Sin(theta)*Cos(phi)*T + nz*Cos(theta)
    !Resassign for nx,ny,nz output
    nx = nx1
    ny = ny1
    nz = nz1
  endif  
  !Number of scatter calls
  j = j+1
end subroutine scatter
!-------------------------------------------------------------------
subroutine initit(npackets,d)
  use bins
  implicit none
  integer, intent(in) :: npackets
  real*8, intent(in) :: d
  Allocate(thetaarray(INt(d)))
  Allocate(xval(npackets))
  Allocate(yval(npackets))
  Allocate(img(npackets,2))
  Allocate(error(Int(d)))
  Allocate(Counts(npackets))
end subroutine initit



  

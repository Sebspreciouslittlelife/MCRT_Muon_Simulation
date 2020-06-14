!------------------------------------------------------
!*********    Borodzin Simulation Program     ********!
!*********      for muons through iron        ********!
!*********         10 March 2020              ********!
!*********           150002081                ********!
!*********Advisors Antje Kohnle and Kenny Wood********!
!------------------------------------------------------
!Bins is a module that holds the lists as global 
!variables
module bins
  implicit none
  real*8, Allocatable :: xval(:)
  real*8, allocatable :: yval(:)
  real*8,allocatable :: img(:,:)
  real*8, Allocatable :: thetaarray(:)
  real*8, Allocatable :: error(:)
  integer, Allocatable :: Counts(:)
end module bins
!Global stores the number of scatterngs as a
!global variable
module global
  implicit none
  integer :: j
end module global
!Constants holds nature constants and 
!the material specificconstants
!as global variables
module Constants
  implicit none
  real*8 :: re,me,Za,a0,hbar,Pi,c
  real*8 :: C1,C2,C3
end module Constants
!------------------------------------------------------
program main
  use bins
  use global
  use Constants
  implicit none
  integer :: i, npackets, iseed, ok, thetabin, k,m
  real*8 :: x, y, z, L, taumax, albedo, ran, d
  real*8 :: xmax, ymax, zmax, t1, t2, xmin, ymin
  real*8 :: theta, nx, ny, nz, thetamax, thetamin, phi
  real*8 :: xim, yim, r, nxim, thetafin, phifin
  integer, Dimension (1) :: Seed
  logical :: inx, iny, inz
  character (len=80) :: name, name2, name3, name4
  call cpu_time(t1)
  !Material Dependent Paramters
!***********************************************!                          
  thetamax = 10.0d-2                             
  npackets=10000                               
  d = 5.d1 !Binning Parameter                
  zmax = 1.d-1 
!***********************************************!
  albedo=1.0 !always a scattering event here
  !Inititialize Random number generator
  Seed(1) = 12345
  Call Random_Seed
  !Dimensions of Material
  xmax = 1.
  ymax = 1.
  xmin = -1.
  ymin= -1.
  !Set scatters to 0
  j=0
  !Initialize the arrays and constants
  call initit(npackets, d)
  !calculate the total optical depth
  taumax = zmax*8.4d28*(4.d0*Pi*C1)/(C2**2 + C2*C3)!n*cross section* zmax
  print*, "taumax is",taumax
  !fire muons
  do i=1, npackets
    print *, "Calling new packet", i
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
      !Exit loop if muon leaves the slab
      if ((z < 0) .or. (z > zmax) .or. (x < xmin) .or. (x > xmax) .or. &
      & (y < ymin) .or. (y > ymax)) exit
      !Otherwise scatter after traveling along random walk
      call Random_Number(ran)
      if(ran<albedo) then
        call scatter(nx,ny,nz,theta,phi)
      else
        exit
      end if
    enddo
    !Project muons onto the bottom of the slab
    if (z <= 0) then
     xim = (x - L*nx) + r*nx
     yim = (y - L*ny) + r*ny
     !Bin the final angles
     thetafin = Pi-ACos(nz)
     phifin = Atan(ny/nx)
     thetabin=Int((thetafin/thetamax)*d) + 1
     If (thetabin < d) then
       thetaarray(thetabin) = thetaarray(thetabin) + 1.
     endif
     !Save the x- and y- values in separate lists
     xval(i) = xim
     yval(i) = yim
    endif
  end do
  !Save img coordinates to file
  name = "ScreenedImg.dat"
  name2 = "ScreenedAngles.dat"
  name3 = "ErrorHistScreen.dat"
  name4 = "CountsError.dat"
  img = Reshape((/ (/xval/),(/yval/)/)  ,(/npackets,2/) ) 
  open(1,file=name,status='replace',iostat=ok)
  do m=1, npackets
    write(1,*) img(m,1),img(m,2)
  enddo
  close(1)
  do k=1, Int(d)
    error(k) = Sqrt(thetaarray(k)/(2.d0*Pi*(k-0.5d0)*(thetamax/d)))
    Counts(k) = thetaarray(k)
    thetaarray(k) = thetaarray(k)/(2.*Pi*(k-0.5)*(thetamax/d))
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
  print*, "the average number of scatterings is,", Real(j)/Real(npackets)
  !print the time to run the progra,
  call cpu_time(t2)
  write (*,fmt="(F8.2)",advance="no") t2-t1
end program main
!----------------------------------------------------
!Set up the initial conditions of the muon
subroutine emit_packet(x,y,z,nx,ny,nz,zmax)
  use Constants
  implicit none
  real*8, intent(in) :: zmax
  real*8, intent(out) :: x,y,z
  real*8, intent(out) :: nx,ny,nz
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
!screened Rutherford scattering
subroutine scatter(nx,ny,nz,theta,phi)
  use global
  use Constants
  implicit none
  real*8, intent(inout) :: nx,ny,nz
  real*8 :: ran
  real*8 :: T,nx1,ny1,nz1
  real*8, intent(out) :: theta,phi
  real*8 :: l
  l=9.9999999d-1
  call Random_Number(ran)
  theta = 2.d0*Asin(Sqrt((C2*ran)/(C2 + C3 - (C3*ran))))
  call Random_Number(ran)
  phi = 2.d0*Pi*ran
  !If the event is close enough to the normal use this
  if (DBLE(Abs(nz)) > l) then
    nx = Sin(theta)*Cos(phi)
    ny = Sin(theta)*Sin(phi)
    nz = Cos(theta)*nz/Abs(nz)
  !Otherwise Change of Frame of Reference
  else
    T = SQRT(1. - nz**2)
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
!Allocate array sizes and assign constants
subroutine initit(npackets,d)
  use bins
  use Constants
  implicit none
  integer, intent(in) :: npackets
  real*8, intent(in) :: d
  real*8 :: lambda,p
  Allocate(thetaarray(INt(d)))
  Allocate(xval(npackets))
  Allocate(yval(npackets))
  Allocate(img(npackets,2))
  Allocate(error(Int(d)))
  Allocate(Counts(npackets))
  !Assigns values
  hbar = 1.055d-34 !Reduced Planck Constant
  Za = 2.6d1 !Atomic number of slab
  a0 = 5.29d-11 !Bohr Radius
  re = 2.82d-15 !Electron radius
  me = 9.11d-31 !Electron mass
  Pi = 4*ATAN(1.d0)
  c = 3.d8 !Speed of light
  !Define Constants for Muon Energy of 3 GeV
  p = Sqrt(((3.d9*1.6d-19)/c)**2-((1.057d8*1.6d-19)/c)**2)
  lambda = Za**(1.d0/3.d0)/a0
  C1 = 4.d0*(p**2)*(Za**2)*(1.d0/137.d0)**2*(hbar**2)
  C2 = (hbar**2)*(lambda**2)
  C3 = 4.d0*(p**2)
end subroutine initit



  

      subroutine ScreenedPot(nxp,nyp,nzp,theta,phi,pi,
     +           twopi,iseed,C1,C2,C3,Z)
      
      implicit none
      
      integer iseed
      real*8 nxp,nyp,nzp,phi,theta
      real*8 pi,twopi
      real*8 C1,C2,C3
      
      integer Z
      real*16 T,nx1,ny1,nz1
      real*8 hbar,a0,re,me,c,p,lambda
      real*16 l
      real ran2
      real*8 random
      
      
      hbar = 1.055d-34
      a0 = 5.29d-11
      re = 2.82d-15
      me = 9.11d-31
      c = 3.d8
      
      p = sqrt(((3.d9*1.6d-19)/c)**2-((1.057d8*1.6d-19)/c)**2)
      lambda = Z**(1.d0/3.d0)/a0
      C1 = 4.d0*(p**2)*(Z**2)*(1.d0/137.d0)**2*(hbar**2)
      C2 = (hbar**2)*(lambda**2)
      C3 = 4.d0*(p**2)

      l=9.9999999d-1
      random=1.d0*ran2(iseed)
      theta=2.d0*asin(sqrt((C2*random)/
     + (C2+C3-(C3*random))))
      phi=2.d0*pi*ran2(iseed)

      if(dble(abs(nzp)).gt.l) then
        nxp=sin(theta)*cos(phi)
        nyp=sin(theta)*sin(phi)
        nzp=cos(theta)*nzp/abs(nzp)
        
      else
        T=sqrt(1.-nzp**2)
        nx1=(sin(theta)*(nxp*nzp*cos(phi)-
     +  nyp*sin(phi))/T)+nxp*cos(theta)
        ny1=(sin(theta)*(nyp*nzp*cos(phi)+
     +  nxp*sin(phi))/T)+nyp*cos(theta)
        nz1=(-sin(theta)*cos(phi)*T+nzp*cos(theta))
        
        nxp=nx1
        nyp=ny1
        nzp=nz1
      endif
      
      return
      end

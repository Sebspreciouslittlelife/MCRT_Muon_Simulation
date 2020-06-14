      subroutine sourceph(xp,yp,zp,nxp,nyp,nzp,
     +                    sint,cost,sinp,cosp,phi,fi,fq,fu,fv,
     +                    xmax,ymax,zmax,twopi,
     +                    xcell,ycell,zcell,nxg,nyg,nzg,iseed,C1,C2,C3)

      implicit none

      include 'photon.txt'

      integer xcell,ycell,zcell,nxg,nyg,nzg,iseed
      real*8 xmax,ymax,zmax,twopi
      real*8 hbar,Za,a0,re,me,c,p,lambda
      real*8 C1,C2,C3
      real ran2
      

c***** emit photon randomly from origin
      xp = 9.d-1*xmax*(1.d0-2.d0*ran2(iseed))
      yp = 9.d-1*ymax*(1.d0-2.d0*ran2(iseed))
            
    
      zp=0.999999*zmax

      cost=2.*ran2(iseed)-1.
      cost=-1.
      sint=(1.-cost*cost)
      if(sint.le.0.)then
        sint=0.
      else
        sint=sqrt(sint)
      endif

      phi=twopi*ran2(iseed)
      cosp=cos(phi)
      sinp=sin(phi)
      

c***** Set photon direction cosines for direction of travel *********
      nxp=sint*cosp  
      nyp=sint*sinp
      nzp=cost

c***** Set Stokes fluxes ********************************************
      fi=1.
      fq=0.
      fu=0.
      fv=0.

c*************** Linear Grid *************************
      xcell=int(nxg*(xp+xmax)/(2.*xmax))+1
      ycell=int(nyg*(yp+ymax)/(2.*ymax))+1
      zcell=int(nzg*(zp+zmax)/(2.*zmax))+1
c*****************************************************
      hbar = 1.055d-34
      Za = 1.2d1
      a0 = 5.29d-11
      re = 2.82d-15
      me = 9.11d-31
      c = 3.d8
      
      p = Sqrt(((3.d9*1.6d-19)/c)**2-((1.057d8*1.6d-19)/c)**2)
      lambda = Za**(1.d0/3.d0)/a0
      C1 = 4.d0*(p**2)*(Za**2)*(1.d0/137.d0)**2*(hbar**2)
      C2 = (hbar**2)*(lambda**2)
      C3 = 4.d0*(p**2)

      return
      end


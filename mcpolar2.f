      program mcpolar
      implicit none
      include 'grid.txt'
      include 'photon.txt'
c***** Parameter declarations ****************************************
      integer nphotons,iseed,j,xcell,ycell,zcell,tflag,k
      real*8 nscatt
      real*8 kappa,albedo,hgg,pl,pc,sc,xmax,ymax,zmax,theta
      real*8 pi,twopi,fourpi,g2,delta,thetamax
      real*8 xexit,yexit,zexit,thetafin,phifin
      real*8 C1,C2,C3
      real*8 thix,thiy,thiz
      real*8 pix,piy,piz
      real*8 mu,lam
      integer voxx,voxy,voxz,k1,k2
      integer voxx1,voxy1,voxz1,errs,errs2
      real*8 xexitlist(100000),yexitlist(100000)
      real*8 zexitlist(100000),thetaarray(50)
      integer thetabin, error(50)
      character (len=80) name1, name2,name3
      real ran2
      !File names
      name1='gridci.dat'
      name2='3DAnglesci.dat'
      name3='imgci.dat'
      !Create files
      open(11,file=name1,status='replace')
      close(11)
      open(12,file=name2,status='replace')
      close(12)
      open(13,file=name3,status='replace')
      close(13)
c**** Read in parameters from the file input.params
      open(10,file='input.params',status='old')
          read(10,*) nphotons
          read(10,*) iseed
          read(10,*) kappa
          read(10,*) albedo
          read(10,*) xmax
          read(10,*) ymax
          read(10,*) zmax
          close(10)
c***** Set up constants, pi and 2*pi  ********************************
      pi=4.*atan(1.)
      twopi=2.*pi
      fourpi=4.*pi
      thetamax = 10.0d-2 
      iseed=-abs(iseed)  ! Random number seed must be negative for ran2
c**** Initialize arrays to zero *************************************
      call iarray(xface,yface,zface,rhokap,vox)
c***** Set up density grid *******************************************
      call gridset(xface,yface,zface,rhokap,xmax,ymax,zmax,
     +            kappa,znumber)
c***** Set small distance for use in optical depth integration routines 
c***** for roundoff effects when crossing cell walls
      delta=1.e-7*(2.*xmax/nxg)
      !Initialize error variables
      errs=0
      errs2=0
c**** Loop over nph photons from each source *************************
        nscatt=0 !scatter counter
        do j=1,nphotons
          !Print every 10 muons
          if(mod(j,10).eq.0)then
             print *, j,' scattered muons completed'
          end if
c***** Release photon from point source *******************************
          call sourceph(xp,yp,zp,nxp,nyp,nzp,sint,cost,sinp,cosp,phi,
     +         fi,fq,fu,fv,xmax,ymax,zmax,twopi,
     +         xcell,ycell,zcell,nxg,nyg,nzg,iseed,C1,C2,C3)
          !Save the incident coordinates and angles
          pix=xp
          piy=yp
          piz=zp
          thix=0.
          thiy=0.
          thiz=-1.
c****** Find scattering location (Designed by K Wood)
          call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +         xface,yface,zface,rhokap,xcell,ycell,zcell,
     +         tflag,iseed,delta,xexit,yexit,zexit)
c******** Photon scatters in grid until it exits (tflag=1) 
c          tflag=0
          dowhile(tflag.eq.0)
                 if( (ran2(iseed).lt.albedo) ) then
!Scatter photon into new direction using screened Rutherford scattering
                   call ScreenedPot(nxp,nyp,nzp,theta,phi,pi,
     +                  twopi,iseed,C1,C2,C3,znumber(xcell,ycell,zcell))
                   nscatt=nscatt+1 !update number of scatters
                else
                   goto 100
                endif
c************ Find next scattering location
              call tauint2(xp,yp,zp,nxp,nyp,nzp,xmax,ymax,zmax,
     +         xface,yface,zface,rhokap,xcell,ycell,zcell,
     +         tflag,iseed,delta,xexit,yexit,zexit)
          end do
100      continue
        !Save the exit coordinates to lists
        xexitlist(j)=xexit
        yexitlist(j)=yexit
        zexitlist(j)=zexit
        !Calculate the exit angle
        thetafin = pi-acos(nzp)
!         Calculate the PoCA
        mu=(nxp*(pix-xexit)+nyp*(piy-yexit))/(nxp**2+nyp**2)
        lam=zexit-piz+(mu*nzp)
        !Assign PoCA vocels for incident vector
        voxx=Int((pix)*50.d0+2.5)+1
        voxy=Int((piy)*50.d0+2.5)+1
        voxz=Int((piz+lam)*50.d0+2.5)+1
        !Assign PoCA vocels for exit vector
        voxx1=Int((xexit+mu*nxp)*50.d0+2.5)+1
        voxy1=Int((yexit+mu*nyp)*50.d0+2.5)+1
        voxz1=Int((zexit+mu*nzp)*50.d0+2.5)+1
        !Record large errors in voxel assignation
        if((abs(voxx-voxx1).gt.0.d0).or.(abs(voxy-voxy1)
     +     .gt.0.d0).or.(abs(voxz-voxz1).gt.0.d0))then
            errs=errs+1
        endif
        !Assign signal to PoCA voxel and 0 to 
        !other candidiates
        if((voxz.gt.0).and.(voxz.lt.6))then
          do k=1,voxz-1
            counterv(voxx,voxy,k)=counterv(voxx,voxy,k)+1
          enddo
!           Assuming output angle is small
          do k=voxz+1,5
            counterv(voxx,voxy,k)=counterv(voxx,voxy,k)+1
          enddo
!           Signal
          vox(voxx,voxy,voxz)=vox(voxx,voxy,voxz)+thetafin**2
          counterv(voxx,voxy,voxz)=counterv(voxx,voxy,voxz)+1
        else
        !Otherwise record if PoCA is outside the slab
          errs2=errs2+1
        endif
        !Bin output angles  
        thetabin = Int((thetafin/thetamax)*50) + 1
        If (thetabin < 50) then
          thetaarray(thetabin) = thetaarray(thetabin) + 1.
        endif
        end do      ! end loop over nph photons
     !Write to file
      open(11,file=name1,position='append',status='old')
      do k=1,nphotons
        write(11,*) xexitlist(k),yexitlist(k),zexitlist(k)
      enddo
      close(11)
      do k=1, 50
        error(k) = sqrt(thetaarray(k)/(2.d0*Pi*(k-0.5d0)*
     +  (thetamax/50.d0)))
        thetaarray(k) = thetaarray(k)/(2.*Pi*(k-0.5d0)*
     +  (thetamax/50.d0))
      enddo 
      open(12,file=name2,position='append',status='old')
      do k=1,50
        write(12,*) thetaarray(k), error(k)
      enddo
      open(13,file=name3,position='append',status='old')
      do k=1,5
        do k1=1,5
          do k2=1,5
            if(counterv(k,k1,k2).eq.0)then
            write(13,*) (k-3)/50.,(k1-3)/50.,
     +              (k2-3)/50.,0
            else
            write(13,*) (k-3)/50.,(k1-3)/50.,
     +              (k2-3)/50.,Real(vox(k,k1,k2)/counterv(k,k1,k2))
            endif
          enddo
        enddo
      enddo
      close(13)
      print*,'Avereage number of scatterings = ',(nscatt/nphotons)
      print*,'Avereage number of misaligned voxels = ',
     +    (real(errs)/nphotons)
      print*,'Avereage number of lost muons is = ', errs2
      stop
      end

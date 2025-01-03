!*************************
! PROGRAM TriffyBloc
! A QuarterSpace with symmetry on three sides
! and absorbing boundaries on the three others.
! for simulation of an infinite medium
! *************************
! Stefan Nielsen September 2000
!
! PILOT
!
include 'triffy.dec'            !declaration file
open (unit=17,file='check',form='formatted')
it=0
call zero                    !initialize variables & tools
call inizia
it=0
do while(it.le.ntime)  !start time iteration
 it = it+1
!print*,'in loop...', it
 call stress                !calculate the stresses
 call drive_PML_stress
 call imposed_source
 !call rupture
 call save               !visualize and/or store variables 
 call extrapolate           !extrapolate in time
 call drive_PML_velo
end do                      !end of time iteration
close (12)
close (26)
stop 
end
! **********************************
! SUBROUTINES:
! **********************************
       subroutine zero
       include 'triffy.dec'
       character*80 fin,fmodel
       character*16 nam
       integer getlen
       real muAll,lamAll,rhoAll
! FIXED PARAMETERS:
       nos1=100
       savefields=1        ! 1 if want to save 
       savetraces=1        ! 1 if want to save 
       rn=9./8.            ! interpolation coefficients
       rnn=-1./24.         !  "  "
       pi=acos(-1.)
       is1=0               ! markers initialization
       is2=0               ! 
       nxs=20
       nys=20
       nzs=nz/2
! MEDIUM PARAMETERS - NO INPUT FILE:
       muAll=25000.0;lamAll=25000.0;rhoAll=.0028 !realistiv values but all divided by 1e6
       do k=1,nz
       do j=1,ny
       do i=1,nx
        mu(i,j,k)  = muAll
        lam(i,j,k) = lamAll
        rho(i,j,k) = rhoAll
       enddo
       enddo
       enddo
      fmodel='model'
      il=getlen(fmodel)
!      call readfiles(fmodel,il)
!      call MediumPMLparameters
      call iniPMLini(muAll,lamAll,rhoAll)
! Fault parameters- NO INPUT FILE:
      nzf=nz/2 ! this also serves as mid-plane for 2D/3d
      nxf1=5 
      nxf2=nx-4
!     nyf1=5      ! CAUTION : this is defined after free surface
!     nyf2=ny-4   ! CAUTION : this "    .   .   .   .   .      "
!     peak=1.
!     mu_s=0.
!     mu_d=0.
!     delta=2.
!     vmax=0.
!     dissipation=0.00
! READ INPUT PARAMETERS:
       open (11,file='fpar',form='formatted',status='old')
       read (11,*)
       read (11,*) dx    
       read (11,*) dt
       read (11,*) ntime
       read (11,*)
       read (11,*) tw  
       read (11,*) vforce
       read (11,*) peak 
       read (11,*) stremin,stremax,struncate
       read (11,*) delta, deltamax
       read (11,*) vmax
       read (11,*) dissipation
       read (11,*) mu_s 
       read (11,*) mu_d 
       read (11,*) FricType
       read (11,*) tau2
       read (11,*) rInitRadius
       read (11,*)
       read (11,*) rxs,rys
       read (11,*)
       nr=1
       do nr=1,6
        read (11,*) rxr,ryr
        nxr(nr)=rxr/dx-npm+3
        nyr(nr)=ryr/dx-npm+3
       enddo
       close (11)
       InitRadius=nint(rInitRadius/dx)
       open (unit=24,file='filpar.dec',form='formatted')
       write (24,'(a19,G8.4,a1)') '      parameter(dt=',dt,')'
       rdis=0.
       broken=0
       incrack=1
       strini23=0.
       nxs=int(rxs/dx)
       nys=int(rys/dx)
! read input prestress file:
       open (11,file='uy_160_80',form='formatted')
        do j=nyf1,nyf2
         read(11,*) (strini23(i,j),i=nxf1,nxf2)
        enddo
       rmaxin=0.
       rminin=strini23(10,10)
       do j=nyf1,nyf2
        do i=nxf1,nxf2
         if (strini23(i,j).gt.rmaxin) rmaxin=strini23(i,j)
         if (strini23(i,j).lt.rminin) rminin=strini23(i,j)
        enddo
       enddo
       strini23=stremin+stremax*(strini23-rminin)/(rmaxin-rminin)
! initial asperity of InitRadius
! and truncate stress above threshold
      do j=1,ny
      do i=1,nx
       if (strini23(i,j).gt.struncate) strini23(i,j)=struncate
       if ((nxs-i)**2+(nys-j)**2.lt.(InitRadius)**2) then
        strini23(i,j)=peak*1.05
       endif
      enddo
      enddo
!cxcx
      open (34,file='strini',form='formatted')
      do j=1,ny;write(34,'(1000F14.7)')(strini23(i,j),i=1,nx);enddo;close(3)
!cxcx
! the dirac functions remove forth order in layers nzf+1 and nzf-1 
! form some terms, eventually add non zero dissipation:
      dirac=0.
      dirac(nzf)=1. 
! EFFECTIVE MEDIUM PARAMETERS: computed in iterations to save mem.
! OUTPUT FILES, OPEN:
       open (unit=15,file='source',form='formatted')
       open (unit=26, file='iter',form='formatted')
!
       numsta=7
       do i=1,numsta !open single trace files for output:
        write (nam,'(a,I3.3)') './RES/tracex',i
        open (60+i,file=nam,form='formatted')
        write (nam,'(a,I3.3)') './RES/tracey',i
        open (60+numsta+i,file=nam,form='formatted')
       enddo
!
       dtx=dt/dx
! Initialisation of the PML layers, we hypothise that the medium is not 
! heterogeneous in the corners
       call InitDist
       return
       end
!************************************
       subroutine stress
       include 'triffy.dec'      
       real fsa
       fsa=3.
! Sij = lam*dij*Ekk + 2*mu*Eij
! 4th order interpolation:
       do k=3,nz-2  ! Mod
       do j=3,ny-2
       do i=3,nx-2  ! Mod
        s11(i,j,k)=  s11(i,j,k) + (dt/dx) * ( &
        (lam(i,j,k)+2*mu(i,j,k)) * ( &
         rn*(v1(i,j,k)-v1(i-1,j,k))+rnn*(v1(i+1,j,k)-v1(i-2,j,k)) ) &
        +lam(i,j,k) * ( &
         rn*(v2(i,j+1,k)-v2(i,j,k))+rnn*(v2(i,j+2,k)-v2(i,j-1,k))+  &
         rn*(v3(i,j,k)-v3(i,j,k-1))+rnn*(v3(i,j,k+1)-v3(i,j,k-2)))) 

        s22(i,j,k)= s22(i,j,k) + (dt/dx) * ( &
        (lam(i,j,k)+2*mu(i,j,k)) * ( &
         rn*(v2(i,j+1,k)-v2(i,j,k))+rnn*(v2(i,j+2,k)-v2(i,j-1,k))) &
        +lam(i,j,k) * ( &
         rn*(v1(i,j,k)-v1(i-1,j,k))+rnn*(v1(i+1,j,k)-v1(i-2,j,k))+  &
         rn*(v3(i,j,k)-v3(i,j,k-1))+rnn*(v3(i,j,k+1)-v3(i,j,k-2)))) 

        s33(i,j,k)= s33(i,j,k) + (dt/dx) * ( &
        (lam(i,j,k)+2*mu(i,j,k)) * ( &
         rn*(v3(i,j,k)-v3(i,j,k-1))+rnn*(v3(i,j,k+1)-v3(i,j,k-2))) &
        +lam(i,j,k) * ( &
         rn*(v1(i,j,k)-v1(i-1,j,k))+rnn*(v1(i+1,j,k)-v1(i-2,j,k))+  &
         rn*(v2(i,j+1,k)-v2(i,j,k))+rnn*(v2(i,j+2,k)-v2(i,j-1,k))))

        s12(i,j,k)=  s12(i,j,k) + &
         (4./(1/mu(i,j,k)+ 1/mu(i+1,j,k)+ &
         1/mu(i+1,j-1,k)+ 1/mu(i,j-1,k)))*(dt/dx)*( &
         rn*(v1(i,j,k)-v1(i,j-1,k))+rnn*(v1(i,j+1,k)-v1(i,j-2,k))+ &
         rn*(v2(i+1,j,k)-v2(i,j,k))+rnn*(v2(i+2,j,k)-v2(i-1,j,k)))

        s23(i,j,k)=  s23(i,j,k) + mu(i,j,k)*(dt/dx)*( &
         rn*(v3(i,j,k)-v3(i,j-1,k))+rnn*(v3(i,j+1,k)-v3(i,j-2,k))+ &
         rn*(v2(i,j,k+1)-v2(i,j,k))+rnn*(v2(i,j,k+2)-v2(i,j,k-1)))
  
        s13(i,j,k)=s13(i,j,k)+mu(i,j,k)*(dt/dx)*( &
         (rn*(v1(i,j,k+1)-v1(i,j,k))+rnn*(v1(i,j,k+2)-v1(i,j,k-1)) &
         +rn*(v3(i+1,j,k)-v3(i,j,k))+rnn*(v3(i+2,j,k)-v3(i-1,j,k))))


        enddo
        enddo
        enddo
!      print *,s11(5,5,5),s13TXTY(5,5,5)

! free surface by "vacuum":
! TODO: Change with symmetry/mirror condition for better match
      do k=1,nz
       do j=1,4
        do i=1,nx
         s11(i,j,k)=0.; s22(i,j,k)=0.; s33(i,j,k)=0.
         s13(i,j,k)=0.; s12(i,j,k)=0.; s23(i,j,k)=0.
        enddo
       enddo
      enddo

       return
       end
! *******************************************
       subroutine extrapolate
       include 'triffy.dec'
!               23........
!             . .      . .  
!           .   .    .   . 
!         .     .  .     . 
!       .       ..       . STAGGERING
!     .        ..        .   OF THE 
!   .        .  3.......31    GRID
! 2.......12  .        .    
! .        ..        .           3
! .       ..       .            /
! .     .  .     .             .---> 1
! .   .    .   .               |
! . .      . .                 |
! kk.......1                   2
!                    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!      3
!      ^        |  3   31 |  3   31 |  3   31 |
!      |        |         |         |       : |
!      o-->1    | kk    1 | kk    1 | kk    1 |  NZF+1 
!...............|.........|.........|..............
!    fault>--------3---31----3---31----3---31---
!               |         |         |       : |  NZF 
!               | kk    1 | kk    1 | kk    1 |
!...............|.........|.........|.................
!               |      31 |      31 |      31 |
!               |         |         |       : |  NZF-1 
!               | kk    1 | kk    1 | kk    1 |
!...............|.........|.........|.................
!               |      31 |      31 |      31 |
!               |         |         |       : |  NZF-2 
!               | kk    1 | kk    1 | kk    1 |
!...............|.........|.........|.................
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      3
!      ^        | 23      | 23      | 23      |
!      |        |         |         |  :      | NZF+1
!      o-->1    |  2   12 |  2   12 |  2   12 |
!...............|.........|.........|..............
!    fault>-------23--------23--------23--------
!               |         |         |  :      | NZF
!               |  2   12 |  2   12 |  2   12 |
!...............|.........|.........|.................
!               | 23      | 23      | 23      |
!               |         |         |  :      | NZF-1
!               |  2   12 |  2   12 |  2   12 |
!...............|.........|.........|.................
!               | 23      | 23      | 23      |
!               |         |         |  :      | NZF-2
!               |  2   12 |  2   12 |  2   12 |
!...............|.........|.........|.................
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                              
! EXTRAPOL SPEED:
!  4th order interpolation velo. in"third" layer:
      do k=3,nz-2
      do j=3,ny-2
      do i=3,nx-2
      v1(i,j,k)= v1(i,j,k)+ &
      (dt/(dx*2./((1/rho(i,j,k)+1/rho(i+1,j,k)))))*( &
      rn*(s11(i+1,j,k)-s11(i,j,k))+rnn*(s11(i+2,j,k)-s11(i-1,j,k))+ &
      rn*(s12(i,j+1,k)-s12(i,j,k))+rnn*(s12(i,j+2,k)-s12(i,j-1,k))+ &
      rn*(s13(i,j,k)-s13(i,j,k-1))+rnn*(s13(i,j,k+1)-s13(i,j,k-2)))
      v2(i,j,k)= v2(i,j,k)+ &
      (dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j-1,k)))))*( &
      rn*(s22(i,j,k)-s22(i,j-1,k))+rnn*(s22(i,j+1,k)-s22(i,j-2,k))+ &
      rn*(s23(i,j,k)-s23(i,j,k-1))+rnn*(s23(i,j,k+1)-s23(i,j,k-2))+ &
      rn*(s12(i,j,k)-s12(i-1,j,k))+rnn*(s12(i+1,j,k)-s12(i-2,j,k)))
      v3(i,j,k)= v3(i,j,k)+ &
      (dt/(dx*2./((1/rho(i,j,k)+1/rho(i,j,k+1)))))*( &
      rn*(s33(i,j,k+1)-s33(i,j,k))+rnn*(s33(i,j,k+2)-s33(i,j,k-1))+ &
      rn*(s23(i,j+1,k)-s23(i,j,k))+rnn*(s23(i,j+2,k)-s23(i,j-1,k))+ &
      rn*(s13(i,j,k)-s13(i-1,j,k))+rnn*(s13(i+1,j,k)-s13(i-2,j,k)))
      enddo
      enddo
      enddo
      write (11,*) v1(nx/4,ny/4,nz/4),v2(nx/4,ny/4,nz/4)
      return         
      end
!cccccccccccccccccccccccccccccccccccccccc
       subroutine save
       include 'triffy.dec'
       real bluff(ny,nz)
       real temp(nx,ny)
       character*16 namx1
       character*16 namx2
       character*16 namy
       character*16 namz
       character*16 namz1
       character*16 namz2
       real*8 r8
       integer*4 it4
       integer*4 ip4
       integer kmi
!      here files trace* are replenished with vertical
!      displacement at the surface. v2 can be replaced by v1
!      or any combination of the 2 components. The convention
!      here is that the velocity is the interpolated value at
!      node s11/s22, so that all velocity components are
!      known at the same spot.
!
       is1=is1+1
       is2=is2+1
!ccccc
       backspace (26)
       write (26,*) 'current iteration ',it
       call flush(26)
!%%% save fields for snapshots:
       if (savefields.eq.1.and.is1.ge.nos1) then
         print*,'current iteration ',it
         is1=0
         numf=numf+1
         write (namx1,'(a,I5.5)') './RES/ss13',it
         write (namx2,'(a,I5.5)') './RES/ss23',it
         write (namy,'(a,I5.5)') './RES/sli1',it
         write (namz,'(a,I5.5)') './RES/sli2',it
         write (namz1,'(a,I5.5)') './RES/rdi1',it
         write (namz2,'(a,I5.5)') './RES/rdi2',it

!c field1:
         open (50,file=namx1,form='formatted')
         open (51,file=namx2,form='formatted')
         open (52,file=namy,form='formatted')
         open (53,file=namz,form='formatted')
         open (41,file=namz1,form='formatted')
         open (42,file=namz2,form='formatted')

         k=nzf
         do j=1,nz
          write(50,79) (s13(i,j,k),i=1,nx)
          write(51,79) (s23(i,j,k),i=1,nx)
         enddo
         do j=1,nz
          write(52,79) (slip_mem1(i,j),i=1,nx)
          write(53,79) (slip_mem2(i,j),i=1,nx)
          write(41,79) (v1(i,j,nzf),i=1,nx)
          write(42,79) (v2(i,j,nzf),i=1,nx)
         enddo
         close (50)
         close (51)
         close (52)
         close (53)
         close (41)
         close (42)
!cccc
       endif
79     format(1000G17.9)
! save receivers
         do ir=1,numsta
          write(60+ir,*) v1(nxr(ir),nyr(ir),nzf)
          write(60+numsta+ir,*) v2(nxr(ir),nyr(ir),nzf)
          call flush(60+ir);call flush(60+numsta+ir)
         enddo
       return
       end 
!cccccccccccccccccccccccccccccccccccccccccc
      subroutine rupture
      include 'triffy.dec'
      real dv,rns,rad,ftw
      k=nzf
!cc   check for rupture:
      do j=nyf1,nyf2
       do i=nxf1,nxf2
        rns=( s33(i,j,k+1)+s33(i,j,k) ) / 2.
        rad=dx*sqrt(float(nxf1-i)**2+float(nyf1-j)**2)
        rstre(i,j)=sqrt((s23(i,j,k)+strini23(i,j))**2+s13(i,j,k)**2)
        if(rstre(i,j).ge.peak+rns*mu_s.and.&
         broken(i,j).eq.0.and.incrack(i,j).eq.1 ) then
          peak_tf(i,j)=rstre(i,j)*.99
          broken(i,j)=1
          rdis(i,j)=0.
          write (37,*) i,j,it
          nbr=nbr+1
        endif
!cxcxc
        deltah=delta
        if (i.gt.nxf1+120) &
       deltah=delta+(deltamax-delta)*float(i-(nxf1+120))/40.
        if (FricType.eq.'stat') then
         if (broken(i,j).eq.0) then
           theta(i,j)=theta(i,j)+(dt*((1.-theta(i,j))/tau2))
           if (theta(i,j).lt.0.) then
              theta(i,j)=0.
           elseif (theta(i,j).gt.1.) then
              theta(i,j)=1.
           endif
         else  ! broken :
           theta(i,j)=theta(i,j)+dt*( ((1.-theta(i,j))/tau2) &
                    -theta(i,j)*abs(rdis(i,j))/deltah)
           if (theta(i,j).lt.0.) then
              theta(i,j)=0.
           elseif (theta(i,j).gt.1.) then
              theta(i,j)=1.
           endif
           friction(i,j)=peak_tf(i,j)*theta(i,j)
           rdis(i,j)=-2.*v2(i,j,nzf)
           slip(i,j)=slip(i,j)+Abs(rdis(i,j))*dt
           slip_mem1(i,j)=slip_mem1(i,j)+Abs(rdis1(i,j))*dt
           slip_mem2(i,j)=slip_mem2(i,j)+Abs(rdis2(i,j))*dt
!          s23(i,j,k)=friction(i,j)-strini23(i,j)
         endif
        elseif (FricType.eq.'weak'.or.FricType.eq.' ') then
          if (broken(i,j).eq.1) then
            if(deltah.eq.0.)then
              dd=0.
            else
              dd=deltah/(deltah+slip(i,j))
!             dd=1. - slip(i,j)/deltah
              if (dd.lt.0.) dd=0.
            endif
            if (vmax.eq.0.) then
              dv=0.
            else
              dv=vmax/(vmax+Abs(rdis(i,j)))
            endif
            if (dv.gt.dd) dd=dv
            friction(i,j)=dd*peak_tf(i,j)
            rdis1(i,j)=v1(i,j,nzf+1)-v1(i,j,nzf)
            rdis2(i,j)=v2(i,j,nzf+1)-v2(i,j,nzf)
            rdis(i,j) =sqrt(rdis1(i,j)**2+rdis2(i,j)**2)
            slip(i,j)=slip(i,j)+rdis(i,j)*dt
            slip_mem1(i,j)=slip_mem1(i,j)+rdis1(i,j)*dt
            slip_mem2(i,j)=slip_mem2(i,j)+rdis2(i,j)*dt
!           s23(i,j,k)=friction(i,j)-strini23(i,j)
         endif
        endif
        enddo
       enddo
       do i=nxf1,nxf2
        do j=nyf1,nyf2
          if (broken(i,j).eq.1) then
!cc   check for healing:
!          if (rdis(i,j).lt.0.) then
           if (friction(i,j).gt.rstre(i,j)) then
             broken(i,j)=0
             rdis(i,j)=0.
             rdis1(i,j)=0.
             rdis2(i,j)=0.
             nbr=nbr-1
!            incrack(i,j)=0 ! avoid repeated fracture if 0
             slip(i,j)=0.
             write(39,*) i,j,it
           endif
          endif
        enddo
       enddo
       write (40,*) nbr
!cc update stress:
       where(broken.eq.1)
           s23(:,:,k)=friction*(s23(:,:,k)+strini23)/rstre-strini23
           s13(:,:,k)=friction*s13(:,:,k)/rstre
       endwhere
       call flush (39)
       call flush (37)
       call flush (40)
       return
       end
!ccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine imposed_source
      include 'triffy.dec'
! STRESS SOURCE:
      rf=0.
      if (float(it)*dt.le.10.*tw) then
        rf=exp(-((float(it)*dt-5.*tw)/tw)**2)
        s11(nxs,nys,nzs)=s11(nxs,nys,nzs)+rf
        s22(nxs,nys,nzs)=s22(nxs,nys,nzs)+rf
        s33(nxs,nys,nzs)=s33(nxs,nys,nzs)+rf
        write (12,'(I4,F17.9)') it,rf
        call flush(15)
      endif
      return
      end
!cccccccccccccccccccccccccccccccc
       subroutine readfiles(fmodel,il)
       include 'triffy.dec'
       integer il
       real rtemp1(nz+2*npm-6)
       real rtemp2(nz+2*npm-6)
       real rtemp3(nz+2*npm-6)
       character*80 fmodel
       open(12,file=fmodel(1:il)//'.mu', form='formatted')
       open(13,file=fmodel(1:il)//'.lam',form='formatted')
       open(14,file=fmodel(1:il)//'.rho',form='formatted')
       do j=1,npm
        read(12,*) (rtemp1(k),k=1,nz+2*npm-6)
        read(13,*) (rtemp2(k),k=1,nz+2*npm-6)
        read(14,*) (rtemp3(k),k=1,nz+2*npm-6)
! FILL jk:
        do i=1,nx
        do k=1,npm
          muBYBZ(i,j,k)=rtemp1(k)
         lamBYBZ(i,j,k)=rtemp2(k)
         rhoBYBZ(i,j,k)=rtemp3(k)
        enddo
        enddo
! FILL j:
        do i=1,nx
        do k=1,nz
          muBY(i,j,k)=rtemp1(k+npm-3)
         lamBY(i,j,k)=rtemp2(k+npm-3)
         rhoBY(i,j,k)=rtemp3(k+npm-3)
        enddo
        enddo
! FILL jk:
        do i=1,nx
        do k=1,npm
          muBYTZ(i,j,k)=rtemp1(k+nz+npm-6)
         lamBYTZ(i,j,k)=rtemp2(k+nz+npm-6)
         rhoBYTZ(i,j,k)=rtemp3(k+nz+npm-6)
        enddo
        enddo
       enddo
!
       do k=1,3*(nz+2*npm-6)
        backspace(12)
        backspace(13)
        backspace(14)
       enddo
!
       do j=1,ny
        read(12,*) (rtemp1(k),k=1,nz+2*npm-6)
        read(13,*) (rtemp2(k),k=1,nz+2*npm-6)
        read(14,*) (rtemp3(k),k=1,nz+2*npm-6)
! FILL k:
        do i=1,nx
        do k=1,npm
          muBZ(i,j,k)=rtemp1(k)
         lamBZ(i,j,k)=rtemp2(k)
         rhoBZ(i,j,k)=rtemp3(k)
        enddo
        enddo
! FILL CENTRAL BLOC:
        do i=1,nx
        do k=1,nz
          mu(i,j,k)=rtemp1(k+npm-3)
         lam(i,j,k)=rtemp2(k+npm-3)
         rho(i,j,k)=rtemp3(k+npm-3)
        enddo
        enddo

! FILL k:
        do i=1,nx
        do k=1,npm
          muTZ(i,j,k)=rtemp1(k+nz+npm-6)
         lamTZ(i,j,k)=rtemp2(k+nz+npm-6)
         rhoTZ(i,j,k)=rtemp3(k+nz+npm-6)
        enddo
        enddo
       enddo
! PER IL TOP 3 BACKSPACE (overlap fondo/centro)
       do k=1,3*(nz+2*npm-6)
        backspace(12)
        backspace(13)
        backspace(14)
       enddo
! TOP:
       do j=1,npm
        read(12,*) (rtemp1(k),k=1,nz+2*npm-6)
        read(13,*) (rtemp2(k),k=1,nz+2*npm-6)
        read(14,*) (rtemp3(k),k=1,nz+2*npm-6)
! FILL jk:
        do i=1,nx
        do k=1,npm
          muTYBZ(i,j,k)=rtemp1(k)
         lamTYBZ(i,j,k)=rtemp2(k)
         rhoTYBZ(i,j,k)=rtemp3(k)
        enddo
        enddo
! FILL j:
        do i=1,nx
        do k=1,nz
          muTY(i,j,k)=rtemp1(k+npm-3) 
         lamTY(i,j,k)=rtemp2(k+npm-3)
         rhoTY(i,j,k)=rtemp3(k+npm-3)
        enddo
        enddo
! FILL jk:
        do i=1,nx
        do k=1,npm
          muTYTZ(i,j,k)=rtemp1(k+nz+npm-6)
         lamTYTZ(i,j,k)=rtemp2(k+nz+npm-6)
         rhoTYTZ(i,j,k)=rtemp3(k+nz+npm-6)
        enddo
        enddo
       enddo
! FILL i:
       do i=1,npm
       do j=1,ny
       do k=1,nz
         muBX(i,j,k)= mu(1,j,k) 
         muTX(i,j,k)= mu(1,j,k) 
        lamBX(i,j,k)=lam(1,j,k) 
        lamTX(i,j,k)=lam(1,j,k) 
        rhoBX(i,j,k)=rho(1,j,k) 
        rhoTX(i,j,k)=rho(1,j,k) 
       enddo 
       enddo 
       enddo 
! FILL ik
       do i=1,npm
       do j=1,ny
       do k=1,npm
         muBXTZ(i,j,k)= muTZ(1,j,k) 
        lamBXTZ(i,j,k)=lamTZ(1,j,k) 
        rhoBXTZ(i,j,k)=rhoTZ(1,j,k) 
         muTXTZ(i,j,k)= muTZ(1,j,k) 
        lamTXTZ(i,j,k)=lamTZ(1,j,k) 
        rhoTXTZ(i,j,k)=rhoTZ(1,j,k) 

         muBXBZ(i,j,k)= muBZ(1,j,k) 
        lamBXBZ(i,j,k)=lamBZ(1,j,k) 
        rhoBXBZ(i,j,k)=rhoBZ(1,j,k) 
         muTXBZ(i,j,k)= muBZ(1,j,k) 
        lamTXBZ(i,j,k)=lamBZ(1,j,k) 
        rhoTXBZ(i,j,k)=rhoBZ(1,j,k) 
       enddo
       enddo
       enddo 
! FILL ij
       do i=1,npm
       do j=1,npm
       do k=1,nz
         muBXTY(i,j,k)= muTY(1,j,k) 
        lamBXTY(i,j,k)=lamTY(1,j,k) 
        rhoBXTY(i,j,k)=rhoTY(1,j,k) 
         muTXTY(i,j,k)= muTY(1,j,k) 
        lamTXTY(i,j,k)=lamTY(1,j,k) 
        rhoTXTY(i,j,k)=rhoTY(1,j,k) 

         muBXBY(i,j,k)= muBY(1,j,k) 
        lamBXBY(i,j,k)=lamBY(1,j,k) 
        rhoBXBY(i,j,k)=rhoBY(1,j,k) 
         muTXBY(i,j,k)= muBY(1,j,k) 
        lamTXBY(i,j,k)=lamBY(1,j,k) 
        rhoTXBY(i,j,k)=rhoBY(1,j,k) 
       enddo
       enddo
       enddo 
! FILL ijk
       do i=1,npm
       do j=1,npm
       do k=1,npm
         muBXBYBZ(i,j,k)= muBYBZ(1,j,k) 
        lamBXBYBZ(i,j,k)=lamBYBZ(1,j,k) 
        rhoBXBYBZ(i,j,k)=rhoBYBZ(1,j,k) 
         muTXBYBZ(i,j,k)= muBYBZ(1,j,k) 
        lamTXBYBZ(i,j,k)=lamBYBZ(1,j,k) 
        rhoTXBYBZ(i,j,k)=rhoBYBZ(1,j,k) 

         muBXBYTZ(i,j,k)= muBYTZ(1,j,k) 
        lamBXBYTZ(i,j,k)=lamBYTZ(1,j,k) 
        rhoBXBYTZ(i,j,k)=rhoBYTZ(1,j,k) 
         muTXBYTZ(i,j,k)= muBYTZ(1,j,k) 
        lamTXBYTZ(i,j,k)=lamBYTZ(1,j,k) 
        rhoTXBYTZ(i,j,k)=rhoBYTZ(1,j,k) 

         muBXTYBZ(i,j,k)= muTYBZ(1,j,k) 
        lamBXTYBZ(i,j,k)=lamTYBZ(1,j,k) 
        rhoBXTYBZ(i,j,k)=rhoTYBZ(1,j,k) 
         muTXTYBZ(i,j,k)= muTYBZ(1,j,k) 
        lamTXTYBZ(i,j,k)=lamTYBZ(1,j,k) 
        rhoTXTYBZ(i,j,k)=rhoTYBZ(1,j,k) 

         muBXTYTZ(i,j,k)= muTYTZ(1,j,k) 
        lamBXTYTZ(i,j,k)=lamTYTZ(1,j,k) 
        rhoBXTYTZ(i,j,k)=rhoTYTZ(1,j,k) 
         muTXTYTZ(i,j,k)= muTYTZ(1,j,k) 
        lamTXTYTZ(i,j,k)=lamTYTZ(1,j,k) 
        rhoTXTYTZ(i,j,k)=rhoTYTZ(1,j,k) 

       enddo
       enddo
       enddo 

       close(12)
       close(13)
       close(14)
!
      return
      end
!cccccccccccccccccccccccccccccccccccccccccccc
       integer function getlen (string)
       character *(*) string
       do i = len(string), 1, -1
          if ( (string(i:i) .ne. ' ') .and.  &
         (ichar(string(i:i)) .ne. 0) ) then 
             getlen = i
             return
          endif
       enddo
       write(6,'(/,"  Warning: blank string passed as name.",/)')
       getlen = 0
       return
       end
!ccccccccccccccccccccccccccccccccccccccc
       subroutine iniPMLini(muAll,lamAll,rhoAll)
       include 'triffy.dec'
       real muAll,lamAll,rhoAll

          muBYBZ=muAll
         lamBYBZ=lamAll
         rhoBYBZ=rhoAll
          muBY=muAll
         lamBY=lamAll
         rhoBY=rhoAll
          muBYTZ=muAll
         lamBYTZ=lamAll
         rhoBYTZ=rhoAll
          muBZ=muAll
         lamBZ=lamAll
         rhoBZ=rhoAll
          mu=muAll
         lam=lamAll
         rho=rhoAll
          muTZ=muAll
         lamTZ=lamAll
         rhoTZ=rhoAll
          muTYBZ=muAll
         lamTYBZ=lamAll
         rhoTYBZ=rhoAll
          muTY=muAll
         lamTY=lamAll
         rhoTY=rhoAll
          muTYTZ=muAll
         lamTYTZ=lamAll
         rhoTYTZ=rhoAll
         muBX=muAll
         muTX=muAll
        lamBX=lamAll
        lamTX=lamAll
        rhoBX=rhoAll
        rhoTX=rhoAll
         muBXTZ=muAll
        lamBXTZ=lamAll
        rhoBXTZ=rhoAll
         muTXTZ=muAll
        lamTXTZ=lamAll
        rhoTXTZ=rhoAll
         muBXBZ=muAll
        lamBXBZ=lamAll
        rhoBXBZ=rhoAll
         muTXBZ=muAll
        lamTXBZ=lamAll
        rhoTXBZ=rhoAll
         muBXTY=muAll
        lamBXTY=lamAll
        rhoBXTY=rhoAll
         muTXTY=muAll
        lamTXTY=lamAll
        rhoTXTY=rhoAll
         muBXBY=muAll
        lamBXBY=lamAll
        rhoBXBY=rhoAll
         muTXBY=muAll
        lamTXBY=lamAll
        rhoTXBY=rhoAll
         muBXBYBZ=muAll
        lamBXBYBZ=lamAll
        rhoBXBYBZ=rhoAll
         muTXBYBZ=muAll
        lamTXBYBZ=lamAll
        rhoTXBYBZ=rhoAll
         muBXBYTZ=muAll
        lamBXBYTZ=lamAll
        rhoBXBYTZ=rhoAll
         muTXBYTZ=muAll
        lamTXBYTZ=lamAll
        rhoTXBYTZ=rhoAll
         muBXTYBZ=muAll
        lamBXTYBZ=lamAll
        rhoBXTYBZ=rhoAll
         muTXTYBZ=muAll
        lamTXTYBZ=lamAll
        rhoTXTYBZ=rhoAll
         muBXTYTZ=muAll
        lamBXTYTZ=lamAll
        rhoBXTYTZ=rhoAll
         muTXTYTZ=muAll
        lamTXTYTZ=lamAll
        rhoTXTYTZ=rhoAll
      return
      end
!cccccccccccccccccccc

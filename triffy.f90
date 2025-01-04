!*************************
! PROGRAM Triffy
! A bloc with absorbing Perfect Matching Layers on all sides.
! Matching layer is overlapping with first 3 nodes on all sides
! so actual medium starts at i=4 and stops at i=nx-3 or (nxf1 to nxf2)  
! and similaly for j and k.
!
! Free surface can also be implemented at j=3
! by calling routine freesurfacesym. 
!
! An extra absorbing cylinder is added at the top edges of the model
! in routine corner_abs
! to attenuate surface wave diffraction from there. 
! 
! *************************
! Stefan Nielsen, Gaetano Festa September 2000
! Revised version by Stefan Nielsen December 2024
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
 call freesurfacesym
 call save                  !visualize and/or store variables 
 call extrapolate           !extrapolate in time
!call corner_abs
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
       integer ux,uy,uz
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
       !nxs=20
       !nys=20
       nzs=nz/2
!      fmodel='model'
!      il=getlen(fmodel)
!      call readfiles(fmodel,il)
!      call MediumPMLparameters
! MEDIUM PARAMETERS - NO INPUT FILE:
      muAll=25000.0;lamAll=25000.0;rhoAll=.0028 !real values but all divided by 1e6!!
      call iniPMLini(muAll,lamAll,rhoAll)
! ADD a spherical inhomogeniety in center of model:
      inrad=5.0
      nxh=float(nx)/2.;nyh=float(ny)/2.;nzh=float(nzs)
      do k=1,nz
       do j=1,ny
        do i=1,nx
         dista = sqrt((float(i)-nxh)**2+ (float(j)-nyh)**2+ (float(k)-nzh)**2)
         if (dista .le. inrad) then
          mu(i,j,k)  = 0.0
          lam(i,j,k) = 0.0
          rho(i,j,k) = 1.0  ! caution: rho not allowed to be zero! take small value instead
         endif
        enddo
       enddo
      enddo
! compute muh : the harmonic mean of mu to match staggered grid.
! If any value is zero then use standard mean instead.
! 
      do k=1,nz
       do j=1,ny
        do i=1,nx
         if (mu(i,j,k) .eq. 0. .or. mu(i+1,j,k) .eq. 0.  &
             .or. mu(i+1,j-1,k) .eq. 0. .or. mu(i,j-1,k) .eq. 0.) then
          muh(i,j,k) = (mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j-1,k)+mu(i,j-1,k))/4.
         else
          muh(i,j,k)=(4./(1./mu(i,j,k)+ 1./mu(i+1,j,k)+1./mu(i+1,j-1,k)+ 1./mu(i,j-1,k)))
         endif
        enddo
       enddo
      enddo

! Fault parameters- NO INPUT FILE:
      nzf=nz/2 ! this also serves as mid-plane for 2D/3d
      nxf1=4 
      nxf2=nx-3
      nyf1=4      ! CAUTION : this is defined after free surface
      nyf2=ny-3   ! CAUTION : this "    .   .   .   .   .      "
      nzf1=4      ! CAUTION : this "    .   .   .   .   .      "
      nzf2=nz-3   ! CAUTION : this "    .   .   .   .   .      "
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
       ux=nxf2-nxf1; uy=nyf2-nyf1; uz=nzf2-nzf1;
       print '(A29,I3.3,A8,I3.3,A3)', 'useful horizontal surface: x=',ux,' dx | z=',uz,' dx.'
       print '(A29,I3.3,A8,I3.3,A3)', 'useful  vertical  section: x=',ux,' dx | y=',uy,' dx.'
       nr=0
       print*,'receivers at:'
       do while (1==1)
        read (11,*,END=51) rxr,ryr,rzr
        nr=nr+1
        nxr(nr)=rxr/dx+nxf1
        nyr(nr)=ryr/dx+nyf1
        nzr(nr)=rzr/dx+nzf1
        print '(A13,I3.3,A5,I3.3,A5,I3.3,A5)', &
          '(x,y,z): ',nxr(nr)-nxf1,' dx, ', nyr(nr)-nyf1,' dx, ',nzr(nr)-nzf1,' dx. '
       enddo
51     close (11)
       print*,'total  nr of receivers =',nr
       !real rux,ruy,ruz
       !rux=dx*float(nxf2-nxf1); ruy=dx*float(nyf2-nyf1); ruz=dx*float(nzf2-nzf1);
       InitRadius=nint(rInitRadius/dx)
       open (unit=24,file='filpar.dec',form='formatted')
       write (24,'(a19,G8.4,a1)') '      parameter(dt=',dt,')'
       rdis=0.
       broken=0
       incrack=1
       strini23=0.
       nxs=int(rxs/dx)+nxf1
       nys=int(rys/dx)+nyf1-1
! the dirac functions remove forth order in layers nzf+1 and nzf-1 
! form some terms, eventually add non zero dissipation:
      dirac=0.
      dirac(nzf)=1. 
! EFFECTIVE MEDIUM PARAMETERS: computed in iterations to save mem.
! OUTPUT FILES, OPEN:
       open (unit=15,file='source',form='formatted')
       open (unit=26, file='iter',form='formatted')
       do i=1,nr !open single trace files for output:
        write (nam,'(a,I3.3)') './RES/tracex',i
        open (60+i,file=nam,form='formatted')
        write (nam,'(a,I3.3)') './RES/tracey',i
        open (60+nr+i,file=nam,form='formatted')
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
         muh(i,j,k)*(dt/dx)*( &
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
       return
       end

! *******************************************
        subroutine freesurfacesym
        include 'triffy.dec'
        ! free surface by "mirror" method:
        ! surface at node j=3 goes through plane of s23&s12
        ! s22(j=2) = -s22(j=3)
        ! s22(j=1) = -s22(j=4)
        ! s23&s12(j=3) = 0
        ! s23&s12(j=2) = -s23&s12(j=4)
        ! s23&s12(j=1) = -s23&s12(j=5)
        ! 
        do k=4,nz-3
          do i=4,nx-3
           s22(i,2,k) = -s22(i,3,k);
           s22(i,1,k) = -s22(i,4,k);
           s23(i,3,k) = 0.; s12(i,3,k)=0.;
           s23(i,2,k) = -s23(i,4,k);s12(i,2,k)= -s12(i,4,k);  
           s23(i,1,k) = -s23(i,5,k);s12(i,1,k)= -s12(i,5,k);  
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
      subroutine corner_abs
      include 'triffy.dec'
      real d1,d2,d3,d4,rlen
      dis=dissipation*dt
      rlen=5.
      do k=3,nz-2
      do j=3,ny-2
      do i=3,nx-2
       d1=dis*exp(-sqrt(float(i-3)**2+float(j-3)**2)/rlen)
       d2=dis*exp(-sqrt(float(nx-2-i)**2+float(j-3)**2)/rlen)
       d3=dis*exp(-sqrt(float(nz-2-k)**2+float(j-3)**2)/rlen)
       d4=dis*exp(-sqrt(float(k-3)**2+float(j-3)**2)/rlen)
       v1(i,j,k)= v1(i,j,k)*(1.-d1)
       v1(i,j,k)= v1(i,j,k)*(1.-d2)
       v1(i,j,k)= v1(i,j,k)*(1.-d3)
       v1(i,j,k)= v1(i,j,k)*(1.-d4)
       v2(i,j,k)= v2(i,j,k)*(1.-d1)
       v2(i,j,k)= v2(i,j,k)*(1.-d2)
       v2(i,j,k)= v2(i,j,k)*(1.-d3)
       v2(i,j,k)= v2(i,j,k)*(1.-d4)
      enddo
      enddo
      enddo
      write (11,*) v1(nx/4,ny/4,nz/4),v2(nx/4,ny/4,nz/4)
      return         
      end
!cccccccccccccccccccccccccccccccccccccccccccccc      
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
         do j=nyf1,nyf2
          write(50,79) (s13(i,j,k),i=nxf1,nxf2)
          write(51,79) (s23(i,j,k),i=nxf1,nxf2)
         enddo
         do j=nyf1,nyf2
          write(52,79) (slip_mem1(i,j),i=nxf1,nxf2)
          write(53,79) (slip_mem2(i,j),i=nxf1,nxf2)
          write(41,79) (v1(i,j,nzf),i=nxf1,nxf2)
          write(42,79) (v2(i,j,nzf),i=nxf1,nxf2)
          !print *, 'writingi i=',nxf1,'->',nxf2,'j=',j,'k=',nzf
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
         do ir=1,nr
          write(60+ir,*) v1(nxr(ir),nyr(ir),nxr(ir))
          write(60+nr+ir,*) v2(nxr(ir),nyr(ir),nzr(ir))
          call flush(60+ir);call flush(60+nr+ir)
         enddo
       return
       end 
!cccccccccccccccccccccccccccccccccccccccccc
      subroutine imposed_source
      include 'triffy.dec'
! STRESS SOURCE:
      rf=0.
      if (float(it)*dt.le.10.*tw) then
        rf=exp(-((float(it)*dt-5.*tw)/tw)**2)
        !s11(nxs,nys,nzs)=s11(nxs,nys,nzs)+rf
        !s22(nxs,nys,nzs)=s22(nxs,nys,nzs)+rf
        !s33(nxs,nys,nzs)=s33(nxs,nys,nzs)+rf
        ! add acceleration to simulate vertical force:
        v2(nxs,nys,nzs) = v2(nxs,nys,nzs) + rf
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

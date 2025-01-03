Subroutine Initdist
include 'triffy.dec'
real lama,mua,rhoa,xx
! An eterogeneous medium is considered
! We define the initilisation for all the layers

! TX
do k=3,nz-2
 do j=3,ny-2
  do i=2,npm

lama=lamTX(i,j,k)
mua=muTX(i,j,k)
rhoa=rhoTX(i,j,k)
vp=sqrt((2*lama+mua)/rhoa) 
dist=4.5*vp/npm/dx
 xx=float(i-2)
 DiTXS0(i,j,k)=dist*(xx/float(npm))**3
 xx=xx+.5
 DiTXS1(i,j,k)=dist*(xx/float(npm))**3

  enddo
 enddo
enddo 
!BX
do k=3,nz-2
 do j=3,ny-2
   do i=1,npm-1

   lama=lamBX(i,j,k)
   mua=muBX(i,j,k)
   rhoa=rhoBX(i,j,k)
   vp=sqrt((2*lama+mua)/rhoa)
   dist=4.5*vp/npm/dx
    xx=float(npm-i)
     DiBXS0(i,j,k)=dist*(xx/float(npm))**3
      xx=xx-.5
       DiBXS1(i,j,k)=dist*(xx/float(npm))**3

    enddo
   enddo
  enddo 

! TY 
do k=3,nz-2
 do j=2,npm
  do i=3,nx-2
  
lama=lamTY(i,j,k)
mua=muTY(i,j,k)
rhoa=rhoTY(i,j,k)
vp=sqrt((2*lama+mua)/rhoa) 
dist=4.5*vp/npm/dx
 xx=float(j-1)
 DiTYS0(i,j,k)=dist*(xx/float(npm))**3
 xx=xx-.5
 DiTYS1(i,j,k)=dist*(xx/float(npm))**3
  enddo
 enddo
enddo

! BY 
do k=3,nz-2
 do j=1,npm-1
   do i=3,nx-2

   lama=lamBY(i,j,k)
   mua=muBY(i,j,k)
   rhoa=rhoBY(i,j,k)
   vp=sqrt((2*lama+mua)/rhoa)
   dist=4.5*vp/npm/dx
    xx=float(npm-1-j)
     DiBYS0(i,j,k)=dist*(xx/float(npm))**3
      xx=xx+.5
       DiBYS1(i,j,k)=dist*(xx/float(npm))**3
         enddo
        enddo
       enddo

! TZ
do k=2,npm
 do j=3,ny-2
  do i=3,nx-2
lama=lamTZ(i,j,k)
mua=muTZ(i,j,k)
rhoa=rhoTZ(i,j,k)
vp=sqrt((2*lama+mua)/rhoa) 
dist=4.5*vp/npm/dx
 xx=float(k-2)
 DiTZS0(i,j,k)=dist*(xx/float(npm))**3
 xx=xx+.5
 DiTZS1(i,j,k)=dist*(xx/float(npm))**3
enddo
enddo
enddo

! BZ
do k=1,npm-1
 do j=3,ny-2
   do i=3,nx-2
   lama=lamBZ(i,j,k)
   mua=muBZ(i,j,k)
   rhoa=rhoBZ(i,j,k)
   vp=sqrt((2*lama+mua)/rhoa)
   dist=4.5*vp/npm/dx
    xx=float(npm-k)
     DiBZS0(i,j,k)=dist*(xx/float(npm))**3
      xx=xx-.5
       DiBZS1(i,j,k)=dist*(xx/float(npm))**3
      enddo
     enddo
    enddo
 


return
end

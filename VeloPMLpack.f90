! ##############################################
Subroutine vPMLTZ
include 'triffy.dec'
do k=2,npm-1 
 do j=3,ny-2
  do i=3,nx-2
   daf3x=DiTZS0(i,j,k)
   c03x=1+daf3x 
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(i,j,k)
   c03z=1+daf3z
!
   bo1=(1/rhoTZ(i,j,k)+1/rhoTZ(i+1,j,k))/2
   bo2=(1/rhoTZ(i,j,k)+1/rhoTZ(i,j-1,k))/2
   bo3=(1/rhoTZ(i,j,k)+1/rhoTZ(i,j,k+1))/2
!
   v1oTZ(i,j,k)=(1-daf3x)/c03x*v1oTZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13TZ(i,j,k)-s13TZ(i,j,k-1))
   v1pTZ(i,j,k)=v1pTZ(i,j,k)+bo1*dtx*(rn*(s11TZ(i+1,j,k)-     &  
    s11TZ(i,j,k)+s12TZ(i,j+1,k)-s12TZ(i,j,k))+        &
    rnn*(s11TZ(i+2,j,k)-s11TZ(i-1,j,k)+         &
    s12TZ(i,j+2,k)-s12TZ(i,j-1,k)))
   v1TZ(i,j,k)=v1oTZ(i,j,k)+v1pTZ(i,j,k)
!
   v2oTZ(i,j,k)=(1-daf3y)/c03y*v2oTZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23TZ(i,j,k)-s23TZ(i,j,k-1))
   v2pTZ(i,j,k)=v2pTZ(i,j,k)+bo2*dtx*(rn*(s12TZ(i,j,k)-       &  
    s12TZ(i-1,j,k)+s22TZ(i,j,k)-s22TZ(i,j-1,k))+  &
    rnn*(s12TZ(i+1,j,k)-s12TZ(i-2,j,k)+         &
    s22TZ(i,j+1,k)-s22TZ(i,j-2,k)))
   v2TZ(i,j,k)=v2oTZ(i,j,k)+v2pTZ(i,j,k)
!
   v3oTZ(i,j,k)=(1-daf3z)/c03z*v3oTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TZ(i,j,k+1)-s33TZ(i,j,k))
   v3pTZ(i,j,k)=v3pTZ(i,j,k)+bo3*dtx*(rn*(s13TZ(i,j,k)-       &
    s13TZ(i-1,j,k)+s23TZ(i,j+1,k)-s23TZ(i,j,k))+   &
    rnn*(s13TZ(i+1,j,k)-s13TZ(i-2,j,k)+           &  
    s23TZ(i,j+2,k)-s23TZ(i,j-1,k)))
   v3TZ(i,j,k)=v3oTZ(i,j,k)+v3pTZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine  
! ##############################################
Subroutine vPMLBZ
include 'triffy.dec'
do k=2,npm-1 
 do j=3,ny-2
  do i=3,nx-2
   daf3x=DiBZS0(i,j,k)
   c03x=1+daf3x 
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(i,j,k)
   c03z=1+daf3z
!
   bo1=(1/rhoBZ(i,j,k)+1/rhoBZ(i+1,j,k))/2
   bo2=(1/rhoBZ(i,j,k)+1/rhoBZ(i,j-1,k))/2
   bo3=(1/rhoBZ(i,j,k)+1/rhoBZ(i,j,k+1))/2
!
   v1oBZ(i,j,k)=(1-daf3x)/c03x*v1oBZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13BZ(i,j,k)-s13BZ(i,j,k-1))
   v1pBZ(i,j,k)=v1pBZ(i,j,k)+bo1*dtx*(rn*(s11BZ(i+1,j,k)-     &  
    s11BZ(i,j,k)+s12BZ(i,j+1,k)-s12BZ(i,j,k))+        &
    rnn*(s11BZ(i+2,j,k)-s11BZ(i-1,j,k)+         &
    s12BZ(i,j+2,k)-s12BZ(i,j-1,k)))
   v1BZ(i,j,k)=v1oBZ(i,j,k)+v1pBZ(i,j,k)
!
   v2oBZ(i,j,k)=(1-daf3y)/c03y*v2oBZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23BZ(i,j,k)-s23BZ(i,j,k-1))
   v2pBZ(i,j,k)=v2pBZ(i,j,k)+bo2*dtx*(rn*(s12BZ(i,j,k)-       &  
    s12BZ(i-1,j,k)+s22BZ(i,j,k)-s22BZ(i,j-1,k))+  &
    rnn*(s12BZ(i+1,j,k)-s12BZ(i-2,j,k)+         &
    s22BZ(i,j+1,k)-s22BZ(i,j-2,k)))
   v2BZ(i,j,k)=v2oBZ(i,j,k)+v2pBZ(i,j,k)
!
   v3oBZ(i,j,k)=(1-daf3z)/c03z*v3oBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BZ(i,j,k+1)-s33BZ(i,j,k))
   v3pBZ(i,j,k)=v3pBZ(i,j,k)+bo3*dtx*(rn*(s13BZ(i,j,k)-       &
    s13BZ(i-1,j,k)+s23BZ(i,j+1,k)-s23BZ(i,j,k))+   &
    rnn*(s13BZ(i+1,j,k)-s13BZ(i-2,j,k)+           &  
    s23BZ(i,j+2,k)-s23BZ(i,j-1,k)))
   v3BZ(i,j,k)=v3oBZ(i,j,k)+v3pBZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine  
! ##################################################
Subroutine vPMLTY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=3,nx-2
!
   daf2x=DiTYS0(i,j,k)
   c02x=1+daf2x  
   daf2y=DiTYS1(i,j,k)
   c02y=1+daf2y  
   daf2z=daf2x
   c02z=c02x
!
   bo1=(1/rhoTY(i,j,k)+1/rhoTY(i+1,j,k))/2
   bo2=(1/rhoTY(i,j,k)+1/rhoTY(i,j-1,k))/2
   bo3=(1/rhoTY(i,j,k)+1/rhoTY(i,j,k+1))/2
!
   v1oTY(i,j,k)=(1-daf2x)/c02x*v1oTY(i,j,k)+bo1*(dtx/c02x)*   &
    (s12TY(i,j+1,k)-s12TY(i,j,k))
   v1pTY(i,j,k)=v1pTY(i,j,k)+bo1*dtx*(rn*(s11TY(i+1,j,k)-     &  
    s11TY(i,j,k)+s13TY(i,j,k)-        &
    s13TY(i,j,k-1))+rnn*(s11TY(i+2,j,k)-s11TY(i-1,j,k)+         &
    s13TY(i,j,k+1)-s13TY(i,j,k-2)))
   v1TY(i,j,k)=v1oTY(i,j,k)+v1pTY(i,j,k)
!
   v2oTY(i,j,k)=(1-daf2y)/c02y*v2oTY(i,j,k)+bo2*(dtx/c02y)*   &
    (s22TY(i,j,k)-s22TY(i,j-1,k))
   v2pTY(i,j,k)=v2pTY(i,j,k)+bo2*dtx*(rn*(s12TY(i,j,k)-       &  
    s12TY(i-1,j,k)+s23TY(i,j,k)-      &
    s23TY(i,j,k-1))+rnn*(s12TY(i+1,j,k)-s12TY(i-2,j,k)+         &
    s23TY(i,j,k+1)-s23TY(i,j,k-2)))
   v2TY(i,j,k)=v2oTY(i,j,k)+v2pTY(i,j,k)
!
   v3oTY(i,j,k)=(1-daf2z)/c02z*v3oTY(i,j,k)+bo3*(dtx/c02z)*   &
    (s23TY(i,j+1,k)-s23TY(i,j,k))
   v3pTY(i,j,k)=v3pTY(i,j,k)+bo3*dtx*(rn*(s13TY(i,j,k)-       &  
    s13TY(i-1,j,k)+s33TY(i,j,k+1)-    &
    s33TY(i,j,k))+rnn*(s13TY(i+1,j,k)-s13TY(i-2,j,k)+           &
    s33TY(i,j,k+2)-s33TY(i,j,k-1)))
   v3TY(i,j,k)=v3oTY(i,j,k)+v3pTY(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ##################################################
Subroutine vPMLBY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=3,nx-2
!
   daf2x=DiBYS0(i,j,k)
   c02x=1+daf2x  
   daf2y=DiBYS1(i,j,k)
   c02y=1+daf2y  
   daf2z=daf2x
   c02z=c02x
!
   bo1=(1/rhoBY(i,j,k)+1/rhoBY(i+1,j,k))/2
   bo2=(1/rhoBY(i,j,k)+1/rhoBY(i,j-1,k))/2
   bo3=(1/rhoBY(i,j,k)+1/rhoBY(i,j,k+1))/2
!
   v1oBY(i,j,k)=(1-daf2x)/c02x*v1oBY(i,j,k)+bo1*(dtx/c02x)*   &
    (s12BY(i,j+1,k)-s12BY(i,j,k))
   v1pBY(i,j,k)=v1pBY(i,j,k)+bo1*dtx*(rn*(s11BY(i+1,j,k)-     &  
    s11BY(i,j,k)+s13BY(i,j,k)-        &
    s13BY(i,j,k-1))+rnn*(s11BY(i+2,j,k)-s11BY(i-1,j,k)+         &
    s13BY(i,j,k+1)-s13BY(i,j,k-2)))
   v1BY(i,j,k)=v1oBY(i,j,k)+v1pBY(i,j,k)
!
   v2oBY(i,j,k)=(1-daf2y)/c02y*v2oBY(i,j,k)+bo2*(dtx/c02y)*   &
    (s22BY(i,j,k)-s22BY(i,j-1,k))
   v2pBY(i,j,k)=v2pBY(i,j,k)+bo2*dtx*(rn*(s12BY(i,j,k)-       &  
    s12BY(i-1,j,k)+s23BY(i,j,k)-      &
    s23BY(i,j,k-1))+rnn*(s12BY(i+1,j,k)-s12BY(i-2,j,k)+         &
    s23BY(i,j,k+1)-s23BY(i,j,k-2)))
   v2BY(i,j,k)=v2oBY(i,j,k)+v2pBY(i,j,k)
!
   v3oBY(i,j,k)=(1-daf2z)/c02z*v3oBY(i,j,k)+bo3*(dtx/c02z)*   &
    (s23BY(i,j+1,k)-s23BY(i,j,k))
   v3pBY(i,j,k)=v3pBY(i,j,k)+bo3*dtx*(rn*(s13BY(i,j,k)-       &  
    s13BY(i-1,j,k)+s33BY(i,j,k+1)-    &
    s33BY(i,j,k))+rnn*(s13BY(i+1,j,k)-s13BY(i-2,j,k)+           &
    s33BY(i,j,k+2)-s33BY(i,j,k-1)))
   v3BY(i,j,k)=v3oBY(i,j,k)+v3pBY(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ################################################
subroutine vPMLTX
include 'triffy.dec'
do k=3,nz-2
 do j=3,ny-2
  do i=2,npm-1
   daf1x=DiTXS1(i,j,k)
   c01x=1+daf1x
   daf1y=DiTXS0(i,j,k)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   bo1=(1/rhoTX(i,j,k)+1/rhoTX(i+1,j,k))/2
   bo2=(1/rhoTX(i,j,k)+1/rhoTX(i,j-1,k))/2
   bo3=(1/rhoTX(i,j,k)+1/rhoTX(i,j,k+1))/2
!
   v1oTX(i,j,k)=(1-daf1x)/c01x*v1oTX(i,j,k)+bo1*(dtx/c01x)*   &
    (s11TX(i+1,j,k)-s11TX(i,j,k))
   v1pTX(i,j,k)=v1pTX(i,j,k)+bo1*dtx*(rn*(     &
    s12TX(i,j+1,k)-s12TX(i,j,k)+s13TX(i,j,k)-        & 
    s13TX(i,j,k-1))+rnn*(         &
    s12TX(i,j+2,k)-s12TX(i,j-1,k)+s13TX(i,j,k+1)-s13TX(i,j,k-2)))
   v1TX(i,j,k)=v1oTX(i,j,k)+v1pTX(i,j,k)
!
   v2oTX(i,j,k)=(1-daf1y)/c01y*v2oTX(i,j,k)+bo2*(dtx/c01y)*   &
    (s12TX(i,j,k)-s12TX(i-1,j,k))
   v2pTX(i,j,k)=v2pTX(i,j,k)+bo2*dtx*(rn*(       &
    s22TX(i,j,k)-s22TX(i,j-1,k)+s23TX(i,j,k)-      &   
    s23TX(i,j,k-1))+rnn*(   &
    s22TX(i,j+1,k)-s22TX(i,j-2,k)+s23TX(i,j,k+1)-s23TX(i,j,k-2)))
   v2TX(i,j,k)=v2oTX(i,j,k)+v2pTX(i,j,k)
!
   v3oTX(i,j,k)=(1-daf1z)/c01z*v3oTX(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TX(i,j,k)-s13TX(i-1,j,k))
   v3pTX(i,j,k)=v3pTX(i,j,k)+bo3*dtx*(rn*(    &
    s23TX(i,j+1,k)-s23TX(i,j,k)+s33TX(i,j,k+1)-    &
    s33TX(i,j,k))+rnn*(     &
    s23TX(i,j+2,k)-s23TX(i,j-1,k)+s33TX(i,j,k+2)-s33TX(i,j,k-1)))
   v3TX(i,j,k)=v3oTX(i,j,k)+v3pTX(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ################################################
subroutine vPMLBX
include 'triffy.dec'
do k=3,nz-2
 do j=3,ny-2
  do i=2,npm-1
   daf1x=DiBXS1(i,j,k)
   c01x=1+daf1x
   daf1y=DiBXS0(i,j,k)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   bo1=(1/rhoBX(i,j,k)+1/rhoBX(i+1,j,k))/2
   bo2=(1/rhoBX(i,j,k)+1/rhoBX(i,j-1,k))/2
   bo3=(1/rhoBX(i,j,k)+1/rhoBX(i,j,k+1))/2
!
   v1oBX(i,j,k)=(1-daf1x)/c01x*v1oBX(i,j,k)+bo1*(dtx/c01x)*   &
    (s11BX(i+1,j,k)-s11BX(i,j,k))
   v1pBX(i,j,k)=v1pBX(i,j,k)+bo1*dtx*(rn*(     &
    s12BX(i,j+1,k)-s12BX(i,j,k)+s13BX(i,j,k)-        & 
    s13BX(i,j,k-1))+rnn*(         &
    s12BX(i,j+2,k)-s12BX(i,j-1,k)+s13BX(i,j,k+1)-s13BX(i,j,k-2)))
   v1BX(i,j,k)=v1oBX(i,j,k)+v1pBX(i,j,k)
!
   v2oBX(i,j,k)=(1-daf1y)/c01y*v2oBX(i,j,k)+bo2*(dtx/c01y)*   &
    (s12BX(i,j,k)-s12BX(i-1,j,k))
   v2pBX(i,j,k)=v2pBX(i,j,k)+bo2*dtx*(rn*(       &
    s22BX(i,j,k)-s22BX(i,j-1,k)+s23BX(i,j,k)-      &   
    s23BX(i,j,k-1))+rnn*(   &
    s22BX(i,j+1,k)-s22BX(i,j-2,k)+s23BX(i,j,k+1)-s23BX(i,j,k-2)))
   v2BX(i,j,k)=v2oBX(i,j,k)+v2pBX(i,j,k)
!
   v3oBX(i,j,k)=(1-daf1z)/c01z*v3oBX(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BX(i,j,k)-s13BX(i-1,j,k))
   v3pBX(i,j,k)=v3pBX(i,j,k)+bo3*dtx*(rn*(    &
    s23BX(i,j+1,k)-s23BX(i,j,k)+s33BX(i,j,k+1)-    &
    s33BX(i,j,k))+rnn*(     &
    s23BX(i,j+2,k)-s23BX(i,j-1,k)+s33BX(i,j,k+2)-s33BX(i,j,k-1)))
   v3BX(i,j,k)=v3oBX(i,j,k)+v3pBX(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ################################################
Subroutine vPMLTXTY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm-1
   daf1x=DiTXS1(i,ny-2,k) 
   c01x=1+daf1x
   daf1y=DiTXS0(i,ny-2,k) 
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!   
   daf2x=DiTYS0(nx-2,j,k)
   c02x=1+daf2x
   daf2y=DiTYS1(nx-2,j,k)
   c02y=1+daf2y
   daf2z=daf2x
   c02z=c02x
      
   bo1=(1/rhoTXTY(i,j,k)+1/rhoTXTY(i+1,j,k))/2
   bo2=(1/rhoTXTY(i,j,k)+1/rhoTXTY(i,j-1,j))/2
   bo3=(1/rhoTXTY(i,j,k)+1/rhoTXTY(i,j,k+1))/2
!     
   v1xTXTY(i,j,k)=(1-daf1x)/c01x*v1xTXTY(i,j,k)+bo1*(dtx/c01x)*   &
    (s11TXTY(i+1,j,k)-s11TXTY(i,j,k))
   v1yTXTY(i,j,k)=(1-daf2x)/c02x*v1yTXTY(i,j,k)+bo1*(dtx/c02x)*   & 
    (s12TXTY(i,j+1,k)-s12TXTY(i,j,k)) 
   v1zTXTY(i,j,k)=v1zTXTY(i,j,k)+bo1*dtx*(rn*(s13TXTY(i,j,k)-     &
    s13TXTY(i,j,k-1))+rnn*(s13TXTY(i,j,k+1)-s13TXTY(i,j,k-2)))         
   v1TXTY(i,j,k)=v1xTXTY(i,j,k)+v1yTXTY(i,j,k)+v1zTXTY(i,j,k)
!
   v2xTXTY(i,j,k)=(1-daf1y)/c01y*v2xTXTY(i,j,k)+bo2*(dtx/c01y)*   &
    (s12TXTY(i,j,k)-s12TXTY(i-1,j,k))
   v2yTXTY(i,j,k)=(1-daf2y)/c02y*v2yTXTY(i,j,k)+bo2*(dtx/c02y)*   & 
    (s22TXTY(i,j,k)-s22TXTY(i,j-1,k)) 
   v2zTXTY(i,j,k)=v2zTXTY(i,j,k)+bo2*dtx*(rn*(s23TXTY(i,j,k)-       &
    s23TXTY(i,j,k-1))+rnn*(s23TXTY(i,j,k+1)-s23TXTY(i,j,k-2)))  
   v2TXTY(i,j,k)=v2xTXTY(i,j,k)+v2yTXTY(i,j,k)+v2zTXTY(i,j,k) 
 
!
   v3xTXTY(i,j,k)=(1-daf1z)/c01z*v3xTXTY(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TXTY(i,j,k)-s13TXTY(i-1,j,k))
   v3yTXTY(i,j,k)=(1-daf2z)/c02z*v3yTXTY(i,j,k)+bo3*(dtx/c02z)*   & 
    (s23TXTY(i,j+1,k)-s23TXTY(i,j,k)) 
   v3zTXTY(i,j,k)=v3zTXTY(i,j,k)+bo3*dtx*(rn*(s33TXTY(i,j,k+1)-    &
    s33TXTY(i,j,k))+rnn*(s33TXTY(i,j,k+2)-s33TXTY(i,j,k-1))) 
   v3TXTY(i,j,k)=v3xTXTY(i,j,k)+v3yTXTY(i,j,k)+v3zTXTY(i,j,k) 
  enddo
 enddo
enddo
return
end subroutine
! ################################################
Subroutine vPMLBXTY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm-1
   daf1x=DiBXS1(i,ny-2,k) 
   c01x=1+daf1x
   daf1y=DiBXS0(i,ny-2,k) 
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!   
   daf2x=DiTYS0(3,j,k)
   c02x=1+daf2x
   daf2y=DiTYS1(3,j,k)
   c02y=1+daf2y
   daf2z=daf2x
   c02z=c02x
      
   bo1=(1/rhoBXTY(i,j,k)+1/rhoBXTY(i+1,j,k))/2
   bo2=(1/rhoBXTY(i,j,k)+1/rhoBXTY(i,j-1,j))/2
   bo3=(1/rhoBXTY(i,j,k)+1/rhoBXTY(i,j,k+1))/2
!     
   v1xBXTY(i,j,k)=(1-daf1x)/c01x*v1xBXTY(i,j,k)+bo1*(dtx/c01x)*   &
    (s11BXTY(i+1,j,k)-s11BXTY(i,j,k))
   v1yBXTY(i,j,k)=(1-daf2x)/c02x*v1yBXTY(i,j,k)+bo1*(dtx/c02x)*   & 
    (s12BXTY(i,j+1,k)-s12BXTY(i,j,k)) 
   v1zBXTY(i,j,k)=v1zBXTY(i,j,k)+bo1*dtx*(rn*(s13BXTY(i,j,k)-     &
    s13BXTY(i,j,k-1))+rnn*(s13BXTY(i,j,k+1)-s13BXTY(i,j,k-2)))         
   v1BXTY(i,j,k)=v1xBXTY(i,j,k)+v1yBXTY(i,j,k)+v1zBXTY(i,j,k)
!
   v2xBXTY(i,j,k)=(1-daf1y)/c01y*v2xBXTY(i,j,k)+bo2*(dtx/c01y)*   &
    (s12BXTY(i,j,k)-s12BXTY(i-1,j,k))
   v2yBXTY(i,j,k)=(1-daf2y)/c02y*v2yBXTY(i,j,k)+bo2*(dtx/c02y)*   & 
    (s22BXTY(i,j,k)-s22BXTY(i,j-1,k)) 
   v2zBXTY(i,j,k)=v2zBXTY(i,j,k)+bo2*dtx*(rn*(s23BXTY(i,j,k)-       &
    s23BXTY(i,j,k-1))+rnn*(s23BXTY(i,j,k+1)-s23BXTY(i,j,k-2)))  
   v2BXTY(i,j,k)=v2xBXTY(i,j,k)+v2yBXTY(i,j,k)+v2zBXTY(i,j,k) 
 
!
   v3xBXTY(i,j,k)=(1-daf1z)/c01z*v3xBXTY(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BXTY(i,j,k)-s13BXTY(i-1,j,k))
   v3yBXTY(i,j,k)=(1-daf2z)/c02z*v3yBXTY(i,j,k)+bo3*(dtx/c02z)*   & 
    (s23BXTY(i,j+1,k)-s23BXTY(i,j,k)) 
   v3zBXTY(i,j,k)=v3zBXTY(i,j,k)+bo3*dtx*(rn*(s33BXTY(i,j,k+1)-    &
    s33BXTY(i,j,k))+rnn*(s33BXTY(i,j,k+2)-s33BXTY(i,j,k-1))) 
   v3BXTY(i,j,k)=v3xBXTY(i,j,k)+v3yBXTY(i,j,k)+v3zBXTY(i,j,k) 
  enddo
 enddo
enddo
return
end subroutine
! ################################################
Subroutine vPMLTXBY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm-1
   daf1x=DiTXS1(i,3,k) 
   c01x=1+daf1x
   daf1y=DiTXS0(i,3,k) 
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!   
   daf2x=DiBYS0(nx-2,j,k)
   c02x=1+daf2x
   daf2y=DiBYS1(nx-2,j,k)
   c02y=1+daf2y
   daf2z=daf2x
   c02z=c02x
      
   bo1=(1/rhoTXBY(i,j,k)+1/rhoTXBY(i+1,j,k))/2
   bo2=(1/rhoTXBY(i,j,k)+1/rhoTXBY(i,j-1,j))/2
   bo3=(1/rhoTXBY(i,j,k)+1/rhoTXBY(i,j,k+1))/2
!     
   v1xTXBY(i,j,k)=(1-daf1x)/c01x*v1xTXBY(i,j,k)+bo1*(dtx/c01x)*   &
    (s11TXBY(i+1,j,k)-s11TXBY(i,j,k))
   v1yTXBY(i,j,k)=(1-daf2x)/c02x*v1yTXBY(i,j,k)+bo1*(dtx/c02x)*   & 
    (s12TXBY(i,j+1,k)-s12TXBY(i,j,k)) 
   v1zTXBY(i,j,k)=v1zTXBY(i,j,k)+bo1*dtx*(rn*(s13TXBY(i,j,k)-     &
    s13TXBY(i,j,k-1))+rnn*(s13TXBY(i,j,k+1)-s13TXBY(i,j,k-2)))         
   v1TXBY(i,j,k)=v1xTXBY(i,j,k)+v1yTXBY(i,j,k)+v1zTXBY(i,j,k)
!
   v2xTXBY(i,j,k)=(1-daf1y)/c01y*v2xTXBY(i,j,k)+bo2*(dtx/c01y)*   &
    (s12TXBY(i,j,k)-s12TXBY(i-1,j,k))
   v2yTXBY(i,j,k)=(1-daf2y)/c02y*v2yTXBY(i,j,k)+bo2*(dtx/c02y)*   & 
    (s22TXBY(i,j,k)-s22TXBY(i,j-1,k)) 
   v2zTXBY(i,j,k)=v2zTXBY(i,j,k)+bo2*dtx*(rn*(s23TXBY(i,j,k)-       &
    s23TXBY(i,j,k-1))+rnn*(s23TXBY(i,j,k+1)-s23TXBY(i,j,k-2)))  
   v2TXBY(i,j,k)=v2xTXBY(i,j,k)+v2yTXBY(i,j,k)+v2zTXBY(i,j,k) 
 
!
   v3xTXBY(i,j,k)=(1-daf1z)/c01z*v3xTXBY(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TXBY(i,j,k)-s13TXBY(i-1,j,k))
   v3yTXBY(i,j,k)=(1-daf2z)/c02z*v3yTXBY(i,j,k)+bo3*(dtx/c02z)*   & 
    (s23TXBY(i,j+1,k)-s23TXBY(i,j,k)) 
   v3zTXBY(i,j,k)=v3zTXBY(i,j,k)+bo3*dtx*(rn*(s33TXBY(i,j,k+1)-    &
    s33TXBY(i,j,k))+rnn*(s33TXBY(i,j,k+2)-s33TXBY(i,j,k-1))) 
   v3TXBY(i,j,k)=v3xTXBY(i,j,k)+v3yTXBY(i,j,k)+v3zTXBY(i,j,k) 
  enddo
 enddo
enddo
return
end subroutine
! ################################################
Subroutine vPMLBXBY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm-1
   daf1x=DiBXS1(i,3,k) 
   c01x=1+daf1x
   daf1y=DiBXS0(i,3,k) 
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!   
   daf2x=DiBYS0(3,j,k)
   c02x=1+daf2x
   daf2y=DiBYS1(3,j,k)
   c02y=1+daf2y
   daf2z=daf2x
   c02z=c02x
      
   bo1=(1/rhoBXBY(i,j,k)+1/rhoBXBY(i+1,j,k))/2
   bo2=(1/rhoBXBY(i,j,k)+1/rhoBXBY(i,j-1,j))/2
   bo3=(1/rhoBXBY(i,j,k)+1/rhoBXBY(i,j,k+1))/2
!     
   v1xBXBY(i,j,k)=(1-daf1x)/c01x*v1xBXBY(i,j,k)+bo1*(dtx/c01x)*   &
    (s11BXBY(i+1,j,k)-s11BXBY(i,j,k))
   v1yBXBY(i,j,k)=(1-daf2x)/c02x*v1yBXBY(i,j,k)+bo1*(dtx/c02x)*   & 
    (s12BXBY(i,j+1,k)-s12BXBY(i,j,k)) 
   v1zBXBY(i,j,k)=v1zBXBY(i,j,k)+bo1*dtx*(rn*(s13BXBY(i,j,k)-     &
    s13BXBY(i,j,k-1))+rnn*(s13BXBY(i,j,k+1)-s13BXBY(i,j,k-2)))         
   v1BXBY(i,j,k)=v1xBXBY(i,j,k)+v1yBXBY(i,j,k)+v1zBXBY(i,j,k)
!
   v2xBXBY(i,j,k)=(1-daf1y)/c01y*v2xBXBY(i,j,k)+bo2*(dtx/c01y)*   &
    (s12BXBY(i,j,k)-s12BXBY(i-1,j,k))
   v2yBXBY(i,j,k)=(1-daf2y)/c02y*v2yBXBY(i,j,k)+bo2*(dtx/c02y)*   & 
    (s22BXBY(i,j,k)-s22BXBY(i,j-1,k)) 
   v2zBXBY(i,j,k)=v2zBXBY(i,j,k)+bo2*dtx*(rn*(s23BXBY(i,j,k)-       &
    s23BXBY(i,j,k-1))+rnn*(s23BXBY(i,j,k+1)-s23BXBY(i,j,k-2)))  
   v2BXBY(i,j,k)=v2xBXBY(i,j,k)+v2yBXBY(i,j,k)+v2zBXBY(i,j,k) 
 
!
   v3xBXBY(i,j,k)=(1-daf1z)/c01z*v3xBXBY(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BXBY(i,j,k)-s13BXBY(i-1,j,k))
   v3yBXBY(i,j,k)=(1-daf2z)/c02z*v3yBXBY(i,j,k)+bo3*(dtx/c02z)*   & 
    (s23BXBY(i,j+1,k)-s23BXBY(i,j,k)) 
   v3zBXBY(i,j,k)=v3zBXBY(i,j,k)+bo3*dtx*(rn*(s33BXBY(i,j,k+1)-    &
    s33BXBY(i,j,k))+rnn*(s33BXBY(i,j,k+2)-s33BXBY(i,j,k-1))) 
   v3BXBY(i,j,k)=v3xBXBY(i,j,k)+v3yBXBY(i,j,k)+v3zBXBY(i,j,k) 
  enddo
 enddo
enddo
return
end subroutine
! ################################################
subroutine vPMLTXTZ
include 'triffy.dec'
real c03x,c03y,c03z,c01x,c01y,c01z,daf1x,daf1y,daf1z
     
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm-1
!
   daf1x=DiTXS1(i,j,nz-2)
   c01x=1+daf1x  
   daf1y=DiTXS0(i,j,nz-2)
   c01y=1+daf1y  
   daf1z=daf1y
   c01z=c01y
!  
   daf3x=DiTZS0(nx-2,j,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(nx-2,j,k)
   c03z=1+daf3z
!
   bo1=(1/rhoTXTZ(i,j,k)+1/rhoTXTZ(i+1,j,k))/2
   bo2=(1/rhoTXTZ(i,j,k)+1/rhoTXTZ(i,j-1,k))/2
   bo3=(1/rhoTXTZ(i,j,k)+1/rhoTXTZ(i,j,k+1))/2
!
   v1xTXTZ(i,j,k)=(1-daf1x)/c01x*v1xTXTZ(i,j,k)+bo1*(dtx/c01x)*   &  
    (s11TXTZ(i+1,j,k)-s11TXTZ(i,j,k))  
   v1yTXTZ(i,j,k)=v1yTXTZ(i,j,k)+bo1*dtx*(rn*(s12TXTZ(i,j+1,k)-    & 
    s12TXTZ(i,j,k))+rnn*(s12TXTZ(i,j+2,k)-s12TXTZ(i,j-1,k))) 
   v1zTXTZ(i,j,k)=(1-daf3x)/c03x*v1zTXTZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13TXTZ(i,j,k)-s13TXTZ(i,j,k-1))    
   v1TXTZ(i,j,k)=v1xTXTZ(i,j,k)+v1yTXTZ(i,j,k)+v1zTXTZ(i,j,k)
!
   v2xTXTZ(i,j,k)=(1-daf1y)/c01y*v2xTXTZ(i,j,k)+bo2*(dtx/c01y)*   & 
    (s12TXTZ(i,j,k)-s12TXTZ(i-1,j,k))    
   v2yTXTZ(i,j,k)=v2yTXTZ(i,j,k)+bo2*dtx*(rn*(s22TXTZ(i,j,k)-   & 
    s22TXTZ(i,j-1,k))+rnn*(s22TXTZ(i,j+1,k)-s22TXTZ(i,j-2,k))) 
   v2zTXTZ(i,j,k)=(1-daf3y)/c03y*v2zTXTZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23TXTZ(i,j,k)-s23TXTZ(i,j,k-1)) 
   v2TXTZ(i,j,k)=v2xTXTZ(i,j,k)+v2yTXTZ(i,j,k)+v2zTXTZ(i,j,k) 
!
   v3xTXTZ(i,j,k)=(1-daf1z)/c01z*v3xTXTZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TXTZ(i,j,k)-s13TXTZ(i-1,j,k))    
   v3yTXTZ(i,j,k)=v3yTXTZ(i,j,k)+bo3*dtx*(rn*(s23TXTZ(i,j+1,k)-     & 
    s23TXTZ(i,j,k))+rnn*(s23TXTZ(i,j+2,k)-s23TXTZ(i,j-1,k)))
   v3zTXTZ(i,j,k)=(1-daf3z)/c03z*v3zTXTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TXTZ(i,j,k+1)-s33TXTZ(i,j,k)) 
   v3TXTZ(i,j,k)=v3xTXTZ(i,j,k)+v3yTXTZ(i,j,k)+v3zTXTZ(i,j,k) 
!
  enddo
 enddo
enddo
return
end subroutine
! ################################################
subroutine vPMLBXTZ
include 'triffy.dec'
real c03x,c03y,c03z,c01x,c01y,c01z,daf1x,daf1y,daf1z
     
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm-1
!
   daf1x=DiBXS1(i,j,nz-2)
   c01x=1+daf1x  
   daf1y=DiBXS0(i,j,nz-2)
   c01y=1+daf1y  
   daf1z=daf1y
   c01z=c01y
!  
   daf3x=DiTZS0(3,j,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(3,j,k)
   c03z=1+daf3z
!
   bo1=(1/rhoBXTZ(i,j,k)+1/rhoBXTZ(i+1,j,k))/2
   bo2=(1/rhoBXTZ(i,j,k)+1/rhoBXTZ(i,j-1,k))/2
   bo3=(1/rhoBXTZ(i,j,k)+1/rhoBXTZ(i,j,k+1))/2
!
   v1xBXTZ(i,j,k)=(1-daf1x)/c01x*v1xBXTZ(i,j,k)+bo1*(dtx/c01x)*   &  
    (s11BXTZ(i+1,j,k)-s11BXTZ(i,j,k))  
   v1yBXTZ(i,j,k)=v1yBXTZ(i,j,k)+bo1*dtx*(rn*(s12BXTZ(i,j+1,k)-    & 
    s12BXTZ(i,j,k))+rnn*(s12BXTZ(i,j+2,k)-s12BXTZ(i,j-1,k))) 
   v1zBXTZ(i,j,k)=(1-daf3x)/c03x*v1zBXTZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13BXTZ(i,j,k)-s13BXTZ(i,j,k-1))    
   v1BXTZ(i,j,k)=v1xBXTZ(i,j,k)+v1yBXTZ(i,j,k)+v1zBXTZ(i,j,k)
!
   v2xBXTZ(i,j,k)=(1-daf1y)/c01y*v2xBXTZ(i,j,k)+bo2*(dtx/c01y)*   & 
    (s12BXTZ(i,j,k)-s12BXTZ(i-1,j,k))    
   v2yBXTZ(i,j,k)=v2yBXTZ(i,j,k)+bo2*dtx*(rn*(s22BXTZ(i,j,k)-   & 
    s22BXTZ(i,j-1,k))+rnn*(s22BXTZ(i,j+1,k)-s22BXTZ(i,j-2,k))) 
   v2zBXTZ(i,j,k)=(1-daf3y)/c03y*v2zBXTZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23BXTZ(i,j,k)-s23BXTZ(i,j,k-1)) 
   v2BXTZ(i,j,k)=v2xBXTZ(i,j,k)+v2yBXTZ(i,j,k)+v2zBXTZ(i,j,k) 
!
   v3xBXTZ(i,j,k)=(1-daf1z)/c01z*v3xBXTZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BXTZ(i,j,k)-s13BXTZ(i-1,j,k))    
   v3yBXTZ(i,j,k)=v3yBXTZ(i,j,k)+bo3*dtx*(rn*(s23BXTZ(i,j+1,k)-     & 
    s23BXTZ(i,j,k))+rnn*(s23BXTZ(i,j+2,k)-s23BXTZ(i,j-1,k)))
   v3zBXTZ(i,j,k)=(1-daf3z)/c03z*v3zBXTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BXTZ(i,j,k+1)-s33BXTZ(i,j,k)) 
   v3BXTZ(i,j,k)=v3xBXTZ(i,j,k)+v3yBXTZ(i,j,k)+v3zBXTZ(i,j,k) 
!
  enddo
 enddo
enddo
return
end subroutine
! ################################################
subroutine vPMLTXBZ
include 'triffy.dec'
real c03x,c03y,c03z,c01x,c01y,c01z,daf1x,daf1y,daf1z
     
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm-1
!
   daf1x=DiTXS1(i,j,3)
   c01x=1+daf1x  
   daf1y=DiTXS0(i,j,3)
   c01y=1+daf1y  
   daf1z=daf1y
   c01z=c01y
!  
   daf3x=DiBZS0(nx-2,j,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(nx-2,j,k)
   c03z=1+daf3z
!
   bo1=(1/rhoTXBZ(i,j,k)+1/rhoTXBZ(i+1,j,k))/2
   bo2=(1/rhoTXBZ(i,j,k)+1/rhoTXBZ(i,j-1,k))/2
   bo3=(1/rhoTXBZ(i,j,k)+1/rhoTXBZ(i,j,k+1))/2
!
   v1xTXBZ(i,j,k)=(1-daf1x)/c01x*v1xTXBZ(i,j,k)+bo1*(dtx/c01x)*   &  
    (s11TXBZ(i+1,j,k)-s11TXBZ(i,j,k))  
   v1yTXBZ(i,j,k)=v1yTXBZ(i,j,k)+bo1*dtx*(rn*(s12TXBZ(i,j+1,k)-    & 
    s12TXBZ(i,j,k))+rnn*(s12TXBZ(i,j+2,k)-s12TXBZ(i,j-1,k))) 
   v1zTXBZ(i,j,k)=(1-daf3x)/c03x*v1zTXBZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13TXBZ(i,j,k)-s13TXBZ(i,j,k-1))    
   v1TXBZ(i,j,k)=v1xTXBZ(i,j,k)+v1yTXBZ(i,j,k)+v1zTXBZ(i,j,k)
!
   v2xTXBZ(i,j,k)=(1-daf1y)/c01y*v2xTXBZ(i,j,k)+bo2*(dtx/c01y)*   & 
    (s12TXBZ(i,j,k)-s12TXBZ(i-1,j,k))    
   v2yTXBZ(i,j,k)=v2yTXBZ(i,j,k)+bo2*dtx*(rn*(s22TXBZ(i,j,k)-   & 
    s22TXBZ(i,j-1,k))+rnn*(s22TXBZ(i,j+1,k)-s22TXBZ(i,j-2,k))) 
   v2zTXBZ(i,j,k)=(1-daf3y)/c03y*v2zTXBZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23TXBZ(i,j,k)-s23TXBZ(i,j,k-1)) 
   v2TXBZ(i,j,k)=v2xTXBZ(i,j,k)+v2yTXBZ(i,j,k)+v2zTXBZ(i,j,k) 
!
   v3xTXBZ(i,j,k)=(1-daf1z)/c01z*v3xTXBZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TXBZ(i,j,k)-s13TXBZ(i-1,j,k))    
   v3yTXBZ(i,j,k)=v3yTXBZ(i,j,k)+bo3*dtx*(rn*(s23TXBZ(i,j+1,k)-     & 
    s23TXBZ(i,j,k))+rnn*(s23TXBZ(i,j+2,k)-s23TXBZ(i,j-1,k)))
   v3zTXBZ(i,j,k)=(1-daf3z)/c03z*v3zTXBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TXBZ(i,j,k+1)-s33TXBZ(i,j,k)) 
   v3TXBZ(i,j,k)=v3xTXBZ(i,j,k)+v3yTXBZ(i,j,k)+v3zTXBZ(i,j,k) 
!
  enddo
 enddo
enddo
!print *, v3zTXBZ(3,ny-2,7)
return
end subroutine
! ################################################
subroutine vPMLBXBZ
include 'triffy.dec'
real c03x,c03y,c03z,c01x,c01y,c01z,daf1x,daf1y,daf1z
     
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm-1
!
   daf1x=DiBXS1(i,j,3)
   c01x=1+daf1x  
   daf1y=DiBXS0(i,j,3)
   c01y=1+daf1y  
   daf1z=daf1y
   c01z=c01y
!  
   daf3x=DiBZS0(3,j,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(3,j,k)
   c03z=1+daf3z
!
   bo1=(1/rhoBXBZ(i,j,k)+1/rhoBXBZ(i+1,j,k))/2
   bo2=(1/rhoBXBZ(i,j,k)+1/rhoBXBZ(i,j-1,k))/2
   bo3=(1/rhoBXBZ(i,j,k)+1/rhoBXBZ(i,j,k+1))/2
!
   v1xBXBZ(i,j,k)=(1-daf1x)/c01x*v1xBXBZ(i,j,k)+bo1*(dtx/c01x)*   &  
    (s11BXBZ(i+1,j,k)-s11BXBZ(i,j,k))  
   v1yBXBZ(i,j,k)=v1yBXBZ(i,j,k)+bo1*dtx*(rn*(s12BXBZ(i,j+1,k)-    & 
    s12BXBZ(i,j,k))+rnn*(s12BXBZ(i,j+2,k)-s12BXBZ(i,j-1,k))) 
   v1zBXBZ(i,j,k)=(1-daf3x)/c03x*v1zBXBZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13BXBZ(i,j,k)-s13BXBZ(i,j,k-1))    
   v1BXBZ(i,j,k)=v1xBXBZ(i,j,k)+v1yBXBZ(i,j,k)+v1zBXBZ(i,j,k)
!
   v2xBXBZ(i,j,k)=(1-daf1y)/c01y*v2xBXBZ(i,j,k)+bo2*(dtx/c01y)*   & 
    (s12BXBZ(i,j,k)-s12BXBZ(i-1,j,k))    
   v2yBXBZ(i,j,k)=v2yBXBZ(i,j,k)+bo2*dtx*(rn*(s22BXBZ(i,j,k)-   & 
    s22BXBZ(i,j-1,k))+rnn*(s22BXBZ(i,j+1,k)-s22BXBZ(i,j-2,k))) 
   v2zBXBZ(i,j,k)=(1-daf3y)/c03y*v2zBXBZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23BXBZ(i,j,k)-s23BXBZ(i,j,k-1)) 
   v2BXBZ(i,j,k)=v2xBXBZ(i,j,k)+v2yBXBZ(i,j,k)+v2zBXBZ(i,j,k) 
!
   v3xBXBZ(i,j,k)=(1-daf1z)/c01z*v3xBXBZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BXBZ(i,j,k)-s13BXBZ(i-1,j,k))    
   v3yBXBZ(i,j,k)=v3yBXBZ(i,j,k)+bo3*dtx*(rn*(s23BXBZ(i,j+1,k)-     & 
    s23BXBZ(i,j,k))+rnn*(s23BXBZ(i,j+2,k)-s23BXBZ(i,j-1,k)))
   v3zBXBZ(i,j,k)=(1-daf3z)/c03z*v3zBXBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BXBZ(i,j,k+1)-s33BXBZ(i,j,k)) 
   v3BXBZ(i,j,k)=v3xBXBZ(i,j,k)+v3yBXBZ(i,j,k)+v3zBXBZ(i,j,k) 
!
  enddo
 enddo
enddo
return
end subroutine
! ############################################### 
subroutine vPMLTYTZ 
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=3,nx-2
! 
   daf2x=DiTYS0(i,j,nz-2)
   c02x=1+daf2x 
   daf2y=DiTYS1(i,j,nz-2)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
!  
   daf3x=DiTZS0(i,ny-2,k)
   c03x=1+daf3x 
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(i,ny-2,k)
   c03z=1+daf3z 
! 
   bo1=(1/rhoTYTZ(i,j,k)+1/rhoTYTZ(i+1,j,k))/2
   bo2=(1/rhoTYTZ(i,j,k)+1/rhoTYTZ(i,j-1,k))/2 
   bo3=(1/rhoTYTZ(i,j,k)+1/rhoTYTZ(i,j,k+1))/2 
! 
   v1xTYTZ(i,j,k)=v1xTYTZ(i,j,k)+bo1*dtx*(rn*(s11TYTZ(i+1,j,k)-    &  
    s11TYTZ(i,j,k))+rnn*(s11TYTZ(i+2,j,k)-s11TYTZ(i-1,j,k)))  
   v1yTYTZ(i,j,k)=(1-daf2x)/c02x*v1yTYTZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12TYTZ(i,j+1,k)-s12TYTZ(i,j,k))   
   v1zTYTZ(i,j,k)=(1-daf3x)/c03x*v1zTYTZ(i,j,k)+bo1*(dtx/c03x)*    & 
    (s13TYTZ(i,j,k)-s13TYTZ(i,j,k-1))      
   v1TYTZ(i,j,k)=v1xTYTZ(i,j,k)+v1yTYTZ(i,j,k)+v1zTYTZ(i,j,k)  
!   
   v2xTYTZ(i,j,k)=v2xTYTZ(i,j,k)+bo2*dtx*(rn*(s12TYTZ(i,j,k)-     & 
    s12TYTZ(i-1,j,k))+rnn*(s12TYTZ(i+1,j,k)-s12TYTZ(i-2,j,k)))  
   v2yTYTZ(i,j,k)=(1-daf2y)/c02y*v2yTYTZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22TYTZ(i,j,k)-s22TYTZ(i,j-1,k))  
   v2zTYTZ(i,j,k)=(1-daf3y)/c03y*v2zTYTZ(i,j,k)+bo2*(dtx/c03y)*  & 
    (s23TYTZ(i,j,k)-s23TYTZ(i,j,k-1))    
   v2TYTZ(i,j,k)=v2xTYTZ(i,j,k)+v2yTYTZ(i,j,k)+v2zTYTZ(i,j,k) 
!   
   v3xTYTZ(i,j,k)=v3xTYTZ(i,j,k)+bo3*dtx*(rn*(s13TYTZ(i,j,k)-       &  
    s13TYTZ(i-1,j,k))+rnn*(s13TYTZ(i+1,j,k)-s13TYTZ(i-2,j,k)))  
   v3yTYTZ(i,j,k)=(1-daf2z)/c02z*v3yTYTZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23TYTZ(i,j+1,k)-s23TYTZ(i,j,k))
   v3zTYTZ(i,j,k)=(1-daf3z)/c03z*v3zTYTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TYTZ(i,j,k+1)-s33TYTZ(i,j,k))
   v3TYTZ(i,j,k)=v3xTYTZ(i,j,k)+v3yTYTZ(i,j,k)+v3zTYTZ(i,j,k)
  enddo 
 enddo  
enddo 
return  
end subroutine 
! ############################################### 
subroutine vPMLBYTZ 
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=3,nx-2
! 
   daf2x=DiBYS0(i,j,nz-2)
   c02x=1+daf2x 
   daf2y=DiBYS1(i,j,nz-2)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
!  
   daf3x=DiTZS0(i,3,k)
   c03x=1+daf3x 
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(i,3,k)
   c03z=1+daf3z 
! 
   bo1=(1/rhoBYTZ(i,j,k)+1/rhoBYTZ(i+1,j,k))/2
   bo2=(1/rhoBYTZ(i,j,k)+1/rhoBYTZ(i,j-1,k))/2 
   bo3=(1/rhoBYTZ(i,j,k)+1/rhoBYTZ(i,j,k+1))/2 
! 
   v1xBYTZ(i,j,k)=v1xBYTZ(i,j,k)+bo1*dtx*(rn*(s11BYTZ(i+1,j,k)-    &  
    s11BYTZ(i,j,k))+rnn*(s11BYTZ(i+2,j,k)-s11BYTZ(i-1,j,k)))  
   v1yBYTZ(i,j,k)=(1-daf2x)/c02x*v1yBYTZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12BYTZ(i,j+1,k)-s12BYTZ(i,j,k))   
   v1zBYTZ(i,j,k)=(1-daf3x)/c03x*v1zBYTZ(i,j,k)+bo1*(dtx/c03x)*    & 
    (s13BYTZ(i,j,k)-s13BYTZ(i,j,k-1))      
   v1BYTZ(i,j,k)=v1xBYTZ(i,j,k)+v1yBYTZ(i,j,k)+v1zBYTZ(i,j,k)  
!   
   v2xBYTZ(i,j,k)=v2xBYTZ(i,j,k)+bo2*dtx*(rn*(s12BYTZ(i,j,k)-     & 
    s12BYTZ(i-1,j,k))+rnn*(s12BYTZ(i+1,j,k)-s12BYTZ(i-2,j,k)))  
   v2yBYTZ(i,j,k)=(1-daf2y)/c02y*v2yBYTZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22BYTZ(i,j,k)-s22BYTZ(i,j-1,k))  
   v2zBYTZ(i,j,k)=(1-daf3y)/c03y*v2zBYTZ(i,j,k)+bo2*(dtx/c03y)*  & 
    (s23BYTZ(i,j,k)-s23BYTZ(i,j,k-1))    
   v2BYTZ(i,j,k)=v2xBYTZ(i,j,k)+v2yBYTZ(i,j,k)+v2zBYTZ(i,j,k) 
!   
   v3xBYTZ(i,j,k)=v3xBYTZ(i,j,k)+bo3*dtx*(rn*(s13BYTZ(i,j,k)-       &  
    s13BYTZ(i-1,j,k))+rnn*(s13BYTZ(i+1,j,k)-s13BYTZ(i-2,j,k)))  
   v3yBYTZ(i,j,k)=(1-daf2z)/c02z*v3yBYTZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23BYTZ(i,j+1,k)-s23BYTZ(i,j,k))
   v3zBYTZ(i,j,k)=(1-daf3z)/c03z*v3zBYTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BYTZ(i,j,k+1)-s33BYTZ(i,j,k))
   v3BYTZ(i,j,k)=v3xBYTZ(i,j,k)+v3yBYTZ(i,j,k)+v3zBYTZ(i,j,k)
  enddo 
 enddo  
enddo 
return  
end subroutine 
! ############################################### 
subroutine vPMLTYBZ 
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=3,nx-2
! 
   daf2x=DiTYS0(i,j,3)
   c02x=1+daf2x 
   daf2y=DiTYS1(i,j,3)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
!  
   daf3x=DiBZS0(i,ny-2,k)
   c03x=1+daf3x 
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(i,ny-2,k)
   c03z=1+daf3z 
! 
   bo1=(1/rhoTYBZ(i,j,k)+1/rhoTYBZ(i+1,j,k))/2
   bo2=(1/rhoTYBZ(i,j,k)+1/rhoTYBZ(i,j-1,k))/2 
   bo3=(1/rhoTYBZ(i,j,k)+1/rhoTYBZ(i,j,k+1))/2 
! 
   v1xTYBZ(i,j,k)=v1xTYBZ(i,j,k)+bo1*dtx*(rn*(s11TYBZ(i+1,j,k)-    &  
    s11TYBZ(i,j,k))+rnn*(s11TYBZ(i+2,j,k)-s11TYBZ(i-1,j,k)))  
   v1yTYBZ(i,j,k)=(1-daf2x)/c02x*v1yTYBZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12TYBZ(i,j+1,k)-s12TYBZ(i,j,k))   
   v1zTYBZ(i,j,k)=(1-daf3x)/c03x*v1zTYBZ(i,j,k)+bo1*(dtx/c03x)*    & 
    (s13TYBZ(i,j,k)-s13TYBZ(i,j,k-1))      
   v1TYBZ(i,j,k)=v1xTYBZ(i,j,k)+v1yTYBZ(i,j,k)+v1zTYBZ(i,j,k)  
!   
   v2xTYBZ(i,j,k)=v2xTYBZ(i,j,k)+bo2*dtx*(rn*(s12TYBZ(i,j,k)-     & 
    s12TYBZ(i-1,j,k))+rnn*(s12TYBZ(i+1,j,k)-s12TYBZ(i-2,j,k)))  
   v2yTYBZ(i,j,k)=(1-daf2y)/c02y*v2yTYBZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22TYBZ(i,j,k)-s22TYBZ(i,j-1,k))  
   v2zTYBZ(i,j,k)=(1-daf3y)/c03y*v2zTYBZ(i,j,k)+bo2*(dtx/c03y)*  & 
    (s23TYBZ(i,j,k)-s23TYBZ(i,j,k-1))    
   v2TYBZ(i,j,k)=v2xTYBZ(i,j,k)+v2yTYBZ(i,j,k)+v2zTYBZ(i,j,k) 
!   
   v3xTYBZ(i,j,k)=v3xTYBZ(i,j,k)+bo3*dtx*(rn*(s13TYBZ(i,j,k)-       &  
    s13TYBZ(i-1,j,k))+rnn*(s13TYBZ(i+1,j,k)-s13TYBZ(i-2,j,k)))  
   v3yTYBZ(i,j,k)=(1-daf2z)/c02z*v3yTYBZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23TYBZ(i,j+1,k)-s23TYBZ(i,j,k))
   v3zTYBZ(i,j,k)=(1-daf3z)/c03z*v3zTYBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TYBZ(i,j,k+1)-s33TYBZ(i,j,k))
   v3TYBZ(i,j,k)=v3xTYBZ(i,j,k)+v3yTYBZ(i,j,k)+v3zTYBZ(i,j,k)
  enddo 
 enddo  
enddo 
return  
end subroutine 
! ############################################### 
subroutine vPMLBYBZ 
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=3,nx-2
! 
   daf2x=DiBYS0(i,j,3)
   c02x=1+daf2x 
   daf2y=DiBYS1(i,j,3)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
!  
   daf3x=DiBZS0(i,3,k)
   c03x=1+daf3x 
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(i,3,k)
   c03z=1+daf3z 
! 
   bo1=(1/rhoBYBZ(i,j,k)+1/rhoBYBZ(i+1,j,k))/2
   bo2=(1/rhoBYBZ(i,j,k)+1/rhoBYBZ(i,j-1,k))/2 
   bo3=(1/rhoBYBZ(i,j,k)+1/rhoBYBZ(i,j,k+1))/2 
! 
   v1xBYBZ(i,j,k)=v1xBYBZ(i,j,k)+bo1*dtx*(rn*(s11BYBZ(i+1,j,k)-    &  
    s11BYBZ(i,j,k))+rnn*(s11BYBZ(i+2,j,k)-s11BYBZ(i-1,j,k)))  
   v1yBYBZ(i,j,k)=(1-daf2x)/c02x*v1yBYBZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12BYBZ(i,j+1,k)-s12BYBZ(i,j,k))   
   v1zBYBZ(i,j,k)=(1-daf3x)/c03x*v1zBYBZ(i,j,k)+bo1*(dtx/c03x)*    & 
    (s13BYBZ(i,j,k)-s13BYBZ(i,j,k-1))      
   v1BYBZ(i,j,k)=v1xBYBZ(i,j,k)+v1yBYBZ(i,j,k)+v1zBYBZ(i,j,k)  
!   
   v2xBYBZ(i,j,k)=v2xBYBZ(i,j,k)+bo2*dtx*(rn*(s12BYBZ(i,j,k)-     & 
    s12BYBZ(i-1,j,k))+rnn*(s12BYBZ(i+1,j,k)-s12BYBZ(i-2,j,k)))  
   v2yBYBZ(i,j,k)=(1-daf2y)/c02y*v2yBYBZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22BYBZ(i,j,k)-s22BYBZ(i,j-1,k))  
   v2zBYBZ(i,j,k)=(1-daf3y)/c03y*v2zBYBZ(i,j,k)+bo2*(dtx/c03y)*  & 
    (s23BYBZ(i,j,k)-s23BYBZ(i,j,k-1))    
   v2BYBZ(i,j,k)=v2xBYBZ(i,j,k)+v2yBYBZ(i,j,k)+v2zBYBZ(i,j,k) 
!   
   v3xBYBZ(i,j,k)=v3xBYBZ(i,j,k)+bo3*dtx*(rn*(s13BYBZ(i,j,k)-       &  
    s13BYBZ(i-1,j,k))+rnn*(s13BYBZ(i+1,j,k)-s13BYBZ(i-2,j,k)))  
   v3yBYBZ(i,j,k)=(1-daf2z)/c02z*v3yBYBZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23BYBZ(i,j+1,k)-s23BYBZ(i,j,k))
   v3zBYBZ(i,j,k)=(1-daf3z)/c03z*v3zBYBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BYBZ(i,j,k+1)-s33BYBZ(i,j,k))
   v3BYBZ(i,j,k)=v3xBYBZ(i,j,k)+v3yBYBZ(i,j,k)+v3zBYBZ(i,j,k)
  enddo 
 enddo  
enddo 
return  
end subroutine 
! ##########################################################
Subroutine vPMLTXTYTZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
!
   daf1x=DiTXS1(i,ny-2,nz-2)
   c01x=1+daf1x
   daf1y=DiTXS0(i,ny-2,nz-2)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   daf2x=DiTYS0(nx-2,j,nz-2)
   c02x=1+daf2x
   daf2y=DiTYS1(nx-2,j,nz-2)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
! 
   daf3x=DiTZS0(nx-2,ny-2,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(nx-2,ny-2,k)
   c03z=1+daf3z
!
   bo1=(1/rhoTXTYTZ(i,j,k)+1/rhoTXTYTZ(i+1,j,k))/2
   bo2=(1/rhoTXTYTZ(i,j,k)+1/rhoTXTYTZ(i,j-1,k))/2
   bo3=(1/rhoTXTYTZ(i,j,k)+1/rhoTXTYTZ(i,j,k+1))/2
!
   v1xTXTYTZ(i,j,k)=(1-daf1x)/c01x*v1xTXTYTZ(i,j,k)+bo1*(dtx/c01x)*   &
    (s11TXTYTZ(i+1,j,k)-s11TXTYTZ(i,j,k))
   v1yTXTYTZ(i,j,k)=(1-daf2x)/c02x*v1yTXTYTZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12TXTYTZ(i,j+1,k)-s12TXTYTZ(i,j,k))   
   v1zTXTYTZ(i,j,k)=(1-daf3x)/c03x*v1zTXTYTZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13TXTYTZ(i,j,k)-s13TXTYTZ(i,j,k-1))  
   v1TXTYTZ(i,j,k)=v1xTXTYTZ(i,j,k)+v1yTXTYTZ(i,j,k)+v1zTXTYTZ(i,j,k)
!
   v2xTXTYTZ(i,j,k)=(1-daf1y)/c01y*v2xTXTYTZ(i,j,k)+bo2*(dtx/c01y)*   &
    (s12TXTYTZ(i,j,k)-s12TXTYTZ(i-1,j,k))    
   v2yTXTYTZ(i,j,k)=(1-daf2y)/c02y*v2yTXTYTZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22TXTYTZ(i,j,k)-s22TXTYTZ(i,j-1,k))    
   v2zTXTYTZ(i,j,k)=(1-daf3y)/c03y*v2zTXTYTZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23TXTYTZ(i,j,k)-s23TXTYTZ(i,j,k-1))   
   v2TXTYTZ(i,j,k)=v2xTXTYTZ(i,j,k)+v2yTXTYTZ(i,j,k)+v2zTXTYTZ(i,j,k)
!
   v3xTXTYTZ(i,j,k)=(1-daf1z)/c01z*v3xTXTYTZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TXTYTZ(i,j,k)-s13TXTYTZ(i-1,j,k))
   v3yTXTYTZ(i,j,k)=(1-daf2z)/c02z*v3yTXTYTZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23TXTYTZ(i,j+1,k)-s23TXTYTZ(i,j,k))  
   v3zTXTYTZ(i,j,k)=(1-daf3z)/c03z*v3zTXTYTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TXTYTZ(i,j,k+1)-s33TXTYTZ(i,j,k))   
   v3TXTYTZ(i,j,k)=v3xTXTYTZ(i,j,k)+v3yTXTYTZ(i,j,k)+v3zTXTYTZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine

! ##########################################################
Subroutine vPMLTXTYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
!
   daf1x=DiTXS1(i,ny-2,3)
   c01x=1+daf1x
   daf1y=DiTXS0(i,ny-2,3)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   daf2x=DiTYS0(nx-2,j,3)
   c02x=1+daf2x
   daf2y=DiTYS1(nx-2,j,3)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
! 
   daf3x=DiBZS0(nx-2,ny-2,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(nx-2,ny-2,k)
   c03z=1+daf3z
!
   bo1=(1/rhoTXTYBZ(i,j,k)+1/rhoTXTYBZ(i+1,j,k))/2
   bo2=(1/rhoTXTYBZ(i,j,k)+1/rhoTXTYBZ(i,j-1,k))/2
   bo3=(1/rhoTXTYBZ(i,j,k)+1/rhoTXTYBZ(i,j,k+1))/2
!
   v1xTXTYBZ(i,j,k)=(1-daf1x)/c01x*v1xTXTYBZ(i,j,k)+bo1*(dtx/c01x)*   &
    (s11TXTYBZ(i+1,j,k)-s11TXTYBZ(i,j,k))
   v1yTXTYBZ(i,j,k)=(1-daf2x)/c02x*v1yTXTYBZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12TXTYBZ(i,j+1,k)-s12TXTYBZ(i,j,k))   
   v1zTXTYBZ(i,j,k)=(1-daf3x)/c03x*v1zTXTYBZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13TXTYBZ(i,j,k)-s13TXTYBZ(i,j,k-1))  
   v1TXTYBZ(i,j,k)=v1xTXTYBZ(i,j,k)+v1yTXTYBZ(i,j,k)+v1zTXTYBZ(i,j,k)
!
   v2xTXTYBZ(i,j,k)=(1-daf1y)/c01y*v2xTXTYBZ(i,j,k)+bo2*(dtx/c01y)*   &
    (s12TXTYBZ(i,j,k)-s12TXTYBZ(i-1,j,k))    
   v2yTXTYBZ(i,j,k)=(1-daf2y)/c02y*v2yTXTYBZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22TXTYBZ(i,j,k)-s22TXTYBZ(i,j-1,k))    
   v2zTXTYBZ(i,j,k)=(1-daf3y)/c03y*v2zTXTYBZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23TXTYBZ(i,j,k)-s23TXTYBZ(i,j,k-1))   
   v2TXTYBZ(i,j,k)=v2xTXTYBZ(i,j,k)+v2yTXTYBZ(i,j,k)+v2zTXTYBZ(i,j,k)
!
   v3xTXTYBZ(i,j,k)=(1-daf1z)/c01z*v3xTXTYBZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TXTYBZ(i,j,k)-s13TXTYBZ(i-1,j,k))
   v3yTXTYBZ(i,j,k)=(1-daf2z)/c02z*v3yTXTYBZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23TXTYBZ(i,j+1,k)-s23TXTYBZ(i,j,k))  
   v3zTXTYBZ(i,j,k)=(1-daf3z)/c03z*v3zTXTYBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TXTYBZ(i,j,k+1)-s33TXTYBZ(i,j,k))   
   v3TXTYBZ(i,j,k)=v3xTXTYBZ(i,j,k)+v3yTXTYBZ(i,j,k)+v3zTXTYBZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ##########################################################
Subroutine vPMLTXBYTZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
!
   daf1x=DiTXS1(i,3,nz-2)
   c01x=1+daf1x
   daf1y=DiTXS0(i,3,nz-2)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   daf2x=DiBYS0(nx-2,j,nz-2)
   c02x=1+daf2x
   daf2y=DiBYS1(nx-2,j,nz-2)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
! 
   daf3x=DiTZS0(nx-2,3,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(nx-2,3,k)
   c03z=1+daf3z
!
   bo1=(1/rhoTXBYTZ(i,j,k)+1/rhoTXBYTZ(i+1,j,k))/2
   bo2=(1/rhoTXBYTZ(i,j,k)+1/rhoTXBYTZ(i,j-1,k))/2
   bo3=(1/rhoTXBYTZ(i,j,k)+1/rhoTXBYTZ(i,j,k+1))/2
!
   v1xTXBYTZ(i,j,k)=(1-daf1x)/c01x*v1xTXBYTZ(i,j,k)+bo1*(dtx/c01x)*   &
    (s11TXBYTZ(i+1,j,k)-s11TXBYTZ(i,j,k))
   v1yTXBYTZ(i,j,k)=(1-daf2x)/c02x*v1yTXBYTZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12TXBYTZ(i,j+1,k)-s12TXBYTZ(i,j,k))   
   v1zTXBYTZ(i,j,k)=(1-daf3x)/c03x*v1zTXBYTZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13TXBYTZ(i,j,k)-s13TXBYTZ(i,j,k-1))  
   v1TXBYTZ(i,j,k)=v1xTXBYTZ(i,j,k)+v1yTXBYTZ(i,j,k)+v1zTXBYTZ(i,j,k)
!
   v2xTXBYTZ(i,j,k)=(1-daf1y)/c01y*v2xTXBYTZ(i,j,k)+bo2*(dtx/c01y)*   &
    (s12TXBYTZ(i,j,k)-s12TXBYTZ(i-1,j,k))    
   v2yTXBYTZ(i,j,k)=(1-daf2y)/c02y*v2yTXBYTZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22TXBYTZ(i,j,k)-s22TXBYTZ(i,j-1,k))    
   v2zTXBYTZ(i,j,k)=(1-daf3y)/c03y*v2zTXBYTZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23TXBYTZ(i,j,k)-s23TXBYTZ(i,j,k-1))   
   v2TXBYTZ(i,j,k)=v2xTXBYTZ(i,j,k)+v2yTXBYTZ(i,j,k)+v2zTXBYTZ(i,j,k)
!
   v3xTXBYTZ(i,j,k)=(1-daf1z)/c01z*v3xTXBYTZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TXBYTZ(i,j,k)-s13TXBYTZ(i-1,j,k))
   v3yTXBYTZ(i,j,k)=(1-daf2z)/c02z*v3yTXBYTZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23TXBYTZ(i,j+1,k)-s23TXBYTZ(i,j,k))  
   v3zTXBYTZ(i,j,k)=(1-daf3z)/c03z*v3zTXBYTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TXBYTZ(i,j,k+1)-s33TXBYTZ(i,j,k))   
   v3TXBYTZ(i,j,k)=v3xTXBYTZ(i,j,k)+v3yTXBYTZ(i,j,k)+v3zTXBYTZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ##########################################################
Subroutine vPMLTXBYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
!
   daf1x=DiTXS1(i,3,3)
   c01x=1+daf1x
   daf1y=DiTXS0(i,3,3)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   daf2x=DiBYS0(nx-2,j,3)
   c02x=1+daf2x
   daf2y=DiBYS1(nx-2,j,3)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
! 
   daf3x=DiBZS0(nx-2,3,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(nx-2,3,k)
   c03z=1+daf3z
!
   bo1=(1/rhoTXBYBZ(i,j,k)+1/rhoTXBYBZ(i+1,j,k))/2
   bo2=(1/rhoTXBYBZ(i,j,k)+1/rhoTXBYBZ(i,j-1,k))/2
   bo3=(1/rhoTXBYBZ(i,j,k)+1/rhoTXBYBZ(i,j,k+1))/2
!
   v1xTXBYBZ(i,j,k)=(1-daf1x)/c01x*v1xTXBYBZ(i,j,k)+bo1*(dtx/c01x)*   &
    (s11TXBYBZ(i+1,j,k)-s11TXBYBZ(i,j,k))
   v1yTXBYBZ(i,j,k)=(1-daf2x)/c02x*v1yTXBYBZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12TXBYBZ(i,j+1,k)-s12TXBYBZ(i,j,k))   
   v1zTXBYBZ(i,j,k)=(1-daf3x)/c03x*v1zTXBYBZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13TXBYBZ(i,j,k)-s13TXBYBZ(i,j,k-1))  
   v1TXBYBZ(i,j,k)=v1xTXBYBZ(i,j,k)+v1yTXBYBZ(i,j,k)+v1zTXBYBZ(i,j,k)
!
   v2xTXBYBZ(i,j,k)=(1-daf1y)/c01y*v2xTXBYBZ(i,j,k)+bo2*(dtx/c01y)*   &
    (s12TXBYBZ(i,j,k)-s12TXBYBZ(i-1,j,k))    
   v2yTXBYBZ(i,j,k)=(1-daf2y)/c02y*v2yTXBYBZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22TXBYBZ(i,j,k)-s22TXBYBZ(i,j-1,k))    
   v2zTXBYBZ(i,j,k)=(1-daf3y)/c03y*v2zTXBYBZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23TXBYBZ(i,j,k)-s23TXBYBZ(i,j,k-1))   
   v2TXBYBZ(i,j,k)=v2xTXBYBZ(i,j,k)+v2yTXBYBZ(i,j,k)+v2zTXBYBZ(i,j,k)
!
   v3xTXBYBZ(i,j,k)=(1-daf1z)/c01z*v3xTXBYBZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13TXBYBZ(i,j,k)-s13TXBYBZ(i-1,j,k))
   v3yTXBYBZ(i,j,k)=(1-daf2z)/c02z*v3yTXBYBZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23TXBYBZ(i,j+1,k)-s23TXBYBZ(i,j,k))  
   v3zTXBYBZ(i,j,k)=(1-daf3z)/c03z*v3zTXBYBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33TXBYBZ(i,j,k+1)-s33TXBYBZ(i,j,k))   
   v3TXBYBZ(i,j,k)=v3xTXBYBZ(i,j,k)+v3yTXBYBZ(i,j,k)+v3zTXBYBZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ##########################################################
Subroutine vPMLBXTYTZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
!
   daf1x=DiBXS1(i,ny-2,nz-2)
   c01x=1+daf1x
   daf1y=DiBXS0(i,ny-2,nz-2)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   daf2x=DiTYS0(3,j,nz-2)
   c02x=1+daf2x
   daf2y=DiTYS1(3,j,nz-2)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
! 
   daf3x=DiTZS0(3,ny-2,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(3,ny-2,k)
   c03z=1+daf3z
!
   bo1=(1/rhoBXTYTZ(i,j,k)+1/rhoBXTYTZ(i+1,j,k))/2
   bo2=(1/rhoBXTYTZ(i,j,k)+1/rhoBXTYTZ(i,j-1,k))/2
   bo3=(1/rhoBXTYTZ(i,j,k)+1/rhoBXTYTZ(i,j,k+1))/2
!
   v1xBXTYTZ(i,j,k)=(1-daf1x)/c01x*v1xBXTYTZ(i,j,k)+bo1*(dtx/c01x)*   &
    (s11BXTYTZ(i+1,j,k)-s11BXTYTZ(i,j,k))
   v1yBXTYTZ(i,j,k)=(1-daf2x)/c02x*v1yBXTYTZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12BXTYTZ(i,j+1,k)-s12BXTYTZ(i,j,k))   
   v1zBXTYTZ(i,j,k)=(1-daf3x)/c03x*v1zBXTYTZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13BXTYTZ(i,j,k)-s13BXTYTZ(i,j,k-1))  
   v1BXTYTZ(i,j,k)=v1xBXTYTZ(i,j,k)+v1yBXTYTZ(i,j,k)+v1zBXTYTZ(i,j,k)
!
   v2xBXTYTZ(i,j,k)=(1-daf1y)/c01y*v2xBXTYTZ(i,j,k)+bo2*(dtx/c01y)*   &
    (s12BXTYTZ(i,j,k)-s12BXTYTZ(i-1,j,k))    
   v2yBXTYTZ(i,j,k)=(1-daf2y)/c02y*v2yBXTYTZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22BXTYTZ(i,j,k)-s22BXTYTZ(i,j-1,k))    
   v2zBXTYTZ(i,j,k)=(1-daf3y)/c03y*v2zBXTYTZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23BXTYTZ(i,j,k)-s23BXTYTZ(i,j,k-1))   
   v2BXTYTZ(i,j,k)=v2xBXTYTZ(i,j,k)+v2yBXTYTZ(i,j,k)+v2zBXTYTZ(i,j,k)
!
   v3xBXTYTZ(i,j,k)=(1-daf1z)/c01z*v3xBXTYTZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BXTYTZ(i,j,k)-s13BXTYTZ(i-1,j,k))
   v3yBXTYTZ(i,j,k)=(1-daf2z)/c02z*v3yBXTYTZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23BXTYTZ(i,j+1,k)-s23BXTYTZ(i,j,k))  
   v3zBXTYTZ(i,j,k)=(1-daf3z)/c03z*v3zBXTYTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BXTYTZ(i,j,k+1)-s33BXTYTZ(i,j,k))   
   v3BXTYTZ(i,j,k)=v3xBXTYTZ(i,j,k)+v3yBXTYTZ(i,j,k)+v3zBXTYTZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine

! ##########################################################
Subroutine vPMLBXTYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
!
   daf1x=DiBXS1(i,ny-2,3)
   c01x=1+daf1x
   daf1y=DiBXS0(i,ny-2,3)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   daf2x=DiTYS0(3,j,3)
   c02x=1+daf2x
   daf2y=DiTYS1(3,j,3)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
! 
   daf3x=DiBZS0(3,ny-2,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(3,ny-2,k)
   c03z=1+daf3z
!
   bo1=(1/rhoBXTYBZ(i,j,k)+1/rhoBXTYBZ(i+1,j,k))/2
   bo2=(1/rhoBXTYBZ(i,j,k)+1/rhoBXTYBZ(i,j-1,k))/2
   bo3=(1/rhoBXTYBZ(i,j,k)+1/rhoBXTYBZ(i,j,k+1))/2
!
   v1xBXTYBZ(i,j,k)=(1-daf1x)/c01x*v1xBXTYBZ(i,j,k)+bo1*(dtx/c01x)*   &
    (s11BXTYBZ(i+1,j,k)-s11BXTYBZ(i,j,k))
   v1yBXTYBZ(i,j,k)=(1-daf2x)/c02x*v1yBXTYBZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12BXTYBZ(i,j+1,k)-s12BXTYBZ(i,j,k))   
   v1zBXTYBZ(i,j,k)=(1-daf3x)/c03x*v1zBXTYBZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13BXTYBZ(i,j,k)-s13BXTYBZ(i,j,k-1))  
   v1BXTYBZ(i,j,k)=v1xBXTYBZ(i,j,k)+v1yBXTYBZ(i,j,k)+v1zBXTYBZ(i,j,k)
!
   v2xBXTYBZ(i,j,k)=(1-daf1y)/c01y*v2xBXTYBZ(i,j,k)+bo2*(dtx/c01y)*   &
    (s12BXTYBZ(i,j,k)-s12BXTYBZ(i-1,j,k))    
   v2yBXTYBZ(i,j,k)=(1-daf2y)/c02y*v2yBXTYBZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22BXTYBZ(i,j,k)-s22BXTYBZ(i,j-1,k))    
   v2zBXTYBZ(i,j,k)=(1-daf3y)/c03y*v2zBXTYBZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23BXTYBZ(i,j,k)-s23BXTYBZ(i,j,k-1))   
   v2BXTYBZ(i,j,k)=v2xBXTYBZ(i,j,k)+v2yBXTYBZ(i,j,k)+v2zBXTYBZ(i,j,k)
!
   v3xBXTYBZ(i,j,k)=(1-daf1z)/c01z*v3xBXTYBZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BXTYBZ(i,j,k)-s13BXTYBZ(i-1,j,k))
   v3yBXTYBZ(i,j,k)=(1-daf2z)/c02z*v3yBXTYBZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23BXTYBZ(i,j+1,k)-s23BXTYBZ(i,j,k))  
   v3zBXTYBZ(i,j,k)=(1-daf3z)/c03z*v3zBXTYBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BXTYBZ(i,j,k+1)-s33BXTYBZ(i,j,k))   
   v3BXTYBZ(i,j,k)=v3xBXTYBZ(i,j,k)+v3yBXTYBZ(i,j,k)+v3zBXTYBZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ##########################################################
Subroutine vPMLBXBYTZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
!
   daf1x=DiBXS1(i,3,nz-2)
   c01x=1+daf1x
   daf1y=DiBXS0(i,3,nz-2)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   daf2x=DiBYS0(3,j,nz-2)
   c02x=1+daf2x
   daf2y=DiBYS1(3,j,nz-2)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
! 
   daf3x=DiTZS0(3,3,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiTZS1(3,3,k)
   c03z=1+daf3z
!
   bo1=(1/rhoBXBYTZ(i,j,k)+1/rhoBXBYTZ(i+1,j,k))/2
   bo2=(1/rhoBXBYTZ(i,j,k)+1/rhoBXBYTZ(i,j-1,k))/2
   bo3=(1/rhoBXBYTZ(i,j,k)+1/rhoBXBYTZ(i,j,k+1))/2
!
   v1xBXBYTZ(i,j,k)=(1-daf1x)/c01x*v1xBXBYTZ(i,j,k)+bo1*(dtx/c01x)*   &
    (s11BXBYTZ(i+1,j,k)-s11BXBYTZ(i,j,k))
   v1yBXBYTZ(i,j,k)=(1-daf2x)/c02x*v1yBXBYTZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12BXBYTZ(i,j+1,k)-s12BXBYTZ(i,j,k))   
   v1zBXBYTZ(i,j,k)=(1-daf3x)/c03x*v1zBXBYTZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13BXBYTZ(i,j,k)-s13BXBYTZ(i,j,k-1))  
   v1BXBYTZ(i,j,k)=v1xBXBYTZ(i,j,k)+v1yBXBYTZ(i,j,k)+v1zBXBYTZ(i,j,k)
!
   v2xBXBYTZ(i,j,k)=(1-daf1y)/c01y*v2xBXBYTZ(i,j,k)+bo2*(dtx/c01y)*   &
    (s12BXBYTZ(i,j,k)-s12BXBYTZ(i-1,j,k))    
   v2yBXBYTZ(i,j,k)=(1-daf2y)/c02y*v2yBXBYTZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22BXBYTZ(i,j,k)-s22BXBYTZ(i,j-1,k))    
   v2zBXBYTZ(i,j,k)=(1-daf3y)/c03y*v2zBXBYTZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23BXBYTZ(i,j,k)-s23BXBYTZ(i,j,k-1))   
   v2BXBYTZ(i,j,k)=v2xBXBYTZ(i,j,k)+v2yBXBYTZ(i,j,k)+v2zBXBYTZ(i,j,k)
!
   v3xBXBYTZ(i,j,k)=(1-daf1z)/c01z*v3xBXBYTZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BXBYTZ(i,j,k)-s13BXBYTZ(i-1,j,k))
   v3yBXBYTZ(i,j,k)=(1-daf2z)/c02z*v3yBXBYTZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23BXBYTZ(i,j+1,k)-s23BXBYTZ(i,j,k))  
   v3zBXBYTZ(i,j,k)=(1-daf3z)/c03z*v3zBXBYTZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BXBYTZ(i,j,k+1)-s33BXBYTZ(i,j,k))   
   v3BXBYTZ(i,j,k)=v3xBXBYTZ(i,j,k)+v3yBXBYTZ(i,j,k)+v3zBXBYTZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ##########################################################
Subroutine vPMLBXBYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
!
   daf1x=DiBXS1(i,3,3)
   c01x=1+daf1x
   daf1y=DiBXS0(i,3,3)
   c01y=1+daf1y
   daf1z=daf1y
   c01z=c01y
!
   daf2x=DiBYS0(3,j,3)
   c02x=1+daf2x
   daf2y=DiBYS1(3,j,3)
   c02y=1+daf2y   
   daf2z=daf2x
   c02z=c02x
! 
   daf3x=DiBZS0(3,3,k)
   c03x=1+daf3x
   daf3y=daf3x
   c03y=c03x
   daf3z=DiBZS1(3,3,k)
   c03z=1+daf3z
!
   bo1=(1/rhoBXBYBZ(i,j,k)+1/rhoBXBYBZ(i+1,j,k))/2
   bo2=(1/rhoBXBYBZ(i,j,k)+1/rhoBXBYBZ(i,j-1,k))/2
   bo3=(1/rhoBXBYBZ(i,j,k)+1/rhoBXBYBZ(i,j,k+1))/2
!
   v1xBXBYBZ(i,j,k)=(1-daf1x)/c01x*v1xBXBYBZ(i,j,k)+bo1*(dtx/c01x)*   &
    (s11BXBYBZ(i+1,j,k)-s11BXBYBZ(i,j,k))
   v1yBXBYBZ(i,j,k)=(1-daf2x)/c02x*v1yBXBYBZ(i,j,k)+bo1*(dtx/c02x)*   &  
    (s12BXBYBZ(i,j+1,k)-s12BXBYBZ(i,j,k))   
   v1zBXBYBZ(i,j,k)=(1-daf3x)/c03x*v1zBXBYBZ(i,j,k)+bo1*(dtx/c03x)*    &
    (s13BXBYBZ(i,j,k)-s13BXBYBZ(i,j,k-1))  
   v1BXBYBZ(i,j,k)=v1xBXBYBZ(i,j,k)+v1yBXBYBZ(i,j,k)+v1zBXBYBZ(i,j,k)
!
   v2xBXBYBZ(i,j,k)=(1-daf1y)/c01y*v2xBXBYBZ(i,j,k)+bo2*(dtx/c01y)*   &
    (s12BXBYBZ(i,j,k)-s12BXBYBZ(i-1,j,k))    
   v2yBXBYBZ(i,j,k)=(1-daf2y)/c02y*v2yBXBYBZ(i,j,k)+bo2*(dtx/c02y)*   &  
    (s22BXBYBZ(i,j,k)-s22BXBYBZ(i,j-1,k))    
   v2zBXBYBZ(i,j,k)=(1-daf3y)/c03y*v2zBXBYBZ(i,j,k)+bo2*(dtx/c03y)*  &
    (s23BXBYBZ(i,j,k)-s23BXBYBZ(i,j,k-1))   
   v2BXBYBZ(i,j,k)=v2xBXBYBZ(i,j,k)+v2yBXBYBZ(i,j,k)+v2zBXBYBZ(i,j,k)
!
   v3xBXBYBZ(i,j,k)=(1-daf1z)/c01z*v3xBXBYBZ(i,j,k)+bo3*(dtx/c01z)*    &
    (s13BXBYBZ(i,j,k)-s13BXBYBZ(i-1,j,k))
   v3yBXBYBZ(i,j,k)=(1-daf2z)/c02z*v3yBXBYBZ(i,j,k)+bo3*(dtx/c02z)*   &  
    (s23BXBYBZ(i,j+1,k)-s23BXBYBZ(i,j,k))  
   v3zBXBYBZ(i,j,k)=(1-daf3z)/c03z*v3zBXBYBZ(i,j,k)+bo3*(dtx/c03z)*  &
    (s33BXBYBZ(i,j,k+1)-s33BXBYBZ(i,j,k))   
   v3BXBYBZ(i,j,k)=v3xBXBYBZ(i,j,k)+v3yBXBYBZ(i,j,k)+v3zBXBYBZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine PasteVeloC_TX
include 'triffy.dec'
 
do k=3,nz-2
 do j=3,ny-2
  v1(nx-1,j,k)=v1TX(2,j,k)
  v1(nx,j,k)=v1TX(3,j,k)
  v1TX(1,j,k)=v1(nx-2,j,k) 
  v2(nx-1,j,k)=v2TX(2,j,k)
  v2(nx,j,k)=v2TX(3,j,k)
  v2TX(1,j,k)=v2(nx-2,j,k) 
  v3(nx-1,j,k)=v3TX(2,j,k)
  v3(nx,j,k)=v3TX(3,j,k)
  v3TX(1,j,k)=v3(nx-2,j,k) 
 enddo 
enddo 
return 
end subroutine   
! ######################################
Subroutine PasteVeloC_BX
include 'triffy.dec'
 
do k=3,nz-2
 do j=3,ny-2
  v1(1,j,k)=v1BX(npm-2,j,k)
  v1(2,j,k)=v1BX(npm-1,j,k)
  v1BX(npm,j,k)=v1(3,j,k) 
  v2(1,j,k)=v2BX(npm-2,j,k)
  v2(2,j,k)=v2BX(npm-1,j,k)
  v2BX(npm,j,k)=v2(3,j,k) 
  v3(1,j,k)=v3BX(npm-2,j,k)
  v3(2,j,k)=v3BX(npm-1,j,k)
  v3BX(npm,j,k)=v3(3,j,k) 
 enddo 
enddo 
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloC_TY
include 'triffy.dec'
!
do k=3,nz-2
 do i=3,nx-2
  v1(i,ny-1,k)=v1TY(i,2,k)
  v1(i,ny,k)=v1TY(i,3,k)
  v1TY(i,1,k)=v1(i,ny-2,k)
  v2(i,ny-1,k)=v2TY(i,2,k)
  v2(i,ny,k)=v2TY(i,3,k)
  v2TY(i,1,k)=v2(i,ny-2,k)
  v3(i,ny-1,k)=v3TY(i,2,k)
  v3(i,ny,k)=v3TY(i,3,k)
  v3TY(i,1,k)=v3(i,ny-2,k)
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloC_BY
include 'triffy.dec'
!
do k=3,nz-2
 do i=3,nx-2
  v1(i,1,k)=v1BY(i,npm-2,k)
  v1(i,2,k)=v1BY(i,npm-1,k)
  v1BY(i,npm,k)=v1(i,3,k)
  v2(i,1,k)=v2BY(i,npm-2,k)
  v2(i,2,k)=v2BY(i,npm-1,k)
  v2BY(i,npm,k)=v2(i,3,k)
  v3(i,1,k)=v3BY(i,npm-2,k)
  v3(i,2,k)=v3BY(i,npm-1,k)
  v3BY(i,npm,k)=v3(i,3,k)
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloC_TZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=3,nx-2
  v1(i,j,nz-1)=v1TZ(i,j,2)
  v1(i,j,nz)=v1TZ(i,j,3)
  v1TZ(i,j,1)=v1(i,j,nz-2)
  v2(i,j,nz-1)=v2TZ(i,j,2) 
  v2(i,j,nz)=v2TZ(i,j,3) 
  v2TZ(i,j,1)=v2(i,j,nz-2) 
  v3(i,j,nz-1)=v3TZ(i,j,2) 
  v3(i,j,nz)=v3TZ(i,j,3) 
  v3TZ(i,j,1)=v3(i,j,nz-2) 
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloC_BZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=3,nx-2
  v1(i,j,1)=v1BZ(i,j,npm-2)
  v1(i,j,2)=v1BZ(i,j,npm-1)
  v1BZ(i,j,npm)=v1(i,j,3)
  v2(i,j,1)=v2BZ(i,j,npm-2) 
  v2(i,j,2)=v2BZ(i,j,npm-1) 
  v2BZ(i,j,npm)=v2(i,j,3) 
  v3(i,j,1)=v3BZ(i,j,npm-2) 
  v3(i,j,2)=v3BZ(i,j,npm-1) 
  v3BZ(i,j,npm)=v3(i,j,3) 
 enddo
enddo
return
end subroutine
! ########################################################### 
Subroutine PasteVeloTX_TXTY 
include 'triffy.dec'
do k=3,nz-2
 do i=1,npm 
  v1TX(i,ny-1,k)=v1TXTY(i,2,k) 
  v1TX(i,ny,k)=v1TXTY(i,3,k) 
  v1TXTY(i,1,k)=v1TX(i,ny-2,k) 
  v2TX(i,ny-1,k)=v2TXTY(i,2,k)  
  v2TX(i,ny,k)=v2TXTY(i,3,k)  
  v2TXTY(i,1,k)=v2TX(i,ny-2,k)   
  v3TX(i,ny-1,k)=v3TXTY(i,2,k)  
  v3TX(i,ny,k)=v3TXTY(i,3,k)  
  v3TXTY(i,1,k)=v3TX(i,ny-2,k)   
 enddo
enddo 
return
end subroutine 
! ########################################################### 
Subroutine PasteVeloBX_BXTY 
include 'triffy.dec'
do k=3,nz-2
 do i=1,npm 
  v1BX(i,ny-1,k)=v1BXTY(i,2,k) 
  v1BX(i,ny,k)=v1BXTY(i,3,k) 
  v1BXTY(i,1,k)=v1BX(i,ny-2,k) 
  v2BX(i,ny-1,k)=v2BXTY(i,2,k)  
  v2BX(i,ny,k)=v2BXTY(i,3,k)  
  v2BXTY(i,1,k)=v2BX(i,ny-2,k)   
  v3BX(i,ny-1,k)=v3BXTY(i,2,k)  
  v3BX(i,ny,k)=v3BXTY(i,3,k)  
  v3BXTY(i,1,k)=v3BX(i,ny-2,k)   
 enddo
enddo 
return
end subroutine 
! ########################################################### 
Subroutine PasteVeloTZ_TYTZ 
include 'triffy.dec'
do k=1,npm
 do i=3,nx-2
  v1TZ(i,ny-1,k)=v1TYTZ(i,2,k) 
  v1TZ(i,ny,k)=v1TYTZ(i,3,k) 
  v1TYTZ(i,1,k)=v1TZ(i,ny-2,k) 
  v2TZ(i,ny-1,k)=v2TYTZ(i,2,k)  
  v2TZ(i,ny,k)=v2TYTZ(i,3,k)  
  v2TYTZ(i,1,k)=v2TZ(i,ny-2,k)   
  v3TZ(i,ny-1,k)=v3TYTZ(i,2,k)  
  v3TZ(i,ny,k)=v3TYTZ(i,3,k)  
  v3TYTZ(i,1,k)=v3TZ(i,ny-2,k)   
 enddo
enddo 
return
end subroutine 
! ########################################################### 
Subroutine PasteVeloBZ_TYBZ 
include 'triffy.dec'
do k=1,npm
 do i=3,nx-2
  v1BZ(i,ny-1,k)=v1TYBZ(i,2,k) 
  v1BZ(i,ny,k)=v1TYBZ(i,3,k) 
  v1TYBZ(i,1,k)=v1BZ(i,ny-2,k) 
  v2BZ(i,ny-1,k)=v2TYBZ(i,2,k)  
  v2BZ(i,ny,k)=v2TYBZ(i,3,k)  
  v2TYBZ(i,1,k)=v2BZ(i,ny-2,k)   
  v3BZ(i,ny-1,k)=v3TYBZ(i,2,k)  
  v3BZ(i,ny,k)=v3TYBZ(i,3,k)  
  v3TYBZ(i,1,k)=v3BZ(i,ny-2,k)   
 enddo
enddo 
return
end subroutine 
! ###########################################################
Subroutine PasteVeloTX_TXBY
include 'triffy.dec'
!
do k=3,nz-2
 do i=1,npm
  v1TX(i,1,k)=v1TXBY(i,npm-2,k)
  v1TX(i,2,k)=v1TXBY(i,npm-1,k)
  v1TXBY(i,npm,k)=v1TX(i,3,k)
  v2TX(i,1,k)=v2TXBY(i,npm-2,k)
  v2TX(i,2,k)=v2TXBY(i,npm-1,k)
  v2TXBY(i,npm,k)=v2TX(i,3,k)
  v3TX(i,1,k)=v3TXBY(i,npm-2,k)
  v3TX(i,2,k)=v3TXBY(i,npm-1,k)
  v3TXBY(i,npm,k)=v3TX(i,3,k)
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloBX_BXBY
include 'triffy.dec'
!
do k=3,nz-2
 do i=1,npm
  v1BX(i,1,k)=v1BXBY(i,npm-2,k)
  v1BX(i,2,k)=v1BXBY(i,npm-1,k)
  v1BXBY(i,npm,k)=v1BX(i,3,k)
  v2BX(i,1,k)=v2BXBY(i,npm-2,k)
  v2BX(i,2,k)=v2BXBY(i,npm-1,k)
  v2BXBY(i,npm,k)=v2BX(i,3,k)
  v3BX(i,1,k)=v3BXBY(i,npm-2,k)
  v3BX(i,2,k)=v3BXBY(i,npm-1,k)
  v3BXBY(i,npm,k)=v3BX(i,3,k)
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloTZ_BYTZ
include 'triffy.dec'
!
do k=1,npm
 do i=3,nx-2
  v1TZ(i,1,k)=v1BYTZ(i,npm-2,k)
  v1TZ(i,2,k)=v1BYTZ(i,npm-1,k)
  v1BYTZ(i,npm,k)=v1TZ(i,3,k)
  v2TZ(i,1,k)=v2BYTZ(i,npm-2,k)
  v2TZ(i,2,k)=v2BYTZ(i,npm-1,k)
  v2BYTZ(i,npm,k)=v2TZ(i,3,k)
  v3TZ(i,1,k)=v3BYTZ(i,npm-2,k)
  v3TZ(i,2,k)=v3BYTZ(i,npm-1,k)
  v3BYTZ(i,npm,k)=v3TZ(i,3,k)
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloBZ_BYBZ
include 'triffy.dec'
!
do k=1,npm
 do i=3,nx-2
  v1BZ(i,1,k)=v1BYBZ(i,npm-2,k)
  v1BZ(i,2,k)=v1BYBZ(i,npm-1,k)
  v1BYBZ(i,npm,k)=v1BZ(i,3,k)
  v2BZ(i,1,k)=v2BYBZ(i,npm-2,k)
  v2BZ(i,2,k)=v2BYBZ(i,npm-1,k)
  v2BYBZ(i,npm,k)=v2BZ(i,3,k)
  v3BZ(i,1,k)=v3BYBZ(i,npm-2,k)
  v3BZ(i,2,k)=v3BYBZ(i,npm-1,k)
  v3BYBZ(i,npm,k)=v3BZ(i,3,k)
 enddo
enddo
return
end subroutine
! #############################################
Subroutine PasteVeloTX_TXTZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=1,npm
  v1TX(i,j,nz-1)=v1TXTZ(i,j,2)
  v1TX(i,j,nz)=v1TXTZ(i,j,3) 
  v1TXTZ(i,j,1)=v1TX(i,j,nz-2) 
  v2TX(i,j,nz-1)=v2TXTZ(i,j,2)
  v2TX(i,j,nz)=v2TXTZ(i,j,3) 
  v2TXTZ(i,j,1)=v2TX(i,j,nz-2)   
  v3TX(i,j,nz-1)=v3TXTZ(i,j,2)
  v3TX(i,j,nz)=v3TXTZ(i,j,3) 
  v3TXTZ(i,j,1)=v3TX(i,j,nz-2)   
 enddo 
enddo  
return 
end subroutine   
! #############################################
Subroutine PasteVeloBX_BXTZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=1,npm
  v1BX(i,j,nz-1)=v1BXTZ(i,j,2)
  v1BX(i,j,nz)=v1BXTZ(i,j,3) 
  v1BXTZ(i,j,1)=v1BX(i,j,nz-2) 
  v2BX(i,j,nz-1)=v2BXTZ(i,j,2)
  v2BX(i,j,nz)=v2BXTZ(i,j,3) 
  v2BXTZ(i,j,1)=v2BX(i,j,nz-2)   
  v3BX(i,j,nz-1)=v3BXTZ(i,j,2)
  v3BX(i,j,nz)=v3BXTZ(i,j,3) 
  v3BXTZ(i,j,1)=v3BX(i,j,nz-2)   
 enddo 
enddo  
return 
end subroutine   
! #############################################
Subroutine PasteVeloTY_TYTZ
include 'triffy.dec'
 
do j=1,npm
 do i=3,nx-2
  v1TY(i,j,nz-1)=v1TYTZ(i,j,2)
  v1TY(i,j,nz)=v1TYTZ(i,j,3) 
  v1TYTZ(i,j,1)=v1TY(i,j,nz-2) 
  v2TY(i,j,nz-1)=v2TYTZ(i,j,2)
  v2TY(i,j,nz)=v2TYTZ(i,j,3) 
  v2TYTZ(i,j,1)=v2TY(i,j,nz-2)   
  v3TY(i,j,nz-1)=v3TYTZ(i,j,2)
  v3TY(i,j,nz)=v3TYTZ(i,j,3) 
  v3TYTZ(i,j,1)=v3TY(i,j,nz-2)   
 enddo 
enddo  
return 
end subroutine   
! #############################################
Subroutine PasteVeloBY_BYTZ
include 'triffy.dec'
 
do j=1,npm
 do i=3,nx-2
  v1BY(i,j,nz-1)=v1BYTZ(i,j,2)
  v1BY(i,j,nz)=v1BYTZ(i,j,3) 
  v1BYTZ(i,j,1)=v1BY(i,j,nz-2) 
  v2BY(i,j,nz-1)=v2BYTZ(i,j,2)
  v2BY(i,j,nz)=v2BYTZ(i,j,3) 
  v2BYTZ(i,j,1)=v2BY(i,j,nz-2)   
  v3BY(i,j,nz-1)=v3BYTZ(i,j,2)
  v3BY(i,j,nz)=v3BYTZ(i,j,3) 
  v3BYTZ(i,j,1)=v3BY(i,j,nz-2)   
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloTX_TXBZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=1,npm
  v1TX(i,j,1)=v1TXBZ(i,j,npm-2)
  v1TX(i,j,2)=v1TXBZ(i,j,npm-1)
  v1TXBZ(i,j,npm)=v1TX(i,j,3)
  v2TX(i,j,1)=v2TXBZ(i,j,npm-2) 
  v2TX(i,j,2)=v2TXBZ(i,j,npm-1) 
  v2TXBZ(i,j,npm)=v2TX(i,j,3) 
  v3TX(i,j,1)=v3TXBZ(i,j,npm-2) 
  v3TX(i,j,2)=v3TXBZ(i,j,npm-1) 
  v3TXBZ(i,j,npm)=v3TX(i,j,3) 
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloBX_BXBZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=1,npm
  v1BX(i,j,1)=v1BXBZ(i,j,npm-2)
  v1BX(i,j,2)=v1BXBZ(i,j,npm-1)
  v1BXBZ(i,j,npm)=v1BX(i,j,3)
  v2BX(i,j,1)=v2BXBZ(i,j,npm-2) 
  v2BX(i,j,2)=v2BXBZ(i,j,npm-1) 
  v2BXBZ(i,j,npm)=v2BX(i,j,3) 
  v3BX(i,j,1)=v3BXBZ(i,j,npm-2) 
  v3BX(i,j,2)=v3BXBZ(i,j,npm-1) 
  v3BXBZ(i,j,npm)=v3BX(i,j,3) 
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloTY_TYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=3,nx-2
  v1TY(i,j,1)=v1TYBZ(i,j,npm-2)
  v1TY(i,j,2)=v1TYBZ(i,j,npm-1)
  v1TYBZ(i,j,npm)=v1TY(i,j,3)
  v2TY(i,j,1)=v2TYBZ(i,j,npm-2) 
  v2TY(i,j,2)=v2TYBZ(i,j,npm-1) 
  v2TYBZ(i,j,npm)=v2TY(i,j,3) 
  v3TY(i,j,1)=v3TYBZ(i,j,npm-2) 
  v3TY(i,j,2)=v3TYBZ(i,j,npm-1) 
  v3TYBZ(i,j,npm)=v3TY(i,j,3) 
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloBY_BYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=3,nx-2
  v1BY(i,j,1)=v1BYBZ(i,j,npm-2)
  v1BY(i,j,2)=v1BYBZ(i,j,npm-1)
  v1BYBZ(i,j,npm)=v1BY(i,j,3)
  v2BY(i,j,1)=v2BYBZ(i,j,npm-2) 
  v2BY(i,j,2)=v2BYBZ(i,j,npm-1) 
  v2BYBZ(i,j,npm)=v2BY(i,j,3) 
  v3BY(i,j,1)=v3BYBZ(i,j,npm-2) 
  v3BY(i,j,2)=v3BYBZ(i,j,npm-1) 
  v3BYBZ(i,j,npm)=v3BY(i,j,3) 
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloTY_TXTY
include 'triffy.dec'
do k=3,nz-2
 do j=1,npm
  v1TY(nx-1,j,k)=v1TXTY(2,j,k)
  v1TY(nx,j,k)=v1TXTY(3,j,k)
  v1TXTY(1,j,k)=v1TY(nx-2,j,k)
  v2TY(nx-1,j,k)=v2TXTY(2,j,k) 
  v2TY(nx,j,k)=v2TXTY(3,j,k) 
  v2TXTY(1,j,k)=v2TY(nx-2,j,k)   
  v3TY(nx-1,j,k)=v3TXTY(2,j,k) 
  v3TY(nx,j,k)=v3TXTY(3,j,k) 
  v3TXTY(1,j,k)=v3TY(nx-2,j,k)   
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloBY_TXBY
include 'triffy.dec'
do k=3,nz-2
 do j=1,npm
  v1BY(nx-1,j,k)=v1TXBY(2,j,k)
  v1BY(nx,j,k)=v1TXBY(3,j,k)
  v1TXBY(1,j,k)=v1BY(nx-2,j,k)
  v2BY(nx-1,j,k)=v2TXBY(2,j,k) 
  v2BY(nx,j,k)=v2TXBY(3,j,k) 
  v2TXBY(1,j,k)=v2BY(nx-2,j,k)   
  v3BY(nx-1,j,k)=v3TXBY(2,j,k) 
  v3BY(nx,j,k)=v3TXBY(3,j,k) 
  v3TXBY(1,j,k)=v3BY(nx-2,j,k)   
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloTZ_TXTZ
include 'triffy.dec'
do k=1,npm
 do j=3,ny-2
  v1TZ(nx-1,j,k)=v1TXTZ(2,j,k)
  v1TZ(nx,j,k)=v1TXTZ(3,j,k)
  v1TXTZ(1,j,k)=v1TZ(nx-2,j,k)
  v2TZ(nx-1,j,k)=v2TXTZ(2,j,k) 
  v2TZ(nx,j,k)=v2TXTZ(3,j,k) 
  v2TXTZ(1,j,k)=v2TZ(nx-2,j,k)   
  v3TZ(nx-1,j,k)=v3TXTZ(2,j,k) 
  v3TZ(nx,j,k)=v3TXTZ(3,j,k) 
  v3TXTZ(1,j,k)=v3TZ(nx-2,j,k)   
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloBZ_TXBZ
include 'triffy.dec'
do k=1,npm
 do j=3,ny-2
  v1BZ(nx-1,j,k)=v1TXBZ(2,j,k)
  v1BZ(nx,j,k)=v1TXBZ(3,j,k)
  v1TXBZ(1,j,k)=v1BZ(nx-2,j,k)
  v2BZ(nx-1,j,k)=v2TXBZ(2,j,k) 
  v2BZ(nx,j,k)=v2TXBZ(3,j,k) 
  v2TXBZ(1,j,k)=v2BZ(nx-2,j,k)   
  v3BZ(nx-1,j,k)=v3TXBZ(2,j,k) 
  v3BZ(nx,j,k)=v3TXBZ(3,j,k) 
  v3TXBZ(1,j,k)=v3BZ(nx-2,j,k)   
 enddo 
enddo  
return 
end subroutine   
! ######################################
Subroutine PasteVeloTY_BXTY
include 'triffy.dec'
 
do k=3,nz-2
 do j=1,npm
  v1TY(1,j,k)=v1BXTY(npm-2,j,k)
  v1TY(2,j,k)=v1BXTY(npm-1,j,k)
  v1BXTY(npm,j,k)=v1TY(3,j,k) 
  v2TY(1,j,k)=v2BXTY(npm-2,j,k)
  v2TY(2,j,k)=v2BXTY(npm-1,j,k)
  v2BXTY(npm,j,k)=v2TY(3,j,k) 
  v3TY(1,j,k)=v3BXTY(npm-2,j,k)
  v3TY(2,j,k)=v3BXTY(npm-1,j,k)
  v3BXTY(npm,j,k)=v3TY(3,j,k) 
 enddo 
enddo 
return 
end subroutine   
! ######################################
Subroutine PasteVeloBY_BXBY
include 'triffy.dec'
 
do k=3,nz-2
 do j=1,npm
  v1BY(1,j,k)=v1BXBY(npm-2,j,k)
  v1BY(2,j,k)=v1BXBY(npm-1,j,k)
  v1BXBY(npm,j,k)=v1BY(3,j,k) 
  v2BY(1,j,k)=v2BXBY(npm-2,j,k)
  v2BY(2,j,k)=v2BXBY(npm-1,j,k)
  v2BXBY(npm,j,k)=v2BY(3,j,k) 
  v3BY(1,j,k)=v3BXBY(npm-2,j,k)
  v3BY(2,j,k)=v3BXBY(npm-1,j,k)
  v3BXBY(npm,j,k)=v3BY(3,j,k) 
 enddo 
enddo 
return 
end subroutine   
! ######################################
Subroutine PasteVeloTZ_BXTZ
include 'triffy.dec'
 
do k=1,npm
 do j=3,ny-2
  v1TZ(1,j,k)=v1BXTZ(npm-2,j,k)
  v1TZ(2,j,k)=v1BXTZ(npm-1,j,k)
  v1BXTZ(npm,j,k)=v1TZ(3,j,k) 
  v2TZ(1,j,k)=v2BXTZ(npm-2,j,k)
  v2TZ(2,j,k)=v2BXTZ(npm-1,j,k)
  v2BXTZ(npm,j,k)=v2TZ(3,j,k) 
  v3TZ(1,j,k)=v3BXTZ(npm-2,j,k)
  v3TZ(2,j,k)=v3BXTZ(npm-1,j,k)
  v3BXTZ(npm,j,k)=v3TZ(3,j,k) 
 enddo 
enddo 
return 
end subroutine   
! ######################################
Subroutine PasteVeloBZ_BXBZ
include 'triffy.dec'
 
do k=1,npm
 do j=3,ny-2
  v1BZ(1,j,k)=v1BXBZ(npm-2,j,k)
  v1BZ(2,j,k)=v1BXBZ(npm-1,j,k)
  v1BXBZ(npm,j,k)=v1BZ(3,j,k) 
  v2BZ(1,j,k)=v2BXBZ(npm-2,j,k)
  v2BZ(2,j,k)=v2BXBZ(npm-1,j,k)
  v2BXBZ(npm,j,k)=v2BZ(3,j,k) 
  v3BZ(1,j,k)=v3BXBZ(npm-2,j,k)
  v3BZ(2,j,k)=v3BXBZ(npm-1,j,k)
  v3BXBZ(npm,j,k)=v3BZ(3,j,k) 
 enddo 
enddo 
return 
end subroutine
! ########################################################### 
Subroutine PasteVeloTXTY_TXTYTZ
include 'triffy.dec'
do j=1,npm 
 do i=1,npm 
  v1TXTY(i,j,nz-1)=v1TXTYTZ(i,j,2) 
  v1TXTY(i,j,nz)=v1TXTYTZ(i,j,3)
  v1TXTYTZ(i,j,1)=v1TXTY(i,j,nz-2) 
  v2TXTY(i,j,nz-1)=v2TXTYTZ(i,j,2) 
  v2TXTY(i,j,nz)=v2TXTYTZ(i,j,3)
  v2TXTYTZ(i,j,1)=v2TXTY(i,j,nz-2) 
  v3TXTY(i,j,nz-1)=v3TXTYTZ(i,j,2) 
  v3TXTY(i,j,nz)=v3TXTYTZ(i,j,3)
  v3TXTYTZ(i,j,1)=v3TXTY(i,j,nz-2) 
 enddo
enddo
return 
end subroutine 
! ########################################################### 
Subroutine PasteVeloTXTZ_TXTYTZ 
include 'triffy.dec'
do k=1,npm
 do i=1,npm 
  v1TXTZ(i,ny-1,k)=v1TXTYTZ(i,2,k) 
  v1TXTZ(i,ny,k)=v1TXTYTZ(i,3,k) 
  v1TXTYTZ(i,1,k)=v1TXTZ(i,ny-2,k) 
  v2TXTZ(i,ny-1,k)=v2TXTYTZ(i,2,k)  
  v2TXTZ(i,ny,k)=v2TXTYTZ(i,3,k)  
  v2TXTYTZ(i,1,k)=v2TXTZ(i,ny-2,k)   
  v3TXTZ(i,ny-1,k)=v3TXTYTZ(i,2,k)  
  v3TXTZ(i,ny,k)=v3TXTYTZ(i,3,k)  
  v3TXTYTZ(i,1,k)=v3TXTZ(i,ny-2,k)   
 enddo
enddo 
return
end subroutine 
! ###########################################################
Subroutine PasteVeloTYTZ_TXTYTZ
include 'triffy.dec'
do k=1,npm
 do j=1,npm
  v1TYTZ(nx-1,j,k)=v1TXTYTZ(2,j,k)
  v1TYTZ(nx,j,k)=v1TXTYTZ(3,j,k)
  v1TXTYTZ(1,j,k)=v1TYTZ(nx-2,j,k)
  v2TYTZ(nx-1,j,k)=v2TXTYTZ(2,j,k) 
  v2TYTZ(nx,j,k)=v2TXTYTZ(3,j,k) 
  v2TXTYTZ(1,j,k)=v2TYTZ(nx-2,j,k)   
  v3TYTZ(nx-1,j,k)=v3TXTYTZ(2,j,k) 
  v3TYTZ(nx,j,k)=v3TXTYTZ(3,j,k) 
  v3TXTYTZ(1,j,k)=v3TYTZ(nx-2,j,k)   
 enddo 
enddo  
return 
end subroutine   
! ########################################################### 
Subroutine PasteVeloTXTY_TXTYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  v1TXTY(i,j,1)=v1TXTYBZ(i,j,npm-2)
  v1TXTY(i,j,2)=v1TXTYBZ(i,j,npm-1)
  v1TXTYBZ(i,j,npm)=v1TXTY(i,j,3)
  v2TXTY(i,j,1)=v2TXTYBZ(i,j,npm-2) 
  v2TXTY(i,j,2)=v2TXTYBZ(i,j,npm-1) 
  v2TXTYBZ(i,j,npm)=v2TXTY(i,j,3) 
  v3TXTY(i,j,1)=v3TXTYBZ(i,j,npm-2) 
  v3TXTY(i,j,2)=v3TXTYBZ(i,j,npm-1) 
  v3TXTYBZ(i,j,npm)=v3TXTY(i,j,3) 
 enddo
enddo
return
end subroutine
! ########################################################### 
Subroutine PasteVeloTXBZ_TXTYBZ 
include 'triffy.dec'
do k=1,npm
 do i=1,npm 
  v1TXBZ(i,ny-1,k)=v1TXTYBZ(i,2,k) 
  v1TXBZ(i,ny,k)=v1TXTYBZ(i,3,k) 
  v1TXTYBZ(i,1,k)=v1TXBZ(i,ny-2,k) 
  v2TXBZ(i,ny-1,k)=v2TXTYBZ(i,2,k)  
  v2TXBZ(i,ny,k)=v2TXTYBZ(i,3,k)  
  v2TXTYBZ(i,1,k)=v2TXBZ(i,ny-2,k)   
  v3TXBZ(i,ny-1,k)=v3TXTYBZ(i,2,k)  
  v3TXBZ(i,ny,k)=v3TXTYBZ(i,3,k)  
  v3TXTYBZ(i,1,k)=v3TXBZ(i,ny-2,k)   
 enddo
enddo 
return
end subroutine 
! ###########################################################
Subroutine PasteVeloTYBZ_TXTYBZ
include 'triffy.dec'
do k=1,npm
 do j=1,npm
  v1TYBZ(nx-1,j,k)=v1TXTYBZ(2,j,k)
  v1TYBZ(nx,j,k)=v1TXTYBZ(3,j,k)
  v1TXTYBZ(1,j,k)=v1TYBZ(nx-2,j,k)
  v2TYBZ(nx-1,j,k)=v2TXTYBZ(2,j,k) 
  v2TYBZ(nx,j,k)=v2TXTYBZ(3,j,k) 
  v2TXTYBZ(1,j,k)=v2TYBZ(nx-2,j,k)   
  v3TYBZ(nx-1,j,k)=v3TXTYBZ(2,j,k) 
  v3TYBZ(nx,j,k)=v3TXTYBZ(3,j,k) 
  v3TXTYBZ(1,j,k)=v3TYBZ(nx-2,j,k)   
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloTXBY_TXBYTZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  v1TXBY(i,j,nz-1)=v1TXBYTZ(i,j,2)
  v1TXBY(i,j,nz)=v1TXBYTZ(i,j,3) 
  v1TXBYTZ(i,j,1)=v1TXBY(i,j,nz-2) 
  v2TXBY(i,j,nz-1)=v2TXBYTZ(i,j,2)
  v2TXBY(i,j,nz)=v2TXBYTZ(i,j,3) 
  v2TXBYTZ(i,j,1)=v2TXBY(i,j,nz-2)   
  v3TXBY(i,j,nz-1)=v3TXBYTZ(i,j,2)
  v3TXBY(i,j,nz)=v3TXBYTZ(i,j,3) 
  v3TXBYTZ(i,j,1)=v3TXBY(i,j,nz-2)   
 enddo 
enddo  
return 
end subroutine
! ###########################################################
Subroutine PasteVeloTXTZ_TXBYTZ
include 'triffy.dec'
!
do k=1,npm
 do i=1,npm
  v1TXTZ(i,1,k)=v1TXBYTZ(i,npm-2,k)
  v1TXTZ(i,2,k)=v1TXBYTZ(i,npm-1,k)
  v1TXBYTZ(i,npm,k)=v1TXTZ(i,3,k)
  v2TXTZ(i,1,k)=v2TXBYTZ(i,npm-2,k)
  v2TXTZ(i,2,k)=v2TXBYTZ(i,npm-1,k)
  v2TXBYTZ(i,npm,k)=v2TXTZ(i,3,k)
  v3TXTZ(i,1,k)=v3TXBYTZ(i,npm-2,k)
  v3TXTZ(i,2,k)=v3TXBYTZ(i,npm-1,k)
  v3TXBYTZ(i,npm,k)=v3TXTZ(i,3,k)
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloBYTZ_TXBYTZ
include 'triffy.dec'
do k=1,npm
 do j=1,npm
  v1BYTZ(nx-1,j,k)=v1TXBYTZ(2,j,k)
  v1BYTZ(nx,j,k)=v1TXBYTZ(3,j,k)
  v1TXBYTZ(1,j,k)=v1BYTZ(nx-2,j,k)
  v2BYTZ(nx-1,j,k)=v2TXBYTZ(2,j,k) 
  v2BYTZ(nx,j,k)=v2TXBYTZ(3,j,k) 
  v2TXBYTZ(1,j,k)=v2BYTZ(nx-2,j,k)   
  v3BYTZ(nx-1,j,k)=v3TXBYTZ(2,j,k) 
  v3BYTZ(nx,j,k)=v3TXBYTZ(3,j,k) 
  v3TXBYTZ(1,j,k)=v3BYTZ(nx-2,j,k)   
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloTXBY_TXBYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  v1TXBY(i,j,1)=v1TXBYBZ(i,j,npm-2)
  v1TXBY(i,j,2)=v1TXBYBZ(i,j,npm-1)
  v1TXBYBZ(i,j,npm)=v1TXBY(i,j,3)
  v2TXBY(i,j,1)=v2TXBYBZ(i,j,npm-2) 
  v2TXBY(i,j,2)=v2TXBYBZ(i,j,npm-1) 
  v2TXBYBZ(i,j,npm)=v2TXBY(i,j,3) 
  v3TXBY(i,j,1)=v3TXBYBZ(i,j,npm-2) 
  v3TXBY(i,j,2)=v3TXBYBZ(i,j,npm-1) 
  v3TXBYBZ(i,j,npm)=v3TXBY(i,j,3) 
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloTXBZ_TXBYBZ
include 'triffy.dec'
!
do k=1,npm
 do i=1,npm
  v1TXBZ(i,1,k)=v1TXBYBZ(i,npm-2,k)
  v1TXBZ(i,2,k)=v1TXBYBZ(i,npm-1,k)
  v1TXBYBZ(i,npm,k)=v1TXBZ(i,3,k)
  v2TXBZ(i,1,k)=v2TXBYBZ(i,npm-2,k)
  v2TXBZ(i,2,k)=v2TXBYBZ(i,npm-1,k)
  v2TXBYBZ(i,npm,k)=v2TXBZ(i,3,k)
  v3TXBZ(i,1,k)=v3TXBYBZ(i,npm-2,k)
  v3TXBZ(i,2,k)=v3TXBYBZ(i,npm-1,k)
  v3TXBYBZ(i,npm,k)=v3TXBZ(i,3,k)
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloBYBZ_TXBYBZ
include 'triffy.dec'
do k=1,npm
 do j=1,npm
  v1BYBZ(nx-1,j,k)=v1TXBYBZ(2,j,k)
  v1BYBZ(nx,j,k)=v1TXBYBZ(3,j,k)
  v1TXBYBZ(1,j,k)=v1BYBZ(nx-2,j,k)
  v2BYBZ(nx-1,j,k)=v2TXBYBZ(2,j,k) 
  v2BYBZ(nx,j,k)=v2TXBYBZ(3,j,k) 
  v2TXBYBZ(1,j,k)=v2BYBZ(nx-2,j,k)   
  v3BYBZ(nx-1,j,k)=v3TXBYBZ(2,j,k) 
  v3BYBZ(nx,j,k)=v3TXBYBZ(3,j,k) 
  v3TXBYBZ(1,j,k)=v3BYBZ(nx-2,j,k)   
 enddo 
enddo  
return 
end subroutine   
! ########################################################### 
Subroutine PasteVeloBXTY_BXTYTZ
include 'triffy.dec'
do j=1,npm 
 do i=1,npm 
  v1BXTY(i,j,nz-1)=v1BXTYTZ(i,j,2) 
  v1BXTY(i,j,nz)=v1BXTYTZ(i,j,3)
  v1BXTYTZ(i,j,1)=v1BXTY(i,j,nz-2) 
  v2BXTY(i,j,nz-1)=v2BXTYTZ(i,j,2) 
  v2BXTY(i,j,nz)=v2BXTYTZ(i,j,3)
  v2BXTYTZ(i,j,1)=v2BXTY(i,j,nz-2) 
  v3BXTY(i,j,nz-1)=v3BXTYTZ(i,j,2) 
  v3BXTY(i,j,nz)=v3BXTYTZ(i,j,3)
  v3BXTYTZ(i,j,1)=v3BXTY(i,j,nz-2) 
 enddo
enddo
return 
end subroutine 
! ########################################################### 
Subroutine PasteVeloBXTZ_BXTYTZ 
include 'triffy.dec'
do k=1,npm
 do i=1,npm 
  v1BXTZ(i,ny-1,k)=v1BXTYTZ(i,2,k) 
  v1BXTZ(i,ny,k)=v1BXTYTZ(i,3,k) 
  v1BXTYTZ(i,1,k)=v1BXTZ(i,ny-2,k) 
  v2BXTZ(i,ny-1,k)=v2BXTYTZ(i,2,k)  
  v2BXTZ(i,ny,k)=v2BXTYTZ(i,3,k)  
  v2BXTYTZ(i,1,k)=v2BXTZ(i,ny-2,k)   
  v3BXTZ(i,ny-1,k)=v3BXTYTZ(i,2,k)  
  v3BXTZ(i,ny,k)=v3BXTYTZ(i,3,k)  
  v3BXTYTZ(i,1,k)=v3BXTZ(i,ny-2,k)   
 enddo
enddo 
return
end subroutine 
! ######################################
Subroutine PasteVeloTYTZ_BXTYTZ
include 'triffy.dec'
 
do k=1,npm
 do j=1,npm
  v1TYTZ(1,j,k)=v1BXTYTZ(npm-2,j,k)
  v1TYTZ(2,j,k)=v1BXTYTZ(npm-1,j,k)
  v1BXTYTZ(npm,j,k)=v1TYTZ(3,j,k) 
  v2TYTZ(1,j,k)=v2BXTYTZ(npm-2,j,k)
  v2TYTZ(2,j,k)=v2BXTYTZ(npm-1,j,k)
  v2BXTYTZ(npm,j,k)=v2TYTZ(3,j,k) 
  v3TYTZ(1,j,k)=v3BXTYTZ(npm-2,j,k)
  v3TYTZ(2,j,k)=v3BXTYTZ(npm-1,j,k)
  v3BXTYTZ(npm,j,k)=v3TYTZ(3,j,k) 
 enddo 
enddo 
return 
end subroutine   
! ########################################################### 
Subroutine PasteVeloBXTY_BXTYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  v1BXTY(i,j,1)=v1BXTYBZ(i,j,npm-2)
  v1BXTY(i,j,2)=v1BXTYBZ(i,j,npm-1)
  v1BXTYBZ(i,j,npm)=v1BXTY(i,j,3)
  v2BXTY(i,j,1)=v2BXTYBZ(i,j,npm-2) 
  v2BXTY(i,j,2)=v2BXTYBZ(i,j,npm-1) 
  v2BXTYBZ(i,j,npm)=v2BXTY(i,j,3) 
  v3BXTY(i,j,1)=v3BXTYBZ(i,j,npm-2) 
  v3BXTY(i,j,2)=v3BXTYBZ(i,j,npm-1) 
  v3BXTYBZ(i,j,npm)=v3BXTY(i,j,3) 
 enddo
enddo
return
end subroutine
! ########################################################### 
Subroutine PasteVeloBXBZ_BXTYBZ 
include 'triffy.dec'
do k=1,npm
 do i=1,npm 
  v1BXBZ(i,ny-1,k)=v1BXTYBZ(i,2,k) 
  v1BXBZ(i,ny,k)=v1BXTYBZ(i,3,k) 
  v1BXTYBZ(i,1,k)=v1BXBZ(i,ny-2,k) 
  v2BXBZ(i,ny-1,k)=v2BXTYBZ(i,2,k)  
  v2BXBZ(i,ny,k)=v2BXTYBZ(i,3,k)  
  v2BXTYBZ(i,1,k)=v2BXBZ(i,ny-2,k)   
  v3BXBZ(i,ny-1,k)=v3BXTYBZ(i,2,k)  
  v3BXBZ(i,ny,k)=v3BXTYBZ(i,3,k)  
  v3BXTYBZ(i,1,k)=v3BXBZ(i,ny-2,k)   
 enddo
enddo 
return
end subroutine 
! ######################################
Subroutine PasteVeloTYBZ_BXTYBZ
include 'triffy.dec'
 
do k=1,npm
 do j=1,npm
  v1TYBZ(1,j,k)=v1BXTYBZ(npm-2,j,k)
  v1TYBZ(2,j,k)=v1BXTYBZ(npm-1,j,k)
  v1BXTYBZ(npm,j,k)=v1TYBZ(3,j,k) 
  v2TYBZ(1,j,k)=v2BXTYBZ(npm-2,j,k)
  v2TYBZ(2,j,k)=v2BXTYBZ(npm-1,j,k)
  v2BXTYBZ(npm,j,k)=v2TYBZ(3,j,k) 
  v3TYBZ(1,j,k)=v3BXTYBZ(npm-2,j,k)
  v3TYBZ(2,j,k)=v3BXTYBZ(npm-1,j,k)
  v3BXTYBZ(npm,j,k)=v3TYBZ(3,j,k) 
 enddo 
enddo 
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloBXBY_BXBYTZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  v1BXBY(i,j,nz-1)=v1BXBYTZ(i,j,2)
  v1BXBY(i,j,nz)=v1BXBYTZ(i,j,3) 
  v1BXBYTZ(i,j,1)=v1BXBY(i,j,nz-2) 
  v2BXBY(i,j,nz-1)=v2BXBYTZ(i,j,2)
  v2BXBY(i,j,nz)=v2BXBYTZ(i,j,3) 
  v2BXBYTZ(i,j,1)=v2BXBY(i,j,nz-2)   
  v3BXBY(i,j,nz-1)=v3BXBYTZ(i,j,2)
  v3BXBY(i,j,nz)=v3BXBYTZ(i,j,3) 
  v3BXBYTZ(i,j,1)=v3BXBY(i,j,nz-2)   
 enddo 
enddo  
return 
end subroutine
! ###########################################################
Subroutine PasteVeloBXTZ_BXBYTZ
include 'triffy.dec'
!
do k=1,npm
 do i=1,npm
  v1BXTZ(i,1,k)=v1BXBYTZ(i,npm-2,k)
  v1BXTZ(i,2,k)=v1BXBYTZ(i,npm-1,k)
  v1BXBYTZ(i,npm,k)=v1BXTZ(i,3,k)
  v2BXTZ(i,1,k)=v2BXBYTZ(i,npm-2,k)
  v2BXTZ(i,2,k)=v2BXBYTZ(i,npm-1,k)
  v2BXBYTZ(i,npm,k)=v2BXTZ(i,3,k)
  v3BXTZ(i,1,k)=v3BXBYTZ(i,npm-2,k)
  v3BXTZ(i,2,k)=v3BXBYTZ(i,npm-1,k)
  v3BXBYTZ(i,npm,k)=v3BXTZ(i,3,k)
 enddo
enddo
return
end subroutine
! ######################################
Subroutine PasteVeloBYTZ_BXBYTZ
include 'triffy.dec'
 
do k=1,npm
 do j=1,npm
  v1BYTZ(1,j,k)=v1BXBYTZ(npm-2,j,k)
  v1BYTZ(2,j,k)=v1BXBYTZ(npm-1,j,k)
  v1BXBYTZ(npm,j,k)=v1BYTZ(3,j,k) 
  v2BYTZ(1,j,k)=v2BXBYTZ(npm-2,j,k)
  v2BYTZ(2,j,k)=v2BXBYTZ(npm-1,j,k)
  v2BXBYTZ(npm,j,k)=v2BYTZ(3,j,k) 
  v3BYTZ(1,j,k)=v3BXBYTZ(npm-2,j,k)
  v3BYTZ(2,j,k)=v3BXBYTZ(npm-1,j,k)
  v3BXBYTZ(npm,j,k)=v3BYTZ(3,j,k) 
 enddo 
enddo 
return 
end subroutine   
! ###########################################################
Subroutine PasteVeloBXBY_BXBYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  v1BXBY(i,j,1)=v1BXBYBZ(i,j,npm-2)
  v1BXBY(i,j,2)=v1BXBYBZ(i,j,npm-1)
  v1BXBYBZ(i,j,npm)=v1BXBY(i,j,3)
  v2BXBY(i,j,1)=v2BXBYBZ(i,j,npm-2) 
  v2BXBY(i,j,2)=v2BXBYBZ(i,j,npm-1) 
  v2BXBYBZ(i,j,npm)=v2BXBY(i,j,3) 
  v3BXBY(i,j,1)=v3BXBYBZ(i,j,npm-2) 
  v3BXBY(i,j,2)=v3BXBYBZ(i,j,npm-1) 
  v3BXBYBZ(i,j,npm)=v3BXBY(i,j,3) 
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteVeloBXBZ_BXBYBZ
include 'triffy.dec'
!
do k=1,npm
 do i=1,npm
  v1BXBZ(i,1,k)=v1BXBYBZ(i,npm-2,k)
  v1BXBZ(i,2,k)=v1BXBYBZ(i,npm-1,k)
  v1BXBYBZ(i,npm,k)=v1BXBZ(i,3,k)
  v2BXBZ(i,1,k)=v2BXBYBZ(i,npm-2,k)
  v2BXBZ(i,2,k)=v2BXBYBZ(i,npm-1,k)
  v2BXBYBZ(i,npm,k)=v2BXBZ(i,3,k)
  v3BXBZ(i,1,k)=v3BXBYBZ(i,npm-2,k)
  v3BXBZ(i,2,k)=v3BXBYBZ(i,npm-1,k)
  v3BXBYBZ(i,npm,k)=v3BXBZ(i,3,k)
 enddo
enddo
return
end subroutine
! ######################################
Subroutine PasteVeloBYBZ_BXBYBZ
include 'triffy.dec'
 
do k=1,npm
 do j=1,npm
  v1BYBZ(1,j,k)=v1BXBYBZ(npm-2,j,k)
  v1BYBZ(2,j,k)=v1BXBYBZ(npm-1,j,k)
  v1BXBYBZ(npm,j,k)=v1BYBZ(3,j,k) 
  v2BYBZ(1,j,k)=v2BXBYBZ(npm-2,j,k)
  v2BYBZ(2,j,k)=v2BXBYBZ(npm-1,j,k)
  v2BXBYBZ(npm,j,k)=v2BYBZ(3,j,k) 
  v3BYBZ(1,j,k)=v3BXBYBZ(npm-2,j,k)
  v3BYBZ(2,j,k)=v3BXBYBZ(npm-1,j,k)
  v3BXBYBZ(npm,j,k)=v3BYBZ(3,j,k) 
 enddo 
enddo 
return 
end subroutine   


! ###################################
Subroutine sPMLTZ
include 'triffy.dec'
do k=2,npm
 do j=3,ny-2
  do i=3,nx-2
   daf=DiTZS0(i,j,k)
   c0=1+daf
   s11oTZ(i,j,k)=(1-daf)/c0*s11oTZ(i,j,k)+dtx/c0*            &
    lamTZ(i,j,k)*(v3TZ(i,j,k)-v3TZ(i,j,k-1))
   s11pTZ(i,j,k)=s11pTZ(i,j,k)+dtx*((lamTZ(i,j,k)             &
    +2.*muTZ(i,j,k))*(rn*(v1TZ(i,j,k)-v1TZ(i-1,j,k))              &
    +rnn*(v1TZ(i+1,j,k)-v1TZ(i-2,j,k)))+lamTZ(i,j,k)*             &
    (rn*(v2TZ(i,j+1,k)-v2TZ(i,j,k))+rnn*(v2TZ(i,j+2,k)-           &
    v2TZ(i,j-1,k))))
   s11TZ(i,j,k)=s11pTZ(i,j,k)+s11oTZ(i,j,k)
! 
   s22oTZ(i,j,k)=(1-daf)/c0*s22oTZ(i,j,k)+dtx/c0*            &  
    lamTZ(i,j,k)*(v3TZ(i,j,k)-v3TZ(i,j,k-1))  
   s22pTZ(i,j,k)=s22pTZ(i,j,k)+dtx*(lamTZ(i,j,k)              & 
    *(rn*(v1TZ(i,j,k)-v1TZ(i-1,j,k))                             & 
    +rnn*(v1TZ(i+1,j,k)-v1TZ(i-2,j,k)))+(lamTZ(i,j,k)+            &     
    2*muTZ(i,j,k))*                                             & 
    (rn*(v2TZ(i,j+1,k)-v2TZ(i,j,k))+rnn*(v2TZ(i,j+2,k)-           & 
    v2TZ(i,j-1,k))) ) 
   s22TZ(i,j,k)=s22oTZ(i,j,k)+s22pTZ(i,j,k)  
! 
   s33oTZ(i,j,k)=(1-daf)/c0*s33oTZ(i,j,k)+dtx/c0*            &   
    (lamTZ(i,j,k)+2*muTZ(i,j,k))*(v3TZ(i,j,k)-v3TZ(i,j,k-1))   
   s33pTZ(i,j,k)=s33pTZ(i,j,k)+dtx*lamTZ(i,j,k)               &  
    *(rn*(v1TZ(i,j,k)-v1TZ(i-1,j,k)+v2TZ(i,j+1,k)-v2TZ(i,j,k))     &
    +rnn*(v1TZ(i+1,j,k)-v1TZ(i-2,j,k)+v2TZ(i,j+2,k)-              &
    v2TZ(i,j-1,k))) 
   s33TZ(i,j,k)=s33pTZ(i,j,k)+s33oTZ(i,j,k)             
!
   mu12=4./(1./muTZ(i,j,k)+1./muTZ(i,j-1,k)+1./muTZ(i+1,j-1,k)+  &
    1./muTZ(i+1,j,k))
   s12TZ(i,j,k)=s12TZ(i,j,k)+dtx*mu12*(rn*(v1TZ(i,j,k)-      &
    v1TZ(i,j-1,k)+v2TZ(i+1,j,k)-v2TZ(i,j,k))+rnn*(v1TZ(i,j+1,k)-  &
    v1TZ(i,j-2,k)+v2TZ(i+2,j,k)-v2TZ(i-1,j,k)))
  enddo 
 enddo 
enddo 
!
do k=2,npm-1
 do j=3,ny-2 
  do i=3,nx-2
   daf=DiTZS1(i,j,k)
   c0=1.+daf
   mu13=4./(1./muTZ(i,j,k)+1./muTZ(i,j,k+1)+1./muTZ(i+1,j,k)+    &
    1./muTZ(i+1,j,k+1))
   mu23=4./(1./muTZ(i,j,k)+1./muTZ(i,j,k+1)+1./muTZ(i,j-1,k)+    &
    1./muTZ(i,j-1,k+1))
!
   s13oTZ(i,j,k)=(1-daf)*s13oTZ(i,j,k)/c0+dtx/c0*mu13       &  
    *(v1TZ(i,j,k+1)-v1TZ(i,j,k))
   s13pTZ(i,j,k)=s13pTZ(i,j,k)+dtx*mu13*                    &  
    (rn*(v3TZ(i+1,j,k)-v3TZ(i,j,k))+rnn*(v3TZ(i+2,j,k)-          &  
    v3TZ(i-1,j,k)))
   s13TZ(i,j,k)=s13oTZ(i,j,k)+s13pTZ(i,j,k)
!
   s23oTZ(i,j,k)=(1-daf)*s23oTZ(i,j,k)/c0+dtx/c0*mu23       &  
    *(v2TZ(i,j,k+1)-v2TZ(i,j,k))
   s23pTZ(i,j,k)=s23pTZ(i,j,k)+dtx*mu23*                    &  
    (rn*(v3TZ(i,j,k)-v3TZ(i,j-1,k))+rnn*(v3TZ(i,j+1,k)-          &  
    v3TZ(i,j-2,k)))
   s23TZ(i,j,k)=s23oTZ(i,j,k)+s23pTZ(i,j,k)
!
  enddo
 enddo
enddo
!
return
end subroutine

! ###################################
Subroutine sPMLBZ
include 'triffy.dec'
do k=2,npm-1
 do j=3,ny-2
  do i=3,nx-2
   daf=DiBZS0(i,j,k)
   c0=1+daf
   s11oBZ(i,j,k)=(1-daf)/c0*s11oBZ(i,j,k)+dtx/c0*            &
    lamBZ(i,j,k)*(v3BZ(i,j,k)-v3BZ(i,j,k-1))
   s11pBZ(i,j,k)=s11pBZ(i,j,k)+dtx*((lamBZ(i,j,k)             &
    +2.*muBZ(i,j,k))*(rn*(v1BZ(i,j,k)-v1BZ(i-1,j,k))              &
    +rnn*(v1BZ(i+1,j,k)-v1BZ(i-2,j,k)))+lamBZ(i,j,k)*             &
    (rn*(v2BZ(i,j+1,k)-v2BZ(i,j,k))+rnn*(v2BZ(i,j+2,k)-           &
    v2BZ(i,j-1,k))))
   s11BZ(i,j,k)=s11pBZ(i,j,k)+s11oBZ(i,j,k)
! 
   s22oBZ(i,j,k)=(1-daf)/c0*s22oBZ(i,j,k)+dtx/c0*            &  
    lamBZ(i,j,k)*(v3BZ(i,j,k)-v3BZ(i,j,k-1))  
   s22pBZ(i,j,k)=s22pBZ(i,j,k)+dtx*(lamBZ(i,j,k)              & 
    *(rn*(v1BZ(i,j,k)-v1BZ(i-1,j,k))                             & 
    +rnn*(v1BZ(i+1,j,k)-v1BZ(i-2,j,k)))+(lamBZ(i,j,k)+            &     
    2*muBZ(i,j,k))*                                             & 
    (rn*(v2BZ(i,j+1,k)-v2BZ(i,j,k))+rnn*(v2BZ(i,j+2,k)-           & 
    v2BZ(i,j-1,k))) ) 
   s22BZ(i,j,k)=s22oBZ(i,j,k)+s22pBZ(i,j,k)  
! 
   s33oBZ(i,j,k)=(1-daf)/c0*s33oBZ(i,j,k)+dtx/c0*            &   
    (lamBZ(i,j,k)+2*muBZ(i,j,k))*(v3BZ(i,j,k)-v3BZ(i,j,k-1))   
   s33pBZ(i,j,k)=s33pBZ(i,j,k)+dtx*lamBZ(i,j,k)               &  
    *(rn*(v1BZ(i,j,k)-v1BZ(i-1,j,k)+v2BZ(i,j+1,k)-v2BZ(i,j,k))     &
    +rnn*(v1BZ(i+1,j,k)-v1BZ(i-2,j,k)+v2BZ(i,j+2,k)-              &
    v2BZ(i,j-1,k))) 
   s33BZ(i,j,k)=s33pBZ(i,j,k)+s33oBZ(i,j,k)             
!
   mu12=4./(1./muBZ(i,j,k)+1./muBZ(i,j-1,k)+1./muBZ(i+1,j-1,k)+  &
    1./muBZ(i+1,j,k))
   s12BZ(i,j,k)=s12BZ(i,j,k)+dtx*mu12*(rn*(v1BZ(i,j,k)-      &
    v1BZ(i,j-1,k)+v2BZ(i+1,j,k)-v2BZ(i,j,k))+rnn*(v1BZ(i,j+1,k)-  &
    v1BZ(i,j-2,k)+v2BZ(i+2,j,k)-v2BZ(i-1,j,k)))
  enddo 
 enddo 
enddo 
!
do k=1,npm-1
 do j=3,ny-2 
  do i=3,nx-2
   daf=DiBZS1(i,j,k)
   c0=1.+daf
   mu13=4./(1./muBZ(i,j,k)+1./muBZ(i,j,k+1)+1./muBZ(i+1,j,k)+    &
    1./muBZ(i+1,j,k+1))
   mu23=4./(1./muBZ(i,j,k)+1./muBZ(i,j,k+1)+1./muBZ(i,j-1,k)+    &
    1./muBZ(i,j-1,k+1))
!
   s13oBZ(i,j,k)=(1-daf)*s13oBZ(i,j,k)/c0+dtx/c0*mu13       &  
    *(v1BZ(i,j,k+1)-v1BZ(i,j,k))
   s13pBZ(i,j,k)=s13pBZ(i,j,k)+dtx*mu13*                    &  
    (rn*(v3BZ(i+1,j,k)-v3BZ(i,j,k))+rnn*(v3BZ(i+2,j,k)-          &  
    v3BZ(i-1,j,k)))
   s13BZ(i,j,k)=s13oBZ(i,j,k)+s13pBZ(i,j,k)
!
   s23oBZ(i,j,k)=(1-daf)*s23oBZ(i,j,k)/c0+dtx/c0*mu23       &  
    *(v2BZ(i,j,k+1)-v2BZ(i,j,k))
   s23pBZ(i,j,k)=s23pBZ(i,j,k)+dtx*mu23*                    &  
    (rn*(v3BZ(i,j,k)-v3BZ(i,j-1,k))+rnn*(v3BZ(i,j+1,k)-          &  
    v3BZ(i,j-2,k)))
   s23BZ(i,j,k)=s23oBZ(i,j,k)+s23pBZ(i,j,k)
!
  enddo
 enddo
enddo
!
return
end subroutine
! ##################################################
Subroutine sPMLTY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=3,nx-2
   daf=DiTYS0(i,j,k)
   c0=1+daf 
   s11oTY(i,j,k)=(1-daf)/c0*s11oTY(i,j,k)+dtx/c0*            & 
    lamTY(i,j,k)*(v2TY(i,j+1,k)-v2TY(i,j,k))
   s11pTY(i,j,k)=s11pTY(i,j,k)+dtx*((lamTY(i,j,k)             & 
    +2.*muTY(i,j,k))*(rn*(v1TY(i,j,k)-v1TY(i-1,j,k))              & 
    +rnn*(v1TY(i+1,j,k)-v1TY(i-2,j,k)))+lamTY(i,j,k)*             & 
    (rn*(v3TY(i,j,k)-v3TY(i,j,k-1))+rnn*(v3TY(i,j,k+1)-           & 
    v3TY(i,j,k-2))))
   s11TY(i,j,k)=s11pTY(i,j,k)+s11oTY(i,j,k)
!
   s22oTY(i,j,k)=(1-daf)/c0*s22oTY(i,j,k)+dtx/c0*            & 
    (lamTY(i,j,k)+2.*muTY(i,j,k))*(v2TY(i,j+1,k)-v2TY(i,j,k))
   s22pTY(i,j,k)=s22pTY(i,j,k)+dtx*(lamTY(i,j,k)              & 
    *(rn*(v1TY(i,j,k)-v1TY(i-1,j,k)+v3TY(i,j,k)-v3TY(i,j,k-1))     & 
    +rnn*(v1TY(i+1,j,k)-v1TY(i-2,j,k)+v3TY(i,j,k+1)               & 
    -v3TY(i,j,k-2))))
   s22TY(i,j,k)=s22oTY(i,j,k)+s22pTY(i,j,k)
!
   s33oTY(i,j,k)=(1-daf)/c0*s33oTY(i,j,k)+dtx/c0*            & 
    lamTY(i,j,k)*(v2TY(i,j+1,k)-v2TY(i,j,k))
   s33pTY(i,j,k)=s33pTY(i,j,k)+dtx*(lamTY(i,j,k)              & 
    *(rn*(v1TY(i,j,k)-v1TY(i-1,j,k))                             & 
    +rnn*(v1TY(i+1,j,k)-v1TY(i-2,j,k)))+(lamTY(i,j,k)+            & 
    2.*muTY(i,j,k))*(rn*(v3TY(i,j,k)-v3TY(i,j,k-1))+              & 
    rnn*(v3TY(i,j,k+1)-v3TY(i,j,k-2))))
   s33TY(i,j,k)=s33pTY(i,j,k)+s33oTY(i,j,k)
!
    mu13=4./(1./muTY(i,j,k)+1./muTY(i,j,k+1)+1./muTY(i+1,j,k)+    &
    1./muTY(i+1,j,k+1))
   s13TY(i,j,k)=s13TY(i,j,k)+dtx*mu13*(rn*(v1TY(i,j,k+1)-    &
    v1TY(i,j,k)+v3TY(i+1,j,k)-v3TY(i,j,k))+rnn*(v1TY(i,j,k+2)-    &
    v1TY(i,j,k-1)+v3TY(i+2,j,k)-v3TY(i-1,j,k)))
!
  enddo  
 enddo
enddo
!
do k=3,nz-2
 do j=2,npm
  do i=3,nx-2
   daf=DiTYS1(i,j,k)
   c0=1.+daf
   mu12=4./(1./muTY(i,j,k)+1./muTY(i,j-1,k)+1./muTY(i+1,j,k)+    &
    1./muTY(i+1,j-1,k))
   mu23=4./(1./muTY(i,j,k)+1./muTY(i,j,k+1)+1./muTY(i,j-1,k)+    &
    1./muTY(i,j-1,k+1))
!
   s12oTY(i,j,k)=(1-daf)*s12oTY(i,j,k)/c0+dtx/c0*mu12       &  
    *(v1TY(i,j,k)-v1TY(i,j-1,k))
   s12pTY(i,j,k)=s12pTY(i,j,k)+dtx*mu12*                    &  
    (rn*(v2TY(i+1,j,k)-v2TY(i,j,k))+rnn*(v2TY(i+2,j,k)-          &  
    v2TY(i-1,j,k)))
   s12TY(i,j,k)=s12oTY(i,j,k)+s12pTY(i,j,k)
!
   s23oTY(i,j,k)=(1-daf)*s23oTY(i,j,k)/c0+dtx/c0*mu23       &  
    *(v3TY(i,j,k)-v3TY(i,j-1,k))
   s23pTY(i,j,k)=s23pTY(i,j,k)+dtx*mu23*                    &  
    (rn*(v2TY(i,j,k+1)-v2TY(i,j,k))+rnn*(v2TY(i,j,k+2)-          &  
    v2TY(i,j,k-1)))
   s23TY(i,j,k)=s23oTY(i,j,k)+s23pTY(i,j,k)
!
  enddo
 enddo
enddo
!
return
end subroutine

! ##################################################
Subroutine sPMLBY
include 'triffy.dec'
do k=3,nz-2
 do j=1,npm-1
  do i=3,nx-2
   daf=DiBYS0(i,j,k)
   c0=1+daf 
   s11oBY(i,j,k)=(1-daf)/c0*s11oBY(i,j,k)+dtx/c0*            & 
    lamBY(i,j,k)*(v2BY(i,j+1,k)-v2BY(i,j,k))
   s11pBY(i,j,k)=s11pBY(i,j,k)+dtx*((lamBY(i,j,k)             & 
    +2.*muBY(i,j,k))*(rn*(v1BY(i,j,k)-v1BY(i-1,j,k))              & 
    +rnn*(v1BY(i+1,j,k)-v1BY(i-2,j,k)))+lamBY(i,j,k)*             & 
    (rn*(v3BY(i,j,k)-v3BY(i,j,k-1))+rnn*(v3BY(i,j,k+1)-           & 
    v3BY(i,j,k-2))))
   s11BY(i,j,k)=s11pBY(i,j,k)+s11oBY(i,j,k)
!
   s22oBY(i,j,k)=(1-daf)/c0*s22oBY(i,j,k)+dtx/c0*            & 
    (lamBY(i,j,k)+2.*muBY(i,j,k))*(v2BY(i,j+1,k)-v2BY(i,j,k))
   s22pBY(i,j,k)=s22pBY(i,j,k)+dtx*(lamBY(i,j,k)              & 
    *(rn*(v1BY(i,j,k)-v1BY(i-1,j,k)+v3BY(i,j,k)-v3BY(i,j,k-1))     & 
    +rnn*(v1BY(i+1,j,k)-v1BY(i-2,j,k)+v3BY(i,j,k+1)               & 
    -v3BY(i,j,k-2))))
   s22BY(i,j,k)=s22oBY(i,j,k)+s22pBY(i,j,k)
!
   s33oBY(i,j,k)=(1-daf)/c0*s33oBY(i,j,k)+dtx/c0*            & 
    lamBY(i,j,k)*(v2BY(i,j+1,k)-v2BY(i,j,k))
   s33pBY(i,j,k)=s33pBY(i,j,k)+dtx*(lamBY(i,j,k)              & 
    *(rn*(v1BY(i,j,k)-v1BY(i-1,j,k))                             & 
    +rnn*(v1BY(i+1,j,k)-v1BY(i-2,j,k)))+(lamBY(i,j,k)+            & 
    2.*muBY(i,j,k))*(rn*(v3BY(i,j,k)-v3BY(i,j,k-1))+              & 
    rnn*(v3BY(i,j,k+1)-v3BY(i,j,k-2))))
   s33BY(i,j,k)=s33pBY(i,j,k)+s33oBY(i,j,k)
!
    mu13=4./(1./muBY(i,j,k)+1./muBY(i,j,k+1)+1./muBY(i+1,j,k)+    &
    1./muBY(i+1,j,k+1))
   s13BY(i,j,k)=s13BY(i,j,k)+dtx*mu13*(rn*(v1BY(i,j,k+1)-    &
    v1BY(i,j,k)+v3BY(i+1,j,k)-v3BY(i,j,k))+rnn*(v1BY(i,j,k+2)-    &
    v1BY(i,j,k-1)+v3BY(i+2,j,k)-v3BY(i-1,j,k)))
!
  enddo  
 enddo
enddo
!
do k=3,nz-2
 do j=2,npm-1
  do i=3,nx-2
   daf=DiBYS1(i,j,k)
   c0=1.+daf
   mu12=4./(1./muBY(i,j,k)+1./muBY(i,j-1,k)+1./muBY(i+1,j,k)+    &
    1./muBY(i+1,j-1,k))
   mu23=4./(1./muBY(i,j,k)+1./muBY(i,j,k+1)+1./muBY(i,j-1,k)+    &
    1./muBY(i,j-1,k+1))
!
   s12oBY(i,j,k)=(1-daf)*s12oBY(i,j,k)/c0+dtx/c0*mu12       &  
    *(v1BY(i,j,k)-v1BY(i,j-1,k))
   s12pBY(i,j,k)=s12pBY(i,j,k)+dtx*mu12*                    &  
    (rn*(v2BY(i+1,j,k)-v2BY(i,j,k))+rnn*(v2BY(i+2,j,k)-          &  
    v2BY(i-1,j,k)))
   s12BY(i,j,k)=s12oBY(i,j,k)+s12pBY(i,j,k)
!
   s23oBY(i,j,k)=(1-daf)*s23oBY(i,j,k)/c0+dtx/c0*mu23       &  
    *(v3BY(i,j,k)-v3BY(i,j-1,k))
   s23pBY(i,j,k)=s23pBY(i,j,k)+dtx*mu23*                    &  
    (rn*(v2BY(i,j,k+1)-v2BY(i,j,k))+rnn*(v2BY(i,j,k+2)-          &  
    v2BY(i,j,k-1)))
   s23BY(i,j,k)=s23oBY(i,j,k)+s23pBY(i,j,k)
!
  enddo
 enddo
enddo
!
return
end subroutine
! ################################################
subroutine sPMLTX
include 'triffy.dec'
do k=3,nz-2
 do j=3,ny-2
  do i=2,npm
   daf=DiTXS0(i,j,k)
   c0=1+daf
   s11oTX(i,j,k)=(1-daf)/c0*s11oTX(i,j,k)+dtx/c0*            & 
    (lamTX(i,j,k)+2*muTX(i,j,k))*(v1TX(i,j,k)-v1TX(i-1,j,k))
   s11pTX(i,j,k)=s11pTX(i,j,k)+dtx*lamTX(i,j,k)*              & 
    (rn*(v2TX(i,j+1,k)-v2TX(i,j,k)+v3TX(i,j,k)-v3TX(i,j,k-1))+     &
    rnn*(v2TX(i,j+2,k)-v2TX(i,j-1,k)+v3TX(i,j,k+1)-v3TX(i,j,k-2)))
   s11TX(i,j,k)=s11pTX(i,j,k)+s11oTX(i,j,k)
!
   s22oTX(i,j,k)=(1-daf)/c0*s22oTX(i,j,k)+dtx/c0*            & 
    lamTX(i,j,k)*(v1TX(i,j,k)-v1TX(i-1,j,k))
   s22pTX(i,j,k)=s22pTX(i,j,k)+dtx*(lamTX(i,j,k)              & 
    *(rn*(v3TX(i,j,k)-v3TX(i,j,k-1))                             & 
    +rnn*(v3TX(i,j,k+1)-v3TX(i,j,k-2)))+(lamTX(i,j,k)+            & 
    2*muTX(i,j,k))*                                             &
    (rn*(v2TX(i,j+1,k)-v2TX(i,j,k))+rnn*(v2TX(i,j+2,k)-           & 
    v2TX(i,j-1,k))) )
   s22TX(i,j,k)=s22oTX(i,j,k)+s22pTX(i,j,k)
!
   s33oTX(i,j,k)=(1-daf)/c0*s33oTX(i,j,k)+dtx/c0*            & 
    lamTX(i,j,k)*(v1TX(i,j,k)-v1TX(i-1,j,k))
   s33pTX(i,j,k)=s33pTX(i,j,k)+dtx*(lamTX(i,j,k)              & 
    *(rn*(v2TX(i,j+1,k)-v2TX(i,j,k))+rnn*(v2TX(i,j+2,k)-          &
    v2TX(i,j-1,k)))+(lamTX(i,j,k)+2*muTX(i,j,k))*(rn*              &
    (v3TX(i,j,k)-v3TX(i,j,k-1))+rnn*(v3TX(i,j,k+1)-               &
    v3TX(i,j,k-2))))
   s33TX(i,j,k)=s33pTX(i,j,k)+s33oTX(i,j,k)
!
   mu23=4./(1./muTX(i,j,k)+1./muTX(i,j,k+1)+1./muTX(i,j-1,k)+    &
    1./muTX(i,j-1,k+1))
   s23TX(i,j,k)=s23TX(i,j,k)+dtx*mu23*(rn*(v2TX(i,j,k+1)-    &
    v2TX(i,j,k)+v3TX(i,j,k)-v3TX(i,j-1,k))+rnn*(v2TX(i,j,k+2)-    &
    v2TX(i,j,k-1)+v3TX(i,j+1,k)-v3TX(i,j-2,k)))
!
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=3,ny-2
  do i=2,npm-1
   daf=DiTXS1(i,j,k)
   c0=1.+daf
   mu13=4./(1./muTX(i,j,k)+1./muTX(i,j,k+1)+1./muTX(i+1,j,k)+    &
    1./muTX(i+1,j,k+1))
   mu12=4./(1./muTX(i,j,k)+1./muTX(i,j-1,k)+1./muTX(i+1,j,k)+    &
    1./muTX(i+1,j-1,k))
!
   s12oTX(i,j,k)=(1-daf)*s12oTX(i,j,k)/c0+dtx/c0*mu12       &  
    *(v2TX(i+1,j,k)-v2TX(i,j,k))
   s12pTX(i,j,k)=s12pTX(i,j,k)+dtx*mu12*                    &  
    (rn*(v1TX(i,j,k)-v1TX(i,j-1,k))+rnn*(v1TX(i,j+1,k)-          &  
    v1TX(i,j-2,k)))
   s12TX(i,j,k)=s12oTX(i,j,k)+s12pTX(i,j,k)
!
   s13oTX(i,j,k)=(1-daf)*s13oTX(i,j,k)/c0+dtx/c0*mu13       &  
    *(v3TX(i+1,j,k)-v3TX(i,j,k))
   s13pTX(i,j,k)=s13pTX(i,j,k)+dtx*mu13*                    &  
    (rn*(v1TX(i,j,k+1)-v1TX(i,j,k))+rnn*(v1TX(i,j,k+2)-          &  
    v1TX(i,j,k-1)))
   s13TX(i,j,k)=s13oTX(i,j,k)+s13pTX(i,j,k)
!
  enddo
 enddo
enddo
!
return
end subroutine
! #############################################
subroutine sPMLBX
include 'triffy.dec'
do k=3,nz-2
 do j=3,ny-2
  do i=2,npm-1
   daf=DiBXS0(i,j,k)
   c0=1+daf
   s11oBX(i,j,k)=(1-daf)/c0*s11oBX(i,j,k)+dtx/c0*            & 
    (lamBX(i,j,k)+2*muBX(i,j,k))*(v1BX(i,j,k)-v1BX(i-1,j,k))
   s11pBX(i,j,k)=s11pBX(i,j,k)+dtx*lamBX(i,j,k)*              & 
    (rn*(v2BX(i,j+1,k)-v2BX(i,j,k)+v3BX(i,j,k)-v3BX(i,j,k-1))+     &
    rnn*(v2BX(i,j+2,k)-v2BX(i,j-1,k)+v3BX(i,j,k+1)-v3BX(i,j,k-2)))
   s11BX(i,j,k)=s11pBX(i,j,k)+s11oBX(i,j,k)
!
   s22oBX(i,j,k)=(1-daf)/c0*s22oBX(i,j,k)+dtx/c0*            & 
    lamBX(i,j,k)*(v1BX(i,j,k)-v1BX(i-1,j,k))
   s22pBX(i,j,k)=s22pBX(i,j,k)+dtx*(lamBX(i,j,k)              & 
    *(rn*(v3BX(i,j,k)-v3BX(i,j,k-1))                             & 
    +rnn*(v3BX(i,j,k+1)-v3BX(i,j,k-2)))+(lamBX(i,j,k)+            & 
    2*muBX(i,j,k))*                                             &
    (rn*(v2BX(i,j+1,k)-v2BX(i,j,k))+rnn*(v2BX(i,j+2,k)-           & 
    v2BX(i,j-1,k))) )
   s22BX(i,j,k)=s22oBX(i,j,k)+s22pBX(i,j,k)
!
   s33oBX(i,j,k)=(1-daf)/c0*s33oBX(i,j,k)+dtx/c0*            & 
    lamBX(i,j,k)*(v1BX(i,j,k)-v1BX(i-1,j,k))
   s33pBX(i,j,k)=s33pBX(i,j,k)+dtx*(lamBX(i,j,k)              & 
    *(rn*(v2BX(i,j+1,k)-v2BX(i,j,k))+rnn*(v2BX(i,j+2,k)-          &
    v2BX(i,j-1,k)))+(lamBX(i,j,k)+2*muBX(i,j,k))*(rn*              &
    (v3BX(i,j,k)-v3BX(i,j,k-1))+rnn*(v3BX(i,j,k+1)-               &
    v3BX(i,j,k-2))))
   s33BX(i,j,k)=s33pBX(i,j,k)+s33oBX(i,j,k)
!
   mu23=4./(1./muBX(i,j,k)+1./muBX(i,j,k+1)+1./muBX(i,j-1,k)+    &
    1./muBX(i,j-1,k+1))
   s23BX(i,j,k)=s23BX(i,j,k)+dtx*mu23*(rn*(v2BX(i,j,k+1)-    &
    v2BX(i,j,k)+v3BX(i,j,k)-v3BX(i,j-1,k))+rnn*(v2BX(i,j,k+2)-    &
    v2BX(i,j,k-1)+v3BX(i,j+1,k)-v3BX(i,j-2,k)))
!
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=3,ny-2
  do i=1,npm-1
   daf=DiBXS1(i,j,k)
   c0=1.+daf
   mu13=4./(1./muBX(i,j,k)+1./muBX(i,j,k+1)+1./muBX(i+1,j,k)+    &
    1./muBX(i+1,j,k+1))
   mu12=4./(1./muBX(i,j,k)+1./muBX(i,j-1,k)+1./muBX(i+1,j,k)+    &
    1./muBX(i+1,j-1,k))
!
   s12oBX(i,j,k)=(1-daf)*s12oBX(i,j,k)/c0+dtx/c0*mu12       &  
    *(v2BX(i+1,j,k)-v2BX(i,j,k))
   s12pBX(i,j,k)=s12pBX(i,j,k)+dtx*mu12*                    &  
    (rn*(v1BX(i,j,k)-v1BX(i,j-1,k))+rnn*(v1BX(i,j+1,k)-          &  
    v1BX(i,j-2,k)))
   s12BX(i,j,k)=s12oBX(i,j,k)+s12pBX(i,j,k)
!
   s13oBX(i,j,k)=(1-daf)*s13oBX(i,j,k)/c0+dtx/c0*mu13       &  
    *(v3BX(i+1,j,k)-v3BX(i,j,k))
   s13pBX(i,j,k)=s13pBX(i,j,k)+dtx*mu13*                    &  
    (rn*(v1BX(i,j,k+1)-v1BX(i,j,k))+rnn*(v1BX(i,j,k+2)-          &  
    v1BX(i,j,k-1)))
   s13BX(i,j,k)=s13oBX(i,j,k)+s13pBX(i,j,k)
!
  enddo
 enddo
enddo
!
return
end subroutine
!###################################################
Subroutine sPMLTXTY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm
   dafx=DiTXS0(i,ny-2,k)
   c0x=1+dafx
   dafy=DiTYS0(nx-2,j,k)
   c0y=1+dafy
   s11xTXTY(i,j,k)=(1-dafx)/c0x*s11xTXTY(i,j,k)+dtx/c0x*       & 
    (lamTXTY(i,j,k)+2*muTXTY(i,j,k))*(v1TXTY(i,j,k)-v1TXTY(i-1,j,k))
   s11yTXTY(i,j,k)=(1-dafy)/c0y*s11yTXTY(i,j,k)+dtx/c0y*       & 
    lamTXTY(i,j,k)*(v2TXTY(i,j+1,k)-v2TXTY(i,j,k))
   s11zTXTY(i,j,k)=s11zTXTY(i,j,k)+dtx*lamTXTY(i,j,k)*           &
    (rn*(v3TXTY(i,j,k)-v3TXTY(i,j,k-1))+rnn*(v3TXTY(i,j,k+1)-        &
    v3TXTY(i,j,k-2)))   
   s11TXTY(i,j,k)=s11xTXTY(i,j,k)+s11yTXTY(i,j,k)+s11zTXTY(i,j,k)
! 
   s22xTXTY(i,j,k)=(1-dafx)/c0x*s22xTXTY(i,j,k)+dtx/c0x*       & 
    lamTXTY(i,j,k)*(v1TXTY(i,j,k)-v1TXTY(i-1,j,k))
   s22yTXTY(i,j,k)=(1-dafy)/c0y*s22yTXTY(i,j,k)+dtx/c0y*       & 
    (lamTXTY(i,j,k)+2.*muTXTY(i,j,k))*(v2TXTY(i,j+1,k)-v2TXTY(i,j,k))
   s22zTXTY(i,j,k)=s22zTXTY(i,j,k)+dtx*lamTXTY(i,j,k)*           &
    (rn*(v3TXTY(i,j,k)-v3TXTY(i,j,k-1))+rnn*(v3TXTY(i,j,k+1)-        & 
    v3TXTY(i,j,k-2)))
   s22TXTY(i,j,k)=s22xTXTY(i,j,k)+s22yTXTY(i,j,k)+s22zTXTY(i,j,k)
!
   s33xTXTY(i,j,k)=(1-dafx)/c0x*s33xTXTY(i,j,k)+dtx/c0x*       & 
    lamTXTY(i,j,k)*(v1TXTY(i,j,k)-v1TXTY(i-1,j,k))
   s33yTXTY(i,j,k)=(1-dafy)/c0y*s33yTXTY(i,j,k)+dtx/c0y*       & 
    lamTXTY(i,j,k)*(v2TXTY(i,j+1,k)-v2TXTY(i,j,k))
   s33zTXTY(i,j,k)=s33zTXTY(i,j,k)+dtx*(lamTXTY(i,j,k)+          &
   2.*muTXTY(i,j,k))*(rn*(v3TXTY(i,j,k)-v3TXTY(i,j,k-1))+            &
   rnn*(v3TXTY(i,j,k+1)-v3TXTY(i,j,k-2)))
   s33TXTY(i,j,k)=s33xTXTY(i,j,k)+s33yTXTY(i,j,k)+s33zTXTY(i,j,k)
!
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=2,npm
  do i=2,npm-1
   dafx=DiTXS1(i,ny-2,k)
   dafy=DiTYS1(nx-2,j,k)
   c0x=1.+dafx
   c0y=1.+dafy
   mu12=4./(1./muTXTY(i,j,k)+1./muTXTY(i,j-1,k)+1./muTXTY(i+1,j,k)+ &
    1./muTXTY(i+1,j-1,k))
!
   s12xTXTY(i,j,k)=(1-dafx)*s12xTXTY(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2TXTY(i+1,j,k)-v2TXTY(i,j,k))
   s12yTXTY(i,j,k)=(1-dafy)*s12yTXTY(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1TXTY(i,j,k)-v1TXTY(i,j-1,k))
   s12TXTY(i,j,k)=s12xTXTY(i,j,k)+s12yTXTY(i,j,k)
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,ny-2,k)
    c0x=1.+dafx
   mu13=4./(1./muTXTY(i,j,k)+1./muTXTY(i,j,k+1)+1./muTXTY(i+1,j,k)+ &
    1./muTXTY(i+1,j,k+1))
   s13xTXTY(i,j,k)=(1-dafx)*s13xTXTY(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3TXTY(i+1,j,k)-v3TXTY(i,j,k))
   s13zTXTY(i,j,k)=s13zTXTY(i,j,k)+dtx*mu13*                  &
    (rn*(v1TXTY(i,j,k+1)-v1TXTY(i,j,k))+rnn*(v1TXTY(i,j,k+2)-       &
    v1TXTY(i,j,k-1)))
   s13TXTY(i,j,k)=s13zTXTY(i,j,k)+s13xTXTY(i,j,k)
  enddo
 enddo
enddo
! 
do k=3,nz-2
 do j=2,npm
  do i=2,npm
   dafy=DiTYS1(nx-2,j,k)
   c0y=1.+dafy
   mu23=4./(1./muTXTY(i,j,k)+1./muTXTY(i,j,k+1)+1./muTXTY(i,j-1,k)+ &
    1./muTXTY(i,j-1,k+1))
   s23yTXTY(i,j,k)=(1-dafy)*s23yTXTY(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3TXTY(i,j,k)-v3TXTY(i,j-1,k))
   s23zTXTY(i,j,k)=s23zTXTY(i,j,k)+dtx*mu23*                  &  
    (rn*(v2TXTY(i,j,k+1)-v2TXTY(i,j,k))+rnn*(v2TXTY(i,j,k+2)-       &  
    v2TXTY(i,j,k-1)))
   s23TXTY(i,j,k)=s23zTXTY(i,j,k)+s23yTXTY(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
!###################################################
Subroutine sPMLTXBY
include 'triffy.dec'
do k=3,nz-2
 do j=1,npm-1
  do i=2,npm
   dafx=DiTXS0(i,3,k)
   c0x=1+dafx
   dafy=DiBYS0(nx-2,j,k)
   c0y=1+dafy
   s11xTXBY(i,j,k)=(1-dafx)/c0x*s11xTXBY(i,j,k)+dtx/c0x*       & 
    (lamTXBY(i,j,k)+2*muTXBY(i,j,k))*(v1TXBY(i,j,k)-v1TXBY(i-1,j,k))
   s11yTXBY(i,j,k)=(1-dafy)/c0y*s11yTXBY(i,j,k)+dtx/c0y*       & 
    lamTXBY(i,j,k)*(v2TXBY(i,j+1,k)-v2TXBY(i,j,k))
   s11zTXBY(i,j,k)=s11zTXBY(i,j,k)+dtx*lamTXBY(i,j,k)*           &
    (rn*(v3TXBY(i,j,k)-v3TXBY(i,j,k-1))+rnn*(v3TXBY(i,j,k+1)-        &
    v3TXBY(i,j,k-2)))   
   s11TXBY(i,j,k)=s11xTXBY(i,j,k)+s11yTXBY(i,j,k)+s11zTXBY(i,j,k)
! 
   s22xTXBY(i,j,k)=(1-dafx)/c0x*s22xTXBY(i,j,k)+dtx/c0x*       & 
    lamTXBY(i,j,k)*(v1TXBY(i,j,k)-v1TXBY(i-1,j,k))
   s22yTXBY(i,j,k)=(1-dafy)/c0y*s22yTXBY(i,j,k)+dtx/c0y*       & 
    (lamTXBY(i,j,k)+2.*muTXBY(i,j,k))*(v2TXBY(i,j+1,k)-v2TXBY(i,j,k))
   s22zTXBY(i,j,k)=s22zTXBY(i,j,k)+dtx*lamTXBY(i,j,k)*           &
    (rn*(v3TXBY(i,j,k)-v3TXBY(i,j,k-1))+rnn*(v3TXBY(i,j,k+1)-        & 
    v3TXBY(i,j,k-2)))
   s22TXBY(i,j,k)=s22xTXBY(i,j,k)+s22yTXBY(i,j,k)+s22zTXBY(i,j,k)
!
   s33xTXBY(i,j,k)=(1-dafx)/c0x*s33xTXBY(i,j,k)+dtx/c0x*       & 
    lamTXBY(i,j,k)*(v1TXBY(i,j,k)-v1TXBY(i-1,j,k))
   s33yTXBY(i,j,k)=(1-dafy)/c0y*s33yTXBY(i,j,k)+dtx/c0y*       & 
    lamTXBY(i,j,k)*(v2TXBY(i,j+1,k)-v2TXBY(i,j,k))
   s33zTXBY(i,j,k)=s33zTXBY(i,j,k)+dtx*(lamTXBY(i,j,k)+          &
   2.*muTXBY(i,j,k))*(rn*(v3TXBY(i,j,k)-v3TXBY(i,j,k-1))+            &
   rnn*(v3TXBY(i,j,k+1)-v3TXBY(i,j,k-2)))
   s33TXBY(i,j,k)=s33xTXBY(i,j,k)+s33yTXBY(i,j,k)+s33zTXBY(i,j,k)
!
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,3,k)
   dafy=DiBYS1(nx-2,j,k)
   c0x=1.+dafx
   c0y=1.+dafy
   mu12=4./(1./muTXBY(i,j,k)+1./muTXBY(i,j-1,k)+1./muTXBY(i+1,j,k)+ &
    1./muTXBY(i+1,j-1,k))
!
   s12xTXBY(i,j,k)=(1-dafx)*s12xTXBY(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2TXBY(i+1,j,k)-v2TXBY(i,j,k))
   s12yTXBY(i,j,k)=(1-dafy)*s12yTXBY(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1TXBY(i,j,k)-v1TXBY(i,j-1,k))
   s12TXBY(i,j,k)=s12xTXBY(i,j,k)+s12yTXBY(i,j,k)
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=1,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,3,k)
    c0x=1.+dafx
   mu13=4./(1./muTXBY(i,j,k)+1./muTXBY(i,j,k+1)+1./muTXBY(i+1,j,k)+ &
    1./muTXBY(i+1,j,k+1))
   s13xTXBY(i,j,k)=(1-dafx)*s13xTXBY(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3TXBY(i+1,j,k)-v3TXBY(i,j,k))
   s13zTXBY(i,j,k)=s13zTXBY(i,j,k)+dtx*mu13*                  &
    (rn*(v1TXBY(i,j,k+1)-v1TXBY(i,j,k))+rnn*(v1TXBY(i,j,k+2)-       &
    v1TXBY(i,j,k-1)))
   s13TXBY(i,j,k)=s13zTXBY(i,j,k)+s13xTXBY(i,j,k)
  enddo
 enddo
enddo
! 
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm
   dafy=DiBYS1(nx-2,j,k)
   c0y=1.+dafy
   mu23=4./(1./muTXBY(i,j,k)+1./muTXBY(i,j,k+1)+1./muTXBY(i,j-1,k)+ &
    1./muTXBY(i,j-1,k+1))
   s23yTXBY(i,j,k)=(1-dafy)*s23yTXBY(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3TXBY(i,j,k)-v3TXBY(i,j-1,k))
   s23zTXBY(i,j,k)=s23zTXBY(i,j,k)+dtx*mu23*                  &  
    (rn*(v2TXBY(i,j,k+1)-v2TXBY(i,j,k))+rnn*(v2TXBY(i,j,k+2)-       &  
    v2TXBY(i,j,k-1)))
   s23TXBY(i,j,k)=s23zTXBY(i,j,k)+s23yTXBY(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
!###################################################
Subroutine sPMLBXTY
include 'triffy.dec'
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiBXS0(i,ny-2,k)
   c0x=1+dafx
   dafy=DiTYS0(3,j,k)
   c0y=1+dafy
   s11xBXTY(i,j,k)=(1-dafx)/c0x*s11xBXTY(i,j,k)+dtx/c0x*       & 
    (lamBXTY(i,j,k)+2*muBXTY(i,j,k))*(v1BXTY(i,j,k)-v1BXTY(i-1,j,k))
   s11yBXTY(i,j,k)=(1-dafy)/c0y*s11yBXTY(i,j,k)+dtx/c0y*       & 
    lamBXTY(i,j,k)*(v2BXTY(i,j+1,k)-v2BXTY(i,j,k))
   s11zBXTY(i,j,k)=s11zBXTY(i,j,k)+dtx*lamBXTY(i,j,k)*           &
    (rn*(v3BXTY(i,j,k)-v3BXTY(i,j,k-1))+rnn*(v3BXTY(i,j,k+1)-        &
    v3BXTY(i,j,k-2)))   
   s11BXTY(i,j,k)=s11xBXTY(i,j,k)+s11yBXTY(i,j,k)+s11zBXTY(i,j,k)
! 
   s22xBXTY(i,j,k)=(1-dafx)/c0x*s22xBXTY(i,j,k)+dtx/c0x*       & 
    lamBXTY(i,j,k)*(v1BXTY(i,j,k)-v1BXTY(i-1,j,k))
   s22yBXTY(i,j,k)=(1-dafy)/c0y*s22yBXTY(i,j,k)+dtx/c0y*       & 
    (lamBXTY(i,j,k)+2.*muBXTY(i,j,k))*(v2BXTY(i,j+1,k)-v2BXTY(i,j,k))
   s22zBXTY(i,j,k)=s22zBXTY(i,j,k)+dtx*lamBXTY(i,j,k)*           &
    (rn*(v3BXTY(i,j,k)-v3BXTY(i,j,k-1))+rnn*(v3BXTY(i,j,k+1)-        & 
    v3BXTY(i,j,k-2)))
   s22BXTY(i,j,k)=s22xBXTY(i,j,k)+s22yBXTY(i,j,k)+s22zBXTY(i,j,k)
!
   s33xBXTY(i,j,k)=(1-dafx)/c0x*s33xBXTY(i,j,k)+dtx/c0x*       & 
    lamBXTY(i,j,k)*(v1BXTY(i,j,k)-v1BXTY(i-1,j,k))
   s33yBXTY(i,j,k)=(1-dafy)/c0y*s33yBXTY(i,j,k)+dtx/c0y*       & 
    lamBXTY(i,j,k)*(v2BXTY(i,j+1,k)-v2BXTY(i,j,k))
   s33zBXTY(i,j,k)=s33zBXTY(i,j,k)+dtx*(lamBXTY(i,j,k)+          &
   2.*muBXTY(i,j,k))*(rn*(v3BXTY(i,j,k)-v3BXTY(i,j,k-1))+            &
   rnn*(v3BXTY(i,j,k+1)-v3BXTY(i,j,k-2)))
   s33BXTY(i,j,k)=s33xBXTY(i,j,k)+s33yBXTY(i,j,k)+s33zBXTY(i,j,k)
!
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=2,npm
  do i=1,npm-1
   dafx=DiBXS1(i,ny-2,k)
   dafy=DiTYS1(3,j,k)
   c0x=1.+dafx
   c0y=1.+dafy
   mu12=4./(1./muBXTY(i,j,k)+1./muBXTY(i,j-1,k)+1./muBXTY(i+1,j,k)+ &
    1./muBXTY(i+1,j-1,k))
!
   s12xBXTY(i,j,k)=(1-dafx)*s12xBXTY(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2BXTY(i+1,j,k)-v2BXTY(i,j,k))
   s12yBXTY(i,j,k)=(1-dafy)*s12yBXTY(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1BXTY(i,j,k)-v1BXTY(i,j-1,k))
   s12BXTY(i,j,k)=s12xBXTY(i,j,k)+s12yBXTY(i,j,k)
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=2,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,ny-2,k)
    c0x=1.+dafx
   mu13=4./(1./muBXTY(i,j,k)+1./muBXTY(i,j,k+1)+1./muBXTY(i+1,j,k)+ &
    1./muBXTY(i+1,j,k+1))
   s13xBXTY(i,j,k)=(1-dafx)*s13xBXTY(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3BXTY(i+1,j,k)-v3BXTY(i,j,k))
   s13zBXTY(i,j,k)=s13zBXTY(i,j,k)+dtx*mu13*                  &
    (rn*(v1BXTY(i,j,k+1)-v1BXTY(i,j,k))+rnn*(v1BXTY(i,j,k+2)-       &
    v1BXTY(i,j,k-1)))
   s13BXTY(i,j,k)=s13zBXTY(i,j,k)+s13xBXTY(i,j,k)
  enddo
 enddo
enddo
! 
do k=3,nz-2
 do j=2,npm
  do i=2,npm-1
   dafy=DiTYS1(nx-2,j,k)
   c0y=1.+dafy
   mu23=4./(1./muBXTY(i,j,k)+1./muBXTY(i,j,k+1)+1./muBXTY(i,j-1,k)+ &
    1./muBXTY(i,j-1,k+1))
   s23yBXTY(i,j,k)=(1-dafy)*s23yBXTY(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3BXTY(i,j,k)-v3BXTY(i,j-1,k))
   s23zBXTY(i,j,k)=s23zBXTY(i,j,k)+dtx*mu23*                  &  
    (rn*(v2BXTY(i,j,k+1)-v2BXTY(i,j,k))+rnn*(v2BXTY(i,j,k+2)-       &  
    v2BXTY(i,j,k-1)))
   s23BXTY(i,j,k)=s23zBXTY(i,j,k)+s23yBXTY(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
!###################################################
Subroutine sPMLBXBY
include 'triffy.dec'
do k=3,nz-2
 do j=1,npm-1
  do i=2,npm-1
   dafx=DiBXS0(i,3,k)
   c0x=1+dafx
   dafy=DiBYS0(3,j,k)
   c0y=1+dafy
   s11xBXBY(i,j,k)=(1-dafx)/c0x*s11xBXBY(i,j,k)+dtx/c0x*       & 
    (lamBXBY(i,j,k)+2*muBXBY(i,j,k))*(v1BXBY(i,j,k)-v1BXBY(i-1,j,k))
   s11yBXBY(i,j,k)=(1-dafy)/c0y*s11yBXBY(i,j,k)+dtx/c0y*       & 
    lamBXBY(i,j,k)*(v2BXBY(i,j+1,k)-v2BXBY(i,j,k))
   s11zBXBY(i,j,k)=s11zBXBY(i,j,k)+dtx*lamBXBY(i,j,k)*           &
    (rn*(v3BXBY(i,j,k)-v3BXBY(i,j,k-1))+rnn*(v3BXBY(i,j,k+1)-        &
    v3BXBY(i,j,k-2)))   
   s11BXBY(i,j,k)=s11xBXBY(i,j,k)+s11yBXBY(i,j,k)+s11zBXBY(i,j,k)
! 
   s22xBXBY(i,j,k)=(1-dafx)/c0x*s22xBXBY(i,j,k)+dtx/c0x*       & 
    lamBXBY(i,j,k)*(v1BXBY(i,j,k)-v1BXBY(i-1,j,k))
   s22yBXBY(i,j,k)=(1-dafy)/c0y*s22yBXBY(i,j,k)+dtx/c0y*       & 
    (lamBXBY(i,j,k)+2.*muBXBY(i,j,k))*(v2BXBY(i,j+1,k)-v2BXBY(i,j,k))
   s22zBXBY(i,j,k)=s22zBXBY(i,j,k)+dtx*lamBXBY(i,j,k)*           &
    (rn*(v3BXBY(i,j,k)-v3BXBY(i,j,k-1))+rnn*(v3BXBY(i,j,k+1)-        & 
    v3BXBY(i,j,k-2)))
   s22BXBY(i,j,k)=s22xBXBY(i,j,k)+s22yBXBY(i,j,k)+s22zBXBY(i,j,k)
!
   s33xBXBY(i,j,k)=(1-dafx)/c0x*s33xBXBY(i,j,k)+dtx/c0x*       & 
    lamBXBY(i,j,k)*(v1BXBY(i,j,k)-v1BXBY(i-1,j,k))
   s33yBXBY(i,j,k)=(1-dafy)/c0y*s33yBXBY(i,j,k)+dtx/c0y*       & 
    lamBXBY(i,j,k)*(v2BXBY(i,j+1,k)-v2BXBY(i,j,k))
   s33zBXBY(i,j,k)=s33zBXBY(i,j,k)+dtx*(lamBXBY(i,j,k)+          &
   2.*muBXBY(i,j,k))*(rn*(v3BXBY(i,j,k)-v3BXBY(i,j,k-1))+            &
   rnn*(v3BXBY(i,j,k+1)-v3BXBY(i,j,k-2)))
   s33BXBY(i,j,k)=s33xBXBY(i,j,k)+s33yBXBY(i,j,k)+s33zBXBY(i,j,k)
!
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=2,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,3,k)
   dafy=DiBYS1(3,j,k)
   c0x=1.+dafx
   c0y=1.+dafy
   mu12=4./(1./muBXBY(i,j,k)+1./muBXBY(i,j-1,k)+1./muBXBY(i+1,j,k)+ &
    1./muBXBY(i+1,j-1,k))
!
   s12xBXBY(i,j,k)=(1-dafx)*s12xBXBY(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2BXBY(i+1,j,k)-v2BXBY(i,j,k))
   s12yBXBY(i,j,k)=(1-dafy)*s12yBXBY(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1BXBY(i,j,k)-v1BXBY(i,j-1,k))
   s12BXBY(i,j,k)=s12xBXBY(i,j,k)+s12yBXBY(i,j,k)
  enddo
 enddo
enddo
!
do k=3,nz-2
 do j=1,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,3,k)
    c0x=1.+dafx
   mu13=4./(1./muBXBY(i,j,k)+1./muBXBY(i,j,k+1)+1./muBXBY(i+1,j,k)+ &
    1./muBXBY(i+1,j,k+1))
   s13xBXBY(i,j,k)=(1-dafx)*s13xBXBY(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3BXBY(i+1,j,k)-v3BXBY(i,j,k))
   s13zBXBY(i,j,k)=s13zBXBY(i,j,k)+dtx*mu13*                  &
    (rn*(v1BXBY(i,j,k+1)-v1BXBY(i,j,k))+rnn*(v1BXBY(i,j,k+2)-       &
    v1BXBY(i,j,k-1)))
   s13BXBY(i,j,k)=s13zBXBY(i,j,k)+s13xBXBY(i,j,k)
  enddo
 enddo
enddo
! 
do k=3,nz-2
 do j=2,npm-1
  do i=2,npm-1
   dafy=DiBYS1(3,j,k)
   c0y=1.+dafy
   mu23=4./(1./muBXBY(i,j,k)+1./muBXBY(i,j,k+1)+1./muBXBY(i,j-1,k)+ &
    1./muBXBY(i,j-1,k+1))
   s23yBXBY(i,j,k)=(1-dafy)*s23yBXBY(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3BXBY(i,j,k)-v3BXBY(i,j-1,k))
   s23zBXBY(i,j,k)=s23zBXBY(i,j,k)+dtx*mu23*                  &  
    (rn*(v2BXBY(i,j,k+1)-v2BXBY(i,j,k))+rnn*(v2BXBY(i,j,k+2)-       &  
    v2BXBY(i,j,k-1)))
   s23BXBY(i,j,k)=s23zBXBY(i,j,k)+s23yBXBY(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! #######################################################
subroutine sPMLTXTZ
include 'triffy.dec'
do k=2,npm
 do j=3,ny-2
  do i=2,npm
   dafx=DiTXS0(i,j,nz-2)
   dafz=DiTZS0(nx-2,j,k)
   c0x=1+dafx
   c0z=1+dafz
   s11xTXTZ(i,j,k)=(1-dafx)/c0x*s11xTXTZ(i,j,k)+dtx/c0x*       & 
    (lamTXTZ(i,j,k)+2*muTXTZ(i,j,k))*(v1TXTZ(i,j,k)-v1TXTZ(i-1,j,k))
   s11yTXTZ(i,j,k)=s11yTXTZ(i,j,k)+dtx*lamTXTZ(i,j,k)*           &
    (rn*(v2TXTZ(i,j+1,k)-v2TXTZ(i,j,k))+rnn*(v2TXTZ(i,j+2,k)-        &
    v2TXTZ(i,j-1,k)))
   s11zTXTZ(i,j,k)=(1-dafz)/c0z*s11zTXTZ(i,j,k)+dtx/c0z*       & 
    lamTXTZ(i,j,k)*(v3TXTZ(i,j,k)-v3TXTZ(i,j,k-1))
   s11TXTZ(i,j,k)=s11xTXTZ(i,j,k)+s11yTXTZ(i,j,k)+s11zTXTZ(i,j,k)
!
   s22xTXTZ(i,j,k)=(1-dafx)/c0x*s22xTXTZ(i,j,k)+dtx/c0x*       & 
    lamTXTZ(i,j,k)*(v1TXTZ(i,j,k)-v1TXTZ(i-1,j,k))
   s22yTXTZ(i,j,k)=s22yTXTZ(i,j,k)+dtx*(lamTXTZ(i,j,k)+          &
    2.*muTXTZ(i,j,k))*                                           &
    (rn*(v2TXTZ(i,j+1,k)-v2TXTZ(i,j,k))+rnn*(v2TXTZ(i,j+2,k)-        & 
    v2TXTZ(i,j-1,k))) 
   s22zTXTZ(i,j,k)=(1-dafz)/c0z*s22zTXTZ(i,j,k)+dtx/c0z*       & 
    lamTXTZ(i,j,k)*(v3TXTZ(i,j,k)-v3TXTZ(i,j,k-1)) 
   s22TXTZ(i,j,k)=s22xTXTZ(i,j,k)+s22yTXTZ(i,j,k)+s22zTXTZ(i,j,k)
!
   s33xTXTZ(i,j,k)=(1-dafx)/c0x*s33xTXTZ(i,j,k)+dtx/c0x*       & 
    lamTXTZ(i,j,k)*(v1TXTZ(i,j,k)-v1TXTZ(i-1,j,k))
   s33yTXTZ(i,j,k)=s33yTXTZ(i,j,k)+dtx*lamTXTZ(i,j,k)*           &
    (rn*(v2TXTZ(i,j+1,k)-v2TXTZ(i,j,k))+rnn*(v2TXTZ(i,j+2,k)-        & 
    v2TXTZ(i,j-1,k))) 
   s33zTXTZ(i,j,k)=(1-dafz)/c0z*s33zTXTZ(i,j,k)+dtx/c0z*       &
   (lamTXTZ(i,j,k)+2*muTXTZ(i,j,k))*(v3TXTZ(i,j,k)-v3TXTZ(i,j,k-1))
   s33TXTZ(i,j,k)=s33xTXTZ(i,j,k)+s33yTXTZ(i,j,k)+s33zTXTZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm
 do j=3,ny-2
  do i=2,npm-1
   dafx=DiTXS1(i,j,nz-2)
   c0x=1.+dafx
   mu12=4./(1./muTXTZ(i,j,k)+1./muTXTZ(i,j-1,k)+1./muTXTZ(i+1,j,k)+ &
    1./muTXTZ(i+1,j-1,k))
  s12xTXTZ(i,j,k)=(1-dafx)*s12xTXTZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2TXTZ(i+1,j,k)-v2TXTZ(i,j,k))
   s12yTXTZ(i,j,k)=s12yTXTZ(i,j,k)+dtx*mu12*  &
    (rn*(v1TXTZ(i,j,k)-v1TXTZ(i,j-1,k))+rnn*(v1TXTZ(i,j+1,k)-       &
    v1TXTZ(i,j-2,k)))
   s12TXTZ(i,j,k)=s12xTXTZ(i,j,k)+s12yTXTZ(i,j,k)
  enddo
 enddo
enddo

do k=2,npm-1
 do j=3,ny-2
  do i=2,npm-1
   dafx=DiTXS1(i,j,nz-2)
   dafz=DiTZS1(nx-2,j,k)
   c0x=1.+dafx
   c0z=1.+dafz
   mu13=4./(1./muTXTZ(i,j,k)+1./muTXTZ(i,j,k+1)+1./muTXTZ(i+1,j,k)+ &
    1./muTXTZ(i+1,j,k+1))
   s13xTXTZ(i,j,k)=(1-dafx)*s13xTXTZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3TXTZ(i+1,j,k)-v3TXTZ(i,j,k))
   s13zTXTZ(i,j,k)=(1-dafz)*s13zTXTZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1TXTZ(i,j,k+1)-v1TXTZ(i,j,k))
   s13TXTZ(i,j,k)=s13xTXTZ(i,j,k)+s13zTXTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm
   dafz=DiTZS1(nx-2,j,k)
   c0z=1.+dafz
   mu23=4./(1./muTXTZ(i,j,k)+1./muTXTZ(i,j,k+1)+1./muTXTZ(i,j-1,k)+ &
    1./muTXTZ(i,j-1,k+1))
   s23yTXTZ(i,j,k)=s23yTXTZ(i,j,k)+dtx*mu23*    &
    (rn*(v3TXTZ(i,j,k)-v3TXTZ(i,j-1,k))+rnn*(v3TXTZ(i,j+1,k)-       & 
    v3TXTZ(i,j-2,k)))
   s23zTXTZ(i,j,k)=(1-dafz)*s23zTXTZ(i,j,k)/c0z+dtx/c0z*mu23  &  
    *(v2TXTZ(i,j,k+1)-v2TXTZ(i,j,k))
   s23TXTZ(i,j,k)=s23yTXTZ(i,j,k)+s23zTXTZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! #######################################################
subroutine sPMLTXBZ
include 'triffy.dec'
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm
   dafx=DiTXS0(i,j,3)
   dafz=DiBZS0(nx-2,j,k)
   c0x=1+dafx
   c0z=1+dafz
   s11xTXBZ(i,j,k)=(1-dafx)/c0x*s11xTXBZ(i,j,k)+dtx/c0x*       & 
    (lamTXBZ(i,j,k)+2*muTXBZ(i,j,k))*(v1TXBZ(i,j,k)-v1TXBZ(i-1,j,k))
   s11yTXBZ(i,j,k)=s11yTXBZ(i,j,k)+dtx*lamTXBZ(i,j,k)*           &
    (rn*(v2TXBZ(i,j+1,k)-v2TXBZ(i,j,k))+rnn*(v2TXBZ(i,j+2,k)-        &
    v2TXBZ(i,j-1,k)))
   s11zTXBZ(i,j,k)=(1-dafz)/c0z*s11zTXBZ(i,j,k)+dtx/c0z*       & 
    lamTXBZ(i,j,k)*(v3TXBZ(i,j,k)-v3TXBZ(i,j,k-1))
   s11TXBZ(i,j,k)=s11xTXBZ(i,j,k)+s11yTXBZ(i,j,k)+s11zTXBZ(i,j,k)
!
   s22xTXBZ(i,j,k)=(1-dafx)/c0x*s22xTXBZ(i,j,k)+dtx/c0x*       & 
    lamTXBZ(i,j,k)*(v1TXBZ(i,j,k)-v1TXBZ(i-1,j,k))
   s22yTXBZ(i,j,k)=s22yTXBZ(i,j,k)+dtx*(lamTXBZ(i,j,k)+          &
    2.*muTXBZ(i,j,k))*                                           &
    (rn*(v2TXBZ(i,j+1,k)-v2TXBZ(i,j,k))+rnn*(v2TXBZ(i,j+2,k)-        & 
    v2TXBZ(i,j-1,k))) 
   s22zTXBZ(i,j,k)=(1-dafz)/c0z*s22zTXBZ(i,j,k)+dtx/c0z*       & 
    lamTXBZ(i,j,k)*(v3TXBZ(i,j,k)-v3TXBZ(i,j,k-1)) 
   s22TXBZ(i,j,k)=s22xTXBZ(i,j,k)+s22yTXBZ(i,j,k)+s22zTXBZ(i,j,k)
!
   s33xTXBZ(i,j,k)=(1-dafx)/c0x*s33xTXBZ(i,j,k)+dtx/c0x*       & 
    lamTXBZ(i,j,k)*(v1TXBZ(i,j,k)-v1TXBZ(i-1,j,k))
   s33yTXBZ(i,j,k)=s33yTXBZ(i,j,k)+dtx*lamTXBZ(i,j,k)*           &
    (rn*(v2TXBZ(i,j+1,k)-v2TXBZ(i,j,k))+rnn*(v2TXBZ(i,j+2,k)-        & 
    v2TXBZ(i,j-1,k))) 
   s33zTXBZ(i,j,k)=(1-dafz)/c0z*s33zTXBZ(i,j,k)+dtx/c0z*       &
   (lamTXBZ(i,j,k)+2*muTXBZ(i,j,k))*(v3TXBZ(i,j,k)-v3TXBZ(i,j,k-1))
   s33TXBZ(i,j,k)=s33xTXBZ(i,j,k)+s33yTXBZ(i,j,k)+s33zTXBZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm-1
   dafx=DiTXS1(i,j,3)
   c0x=1.+dafx
   mu12=4./(1./muTXBZ(i,j,k)+1./muTXBZ(i,j-1,k)+1./muTXBZ(i+1,j,k)+ &
    1./muTXBZ(i+1,j-1,k))
  s12xTXBZ(i,j,k)=(1-dafx)*s12xTXBZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2TXBZ(i+1,j,k)-v2TXBZ(i,j,k))
   s12yTXBZ(i,j,k)=s12yTXBZ(i,j,k)+dtx*mu12*  &
    (rn*(v1TXBZ(i,j,k)-v1TXBZ(i,j-1,k))+rnn*(v1TXBZ(i,j+1,k)-       &
    v1TXBZ(i,j-2,k)))
   s12TXBZ(i,j,k)=s12xTXBZ(i,j,k)+s12yTXBZ(i,j,k)
  enddo
 enddo
enddo

do k=1,npm-1
 do j=3,ny-2
  do i=2,npm-1
   dafx=DiTXS1(i,j,3)
   dafz=DiBZS1(nx-2,j,k)
   c0x=1.+dafx
   c0z=1.+dafz
   mu13=4./(1./muTXBZ(i,j,k)+1./muTXBZ(i,j,k+1)+1./muTXBZ(i+1,j,k)+ &
    1./muTXBZ(i+1,j,k+1))
   s13xTXBZ(i,j,k)=(1-dafx)*s13xTXBZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3TXBZ(i+1,j,k)-v3TXBZ(i,j,k))
   s13zTXBZ(i,j,k)=(1-dafz)*s13zTXBZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1TXBZ(i,j,k+1)-v1TXBZ(i,j,k))
   s13TXBZ(i,j,k)=s13xTXBZ(i,j,k)+s13zTXBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=3,ny-2
  do i=2,npm
   dafz=DiBZS1(nx-2,j,k)
   c0z=1.+dafz
   mu23=4./(1./muTXBZ(i,j,k)+1./muTXBZ(i,j,k+1)+1./muTXBZ(i,j-1,k)+ &
    1./muTXBZ(i,j-1,k+1))
   s23yTXBZ(i,j,k)=s23yTXBZ(i,j,k)+dtx*mu23*    &
    (rn*(v3TXBZ(i,j,k)-v3TXBZ(i,j-1,k))+rnn*(v3TXBZ(i,j+1,k)-       & 
    v3TXBZ(i,j-2,k)))
   s23zTXBZ(i,j,k)=(1-dafz)*s23zTXBZ(i,j,k)/c0z+dtx/c0z*mu23  &  
    *(v2TXBZ(i,j,k+1)-v2TXBZ(i,j,k))
   s23TXBZ(i,j,k)=s23yTXBZ(i,j,k)+s23zTXBZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! #######################################################
subroutine sPMLBXTZ
include 'triffy.dec'
do k=2,npm
 do j=3,ny-2
  do i=2,npm-1
   dafx=DiBXS0(i,j,nz-2)
   dafz=DiTZS0(3,j,k)
   c0x=1+dafx
   c0z=1+dafz
   s11xBXTZ(i,j,k)=(1-dafx)/c0x*s11xBXTZ(i,j,k)+dtx/c0x*       & 
    (lamBXTZ(i,j,k)+2*muBXTZ(i,j,k))*(v1BXTZ(i,j,k)-v1BXTZ(i-1,j,k))
   s11yBXTZ(i,j,k)=s11yBXTZ(i,j,k)+dtx*lamBXTZ(i,j,k)*           &
    (rn*(v2BXTZ(i,j+1,k)-v2BXTZ(i,j,k))+rnn*(v2BXTZ(i,j+2,k)-        &
    v2BXTZ(i,j-1,k)))
   s11zBXTZ(i,j,k)=(1-dafz)/c0z*s11zBXTZ(i,j,k)+dtx/c0z*       & 
    lamBXTZ(i,j,k)*(v3BXTZ(i,j,k)-v3BXTZ(i,j,k-1))
   s11BXTZ(i,j,k)=s11xBXTZ(i,j,k)+s11yBXTZ(i,j,k)+s11zBXTZ(i,j,k)
!
   s22xBXTZ(i,j,k)=(1-dafx)/c0x*s22xBXTZ(i,j,k)+dtx/c0x*       & 
    lamBXTZ(i,j,k)*(v1BXTZ(i,j,k)-v1BXTZ(i-1,j,k))
   s22yBXTZ(i,j,k)=s22yBXTZ(i,j,k)+dtx*(lamBXTZ(i,j,k)+          &
    2.*muBXTZ(i,j,k))*                                           &
    (rn*(v2BXTZ(i,j+1,k)-v2BXTZ(i,j,k))+rnn*(v2BXTZ(i,j+2,k)-        & 
    v2BXTZ(i,j-1,k))) 
   s22zBXTZ(i,j,k)=(1-dafz)/c0z*s22zBXTZ(i,j,k)+dtx/c0z*       & 
    lamBXTZ(i,j,k)*(v3BXTZ(i,j,k)-v3BXTZ(i,j,k-1)) 
   s22BXTZ(i,j,k)=s22xBXTZ(i,j,k)+s22yBXTZ(i,j,k)+s22zBXTZ(i,j,k)
!
   s33xBXTZ(i,j,k)=(1-dafx)/c0x*s33xBXTZ(i,j,k)+dtx/c0x*       & 
    lamBXTZ(i,j,k)*(v1BXTZ(i,j,k)-v1BXTZ(i-1,j,k))
   s33yBXTZ(i,j,k)=s33yBXTZ(i,j,k)+dtx*lamBXTZ(i,j,k)*           &
    (rn*(v2BXTZ(i,j+1,k)-v2BXTZ(i,j,k))+rnn*(v2BXTZ(i,j+2,k)-        & 
    v2BXTZ(i,j-1,k))) 
   s33zBXTZ(i,j,k)=(1-dafz)/c0z*s33zBXTZ(i,j,k)+dtx/c0z*       &
   (lamBXTZ(i,j,k)+2*muBXTZ(i,j,k))*(v3BXTZ(i,j,k)-v3BXTZ(i,j,k-1))
   s33BXTZ(i,j,k)=s33xBXTZ(i,j,k)+s33yBXTZ(i,j,k)+s33zBXTZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm
 do j=3,ny-2
  do i=1,npm-1
   dafx=DiBXS1(i,j,nz-2)
   c0x=1.+dafx
   mu12=4./(1./muBXTZ(i,j,k)+1./muBXTZ(i,j-1,k)+1./muBXTZ(i+1,j,k)+ &
    1./muBXTZ(i+1,j-1,k))
  s12xBXTZ(i,j,k)=(1-dafx)*s12xBXTZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2BXTZ(i+1,j,k)-v2BXTZ(i,j,k))
   s12yBXTZ(i,j,k)=s12yBXTZ(i,j,k)+dtx*mu12*  &
    (rn*(v1BXTZ(i,j,k)-v1BXTZ(i,j-1,k))+rnn*(v1BXTZ(i,j+1,k)-       &
    v1BXTZ(i,j-2,k)))
   s12BXTZ(i,j,k)=s12xBXTZ(i,j,k)+s12yBXTZ(i,j,k)
  enddo
 enddo
enddo

do k=2,npm-1
 do j=3,ny-2
  do i=1,npm-1
   dafx=DiBXS1(i,j,nz-2)
   dafz=DiTZS1(3,j,k)
   c0x=1.+dafx
   c0z=1.+dafz
   mu13=4./(1./muBXTZ(i,j,k)+1./muBXTZ(i,j,k+1)+1./muBXTZ(i+1,j,k)+ &
    1./muBXTZ(i+1,j,k+1))
   s13xBXTZ(i,j,k)=(1-dafx)*s13xBXTZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3BXTZ(i+1,j,k)-v3BXTZ(i,j,k))
   s13zBXTZ(i,j,k)=(1-dafz)*s13zBXTZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1BXTZ(i,j,k+1)-v1BXTZ(i,j,k))
   s13BXTZ(i,j,k)=s13xBXTZ(i,j,k)+s13zBXTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm-1
   dafz=DiTZS1(3,j,k)
   c0z=1.+dafz
   mu23=4./(1./muBXTZ(i,j,k)+1./muBXTZ(i,j,k+1)+1./muBXTZ(i,j-1,k)+ &
    1./muBXTZ(i,j-1,k+1))
   s23yBXTZ(i,j,k)=s23yBXTZ(i,j,k)+dtx*mu23*    &
    (rn*(v3BXTZ(i,j,k)-v3BXTZ(i,j-1,k))+rnn*(v3BXTZ(i,j+1,k)-       & 
    v3BXTZ(i,j-2,k)))
   s23zBXTZ(i,j,k)=(1-dafz)*s23zBXTZ(i,j,k)/c0z+dtx/c0z*mu23  &  
    *(v2BXTZ(i,j,k+1)-v2BXTZ(i,j,k))
   s23BXTZ(i,j,k)=s23yBXTZ(i,j,k)+s23zBXTZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! #######################################################
subroutine sPMLBXBZ
include 'triffy.dec'
do k=2,npm-1
 do j=3,ny-2
  do i=2,npm-1
   dafx=DiBXS0(i,j,3)
   dafz=DiBZS0(3,j,k)
   c0x=1+dafx
   c0z=1+dafz
   s11xBXBZ(i,j,k)=(1-dafx)/c0x*s11xBXBZ(i,j,k)+dtx/c0x*       & 
    (lamBXBZ(i,j,k)+2*muBXBZ(i,j,k))*(v1BXBZ(i,j,k)-v1BXBZ(i-1,j,k))
   s11yBXBZ(i,j,k)=s11yBXBZ(i,j,k)+dtx*lamBXBZ(i,j,k)*           &
    (rn*(v2BXBZ(i,j+1,k)-v2BXBZ(i,j,k))+rnn*(v2BXBZ(i,j+2,k)-        &
    v2BXBZ(i,j-1,k)))
   s11zBXBZ(i,j,k)=(1-dafz)/c0z*s11zBXBZ(i,j,k)+dtx/c0z*       & 
    lamBXBZ(i,j,k)*(v3BXBZ(i,j,k)-v3BXBZ(i,j,k-1))
   s11BXBZ(i,j,k)=s11xBXBZ(i,j,k)+s11yBXBZ(i,j,k)+s11zBXBZ(i,j,k)
!
   s22xBXBZ(i,j,k)=(1-dafx)/c0x*s22xBXBZ(i,j,k)+dtx/c0x*       & 
    lamBXBZ(i,j,k)*(v1BXBZ(i,j,k)-v1BXBZ(i-1,j,k))
   s22yBXBZ(i,j,k)=s22yBXBZ(i,j,k)+dtx*(lamBXBZ(i,j,k)+          &
    2.*muBXBZ(i,j,k))*                                           &
    (rn*(v2BXBZ(i,j+1,k)-v2BXBZ(i,j,k))+rnn*(v2BXBZ(i,j+2,k)-        & 
    v2BXBZ(i,j-1,k))) 
   s22zBXBZ(i,j,k)=(1-dafz)/c0z*s22zBXBZ(i,j,k)+dtx/c0z*       & 
    lamBXBZ(i,j,k)*(v3BXBZ(i,j,k)-v3BXBZ(i,j,k-1)) 
   s22BXBZ(i,j,k)=s22xBXBZ(i,j,k)+s22yBXBZ(i,j,k)+s22zBXBZ(i,j,k)
!
   s33xBXBZ(i,j,k)=(1-dafx)/c0x*s33xBXBZ(i,j,k)+dtx/c0x*       & 
    lamBXBZ(i,j,k)*(v1BXBZ(i,j,k)-v1BXBZ(i-1,j,k))
   s33yBXBZ(i,j,k)=s33yBXBZ(i,j,k)+dtx*lamBXBZ(i,j,k)*           &
    (rn*(v2BXBZ(i,j+1,k)-v2BXBZ(i,j,k))+rnn*(v2BXBZ(i,j+2,k)-        & 
    v2BXBZ(i,j-1,k))) 
   s33zBXBZ(i,j,k)=(1-dafz)/c0z*s33zBXBZ(i,j,k)+dtx/c0z*       &
   (lamBXBZ(i,j,k)+2*muBXBZ(i,j,k))*(v3BXBZ(i,j,k)-v3BXBZ(i,j,k-1))
   s33BXBZ(i,j,k)=s33xBXBZ(i,j,k)+s33yBXBZ(i,j,k)+s33zBXBZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=3,ny-2
  do i=1,npm-1
   dafx=DiBXS1(i,j,3)
   c0x=1.+dafx
   mu12=4./(1./muBXBZ(i,j,k)+1./muBXBZ(i,j-1,k)+1./muBXBZ(i+1,j,k)+ &
    1./muBXBZ(i+1,j-1,k))
  s12xBXBZ(i,j,k)=(1-dafx)*s12xBXBZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2BXBZ(i+1,j,k)-v2BXBZ(i,j,k))
   s12yBXBZ(i,j,k)=s12yBXBZ(i,j,k)+dtx*mu12*  &
    (rn*(v1BXBZ(i,j,k)-v1BXBZ(i,j-1,k))+rnn*(v1BXBZ(i,j+1,k)-       &
    v1BXBZ(i,j-2,k)))
   s12BXBZ(i,j,k)=s12xBXBZ(i,j,k)+s12yBXBZ(i,j,k)
  enddo
 enddo
enddo

do k=1,npm-1
 do j=3,ny-2
  do i=1,npm-1
   dafx=DiBXS1(i,j,3)
   dafz=DiBZS1(3,j,k)
   c0x=1.+dafx
   c0z=1.+dafz
   mu13=4./(1./muBXBZ(i,j,k)+1./muBXBZ(i,j,k+1)+1./muBXBZ(i+1,j,k)+ &
    1./muBXBZ(i+1,j,k+1))
   s13xBXBZ(i,j,k)=(1-dafx)*s13xBXBZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3BXBZ(i+1,j,k)-v3BXBZ(i,j,k))
   s13zBXBZ(i,j,k)=(1-dafz)*s13zBXBZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1BXBZ(i,j,k+1)-v1BXBZ(i,j,k))
   s13BXBZ(i,j,k)=s13xBXBZ(i,j,k)+s13zBXBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=3,ny-2
  do i=2,npm-1
   dafz=DiBZS1(3,j,k)
   c0z=1.+dafz
   mu23=4./(1./muBXBZ(i,j,k)+1./muBXBZ(i,j,k+1)+1./muBXBZ(i,j-1,k)+ &
    1./muBXBZ(i,j-1,k+1))
   s23yBXBZ(i,j,k)=s23yBXBZ(i,j,k)+dtx*mu23*    &
    (rn*(v3BXBZ(i,j,k)-v3BXBZ(i,j-1,k))+rnn*(v3BXBZ(i,j+1,k)-       & 
    v3BXBZ(i,j-2,k)))
   s23zBXBZ(i,j,k)=(1-dafz)*s23zBXBZ(i,j,k)/c0z+dtx/c0z*mu23  &  
    *(v2BXBZ(i,j,k+1)-v2BXBZ(i,j,k))
   s23BXBZ(i,j,k)=s23yBXBZ(i,j,k)+s23zBXBZ(i,j,k)
!
  enddo
 enddo
enddo
return
end subroutine
! #####################################################
subroutine sPMLTYTZ
include 'triffy.dec'
do k=2,npm
 do j=2,npm-1
  do i=3,nx-2
   dafy=DiTYS0(i,j,nz-2)
   c0y=1+dafy
   dafz=DiTZS0(i,ny-2,k)
   c0z=1+dafz
   s11xTYTZ(i,j,k)=s11xTYTZ(i,j,k)+dtx*(lamTYTZ(i,j,k)+          &
    2.*muTYTZ(i,j,k))*(rn*(v1TYTZ(i,j,k)-v1TYTZ(i-1,j,k))+           &
    rnn*(v1TYTZ(i+1,j,k)-v1TYTZ(i-2,j,k)))
   s11yTYTZ(i,j,k)=(1-dafy)/c0y*s11yTYTZ(i,j,k)+dtx/c0y*       & 
    lamTYTZ(i,j,k)*(v2TYTZ(i,j+1,k)-v2TYTZ(i,j,k))
   s11zTYTZ(i,j,k)=(1-dafz)/c0z*s11zTYTZ(i,j,k)+dtx/c0z*       & 
    lamTYTZ(i,j,k)*(v3TYTZ(i,j,k)-v3TYTZ(i,j,k-1))
   s11TYTZ(i,j,k)=s11xTYTZ(i,j,k)+s11yTYTZ(i,j,k)+s11zTYTZ(i,j,k)
!
   s22xTYTZ(i,j,k)=s22xTYTZ(i,j,k)+dtx*lamTYTZ(i,j,k)*           &
    (rn*(v1TYTZ(i,j,k)-v1TYTZ(i-1,j,k))+rnn*                       &   
    (v1TYTZ(i+1,j,k)-v1TYTZ(i-2,j,k)))
   s22yTYTZ(i,j,k)=(1-dafy)/c0y*s22yTYTZ(i,j,k)+dtx/c0y*       & 
    (lamTYTZ(i,j,k)+2.*muTYTZ(i,j,k))*(v2TYTZ(i,j+1,k)-v2TYTZ(i,j,k))
   s22zTYTZ(i,j,k)=(1-dafz)/c0z*s22zTYTZ(i,j,k)+dtx/c0z*       & 
    lamTYTZ(i,j,k)*(v3TYTZ(i,j,k)-v3TYTZ(i,j,k-1))
   s22TYTZ(i,j,k)=s22xTYTZ(i,j,k)+s22yTYTZ(i,j,k)+s22zTYTZ(i,j,k)
!
   s33xTYTZ(i,j,k)=s33xTYTZ(i,j,k)+dtx*lamTYTZ(i,j,k)*           &
    (rn*(v1TYTZ(i,j,k)-v1TYTZ(i-1,j,k))+rnn*                       & 
    (v1TYTZ(i+1,j,k)-v1TYTZ(i-2,j,k)))
   s33yTYTZ(i,j,k)=(1-dafy)/c0y*s33yTYTZ(i,j,k)+dtx/c0y*       & 
    lamTYTZ(i,j,k)*(v2TYTZ(i,j+1,k)-v2TYTZ(i,j,k))
   s33zTYTZ(i,j,k)=(1-dafz)/c0z*s33zTYTZ(i,j,k)+dtx/c0z*       & 
   (lamTYTZ(i,j,k)+2*muTYTZ(i,j,k))*(v3TYTZ(i,j,k)-v3TYTZ(i,j,k-1))
   s33TYTZ(i,j,k)=s33xTYTZ(i,j,k)+s33yTYTZ(i,j,k)+s33zTYTZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm
 do j=2,npm
  do i=3,nx-2
   dafy=DiTYS1(i,j,nz-2)
   c0y=1+dafy
   mu12=4./(1./muTYTZ(i,j,k)+1./muTYTZ(i,j-1,k)+1./muTYTZ(i+1,j,k)+ &
    1./muTYTZ(i+1,j-1,k))
   s12xTYTZ(i,j,k)=s12xTYTZ(i,j,k)+dtx*mu12*                  &
    (rn*(v2TYTZ(i+1,j,k)-v2TYTZ(i,j,k))+rnn*                      &
    (v2TYTZ(i+2,j,k)-v2TYTZ(i-1,j,k)))
   s12yTYTZ(i,j,k)=(1-dafy)*s12yTYTZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1TYTZ(i,j,k)-v1TYTZ(i,j-1,k))
   s12TYTZ(i,j,k)=s12xTYTZ(i,j,k)+s12yTYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=3,nx-2
   dafz=DiTZS1(i,ny-2,k)
   c0z=1.+dafz
   mu13=4./(1./muTYTZ(i,j,k)+1./muTYTZ(i,j,k+1)+1./muTYTZ(i+1,j,k)+ &
    1./muTYTZ(i+1,j,k+1))
   s13xTYTZ(i,j,k)=s13xTYTZ(i,j,k)+dtx*mu13*                  &
    (rn*(v3TYTZ(i+1,j,k)-v3TYTZ(i,j,k))+rnn*                      &
    (v3TYTZ(i+2,j,k)-v3TYTZ(i-1,j,k)))
   s13zTYTZ(i,j,k)=(1-dafz)*s13zTYTZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1TYTZ(i,j,k+1)-v1TYTZ(i,j,k))
   s13TYTZ(i,j,k)=s13xTYTZ(i,j,k)+s13zTYTZ(i,j,k)
  enddo
 enddo
enddo
do k=2,npm-1
 do j=2,npm
  do i=3,nx-2
   dafy=DiTYS1(i,j,nz-2)
   c0y=1+dafy
   dafz=DiTZS1(i,ny-2,k)
   c0z=1.+dafz
   mu23=4./(1./muTYTZ(i,j,k)+1./muTYTZ(i,j,k+1)+1./muTYTZ(i,j-1,k)+ &
    1./muTYTZ(i,j-1,k+1))
   s23yTYTZ(i,j,k)=(1-dafy)*s23yTYTZ(i,j,k)/c0y+dtx/c0y*mu23  &
    *(v3TYTZ(i,j,k)-v3TYTZ(i,j-1,k))
   s23zTYTZ(i,j,k)=(1-dafz)*s23zTYTZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2TYTZ(i,j,k+1)-v2TYTZ(i,j,k))
   s23TYTZ(i,j,k)=s23yTYTZ(i,j,k)+s23zTYTZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine 
! #####################################################
subroutine sPMLTYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=3,nx-2
   dafy=DiTYS0(i,j,3)
   c0y=1+dafy
   dafz=DiBZS0(i,ny-2,k)
   c0z=1+dafz
   s11xTYBZ(i,j,k)=s11xTYBZ(i,j,k)+dtx*(lamTYBZ(i,j,k)+          &
    2.*muTYBZ(i,j,k))*(rn*(v1TYBZ(i,j,k)-v1TYBZ(i-1,j,k))+           &
    rnn*(v1TYBZ(i+1,j,k)-v1TYBZ(i-2,j,k)))
   s11yTYBZ(i,j,k)=(1-dafy)/c0y*s11yTYBZ(i,j,k)+dtx/c0y*       & 
    lamTYBZ(i,j,k)*(v2TYBZ(i,j+1,k)-v2TYBZ(i,j,k))
   s11zTYBZ(i,j,k)=(1-dafz)/c0z*s11zTYBZ(i,j,k)+dtx/c0z*       & 
    lamTYBZ(i,j,k)*(v3TYBZ(i,j,k)-v3TYBZ(i,j,k-1))
   s11TYBZ(i,j,k)=s11xTYBZ(i,j,k)+s11yTYBZ(i,j,k)+s11zTYBZ(i,j,k)
!
   s22xTYBZ(i,j,k)=s22xTYBZ(i,j,k)+dtx*lamTYBZ(i,j,k)*           &
    (rn*(v1TYBZ(i,j,k)-v1TYBZ(i-1,j,k))+rnn*                       &   
    (v1TYBZ(i+1,j,k)-v1TYBZ(i-2,j,k)))
   s22yTYBZ(i,j,k)=(1-dafy)/c0y*s22yTYBZ(i,j,k)+dtx/c0y*       & 
    (lamTYBZ(i,j,k)+2.*muTYBZ(i,j,k))*(v2TYBZ(i,j+1,k)-v2TYBZ(i,j,k))
   s22zTYBZ(i,j,k)=(1-dafz)/c0z*s22zTYBZ(i,j,k)+dtx/c0z*       & 
    lamTYBZ(i,j,k)*(v3TYBZ(i,j,k)-v3TYBZ(i,j,k-1))
   s22TYBZ(i,j,k)=s22xTYBZ(i,j,k)+s22yTYBZ(i,j,k)+s22zTYBZ(i,j,k)
!
   s33xTYBZ(i,j,k)=s33xTYBZ(i,j,k)+dtx*lamTYBZ(i,j,k)*           &
    (rn*(v1TYBZ(i,j,k)-v1TYBZ(i-1,j,k))+rnn*                       & 
    (v1TYBZ(i+1,j,k)-v1TYBZ(i-2,j,k)))
   s33yTYBZ(i,j,k)=(1-dafy)/c0y*s33yTYBZ(i,j,k)+dtx/c0y*       & 
    lamTYBZ(i,j,k)*(v2TYBZ(i,j+1,k)-v2TYBZ(i,j,k))
   s33zTYBZ(i,j,k)=(1-dafz)/c0z*s33zTYBZ(i,j,k)+dtx/c0z*       & 
   (lamTYBZ(i,j,k)+2*muTYBZ(i,j,k))*(v3TYBZ(i,j,k)-v3TYBZ(i,j,k-1))
   s33TYBZ(i,j,k)=s33xTYBZ(i,j,k)+s33yTYBZ(i,j,k)+s33zTYBZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm
  do i=3,nx-2
   dafy=DiTYS1(i,j,3)
   c0y=1+dafy
   mu12=4./(1./muTYBZ(i,j,k)+1./muTYBZ(i,j-1,k)+1./muTYBZ(i+1,j,k)+ &
    1./muTYBZ(i+1,j-1,k))
   s12xTYBZ(i,j,k)=s12xTYBZ(i,j,k)+dtx*mu12*                  &
    (rn*(v2TYBZ(i+1,j,k)-v2TYBZ(i,j,k))+rnn*                      &
    (v2TYBZ(i+2,j,k)-v2TYBZ(i-1,j,k)))
   s12yTYBZ(i,j,k)=(1-dafy)*s12yTYBZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1TYBZ(i,j,k)-v1TYBZ(i,j-1,k))
   s12TYBZ(i,j,k)=s12xTYBZ(i,j,k)+s12yTYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=2,npm-1
  do i=3,nx-2
   dafz=DiBZS1(i,ny-2,k)
   c0z=1.+dafz
   mu13=4./(1./muTYBZ(i,j,k)+1./muTYBZ(i,j,k+1)+1./muTYBZ(i+1,j,k)+ &
    1./muTYBZ(i+1,j,k+1))
   s13xTYBZ(i,j,k)=s13xTYBZ(i,j,k)+dtx*mu13*                  &
    (rn*(v3TYBZ(i+1,j,k)-v3TYBZ(i,j,k))+rnn*                      &
    (v3TYBZ(i+2,j,k)-v3TYBZ(i-1,j,k)))
   s13zTYBZ(i,j,k)=(1-dafz)*s13zTYBZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1TYBZ(i,j,k+1)-v1TYBZ(i,j,k))
   s13TYBZ(i,j,k)=s13xTYBZ(i,j,k)+s13zTYBZ(i,j,k)
  enddo
 enddo
enddo
do k=1,npm-1
 do j=2,npm
  do i=3,nx-2
   dafy=DiTYS1(i,j,3)
   c0y=1+dafy
   dafz=DiBZS1(i,ny-2,k)
   c0z=1.+dafz
   mu23=4./(1./muTYBZ(i,j,k)+1./muTYBZ(i,j,k+1)+1./muTYBZ(i,j-1,k)+ &
    1./muTYBZ(i,j-1,k+1))
   s23yTYBZ(i,j,k)=(1-dafy)*s23yTYBZ(i,j,k)/c0y+dtx/c0y*mu23  &
    *(v3TYBZ(i,j,k)-v3TYBZ(i,j-1,k))
   s23zTYBZ(i,j,k)=(1-dafz)*s23zTYBZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2TYBZ(i,j,k+1)-v2TYBZ(i,j,k))
   s23TYBZ(i,j,k)=s23yTYBZ(i,j,k)+s23zTYBZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine 
! #####################################################
subroutine sPMLBYTZ
include 'triffy.dec'
do k=2,npm
 do j=1,npm-1
  do i=3,nx-2
   dafy=DiBYS0(i,j,nz-2)
   c0y=1+dafy
   dafz=DiTZS0(i,3,k)
   c0z=1+dafz
   s11xBYTZ(i,j,k)=s11xBYTZ(i,j,k)+dtx*(lamBYTZ(i,j,k)+          &
    2.*muBYTZ(i,j,k))*(rn*(v1BYTZ(i,j,k)-v1BYTZ(i-1,j,k))+           &
    rnn*(v1BYTZ(i+1,j,k)-v1BYTZ(i-2,j,k)))
   s11yBYTZ(i,j,k)=(1-dafy)/c0y*s11yBYTZ(i,j,k)+dtx/c0y*       & 
    lamBYTZ(i,j,k)*(v2BYTZ(i,j+1,k)-v2BYTZ(i,j,k))
   s11zBYTZ(i,j,k)=(1-dafz)/c0z*s11zBYTZ(i,j,k)+dtx/c0z*       & 
    lamBYTZ(i,j,k)*(v3BYTZ(i,j,k)-v3BYTZ(i,j,k-1))
   s11BYTZ(i,j,k)=s11xBYTZ(i,j,k)+s11yBYTZ(i,j,k)+s11zBYTZ(i,j,k)
!
   s22xBYTZ(i,j,k)=s22xBYTZ(i,j,k)+dtx*lamBYTZ(i,j,k)*           &
    (rn*(v1BYTZ(i,j,k)-v1BYTZ(i-1,j,k))+rnn*                       &   
    (v1BYTZ(i+1,j,k)-v1BYTZ(i-2,j,k)))
   s22yBYTZ(i,j,k)=(1-dafy)/c0y*s22yBYTZ(i,j,k)+dtx/c0y*       & 
    (lamBYTZ(i,j,k)+2.*muBYTZ(i,j,k))*(v2BYTZ(i,j+1,k)-v2BYTZ(i,j,k))
   s22zBYTZ(i,j,k)=(1-dafz)/c0z*s22zBYTZ(i,j,k)+dtx/c0z*       & 
    lamBYTZ(i,j,k)*(v3BYTZ(i,j,k)-v3BYTZ(i,j,k-1))
   s22BYTZ(i,j,k)=s22xBYTZ(i,j,k)+s22yBYTZ(i,j,k)+s22zBYTZ(i,j,k)
!
   s33xBYTZ(i,j,k)=s33xBYTZ(i,j,k)+dtx*lamBYTZ(i,j,k)*           &
    (rn*(v1BYTZ(i,j,k)-v1BYTZ(i-1,j,k))+rnn*                       & 
    (v1BYTZ(i+1,j,k)-v1BYTZ(i-2,j,k)))
   s33yBYTZ(i,j,k)=(1-dafy)/c0y*s33yBYTZ(i,j,k)+dtx/c0y*       & 
    lamBYTZ(i,j,k)*(v2BYTZ(i,j+1,k)-v2BYTZ(i,j,k))
   s33zBYTZ(i,j,k)=(1-dafz)/c0z*s33zBYTZ(i,j,k)+dtx/c0z*       & 
   (lamBYTZ(i,j,k)+2*muBYTZ(i,j,k))*(v3BYTZ(i,j,k)-v3BYTZ(i,j,k-1))
   s33BYTZ(i,j,k)=s33xBYTZ(i,j,k)+s33yBYTZ(i,j,k)+s33zBYTZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm
 do j=2,npm-1
  do i=3,nx-2
   dafy=DiBYS1(i,j,nz-2)
   c0y=1+dafy
   mu12=4./(1./muBYTZ(i,j,k)+1./muBYTZ(i,j-1,k)+1./muBYTZ(i+1,j,k)+ &
    1./muBYTZ(i+1,j-1,k))
   s12xBYTZ(i,j,k)=s12xBYTZ(i,j,k)+dtx*mu12*                  &
    (rn*(v2BYTZ(i+1,j,k)-v2BYTZ(i,j,k))+rnn*                      &
    (v2BYTZ(i+2,j,k)-v2BYTZ(i-1,j,k)))
   s12yBYTZ(i,j,k)=(1-dafy)*s12yBYTZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1BYTZ(i,j,k)-v1BYTZ(i,j-1,k))
   s12BYTZ(i,j,k)=s12xBYTZ(i,j,k)+s12yBYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=1,npm-1
  do i=3,nx-2
   dafz=DiTZS1(i,3,k)
   c0z=1.+dafz
   mu13=4./(1./muBYTZ(i,j,k)+1./muBYTZ(i,j,k+1)+1./muBYTZ(i+1,j,k)+ &
    1./muBYTZ(i+1,j,k+1))
   s13xBYTZ(i,j,k)=s13xBYTZ(i,j,k)+dtx*mu13*                  &
    (rn*(v3BYTZ(i+1,j,k)-v3BYTZ(i,j,k))+rnn*                      &
    (v3BYTZ(i+2,j,k)-v3BYTZ(i-1,j,k)))
   s13zBYTZ(i,j,k)=(1-dafz)*s13zBYTZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1BYTZ(i,j,k+1)-v1BYTZ(i,j,k))
   s13BYTZ(i,j,k)=s13xBYTZ(i,j,k)+s13zBYTZ(i,j,k)
  enddo
 enddo
enddo
do k=2,npm-1
 do j=2,npm-1
  do i=3,nx-2
   dafy=DiBYS1(i,j,nz-2)
   c0y=1+dafy
   dafz=DiTZS1(i,3,k)
   c0z=1.+dafz
   mu23=4./(1./muBYTZ(i,j,k)+1./muBYTZ(i,j,k+1)+1./muBYTZ(i,j-1,k)+ &
    1./muBYTZ(i,j-1,k+1))
   s23yBYTZ(i,j,k)=(1-dafy)*s23yBYTZ(i,j,k)/c0y+dtx/c0y*mu23  &
    *(v3BYTZ(i,j,k)-v3BYTZ(i,j-1,k))
   s23zBYTZ(i,j,k)=(1-dafz)*s23zBYTZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2BYTZ(i,j,k+1)-v2BYTZ(i,j,k))
   s23BYTZ(i,j,k)=s23yBYTZ(i,j,k)+s23zBYTZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine 
! #####################################################
subroutine sPMLBYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=1,npm-1
  do i=3,nx-2
   dafy=DiBYS0(i,j,3)
   c0y=1+dafy
   dafz=DiBZS0(i,3,k)
   c0z=1+dafz
   s11xBYBZ(i,j,k)=s11xBYBZ(i,j,k)+dtx*(lamBYBZ(i,j,k)+          &
    2.*muBYBZ(i,j,k))*(rn*(v1BYBZ(i,j,k)-v1BYBZ(i-1,j,k))+           &
    rnn*(v1BYBZ(i+1,j,k)-v1BYBZ(i-2,j,k)))
   s11yBYBZ(i,j,k)=(1-dafy)/c0y*s11yBYBZ(i,j,k)+dtx/c0y*       & 
    lamBYBZ(i,j,k)*(v2BYBZ(i,j+1,k)-v2BYBZ(i,j,k))
   s11zBYBZ(i,j,k)=(1-dafz)/c0z*s11zBYBZ(i,j,k)+dtx/c0z*       & 
    lamBYBZ(i,j,k)*(v3BYBZ(i,j,k)-v3BYBZ(i,j,k-1))
   s11BYBZ(i,j,k)=s11xBYBZ(i,j,k)+s11yBYBZ(i,j,k)+s11zBYBZ(i,j,k)
!
   s22xBYBZ(i,j,k)=s22xBYBZ(i,j,k)+dtx*lamBYBZ(i,j,k)*           &
    (rn*(v1BYBZ(i,j,k)-v1BYBZ(i-1,j,k))+rnn*                       &   
    (v1BYBZ(i+1,j,k)-v1BYBZ(i-2,j,k)))
   s22yBYBZ(i,j,k)=(1-dafy)/c0y*s22yBYBZ(i,j,k)+dtx/c0y*       & 
    (lamBYBZ(i,j,k)+2.*muBYBZ(i,j,k))*(v2BYBZ(i,j+1,k)-v2BYBZ(i,j,k))
   s22zBYBZ(i,j,k)=(1-dafz)/c0z*s22zBYBZ(i,j,k)+dtx/c0z*       & 
    lamBYBZ(i,j,k)*(v3BYBZ(i,j,k)-v3BYBZ(i,j,k-1))
   s22BYBZ(i,j,k)=s22xBYBZ(i,j,k)+s22yBYBZ(i,j,k)+s22zBYBZ(i,j,k)
!
   s33xBYBZ(i,j,k)=s33xBYBZ(i,j,k)+dtx*lamBYBZ(i,j,k)*           &
    (rn*(v1BYBZ(i,j,k)-v1BYBZ(i-1,j,k))+rnn*                       & 
    (v1BYBZ(i+1,j,k)-v1BYBZ(i-2,j,k)))
   s33yBYBZ(i,j,k)=(1-dafy)/c0y*s33yBYBZ(i,j,k)+dtx/c0y*       & 
    lamBYBZ(i,j,k)*(v2BYBZ(i,j+1,k)-v2BYBZ(i,j,k))
   s33zBYBZ(i,j,k)=(1-dafz)/c0z*s33zBYBZ(i,j,k)+dtx/c0z*       & 
   (lamBYBZ(i,j,k)+2*muBYBZ(i,j,k))*(v3BYBZ(i,j,k)-v3BYBZ(i,j,k-1))
   s33BYBZ(i,j,k)=s33xBYBZ(i,j,k)+s33yBYBZ(i,j,k)+s33zBYBZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=3,nx-2
   dafy=DiBYS1(i,j,3)
   c0y=1+dafy
   mu12=4./(1./muBYBZ(i,j,k)+1./muBYBZ(i,j-1,k)+1./muBYBZ(i+1,j,k)+ &
    1./muBYBZ(i+1,j-1,k))
   s12xBYBZ(i,j,k)=s12xBYBZ(i,j,k)+dtx*mu12*                  &
    (rn*(v2BYBZ(i+1,j,k)-v2BYBZ(i,j,k))+rnn*                      &
    (v2BYBZ(i+2,j,k)-v2BYBZ(i-1,j,k)))
   s12yBYBZ(i,j,k)=(1-dafy)*s12yBYBZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1BYBZ(i,j,k)-v1BYBZ(i,j-1,k))
   s12BYBZ(i,j,k)=s12xBYBZ(i,j,k)+s12yBYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=1,npm-1
  do i=3,nx-2
   dafz=DiBZS1(i,3,k)
   c0z=1.+dafz
   mu13=4./(1./muBYBZ(i,j,k)+1./muBYBZ(i,j,k+1)+1./muBYBZ(i+1,j,k)+ &
    1./muBYBZ(i+1,j,k+1))
   s13xBYBZ(i,j,k)=s13xBYBZ(i,j,k)+dtx*mu13*                  &
    (rn*(v3BYBZ(i+1,j,k)-v3BYBZ(i,j,k))+rnn*                      &
    (v3BYBZ(i+2,j,k)-v3BYBZ(i-1,j,k)))
   s13zBYBZ(i,j,k)=(1-dafz)*s13zBYBZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1BYBZ(i,j,k+1)-v1BYBZ(i,j,k))
   s13BYBZ(i,j,k)=s13xBYBZ(i,j,k)+s13zBYBZ(i,j,k)
  enddo
 enddo
enddo
do k=1,npm-1
 do j=2,npm-1
  do i=3,nx-2
   dafy=DiBYS1(i,j,3)
   c0y=1+dafy
   dafz=DiBZS1(i,3,k)
   c0z=1.+dafz
   mu23=4./(1./muBYBZ(i,j,k)+1./muBYBZ(i,j,k+1)+1./muBYBZ(i,j-1,k)+ &
    1./muBYBZ(i,j-1,k+1))
   s23yBYBZ(i,j,k)=(1-dafy)*s23yBYBZ(i,j,k)/c0y+dtx/c0y*mu23  &
    *(v3BYBZ(i,j,k)-v3BYBZ(i,j-1,k))
   s23zBYBZ(i,j,k)=(1-dafz)*s23zBYBZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2BYBZ(i,j,k+1)-v2BYBZ(i,j,k))
   s23BYBZ(i,j,k)=s23yBYBZ(i,j,k)+s23zBYBZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine 
! ######################################
Subroutine sPMLTXTYTZ
include 'triffy.dec'
do k=2,npm
 do j=2,npm-1
  do i=2,npm
   dafx=DiTXS0(i,ny-2,nz-2)
   c0x=1+dafx
   dafy=DiTYS0(nx-2,j,nz-2)
   c0y=1+dafy
   dafz=DiTZS0(nx-2,ny-2,k)
   c0z=1+dafz
   s11xTXTYTZ(i,j,k)=(1-dafx)/c0x*s11xTXTYTZ(i,j,k)+dtx/c0x*       & 
    (lamTXTYTZ(i,j,k)+2*muTXTYTZ(i,j,k))*(v1TXTYTZ(i,j,k)-v1TXTYTZ(i-1,j,k))
   s11yTXTYTZ(i,j,k)=(1-dafy)/c0y*s11yTXTYTZ(i,j,k)+dtx/c0y*       & 
    lamTXTYTZ(i,j,k)*(v2TXTYTZ(i,j+1,k)-v2TXTYTZ(i,j,k))
   s11zTXTYTZ(i,j,k)=(1-dafz)/c0z*s11zTXTYTZ(i,j,k)+dtx/c0z*       & 
    lamTXTYTZ(i,j,k)*(v3TXTYTZ(i,j,k)-v3TXTYTZ(i,j,k-1))
   s11TXTYTZ(i,j,k)=s11xTXTYTZ(i,j,k)+s11yTXTYTZ(i,j,k)+s11zTXTYTZ(i,j,k)
! 
   s22xTXTYTZ(i,j,k)=(1-dafx)/c0x*s22xTXTYTZ(i,j,k)+dtx/c0x*       & 
    lamTXTYTZ(i,j,k)*(v1TXTYTZ(i,j,k)-v1TXTYTZ(i-1,j,k))
   s22yTXTYTZ(i,j,k)=(1-dafy)/c0y*s22yTXTYTZ(i,j,k)+dtx/c0y*       & 
    (lamTXTYTZ(i,j,k)+2.*muTXTYTZ(i,j,k))*(v2TXTYTZ(i,j+1,k)-v2TXTYTZ(i,j,k))
   s22zTXTYTZ(i,j,k)=(1-dafz)/c0z*s22zTXTYTZ(i,j,k)+dtx/c0z*       & 
    lamTXTYTZ(i,j,k)*(v3TXTYTZ(i,j,k)-v3TXTYTZ(i,j,k-1))
   s22TXTYTZ(i,j,k)=s22xTXTYTZ(i,j,k)+s22yTXTYTZ(i,j,k)+s22zTXTYTZ(i,j,k)
!
   s33xTXTYTZ(i,j,k)=(1-dafx)/c0x*s33xTXTYTZ(i,j,k)+dtx/c0x*       & 
    lamTXTYTZ(i,j,k)*(v1TXTYTZ(i,j,k)-v1TXTYTZ(i-1,j,k))
   s33yTXTYTZ(i,j,k)=(1-dafy)/c0y*s33yTXTYTZ(i,j,k)+dtx/c0y*       & 
    lamTXTYTZ(i,j,k)*(v2TXTYTZ(i,j+1,k)-v2TXTYTZ(i,j,k))
   s33zTXTYTZ(i,j,k)=(1-dafz)/c0z*s33zTXTYTZ(i,j,k)+dtx/c0z*       & 
   (lamTXTYTZ(i,j,k)+2*muTXTYTZ(i,j,k))*(v3TXTYTZ(i,j,k)-v3TXTYTZ(i,j,k-1))
   s33TXTYTZ(i,j,k)=s33xTXTYTZ(i,j,k)+s33yTXTYTZ(i,j,k)+s33zTXTYTZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm
 do j=2,npm
  do i=2,npm-1
   dafx=DiTXS1(i,ny-2,nz-2)
   c0x=1.+dafx
   dafy=DiTYS1(nx-2,j,nz-2)
   c0y=1.+dafy
   mu12=4./(1./muTXTYTZ(i,j,k)+1./muTXTYTZ(i,j-1,k)+1./muTXTYTZ(i+1,j,k)+ &
    1./muTXTYTZ(i+1,j-1,k))
   s12xTXTYTZ(i,j,k)=(1-dafx)*s12xTXTYTZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2TXTYTZ(i+1,j,k)-v2TXTYTZ(i,j,k))
   s12yTXTYTZ(i,j,k)=(1-dafy)*s12yTXTYTZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1TXTYTZ(i,j,k)-v1TXTYTZ(i,j-1,k))
   s12TXTYTZ(i,j,k)=s12xTXTYTZ(i,j,k)+s12yTXTYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,ny-2,nz-2)
   c0x=1.+dafx
   dafz=DiTZS1(nx-2,ny-2,k)
   c0z=1.+dafz
   mu13=4./(1./muTXTYTZ(i,j,k)+1./muTXTYTZ(i,j,k+1)+1./muTXTYTZ(i+1,j,k)+ &
    1./muTXTYTZ(i+1,j,k+1))
   s13xTXTYTZ(i,j,k)=(1-dafx)*s13xTXTYTZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3TXTYTZ(i+1,j,k)-v3TXTYTZ(i,j,k))
   s13zTXTYTZ(i,j,k)=(1-dafz)*s13zTXTYTZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1TXTYTZ(i,j,k+1)-v1TXTYTZ(i,j,k))
   s13TXTYTZ(i,j,k)=s13xTXTYTZ(i,j,k)+s13zTXTYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm
  do i=2,npm
   dafy=DiTYS1(nx-2,j,nz-2)
   c0y=1.+dafy
   dafz=DiTZS1(nx-2,ny-2,k)
   c0z=1.+dafz
   mu23=4./(1./muTXTYTZ(i,j,k)+1./muTXTYTZ(i,j,k+1)+1./muTXTYTZ(i,j-1,k)+ &
    1./muTXTYTZ(i,j-1,k+1))
   s23yTXTYTZ(i,j,k)=(1-dafy)*s23yTXTYTZ(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3TXTYTZ(i,j,k)-v3TXTYTZ(i,j-1,k))
   s23zTXTYTZ(i,j,k)=(1-dafz)*s23zTXTYTZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2TXTYTZ(i,j,k+1)-v2TXTYTZ(i,j,k))
   s23TXTYTZ(i,j,k)=s23yTXTYTZ(i,j,k)+s23zTXTYTZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine sPMLTXTYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm
   dafx=DiTXS0(i,ny-2,3)
   c0x=1+dafx
   dafy=DiTYS0(nx-2,j,3)
   c0y=1+dafy
   dafz=DiBZS0(nx-2,ny-2,k)
   c0z=1+dafz
   s11xTXTYBZ(i,j,k)=(1-dafx)/c0x*s11xTXTYBZ(i,j,k)+dtx/c0x*       & 
    (lamTXTYBZ(i,j,k)+2*muTXTYBZ(i,j,k))*(v1TXTYBZ(i,j,k)-v1TXTYBZ(i-1,j,k))
   s11yTXTYBZ(i,j,k)=(1-dafy)/c0y*s11yTXTYBZ(i,j,k)+dtx/c0y*       & 
    lamTXTYBZ(i,j,k)*(v2TXTYBZ(i,j+1,k)-v2TXTYBZ(i,j,k))
   s11zTXTYBZ(i,j,k)=(1-dafz)/c0z*s11zTXTYBZ(i,j,k)+dtx/c0z*       & 
    lamTXTYBZ(i,j,k)*(v3TXTYBZ(i,j,k)-v3TXTYBZ(i,j,k-1))
   s11TXTYBZ(i,j,k)=s11xTXTYBZ(i,j,k)+s11yTXTYBZ(i,j,k)+s11zTXTYBZ(i,j,k)
! 
   s22xTXTYBZ(i,j,k)=(1-dafx)/c0x*s22xTXTYBZ(i,j,k)+dtx/c0x*       & 
    lamTXTYBZ(i,j,k)*(v1TXTYBZ(i,j,k)-v1TXTYBZ(i-1,j,k))
   s22yTXTYBZ(i,j,k)=(1-dafy)/c0y*s22yTXTYBZ(i,j,k)+dtx/c0y*       & 
    (lamTXTYBZ(i,j,k)+2.*muTXTYBZ(i,j,k))*(v2TXTYBZ(i,j+1,k)-v2TXTYBZ(i,j,k))
   s22zTXTYBZ(i,j,k)=(1-dafz)/c0z*s22zTXTYBZ(i,j,k)+dtx/c0z*       & 
    lamTXTYBZ(i,j,k)*(v3TXTYBZ(i,j,k)-v3TXTYBZ(i,j,k-1))
   s22TXTYBZ(i,j,k)=s22xTXTYBZ(i,j,k)+s22yTXTYBZ(i,j,k)+s22zTXTYBZ(i,j,k)
!
   s33xTXTYBZ(i,j,k)=(1-dafx)/c0x*s33xTXTYBZ(i,j,k)+dtx/c0x*       & 
    lamTXTYBZ(i,j,k)*(v1TXTYBZ(i,j,k)-v1TXTYBZ(i-1,j,k))
   s33yTXTYBZ(i,j,k)=(1-dafy)/c0y*s33yTXTYBZ(i,j,k)+dtx/c0y*       & 
    lamTXTYBZ(i,j,k)*(v2TXTYBZ(i,j+1,k)-v2TXTYBZ(i,j,k))
   s33zTXTYBZ(i,j,k)=(1-dafz)/c0z*s33zTXTYBZ(i,j,k)+dtx/c0z*       & 
   (lamTXTYBZ(i,j,k)+2*muTXTYBZ(i,j,k))*(v3TXTYBZ(i,j,k)-v3TXTYBZ(i,j,k-1))
   s33TXTYBZ(i,j,k)=s33xTXTYBZ(i,j,k)+s33yTXTYBZ(i,j,k)+s33zTXTYBZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm
  do i=2,npm-1
   dafx=DiTXS1(i,ny-2,3)
   c0x=1.+dafx
   dafy=DiTYS1(nx-2,j,3)
   c0y=1.+dafy
   mu12=4./(1./muTXTYBZ(i,j,k)+1./muTXTYBZ(i,j-1,k)+1./muTXTYBZ(i+1,j,k)+ &
    1./muTXTYBZ(i+1,j-1,k))
   s12xTXTYBZ(i,j,k)=(1-dafx)*s12xTXTYBZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2TXTYBZ(i+1,j,k)-v2TXTYBZ(i,j,k))
   s12yTXTYBZ(i,j,k)=(1-dafy)*s12yTXTYBZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1TXTYBZ(i,j,k)-v1TXTYBZ(i,j-1,k))
   s12TXTYBZ(i,j,k)=s12xTXTYBZ(i,j,k)+s12yTXTYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,ny-2,3)
   c0x=1.+dafx
   dafz=DiBZS1(nx-2,ny-2,k)
   c0z=1.+dafz
   mu13=4./(1./muTXTYBZ(i,j,k)+1./muTXTYBZ(i,j,k+1)+1./muTXTYBZ(i+1,j,k)+ &
    1./muTXTYBZ(i+1,j,k+1))
   s13xTXTYBZ(i,j,k)=(1-dafx)*s13xTXTYBZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3TXTYBZ(i+1,j,k)-v3TXTYBZ(i,j,k))
   s13zTXTYBZ(i,j,k)=(1-dafz)*s13zTXTYBZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1TXTYBZ(i,j,k+1)-v1TXTYBZ(i,j,k))
   s13TXTYBZ(i,j,k)=s13xTXTYBZ(i,j,k)+s13zTXTYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=2,npm
  do i=2,npm
   dafy=DiTYS1(nx-2,j,3)
   c0y=1.+dafy
   dafz=DiBZS1(nx-2,ny-2,k)
   c0z=1.+dafz
   mu23=4./(1./muTXTYBZ(i,j,k)+1./muTXTYBZ(i,j,k+1)+1./muTXTYBZ(i,j-1,k)+ &
    1./muTXTYBZ(i,j-1,k+1))
!    print *, s23yTXTYBZ(i,j,k),100
   s23yTXTYBZ(i,j,k)=(1-dafy)*s23yTXTYBZ(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3TXTYBZ(i,j,k)-v3TXTYBZ(i,j-1,k))
!    print *, s23yTXTYBZ(i,j,k),v3TXTYBZ(i,j,k),v3TXTYBZ(i,j-1,k),mu23,c0y
   s23zTXTYBZ(i,j,k)=(1-dafz)*s23zTXTYBZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2TXTYBZ(i,j,k+1)-v2TXTYBZ(i,j,k))
   s23TXTYBZ(i,j,k)=s23yTXTYBZ(i,j,k)+s23zTXTYBZ(i,j,k)
!   print *,s23yTXTYBZ(i,j,k),s23TXTYBZ(i,j,k),v3TXTYBZ(i,j,k),i,j,k,it
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine sPMLTXBYTZ
include 'triffy.dec'
do k=2,npm
 do j=1,npm-1
  do i=2,npm
   dafx=DiTXS0(i,3,nz-2)
   c0x=1+dafx
   dafy=DiBYS0(nx-2,j,nz-2)
   c0y=1+dafy
   dafz=DiTZS0(nx-2,3,k)
   c0z=1+dafz
   s11xTXBYTZ(i,j,k)=(1-dafx)/c0x*s11xTXBYTZ(i,j,k)+dtx/c0x*       & 
    (lamTXBYTZ(i,j,k)+2*muTXBYTZ(i,j,k))*(v1TXBYTZ(i,j,k)-v1TXBYTZ(i-1,j,k))
   s11yTXBYTZ(i,j,k)=(1-dafy)/c0y*s11yTXBYTZ(i,j,k)+dtx/c0y*       & 
    lamTXBYTZ(i,j,k)*(v2TXBYTZ(i,j+1,k)-v2TXBYTZ(i,j,k))
   s11zTXBYTZ(i,j,k)=(1-dafz)/c0z*s11zTXBYTZ(i,j,k)+dtx/c0z*       & 
    lamTXBYTZ(i,j,k)*(v3TXBYTZ(i,j,k)-v3TXBYTZ(i,j,k-1))
   s11TXBYTZ(i,j,k)=s11xTXBYTZ(i,j,k)+s11yTXBYTZ(i,j,k)+s11zTXBYTZ(i,j,k)
! 
   s22xTXBYTZ(i,j,k)=(1-dafx)/c0x*s22xTXBYTZ(i,j,k)+dtx/c0x*       & 
    lamTXBYTZ(i,j,k)*(v1TXBYTZ(i,j,k)-v1TXBYTZ(i-1,j,k))
   s22yTXBYTZ(i,j,k)=(1-dafy)/c0y*s22yTXBYTZ(i,j,k)+dtx/c0y*       & 
    (lamTXBYTZ(i,j,k)+2.*muTXBYTZ(i,j,k))*(v2TXBYTZ(i,j+1,k)-v2TXBYTZ(i,j,k))
   s22zTXBYTZ(i,j,k)=(1-dafz)/c0z*s22zTXBYTZ(i,j,k)+dtx/c0z*       & 
    lamTXBYTZ(i,j,k)*(v3TXBYTZ(i,j,k)-v3TXBYTZ(i,j,k-1))
   s22TXBYTZ(i,j,k)=s22xTXBYTZ(i,j,k)+s22yTXBYTZ(i,j,k)+s22zTXBYTZ(i,j,k)
!
   s33xTXBYTZ(i,j,k)=(1-dafx)/c0x*s33xTXBYTZ(i,j,k)+dtx/c0x*       & 
    lamTXBYTZ(i,j,k)*(v1TXBYTZ(i,j,k)-v1TXBYTZ(i-1,j,k))
   s33yTXBYTZ(i,j,k)=(1-dafy)/c0y*s33yTXBYTZ(i,j,k)+dtx/c0y*       & 
    lamTXBYTZ(i,j,k)*(v2TXBYTZ(i,j+1,k)-v2TXBYTZ(i,j,k))
   s33zTXBYTZ(i,j,k)=(1-dafz)/c0z*s33zTXBYTZ(i,j,k)+dtx/c0z*       & 
   (lamTXBYTZ(i,j,k)+2*muTXBYTZ(i,j,k))*(v3TXBYTZ(i,j,k)-v3TXBYTZ(i,j,k-1))
   s33TXBYTZ(i,j,k)=s33xTXBYTZ(i,j,k)+s33yTXBYTZ(i,j,k)+s33zTXBYTZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,3,nz-2)
   c0x=1.+dafx
   dafy=DiBYS1(nx-2,j,nz-2)
   c0y=1.+dafy
   mu12=4./(1./muTXBYTZ(i,j,k)+1./muTXBYTZ(i,j-1,k)+1./muTXBYTZ(i+1,j,k)+ &
    1./muTXBYTZ(i+1,j-1,k))
   s12xTXBYTZ(i,j,k)=(1-dafx)*s12xTXBYTZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2TXBYTZ(i+1,j,k)-v2TXBYTZ(i,j,k))
   s12yTXBYTZ(i,j,k)=(1-dafy)*s12yTXBYTZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1TXBYTZ(i,j,k)-v1TXBYTZ(i,j-1,k))
   s12TXBYTZ(i,j,k)=s12xTXBYTZ(i,j,k)+s12yTXBYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=1,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,3,nz-2)
   c0x=1.+dafx
   dafz=DiTZS1(nx-2,3,k)
   c0z=1.+dafz
   mu13=4./(1./muTXBYTZ(i,j,k)+1./muTXBYTZ(i,j,k+1)+1./muTXBYTZ(i+1,j,k)+ &
    1./muTXBYTZ(i+1,j,k+1))
   s13xTXBYTZ(i,j,k)=(1-dafx)*s13xTXBYTZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3TXBYTZ(i+1,j,k)-v3TXBYTZ(i,j,k))
   s13zTXBYTZ(i,j,k)=(1-dafz)*s13zTXBYTZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1TXBYTZ(i,j,k+1)-v1TXBYTZ(i,j,k))
   s13TXBYTZ(i,j,k)=s13xTXBYTZ(i,j,k)+s13zTXBYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm
   dafy=DiBYS1(nx-2,j,nz-2)
   c0y=1.+dafy
   dafz=DiTZS1(nx-2,3,k)
   c0z=1.+dafz
   mu23=4./(1./muTXBYTZ(i,j,k)+1./muTXBYTZ(i,j,k+1)+1./muTXBYTZ(i,j-1,k)+ &
    1./muTXBYTZ(i,j-1,k+1))
   s23yTXBYTZ(i,j,k)=(1-dafy)*s23yTXBYTZ(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3TXBYTZ(i,j,k)-v3TXBYTZ(i,j-1,k))
   s23zTXBYTZ(i,j,k)=(1-dafz)*s23zTXBYTZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2TXBYTZ(i,j,k+1)-v2TXBYTZ(i,j,k))
   s23TXBYTZ(i,j,k)=s23yTXBYTZ(i,j,k)+s23zTXBYTZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine sPMLTXBYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=1,npm-1
  do i=2,npm
   dafx=DiTXS0(i,3,3)
   c0x=1+dafx
   dafy=DiBYS0(nx-2,j,3)
   c0y=1+dafy
   dafz=DiBZS0(nx-2,3,k)
   c0z=1+dafz
   s11xTXBYBZ(i,j,k)=(1-dafx)/c0x*s11xTXBYBZ(i,j,k)+dtx/c0x*       & 
    (lamTXBYBZ(i,j,k)+2*muTXBYBZ(i,j,k))*(v1TXBYBZ(i,j,k)-v1TXBYBZ(i-1,j,k))
   s11yTXBYBZ(i,j,k)=(1-dafy)/c0y*s11yTXBYBZ(i,j,k)+dtx/c0y*       & 
    lamTXBYBZ(i,j,k)*(v2TXBYBZ(i,j+1,k)-v2TXBYBZ(i,j,k))
   s11zTXBYBZ(i,j,k)=(1-dafz)/c0z*s11zTXBYBZ(i,j,k)+dtx/c0z*       & 
    lamTXBYBZ(i,j,k)*(v3TXBYBZ(i,j,k)-v3TXBYBZ(i,j,k-1))
   s11TXBYBZ(i,j,k)=s11xTXBYBZ(i,j,k)+s11yTXBYBZ(i,j,k)+s11zTXBYBZ(i,j,k)
! 
   s22xTXBYBZ(i,j,k)=(1-dafx)/c0x*s22xTXBYBZ(i,j,k)+dtx/c0x*       & 
    lamTXBYBZ(i,j,k)*(v1TXBYBZ(i,j,k)-v1TXBYBZ(i-1,j,k))
   s22yTXBYBZ(i,j,k)=(1-dafy)/c0y*s22yTXBYBZ(i,j,k)+dtx/c0y*       & 
    (lamTXBYBZ(i,j,k)+2.*muTXBYBZ(i,j,k))*(v2TXBYBZ(i,j+1,k)-v2TXBYBZ(i,j,k))
   s22zTXBYBZ(i,j,k)=(1-dafz)/c0z*s22zTXBYBZ(i,j,k)+dtx/c0z*       & 
    lamTXBYBZ(i,j,k)*(v3TXBYBZ(i,j,k)-v3TXBYBZ(i,j,k-1))
   s22TXBYBZ(i,j,k)=s22xTXBYBZ(i,j,k)+s22yTXBYBZ(i,j,k)+s22zTXBYBZ(i,j,k)
!
   s33xTXBYBZ(i,j,k)=(1-dafx)/c0x*s33xTXBYBZ(i,j,k)+dtx/c0x*       & 
    lamTXBYBZ(i,j,k)*(v1TXBYBZ(i,j,k)-v1TXBYBZ(i-1,j,k))
   s33yTXBYBZ(i,j,k)=(1-dafy)/c0y*s33yTXBYBZ(i,j,k)+dtx/c0y*       & 
    lamTXBYBZ(i,j,k)*(v2TXBYBZ(i,j+1,k)-v2TXBYBZ(i,j,k))
   s33zTXBYBZ(i,j,k)=(1-dafz)/c0z*s33zTXBYBZ(i,j,k)+dtx/c0z*       & 
   (lamTXBYBZ(i,j,k)+2*muTXBYBZ(i,j,k))*(v3TXBYBZ(i,j,k)-v3TXBYBZ(i,j,k-1))
   s33TXBYBZ(i,j,k)=s33xTXBYBZ(i,j,k)+s33yTXBYBZ(i,j,k)+s33zTXBYBZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,3,3)
   c0x=1.+dafx
   dafy=DiBYS1(nx-2,j,3)
   c0y=1.+dafy
   mu12=4./(1./muTXBYBZ(i,j,k)+1./muTXBYBZ(i,j-1,k)+1./muTXBYBZ(i+1,j,k)+ &
    1./muTXBYBZ(i+1,j-1,k))
   s12xTXBYBZ(i,j,k)=(1-dafx)*s12xTXBYBZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2TXBYBZ(i+1,j,k)-v2TXBYBZ(i,j,k))
   s12yTXBYBZ(i,j,k)=(1-dafy)*s12yTXBYBZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1TXBYBZ(i,j,k)-v1TXBYBZ(i,j-1,k))
   s12TXBYBZ(i,j,k)=s12xTXBYBZ(i,j,k)+s12yTXBYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=1,npm-1
  do i=2,npm-1
   dafx=DiTXS1(i,3,3)
   c0x=1.+dafx
   dafz=DiBZS1(nx-2,3,k)
   c0z=1.+dafz
   mu13=4./(1./muTXBYBZ(i,j,k)+1./muTXBYBZ(i,j,k+1)+1./muTXBYBZ(i+1,j,k)+ &
    1./muTXBYBZ(i+1,j,k+1))
   s13xTXBYBZ(i,j,k)=(1-dafx)*s13xTXBYBZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3TXBYBZ(i+1,j,k)-v3TXBYBZ(i,j,k))
   s13zTXBYBZ(i,j,k)=(1-dafz)*s13zTXBYBZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1TXBYBZ(i,j,k+1)-v1TXBYBZ(i,j,k))
   s13TXBYBZ(i,j,k)=s13xTXBYBZ(i,j,k)+s13zTXBYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=2,npm-1
  do i=2,npm
   dafy=DiBYS1(nx-2,j,3)
   c0y=1.+dafy
   dafz=DiBZS1(nx-2,3,k)
   c0z=1.+dafz
   mu23=4./(1./muTXBYBZ(i,j,k)+1./muTXBYBZ(i,j,k+1)+1./muTXBYBZ(i,j-1,k)+ &
    1./muTXBYBZ(i,j-1,k+1))
   s23yTXBYBZ(i,j,k)=(1-dafy)*s23yTXBYBZ(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3TXBYBZ(i,j,k)-v3TXBYBZ(i,j-1,k))
   s23zTXBYBZ(i,j,k)=(1-dafz)*s23zTXBYBZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2TXBYBZ(i,j,k+1)-v2TXBYBZ(i,j,k))
   s23TXBYBZ(i,j,k)=s23yTXBYBZ(i,j,k)+s23zTXBYBZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine sPMLBXTYTZ
include 'triffy.dec'
do k=2,npm
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiBXS0(i,ny-2,nz-2)
   c0x=1+dafx
   dafy=DiTYS0(3,j,nz-2)
   c0y=1+dafy
   dafz=DiTZS0(3,ny-2,k)
   c0z=1+dafz
   s11xBXTYTZ(i,j,k)=(1-dafx)/c0x*s11xBXTYTZ(i,j,k)+dtx/c0x*       & 
    (lamBXTYTZ(i,j,k)+2*muBXTYTZ(i,j,k))*(v1BXTYTZ(i,j,k)-v1BXTYTZ(i-1,j,k))
   s11yBXTYTZ(i,j,k)=(1-dafy)/c0y*s11yBXTYTZ(i,j,k)+dtx/c0y*       & 
    lamBXTYTZ(i,j,k)*(v2BXTYTZ(i,j+1,k)-v2BXTYTZ(i,j,k))
   s11zBXTYTZ(i,j,k)=(1-dafz)/c0z*s11zBXTYTZ(i,j,k)+dtx/c0z*       & 
    lamBXTYTZ(i,j,k)*(v3BXTYTZ(i,j,k)-v3BXTYTZ(i,j,k-1))
   s11BXTYTZ(i,j,k)=s11xBXTYTZ(i,j,k)+s11yBXTYTZ(i,j,k)+s11zBXTYTZ(i,j,k)
! 
   s22xBXTYTZ(i,j,k)=(1-dafx)/c0x*s22xBXTYTZ(i,j,k)+dtx/c0x*       & 
    lamBXTYTZ(i,j,k)*(v1BXTYTZ(i,j,k)-v1BXTYTZ(i-1,j,k))
   s22yBXTYTZ(i,j,k)=(1-dafy)/c0y*s22yBXTYTZ(i,j,k)+dtx/c0y*       & 
    (lamBXTYTZ(i,j,k)+2.*muBXTYTZ(i,j,k))*(v2BXTYTZ(i,j+1,k)-v2BXTYTZ(i,j,k))
   s22zBXTYTZ(i,j,k)=(1-dafz)/c0z*s22zBXTYTZ(i,j,k)+dtx/c0z*       & 
    lamBXTYTZ(i,j,k)*(v3BXTYTZ(i,j,k)-v3BXTYTZ(i,j,k-1))
   s22BXTYTZ(i,j,k)=s22xBXTYTZ(i,j,k)+s22yBXTYTZ(i,j,k)+s22zBXTYTZ(i,j,k)
!
   s33xBXTYTZ(i,j,k)=(1-dafx)/c0x*s33xBXTYTZ(i,j,k)+dtx/c0x*       & 
    lamBXTYTZ(i,j,k)*(v1BXTYTZ(i,j,k)-v1BXTYTZ(i-1,j,k))
   s33yBXTYTZ(i,j,k)=(1-dafy)/c0y*s33yBXTYTZ(i,j,k)+dtx/c0y*       & 
    lamBXTYTZ(i,j,k)*(v2BXTYTZ(i,j+1,k)-v2BXTYTZ(i,j,k))
   s33zBXTYTZ(i,j,k)=(1-dafz)/c0z*s33zBXTYTZ(i,j,k)+dtx/c0z*       & 
   (lamBXTYTZ(i,j,k)+2*muBXTYTZ(i,j,k))*(v3BXTYTZ(i,j,k)-v3BXTYTZ(i,j,k-1))
   s33BXTYTZ(i,j,k)=s33xBXTYTZ(i,j,k)+s33yBXTYTZ(i,j,k)+s33zBXTYTZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm
 do j=2,npm
  do i=1,npm-1
   dafx=DiBXS1(i,ny-2,nz-2)
   c0x=1.+dafx
   dafy=DiTYS1(3,j,nz-2)
   c0y=1.+dafy
   mu12=4./(1./muBXTYTZ(i,j,k)+1./muBXTYTZ(i,j-1,k)+1./muBXTYTZ(i+1,j,k)+ &
    1./muBXTYTZ(i+1,j-1,k))
   s12xBXTYTZ(i,j,k)=(1-dafx)*s12xBXTYTZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2BXTYTZ(i+1,j,k)-v2BXTYTZ(i,j,k))
   s12yBXTYTZ(i,j,k)=(1-dafy)*s12yBXTYTZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1BXTYTZ(i,j,k)-v1BXTYTZ(i,j-1,k))
   s12BXTYTZ(i,j,k)=s12xBXTYTZ(i,j,k)+s12yBXTYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,ny-2,nz-2)
   c0x=1.+dafx
   dafz=DiTZS1(3,ny-2,k)
   c0z=1.+dafz
   mu13=4./(1./muBXTYTZ(i,j,k)+1./muBXTYTZ(i,j,k+1)+1./muBXTYTZ(i+1,j,k)+ &
    1./muBXTYTZ(i+1,j,k+1))
   s13xBXTYTZ(i,j,k)=(1-dafx)*s13xBXTYTZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3BXTYTZ(i+1,j,k)-v3BXTYTZ(i,j,k))
   s13zBXTYTZ(i,j,k)=(1-dafz)*s13zBXTYTZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1BXTYTZ(i,j,k+1)-v1BXTYTZ(i,j,k))
   s13BXTYTZ(i,j,k)=s13xBXTYTZ(i,j,k)+s13zBXTYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm
  do i=2,npm-1
   dafy=DiTYS1(3,j,nz-2)
   c0y=1.+dafy
   dafz=DiTZS1(3,ny-2,k)
   c0z=1.+dafz
   mu23=4./(1./muBXTYTZ(i,j,k)+1./muBXTYTZ(i,j,k+1)+1./muBXTYTZ(i,j-1,k)+ &
    1./muBXTYTZ(i,j-1,k+1))
   s23yBXTYTZ(i,j,k)=(1-dafy)*s23yBXTYTZ(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3BXTYTZ(i,j,k)-v3BXTYTZ(i,j-1,k))
   s23zBXTYTZ(i,j,k)=(1-dafz)*s23zBXTYTZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2BXTYTZ(i,j,k+1)-v2BXTYTZ(i,j,k))
   s23BXTYTZ(i,j,k)=s23yBXTYTZ(i,j,k)+s23zBXTYTZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine sPMLBXTYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
   dafx=DiBXS0(i,ny-2,3)
   c0x=1+dafx
   dafy=DiTYS0(3,j,3)
   c0y=1+dafy
   dafz=DiBZS0(3,ny-2,k)
   c0z=1+dafz
   s11xBXTYBZ(i,j,k)=(1-dafx)/c0x*s11xBXTYBZ(i,j,k)+dtx/c0x*       & 
    (lamBXTYBZ(i,j,k)+2*muBXTYBZ(i,j,k))*(v1BXTYBZ(i,j,k)-v1BXTYBZ(i-1,j,k))
   s11yBXTYBZ(i,j,k)=(1-dafy)/c0y*s11yBXTYBZ(i,j,k)+dtx/c0y*       & 
    lamBXTYBZ(i,j,k)*(v2BXTYBZ(i,j+1,k)-v2BXTYBZ(i,j,k))
   s11zBXTYBZ(i,j,k)=(1-dafz)/c0z*s11zBXTYBZ(i,j,k)+dtx/c0z*       & 
    lamBXTYBZ(i,j,k)*(v3BXTYBZ(i,j,k)-v3BXTYBZ(i,j,k-1))
   s11BXTYBZ(i,j,k)=s11xBXTYBZ(i,j,k)+s11yBXTYBZ(i,j,k)+s11zBXTYBZ(i,j,k)
! 
   s22xBXTYBZ(i,j,k)=(1-dafx)/c0x*s22xBXTYBZ(i,j,k)+dtx/c0x*       & 
    lamBXTYBZ(i,j,k)*(v1BXTYBZ(i,j,k)-v1BXTYBZ(i-1,j,k))
   s22yBXTYBZ(i,j,k)=(1-dafy)/c0y*s22yBXTYBZ(i,j,k)+dtx/c0y*       & 
    (lamBXTYBZ(i,j,k)+2.*muBXTYBZ(i,j,k))*(v2BXTYBZ(i,j+1,k)-v2BXTYBZ(i,j,k))
   s22zBXTYBZ(i,j,k)=(1-dafz)/c0z*s22zBXTYBZ(i,j,k)+dtx/c0z*       & 
    lamBXTYBZ(i,j,k)*(v3BXTYBZ(i,j,k)-v3BXTYBZ(i,j,k-1))
   s22BXTYBZ(i,j,k)=s22xBXTYBZ(i,j,k)+s22yBXTYBZ(i,j,k)+s22zBXTYBZ(i,j,k)
!
   s33xBXTYBZ(i,j,k)=(1-dafx)/c0x*s33xBXTYBZ(i,j,k)+dtx/c0x*       & 
    lamBXTYBZ(i,j,k)*(v1BXTYBZ(i,j,k)-v1BXTYBZ(i-1,j,k))
   s33yBXTYBZ(i,j,k)=(1-dafy)/c0y*s33yBXTYBZ(i,j,k)+dtx/c0y*       & 
    lamBXTYBZ(i,j,k)*(v2BXTYBZ(i,j+1,k)-v2BXTYBZ(i,j,k))
   s33zBXTYBZ(i,j,k)=(1-dafz)/c0z*s33zBXTYBZ(i,j,k)+dtx/c0z*       & 
   (lamBXTYBZ(i,j,k)+2*muBXTYBZ(i,j,k))*(v3BXTYBZ(i,j,k)-v3BXTYBZ(i,j,k-1))
   s33BXTYBZ(i,j,k)=s33xBXTYBZ(i,j,k)+s33yBXTYBZ(i,j,k)+s33zBXTYBZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,ny-2,3)
   c0x=1.+dafx
   dafy=DiTYS1(3,j,3)
   c0y=1.+dafy
   mu12=4./(1./muBXTYBZ(i,j,k)+1./muBXTYBZ(i,j-1,k)+1./muBXTYBZ(i+1,j,k)+ &
    1./muBXTYBZ(i+1,j-1,k))
   s12xBXTYBZ(i,j,k)=(1-dafx)*s12xBXTYBZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2BXTYBZ(i+1,j,k)-v2BXTYBZ(i,j,k))
   s12yBXTYBZ(i,j,k)=(1-dafy)*s12yBXTYBZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1BXTYBZ(i,j,k)-v1BXTYBZ(i,j-1,k))
   s12BXTYBZ(i,j,k)=s12xBXTYBZ(i,j,k)+s12yBXTYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=2,npm
  do i=1,npm-1
   dafx=DiBXS1(i,ny-2,3)
   c0x=1.+dafx
   dafz=DiBZS1(3,ny-2,k)
   c0z=1.+dafz
   mu13=4./(1./muBXTYBZ(i,j,k)+1./muBXTYBZ(i,j,k+1)+1./muBXTYBZ(i+1,j,k)+ &
    1./muBXTYBZ(i+1,j,k+1))
   s13xBXTYBZ(i,j,k)=(1-dafx)*s13xBXTYBZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3BXTYBZ(i+1,j,k)-v3BXTYBZ(i,j,k))
   s13zBXTYBZ(i,j,k)=(1-dafz)*s13zBXTYBZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1BXTYBZ(i,j,k+1)-v1BXTYBZ(i,j,k))
   s13BXTYBZ(i,j,k)=s13xBXTYBZ(i,j,k)+s13zBXTYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=2,npm
  do i=2,npm-1
   dafy=DiTYS1(3,j,3)
   c0y=1.+dafy
   dafz=DiBZS1(3,ny-2,k)
   c0z=1.+dafz
   mu23=4./(1./muBXTYBZ(i,j,k)+1./muBXTYBZ(i,j,k+1)+1./muBXTYBZ(i,j-1,k)+ &
    1./muBXTYBZ(i,j-1,k+1))
   s23yBXTYBZ(i,j,k)=(1-dafy)*s23yBXTYBZ(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3BXTYBZ(i,j,k)-v3BXTYBZ(i,j-1,k))
   s23zBXTYBZ(i,j,k)=(1-dafz)*s23zBXTYBZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2BXTYBZ(i,j,k+1)-v2BXTYBZ(i,j,k))
   s23BXTYBZ(i,j,k)=s23yBXTYBZ(i,j,k)+s23zBXTYBZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine sPMLBXBYTZ
include 'triffy.dec'
do k=2,npm
 do j=1,npm-1
  do i=2,npm-1
   dafx=DiBXS0(i,3,nz-2)
   c0x=1+dafx
   dafy=DiBYS0(3,j,nz-2)
   c0y=1+dafy
   dafz=DiTZS0(3,3,k)
   c0z=1+dafz
   s11xBXBYTZ(i,j,k)=(1-dafx)/c0x*s11xBXBYTZ(i,j,k)+dtx/c0x*       & 
    (lamBXBYTZ(i,j,k)+2*muBXBYTZ(i,j,k))*(v1BXBYTZ(i,j,k)-v1BXBYTZ(i-1,j,k))
   s11yBXBYTZ(i,j,k)=(1-dafy)/c0y*s11yBXBYTZ(i,j,k)+dtx/c0y*       & 
    lamBXBYTZ(i,j,k)*(v2BXBYTZ(i,j+1,k)-v2BXBYTZ(i,j,k))
   s11zBXBYTZ(i,j,k)=(1-dafz)/c0z*s11zBXBYTZ(i,j,k)+dtx/c0z*       & 
    lamBXBYTZ(i,j,k)*(v3BXBYTZ(i,j,k)-v3BXBYTZ(i,j,k-1))
   s11BXBYTZ(i,j,k)=s11xBXBYTZ(i,j,k)+s11yBXBYTZ(i,j,k)+s11zBXBYTZ(i,j,k)
! 
   s22xBXBYTZ(i,j,k)=(1-dafx)/c0x*s22xBXBYTZ(i,j,k)+dtx/c0x*       & 
    lamBXBYTZ(i,j,k)*(v1BXBYTZ(i,j,k)-v1BXBYTZ(i-1,j,k))
   s22yBXBYTZ(i,j,k)=(1-dafy)/c0y*s22yBXBYTZ(i,j,k)+dtx/c0y*       & 
    (lamBXBYTZ(i,j,k)+2.*muBXBYTZ(i,j,k))*(v2BXBYTZ(i,j+1,k)-v2BXBYTZ(i,j,k))
   s22zBXBYTZ(i,j,k)=(1-dafz)/c0z*s22zBXBYTZ(i,j,k)+dtx/c0z*       & 
    lamBXBYTZ(i,j,k)*(v3BXBYTZ(i,j,k)-v3BXBYTZ(i,j,k-1))
   s22BXBYTZ(i,j,k)=s22xBXBYTZ(i,j,k)+s22yBXBYTZ(i,j,k)+s22zBXBYTZ(i,j,k)
!
   s33xBXBYTZ(i,j,k)=(1-dafx)/c0x*s33xBXBYTZ(i,j,k)+dtx/c0x*       & 
    lamBXBYTZ(i,j,k)*(v1BXBYTZ(i,j,k)-v1BXBYTZ(i-1,j,k))
   s33yBXBYTZ(i,j,k)=(1-dafy)/c0y*s33yBXBYTZ(i,j,k)+dtx/c0y*       & 
    lamBXBYTZ(i,j,k)*(v2BXBYTZ(i,j+1,k)-v2BXBYTZ(i,j,k))
   s33zBXBYTZ(i,j,k)=(1-dafz)/c0z*s33zBXBYTZ(i,j,k)+dtx/c0z*       & 
   (lamBXBYTZ(i,j,k)+2*muBXBYTZ(i,j,k))*(v3BXBYTZ(i,j,k)-v3BXBYTZ(i,j,k-1))
   s33BXBYTZ(i,j,k)=s33xBXBYTZ(i,j,k)+s33yBXBYTZ(i,j,k)+s33zBXBYTZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm
 do j=2,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,3,nz-2)
   c0x=1.+dafx
   dafy=DiBYS1(3,j,nz-2)
   c0y=1.+dafy
   mu12=4./(1./muBXBYTZ(i,j,k)+1./muBXBYTZ(i,j-1,k)+1./muBXBYTZ(i+1,j,k)+ &
    1./muBXBYTZ(i+1,j-1,k))
   s12xBXBYTZ(i,j,k)=(1-dafx)*s12xBXBYTZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2BXBYTZ(i+1,j,k)-v2BXBYTZ(i,j,k))
   s12yBXBYTZ(i,j,k)=(1-dafy)*s12yBXBYTZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1BXBYTZ(i,j,k)-v1BXBYTZ(i,j-1,k))
   s12BXBYTZ(i,j,k)=s12xBXBYTZ(i,j,k)+s12yBXBYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=1,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,3,nz-2)
   c0x=1.+dafx
   dafz=DiTZS1(3,3,k)
   c0z=1.+dafz
   mu13=4./(1./muBXBYTZ(i,j,k)+1./muBXBYTZ(i,j,k+1)+1./muBXBYTZ(i+1,j,k)+ &
    1./muBXBYTZ(i+1,j,k+1))
   s13xBXBYTZ(i,j,k)=(1-dafx)*s13xBXBYTZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3BXBYTZ(i+1,j,k)-v3BXBYTZ(i,j,k))
   s13zBXBYTZ(i,j,k)=(1-dafz)*s13zBXBYTZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1BXBYTZ(i,j,k+1)-v1BXBYTZ(i,j,k))
   s13BXBYTZ(i,j,k)=s13xBXBYTZ(i,j,k)+s13zBXBYTZ(i,j,k)
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=2,npm-1
   dafy=DiBYS1(3,j,nz-2)
   c0y=1.+dafy
   dafz=DiTZS1(3,3,k)
   c0z=1.+dafz
   mu23=4./(1./muBXBYTZ(i,j,k)+1./muBXBYTZ(i,j,k+1)+1./muBXBYTZ(i,j-1,k)+ &
    1./muBXBYTZ(i,j-1,k+1))
   s23yBXBYTZ(i,j,k)=(1-dafy)*s23yBXBYTZ(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3BXBYTZ(i,j,k)-v3BXBYTZ(i,j-1,k))
   s23zBXBYTZ(i,j,k)=(1-dafz)*s23zBXBYTZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2BXBYTZ(i,j,k+1)-v2BXBYTZ(i,j,k))
   s23BXBYTZ(i,j,k)=s23yBXBYTZ(i,j,k)+s23zBXBYTZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine sPMLBXBYBZ
include 'triffy.dec'
do k=2,npm-1
 do j=1,npm-1
  do i=2,npm-1
   dafx=DiBXS0(i,3,3)
   c0x=1+dafx
   dafy=DiBYS0(3,j,3)
   c0y=1+dafy
   dafz=DiBZS0(3,3,k)
   c0z=1+dafz
   s11xBXBYBZ(i,j,k)=(1-dafx)/c0x*s11xBXBYBZ(i,j,k)+dtx/c0x*       & 
    (lamBXBYBZ(i,j,k)+2*muBXBYBZ(i,j,k))*(v1BXBYBZ(i,j,k)-v1BXBYBZ(i-1,j,k))
   s11yBXBYBZ(i,j,k)=(1-dafy)/c0y*s11yBXBYBZ(i,j,k)+dtx/c0y*       & 
    lamBXBYBZ(i,j,k)*(v2BXBYBZ(i,j+1,k)-v2BXBYBZ(i,j,k))
   s11zBXBYBZ(i,j,k)=(1-dafz)/c0z*s11zBXBYBZ(i,j,k)+dtx/c0z*       & 
    lamBXBYBZ(i,j,k)*(v3BXBYBZ(i,j,k)-v3BXBYBZ(i,j,k-1))
   s11BXBYBZ(i,j,k)=s11xBXBYBZ(i,j,k)+s11yBXBYBZ(i,j,k)+s11zBXBYBZ(i,j,k)
! 
   s22xBXBYBZ(i,j,k)=(1-dafx)/c0x*s22xBXBYBZ(i,j,k)+dtx/c0x*       & 
    lamBXBYBZ(i,j,k)*(v1BXBYBZ(i,j,k)-v1BXBYBZ(i-1,j,k))
   s22yBXBYBZ(i,j,k)=(1-dafy)/c0y*s22yBXBYBZ(i,j,k)+dtx/c0y*       & 
    (lamBXBYBZ(i,j,k)+2.*muBXBYBZ(i,j,k))*(v2BXBYBZ(i,j+1,k)-v2BXBYBZ(i,j,k))
   s22zBXBYBZ(i,j,k)=(1-dafz)/c0z*s22zBXBYBZ(i,j,k)+dtx/c0z*       & 
    lamBXBYBZ(i,j,k)*(v3BXBYBZ(i,j,k)-v3BXBYBZ(i,j,k-1))
   s22BXBYBZ(i,j,k)=s22xBXBYBZ(i,j,k)+s22yBXBYBZ(i,j,k)+s22zBXBYBZ(i,j,k)
!
   s33xBXBYBZ(i,j,k)=(1-dafx)/c0x*s33xBXBYBZ(i,j,k)+dtx/c0x*       & 
    lamBXBYBZ(i,j,k)*(v1BXBYBZ(i,j,k)-v1BXBYBZ(i-1,j,k))
   s33yBXBYBZ(i,j,k)=(1-dafy)/c0y*s33yBXBYBZ(i,j,k)+dtx/c0y*       & 
    lamBXBYBZ(i,j,k)*(v2BXBYBZ(i,j+1,k)-v2BXBYBZ(i,j,k))
   s33zBXBYBZ(i,j,k)=(1-dafz)/c0z*s33zBXBYBZ(i,j,k)+dtx/c0z*       & 
   (lamBXBYBZ(i,j,k)+2*muBXBYBZ(i,j,k))*(v3BXBYBZ(i,j,k)-v3BXBYBZ(i,j,k-1))
   s33BXBYBZ(i,j,k)=s33xBXBYBZ(i,j,k)+s33yBXBYBZ(i,j,k)+s33zBXBYBZ(i,j,k)
!
  enddo
 enddo
enddo
!
do k=2,npm-1
 do j=2,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,3,3)
   c0x=1.+dafx
   dafy=DiBYS1(3,j,3)
   c0y=1.+dafy
   mu12=4./(1./muBXBYBZ(i,j,k)+1./muBXBYBZ(i,j-1,k)+1./muBXBYBZ(i+1,j,k)+ &
    1./muBXBYBZ(i+1,j-1,k))
   s12xBXBYBZ(i,j,k)=(1-dafx)*s12xBXBYBZ(i,j,k)/c0x+dtx/c0x*mu12  &
    *(v2BXBYBZ(i+1,j,k)-v2BXBYBZ(i,j,k))
   s12yBXBYBZ(i,j,k)=(1-dafy)*s12yBXBYBZ(i,j,k)/c0y+dtx/c0y*mu12  &
    *(v1BXBYBZ(i,j,k)-v1BXBYBZ(i,j-1,k))
   s12BXBYBZ(i,j,k)=s12xBXBYBZ(i,j,k)+s12yBXBYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=1,npm-1
  do i=1,npm-1
   dafx=DiBXS1(i,3,3)
   c0x=1.+dafx
   dafz=DiBZS1(3,3,k)
   c0z=1.+dafz
   mu13=4./(1./muBXBYBZ(i,j,k)+1./muBXBYBZ(i,j,k+1)+1./muBXBYBZ(i+1,j,k)+ &
    1./muBXBYBZ(i+1,j,k+1))
   s13xBXBYBZ(i,j,k)=(1-dafx)*s13xBXBYBZ(i,j,k)/c0x+dtx/c0x*mu13  &
    *(v3BXBYBZ(i+1,j,k)-v3BXBYBZ(i,j,k))
   s13zBXBYBZ(i,j,k)=(1-dafz)*s13zBXBYBZ(i,j,k)/c0z+dtx/c0z*mu13  &
    *(v1BXBYBZ(i,j,k+1)-v1BXBYBZ(i,j,k))
   s13BXBYBZ(i,j,k)=s13xBXBYBZ(i,j,k)+s13zBXBYBZ(i,j,k)
  enddo
 enddo
enddo
!
do k=1,npm-1
 do j=2,npm-1
  do i=2,npm-1
   dafy=DiBYS1(3,j,3)
   c0y=1.+dafy
   dafz=DiBZS1(3,3,k)
   c0z=1.+dafz
   mu23=4./(1./muBXBYBZ(i,j,k)+1./muBXBYBZ(i,j,k+1)+1./muBXBYBZ(i,j-1,k)+ &
    1./muBXBYBZ(i,j-1,k+1))
   s23yBXBYBZ(i,j,k)=(1-dafy)*s23yBXBYBZ(i,j,k)/c0y+dtx/c0y*mu23  &  
    *(v3BXBYBZ(i,j,k)-v3BXBYBZ(i,j-1,k))
   s23zBXBYBZ(i,j,k)=(1-dafz)*s23zBXBYBZ(i,j,k)/c0z+dtx/c0z*mu23  &
    *(v2BXBYBZ(i,j,k+1)-v2BXBYBZ(i,j,k))
   s23BXBYBZ(i,j,k)=s23yBXBYBZ(i,j,k)+s23zBXBYBZ(i,j,k)
  enddo
 enddo
enddo
return
end subroutine
! ######################################
Subroutine PasteStressC_TX
include 'triffy.dec'
 
do k=3,nz-2
 do j=3,ny-2
  s11(nx-1,j,k)=s11TX(2,j,k)
  s11(nx,j,k)=s11TX(3,j,k)
  s11TX(1,j,k)=s11(nx-2,j,k) 
  s12(nx-1,j,k)=s12TX(2,j,k)
  s12(nx,j,k)=s12TX(3,j,k)
  s12TX(1,j,k)=s12(nx-2,j,k) 
  s13(nx-1,j,k)=s13TX(2,j,k)
  s13(nx,j,k)=s13TX(3,j,k)
  s13TX(1,j,k)=s13(nx-2,j,k) 
  s22(nx-1,j,k)=s22TX(2,j,k) 
  s22(nx,j,k)=s22TX(3,j,k) 
  s22TX(1,j,k)=s22(nx-2,j,k) 
  s23(nx-1,j,k)=s23TX(2,j,k)
  s23(nx,j,k)=s23TX(3,j,k) 
  s23TX(1,j,k)=s23(nx-2,j,k)  
  s33(nx-1,j,k)=s33TX(2,j,k) 
  s33(nx,j,k)=s33TX(3,j,k)  
  s33TX(1,j,k)=s33(nx-2,j,k)   
 enddo 
enddo 
return 
end subroutine   
! ######################################
Subroutine PasteStressC_BX
include 'triffy.dec'
 
do k=3,nz-2
 do j=3,ny-2
  s11(1,j,k)=s11BX(npm-2,j,k)
  s11(2,j,k)=s11BX(npm-1,j,k)
  s11BX(npm,j,k)=s11(3,j,k) 
  s12(1,j,k)=s12BX(npm-2,j,k)
  s12(2,j,k)=s12BX(npm-1,j,k)
  s12BX(npm,j,k)=s12(3,j,k) 
  s13(1,j,k)=s13BX(npm-2,j,k)
  s13(2,j,k)=s13BX(npm-1,j,k)
  s13BX(npm,j,k)=s13(3,j,k) 
  s22(1,j,k)=s22BX(npm-2,j,k) 
  s22(2,j,k)=s22BX(npm-1,j,k) 
  s22BX(npm,j,k)=s22(3,j,k) 
  s23(1,j,k)=s23BX(npm-2,j,k)
  s23(2,j,k)=s23BX(npm-1,j,k) 
  s23BX(npm,j,k)=s23(3,j,k)  
  s33(1,j,k)=s33BX(npm-2,j,k) 
  s33(2,j,k)=s33BX(npm-1,j,k)  
  s33BX(npm,j,k)=s33(3,j,k)   
 enddo 
enddo 
return 
end subroutine   
! ###########################################################
Subroutine PasteStressC_TY
include 'triffy.dec'
!
do k=3,nz-2
 do i=3,nx-2
  s11(i,ny-1,k)=s11TY(i,2,k)
  s11(i,ny,k)=s11TY(i,3,k)
  s11TY(i,1,k)=s11(i,ny-2,k)
  s12(i,ny-1,k)=s12TY(i,2,k)
  s12(i,ny,k)=s12TY(i,3,k)
  s12TY(i,1,k)=s12(i,ny-2,k)
  s13(i,ny-1,k)=s13TY(i,2,k)
  s13(i,ny,k)=s13TY(i,3,k)
  s13TY(i,1,k)=s13(i,ny-2,k)
  s22(i,ny-1,k)=s22TY(i,2,k) 
  s22(i,ny,k)=s22TY(i,3,k) 
  s22TY(i,1,k)=s22(i,ny-2,k)
  s23(i,ny-1,k)=s23TY(i,2,k) 
  s23(i,ny,k)=s23TY(i,3,k) 
  s23TY(i,1,k)=s23(i,ny-2,k) 
  s33(i,ny-1,k)=s33TY(i,2,k) 
  s33(i,ny,k)=s33TY(i,3,k)  
  s33TY(i,1,k)=s33(i,ny-2,k)  
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressC_BY
include 'triffy.dec'
!
do k=3,nz-2
 do i=3,nx-2
  s11(i,1,k)=s11BY(i,npm-2,k)
  s11(i,2,k)=s11BY(i,npm-1,k)
  s11BY(i,npm,k)=s11(i,3,k)
  s12(i,1,k)=s12BY(i,npm-2,k)
  s12(i,2,k)=s12BY(i,npm-1,k)
  s12BY(i,npm,k)=s12(i,3,k)
  s13(i,1,k)=s13BY(i,npm-2,k)
  s13(i,2,k)=s13BY(i,npm-1,k)
  s13BY(i,npm,k)=s13(i,3,k)
  s22(i,1,k)=s22BY(i,npm-2,k) 
  s22(i,2,k)=s22BY(i,npm-1,k) 
  s22BY(i,npm,k)=s22(i,3,k)
  s23(i,1,k)=s23BY(i,npm-2,k) 
  s23(i,2,k)=s23BY(i,npm-1,k) 
  s23BY(i,npm,k)=s23(i,3,k) 
  s33(i,1,k)=s33BY(i,npm-2,k) 
  s33(i,2,k)=s33BY(i,npm-1,k)  
  s33BY(i,npm,k)=s33(i,3,k)  
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressC_TZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=3,nx-2
  s11(i,j,nz-1)=s11TZ(i,j,2)
  s11(i,j,nz)=s11TZ(i,j,3)
  s11TZ(i,j,1)=s11(i,j,nz-2)
  s12(i,j,nz-1)=s12TZ(i,j,2) 
  s12(i,j,nz)=s12TZ(i,j,3) 
  s12TZ(i,j,1)=s12(i,j,nz-2) 
  s13(i,j,nz-1)=s13TZ(i,j,2) 
  s13(i,j,nz)=s13TZ(i,j,3) 
  s13TZ(i,j,1)=s13(i,j,nz-2) 
  s22(i,j,nz-1)=s22TZ(i,j,2) 
  s22(i,j,nz)=s22TZ(i,j,3)  
  s22TZ(i,j,1)=s22(i,j,nz-2)  
  s23(i,j,nz-1)=s23TZ(i,j,2)  
  s23(i,j,nz)=s23TZ(i,j,3)  
  s23TZ(i,j,1)=s23(i,j,nz-2)  
  s33(i,j,nz-1)=s33TZ(i,j,2)   
  s33(i,j,nz)=s33TZ(i,j,3)   
  s33TZ(i,j,1)=s33(i,j,nz-2)   
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressC_BZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=3,nx-2
  s11(i,j,1)=s11BZ(i,j,npm-2)
  s11(i,j,2)=s11BZ(i,j,npm-1)
  s11BZ(i,j,npm)=s11(i,j,3)
  s12(i,j,1)=s12BZ(i,j,npm-2) 
  s12(i,j,2)=s12BZ(i,j,npm-1) 
  s12BZ(i,j,npm)=s12(i,j,3) 
  s13(i,j,1)=s13BZ(i,j,npm-2) 
  s13(i,j,2)=s13BZ(i,j,npm-1) 
  s13BZ(i,j,npm)=s13(i,j,3) 
  s22(i,j,1)=s22BZ(i,j,npm-2) 
  s22(i,j,2)=s22BZ(i,j,npm-1)  
  s22BZ(i,j,npm)=s22(i,j,3)  
  s23(i,j,1)=s23BZ(i,j,npm-2)  
  s23(i,j,2)=s23BZ(i,j,npm-1)  
  s23BZ(i,j,npm)=s23(i,j,3)  
  s33(i,j,1)=s33BZ(i,j,npm-2)   
  s33(i,j,2)=s33BZ(i,j,npm-1)   
  s33BZ(i,j,npm)=s33(i,j,3)   
 enddo
enddo
return
end subroutine
! ########################################################### 
Subroutine PasteStressTX_TXTY 
include 'triffy.dec'
do k=3,nz-2
 do i=1,npm 
  s11TX(i,ny-1,k)=s11TXTY(i,2,k) 
  s11TX(i,ny,k)=s11TXTY(i,3,k) 
  s11TXTY(i,1,k)=s11TX(i,ny-2,k) 
  s12TX(i,ny-1,k)=s12TXTY(i,2,k)  
  s12TX(i,ny,k)=s12TXTY(i,3,k)  
  s12TXTY(i,1,k)=s12TX(i,ny-2,k)   
  s13TX(i,ny-1,k)=s13TXTY(i,2,k)  
  s13TX(i,ny,k)=s13TXTY(i,3,k)  
  s13TXTY(i,1,k)=s13TX(i,ny-2,k)   
  s22TX(i,ny-1,k)=s22TXTY(i,2,k)  
  s22TX(i,ny,k)=s22TXTY(i,3,k)  
  s22TXTY(i,1,k)=s22TX(i,ny-2,k)    
  s23TX(i,ny-1,k)=s23TXTY(i,2,k)   
  s23TX(i,ny,k)=s23TXTY(i,3,k)   
  s23TXTY(i,1,k)=s23TX(i,ny-2,k)   
  s33TX(i,ny-1,k)=s33TXTY(i,2,k)    
  s33TX(i,ny,k)=s33TXTY(i,3,k)    
  s33TXTY(i,1,k)=s33TX(i,ny-2,k)    
 enddo
enddo 
return
end subroutine 
! ########################################################### 
Subroutine PasteStressBX_BXTY 
include 'triffy.dec'
do k=3,nz-2
 do i=1,npm 
  s11BX(i,ny-1,k)=s11BXTY(i,2,k) 
  s11BX(i,ny,k)=s11BXTY(i,3,k) 
  s11BXTY(i,1,k)=s11BX(i,ny-2,k) 
  s12BX(i,ny-1,k)=s12BXTY(i,2,k)  
  s12BX(i,ny,k)=s12BXTY(i,3,k)  
  s12BXTY(i,1,k)=s12BX(i,ny-2,k)   
  s13BX(i,ny-1,k)=s13BXTY(i,2,k)  
  s13BX(i,ny,k)=s13BXTY(i,3,k)  
  s13BXTY(i,1,k)=s13BX(i,ny-2,k)   
  s22BX(i,ny-1,k)=s22BXTY(i,2,k)  
  s22BX(i,ny,k)=s22BXTY(i,3,k)  
  s22BXTY(i,1,k)=s22BX(i,ny-2,k)    
  s23BX(i,ny-1,k)=s23BXTY(i,2,k)   
  s23BX(i,ny,k)=s23BXTY(i,3,k)   
  s23BXTY(i,1,k)=s23BX(i,ny-2,k)   
  s33BX(i,ny-1,k)=s33BXTY(i,2,k)    
  s33BX(i,ny,k)=s33BXTY(i,3,k)    
  s33BXTY(i,1,k)=s33BX(i,ny-2,k)    
 enddo
enddo 
return
end subroutine 
! ########################################################### 
Subroutine PasteStressTZ_TYTZ 
include 'triffy.dec'
do k=1,npm
 do i=3,nx-2
  s11TZ(i,ny-1,k)=s11TYTZ(i,2,k) 
  s11TZ(i,ny,k)=s11TYTZ(i,3,k) 
  s11TYTZ(i,1,k)=s11TZ(i,ny-2,k) 
  s12TZ(i,ny-1,k)=s12TYTZ(i,2,k)  
  s12TZ(i,ny,k)=s12TYTZ(i,3,k)  
  s12TYTZ(i,1,k)=s12TZ(i,ny-2,k)   
  s13TZ(i,ny-1,k)=s13TYTZ(i,2,k)  
  s13TZ(i,ny,k)=s13TYTZ(i,3,k)  
  s13TYTZ(i,1,k)=s13TZ(i,ny-2,k)   
  s22TZ(i,ny-1,k)=s22TYTZ(i,2,k)  
  s22TZ(i,ny,k)=s22TYTZ(i,3,k)  
  s22TYTZ(i,1,k)=s22TZ(i,ny-2,k)    
  s23TZ(i,ny-1,k)=s23TYTZ(i,2,k)   
  s23TZ(i,ny,k)=s23TYTZ(i,3,k)   
  s23TYTZ(i,1,k)=s23TZ(i,ny-2,k)   
  s33TZ(i,ny-1,k)=s33TYTZ(i,2,k)    
  s33TZ(i,ny,k)=s33TYTZ(i,3,k)    
  s33TYTZ(i,1,k)=s33TZ(i,ny-2,k)    
 enddo
enddo 
return
end subroutine 
! ########################################################### 
Subroutine PasteStressBZ_TYBZ 
include 'triffy.dec'
do k=1,npm
 do i=3,nx-2
  s11BZ(i,ny-1,k)=s11TYBZ(i,2,k) 
  s11BZ(i,ny,k)=s11TYBZ(i,3,k) 
  s11TYBZ(i,1,k)=s11BZ(i,ny-2,k) 
  s12BZ(i,ny-1,k)=s12TYBZ(i,2,k)  
  s12BZ(i,ny,k)=s12TYBZ(i,3,k)  
  s12TYBZ(i,1,k)=s12BZ(i,ny-2,k)   
  s13BZ(i,ny-1,k)=s13TYBZ(i,2,k)  
  s13BZ(i,ny,k)=s13TYBZ(i,3,k)  
  s13TYBZ(i,1,k)=s13BZ(i,ny-2,k)   
  s22BZ(i,ny-1,k)=s22TYBZ(i,2,k)  
  s22BZ(i,ny,k)=s22TYBZ(i,3,k)  
  s22TYBZ(i,1,k)=s22BZ(i,ny-2,k)    
  s23BZ(i,ny-1,k)=s23TYBZ(i,2,k)   
  s23BZ(i,ny,k)=s23TYBZ(i,3,k)   
  s23TYBZ(i,1,k)=s23BZ(i,ny-2,k)   
  s33BZ(i,ny-1,k)=s33TYBZ(i,2,k)    
  s33BZ(i,ny,k)=s33TYBZ(i,3,k)    
  s33TYBZ(i,1,k)=s33BZ(i,ny-2,k)    
 enddo
enddo 
return
end subroutine 
! ###########################################################
Subroutine PasteStressTX_TXBY
include 'triffy.dec'
!
do k=3,nz-2
 do i=1,npm
  s11TX(i,1,k)=s11TXBY(i,npm-2,k)
  s11TX(i,2,k)=s11TXBY(i,npm-1,k)
  s11TXBY(i,npm,k)=s11TX(i,3,k)
  s12TX(i,1,k)=s12TXBY(i,npm-2,k)
  s12TX(i,2,k)=s12TXBY(i,npm-1,k)
  s12TXBY(i,npm,k)=s12TX(i,3,k)
  s13TX(i,1,k)=s13TXBY(i,npm-2,k)
  s13TX(i,2,k)=s13TXBY(i,npm-1,k)
  s13TXBY(i,npm,k)=s13TX(i,3,k)
  s22TX(i,1,k)=s22TXBY(i,npm-2,k) 
  s22TX(i,2,k)=s22TXBY(i,npm-1,k) 
  s22TXBY(i,npm,k)=s22TX(i,3,k)
  s23TX(i,1,k)=s23TXBY(i,npm-2,k) 
  s23TX(i,2,k)=s23TXBY(i,npm-1,k) 
  s23TXBY(i,npm,k)=s23TX(i,3,k) 
  s33TX(i,1,k)=s33TXBY(i,npm-2,k) 
  s33TX(i,2,k)=s33TXBY(i,npm-1,k)  
  s33TXBY(i,npm,k)=s33TX(i,3,k)  
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressBX_BXBY
include 'triffy.dec'
!
do k=3,nz-2
 do i=1,npm
  s11BX(i,1,k)=s11BXBY(i,npm-2,k)
  s11BX(i,2,k)=s11BXBY(i,npm-1,k)
  s11BXBY(i,npm,k)=s11BX(i,3,k)
  s12BX(i,1,k)=s12BXBY(i,npm-2,k)
  s12BX(i,2,k)=s12BXBY(i,npm-1,k)
  s12BXBY(i,npm,k)=s12BX(i,3,k)
  s13BX(i,1,k)=s13BXBY(i,npm-2,k)
  s13BX(i,2,k)=s13BXBY(i,npm-1,k)
  s13BXBY(i,npm,k)=s13BX(i,3,k)
  s22BX(i,1,k)=s22BXBY(i,npm-2,k) 
  s22BX(i,2,k)=s22BXBY(i,npm-1,k) 
  s22BXBY(i,npm,k)=s22BX(i,3,k)
  s23BX(i,1,k)=s23BXBY(i,npm-2,k) 
  s23BX(i,2,k)=s23BXBY(i,npm-1,k) 
  s23BXBY(i,npm,k)=s23BX(i,3,k) 
  s33BX(i,1,k)=s33BXBY(i,npm-2,k) 
  s33BX(i,2,k)=s33BXBY(i,npm-1,k)  
  s33BXBY(i,npm,k)=s33BX(i,3,k)  
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressTZ_BYTZ
include 'triffy.dec'
!
do k=1,npm
 do i=3,nx-2
  s11TZ(i,1,k)=s11BYTZ(i,npm-2,k)
  s11TZ(i,2,k)=s11BYTZ(i,npm-1,k)
  s11BYTZ(i,npm,k)=s11TZ(i,3,k)
  s12TZ(i,1,k)=s12BYTZ(i,npm-2,k)
  s12TZ(i,2,k)=s12BYTZ(i,npm-1,k)
  s12BYTZ(i,npm,k)=s12TZ(i,3,k)
  s13TZ(i,1,k)=s13BYTZ(i,npm-2,k)
  s13TZ(i,2,k)=s13BYTZ(i,npm-1,k)
  s13BYTZ(i,npm,k)=s13TZ(i,3,k)
  s22TZ(i,1,k)=s22BYTZ(i,npm-2,k) 
  s22TZ(i,2,k)=s22BYTZ(i,npm-1,k) 
  s22BYTZ(i,npm,k)=s22TZ(i,3,k)
  s23TZ(i,1,k)=s23BYTZ(i,npm-2,k) 
  s23TZ(i,2,k)=s23BYTZ(i,npm-1,k) 
  s23BYTZ(i,npm,k)=s23TZ(i,3,k) 
  s33TZ(i,1,k)=s33BYTZ(i,npm-2,k) 
  s33TZ(i,2,k)=s33BYTZ(i,npm-1,k)  
  s33BYTZ(i,npm,k)=s33TZ(i,3,k)  
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressBZ_BYBZ
include 'triffy.dec'
!
do k=1,npm
 do i=3,nx-2
  s11BZ(i,1,k)=s11BYBZ(i,npm-2,k)
  s11BZ(i,2,k)=s11BYBZ(i,npm-1,k)
  s11BYBZ(i,npm,k)=s11BZ(i,3,k)
  s12BZ(i,1,k)=s12BYBZ(i,npm-2,k)
  s12BZ(i,2,k)=s12BYBZ(i,npm-1,k)
  s12BYBZ(i,npm,k)=s12BZ(i,3,k)
  s13BZ(i,1,k)=s13BYBZ(i,npm-2,k)
  s13BZ(i,2,k)=s13BYBZ(i,npm-1,k)
  s13BYBZ(i,npm,k)=s13BZ(i,3,k)
  s22BZ(i,1,k)=s22BYBZ(i,npm-2,k) 
  s22BZ(i,2,k)=s22BYBZ(i,npm-1,k) 
  s22BYBZ(i,npm,k)=s22BZ(i,3,k)
  s23BZ(i,1,k)=s23BYBZ(i,npm-2,k) 
  s23BZ(i,2,k)=s23BYBZ(i,npm-1,k) 
  s23BYBZ(i,npm,k)=s23BZ(i,3,k) 
  s33BZ(i,1,k)=s33BYBZ(i,npm-2,k) 
  s33BZ(i,2,k)=s33BYBZ(i,npm-1,k)  
  s33BYBZ(i,npm,k)=s33BZ(i,3,k)  
 enddo
enddo
return
end subroutine
! #############################################
Subroutine PasteStressTX_TXTZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=1,npm
  s11TX(i,j,nz-1)=s11TXTZ(i,j,2)
  s11TX(i,j,nz)=s11TXTZ(i,j,3) 
  s11TXTZ(i,j,1)=s11TX(i,j,nz-2) 
  s12TX(i,j,nz-1)=s12TXTZ(i,j,2)
  s12TX(i,j,nz)=s12TXTZ(i,j,3) 
  s12TXTZ(i,j,1)=s12TX(i,j,nz-2)   
  s13TX(i,j,nz-1)=s13TXTZ(i,j,2)
  s13TX(i,j,nz)=s13TXTZ(i,j,3) 
  s13TXTZ(i,j,1)=s13TX(i,j,nz-2)   
  s22TX(i,j,nz-1)=s22TXTZ(i,j,2)
  s22TX(i,j,nz)=s22TXTZ(i,j,3)  
  s22TXTZ(i,j,1)=s22TX(i,j,nz-2)    
  s23TX(i,j,nz-1)=s23TXTZ(i,j,2) 
  s23TX(i,j,nz)=s23TXTZ(i,j,3)  
  s23TXTZ(i,j,1)=s23TX(i,j,nz-2)    
  s33TX(i,j,nz-1)=s33TXTZ(i,j,2)   
  s33TX(i,j,nz)=s33TXTZ(i,j,3)   
  s33TXTZ(i,j,1)=s33TX(i,j,nz-2) 
 enddo 
enddo  
return 
end subroutine   
! #############################################
Subroutine PasteStressBX_BXTZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=1,npm
  s11BX(i,j,nz-1)=s11BXTZ(i,j,2)
  s11BX(i,j,nz)=s11BXTZ(i,j,3) 
  s11BXTZ(i,j,1)=s11BX(i,j,nz-2) 
  s12BX(i,j,nz-1)=s12BXTZ(i,j,2)
  s12BX(i,j,nz)=s12BXTZ(i,j,3) 
  s12BXTZ(i,j,1)=s12BX(i,j,nz-2)   
  s13BX(i,j,nz-1)=s13BXTZ(i,j,2)
  s13BX(i,j,nz)=s13BXTZ(i,j,3) 
  s13BXTZ(i,j,1)=s13BX(i,j,nz-2)   
  s22BX(i,j,nz-1)=s22BXTZ(i,j,2)
  s22BX(i,j,nz)=s22BXTZ(i,j,3)  
  s22BXTZ(i,j,1)=s22BX(i,j,nz-2)    
  s23BX(i,j,nz-1)=s23BXTZ(i,j,2) 
  s23BX(i,j,nz)=s23BXTZ(i,j,3)  
  s23BXTZ(i,j,1)=s23BX(i,j,nz-2)    
  s33BX(i,j,nz-1)=s33BXTZ(i,j,2)   
  s33BX(i,j,nz)=s33BXTZ(i,j,3)   
  s33BXTZ(i,j,1)=s33BX(i,j,nz-2) 
 enddo 
enddo  
return 
end subroutine   
! #############################################
Subroutine PasteStressTY_TYTZ
include 'triffy.dec'
 
do j=1,npm
 do i=3,nx-2
  s11TY(i,j,nz-1)=s11TYTZ(i,j,2)
  s11TY(i,j,nz)=s11TYTZ(i,j,3) 
  s11TYTZ(i,j,1)=s11TY(i,j,nz-2) 
  s12TY(i,j,nz-1)=s12TYTZ(i,j,2)
  s12TY(i,j,nz)=s12TYTZ(i,j,3) 
  s12TYTZ(i,j,1)=s12TY(i,j,nz-2)   
  s13TY(i,j,nz-1)=s13TYTZ(i,j,2)
  s13TY(i,j,nz)=s13TYTZ(i,j,3) 
  s13TYTZ(i,j,1)=s13TY(i,j,nz-2)   
  s22TY(i,j,nz-1)=s22TYTZ(i,j,2)
  s22TY(i,j,nz)=s22TYTZ(i,j,3)  
  s22TYTZ(i,j,1)=s22TY(i,j,nz-2)    
  s23TY(i,j,nz-1)=s23TYTZ(i,j,2) 
  s23TY(i,j,nz)=s23TYTZ(i,j,3)  
  s23TYTZ(i,j,1)=s23TY(i,j,nz-2)    
  s33TY(i,j,nz-1)=s33TYTZ(i,j,2)   
  s33TY(i,j,nz)=s33TYTZ(i,j,3)   
  s33TYTZ(i,j,1)=s33TY(i,j,nz-2) 
 enddo 
enddo  
return 
end subroutine   
! #############################################
Subroutine PasteStressBY_BYTZ
include 'triffy.dec'
 
do j=1,npm
 do i=3,nx-2
  s11BY(i,j,nz-1)=s11BYTZ(i,j,2)
  s11BY(i,j,nz)=s11BYTZ(i,j,3) 
  s11BYTZ(i,j,1)=s11BY(i,j,nz-2) 
  s12BY(i,j,nz-1)=s12BYTZ(i,j,2)
  s12BY(i,j,nz)=s12BYTZ(i,j,3) 
  s12BYTZ(i,j,1)=s12BY(i,j,nz-2)   
  s13BY(i,j,nz-1)=s13BYTZ(i,j,2)
  s13BY(i,j,nz)=s13BYTZ(i,j,3) 
  s13BYTZ(i,j,1)=s13BY(i,j,nz-2)   
  s22BY(i,j,nz-1)=s22BYTZ(i,j,2)
  s22BY(i,j,nz)=s22BYTZ(i,j,3)  
  s22BYTZ(i,j,1)=s22BY(i,j,nz-2)    
  s23BY(i,j,nz-1)=s23BYTZ(i,j,2) 
  s23BY(i,j,nz)=s23BYTZ(i,j,3)  
  s23BYTZ(i,j,1)=s23BY(i,j,nz-2)    
  s33BY(i,j,nz-1)=s33BYTZ(i,j,2)   
  s33BY(i,j,nz)=s33BYTZ(i,j,3)   
  s33BYTZ(i,j,1)=s33BY(i,j,nz-2) 
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteStressTX_TXBZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=1,npm
  s11TX(i,j,1)=s11TXBZ(i,j,npm-2)
  s11TX(i,j,2)=s11TXBZ(i,j,npm-1)
  s11TXBZ(i,j,npm)=s11TX(i,j,3)
  s12TX(i,j,1)=s12TXBZ(i,j,npm-2) 
  s12TX(i,j,2)=s12TXBZ(i,j,npm-1) 
  s12TXBZ(i,j,npm)=s12TX(i,j,3) 
  s13TX(i,j,1)=s13TXBZ(i,j,npm-2) 
  s13TX(i,j,2)=s13TXBZ(i,j,npm-1) 
  s13TXBZ(i,j,npm)=s13TX(i,j,3) 
  s22TX(i,j,1)=s22TXBZ(i,j,npm-2) 
  s22TX(i,j,2)=s22TXBZ(i,j,npm-1)  
  s22TXBZ(i,j,npm)=s22TX(i,j,3)  
  s23TX(i,j,1)=s23TXBZ(i,j,npm-2)  
  s23TX(i,j,2)=s23TXBZ(i,j,npm-1)  
  s23TXBZ(i,j,npm)=s23TX(i,j,3)  
  s33TX(i,j,1)=s33TXBZ(i,j,npm-2)   
  s33TX(i,j,2)=s33TXBZ(i,j,npm-1)   
  s33TXBZ(i,j,npm)=s33TX(i,j,3)   
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressBX_BXBZ
include 'triffy.dec'
 
do j=3,ny-2
 do i=1,npm
  s11BX(i,j,1)=s11BXBZ(i,j,npm-2)
  s11BX(i,j,2)=s11BXBZ(i,j,npm-1)
  s11BXBZ(i,j,npm)=s11BX(i,j,3)
  s12BX(i,j,1)=s12BXBZ(i,j,npm-2) 
  s12BX(i,j,2)=s12BXBZ(i,j,npm-1) 
  s12BXBZ(i,j,npm)=s12BX(i,j,3) 
  s13BX(i,j,1)=s13BXBZ(i,j,npm-2) 
  s13BX(i,j,2)=s13BXBZ(i,j,npm-1) 
  s13BXBZ(i,j,npm)=s13BX(i,j,3) 
  s22BX(i,j,1)=s22BXBZ(i,j,npm-2) 
  s22BX(i,j,2)=s22BXBZ(i,j,npm-1)  
  s22BXBZ(i,j,npm)=s22BX(i,j,3)  
  s23BX(i,j,1)=s23BXBZ(i,j,npm-2)  
  s23BX(i,j,2)=s23BXBZ(i,j,npm-1)  
  s23BXBZ(i,j,npm)=s23BX(i,j,3)  
  s33BX(i,j,1)=s33BXBZ(i,j,npm-2)   
  s33BX(i,j,2)=s33BXBZ(i,j,npm-1)   
  s33BXBZ(i,j,npm)=s33BX(i,j,3)   
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressTY_TYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=3,nx-2
  s11TY(i,j,1)=s11TYBZ(i,j,npm-2)
  s11TY(i,j,2)=s11TYBZ(i,j,npm-1)
  s11TYBZ(i,j,npm)=s11TY(i,j,3)
  s12TY(i,j,1)=s12TYBZ(i,j,npm-2) 
  s12TY(i,j,2)=s12TYBZ(i,j,npm-1) 
  s12TYBZ(i,j,npm)=s12TY(i,j,3) 
  s13TY(i,j,1)=s13TYBZ(i,j,npm-2) 
  s13TY(i,j,2)=s13TYBZ(i,j,npm-1) 
  s13TYBZ(i,j,npm)=s13TY(i,j,3) 
  s22TY(i,j,1)=s22TYBZ(i,j,npm-2) 
  s22TY(i,j,2)=s22TYBZ(i,j,npm-1)  
  s22TYBZ(i,j,npm)=s22TY(i,j,3)  
  s23TY(i,j,1)=s23TYBZ(i,j,npm-2)  
  s23TY(i,j,2)=s23TYBZ(i,j,npm-1)  
  s23TYBZ(i,j,npm)=s23TY(i,j,3)  
  s33TY(i,j,1)=s33TYBZ(i,j,npm-2)   
  s33TY(i,j,2)=s33TYBZ(i,j,npm-1)   
  s33TYBZ(i,j,npm)=s33TY(i,j,3)   
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressBY_BYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=3,nx-2
  s11BY(i,j,1)=s11BYBZ(i,j,npm-2)
  s11BY(i,j,2)=s11BYBZ(i,j,npm-1)
  s11BYBZ(i,j,npm)=s11BY(i,j,3)
  s12BY(i,j,1)=s12BYBZ(i,j,npm-2) 
  s12BY(i,j,2)=s12BYBZ(i,j,npm-1) 
  s12BYBZ(i,j,npm)=s12BY(i,j,3) 
  s13BY(i,j,1)=s13BYBZ(i,j,npm-2) 
  s13BY(i,j,2)=s13BYBZ(i,j,npm-1) 
  s13BYBZ(i,j,npm)=s13BY(i,j,3) 
  s22BY(i,j,1)=s22BYBZ(i,j,npm-2) 
  s22BY(i,j,2)=s22BYBZ(i,j,npm-1)  
  s22BYBZ(i,j,npm)=s22BY(i,j,3)  
  s23BY(i,j,1)=s23BYBZ(i,j,npm-2)  
  s23BY(i,j,2)=s23BYBZ(i,j,npm-1)  
  s23BYBZ(i,j,npm)=s23BY(i,j,3)  
  s33BY(i,j,1)=s33BYBZ(i,j,npm-2)   
  s33BY(i,j,2)=s33BYBZ(i,j,npm-1)   
  s33BYBZ(i,j,npm)=s33BY(i,j,3)   
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressTY_TXTY
include 'triffy.dec'
do k=3,nz-2
 do j=1,npm
  s11TY(nx-1,j,k)=s11TXTY(2,j,k)
  s11TY(nx,j,k)=s11TXTY(3,j,k)
  s11TXTY(1,j,k)=s11TY(nx-2,j,k)
  s12TY(nx-1,j,k)=s12TXTY(2,j,k) 
  s12TY(nx,j,k)=s12TXTY(3,j,k) 
  s12TXTY(1,j,k)=s12TY(nx-2,j,k)   
  s13TY(nx-1,j,k)=s13TXTY(2,j,k) 
  s13TY(nx,j,k)=s13TXTY(3,j,k) 
  s13TXTY(1,j,k)=s13TY(nx-2,j,k)   
  s22TY(nx-1,j,k)=s22TXTY(2,j,k)  
  s22TY(nx,j,k)=s22TXTY(3,j,k)  
  s22TXTY(1,j,k)=s22TY(nx-2,j,k)    
  s23TY(nx-1,j,k)=s23TXTY(2,j,k)  
  s23TY(nx,j,k)=s23TXTY(3,j,k)  
  s23TXTY(1,j,k)=s23TY(nx-2,j,k)    
  s33TY(nx-1,j,k)=s33TXTY(2,j,k)   
  s33TY(nx,j,k)=s33TXTY(3,j,k)   
  s33TXTY(1,j,k)=s33TY(nx-2,j,k) 
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteStressBY_TXBY
include 'triffy.dec'
do k=3,nz-2
 do j=1,npm
  s11BY(nx-1,j,k)=s11TXBY(2,j,k)
  s11BY(nx,j,k)=s11TXBY(3,j,k)
  s11TXBY(1,j,k)=s11BY(nx-2,j,k)
  s12BY(nx-1,j,k)=s12TXBY(2,j,k) 
  s12BY(nx,j,k)=s12TXBY(3,j,k) 
  s12TXBY(1,j,k)=s12BY(nx-2,j,k)   
  s13BY(nx-1,j,k)=s13TXBY(2,j,k) 
  s13BY(nx,j,k)=s13TXBY(3,j,k) 
  s13TXBY(1,j,k)=s13BY(nx-2,j,k)   
  s22BY(nx-1,j,k)=s22TXBY(2,j,k)  
  s22BY(nx,j,k)=s22TXBY(3,j,k)  
  s22TXBY(1,j,k)=s22BY(nx-2,j,k)    
  s23BY(nx-1,j,k)=s23TXBY(2,j,k)  
  s23BY(nx,j,k)=s23TXBY(3,j,k)  
  s23TXBY(1,j,k)=s23BY(nx-2,j,k)    
  s33BY(nx-1,j,k)=s33TXBY(2,j,k)   
  s33BY(nx,j,k)=s33TXBY(3,j,k)   
  s33TXBY(1,j,k)=s33BY(nx-2,j,k) 
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteStressTZ_TXTZ
include 'triffy.dec'
do k=1,npm
 do j=3,ny-2
  s11TZ(nx-1,j,k)=s11TXTZ(2,j,k)
  s11TZ(nx,j,k)=s11TXTZ(3,j,k)
  s11TXTZ(1,j,k)=s11TZ(nx-2,j,k)
  s12TZ(nx-1,j,k)=s12TXTZ(2,j,k) 
  s12TZ(nx,j,k)=s12TXTZ(3,j,k) 
  s12TXTZ(1,j,k)=s12TZ(nx-2,j,k)   
  s13TZ(nx-1,j,k)=s13TXTZ(2,j,k) 
  s13TZ(nx,j,k)=s13TXTZ(3,j,k) 
  s13TXTZ(1,j,k)=s13TZ(nx-2,j,k)   
  s22TZ(nx-1,j,k)=s22TXTZ(2,j,k)  
  s22TZ(nx,j,k)=s22TXTZ(3,j,k)  
  s22TXTZ(1,j,k)=s22TZ(nx-2,j,k)    
  s23TZ(nx-1,j,k)=s23TXTZ(2,j,k)  
  s23TZ(nx,j,k)=s23TXTZ(3,j,k)  
  s23TXTZ(1,j,k)=s23TZ(nx-2,j,k)    
  s33TZ(nx-1,j,k)=s33TXTZ(2,j,k)   
  s33TZ(nx,j,k)=s33TXTZ(3,j,k)   
  s33TXTZ(1,j,k)=s33TZ(nx-2,j,k) 
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteStressBZ_TXBZ
include 'triffy.dec'
do k=1,npm
 do j=3,ny-2
  s11BZ(nx-1,j,k)=s11TXBZ(2,j,k)
  s11BZ(nx,j,k)=s11TXBZ(3,j,k)
  s11TXBZ(1,j,k)=s11BZ(nx-2,j,k)
  s12BZ(nx-1,j,k)=s12TXBZ(2,j,k) 
  s12BZ(nx,j,k)=s12TXBZ(3,j,k) 
  s12TXBZ(1,j,k)=s12BZ(nx-2,j,k)   
  s13BZ(nx-1,j,k)=s13TXBZ(2,j,k) 
  s13BZ(nx,j,k)=s13TXBZ(3,j,k) 
  s13TXBZ(1,j,k)=s13BZ(nx-2,j,k)   
  s22BZ(nx-1,j,k)=s22TXBZ(2,j,k)  
  s22BZ(nx,j,k)=s22TXBZ(3,j,k)  
  s22TXBZ(1,j,k)=s22BZ(nx-2,j,k)    
  s23BZ(nx-1,j,k)=s23TXBZ(2,j,k)  
  s23BZ(nx,j,k)=s23TXBZ(3,j,k)  
  s23TXBZ(1,j,k)=s23BZ(nx-2,j,k)    
  s33BZ(nx-1,j,k)=s33TXBZ(2,j,k)   
  s33BZ(nx,j,k)=s33TXBZ(3,j,k)   
  s33TXBZ(1,j,k)=s33BZ(nx-2,j,k) 
 enddo 
enddo  
return 
end subroutine   
! ######################################
Subroutine PasteStressTY_BXTY
include 'triffy.dec'
 
do k=3,nz-2
 do j=1,npm
  s11TY(1,j,k)=s11BXTY(npm-2,j,k)
  s11TY(2,j,k)=s11BXTY(npm-1,j,k)
  s11BXTY(npm,j,k)=s11TY(3,j,k) 
  s12TY(1,j,k)=s12BXTY(npm-2,j,k)
  s12TY(2,j,k)=s12BXTY(npm-1,j,k)
  s12BXTY(npm,j,k)=s12TY(3,j,k) 
  s13TY(1,j,k)=s13BXTY(npm-2,j,k)
  s13TY(2,j,k)=s13BXTY(npm-1,j,k)
  s13BXTY(npm,j,k)=s13TY(3,j,k) 
  s22TY(1,j,k)=s22BXTY(npm-2,j,k) 
  s22TY(2,j,k)=s22BXTY(npm-1,j,k) 
  s22BXTY(npm,j,k)=s22TY(3,j,k) 
  s23TY(1,j,k)=s23BXTY(npm-2,j,k)
  s23TY(2,j,k)=s23BXTY(npm-1,j,k) 
  s23BXTY(npm,j,k)=s23TY(3,j,k)  
  s33TY(1,j,k)=s33BXTY(npm-2,j,k) 
  s33TY(2,j,k)=s33BXTY(npm-1,j,k)  
  s33BXTY(npm,j,k)=s33TY(3,j,k)   
 enddo 
enddo 
return 
end subroutine   
! ######################################
Subroutine PasteStressBY_BXBY
include 'triffy.dec'
 
do k=3,nz-2
 do j=1,npm
  s11BY(1,j,k)=s11BXBY(npm-2,j,k)
  s11BY(2,j,k)=s11BXBY(npm-1,j,k)
  s11BXBY(npm,j,k)=s11BY(3,j,k) 
  s12BY(1,j,k)=s12BXBY(npm-2,j,k)
  s12BY(2,j,k)=s12BXBY(npm-1,j,k)
  s12BXBY(npm,j,k)=s12BY(3,j,k) 
  s13BY(1,j,k)=s13BXBY(npm-2,j,k)
  s13BY(2,j,k)=s13BXBY(npm-1,j,k)
  s13BXBY(npm,j,k)=s13BY(3,j,k) 
  s22BY(1,j,k)=s22BXBY(npm-2,j,k) 
  s22BY(2,j,k)=s22BXBY(npm-1,j,k) 
  s22BXBY(npm,j,k)=s22BY(3,j,k) 
  s23BY(1,j,k)=s23BXBY(npm-2,j,k)
  s23BY(2,j,k)=s23BXBY(npm-1,j,k) 
  s23BXBY(npm,j,k)=s23BY(3,j,k)  
  s33BY(1,j,k)=s33BXBY(npm-2,j,k) 
  s33BY(2,j,k)=s33BXBY(npm-1,j,k)  
  s33BXBY(npm,j,k)=s33BY(3,j,k)   
 enddo 
enddo 
return 
end subroutine   
! ######################################
Subroutine PasteStressTZ_BXTZ
include 'triffy.dec'
 
do k=1,npm
 do j=3,ny-2
  s11TZ(1,j,k)=s11BXTZ(npm-2,j,k)
  s11TZ(2,j,k)=s11BXTZ(npm-1,j,k)
  s11BXTZ(npm,j,k)=s11TZ(3,j,k) 
  s12TZ(1,j,k)=s12BXTZ(npm-2,j,k)
  s12TZ(2,j,k)=s12BXTZ(npm-1,j,k)
  s12BXTZ(npm,j,k)=s12TZ(3,j,k) 
  s13TZ(1,j,k)=s13BXTZ(npm-2,j,k)
  s13TZ(2,j,k)=s13BXTZ(npm-1,j,k)
  s13BXTZ(npm,j,k)=s13TZ(3,j,k) 
  s22TZ(1,j,k)=s22BXTZ(npm-2,j,k) 
  s22TZ(2,j,k)=s22BXTZ(npm-1,j,k) 
  s22BXTZ(npm,j,k)=s22TZ(3,j,k) 
  s23TZ(1,j,k)=s23BXTZ(npm-2,j,k)
  s23TZ(2,j,k)=s23BXTZ(npm-1,j,k) 
  s23BXTZ(npm,j,k)=s23TZ(3,j,k)  
  s33TZ(1,j,k)=s33BXTZ(npm-2,j,k) 
  s33TZ(2,j,k)=s33BXTZ(npm-1,j,k)  
  s33BXTZ(npm,j,k)=s33TZ(3,j,k)   
 enddo 
enddo 
return 
end subroutine   
! ######################################
Subroutine PasteStressBZ_BXBZ
include 'triffy.dec'
 
do k=1,npm
 do j=3,ny-2
  s11BZ(1,j,k)=s11BXBZ(npm-2,j,k)
  s11BZ(2,j,k)=s11BXBZ(npm-1,j,k)
  s11BXBZ(npm,j,k)=s11BZ(3,j,k) 
  s12BZ(1,j,k)=s12BXBZ(npm-2,j,k)
  s12BZ(2,j,k)=s12BXBZ(npm-1,j,k)
  s12BXBZ(npm,j,k)=s12BZ(3,j,k) 
  s13BZ(1,j,k)=s13BXBZ(npm-2,j,k)
  s13BZ(2,j,k)=s13BXBZ(npm-1,j,k)
  s13BXBZ(npm,j,k)=s13BZ(3,j,k) 
  s22BZ(1,j,k)=s22BXBZ(npm-2,j,k) 
  s22BZ(2,j,k)=s22BXBZ(npm-1,j,k) 
  s22BXBZ(npm,j,k)=s22BZ(3,j,k) 
  s23BZ(1,j,k)=s23BXBZ(npm-2,j,k)
  s23BZ(2,j,k)=s23BXBZ(npm-1,j,k) 
  s23BXBZ(npm,j,k)=s23BZ(3,j,k)  
  s33BZ(1,j,k)=s33BXBZ(npm-2,j,k) 
  s33BZ(2,j,k)=s33BXBZ(npm-1,j,k)  
  s33BXBZ(npm,j,k)=s33BZ(3,j,k)   
 enddo 
enddo 
return 
end subroutine
! ########################################################### 
Subroutine PasteStressTXTY_TXTYTZ
include 'triffy.dec'
do j=1,npm 
 do i=1,npm 
  s11TXTY(i,j,nz-1)=s11TXTYTZ(i,j,2) 
  s11TXTY(i,j,nz)=s11TXTYTZ(i,j,3)
  s11TXTYTZ(i,j,1)=s11TXTY(i,j,nz-2) 
  s12TXTY(i,j,nz-1)=s12TXTYTZ(i,j,2) 
  s12TXTY(i,j,nz)=s12TXTYTZ(i,j,3)
  s12TXTYTZ(i,j,1)=s12TXTY(i,j,nz-2) 
  s13TXTY(i,j,nz-1)=s13TXTYTZ(i,j,2) 
  s13TXTY(i,j,nz)=s13TXTYTZ(i,j,3)
  s13TXTYTZ(i,j,1)=s13TXTY(i,j,nz-2) 
  s22TXTY(i,j,nz-1)=s22TXTYTZ(i,j,2) 
  s22TXTY(i,j,nz)=s22TXTYTZ(i,j,3) 
  s22TXTYTZ(i,j,1)=s22TXTY(i,j,nz-2)  
  s23TXTY(i,j,nz-1)=s23TXTYTZ(i,j,2)  
  s23TXTY(i,j,nz)=s23TXTYTZ(i,j,3) 
  s23TXTYTZ(i,j,1)=s23TXTY(i,j,nz-2)  
  s33TXTY(i,j,nz-1)=s33TXTYTZ(i,j,2)   
  s33TXTY(i,j,nz)=s33TXTYTZ(i,j,3)  
  s33TXTYTZ(i,j,1)=s33TXTY(i,j,nz-2)
 enddo
enddo
return 
end subroutine 
! ########################################################### 
Subroutine PasteStressTXTZ_TXTYTZ 
include 'triffy.dec'
do k=1,npm
 do i=1,npm 
  s11TXTZ(i,ny-1,k)=s11TXTYTZ(i,2,k) 
  s11TXTZ(i,ny,k)=s11TXTYTZ(i,3,k) 
  s11TXTYTZ(i,1,k)=s11TXTZ(i,ny-2,k) 
  s12TXTZ(i,ny-1,k)=s12TXTYTZ(i,2,k)  
  s12TXTZ(i,ny,k)=s12TXTYTZ(i,3,k)  
  s12TXTYTZ(i,1,k)=s12TXTZ(i,ny-2,k)   
  s13TXTZ(i,ny-1,k)=s13TXTYTZ(i,2,k)  
  s13TXTZ(i,ny,k)=s13TXTYTZ(i,3,k)  
  s13TXTYTZ(i,1,k)=s13TXTZ(i,ny-2,k)   
  s22TXTZ(i,ny-1,k)=s22TXTYTZ(i,2,k)  
  s22TXTZ(i,ny,k)=s22TXTYTZ(i,3,k)  
  s22TXTYTZ(i,1,k)=s22TXTZ(i,ny-2,k)    
  s23TXTZ(i,ny-1,k)=s23TXTYTZ(i,2,k)   
  s23TXTZ(i,ny,k)=s23TXTYTZ(i,3,k)   
  s23TXTYTZ(i,1,k)=s23TXTZ(i,ny-2,k)   
  s33TXTZ(i,ny-1,k)=s33TXTYTZ(i,2,k)    
  s33TXTZ(i,ny,k)=s33TXTYTZ(i,3,k)    
  s33TXTYTZ(i,1,k)=s33TXTZ(i,ny-2,k)    
 enddo
enddo 
return
end subroutine 
! ###########################################################
Subroutine PasteStressTYTZ_TXTYTZ
include 'triffy.dec'
do k=1,npm
 do j=1,npm
  s11TYTZ(nx-1,j,k)=s11TXTYTZ(2,j,k)
  s11TYTZ(nx,j,k)=s11TXTYTZ(3,j,k)
  s11TXTYTZ(1,j,k)=s11TYTZ(nx-2,j,k)
  s12TYTZ(nx-1,j,k)=s12TXTYTZ(2,j,k) 
  s12TYTZ(nx,j,k)=s12TXTYTZ(3,j,k) 
  s12TXTYTZ(1,j,k)=s12TYTZ(nx-2,j,k)   
  s13TYTZ(nx-1,j,k)=s13TXTYTZ(2,j,k) 
  s13TYTZ(nx,j,k)=s13TXTYTZ(3,j,k) 
  s13TXTYTZ(1,j,k)=s13TYTZ(nx-2,j,k)   
  s22TYTZ(nx-1,j,k)=s22TXTYTZ(2,j,k)  
  s22TYTZ(nx,j,k)=s22TXTYTZ(3,j,k)  
  s22TXTYTZ(1,j,k)=s22TYTZ(nx-2,j,k)    
  s23TYTZ(nx-1,j,k)=s23TXTYTZ(2,j,k)  
  s23TYTZ(nx,j,k)=s23TXTYTZ(3,j,k)  
  s23TXTYTZ(1,j,k)=s23TYTZ(nx-2,j,k)    
  s33TYTZ(nx-1,j,k)=s33TXTYTZ(2,j,k)   
  s33TYTZ(nx,j,k)=s33TXTYTZ(3,j,k)   
  s33TXTYTZ(1,j,k)=s33TYTZ(nx-2,j,k) 
 enddo 
enddo  
return 
end subroutine   
! ########################################################### 
Subroutine PasteStressTXTY_TXTYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  s11TXTY(i,j,1)=s11TXTYBZ(i,j,npm-2)
  s11TXTY(i,j,2)=s11TXTYBZ(i,j,npm-1)
  s11TXTYBZ(i,j,npm)=s11TXTY(i,j,3)
  s12TXTY(i,j,1)=s12TXTYBZ(i,j,npm-2) 
  s12TXTY(i,j,2)=s12TXTYBZ(i,j,npm-1) 
  s12TXTYBZ(i,j,npm)=s12TXTY(i,j,3) 
  s13TXTY(i,j,1)=s13TXTYBZ(i,j,npm-2) 
  s13TXTY(i,j,2)=s13TXTYBZ(i,j,npm-1) 
  s13TXTYBZ(i,j,npm)=s13TXTY(i,j,3) 
  s22TXTY(i,j,1)=s22TXTYBZ(i,j,npm-2) 
  s22TXTY(i,j,2)=s22TXTYBZ(i,j,npm-1)  
  s22TXTYBZ(i,j,npm)=s22TXTY(i,j,3)  
  s23TXTY(i,j,1)=s23TXTYBZ(i,j,npm-2)  
  s23TXTY(i,j,2)=s23TXTYBZ(i,j,npm-1)  
  s23TXTYBZ(i,j,npm)=s23TXTY(i,j,3)  
  s33TXTY(i,j,1)=s33TXTYBZ(i,j,npm-2)   
  s33TXTY(i,j,2)=s33TXTYBZ(i,j,npm-1)   
  s33TXTYBZ(i,j,npm)=s33TXTY(i,j,3)   
 enddo
enddo
return
end subroutine
! ########################################################### 
Subroutine PasteStressTXBZ_TXTYBZ 
include 'triffy.dec'
do k=1,npm
 do i=1,npm 
  s11TXBZ(i,ny-1,k)=s11TXTYBZ(i,2,k) 
  s11TXBZ(i,ny,k)=s11TXTYBZ(i,3,k) 
  s11TXTYBZ(i,1,k)=s11TXBZ(i,ny-2,k) 
  s12TXBZ(i,ny-1,k)=s12TXTYBZ(i,2,k)  
  s12TXBZ(i,ny,k)=s12TXTYBZ(i,3,k)  
  s12TXTYBZ(i,1,k)=s12TXBZ(i,ny-2,k)   
  s13TXBZ(i,ny-1,k)=s13TXTYBZ(i,2,k)  
  s13TXBZ(i,ny,k)=s13TXTYBZ(i,3,k)  
  s13TXTYBZ(i,1,k)=s13TXBZ(i,ny-2,k)   
  s22TXBZ(i,ny-1,k)=s22TXTYBZ(i,2,k)  
  s22TXBZ(i,ny,k)=s22TXTYBZ(i,3,k)  
  s22TXTYBZ(i,1,k)=s22TXBZ(i,ny-2,k)    
  s23TXBZ(i,ny-1,k)=s23TXTYBZ(i,2,k)   
  s23TXBZ(i,ny,k)=s23TXTYBZ(i,3,k)   
  s23TXTYBZ(i,1,k)=s23TXBZ(i,ny-2,k)   
  s33TXBZ(i,ny-1,k)=s33TXTYBZ(i,2,k)    
  s33TXBZ(i,ny,k)=s33TXTYBZ(i,3,k)    
  s33TXTYBZ(i,1,k)=s33TXBZ(i,ny-2,k)    
 enddo
enddo 
return
end subroutine 
! ###########################################################
Subroutine PasteStressTYBZ_TXTYBZ
include 'triffy.dec'
do k=1,npm
 do j=1,npm
  s11TYBZ(nx-1,j,k)=s11TXTYBZ(2,j,k)
  s11TYBZ(nx,j,k)=s11TXTYBZ(3,j,k)
  s11TXTYBZ(1,j,k)=s11TYBZ(nx-2,j,k)
  s12TYBZ(nx-1,j,k)=s12TXTYBZ(2,j,k) 
  s12TYBZ(nx,j,k)=s12TXTYBZ(3,j,k) 
  s12TXTYBZ(1,j,k)=s12TYBZ(nx-2,j,k)   
  s13TYBZ(nx-1,j,k)=s13TXTYBZ(2,j,k) 
  s13TYBZ(nx,j,k)=s13TXTYBZ(3,j,k) 
  s13TXTYBZ(1,j,k)=s13TYBZ(nx-2,j,k)   
  s22TYBZ(nx-1,j,k)=s22TXTYBZ(2,j,k)  
  s22TYBZ(nx,j,k)=s22TXTYBZ(3,j,k)  
  s22TXTYBZ(1,j,k)=s22TYBZ(nx-2,j,k)    
  s23TYBZ(nx-1,j,k)=s23TXTYBZ(2,j,k)  
  s23TYBZ(nx,j,k)=s23TXTYBZ(3,j,k)  
  s23TXTYBZ(1,j,k)=s23TYBZ(nx-2,j,k)    
  s33TYBZ(nx-1,j,k)=s33TXTYBZ(2,j,k)   
  s33TYBZ(nx,j,k)=s33TXTYBZ(3,j,k)   
  s33TXTYBZ(1,j,k)=s33TYBZ(nx-2,j,k) 
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteStressTXBY_TXBYTZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  s11TXBY(i,j,nz-1)=s11TXBYTZ(i,j,2)
  s11TXBY(i,j,nz)=s11TXBYTZ(i,j,3) 
  s11TXBYTZ(i,j,1)=s11TXBY(i,j,nz-2) 
  s12TXBY(i,j,nz-1)=s12TXBYTZ(i,j,2)
  s12TXBY(i,j,nz)=s12TXBYTZ(i,j,3) 
  s12TXBYTZ(i,j,1)=s12TXBY(i,j,nz-2)   
  s13TXBY(i,j,nz-1)=s13TXBYTZ(i,j,2)
  s13TXBY(i,j,nz)=s13TXBYTZ(i,j,3) 
  s13TXBYTZ(i,j,1)=s13TXBY(i,j,nz-2)   
  s22TXBY(i,j,nz-1)=s22TXBYTZ(i,j,2)
  s22TXBY(i,j,nz)=s22TXBYTZ(i,j,3)  
  s22TXBYTZ(i,j,1)=s22TXBY(i,j,nz-2)    
  s23TXBY(i,j,nz-1)=s23TXBYTZ(i,j,2) 
  s23TXBY(i,j,nz)=s23TXBYTZ(i,j,3)  
  s23TXBYTZ(i,j,1)=s23TXBY(i,j,nz-2)    
  s33TXBY(i,j,nz-1)=s33TXBYTZ(i,j,2)   
  s33TXBY(i,j,nz)=s33TXBYTZ(i,j,3)   
  s33TXBYTZ(i,j,1)=s33TXBY(i,j,nz-2) 
 enddo 
enddo  
return 
end subroutine
! ###########################################################
Subroutine PasteStressTXTZ_TXBYTZ
include 'triffy.dec'
!
do k=1,npm
 do i=1,npm
  s11TXTZ(i,1,k)=s11TXBYTZ(i,npm-2,k)
  s11TXTZ(i,2,k)=s11TXBYTZ(i,npm-1,k)
  s11TXBYTZ(i,npm,k)=s11TXTZ(i,3,k)
  s12TXTZ(i,1,k)=s12TXBYTZ(i,npm-2,k)
  s12TXTZ(i,2,k)=s12TXBYTZ(i,npm-1,k)
  s12TXBYTZ(i,npm,k)=s12TXTZ(i,3,k)
  s13TXTZ(i,1,k)=s13TXBYTZ(i,npm-2,k)
  s13TXTZ(i,2,k)=s13TXBYTZ(i,npm-1,k)
  s13TXBYTZ(i,npm,k)=s13TXTZ(i,3,k)
  s22TXTZ(i,1,k)=s22TXBYTZ(i,npm-2,k) 
  s22TXTZ(i,2,k)=s22TXBYTZ(i,npm-1,k) 
  s22TXBYTZ(i,npm,k)=s22TXTZ(i,3,k)
  s23TXTZ(i,1,k)=s23TXBYTZ(i,npm-2,k) 
  s23TXTZ(i,2,k)=s23TXBYTZ(i,npm-1,k) 
  s23TXBYTZ(i,npm,k)=s23TXTZ(i,3,k) 
  s33TXTZ(i,1,k)=s33TXBYTZ(i,npm-2,k) 
  s33TXTZ(i,2,k)=s33TXBYTZ(i,npm-1,k)  
  s33TXBYTZ(i,npm,k)=s33TXTZ(i,3,k)  
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressBYTZ_TXBYTZ
include 'triffy.dec'
do k=1,npm
 do j=1,npm
  s11BYTZ(nx-1,j,k)=s11TXBYTZ(2,j,k)
  s11BYTZ(nx,j,k)=s11TXBYTZ(3,j,k)
  s11TXBYTZ(1,j,k)=s11BYTZ(nx-2,j,k)
  s12BYTZ(nx-1,j,k)=s12TXBYTZ(2,j,k) 
  s12BYTZ(nx,j,k)=s12TXBYTZ(3,j,k) 
  s12TXBYTZ(1,j,k)=s12BYTZ(nx-2,j,k)   
  s13BYTZ(nx-1,j,k)=s13TXBYTZ(2,j,k) 
  s13BYTZ(nx,j,k)=s13TXBYTZ(3,j,k) 
  s13TXBYTZ(1,j,k)=s13BYTZ(nx-2,j,k)   
  s22BYTZ(nx-1,j,k)=s22TXBYTZ(2,j,k)  
  s22BYTZ(nx,j,k)=s22TXBYTZ(3,j,k)  
  s22TXBYTZ(1,j,k)=s22BYTZ(nx-2,j,k)    
  s23BYTZ(nx-1,j,k)=s23TXBYTZ(2,j,k)  
  s23BYTZ(nx,j,k)=s23TXBYTZ(3,j,k)  
  s23TXBYTZ(1,j,k)=s23BYTZ(nx-2,j,k)    
  s33BYTZ(nx-1,j,k)=s33TXBYTZ(2,j,k)   
  s33BYTZ(nx,j,k)=s33TXBYTZ(3,j,k)   
  s33TXBYTZ(1,j,k)=s33BYTZ(nx-2,j,k) 
 enddo 
enddo  
return 
end subroutine   
! ###########################################################
Subroutine PasteStressTXBY_TXBYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  s11TXBY(i,j,1)=s11TXBYBZ(i,j,npm-2)
  s11TXBY(i,j,2)=s11TXBYBZ(i,j,npm-1)
  s11TXBYBZ(i,j,npm)=s11TXBY(i,j,3)
  s12TXBY(i,j,1)=s12TXBYBZ(i,j,npm-2) 
  s12TXBY(i,j,2)=s12TXBYBZ(i,j,npm-1) 
  s12TXBYBZ(i,j,npm)=s12TXBY(i,j,3) 
  s13TXBY(i,j,1)=s13TXBYBZ(i,j,npm-2) 
  s13TXBY(i,j,2)=s13TXBYBZ(i,j,npm-1) 
  s13TXBYBZ(i,j,npm)=s13TXBY(i,j,3) 
  s22TXBY(i,j,1)=s22TXBYBZ(i,j,npm-2) 
  s22TXBY(i,j,2)=s22TXBYBZ(i,j,npm-1)  
  s22TXBYBZ(i,j,npm)=s22TXBY(i,j,3)  
  s23TXBY(i,j,1)=s23TXBYBZ(i,j,npm-2)  
  s23TXBY(i,j,2)=s23TXBYBZ(i,j,npm-1)  
  s23TXBYBZ(i,j,npm)=s23TXBY(i,j,3)  
  s33TXBY(i,j,1)=s33TXBYBZ(i,j,npm-2)   
  s33TXBY(i,j,2)=s33TXBYBZ(i,j,npm-1)   
  s33TXBYBZ(i,j,npm)=s33TXBY(i,j,3)   
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressTXBZ_TXBYBZ
include 'triffy.dec'
!
do k=1,npm
 do i=1,npm
  s11TXBZ(i,1,k)=s11TXBYBZ(i,npm-2,k)
  s11TXBZ(i,2,k)=s11TXBYBZ(i,npm-1,k)
  s11TXBYBZ(i,npm,k)=s11TXBZ(i,3,k)
  s12TXBZ(i,1,k)=s12TXBYBZ(i,npm-2,k)
  s12TXBZ(i,2,k)=s12TXBYBZ(i,npm-1,k)
  s12TXBYBZ(i,npm,k)=s12TXBZ(i,3,k)
  s13TXBZ(i,1,k)=s13TXBYBZ(i,npm-2,k)
  s13TXBZ(i,2,k)=s13TXBYBZ(i,npm-1,k)
  s13TXBYBZ(i,npm,k)=s13TXBZ(i,3,k)
  s22TXBZ(i,1,k)=s22TXBYBZ(i,npm-2,k) 
  s22TXBZ(i,2,k)=s22TXBYBZ(i,npm-1,k) 
  s22TXBYBZ(i,npm,k)=s22TXBZ(i,3,k)
  s23TXBZ(i,1,k)=s23TXBYBZ(i,npm-2,k) 
  s23TXBZ(i,2,k)=s23TXBYBZ(i,npm-1,k) 
  s23TXBYBZ(i,npm,k)=s23TXBZ(i,3,k) 
  s33TXBZ(i,1,k)=s33TXBYBZ(i,npm-2,k) 
  s33TXBZ(i,2,k)=s33TXBYBZ(i,npm-1,k)  
  s33TXBYBZ(i,npm,k)=s33TXBZ(i,3,k)  
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressBYBZ_TXBYBZ
include 'triffy.dec'
do k=1,npm
 do j=1,npm
  s11BYBZ(nx-1,j,k)=s11TXBYBZ(2,j,k)
  s11BYBZ(nx,j,k)=s11TXBYBZ(3,j,k)
  s11TXBYBZ(1,j,k)=s11BYBZ(nx-2,j,k)
  s12BYBZ(nx-1,j,k)=s12TXBYBZ(2,j,k) 
  s12BYBZ(nx,j,k)=s12TXBYBZ(3,j,k) 
  s12TXBYBZ(1,j,k)=s12BYBZ(nx-2,j,k)   
  s13BYBZ(nx-1,j,k)=s13TXBYBZ(2,j,k) 
  s13BYBZ(nx,j,k)=s13TXBYBZ(3,j,k) 
  s13TXBYBZ(1,j,k)=s13BYBZ(nx-2,j,k)   
  s22BYBZ(nx-1,j,k)=s22TXBYBZ(2,j,k)  
  s22BYBZ(nx,j,k)=s22TXBYBZ(3,j,k)  
  s22TXBYBZ(1,j,k)=s22BYBZ(nx-2,j,k)    
  s23BYBZ(nx-1,j,k)=s23TXBYBZ(2,j,k)  
  s23BYBZ(nx,j,k)=s23TXBYBZ(3,j,k)  
  s23TXBYBZ(1,j,k)=s23BYBZ(nx-2,j,k)    
  s33BYBZ(nx-1,j,k)=s33TXBYBZ(2,j,k)   
  s33BYBZ(nx,j,k)=s33TXBYBZ(3,j,k)   
  s33TXBYBZ(1,j,k)=s33BYBZ(nx-2,j,k) 
 enddo 
enddo  
return 
end subroutine   
! ########################################################### 
Subroutine PasteStressBXTY_BXTYTZ
include 'triffy.dec'
do j=1,npm 
 do i=1,npm 
  s11BXTY(i,j,nz-1)=s11BXTYTZ(i,j,2) 
  s11BXTY(i,j,nz)=s11BXTYTZ(i,j,3)
  s11BXTYTZ(i,j,1)=s11BXTY(i,j,nz-2) 
  s12BXTY(i,j,nz-1)=s12BXTYTZ(i,j,2) 
  s12BXTY(i,j,nz)=s12BXTYTZ(i,j,3)
  s12BXTYTZ(i,j,1)=s12BXTY(i,j,nz-2) 
  s13BXTY(i,j,nz-1)=s13BXTYTZ(i,j,2) 
  s13BXTY(i,j,nz)=s13BXTYTZ(i,j,3)
  s13BXTYTZ(i,j,1)=s13BXTY(i,j,nz-2) 
  s22BXTY(i,j,nz-1)=s22BXTYTZ(i,j,2) 
  s22BXTY(i,j,nz)=s22BXTYTZ(i,j,3) 
  s22BXTYTZ(i,j,1)=s22BXTY(i,j,nz-2)  
  s23BXTY(i,j,nz-1)=s23BXTYTZ(i,j,2)  
  s23BXTY(i,j,nz)=s23BXTYTZ(i,j,3) 
  s23BXTYTZ(i,j,1)=s23BXTY(i,j,nz-2)  
  s33BXTY(i,j,nz-1)=s33BXTYTZ(i,j,2)   
  s33BXTY(i,j,nz)=s33BXTYTZ(i,j,3)  
  s33BXTYTZ(i,j,1)=s33BXTY(i,j,nz-2)
 enddo
enddo
return 
end subroutine 
! ########################################################### 
Subroutine PasteStressBXTZ_BXTYTZ 
include 'triffy.dec'
do k=1,npm
 do i=1,npm 
  s11BXTZ(i,ny-1,k)=s11BXTYTZ(i,2,k) 
  s11BXTZ(i,ny,k)=s11BXTYTZ(i,3,k) 
  s11BXTYTZ(i,1,k)=s11BXTZ(i,ny-2,k) 
  s12BXTZ(i,ny-1,k)=s12BXTYTZ(i,2,k)  
  s12BXTZ(i,ny,k)=s12BXTYTZ(i,3,k)  
  s12BXTYTZ(i,1,k)=s12BXTZ(i,ny-2,k)   
  s13BXTZ(i,ny-1,k)=s13BXTYTZ(i,2,k)  
  s13BXTZ(i,ny,k)=s13BXTYTZ(i,3,k)  
  s13BXTYTZ(i,1,k)=s13BXTZ(i,ny-2,k)   
  s22BXTZ(i,ny-1,k)=s22BXTYTZ(i,2,k)  
  s22BXTZ(i,ny,k)=s22BXTYTZ(i,3,k)  
  s22BXTYTZ(i,1,k)=s22BXTZ(i,ny-2,k)    
  s23BXTZ(i,ny-1,k)=s23BXTYTZ(i,2,k)   
  s23BXTZ(i,ny,k)=s23BXTYTZ(i,3,k)   
  s23BXTYTZ(i,1,k)=s23BXTZ(i,ny-2,k)   
  s33BXTZ(i,ny-1,k)=s33BXTYTZ(i,2,k)    
  s33BXTZ(i,ny,k)=s33BXTYTZ(i,3,k)    
  s33BXTYTZ(i,1,k)=s33BXTZ(i,ny-2,k)    
 enddo
enddo 
return
end subroutine 
! ######################################
Subroutine PasteStressTYTZ_BXTYTZ
include 'triffy.dec'
 
do k=1,npm
 do j=1,npm
  s11TYTZ(1,j,k)=s11BXTYTZ(npm-2,j,k)
  s11TYTZ(2,j,k)=s11BXTYTZ(npm-1,j,k)
  s11BXTYTZ(npm,j,k)=s11TYTZ(3,j,k) 
  s12TYTZ(1,j,k)=s12BXTYTZ(npm-2,j,k)
  s12TYTZ(2,j,k)=s12BXTYTZ(npm-1,j,k)
  s12BXTYTZ(npm,j,k)=s12TYTZ(3,j,k) 
  s13TYTZ(1,j,k)=s13BXTYTZ(npm-2,j,k)
  s13TYTZ(2,j,k)=s13BXTYTZ(npm-1,j,k)
  s13BXTYTZ(npm,j,k)=s13TYTZ(3,j,k) 
  s22TYTZ(1,j,k)=s22BXTYTZ(npm-2,j,k) 
  s22TYTZ(2,j,k)=s22BXTYTZ(npm-1,j,k) 
  s22BXTYTZ(npm,j,k)=s22TYTZ(3,j,k) 
  s23TYTZ(1,j,k)=s23BXTYTZ(npm-2,j,k)
  s23TYTZ(2,j,k)=s23BXTYTZ(npm-1,j,k) 
  s23BXTYTZ(npm,j,k)=s23TYTZ(3,j,k)  
  s33TYTZ(1,j,k)=s33BXTYTZ(npm-2,j,k) 
  s33TYTZ(2,j,k)=s33BXTYTZ(npm-1,j,k)  
  s33BXTYTZ(npm,j,k)=s33TYTZ(3,j,k)   
 enddo 
enddo 
return 
end subroutine   
! ########################################################### 
Subroutine PasteStressBXTY_BXTYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  s11BXTY(i,j,1)=s11BXTYBZ(i,j,npm-2)
  s11BXTY(i,j,2)=s11BXTYBZ(i,j,npm-1)
  s11BXTYBZ(i,j,npm)=s11BXTY(i,j,3)
  s12BXTY(i,j,1)=s12BXTYBZ(i,j,npm-2) 
  s12BXTY(i,j,2)=s12BXTYBZ(i,j,npm-1) 
  s12BXTYBZ(i,j,npm)=s12BXTY(i,j,3) 
  s13BXTY(i,j,1)=s13BXTYBZ(i,j,npm-2) 
  s13BXTY(i,j,2)=s13BXTYBZ(i,j,npm-1) 
  s13BXTYBZ(i,j,npm)=s13BXTY(i,j,3) 
  s22BXTY(i,j,1)=s22BXTYBZ(i,j,npm-2) 
  s22BXTY(i,j,2)=s22BXTYBZ(i,j,npm-1)  
  s22BXTYBZ(i,j,npm)=s22BXTY(i,j,3)  
  s23BXTY(i,j,1)=s23BXTYBZ(i,j,npm-2)  
  s23BXTY(i,j,2)=s23BXTYBZ(i,j,npm-1)  
  s23BXTYBZ(i,j,npm)=s23BXTY(i,j,3)  
  s33BXTY(i,j,1)=s33BXTYBZ(i,j,npm-2)   
  s33BXTY(i,j,2)=s33BXTYBZ(i,j,npm-1)   
  s33BXTYBZ(i,j,npm)=s33BXTY(i,j,3)   
 enddo
enddo
return
end subroutine
! ########################################################### 
Subroutine PasteStressBXBZ_BXTYBZ 
include 'triffy.dec'
do k=1,npm
 do i=1,npm 
  s11BXBZ(i,ny-1,k)=s11BXTYBZ(i,2,k) 
  s11BXBZ(i,ny,k)=s11BXTYBZ(i,3,k) 
  s11BXTYBZ(i,1,k)=s11BXBZ(i,ny-2,k) 
  s12BXBZ(i,ny-1,k)=s12BXTYBZ(i,2,k)  
  s12BXBZ(i,ny,k)=s12BXTYBZ(i,3,k)  
  s12BXTYBZ(i,1,k)=s12BXBZ(i,ny-2,k)   
  s13BXBZ(i,ny-1,k)=s13BXTYBZ(i,2,k)  
  s13BXBZ(i,ny,k)=s13BXTYBZ(i,3,k)  
  s13BXTYBZ(i,1,k)=s13BXBZ(i,ny-2,k)   
  s22BXBZ(i,ny-1,k)=s22BXTYBZ(i,2,k)  
  s22BXBZ(i,ny,k)=s22BXTYBZ(i,3,k)  
  s22BXTYBZ(i,1,k)=s22BXBZ(i,ny-2,k)    
  s23BXBZ(i,ny-1,k)=s23BXTYBZ(i,2,k)   
  s23BXBZ(i,ny,k)=s23BXTYBZ(i,3,k)   
  s23BXTYBZ(i,1,k)=s23BXBZ(i,ny-2,k)   
  s33BXBZ(i,ny-1,k)=s33BXTYBZ(i,2,k)    
  s33BXBZ(i,ny,k)=s33BXTYBZ(i,3,k)    
  s33BXTYBZ(i,1,k)=s33BXBZ(i,ny-2,k)    
 enddo
enddo 
return
end subroutine 
! ######################################
Subroutine PasteStressTYBZ_BXTYBZ
include 'triffy.dec'
 
do k=1,npm
 do j=1,npm
  s11TYBZ(1,j,k)=s11BXTYBZ(npm-2,j,k)
  s11TYBZ(2,j,k)=s11BXTYBZ(npm-1,j,k)
  s11BXTYBZ(npm,j,k)=s11TYBZ(3,j,k) 
  s12TYBZ(1,j,k)=s12BXTYBZ(npm-2,j,k)
  s12TYBZ(2,j,k)=s12BXTYBZ(npm-1,j,k)
  s12BXTYBZ(npm,j,k)=s12TYBZ(3,j,k) 
  s13TYBZ(1,j,k)=s13BXTYBZ(npm-2,j,k)
  s13TYBZ(2,j,k)=s13BXTYBZ(npm-1,j,k)
  s13BXTYBZ(npm,j,k)=s13TYBZ(3,j,k) 
  s22TYBZ(1,j,k)=s22BXTYBZ(npm-2,j,k) 
  s22TYBZ(2,j,k)=s22BXTYBZ(npm-1,j,k) 
  s22BXTYBZ(npm,j,k)=s22TYBZ(3,j,k) 
  s23TYBZ(1,j,k)=s23BXTYBZ(npm-2,j,k)
  s23TYBZ(2,j,k)=s23BXTYBZ(npm-1,j,k) 
  s23BXTYBZ(npm,j,k)=s23TYBZ(3,j,k)  
  s33TYBZ(1,j,k)=s33BXTYBZ(npm-2,j,k) 
  s33TYBZ(2,j,k)=s33BXTYBZ(npm-1,j,k)  
  s33BXTYBZ(npm,j,k)=s33TYBZ(3,j,k)   
 enddo 
enddo 
return 
end subroutine   
! ###########################################################
Subroutine PasteStressBXBY_BXBYTZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  s11BXBY(i,j,nz-1)=s11BXBYTZ(i,j,2)
  s11BXBY(i,j,nz)=s11BXBYTZ(i,j,3) 
  s11BXBYTZ(i,j,1)=s11BXBY(i,j,nz-2) 
  s12BXBY(i,j,nz-1)=s12BXBYTZ(i,j,2)
  s12BXBY(i,j,nz)=s12BXBYTZ(i,j,3) 
  s12BXBYTZ(i,j,1)=s12BXBY(i,j,nz-2)   
  s13BXBY(i,j,nz-1)=s13BXBYTZ(i,j,2)
  s13BXBY(i,j,nz)=s13BXBYTZ(i,j,3) 
  s13BXBYTZ(i,j,1)=s13BXBY(i,j,nz-2)   
  s22BXBY(i,j,nz-1)=s22BXBYTZ(i,j,2)
  s22BXBY(i,j,nz)=s22BXBYTZ(i,j,3)  
  s22BXBYTZ(i,j,1)=s22BXBY(i,j,nz-2)    
  s23BXBY(i,j,nz-1)=s23BXBYTZ(i,j,2) 
  s23BXBY(i,j,nz)=s23BXBYTZ(i,j,3)  
  s23BXBYTZ(i,j,1)=s23BXBY(i,j,nz-2)    
  s33BXBY(i,j,nz-1)=s33BXBYTZ(i,j,2)   
  s33BXBY(i,j,nz)=s33BXBYTZ(i,j,3)   
  s33BXBYTZ(i,j,1)=s33BXBY(i,j,nz-2) 
 enddo 
enddo  
return 
end subroutine
! ###########################################################
Subroutine PasteStressBXTZ_BXBYTZ
include 'triffy.dec'
!
do k=1,npm
 do i=1,npm
  s11BXTZ(i,1,k)=s11BXBYTZ(i,npm-2,k)
  s11BXTZ(i,2,k)=s11BXBYTZ(i,npm-1,k)
  s11BXBYTZ(i,npm,k)=s11BXTZ(i,3,k)
  s12BXTZ(i,1,k)=s12BXBYTZ(i,npm-2,k)
  s12BXTZ(i,2,k)=s12BXBYTZ(i,npm-1,k)
  s12BXBYTZ(i,npm,k)=s12BXTZ(i,3,k)
  s13BXTZ(i,1,k)=s13BXBYTZ(i,npm-2,k)
  s13BXTZ(i,2,k)=s13BXBYTZ(i,npm-1,k)
  s13BXBYTZ(i,npm,k)=s13BXTZ(i,3,k)
  s22BXTZ(i,1,k)=s22BXBYTZ(i,npm-2,k) 
  s22BXTZ(i,2,k)=s22BXBYTZ(i,npm-1,k) 
  s22BXBYTZ(i,npm,k)=s22BXTZ(i,3,k)
  s23BXTZ(i,1,k)=s23BXBYTZ(i,npm-2,k) 
  s23BXTZ(i,2,k)=s23BXBYTZ(i,npm-1,k) 
  s23BXBYTZ(i,npm,k)=s23BXTZ(i,3,k) 
  s33BXTZ(i,1,k)=s33BXBYTZ(i,npm-2,k) 
  s33BXTZ(i,2,k)=s33BXBYTZ(i,npm-1,k)  
  s33BXBYTZ(i,npm,k)=s33BXTZ(i,3,k)  
 enddo
enddo
return
end subroutine
! ######################################
Subroutine PasteStressBYTZ_BXBYTZ
include 'triffy.dec'
 
do k=1,npm
 do j=1,npm
  s11BYTZ(1,j,k)=s11BXBYTZ(npm-2,j,k)
  s11BYTZ(2,j,k)=s11BXBYTZ(npm-1,j,k)
  s11BXBYTZ(npm,j,k)=s11BYTZ(3,j,k) 
  s12BYTZ(1,j,k)=s12BXBYTZ(npm-2,j,k)
  s12BYTZ(2,j,k)=s12BXBYTZ(npm-1,j,k)
  s12BXBYTZ(npm,j,k)=s12BYTZ(3,j,k) 
  s13BYTZ(1,j,k)=s13BXBYTZ(npm-2,j,k)
  s13BYTZ(2,j,k)=s13BXBYTZ(npm-1,j,k)
  s13BXBYTZ(npm,j,k)=s13BYTZ(3,j,k) 
  s22BYTZ(1,j,k)=s22BXBYTZ(npm-2,j,k) 
  s22BYTZ(2,j,k)=s22BXBYTZ(npm-1,j,k) 
  s22BXBYTZ(npm,j,k)=s22BYTZ(3,j,k) 
  s23BYTZ(1,j,k)=s23BXBYTZ(npm-2,j,k)
  s23BYTZ(2,j,k)=s23BXBYTZ(npm-1,j,k) 
  s23BXBYTZ(npm,j,k)=s23BYTZ(3,j,k)  
  s33BYTZ(1,j,k)=s33BXBYTZ(npm-2,j,k) 
  s33BYTZ(2,j,k)=s33BXBYTZ(npm-1,j,k)  
  s33BXBYTZ(npm,j,k)=s33BYTZ(3,j,k)   
 enddo 
enddo 
return 
end subroutine   
! ###########################################################
Subroutine PasteStressBXBY_BXBYBZ
include 'triffy.dec'
 
do j=1,npm
 do i=1,npm
  s11BXBY(i,j,1)=s11BXBYBZ(i,j,npm-2)
  s11BXBY(i,j,2)=s11BXBYBZ(i,j,npm-1)
  s11BXBYBZ(i,j,npm)=s11BXBY(i,j,3)
  s12BXBY(i,j,1)=s12BXBYBZ(i,j,npm-2) 
  s12BXBY(i,j,2)=s12BXBYBZ(i,j,npm-1) 
  s12BXBYBZ(i,j,npm)=s12BXBY(i,j,3) 
  s13BXBY(i,j,1)=s13BXBYBZ(i,j,npm-2) 
  s13BXBY(i,j,2)=s13BXBYBZ(i,j,npm-1) 
  s13BXBYBZ(i,j,npm)=s13BXBY(i,j,3) 
  s22BXBY(i,j,1)=s22BXBYBZ(i,j,npm-2) 
  s22BXBY(i,j,2)=s22BXBYBZ(i,j,npm-1)  
  s22BXBYBZ(i,j,npm)=s22BXBY(i,j,3)  
  s23BXBY(i,j,1)=s23BXBYBZ(i,j,npm-2)  
  s23BXBY(i,j,2)=s23BXBYBZ(i,j,npm-1)  
  s23BXBYBZ(i,j,npm)=s23BXBY(i,j,3)  
  s33BXBY(i,j,1)=s33BXBYBZ(i,j,npm-2)   
  s33BXBY(i,j,2)=s33BXBYBZ(i,j,npm-1)   
  s33BXBYBZ(i,j,npm)=s33BXBY(i,j,3)   
 enddo
enddo
return
end subroutine
! ###########################################################
Subroutine PasteStressBXBZ_BXBYBZ
include 'triffy.dec'
!
do k=1,npm
 do i=1,npm
  s11BXBZ(i,1,k)=s11BXBYBZ(i,npm-2,k)
  s11BXBZ(i,2,k)=s11BXBYBZ(i,npm-1,k)
  s11BXBYBZ(i,npm,k)=s11BXBZ(i,3,k)
  s12BXBZ(i,1,k)=s12BXBYBZ(i,npm-2,k)
  s12BXBZ(i,2,k)=s12BXBYBZ(i,npm-1,k)
  s12BXBYBZ(i,npm,k)=s12BXBZ(i,3,k)
  s13BXBZ(i,1,k)=s13BXBYBZ(i,npm-2,k)
  s13BXBZ(i,2,k)=s13BXBYBZ(i,npm-1,k)
  s13BXBYBZ(i,npm,k)=s13BXBZ(i,3,k)
  s22BXBZ(i,1,k)=s22BXBYBZ(i,npm-2,k) 
  s22BXBZ(i,2,k)=s22BXBYBZ(i,npm-1,k) 
  s22BXBYBZ(i,npm,k)=s22BXBZ(i,3,k)
  s23BXBZ(i,1,k)=s23BXBYBZ(i,npm-2,k) 
  s23BXBZ(i,2,k)=s23BXBYBZ(i,npm-1,k) 
  s23BXBYBZ(i,npm,k)=s23BXBZ(i,3,k) 
  s33BXBZ(i,1,k)=s33BXBYBZ(i,npm-2,k) 
  s33BXBZ(i,2,k)=s33BXBYBZ(i,npm-1,k)  
  s33BXBYBZ(i,npm,k)=s33BXBZ(i,3,k)  
 enddo
enddo
return
end subroutine
! ######################################
Subroutine PasteStressBYBZ_BXBYBZ
include 'triffy.dec'
 
do k=1,npm
 do j=1,npm
  s11BYBZ(1,j,k)=s11BXBYBZ(npm-2,j,k)
  s11BYBZ(2,j,k)=s11BXBYBZ(npm-1,j,k)
  s11BXBYBZ(npm,j,k)=s11BYBZ(3,j,k) 
  s12BYBZ(1,j,k)=s12BXBYBZ(npm-2,j,k)
  s12BYBZ(2,j,k)=s12BXBYBZ(npm-1,j,k)
  s12BXBYBZ(npm,j,k)=s12BYBZ(3,j,k) 
  s13BYBZ(1,j,k)=s13BXBYBZ(npm-2,j,k)
  s13BYBZ(2,j,k)=s13BXBYBZ(npm-1,j,k)
  s13BXBYBZ(npm,j,k)=s13BYBZ(3,j,k) 
  s22BYBZ(1,j,k)=s22BXBYBZ(npm-2,j,k) 
  s22BYBZ(2,j,k)=s22BXBYBZ(npm-1,j,k) 
  s22BXBYBZ(npm,j,k)=s22BYBZ(3,j,k) 
  s23BYBZ(1,j,k)=s23BXBYBZ(npm-2,j,k)
  s23BYBZ(2,j,k)=s23BXBYBZ(npm-1,j,k) 
  s23BXBYBZ(npm,j,k)=s23BYBZ(3,j,k)  
  s33BYBZ(1,j,k)=s33BXBYBZ(npm-2,j,k) 
  s33BYBZ(2,j,k)=s33BXBYBZ(npm-1,j,k)  
  s33BXBYBZ(npm,j,k)=s33BYBZ(3,j,k)   
 enddo 
enddo 
return 
end subroutine   

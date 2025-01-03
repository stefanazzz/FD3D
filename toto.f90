       open (11,file='fpar',form='formatted',status='old')
       read (11,'(a)')
       read (11,*) dx
       read (11,*) dt
       read (11,*) ntime
       read (11,'(a)')
       read (11,*) tw
       read (11,*) vforce
       read (11,*) peak
       read (11,*) ssini
       read (11,*) delta
       read (11,*) vmax
       read (11,*) dissipation
       read (11,*) initradius
       read (11,*) mu_s
       read (11,*) mu_d
       read (11,'(a)')
       read (11,*) nxs,nys
       read (11,'(a)')
       nr=1
       do while (0.eq.0)
        read (11,*,err=30) nxr,nyr
        nr=nr+1
       enddo
30     close (11)
       print*, nr
       end

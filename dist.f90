function dist(xx,lama,mua,rhoa)
       include 'triffy.dec'
       real Rpp,xx,lama,mua,rhoa
       parameter (Rpp=1.e-3)

       vp=sqrt((2*lama+mua)/rhoa)
       dist=3*3.*vp/2./npm/dx
       dist=dist*(xx/npm)**6
       return
       end

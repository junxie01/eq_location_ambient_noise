       subroutine cal_dist(lon1,lat1,lon2,lat2,dist)
       real lon1,lat1,lon2,lat2
       real dist,cva,cvb,cv
       real l,tpdist,u1,u2,numda,pi,numda1
       real cv1,cv2,cv3,cv4,cv5,cvc,mius,cv6,cv7,deltacv
       real az,evla,evlo,stla,stlo,delta_o,sinp,cosp,cosaz
       cva=6378.137
       cvb= 6356.7523142
       f = 1/298.257223563
       !pi=atan(1)*4.0
       pi=3.14159265
       L = lon1-lon2
       if (L > 180.000)  L =360.000000 - L
       if (L < -180.000) L =  360.000 - abs(L)
       L = abs(L); 
 
       if(lat1 .eq.lat2 .and. lon1.eq.lon2 ) then
             dist = 0.0
             return
       endif
       U1 = 0
       U1 = atan((1-f)*tan(lat1/180*pi))+ 0.001
       U2 = 0
       U2 = atan((1-f)*tan(lat2/180*pi))+ 0.001

       L = L*pi/180
       numda=L
       numda1 = numda
       do  
          numda = numda1;
          cv1 =  sqrt( (cos(U2)*sin(numda))*(cos(U2)*sin(numda))
     1        + (cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda))
     1        *(cos(U1)*sin(U2)-sin(U1)*cos(U2)*cos(numda)) )
          cv2 = sin(U1)*sin(U2)+ cos(U1)*cos(U2)*cos(numda)   
          cv = atan2(cv1,cv2)
          cv3 = cos(U1)*cos(U2)*sin(numda)/sin(cv)   
          cv4 = 1 - cv3*cv3
          cv5 = cos(cv) - 2*sin(U1)*sin(U2)/cv4
          cvC = f/16*cv4*(4 + f*(4 - 3*cv4));
         numda1 = L + (1-cvC)*f*cv3*(cv + cvC*cv1*
     1           (cv5 + cvC*cv2*(-1 +2*cv5*cv5)))
        if( abs(numda - numda1) < 0.0000001) goto 100
       enddo

100    mius = cv4*(cva*cva - cvb*cvb)/(cvb*cvb)
       cv6 = 1+mius/16384*(4096 + mius*(-768 + mius*(320 - 175*mius)))
       cv7 = mius/1024*(256+ mius*(-128 + mius*(74 - 47*mius)))
       deltacv = cv7*cv1*(cv5 +cv7/4*(cv2*(-1 + 2*cv5*cv5)
     1         -cv7/6*cv5*(-3+4*cv1*cv1)*(-3+4*cv5*cv5) ))

      dist = cvb * cv6 *(cv - deltacv)
!     calculate the azimuth 
!      evlo=lon1*pi/180 
!      stlo=lon2*pi/180
!      evla=(90-lat1)*pi/180
!      stla=(90-lat2)*pi/180
!      delta_o=stlo-evlo
!      cosp=cos(evla)*cos(stla)+sin(evla)*sin(stla)*cos(delta_o)
!      sinp=sqrt(1-cosp**2)
!      cosaz=(cos(stla)-cosp*cos(evla))/sinp/sin(evla)
!      az=acos(cosaz)*180/pi
!      if(delta_o<0)az=360-az
      return
      end subroutine

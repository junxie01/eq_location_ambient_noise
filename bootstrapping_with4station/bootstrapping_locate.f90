! program to locate "earthquake" using grid search
! 2019/07/09 fix t00 and use array
! change to bootstrapping
! 2019/07/27 every quadrant has one random station
! 2019/07/28 use another way
! 2019/07/30
program bootstrapping
integer,parameter :: maxnp=100,maxnsta=100
real,dimension(maxnp,maxnsta) :: period_reference,disp_reference,period_target,disp_target
real,dimension(maxnsta) :: stla,stlo
real :: stlaa,stloo,disp,stla_reference,stlo_reference
real :: stlo_b,stla_b,stlo_e,stla_e,stplo,stpla,tt,t_step,begin_t,end_t,t,tt1,tt2
real :: stlo_g,stla_g,azall(maxnsta),t0(maxnsta),t00,error(maxnsta),error_all
real :: vest(maxnp,maxnsta),az
real :: dist(maxnsta),dist0(maxnsta),distg(maxnsta)
real :: errorabs(maxnsta),errorabsall
real :: llo,lla,tt0,ddt,ddist,minerror
real :: num
real :: az_range(4,2)
integer :: azid(4,maxnsta),numaz(4),id
integer :: iis,qmark(4),iaz,dmark,idsta(4),idds,picked_idsta(4)
character(7)  :: sta_reference,sta_target
character(20) :: filen1,filen2,station(maxnsta)
character(80) :: file_reference(maxnsta),file_target(maxnsta)
character(80) :: para,sta_list,list_reference,list_target,output,ref,target
integer i,j,n,k,np,nsta,nstp_la,nstp_lo,nref,ntarget,nstp_t
integer ilo,ila,ip,if,is,ip1,ip2,it,imark,ii
integer no_p_reference(maxnsta),no_p_target(maxnsta)
integer nboot
logical ext1,ext2
if(iargc().ne.1)then
   write(*,*)'Usage: locat locate.par'
   write(*,*)'locate.par:  '
   write(*,*)'station.list '
   write(*,*)'ref_sta'
   write(*,*)'stlo_g,stla_g'
   write(*,*)'output file' 
   write(*,*)'stlo_b stlo_e steplo'
   write(*,*)'stla_b stla_e stepla'
   call exit(-1)
endif
call getarg(1,para)
open(10,file=para)
read(10,*) sta_list
read(10,*) sta_reference,stlo_reference,stla_reference
read(10,*) stlo_g,stla_g
read(10,*) output
read(10,*) stlo_b,stlo_e,stplo
read(10,*) stla_b,stla_e,stpla
read(10,*) nboot
close(10)
az_range(1,1)=0
az_range(1,2)=90
az_range(2,1)=90
az_range(2,2)=180
az_range(3,1)=180
az_range(3,2)=270
az_range(4,1)=270
az_range(4,2)=360
no_p_reference=0
no_p_target=0
is=1
write(*,'("from",f10.4,"to",f10.4,", from",f10.4,"to",f10.4,",")')stlo_b,stlo_e,stla_b,stla_e
write(*,'("taking station",1x,1a,1x,"as reference.")')sta_reference
! read the station list
! the last one is the one to be found 
! and here it is the reference location assume it is the first guess location
open(7,file=sta_list)
numaz=0;azid=0
do is=1,maxnsta
   read(7,*,end=20,err=20)station(is), stlo(is),stla(is)
   call cal_dist(stlo(is),stla(is),stlo_reference,stla_reference,dist0(is),azall(is))    ! distance between tele station and the reference location
   call cal_dist(stlo(is),stla(is),stlo_g,stla_g,distg(is),az)    ! distance between tele station and the first guess location
   do iaz=1,4
      if(az_range(iaz,1).le.azall(is).and.az_range(iaz,2).gt.azall(is))then
         numaz(iaz)=numaz(iaz)+1          ! number of stations at quadrant iaz
         azid(iaz,numaz(iaz))=is          ! id of the station at quadrant iaz
      endif
   enddo
enddo
20 close(7)
nsta=is-1                                 ! number of station
nstp_lo=int((stlo_e-stlo_b)/stplo)+1     ! number of longitude step
nstp_la=int((stla_e-stla_b)/stpla)+1     ! number of lanitude step
! read in the dispersion curve
do is=1,nsta                             
   ref=trim(sta_reference)//'_'//trim(station(is))//'.dsp'
   ref='dsp_ncf_'//trim(station(is))//'.dat'
   target=trim(sta_target)//'_'//trim(station(is))//'.dsp'
   target='dsp_eqe_'//trim(station(is))//'.dat'
   inquire(file=ref,exist=ext1)
   inquire(file=target,exist=ext2)
   if(ext1.and.ext2)then
      !write(*,'("Read in the reference dispersion curve",1x,1a)')trim(ref)
      ! read in reference group velocity disperison curve 
      ip=1
      open(10,file=ref) 
 12   read(10,*,err=11,end=11)period_reference(ip,is),disp_reference(ip,is)
      ip=ip+1
      goto 12
 11   continue
      no_p_reference(is)=ip-1
      call sort(period_reference(:,is),disp_reference(:,is),no_p_reference(is),1)
      ! read in earthquake group velocity disperison curve 
      ip=1
      open(13,file=target) 
      !write(*,'("Read in dispersion curve",1x,1a)')trim(target)
  15  read(13,*,err=14,end=14)period_target(ip,is),disp_target(ip,is)
      ip=ip+1
      goto 15
 14   continue
      no_p_target(is)=ip-1
      call sort(period_target(:,is),disp_target(:,is),no_p_target(is),1)
   endif
enddo
call init_random_seed(i)
open(20,file=output)
do ibt=1,nboot                                              ! begin bootstrap
   llo=0;lla=0;tt0=0;ddt=0;ddist=0;minerror=10000
   if(mod(ibt,nboot/10).eq.0)write(*,*)'Bootstrapping completed: ',real(ibt*10)/real(nboot)*10,"%"
   ! Begin grid search
   do i=1,nstp_lo                                           ! loop over lontitude
      stloo=stlo_b+(i-1)*stplo
      do j=1,nstp_la                                        ! loop over lantitude
         stlaa=stla_b+(j-1)*stpla
         t00=0;imark=0;vest=0;idsta=0
         do iaz=1,4                                         ! only four stations are used
            call random_number(num)
            id=int(num*numaz(iaz))+1                        ! id of the station at quadrant iaz
            is=azid(iaz,id)                                 ! is is the station id
            idsta(iaz)=is                                   ! restore the station id
            t0(iaz)=0;dist(iaz)=0
         !   write(*,*)'quadrant:',iaz,'id=',id,'is=',is,station(is)
            if(no_p_target(is).eq.0.or.no_p_reference(is).eq.0)cycle
            call cal_dist(stloo,stlaa,stlo(is),stla(is),dist(iaz),az)      ! distance between tele station and the grid
            do ip1=1,no_p_target(is)                      ! loop over period of "earthquake"
               do ip2=2,no_p_reference(is)                ! loop over period of reference dispersion curve
                  if(period_reference(ip2-1,is)<=period_target(ip1,is)&
                  .and.period_reference(ip2,is)>period_target(ip1,is))then
                     vest(ip1,iaz)=disp_reference(ip2-1,is)+(disp_reference(ip2,is)-&
                     disp_reference(ip2-1,is))*(period_target(ip1,is)-period_reference(ip2-1,is))&
                     /(period_reference(ip2,is)-period_reference(ip2-1,is))
                     imark=imark+1
                     exit
                  endif
               enddo                                      ! end loop over period of reference dispersion curve
               if(ip2.ne.no_p_reference(is)+1)then
                  tt=distg(is)/disp_target(ip1,is)        ! travel time of the earthquake signal
                  t0(iaz)=t0(iaz)+tt-dist(iaz)/vest(ip1,iaz)
               endif
            enddo                                         ! end loop over period of "earthquake"
            t00=t00+t0(iaz)
         enddo                                            ! end loop over stations
         t00=t00/imark                                    ! the original time
         error=0;error_all=0;errorabs=0;errorabsall=0
         !write(*,*)'t00=',t00
         do iaz=1,4                                        ! loop over tele stations
            is=idsta(iaz)
            if(no_p_target(is).eq.0.or.no_p_reference(is).eq.0)cycle
            do ip=1,no_p_target(is)                     ! loop over period of "earthquake"
               if(abs(vest(ip,iaz)).lt.1e-7)cycle
               tt=distg(is)/disp_target(ip,is)        ! travel time of the earthquake signal
               error(iaz)=error(iaz)+abs(tt-dist(iaz)/vest(ip,iaz)-t00)**2  !error1=error1+abs(dist/tt-v)**2
               errorabs(iaz)=errorabs(iaz)+abs(tt-dist(iaz)/vest(ip,iaz)-t00)
            enddo                                         ! end loop over period of earthquake
            error_all=error_all+error(iaz)/imark
            errorabsall=errorabsall+errorabs(iaz)/imark
         enddo                                            ! end loop over stations
         error_all=sqrt(error_all)                        ! travel time error
         !write(*,*)'errorall=',error_all
         if(error_all.lt.minerror)then
            llo=stloo;lla=stlaa;minerror=error_all;tt0=t00;ddt=errorabsall;picked_idsta=idsta
            call cal_dist(llo,lla,stlo_g,stla_g,ddist,az) ! location variation
         endif
         !write(*,'(6f15.5,1x,4i2.2)') llo,lla,minerror,tt0,ddt,ddist,(idsta(ii),ii=1,4)
      enddo                                               ! end loop over lantitude
   enddo                                                  ! end loop over longtitude
   !stop
   write(20,'(6f15.5,1x,4i2.2)') llo,lla,minerror,tt0,ddt,ddist,(picked_idsta(ii),ii=1,4)
enddo                                                     ! end loop over bootstrapping
close(20)
end program
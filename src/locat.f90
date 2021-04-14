! program to locate "earthquake" using grid search
! 2019/07/09 fix t00 and use array
! 2019/09/16 fix the input style
program main
integer,parameter :: maxnp=100,maxnsta=100
real,dimension(maxnp,maxnsta) :: period_reference,disp_reference,period_target,disp_target
real,dimension(maxnsta) :: stla,stlo
real :: stlaa,stloo,disp,stla_reference,stlo_reference
real :: stlo_b,stla_b,stlo_e,stla_e,stplo,stpla,tt,t_step,begin_t,end_t,t,tt1,tt2
real :: stlo_g,stla_g,az,t0(maxnsta),t00,error(maxnsta),error_all
real :: vest(maxnp,maxnsta)
real :: dist(maxnsta),dist0(maxnsta)
real :: errorabs,tmp
character(7)  :: sta_reference,sta_target
character(20) :: filen1,filen2,id,station(maxnsta)
character(80) :: file_reference(maxnsta),file_target(maxnsta)
character(80) :: para,sta_list,list_reference,list_target,output,ref,target
integer i,j,n,k,np,nsta,nstp_la,nstp_lo,nref,ntarget,nstp_t
integer ilo,ila,ip,if,is,ip1,ip2,it,imark,ista
integer no_p_reference(maxnsta),no_p_target(maxnsta)
logical ext1,ext2
if(iargc().ne.1)then
   write(*,*)'Usage: locat locate.par'
   write(*,*)'locate.par:  '
   write(*,*)'station.list '
   write(*,*)'ref_sta ref_stlo ref_stla'
   write(*,*)'stlo_g,stla_g'
   write(*,*)'output file' 
   write(*,*)'stlo_b stlo_e steplo'
   write(*,*)'stla_b stla_e stepla'
   write(*,*)'begin_t end_t step_t'
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
close(10)
no_p_reference=0
no_p_target=0
is=1
write(*,'("from",f10.4,"to",f10.4,", from",f10.4,"to",f10.4,",")')stlo_b,stlo_e,stla_b,stla_e
write(*,'("taking station",1x,1a,1x,"as reference.")')sta_reference
! read the station list
! the last one is the one to be found 
! and here it is the reference location assume it is the first guess location
open(7,file=sta_list)
do is=1,maxnsta
   read(7,*,end=20,err=20)station(is), stlo(is),stla(is)
enddo
20 close(7)
nsta=is-1                                ! number of station
nstp_lo=int((stlo_e-stlo_b)/stplo)+1     ! number of longitude step
nstp_la=int((stla_e-stla_b)/stpla)+1     ! number of lanitude step
!write(*,*)nstp_lo,nstp_la,stlo_b,stla_b,stlo_e,stla_e,stplo,stpla
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
! Begin grid search
open(20,file=output)
do i=1,nstp_lo                                                       ! loop over lontitude
   stloo=stlo_b+(i-1)*stplo
   do j=1,nstp_la                                                    ! loop over lantitude
      stlaa=stla_b+(j-1)*stpla
      ista=0;t00=0;imark=0;vest=0
      do is=1,nsta                                                   ! loop over tele stations
         t0(is)=0;dist(is)=0;dist0(is)=0
         if(no_p_target(is).eq.0.or.no_p_reference(is).eq.0)cycle
         write(id,'("File id: ",i5)')is
         call cal_dist(stloo,stlaa,stlo(is),stla(is),dist(is),tmp)       ! distance between tele station and the grid
         call cal_dist(stlo(is),stla(is),stlo_g,stla_g,dist0(is),tmp)    ! distance between tele station and the first guess location
         do ip1=1,no_p_target(is)                      ! loop over period of "earthquake"
            do ip2=2,no_p_reference(is)                ! loop over period of reference dispersion curve
               if(period_reference(ip2-1,is)<=period_target(ip1,is)&
                  !.and.period_reference(ip2,is)>period_target(ip1,is)&
                  .and.period_reference(ip2,is)>period_target(ip1,is))then
                  !.and.period_target(ip1,is).ge.10.and.period_target(ip1,is).le.20)then
                  vest(ip1,is)=disp_reference(ip2-1,is)+(disp_reference(ip2,is)-disp_reference(ip2-1,is))&
                  *(period_target(ip1,is)-period_reference(ip2-1,is))&
                  /(period_reference(ip2,is)-period_reference(ip2-1,is))
                  imark=imark+1
                  exit
               endif
            enddo                                      ! end loop over period of reference dispersion curve
            if(ip2.ne.no_p_reference(is)+1)then
               tt=dist0(is)/disp_target(ip1,is)        ! travel time of the earthquake signal
               t0(is)=t0(is)+tt-dist(is)/vest(ip1,is)
            endif
         enddo                                         ! end loop over period of "earthquake"
         t00=t00+t0(is)
      enddo                                            ! end loop over stations
      t00=t00/imark                                    ! the original time
     ! write(*,*)imark
      ista=0;error_all=0;errorabs=0
      do is=1,nsta                                     ! loop over tele stations
         error(is)=0
         if(no_p_target(is).eq.0.or.no_p_reference(is).eq.0)cycle
         !write(*,*)'hello nsta=',nsta,' t00=',t00
         do ip=1,no_p_target(is)                  ! loop over period of "earthquake"
            if(abs(vest(ip,is)).lt.1e-7)cycle
            tt=dist0(is)/disp_target(ip,is)        ! travel time of the earthquake signal
            error(is)=error(is)+abs(tt-dist(is)/vest(ip,is)-t00)**2  !error1=error1+abs(dist/tt-v)**2
            error_all=error_all+error(is)
            errorabs=errorabs+abs(tt-dist(is)/vest(ip,is)-t00)
         enddo                                     ! end loop over period of earthquake
      enddo                                        ! end loop over stations
      error_all=error_all/imark                    ! travel time error
      write(20,'(5f15.5)')stloo,stlaa,sqrt(error_all),t00,errorabs/imark
   enddo
enddo
close(20)
end program main

! program to locate "earthquake" using grid search
! 2019/07/09 fix t00 and use array
program main
integer,parameter :: maxnp=100,maxnsta=100
real,dimension(maxnp) :: period_reference,disp_reference,period_target,disp_target
real,dimension(maxnsta) :: stla,stlo
real :: stlaa,stloo,disp,stla_reference,stlo_reference
real :: stlo_b,stla_b,stlo_e,stla_e,stplo,stpla,t_step,begin_t,end_t,t,tt1,tt2
real :: vb,ve,stp_v,tp,vv,tt(maxnsta),sto,sta
integer :: iv,nstpv
real :: stlo_g,stla_g,az,t0(maxnsta),t00,error(maxnsta),error_all
real :: vest(maxnsta)
real :: dist(maxnsta),dist0(maxnsta)
real :: errorabs,ddist,serr,serrabs,ft0
character(7)  :: sta_reference,sta_target
character(20) :: filen1,filen2,id,station(maxnsta)
character(80) :: file_reference(maxnsta),file_target(maxnsta)
character(80) :: para,sta_list,list_reference,list_target,output1,output2,ref,target
integer i,j,n,k,np,nsta,nstp_la,nstp_lo,nref,ntarget,nstp_t
integer ilo,ila,ip,if,is,ip1,ip2,it,imark,ista
logical ext1,ext2
if(iargc().ne.1)then
   write(*,*)'Usage: locat locate.par'
   write(*,*)'locate.par:  '
   write(*,*)'station.list '
   write(*,*)'ref_sta'
   write(*,*)'stlo_g,stla_g'
   write(*,*)'output file1/2' 
   write(*,*)'stlo_b stlo_e steplo'
   write(*,*)'stla_b stla_e stepla'
   write(*,*)'begin_v end_v step_v'
   write(*,*)'period'
   call exit(-1)
endif
call getarg(1,para)
open( 11,file=para)
read( 11,*) sta_list
read( 11,*) sta_reference
read( 11,*) stlo_g,stla_g
read( 11,*) output1,output2
read( 11,*) stlo_b,stlo_e,stplo
read( 11,*) stla_b,stla_e,stpla
read( 11,*) vb,ve,stpv
read( 11,*) tp
close(11)
write(*,*)'target period is ', tp
is=1
write(*,'("from",f10.4,"to",f10.4,", from",f10.4,"to",f10.4,",")')stlo_b,stlo_e,stla_b,stla_e
write(*,'("taking station",1x,1a,1x,"as reference.")')sta_reference
! read the station list
! the last one is the one to be found 
! and here it is the reference location assume it is the first guess location
open(7,file=sta_list)
do is=1,maxnsta
   read(7,*,end=20,err=20)station(is), stlo(is),stla(is)
   call cal_dist(stlo_g,stla_g,stlo(is),stla(is),dist0(is))      ! distance between tele station and the grid
enddo
20 close(7)
nsta=is-1                                ! number of station
nstp_lo=int((stlo_e-stlo_b)/stplo)+1     ! number of longitude step
nstp_la=int((stla_e-stla_b)/stpla)+1     ! number of lanitude step
imark=0;vest=0
do is=1,nsta                              ! read in the dispersion curve
   target='dsp_eqe_'//trim(station(is))//'.dat'
   inquire(file=target,exist=ext1)
   if(.not.ext1)cycle
   open(13,file=target) 
   do ip=1,maxnp
      read(13,*,err=14,end=14)period_target(ip),disp_target(ip)
      !write(*,*)period_target(ip),disp_target(ip)
      if(ip.gt.1.and.(period_target(ip)-tp)*(period_target(ip-1)-tp).le.0)then
         vest(is)=disp_target(ip-1)+(disp_target(ip)-disp_target(ip-1))&
         *(tp-period_target(ip-1))/(period_target(ip)-period_target(ip-1))
         tt(is)=dist0(is)/vest(is)    ! travel time of the earthuqake signal
         imark=imark+1
         !write(*,*)tt(is)
         exit
      endif
   enddo
14 continue
   close(13)
enddo
! Begin grid search

nstpv=int((ve-vb)/stpv)+1
open(20,file=output1)
open(21,file=output2)
do iv=1,nstpv
   serr=100000;serrabs=0
   vv=(iv-1)*stpv+vb
   do i=1,nstp_lo                                                       ! loop over lontitude
      stloo=stlo_b+(i-1)*stplo
      do j=1,nstp_la                                                    ! loop over lantitude
         stlaa=stla_b+(j-1)*stpla
         t00=0
         do is=1,nsta                                                   ! loop over tele stations
            dist(is)=0
            call cal_dist(stloo,stlaa,stlo(is),stla(is),dist(is))       ! distance between tele station and the grid
            if(tt(is).ne.0)&
            t00=t00+tt(is)-dist(is)/vv
         enddo                                                          ! end loop over stations
         t00=t00/imark                                                  ! the original time
         error_all=0;errorabs=0
         do is=1,nsta                                                   ! loop over tele stations
            if(tt(is).ne.0)then
               error_all=error_all+(tt(is)-dist(is)/vv-t00)**2/imark     !error1=error1+abs(dist/tt-v)**2
               errorabs =errorabs+ abs(tt(is)-dist(is)/vv-t00)/imark
            endif
         enddo                                                          ! end loop over stations
         !if(vv.eq.3.16)write(20,'(5f15.5)')stloo,stlaa,sqrt(error_all),t00,errorabs
         if(vv.eq.3.06)write(20,'(5f15.5)')stloo,stlaa,sqrt(error_all),t00,errorabs
         if(serr.gt.sqrt(error_all))then
            serr=sqrt(error_all)
            serrabs=errorabs
            ft0=t00 
            sta=stlaa 
            sto=stloo 
            call cal_dist(stloo,stlaa,stlo_g,stla_g,ddist)              ! distance between tele station and the grid
         endif
      enddo                                                             ! end loop over lantitude
   enddo                                                                ! loop over longtitude
   write(21,'(7f15.5)')vv,sto,sta,ddist,serr,serrabs,ft0
enddo
close(20)
close(21)
end program main

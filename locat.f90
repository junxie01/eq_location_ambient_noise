program main
! program to locate "earthquake" using grid search
! virtual earthquake is represented by station
parameter (maxnp=100,maxnsta=100)
real,dimension(maxnp,maxnsta) :: period_ref,disp_ref,period_tgt,disp_tgt
real,dimension(maxnsta) :: stla,stlo
real :: dist,dist1stlaa,stloo,error,disp,stla_ref,stlo_ref
real :: stlo_b,stla_b,stlo_e,stla_e,stplo,stpla,tt,v,t_step,begin_t,end_t,t,tt1,tt2
real :: stlo_g,stla_g,az
character(7)  :: sta_ref,sta_tgt
character(20) :: filen1,filen2,id,station(maxnsta)
character(80) :: file_ref(maxnsta),file_tgt(maxnsta)
character(80) :: para,sta_list,list_ref,list_tgt,output,ref,tgt
integer i,j,n,k,np,nsta,nstp_la,nstp_lo,nref,ntgt,nstp_t
integer ilo,ila,ip,if,is,ip1,ip2,it
integer nump_ref(maxnsta),nump_tgt(maxnsta)
logical ext1,ext2
if(iargc().ne.1)then
        write(*,*)'Usage: locat locate.par'
        write(*,*)'locate.par:  '
        write(*,*)'station.list '
        write(*,*)'reference station'
        write(*,*)'target station'
        write(*,*)'output file' 
        write(*,*)'stlo_b stlo_e steplo'
        write(*,*)'stla_b stla_e stepla'
        write(*,*)'begin_t end_t step_t'
        call exit(-1)
endif
call getarg(1,para)
open(10,file=para)
read(10,*) sta_list
read(10,*) sta_ref
read(10,*) sta_tgt
read(10,*) output
read(10,*) stlo_b,stlo_e,stplo
read(10,*) stla_b,stla_e,stpla
read(10,*) begin_t,end_t,t_step
close(10)
nump_ref=0
nump_tgt=0
is=1
write(*,'("Search station",1x,1a,",")')sta_tgt
write(*,'("from",f10.4,"to",f10.4,", from",f10.4,"to",f10.4,",")')stlo_b,stlo_e,stla_b,stla_e
write(*,'("taking station",1x,1a,1x,"as reference.")')sta_ref
! read the station list
! the last one is the one to be found 
! and here it is the reference location assume it is the first guess location
open(7,file=sta_list)
10 read(7,*,end=20,err=20)station(is), stlo(is),stla(is)
      is=is+1
      goto 10
20   continue
close(7)
nsta=is-1                                ! number of station
nstp_lo=int((stlo_e-stlo_b)/stplo)+1     ! number of longitude step
nstp_la=int((stla_e-stla_b)/stpla)+1     ! number of lanitude step
if(end_t.eq.begin_t)then
       nstp_t=0
else
       nstp_t=int((end_t-begin_t)/t_step)+1     ! number of time step
endif
!write(*,*)nstp_lo,nstp_la,stlo_b,stla_b,stlo_e,stla_e,stplo,stpla
do is=1,nsta
       if(trim(station(is)).eq.trim(sta_ref))then
              stla_ref=stla(is)
              stlo_ref=stlo(is)
       endif
enddo
do is=1,nsta
       if(trim(station(is)).eq.trim(sta_tgt))then
              stla_g=stla(is)
              stlo_g=stlo(is)
       endif
enddo
do is=1,nsta                             ! read in the dispersion curve
!       ref='COR_'//trim(sta_ref)//'_'//trim(station(is))//'.dsp'
!       tgt='COR_'//trim(sta_tgt)//'_'//trim(station(is))//'.dsp'
       ref=trim(sta_ref)//'_'//trim(station(is))//'.dsp'
       tgt=trim(sta_tgt)//'_'//trim(station(is))//'.dsp'
       inquire(file=ref,exist=ext1)
       inquire(file=tgt,exist=ext2)
       if(ext1.and.ext2)then
              !write(*,'("Read in the reference dispersion curve",1x,1a)')trim(ref)
                   ! read in reference group velocity disperison curve 
              ip=1
              open(10,file=ref) 
         12   read(10,*,err=11,end=11)period_ref(ip,is),disp_ref(ip,is)
              ip=ip+1
              goto 12
         11   continue
              nump_ref(is)=ip-1
              call sort(period_ref(:,is),disp_ref(:,is),nump_ref(is),1)
                   ! read in earthquake group velocity disperison curve 
              ip=1
              open(13,file=tgt) 
              !write(*,'("Read in dispersion curve",1x,1a)')trim(tgt)
         15   read(13,*,err=14,end=14)period_tgt(ip,is),disp_tgt(ip,is)
              ip=ip+1
              goto 15
         14   continue
              nump_tgt(is)=ip-1
              call sort(period_tgt(:,is),disp_tgt(:,is),nump_tgt(is),1)
       endif
enddo
! Begin grid search
open(20,file=output)
if (nstp_t.eq.0)then                                                                    ! here we suppose the origional time is known for sure
       do i=1,nstp_lo                                                                   ! loop over lontitude
               stloo=stlo_b+(i-1)*stplo
               do j=1,nstp_la                                                           ! loop over lantitude
                      stlaa=stla_b+(j-1)*stpla
                      error1=0
                      error2=0 
                      do is=1,nsta                                                      ! loop over tele stations
                              if(nump_tgt(is).ne.0.and.nump_ref(is).ne.0)then
                                      write(id,'("File id: ",i5)')is
                                      call cal_dist(stloo,stlaa,stlo(is),stla(is),dist)         ! distance between tele station and the grid
                                      call cal_dist(stlo(is),stla(is),stlo_g,stla_g,dist0)      ! distance between tele station and the first guess location
                                      call cal_dist(stlo_ref,stla_ref,stla(is),stlo(is),dist1)  ! distance between tele station and the reference station
                              !write(*,*)dist,dist0,dist1
                                      do ip1=1,nump_tgt(is)                             ! loop over period of "earthquake"
                                               do ip2=2,nump_ref(is)                    ! loop over period of reference dispersion curve
                                                        if(period_ref(ip2-1,is)<=period_tgt(ip1,is)&
                                                                     .and.period_ref(ip2,is)>period_tgt(ip1,is))then
                                                                      v=disp_ref(ip2-1,is)+(disp_ref(ip2,is)-disp_ref(ip2-1,is))&
                                                                      *(period_tgt(ip1,is)-period_ref(ip2-1,is))&
                                                                      /(period_ref(ip2,is)-period_ref(ip2-1,is))
                                                                      !write(*,*)v
                                                                      goto 16
                                                        endif
                                               enddo
                                          16   continue
                                               if(ip2.ne.nump_ref(is)+1)then
                                                        tt=dist0/disp_tgt(ip1,is)       ! travel time of the earthquake signal
                                                        !error1=error1+abs(dist/tt-v)**2
                                                        error1=error1+abs(dist/v-tt)**2
                                               endif
                                      enddo
                             endif
                     enddo
                     write(20,'(3f15.5)')stloo,stlaa,error1
              enddo
       enddo
else
! here we suppose the origional time is unknown
       !output=trim(output)//'_1'
       ! open(20,file=output)
       do it=1, nstp_t                                                                       ! loop over time
              t=begin_t+(it-1)*t_step
              do i=1,nstp_lo                                                                 ! loop over lontitude
                     stloo=stlo_b+(i-1)*stplo
                     do j=1,nstp_la                                                          ! loop over lantitude
                            stlaa=stla_b+(j-1)*stpla
                            error1=0
                            do is=1,nsta                                                     ! loop over tele stations
                                   write(id,'("File id: ",i5)')is
                                   call cal_dist(stloo,stlaa,stlo(is),stla(is),dist)         ! distance between tele station and the grid
                                   call cal_dist(stlo(is),stla(is),stlo_g,stla_g,dist0)      ! distance between tele station and the first guess location
                                   call cal_dist(stlo_ref,stla_ref,stla(is),stlo(is),dist1)  ! distance between tele station and the reference station
                               !write(*,*)dist,dist0,dist1
                                   ref='COR_'//trim(sta_ref)//'_'//trim(station(is))//'.dsp'
                                   tgt='COR_'//trim(sta_tgt)//'_'//trim(station(is))//'.dsp'
                                   !write(ref,'(1i0,".dsp")')is
                                   !write(tgt,'(1i0,".dsp2")')is
                                   inquire(file=ref,exist=ext1)
                                   inquire(file=tgt,exist=ext2)
                                   if(ext1.and.ext2)then
                                         do ip1=1,nump_tgt(is)           ! loop over period of earthquakes      
                                                  do ip2=2,nump_ref(is)  ! loop over period of reference
                                                            if(period_ref(ip2-1,is)<=period_tgt(ip1,is)&
                                                               .and.period_ref(ip2,is)>period_tgt(ip1,is))then
                                                                      v=disp_ref(ip2-1,is)+(disp_ref(ip2,is)-disp_ref(ip2-1,is))&
                                                                      *(period_tgt(ip1,is)-period_ref(ip2-1,is))&
                                                                      /(period_ref(ip2,is)-period_ref(ip2-1,is))
                                                                      !write(*,*)v
                                                                      goto 17
                                                            endif
                                                  enddo
                                       17         continue
                                                  if(ip2.ne.nump_ref(is)+1)then
                                                          tt1=dist0/disp_tgt(ip1,is)+t ! real group travel time of the earthquake
                                                          !tt2=dist1/v                  ! real travel time of ccf the refernece station
                                                          !error1=error1+abs(dist/tt1-dist1/tt2)**2
                                                          error1=error1+abs(dist/v-tt1)
                                                          !error1=error1+abs(dist/v-tt1-t)
                                                  endif
                                         enddo
                                  endif
                            enddo                                       ! loop over remote earthquakes
                            write(20,'(4f16.5)')t,stloo,stlaa,error1
                    enddo                                               ! loop over lantitude
            enddo                                                       ! loop over lontitude
       enddo                                                            ! loop over time
endif
close(20)
end program main

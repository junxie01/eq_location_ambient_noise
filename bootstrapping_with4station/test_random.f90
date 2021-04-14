program main
integer i,m
real :: h
i=1;m=1
call init_random_seed(i)
open(11,file='num.dat')
do i=1,1000
call random_number(h)
write(11,*)int(h*m)+1
enddo
close(11)
end program

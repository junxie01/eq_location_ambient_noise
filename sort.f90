!subroutine sort(array1,array2,n,order)
subroutine sort(array1,array2,n,order)
real :: array1(n),array2(n),array(n),tmp1(n),temp
integer :: n,order,i,j,index(n)
tmp1=0
do i=1,n-1
         do j=i+1,n
                 if(order>0)then
                         if(array1(i)>array1(j))then
                                 temp=array1(i)
                                 array1(i)=array1(j)
                                 array1(j)=temp
                                 temp=array2(i)
                                 array2(i)=array2(j)
                                 array2(j)=temp
                         endif
                 else
                         if(array1(i)<array1(j))then
                                 temp=array1(i)
                                 array1(i)=array1(j)
                                 array1(j)=temp
                                 temp=array2(i)
                                 array2(i)=array2(j)
                                 array2(j)=temp
                         endif
                 endif
         enddo
!         write(*,*)temp
enddo
return
end subroutine

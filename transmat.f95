program transmat
implicit none
integer::n,i,j
real*8,dimension(:,:),allocatable::a,vx
print*,"enter dimension"
read*,n
allocate (a(n,n))
allocate (vx(n,n))
print*,"enter matrix"
do i=1,n
  read*,(a(i,j),j=1,n)
enddo
vx=transpose(a)
print*,"original matrix:"
do i=1,n
  print*,(a(i,j),j=1,n)
enddo
print*,"transposed matrix:"
do i=1,n
  print*,(vx(i,j),j=1,n)
enddo
end program transmat

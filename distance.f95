program distance
implicit none
integer::i,j,k,l,m,natom
real,dimension(7,3)::A
real,dimension(7,7):: R
integer,dimension(7):: atomic_no
real,dimension(7,7,7)::phi
real::xcord,ycord,zcord,angle,angleop,angle_oop,eklm_x,eklm_y,eklm_z,exx,eyy,ezz
!reading data
open (unit=1,file='acetaldehyde.txt',status="old",action='read')
read(1,*) natom
do i=1,7
  read(1,*) atomic_no(i),(A(i,j),j=1,3)
end do
close(1)
!calculating the bond lengths
do i=1,7
  do k=1,i-1
    R(i,k)=sqrt((A(i,1)-A(k,1))**2 + (A(i,2)-A(k,2))**2+(A(i,3)-A(k,3))**2)
    R(k,i)=R(i,k)
    Write(*,*)R(i,k),i,k
  enddo
enddo

!calculating bond angles
 do i=1,7
   do k=1,i-1
     do l=1,k-1
       if(i/=k .and. k/=l)then
        xcord=((A(k,1)-A(i,1))*(A(k,1)-A(l,1)))
        ycord=((A(k,2)-A(i,2))*(A(k,2)-A(l,2)))
        zcord=((A(k,3)-A(i,3))*(A(k,3)-A(l,3)))

        phi(i,k,l)=(1.0/R(i,k))*(1.0/R(k,l))*(xcord+ycord+zcord)

        angle=acos(phi(i,k,l))*(180.0/acos(-1.0))
        write(*,*)angle,i,k,l
       endif
      end do
   end do
 end do

!calculating out of plane angle
do i=1,4
 do k=2,5
  do l=3,6
   do m=4,7
     eklm_x=(((A(l,2)-A(k,2))*(A(l,3)-A(m,3)))-((A(l,3)-A(k,3))*(A(l,2)-A(m,2))))*(1/R(l,k))*(1/R(l,m))
     !write(*,*)eklm_x
      eklm_y=(((A(l,1)-A(k,1))*(A(l,3)-A(m,3)))-((A(l,1)-A(m,1))*(A(l,3)-A(m,3))))*(1/R(l,k))*(1/R(l,m))
     eklm_z=(((A(l,2)-A(k,2))*(A(l,3)-A(m,3)))-((A(l,2)-A(m,2))*(A(l,3)-A(k,3))))*(1/R(l,k))*(1/R(l,m))

     exx=eklm_x*(A(l,1)-A(i,1))/R(l,i)

     eyy=eklm_y*(A(l,2)-A(i,2))/R(l,i)
     ezz=eklm_y*(A(l,3)-A(i,3))/R(l,i)

     angleop=(exx-eyy+ezz)/sin(phi(k,l,m))*(180.0/3.14)
     angle_oop=asin(angleop)
    end do
   end do
 end do
end do
end program distance
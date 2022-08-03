program HF
implicit none
integer::i,j,ierr,k,l,ij,kl,ijkl,indx,ik,jl,ikjl,c
double precision,dimension(7,7)::ovlp,ke,enu,hcore,ls,lamda,dummy,lamda_inv,lamda_halfinv,s_halfinv,f_prime,c0,c0_prime,epsilon,d0
double precision,dimension(7,7)::f,d1,fdcomm
double precision,dimension(665)::tei
double precision:: nuc,val,e0,e1,deltae,trace
double precision, parameter:: abserr=1.0e-12
e0=0.
tei=0.
d0=0.
c=1
trace=0.
open (unit=13, file="nuc.txt", status="old", action="read")
read(13,*)nuc
close(13)
open (unit=14, file="overlap.txt", status="old", action="read")
do i=1,7
  do j=1,i
    read(14,*)ovlp(i,j)
    ovlp(j,i)=ovlp(i,j)
  enddo
enddo
close(14)
open (unit=15, file="ke.txt", status="old", action="read")
do i=1,7
  do j=1,i
    read(15,*)ke(i,j)
    ke(j,i)=ke(i,j)
  enddo
enddo
close(15)
open (unit=16, file="enu.txt", status="old", action="read")
do i=1,7
  do j=1,i
    read(16,*)enu(i,j)
    enu(j,i)=enu(i,j)
  enddo
enddo
close(16)
hcore=ke+enu
print*,"Core Hamiltonian:"
do i=1,7
  print"(7(f12.8,3x))",(hcore(i,j),j=1,7)
enddo
open (unit=17, file="eri.txt", status="old", action="read", iostat=ierr)
if ( ierr /= 0 ) then
print*, "Failed to open!"
stop
end if
do
  read(17,*,iostat=ierr)i,j,k,l,val
  if(ierr==0)then
    ij=indx(i,j)
    kl=indx(k,l)
    ijkl=indx(ij,kl)
    tei(ijkl)=val
  else
    exit
  endif
enddo
close(16)
lamda(:,:)=ovlp(:,:)
call Jacobi(lamda,ls,abserr,7)
do i=1,7
  do j=1,7
    if(i/=j)then
      lamda(i,j)=0.
    endif
  enddo
enddo
call eigsort(lamda,ls,7)
dummy(:,:)=lamda(:,:)
call inverse(dummy,lamda_inv,7)
lamda_halfinv=sqrt(lamda_inv)
s_halfinv=matmul(matmul(ls,lamda_halfinv),transpose(ls))
print*,"S-1/2 matrix:"
do i=1,7
  print"(7(f12.8,3x))",(s_halfinv(i,j),j=1,7)
enddo
f_prime=matmul(matmul(transpose(s_halfinv),hcore),s_halfinv)
print*,"F0prime matrix:"
do i=1,7
  print"(7(f12.8,3x))",(f_prime(i,j),j=1,7)
enddo
epsilon(:,:)=f_prime(:,:)
call Jacobi(epsilon,c0_prime,abserr,7)
call eigsort(epsilon,c0_prime,7)
c0=matmul(s_halfinv,c0_prime)
print*,"Initial mo coeff:"
do i=1,7
  print"(7(f12.8,3x))",(c0(i,j),j=1,7)
enddo
do i=1,7
  do j=1,7
    do k=1,5
      d0(i,j)= d0(i,j)+c0(i,k)*c0(j,k)
    enddo
  enddo
enddo
print*,"Initial density:"
do i=1,7
  print"(7(f12.8,3x))",(d0(i,j),j=1,7)
enddo
do i=1,7
  do j=1,7
    e0=e0+d0(i,j)*(hcore(i,j)+hcore(i,j))
  enddo
enddo
e0=e0+nuc
print*,"initial electronic energy:",e0
do
  do i=1,7
    do j=1,7
      f(i,j) = hcore(i,j)
      do k=1,7
        do l=1,7
          ij = INDX(i,j)
          kl = INDX(k,l)
          ijkl = INDX(ij,kl)
          ik = INDX(i,k)
          jl = INDX(j,l)
          ikjl = INDX(ik,jl)
          f(i,j) = f(i,j) + d0(k,l) * (2.0 * TEI(ijkl) - TEI(ikjl))
        enddo
      enddo
    enddo
  enddo
  f_prime=matmul(matmul(transpose(s_halfinv),f),s_halfinv)
  epsilon(:,:)=f_prime(:,:)
  call Jacobi(epsilon,c0_prime,abserr,7)
  call eigsort(epsilon,c0_prime,7)
  c0=matmul(s_halfinv,c0_prime)
  d1=0.
  do i=1,7
    do j=1,7
      do k=1,5
        d1(i,j)= d1(i,j)+c0(i,k)*c0(j,k)
      enddo
    enddo
  enddo
  e1=0.
  do i=1,7
    do j=1,7
      e1=e1+d1(i,j)*(hcore(i,j)+f(i,j))
    enddo
  enddo
  e1=e1+nuc
  deltae=e1-e0
  print*,"iteration no.,e0,e1,deltaE:",c,e0,e1,deltae
  if(abs(deltae)<abserr)then
    print*,"Converged!"
    exit
  else
    e0=e1
    d0=d1
  endif
  c=c+1
enddo
print*,"HF energy of water:",e1
end program HF

function indx(i,j)
integer::i,j,indx
if(i>=j)then
  indx=((i*(i+1))/2)+j
else
  indx=((j*(j+1))/2)+i
endif
end function indx

subroutine eigsort(a,x,n)
implicit none
double precision,dimension(n,n)::a,x
integer::i,j,k,n
double precision::temp
do i=1,n-1
  do j=i+1,n
    if (a(i,i)>a(j,j))then
      temp=a(j,j)
      a(j,j)=a(i,i)
      a(i,i)=temp
      do k=1,n
        temp=x(k,j)
        x(k,j)=x(k,i)
        x(k,i)=temp
      enddo
    endif
  enddo
enddo
end subroutine eigsort

subroutine Jacobi(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer i, j, k, n
double precision a(n,n),x(n,n)
double precision abserr, b2, bar
double precision beta, coeff, c, s, cs, sc

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do
return
end

subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse
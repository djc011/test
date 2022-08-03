program diagmat
implicit none
integer::n,i,j,m,it
real*8,dimension(:,:),allocatable::a,vx,r
!real*8,dimension(:),allocatable::r
real*8 eps,d1,d2
double precision, parameter:: abserr=1.0e-15
eps=1.d-10
d1=1.d-8
d2=1.d-8
m=600
print*,"enter dimension"
read*,n
allocate (a(n,n))
allocate (vx(n,n))
allocate (r(n,n))
print*,"enter matrix"
do i=1,n
  read*,(a(i,j),j=1,n)
enddo
r(:,:)=a(:,:)
call Jacobi(r,vx,abserr,n)
!call EPMRI(eps,d1,d2,m,n,a,it,R,VX)
!if (it==-1) then
!      print *,'  No convergence!'
!	  print *,' '
!else
 print*,"eigen values"
  do i=1,n
    print"(5(f12.8,3x))",(r(i,j),j=1,n)
  enddo
  print*,"eigen vectors"
  do i=1,n
    print"(5(f12.8,3x))",(VX(i,j),j=1,n)
  enddo
!endif
end program diagmat

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

 !***************************************************************
  !* Subroutine DECCRM determines the lower triangular matrix and*
  !* the upper triangukar matrix of Crout's decomposition of a   *
  !* given square real matrix, A.                                *
  !* ----------------------------------------------------------- *
  !* INPUTS:                                                     *
  !*         eps: required precision (double)                    *
  !*          n : size of matrix A (integer)                     *
  !*          A : input matrix (n x n)                           *
  !* OUTPUTS:                                                    *
  !*          it: flag, =0 if method does not apply              *
  !*                    =1 if method is ok.                      *
  !*           U: lower triangular matrix.                       *
  !*           V: upper triangular matrix.                       *
  !***************************************************************
  Subroutine DECCRM(eps, n, A, it, U, V)
    real*8 A(n,n), U(n,n), V(n,n)
    real*8 eps, s
    integer n,it

    if (dabs(A(1,1)) < eps)  then
	  it=0
    else
      do i=1, n
	    U(i,1)=A(i,1)
      end do 
      V(1,1)=1.d0
      do j=2, n
	    V(1,j)=A(1,j)/U(1,1)
      end do
      it=1; k=2
      do while (it.ne.0.and.k<=n)
	    do i=1, n
	      if (i < k) then  
		    U(i,k)=0.d0
		  else
	        s=0.d0
	        do j=1, k-1
			  s = s + U(i,j)*V(j,k)
            end do
	        U(i,k) = A(i,k) - s
		  end if
        end do
	    if (dabs(U(k,k)) < eps) then 
          it=0
		else
	      do j=1, n
	        if (j < k) then 
			  V(k,j)=0.d0
	        else if (j==k) then
			  V(k,j)=1.d0
            else
	          s=0.d0
	          do i=1, k-1
			    s = s + U(k,i)*V(i,j)
              end do
	          V(k,j) = A(k,j)/U(k,k);
			end if
          end do
	      k = k + 1
		end if
      end do !while
    end if
	return
  End

!  void MatPrint(char *title, MAT A, int n) {
!  int i,j;
!   printf("\n %s\n", title);
!for (i=1; i<=n; i++) {
!     for (j=1; j<=n; j++)
!       printf(" %f", A[i][j]);
!	  printf("\n");
!   }
! }

! void VecPrint(char *title, VEC X, int n) {
!   int i;
!   printf("\n %s\n", title);
!   for (i=1; i<=n; i++) printf(" %f", X[i]);
!	printf("\n");
! }

  !*********************************************************
  !* Calculate the eigenvalues of a real square matrix by  *
  !* Rutishauser's Method.                                 *
  !* ----------------------------------------------------- *
  !* INPUTS:                                               *
  !*        eps: absolute precision (double)               *
  !*        dta: relative precision (double)               *
  !*         m : maximum number of iterations (integer)    *
  !*         n : size of matrix A (integer)                *
  !*         A : input real square matrix (n x n)          *
  !* OUTPUTS:                                              *
  !*         it: flag, =-1 if convergence is not obtained  *
  !*                   =1 if convergence is ok.            *
  !*         R : contains in output the n eigenvalues of A *
  !*                                                       *         
  !*********************************************************
  Subroutine VAMR(eps, dta, m, n, A, it, R)
    real*8 A(n,n), R(n)
    real*8 eps, dta
    real*8 phi,s,t0
    real*8 U(n,n), V(n,n)
    integer m,n,it

    t0=0.d0
    l=1
    it=0
    do while (l<=m.and.it.ne.1)
      do i=1, n
	    R(i)=A(i,i)
      end do

      call DECCRM(eps, n, A, it, U, V)

      if (it==0) then
	    do i=1, n
		  A(i,i)=A(i,i) + 1.d0
        end do
	    t0 = t0 + 1.d0
      else
	    do i=1, n
		  do j=1, n
	        s=0.d0
	        do k=1, n
			  s = s + V(i,k) * U(k,j)
            end do
	        A(i,j) = s
		  end do
        end do
	    phi=0.d0
        do i=1, n
	      s= dabs(A(i,i)-R(i))
	      if (s > phi)  phi=s
		end do
	    if (phi < dta) then
	      do i=1, n
		    R(i) = A(i,i) - t0
          end do
		else
	      l = l + 1
	      it=-1
		end if
      end if
    end do !while
	return
  End

  !************************************************************
  !* Procedure IIM calculates a real eigenvalue and the asso- *
  !* ciated eigenvector of a real square matrix the inverse   *
  !* iteration method.                                        *
  !* -------------------------------------------------------- *
  !* INPUTS:                                                  *
  !*         eps : absolute precision (double)                *
  !*         dta : relative precision (double)                *
  !*          m  : maximum number of iterations (integer)     *
  !*          n  : size of matrix A                           *
  !*          A  : input real square matrix (n x n)           *
  !* OUTPUTS:                                                 *
  !*          it : flag, =-1 if convergence is not obtained   *
  !*                     =1 if convergence is ok.             *
  !*        Gamma: starting value for the eigenvalue as input *
  !*               approximation of the eigenvalue with preci-*
  !*               sion dta in output.                        *
  !*          X1 : contains in output the associated eigen-   *
  !*               vector.                                    *
  !*                                                          *
  !************************************************************
  Subroutine IIM(eps, dta, m, n, A, it, gamma, X1)
    real*8 A(n,n), X1(n)
    real*8 eps, dta, gamma
    real*8  p0,phi,s,t0
    real*8 W(n), X0(n)
    integer LP(n)

    do i=1, n
	  A(i,i) = A(i,i) - gamma
    end do
    do k=1, n-1
      p0=A(k,k); l0=k
      do i=k+1, n
		if (dabs(A(i,k)) > dabs(p0)) then
	      p0=A(i,k); l0=i
		end if
      end do 
      LP(k)=l0
      if (dabs(p0) < eps) then
	    p0=eps; A(l0,k)=eps
      end if
      if (l0.ne.k) then
	    do j=k, n
	      t0=A(k,j); A(k,j)=A(l0,j); A(l0,j)=t0
		end do
      end if
	  do i=k+1, n
	    A(i,k)=A(i,k)/p0
	    do j=k+1, n
	      A(i,j)=A(i,j)-A(i,k)*A(k,j)
        end do
      end do
	end do !k loop

    if (dabs(A(n,n)) < eps)  A(n,n)=eps
    do i=1, n
	  X0(i)=1.d0/sqrt(1.d0*i)
    end do 
    it=-1; l=1
    do while (it==-1.and.l<=m)
      do i=1, n
	    W(i)=X0(i)
      end do
      do k=1, n-1
        l0=LP(k)
        if (l0.ne.k) then
          t0=W(k); W(k)=W(l0); W(l0)=t0
        end if
        do i=k+1, n
		  W(i)=W(i)-A(i,k)*W(k)
        end do
      end do
      X1(n)=W(n)/A(n,n)
      do i=n-1, 1, -1
        s=0.d0
        do j=i+1, n
		  s = s + A(i,j)*X1(j)
        end do
	    X1(i)=(W(i)-s)/A(i,i)
      end do
      p0=0.d0
      do i=1, n
        if (dabs(X1(i)) > dabs(p0))  p0=X1(i)
      end do
      do i=1, n
        X1(i)=X1(i)/p0
      end do
      phi=0.d0
      do i=1, n
        s=dabs(X1(i)-X0(i))
        if (s > phi)  phi=s
      end do

      if (phi < dta) then
        gamma = gamma + 1.d0/p0
        it=1
      else
        do i=1, n
		  X0(i)=X1(i)
        end do 
        l = l + 1
      end if
	end do !while
	return
  End  !IIM  	  						    


  !******************************************************
  !* INPUTS:                                            *
  !* EPS : precision (Double)                           *
  !* D1  : precision d1 (Double)                        *
  !* D2  : precision d2 (Double)                        *
  !* M   : maximum number of iterations (integer)       *
  !* N   : order of matrix A (integer)                  *
  !* A   : input matrix to study (of MAT type)          *
  !* -------------------------------------------------- *
  !* OUTPUTS:                                           *
  !* IT  : -1 if no convergence is obtained (integer)   *
  !* R   : table of eigenvalues (of VEC type)           *
  !* VX  : table of eigenvectors (of MAT type)          *
  !******************************************************
  Subroutine EPMRI(eps, d1, d2, m, n, A, it, R, VX)
    real*8 A(n,n), R(n), VX(n,n)
    real*8 X(n), A1(n,n)
    real*8 eps,d1,d2
    integer m,n,it

    do i=1, n
      do j=1, n
        A1(i,j) = A(i,j)
      end do
    end do

    call VAMR(eps,d2,m,n,A1,it,R)

! restore A1 after VAMR
    do i=1, n
      do j=1, n
        A1(i,j) = A(i,j)
      end do
	end do   

    j=1
    do while (it==1.and.j<=n)

    call IIM(eps,d1,m,n,A1,it,R(j),X)
  
! restore A1 after IIM
      do i=1, n
        do k=1, n
          A1(i,k) = A(i,k)
        end do
        VX(i,j) = X(i)
      end do
      j = j + 1     
    end do !while
	
	return
  End
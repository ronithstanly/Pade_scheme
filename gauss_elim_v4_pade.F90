MODULE matrices
    INTEGER, PARAMETER :: DP = KIND(1.0d0) 
  REAL, DIMENSION(:,:), ALLOCATABLE :: a
  REAL, DIMENSION(:), ALLOCATABLE :: b,sol
  REAL(DP), DIMENSION(:), ALLOCATABLE :: f,x,exact_first,exact_second
  REAL(DP) :: h
  INTEGER :: n, n_grid !Size of matrix n x n
    
    CONTAINS

    SUBROUTINE A_matrix
        IMPLICIT NONE
        INTEGER :: i,j

        ALLOCATE(a(1:n,1:n))

        a(1:n,1:n) = 0.0d0  

        a(1,1)=1.0d0; a(1,2)=2.0d0;

        j=1;
        DO i=2,n-1
                a(i,j)=1.0d0; a(i,j+1)=4.0d0; a(i,j+2)=1.0d0; 
                j = j+1
        END DO

        a(n,n-1) = 2.0d0; a(n,n) = 1.0d0;

    OPEN(20,file='test_matrix_a.dat',status='replace')
            DO i=1, n
                WRITE(20,'(50f5.1)') a(i,:) !Printing out array assuming max of 50 elements in each row
            END DO
    CLOSE(20)
    
    END SUBROUTINE A_Matrix


    SUBROUTINE B_Matrix
        IMPLICIT NONE
        INTEGER :: i,j

        ALLOCATE(b(1:n))
        ALLOCATE(x(1:n_grid),f(1:n_grid),exact_first(1:n_grid),exact_second(1:n_grid))

        b(1:n) = 0.0d0
    
        x(1:n_grid) = 0.0d0
        f(1:n_grid) = 0.0d0
        exact_first(1:n_grid)   = 0.0d0
        exact_second(1:n_grid)  = 0.0d0
       
        DO i=1, n_grid
            x(i)            = 2.0d0 + FLOAT(i)*h    ! Grid
            f(i)            = x(i)**3  ! Function
            exact_first(i)  = 3.0d0*x(i)**2 ! Exact first derivative
            exact_second(i) = 6.0d0*x(i) 
        END DO

        !B Matrix
        b(1) = (-5.0d0/2.0d0)*f(1) + 2.0d0*f(2) + 0.5d0*f(3)

        DO i=2, n-1
                b(i) = 3.0d0*(f(i+1)-f(i-1))
        END DO

        b(n)    = (5.0d0/2.0d0)*f(n) - 2.0d0*f(n-1) - 0.5d0*f(n-2)   
        
        b(:)    = b(:)/h
    OPEN(24,file='test_matrix_b.dat',status='replace')
            DO i=1, n
                WRITE(24,'(50f5.1)') b(i) !Printing out array assuming max of 50 elements in each row
            END DO
    CLOSE(24)
              
    OPEN(21,file='test_out.dat',status='replace')
            WRITE(21,*) 'x, f, exact_first, exact_second'
            DO i=1, n_grid
                WRITE(21,'(100g15.5)') x(i),f(i),exact_first(i),exact_second(i)
            END DO
    CLOSE(21)

    END SUBROUTINE B_Matrix


    SUBROUTINE invert
        IMPLICIT NONE
        INTEGER :: i,j

    OPEN(22,file='test_sol_first_der.dat',status='replace')
            DO i=1, n
                IF (MOD(i,2).NE.0) THEN ! Odd
                    WRITE(22,*) sol(i)
                END IF 
            END DO
    CLOSE(22)

    OPEN(23,file='test_sol_sec_der.dat',status='replace')
            DO i=1, n
                IF (MOD(i,2)==0) THEN ! Even
                    WRITE(23,*) sol(i) ! change this
                END IF 
            END DO
    CLOSE(23)

    END SUBROUTINE invert

END MODULE matrices

!###########################
PROGRAM gauss_elim
  USE matrices

    n = 10 ! Double of grid points (Matrices are n/2 X n/2)
    n_grid = n ! Grid points
    h = 1/(n/2.0d0)!4.0d0/FLOAT(n_grid)!0.1d0 ! Grid spacing

    CALL A_Matrix
    CALL B_Matrix

        ALLOCATE(sol(1:n))
        sol(1:n) = 0.0d0

        sol = solve_wbs(ge_wpp(a,b))
    
        PRINT'(f15.7)',sol

    CALL invert
 
	CONTAINS
 
	FUNCTION solve_wbs(u) result(x) ! solve with backward substitution
      REAL                 :: u(:,:)
      INTEGER              :: i,n
      REAL   , ALLOCATABLE :: x(:)
      n = SIZE(u,1)
      ALLOCATE(x(n))
        FORALL (i=n:1:-1) x(i) = ( u(i,n+1) - SUM(u(i,i+1:n)*x(i+1:n)) ) / u(i,i)
    END FUNCTION
 
    FUNCTION  ge_wpp(a,b) result(u) ! gaussian eliminate with partial pivoting
      real                 :: a(:,:),b(:),upi
      integer              :: i,j,n,p
      real   , allocatable :: u(:,:)
      n = SIZE(a,1)
      u = reshape( [a,b], [n,n+1] )
      DO j=1,n
        p = MAXLOC(ABS(u(j:n,j)),1) + j-1 ! maxloc returns indices between (1,n-j+1)
        IF (p /= j) u([p,j],j) = u([j,p],j)
        u(j+1:,j) = u(j+1:,j)/u(j,j)
        DO i=j+1,n+1
          upi = u(p,i)
          IF (p /= j) u([p,j],i) = u([j,p],i)
          u(j+1:n,i) = u(j+1:n,i) - upi*u(j+1:n,j)
        END DO
      END DO
    END FUNCTION
 
END PROGRAM gauss_elim
 

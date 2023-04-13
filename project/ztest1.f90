! Arnoldi iteration to solve for eigenvalue problem
! solve A*x=lamda*M*x
! parallel LU decompostion
! All the matrix info are stored in 0 proc 
! Lu Chen 2023/4/13

Program Ztest1

IMPLICIT NONE

INCLUDE 'mpif.h'
INCLUDE 'zmumps_struc.h'

	
TYPE (ZMUMPS_STRUC) mumps_par   !< MUMPS solver  Matrix A-sigma*M

INTEGER IERR    !< error code of the MUMPS solver
INTEGER,parameter::nsize=50958  !<nsize
INTEGER,parameter::nev=40  !<nev 
INTEGER,parameter::ncv=120  !<ncv

INTEGER::ido,iparam(11),ipntr(14),lworkl,info
COMPLEX*16::resid(nsize),v(nsize,ncv),workd(3*nsize),workl(3*ncv**2 + 5*ncv)
REAL*16:: rwork(ncv)
CHARACTER(2):: which
CHARACTER:: bmat
REAL*16::tol
LOGICAL::select(ncv)
COMPLEX*16::sigma,workev(3*ncv)
INTEGER::ierrA,ierrM
COMPLEX*16:: eigen_vec(nsize,ncv) !<egienvector
COMPLEX*16 :: d(nsize)    !<egienvalue
COMPLEX*16 ::w(nsize)     !<RHS vector
COMPLEX*16 ::u(nsize)     !<solution vector u 
REAL*16 :: time_begin,time_end,total_time

INTEGER:: I,J

INTEGER:: Nm,NNZm
INTEGER,ALLOCATABLE::MIRN(:)
INTEGER,ALLOCATABLE::MJCN(:)
COMPLEX*16,ALLOCATABLE::MA(:)

CALL MPI_INIT(IERR)
time_begin = MPI_Wtime()
! Define a communicator for the package.
mumps_par%COMM = MPI_COMM_WORLD
!  Initialize an instance of the package
!  for L U factorization (sym = 0, with working host)
mumps_par%JOB = -1
mumps_par%SYM = 0
mumps_par%PAR = 1
CALL ZMUMPS(mumps_par)

IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'error1 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
    ERROR STOP 200
ENDIF    


!locd Matrix M and (A-sigma*M)
! M first
IF ( mumps_par%MYID .eq. 0 ) THEN
    OPEN(5,file="input_M")
    READ(5,*) mumps_par%N
    READ(5,*) mumps_par%NNZ
    ALLOCATE( mumps_par%IRN ( mumps_par%NNZ ) )
    ALLOCATE( mumps_par%JCN ( mumps_par%NNZ ) )
    ALLOCATE( mumps_par%A( mumps_par%NNZ ) )
    ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )
    DO I = 1, mumps_par%NNZ
       READ(5,*) mumps_par%IRN(I),mumps_par%JCN(I),mumps_par%A(I)
    END DO
    Nm = mumps_par%N
    NNZm = mumps_par%NNZ
    ALLOCATE(MIRN(NNZm))
    ALLOCATE(MJCN(NNZm))   
    ALLOCATE(MA(NNZm))
    MIRN = mumps_par%IRN
    MJCN = mumps_par%JCN
    MA = mumps_par%A
END IF

IF ( mumps_par%MYID .eq. 0 )THEN
	DEALLOCATE( mumps_par%IRN )
	DEALLOCATE( mumps_par%JCN )
	DEALLOCATE( mumps_par%A )
	DEALLOCATE( mumps_par%RHS )
ENDIF

! (A-sigma*M) second
IF ( mumps_par%MYID .eq. 0 ) THEN
    OPEN(5,file="input_L")
    READ(5,*) mumps_par%N
    READ(5,*) mumps_par%NNZ
    ALLOCATE( mumps_par%IRN ( mumps_par%NNZ ) )
    ALLOCATE( mumps_par%JCN ( mumps_par%NNZ ) )
    ALLOCATE( mumps_par%A( mumps_par%NNZ ) )
    ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )    
    DO I = 1, mumps_par%NNZ
       READ(5,*) mumps_par%IRN(I),mumps_par%JCN(I),mumps_par%A(I)
    END DO
END IF


! Analysis and factorisation 
mumps_par%JOB = 4
CALL ZMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'error2 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
        ERROR STOP 200
ENDIF   



!Using ARPACK to find the needed eigen values and eigen vectors
!-----------------------------------------------------------------------------
ido=0
tol=0.0
iparam(1)=1     ! exact shifts
iparam(3) = 1000 ! max iterations
iparam(7) = 3 
sigma = (0.0,0.74)
lworkl=3*ncv**2 + 5*ncv
info=0

! define the spectral 
! the closet eigenvalues to the sigma
 which = 'LM'
 bmat  = 'G'

!******************************************************************************   
DO 
	! %---------------------------------------------%
	! | Repeatedly call the routine ZNAUPD and take |
	! | actions indicated by parameter IDO until    |
	! | either convergence is indicated or maxitr   |
	! | has been exceeded.                          |
	! %---------------------------------------------%


        IF( mumps_par%MYID .eq. 0 ) THEN	
		CALL znaupd ( ido, bmat, nsize, which, nev, tol, resid, ncv, &
				&     v, nsize, iparam, ipntr, workd, workl, lworkl,& 
				&     rwork, info )
		write(0,*)'ido is',ido
	ENDIF
			
	! important (ido for other mumps_par%MYID .nq. 0)	
	call MPI_BCAST(ido,1,MPI_INT,0,MPI_COMM_WORLD,ierrM)
	
        
        IF (ido .eq. -1)  THEN	

          ! %-------------------------------------------%
          ! | Perform  y <--- OP*x = inv[A-SIGMA*M]*M*x |
          ! | to force starting vector into the range   |
          ! | of OP.   The user should supply his/her   |
          ! | own matrix vector multiplication routine  |
          ! | and a linear system solver.  The matrix   |
          ! | vector multiplication routine should take |
          ! | workd(ipntr(1)) as the input. The final   |
          ! | result should be returned to              |
          ! | workd(ipntr(2)).                          |
          ! %-------------------------------------------%
          
          	! Call matrix vector multiplication
         	
		CALL MV(workd(ipntr(2)),nsize,NNZm,MIRN,MJCN,MA,workd(ipntr(1)))

		
		CALL TRANSFORM(w,nsize,workd(ipntr(2))) 
						
		IF ( mumps_par%MYID .eq. 0 ) THEN
			mumps_par%RHS=(0.0,0.0)
			DO i = 1,nsize
			   mumps_par%RHS(i) = w(i)
			ENDDO
		ENDIF

		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error3 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF    

		
		!  Call package for solution of (A-sigma*M)*x=b
		mumps_par%JOB = 3		
		CALL ZMUMPS(mumps_par)
		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error4 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF    

		!  Solution has been assembled on the host
		IF ( mumps_par%MYID .eq. 0 ) THEN
		    u=mumps_par%RHS
		END IF
				
		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error5 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF   	
			
        	CALL TRANSFORM(workd(ipntr(2)),nsize,u)      
        	
		! %-----------------------------------------%
		! | L O O P   B A C K to call ZNAUPD again. |
		! %-----------------------------------------%
	
        ELSEIF ( ido .eq. 1)  THEN	
        
          ! %-----------------------------------------%
          ! | Perform y <-- OP*x = inv[A-sigma*M]*M*x |
          ! | M*x has been saved in workd(ipntr(3)).  |
          ! | The user only need the linear system    |
          ! | solver here that takes workd(ipntr(3))  |
          ! | as input, and returns the result to     |
          ! | workd(ipntr(2)).                        |
          ! %-----------------------------------------%
                
		!  Call package for solution of (A-sigma*M)*x=b
		CALL TRANSFORM(w,nsize,workd(ipntr(3)))      	

		IF ( mumps_par%MYID .eq. 0 ) THEN
			mumps_par%RHS=(0.0,0.0)
			DO i = 1,nsize
			   mumps_par%RHS(i) = w(i)
			ENDDO
		ENDIF

		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error3 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF    

		
		!  Call package for solution
		mumps_par%JOB = 3		
		CALL ZMUMPS(mumps_par)
		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error4 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF    

		!  Solution has been assembled on the host
		IF ( mumps_par%MYID .eq. 0 ) THEN
		    u=mumps_par%RHS
		END IF
	
		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error5 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF   	
			
        	CALL TRANSFORM(workd(ipntr(2)),nsize,u) 
        	
		! %-----------------------------------------%
		! | L O O P   B A C K to call ZNAUPD again. |
		! %-----------------------------------------%
        
        
        ELSEIF ( ido .eq. 2)  THEN	       	

          ! %---------------------------------------------%
          ! |          Perform  y <--- M*x                |
          ! | Need matrix vector multiplication routine   |
          ! | here that takes workd(ipntr(1)) as input    |
          ! | and returns the result to workd(ipntr(2)).  |
          ! %---------------------------------------------%
		
          	! Call matrix vector multiplication
      	
		CALL MV(u,nsize,NNZm,MIRN,MJCN,MA,workd(ipntr(1)))
		
        	CALL TRANSFORM(workd(ipntr(2)),nsize,u)		

		! %-----------------------------------------%
		! | L O O P   B A C K to call ZNAUPD again. |
		! %-----------------------------------------%	
	
	ELSE
			
		EXIT
	ENDIF 
	
ENDDO
!********************************************************************************
	
IF ( mumps_par%MYID .eq. 0 ) THEN	
CALL zneupd (.TRUE., 'A', select, d, v, nsize, sigma, workev, &
           & bmat, nsize, which, nev, tol, resid, ncv, v, &
           & nsize, iparam, ipntr, workd, workl, lworkl, rwork, ierrA )
ENDIF                        
call MPI_BCAST(ido,1,MPI_INT,0,MPI_COMM_WORLD,ierrA)
   
IF( ierrA < 0 ) THEN
    ERROR STOP 106
ELSE
	eigen_vec(:,1:nev)=v(:,1:nev)
ENDIF


IF ( mumps_par%MYID .eq. 0 ) THEN
	write(0,*) 'eigenvalue',d(1),d(2),d(3),d(4)
ENDIF

time_end = MPI_Wtime()
total_time=time_end - time_begin
IF  ( mumps_par%MYID .eq. 0 ) THEN
	write(0,*) 'total_time',total_time
ENDIF


!  Deallocate user data
IF ( mumps_par%MYID .eq. 0 )THEN
	DEALLOCATE( mumps_par%IRN )
	DEALLOCATE( mumps_par%JCN )
	DEALLOCATE( mumps_par%A )
	DEALLOCATE( mumps_par%RHS )
	DEALLOCATE( MIRN )
	DEALLOCATE( MJCN )
	DEALLOCATE( MA )	
ENDIF

! Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
CALL ZMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'error6 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
    ERROR STOP 200
ENDIF    
CALL MPI_FINALIZE(IERR)


END 
!=========================================================

!*********************************************************
SUBROUTINE TRANSFORM(x,nsize,u)

IMPLICIT NONE

INTEGER::nsize 
COMPLEX*16 ::u(nsize)    !< input
COMPLEX*16 ::x(nsize)    !< solution vector 

INTEGER::i

DO i =1,nsize
   x(i) = u(i)
ENDDO
   

END SUBROUTINE TRANSFORM
!*********************************************************



!*********************************************************
SUBROUTINE MV(x,nsize,NNZm,MIRN,MJCN,MA,u)

IMPLICIT NONE

INTEGER::nsize,NNZm
COMPLEX*16 ::u(nsize)    !< input
COMPLEX*16 ::x(nsize)    !< solution vector 
INTEGER::MIRN(NNZm)
INTEGER::MJCN(NNZm)
COMPLEX*16::MA(NNZm)

INTEGER::i

DO i =1,nsize
   x(i) = (0.0,0.0)
ENDDO

DO i =1, NNZm
   x(MIRN(i)) = x(MIRN(i)) + MA(i) * u(MJCN(i))
ENDDO


END SUBROUTINE MV
!*********************************************************


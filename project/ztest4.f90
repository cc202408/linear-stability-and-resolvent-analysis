! Arnoldi iteration to solve for eigenvalue problem
! solve A*x=lamda*M*x
! parallel LU decompostion
! All the matrix info are stored in 0 proc 
! Lu Chen 2023/4/13

! Load A and M (instead of (A-sigma*M) and M)
! Lu Chen 2023/10/10

! M is not Hermitian Positive Semi-Definite
! standard eigenvalue problem C*x=lamda*x instead of A*x=lamda*M*x
! C=inv(A-sigma*M)*M
! two steps for the v<-C*w
! (1) z<-M*w
! (2) v<-inv(A-sigma*M)*w
! Lu Chen 2023/10/31


Program Ztest4

IMPLICIT NONE

INCLUDE 'mpif.h'
INCLUDE 'zmumps_struc.h'

	
TYPE (ZMUMPS_STRUC) mumps_par   !< MUMPS solver  Matrix A-sigma*M

INTEGER IERR    !< error code of the MUMPS solver
INTEGER::nsize
INTEGER,parameter::nev=5  !<nev 
INTEGER,parameter::ncv=25  !<ncv

INTEGER::ido,iparam(11),ipntr(14),lworkl,info
COMPLEX*16,ALLOCATABLE::resid(:),v(:,:),workd(:)
COMPLEX*16::workl(3*ncv**2 + 5*ncv)
REAL*16:: rwork(ncv)
CHARACTER(2):: which
CHARACTER:: bmat
REAL*16::tol
LOGICAL::select(ncv)
COMPLEX*16::sigma,workev(3*ncv)
INTEGER::ierrA,ierrM
COMPLEX*16,ALLOCATABLE::eigen_vec(:,:),d(:),w(:),u(:)
REAL*16 :: time_begin,time_end,total_time

REAL*16 :: shiftOP

CHARACTER(len=100)::filename

INTEGER:: I,J
INTEGER:: ITE
INTEGER:: indexA,indexM,indexL
INTEGER:: INTtemp1,INTtemp2

INTEGER:: NM,NNZM,NA,NNZA,NNZL
INTEGER:: NMAX
INTEGER,ALLOCATABLE::MIRN(:),AIRN(:),LIRN(:),LIRN_TMP(:)
INTEGER,ALLOCATABLE::MJCN(:),AJCN(:),LJCN(:),LJCN_TMP(:)
COMPLEX*16,ALLOCATABLE::MA(:),AA(:),LA(:),LA_TMP(:)

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
IF  (mumps_par%MYID .eq. 0 ) THEN
! M first
    OPEN(5,file="B.txt")
    READ(5,*)
    READ(5,*)
    READ(5,*)
    READ(5,*) NM,INTtemp1,INTtemp2,NNZM
    ALLOCATE(MIRN(NNZM))
    ALLOCATE(MJCN(NNZM))   
    ALLOCATE(MA(NNZM))
    DO I = 1, NNZM
       READ(5,*) MIRN(I),MJCN(I),MA(I)
    END DO
    
! A second   
    OPEN(5,file="A.txt")
    READ(5,*)
    READ(5,*)
    READ(5,*)
    READ(5,*) NA,INTtemp1,INTtemp2,NNZA
    ALLOCATE(AIRN(NNZA))
    ALLOCATE(AJCN(NNZA))   
    ALLOCATE(AA(NNZA))
    DO I = 1, NNZA
       READ(5,*) AIRN(I),AJCN(I),AA(I)
    END DO   
ENDIF

! calculate (A-sigma*M)
! loal shift value
IF ( mumps_par%MYID .eq. 0 ) THEN
    OPEN(5,file="shiftOP.txt")
    READ(5,*) shiftOP
ENDIF

sigma = cmplx(0.0,shiftOP)

IF ( mumps_par%MYID .eq. 0 )THEN

	!A+(-sigma*M)
	indexA=1
	indexM=1
	indexL=1
	NMAX = NNZM + NNZA
	ALLOCATE(LIRN_TMP(NMAX))
	ALLOCATE(LJCN_TMP(NMAX))   
	ALLOCATE(LA_TMP(NMAX))
       
       
	DO
	     IF (indexA .LE. NNZA .or. indexM .LE. NNZM) THEN
		  
		  IF (AIRN(indexA) .LT. MIRN(indexM)) THEN
		  	
		  	LIRN_TMP(indexL) = AIRN(indexA)
		  	LJCN_TMP(indexL) = AJCN(indexA)
		  	LA_TMP(indexL) = AA(indexA)
		  	indexA = indexA + 1
		  	indexL = indexL + 1
		  
		  ELSEIF (AIRN(indexA) .GT. MIRN(indexM)) THEN
		  
		  	LIRN_TMP(indexL) = MIRN(indexM)
		  	LJCN_TMP(indexL) = MJCN(indexM)
		  	LA_TMP(indexL) = MA(indexM)*(-sigma)
		  	indexM = indexM + 1
		  	indexL = indexL + 1
		  			 
		  ELSE
		  	IF (AJCN(indexA) .LT. MJCN(indexM)) THEN
		  	
			  	LIRN_TMP(indexL) = AIRN(indexA)
			  	LJCN_TMP(indexL) = AJCN(indexA)
			  	LA_TMP(indexL) = AA(indexA)
			  	indexA = indexA + 1
			  	indexL = indexL + 1
		  	
		        ELSEIF (AJCN(indexA) .GT. MJCN(indexM)) THEN	
		        
			  	LIRN_TMP(indexL) = MIRN(indexM)
			  	LJCN_TMP(indexL) = MJCN(indexM)
			  	LA_TMP(indexL) = MA(indexM)*(-sigma)
			  	indexM = indexM + 1
			  	indexL = indexL + 1	
			  		        
		        ELSE
		        
			  	LIRN_TMP(indexL) = AIRN(indexA)
			  	LJCN_TMP(indexL) = AJCN(indexA)
			  	LA_TMP(indexL) = AA(indexA) + MA(indexM)*(-sigma)
			  	indexA = indexA + 1
			  	indexM = indexM + 1
			  	indexL = indexL + 1		        
		        
		        ENDIF	  	
		  	
		  
		  ENDIF
	          
	     
	     ELSE
	     
		  EXIT
		  
	     ENDIF

	ENDDO	
	
	indexL = indexL - 1
	
	ALLOCATE(LIRN(indexL))
	ALLOCATE(LJCN(indexL))   
	ALLOCATE(LA(indexL))
	
	
	DO I = 1, indexL
		LIRN(I) = LIRN_TMP(I)
		LJCN(I) = LJCN_TMP(I)
		LA(I) = LA_TMP(I)	
	ENDDO
	DEALLOCATE(LIRN_TMP)
	DEALLOCATE(LJCN_TMP)   
	DEALLOCATE(LA_TMP)
	
			    	
ENDIF


! load matrix (A-sigma*M)        		    		 
IF ( mumps_par%MYID .eq. 0 )THEN      		
	mumps_par%N = NA
	mumps_par%NNZ = indexL
	ALLOCATE( mumps_par%IRN ( mumps_par%NNZ ) )
	ALLOCATE( mumps_par%JCN ( mumps_par%NNZ ) )
	ALLOCATE( mumps_par%A( mumps_par%NNZ ) )
	ALLOCATE( mumps_par%RHS ( mumps_par%N  ) )    
	DO I = 1, mumps_par%NNZ
		mumps_par%IRN(I) = LIRN(I)
		mumps_par%JCN(I) = LJCN(I)
		mumps_par%A(I) = LA(I)
	ENDDO							       		
ENDIF 

! Analysis and factorisation 
!mumps_par%ICNTL(14) = 30
mumps_par%JOB = 4
CALL ZMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'error2 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
        ERROR STOP 200
ENDIF   

!allocate varabie
nsize=NA

ALLOCATE(eigen_vec(nsize,ncv)) !<egienvector
ALLOCATE(d(nsize))    !<egienvalue
ALLOCATE(w(nsize))     !<RHS vector
ALLOCATE(u(nsize))     !<solution vector u 

ALLOCATE(resid(nsize))
ALLOCATE(v(nsize,ncv))
ALLOCATE(workd(3*nsize))


!Using ARPACK to find the needed eigen values and eigen vectors
!-----------------------------------------------------------------------------
ido=0
tol=0.0
iparam(1)= 1     ! exact shifts
iparam(3) = 100000 ! max iterations
iparam(7) = 1 
lworkl=3*ncv**2 + 5*ncv
info=0

! define the spectral 
! the closet eigenvalues to the sigma
 which = 'LM'
 bmat  = 'I'

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
	
        
        if (ido .eq. -1 .or. ido .eq. 1) then
	

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


! transform eigenvalue: lamda=sigma+1/miu
IF ( mumps_par%MYID .eq. 0 ) THEN
    DO I = 1, nev
       d(I)=sigma+1.0/d(I)
    ENDDO
ENDIF
   
IF( ierrA < 0 ) THEN
    ERROR STOP 106
ELSE
	eigen_vec(:,1:nev)=v(:,1:nev)
ENDIF

! output eigenvector
IF ( mumps_par%MYID .eq. 0 ) THEN
    DO I = 1, nev
        WRITE(filename,'(a,i0,a)')'eigenvector',I,'.txt'
	OPEN(1,file=trim(filename),status='replace')
	DO J = 1, NA
	    write(1,*) eigen_vec(J,I)
	ENDDO
    ENDDO
END IF

IF ( mumps_par%MYID .eq. 0 ) THEN
	write(0,*) 'eigenvalue',d(1),d(2)
ENDIF

IF ( mumps_par%MYID .eq. 0 ) THEN
    OPEN(1,file="result_eigenvalue.txt")
    DO I = 1, nev
    	write(1,*) d(I)
    ENDDO
END IF


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


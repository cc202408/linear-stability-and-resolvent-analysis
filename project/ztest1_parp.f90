! Arnoldi iteration to solve for eigenvalue problem
! solve A*x=lamda*M*x
! parallel LU decompostion with MUMPS
! All the matrix info are stored in 0 proc 
! Lu Chen 2023/4/13

! M is not Hermitian Positive Semi-Definite
! standard eigenvalue problem C*x=miu*x instead of A*x=lamda*M*x
! C=inv(A-sigma*M)*M
! two steps for the v<-C*w
! (1) z<-M*w
! (2) v<-inv(A-sigma*M)*w
! Lu Chen 2023/7/1

! input: matrix A M
! calculate (A-sigma*M) and solve for inv(A-sigma*M)*M*x=miu*x
! Lu Chen 2023/7/3

! resolvent analysis
! input: matrix A M P Mf Mq
! calculate (P'*R'*Mq*R*P)*f = (sigma^2)*Mf*f
! Lu Chen 2023/7/6

! update (A^H)*x=b
! Lu Chen 2023/7/7

! solve for inv(Mf) with SuperLU
! MF/MQ/P -> real*16 instead of complex*16
! Lu Chen 2023/7/9

! update I/O for FreeFem++ matrix in COO directly
! Lu Chen 2023/7/10

! resolvent forcing and responding mode
! Lu Chen 2023/7/16

Program Ztest1

IMPLICIT NONE

INCLUDE 'mpif.h'
INCLUDE 'zmumps_struc.h'

	
TYPE (ZMUMPS_STRUC) mumps_par   !< MUMPS solver  Matrix A-sigma*M

INTEGER IERR    !< error code of the MUMPS solver
INTEGER,parameter::nev=10   !<nev 
INTEGER,parameter::ncv=30   !<ncv

INTEGER::ido,iparam(11),ipntr(14),lworkl,info
COMPLEX*16,ALLOCATABLE::resid(:),v(:,:),workd(:)
COMPLEX*16::workl(3*ncv**2 + 5*ncv)
REAL*16:: rwork(ncv)
CHARACTER(2):: which
CHARACTER:: bmat
REAL*16::tol
LOGICAL::select(ncv)
COMPLEX*16::sigma,workev(3*ncv),omega
INTEGER::ierrA,ierrM
COMPLEX*16,ALLOCATABLE::eigen_vec(:,:),d(:),w1(:),w2(:),u(:)
REAL*16 :: time_begin,time_end,total_time


INTEGER:: N_omega
REAL*16::Omega_MAX,Omega_begin
COMPLEX*16,ALLOCATABLE::singular_value(:,:)
COMPLEX*16,ALLOCATABLE::forcing_mode(:)
COMPLEX*16,ALLOCATABLE::responding_mode(:)

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

!matrix for P,MF,MQ,P'
INTEGER:: NP,NNZP,NMF,NNZMF,NMQ,NNZMQ
INTEGER,ALLOCATABLE::PIRN(:),MFIRN(:),MQIRN(:),PHIRN(:)
INTEGER,ALLOCATABLE::PJCN(:),MFJCN(:),MQJCN(:),PHJCN(:)
REAL*16,ALLOCATABLE::PA(:),MQA(:),MFA(:),PHA(:)

!SuperLU
INTEGER,ALLOCATABLE::rowind(:),colptr(:)
REAL*16,ALLOCATABLE::values0(:)
COMPLEX*16,ALLOCATABLE:: values(:),b(:)
INTEGER nrhs,ldb,infoSuperLu,iopt
INTEGER*8 factors

!transfer st_to_ccs
INTEGER i_max,i_min,i4vec_max,i4vec_min,J_max,j_min,ncc


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


!--------------------------------------------------------------------------------
!load Matrix 
IF ( mumps_par%MYID .eq. 0 ) THEN
! M first
    OPEN(5,file="M.txt")
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
    
!Load MF   
    OPEN(5,file="MF.txt")
    READ(5,*)
    READ(5,*)
    READ(5,*)
    READ(5,*) NMF,INTtemp1,INTtemp2,NNZMF
    ALLOCATE(MFIRN(NNZMF))
    ALLOCATE(MFJCN(NNZMF))   
    ALLOCATE(MFA(NNZMF))
    DO I = 1, NNZMF
       READ(5,*) MFIRN(I),MFJCN(I),MFA(I)
    END DO   

!Load MQ 
    OPEN(5,file="MQ.txt")
    READ(5,*)
    READ(5,*)
    READ(5,*)
    READ(5,*) NMQ,INTtemp1,INTtemp2,NNZMQ
    ALLOCATE(MQIRN(NNZMQ))
    ALLOCATE(MQJCN(NNZMQ))   
    ALLOCATE(MQA(NNZMQ))
    DO I = 1, NNZMQ
       READ(5,*) MQIRN(I),MQJCN(I),MQA(I)
    END DO  

!Load P
    OPEN(5,file="P.txt")
    READ(5,*)
    READ(5,*)
    READ(5,*)
    READ(5,*) NP,INTtemp1,INTtemp2,NNZP
    ALLOCATE(PIRN(NNZP))
    ALLOCATE(PJCN(NNZP))   
    ALLOCATE(PA(NNZP))
    DO I = 1, NNZP
       READ(5,*) PIRN(I),PJCN(I),PA(I)
    END DO  
    ALLOCATE(PHIRN(NNZP))
    ALLOCATE(PHJCN(NNZP))   
    ALLOCATE(PHA(NNZP))    
    DO I = 1, NNZP
       PHIRN(I) = PJCN(I)
       PHJCN(I) = PIRN(I)
       PHA(I) = PA(I)
    ENDDO

END IF

!--------------------------------------------------------------------------------
!factorizae matrix MF with SuperLU
!Step1 Transform Matrix MF in SuperLU format (CCS)

i_min = i4vec_min ( NNZMF, MFIRN )
i_max = i4vec_max ( NNZMF, MFIRN )
j_min = i4vec_min ( NNZMF, MFJCN )
j_max = i4vec_max ( NNZMF, MFJCN )

!Get the CCS size.
CALL st_to_ccs_size ( NNZMF, MFIRN, MFJCN, ncc )


!Create the CCS indices.
ALLOCATE ( rowind(1:ncc) )
ALLOCATE ( colptr(1:NMF+1) )
CALL st_to_ccs_index ( NNZMF, MFIRN, MFJCN, ncc, NMF, rowind, colptr )
  
!Create the CCS values.
ALLOCATE ( values0(1:ncc) )
CALL st_to_ccs_values ( NNZMF, MFIRN, MFJCN, MFA, ncc, NMF, rowind, colptr, values0 )


ALLOCATE ( values(1:NNZMF) )
DO I = 1,NNZMF
    values(I)=cmplx(values0(I),0.0)
ENDDO
      


!Step2 The factors are stored in *factors* handle.

nrhs = 1
ldb  = NMF

ALLOCATE (b(1:NMF))
DO I =1, NMF
   b(I) = cmplx(0.0,0.0)
ENDDO

IF ( mumps_par%MYID .eq. 0 )THEN
	iopt = 1
	CALL c_fortran_zgssv( iopt, NMF, NNZMF, nrhs, values, rowind, colptr, &
                      &    b, ldb, factors, infoSuperLu )
ENDIF                      

!---------------------------------------------------------------------------------
!Setting Omega
!loal input_file
IF ( mumps_par%MYID .eq. 0 ) THEN
    OPEN(5,file="input_omega.txt")
    READ(5,*) Omega_begin,Omega_MAX
ENDIF
N_omega = 20
ALLOCATE(singular_value(N_omega,nev))

DO ITE = 1, N_omega

!calculate A-omega*M
!--------------------------------------------------------------------------------

omega=cmplx(0.0,Omega_begin+ITE*(Omega_MAX-Omega_begin)/N_omega*1.0)

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
		  	LA_TMP(indexL) = MA(indexM)*(-omega)
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
			  	LA_TMP(indexL) = MA(indexM)*(-omega)
			  	indexM = indexM + 1
			  	indexL = indexL + 1	
			  		        
		        ELSE
		        
			  	LIRN_TMP(indexL) = AIRN(indexA)
			  	LJCN_TMP(indexL) = AJCN(indexA)
			  	LA_TMP(indexL) = AA(indexA) + MA(indexM)*(-omega)
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


! load matrix (A-omega*M)        		 
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
mumps_par%JOB = 4
mumps_par%ICNTL(4) = 0
CALL ZMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
        WRITE (0,*) 'error2 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
	ERROR STOP 200
ENDIF      	
			 


!allocate varabie
!---------------------------------------------------------------------------------
ALLOCATE(eigen_vec(NMF,ncv)) !<egienvector
ALLOCATE(d(NMF))     !<egienvalue
ALLOCATE(w1(NMQ))     !<RHS vector
ALLOCATE(w2(NMQ))     !<RHS vector
ALLOCATE(u(NMF))     !<solution vector u 

ALLOCATE(resid(NMF))
ALLOCATE(v(NMF,ncv))
ALLOCATE(workd(3*NMF))


!Using ARPACK to find the needed eigen values and eigen vectors
!-----------------------------------------------------------------------------
ido=0
tol=0.0
iparam(1)=1        ! exact shifts
iparam(3) = 1000000! max iterations
iparam(7) = 2      ! mode 2  regular inverse mode
sigma=(0.0,0.0)
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
		CALL znaupd ( ido, bmat, NMF, which, nev, tol, resid, ncv, &
				&     v, NMF, iparam, ipntr, workd, workl, lworkl,& 
				&     rwork, info )
!		write(0,*)'ido is',ido
	ENDIF
			
	! important (ido for other mumps_par%MYID .nq. 0)	
	call MPI_BCAST(ido,1,MPI_INT,0,MPI_COMM_WORLD,ierrM)
	
        
        IF (ido .eq. -1 .or. ido .eq. 1) then

!           %----------------------------------------%
!           | Perform  y <--- OP*x = inv[MF]*OA*x    |
!           | The user should supply his/her own     |
!           | matrix vector routine and a linear     |
!           | system solver.  The matrix-vector      |
!           | subroutine should take workd(ipntr(1)) |
!           | as input, and the final result should  |
!           | be returned to workd(ipntr(2)).        |
!           | OA=P'*R'*MQ*R*P, R=inv(A-omega*M)      |
!           %----------------------------------------%
          
               ! step1 
               ! Call matrix vector multiplication for P*v
         	
		CALL PV(w1,NMF,NMQ,NNZP,PIRN,PJCN,PA,workd(ipntr(1)))

               ! step2
               ! Call MUMPS for solution of (A-omega*M)*x=b for R*(P*v) 	
               	              		        
               ! load RHS		
		IF ( mumps_par%MYID .eq. 0 ) THEN
			mumps_par%RHS=(0.0,0.0)
			DO i = 1,NA
			   mumps_par%RHS(i) = w1(i)
			ENDDO
		ENDIF
		
		! solve Ax=b
		mumps_par%JOB = 3
		mumps_par%ICNTL(9) = 1		
		CALL ZMUMPS(mumps_par)
		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error4 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF    


		!  Solution has been assembled on the host
		IF ( mumps_par%MYID .eq. 0 ) THEN
		    w2=mumps_par%RHS
		END IF
				
		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error5 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF   
		

        	! step3
           	! Call matrix vector multiplication for MQ*(R*(P*v))
         	
		CALL MQV(w1,NMQ,NNZMQ,MQIRN,MQJCN,MQA,w2)
		
		
               ! step4       	
               ! Call MUMPS for solution of (A-omega*M)'*x=b for R'*(MQ*(R*(P*v))) 
                

               ! load RHS (x=INV(A'^T)*b-> x'=INV(A^T)*b'->x=(x')')	
		IF ( mumps_par%MYID .eq. 0 ) THEN
			mumps_par%RHS=(0.0,0.0)
			DO i = 1,NA
			   mumps_par%RHS(i) = CONJG(w1(i))
			ENDDO
		ENDIF
		
		! solve (A^T)*x'=b'		
		mumps_par%JOB = 3
		!(INCTL(9) = 0 ) for A^T
		mumps_par%ICNTL(9) = 0		
		CALL ZMUMPS(mumps_par)
		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error4 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF    

		!  Solution has been assembled on the host and x=(x')'
		IF ( mumps_par%MYID .eq. 0 ) THEN
		    w2=mumps_par%RHS
		    DO i = 1,NA
			   w2(i) = CONJG(w2(i))
		    ENDDO
		END IF
				
		IF (mumps_par%INFOG(1) <0 ) THEN
		    WRITE (0,*) ' error5 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
		    ERROR STOP 200
		ENDIF   
		      	                                                                                            
               ! step5
               ! Call matrix vector multiplication for P'*(R'*(MQ*(R*(P*v))))
         	
		CALL PHV(u,NMF,NMQ,NNZP,PHIRN,PHJCN,PHA,w2)               
                
                                
               ! step6
               ! Call SuperLU for solution of Mf*x=b for inv(Mf)*(R'*(MQ*(R*(P*v))))    
               
               
               ! Solve the system using the existing factors in SuperLU
		IF ( mumps_par%MYID .eq. 0 )THEN             
               	iopt = 2
               	CALL c_fortran_zgssv( iopt, NMF, NNZMF, nrhs, values, rowind, colptr, &
                        	 &   u, ldb, factors, infoSuperLu ) 
               ENDIF                  
               

        	CALL TRANSFORM(workd(ipntr(2)),NMF,u)        	
	        	
        	
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
      	
		CALL MFV(u,NMF,NNZMF,MFIRN,MFJCN,MFA,workd(ipntr(1)))
		
        	CALL TRANSFORM(workd(ipntr(2)),NMF,u)		

		! %-----------------------------------------%
		! | L O O P   B A C K to call ZNAUPD again. |
		! %-----------------------------------------%	
	
	ELSE
			
		EXIT
	ENDIF 
	
ENDDO
!********************************************************************************

	
IF ( mumps_par%MYID .eq. 0 ) THEN	
CALL zneupd (.TRUE., 'A', select, d, v, NMF, sigma, workev, &
           & bmat, NMF, which, nev, tol, resid, ncv, v, &
           & NMF, iparam, ipntr, workd, workl, lworkl, rwork, ierrA )
ENDIF                        
call MPI_BCAST(ido,1,MPI_INT,0,MPI_COMM_WORLD,ierrA)
   
IF( ierrA < 0 ) THEN
    ERROR STOP 106
ELSE
	eigen_vec(:,1:nev)=v(:,1:nev)
ENDIF

ALLOCATE(forcing_mode(1:NMF))
ALLOCATE(responding_mode(1:NMQ))

! output singular_vector (forcing mode)
IF ( mumps_par%MYID .eq. 0 ) THEN
    DO I = 1, NMF
        forcing_mode(I) = eigen_vec(I,1)
    ENDDO
    WRITE(filename,'(a,i0,a)')'forcing',ITE,'.txt'
    OPEN(10,file=trim(filename),status='replace')
    DO I = 1,NMF
        WRITE(10,*) forcing_mode(I)
    ENDDO
    CLOSE(10)    
END IF

! calculating responding mode q=R*P*f
! step 1 : P*f
CALL PV(responding_mode,NMF,NMQ,NNZP,PIRN,PJCN,PA,forcing_mode)
! step 2 : R*(P*f)
! Call MUMPS for solution of (A-omega*M)*x=b for R*(P*v) 	
! load RHS		
IF ( mumps_par%MYID .eq. 0 ) THEN
	mumps_par%RHS=(0.0,0.0)
	DO i = 1,NA
	     mumps_par%RHS(i) = responding_mode(i)
	ENDDO
ENDIF
! solve Ax=b
mumps_par%JOB = 3
mumps_par%ICNTL(9) = 1		
CALL ZMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
	WRITE (0,*) ' error4 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
	ERROR STOP 200
ENDIF    
!  Solution has been assembled on the host
IF ( mumps_par%MYID .eq. 0 ) THEN
	responding_mode=mumps_par%RHS
END IF

! output responding mode
IF ( mumps_par%MYID .eq. 0 ) THEN
    WRITE(filename,'(a,i0,a)')'responding',ITE,'.txt'
    OPEN(10,file=trim(filename),status='replace')
    DO I = 1,NMQ
        WRITE(10,*) responding_mode(I)
    ENDDO
    CLOSE(10)    
END IF
	
DEALLOCATE(forcing_mode)
DEALLOCATE(responding_mode)

  
IF ( mumps_par%MYID .eq. 0 ) THEN
     DO I = 1, nev
         singular_value(ITE,I) = d(I)
     ENDDO
ENDIF  
   

IF ( mumps_par%MYID .eq. 0 ) THEN
	write(0,*) 'singular_value',d(1),d(2)
	write(0,*) 'ITE',ITE	
	write(0,*) 'omega',omega
ENDIF

time_end = MPI_Wtime()
total_time=time_end - time_begin
IF  ( mumps_par%MYID .eq. 0 ) THEN
	write(0,*) 'total_time',total_time
ENDIF

!deallocate varabie
DEALLOCATE(eigen_vec) 
DEALLOCATE(d)     
DEALLOCATE(w1)
DEALLOCATE(w2)     
DEALLOCATE(u)     
DEALLOCATE(resid)
DEALLOCATE(v)
DEALLOCATE(workd)
                

!  Deallocate matrix
IF ( mumps_par%MYID .eq. 0 )THEN
	DEALLOCATE(LIRN)
	DEALLOCATE(LJCN)   
	DEALLOCATE(LA) 				
ENDIF

!deallocate matrix	
IF ( mumps_par%MYID .eq. 0 )THEN
	DEALLOCATE( mumps_par%IRN )
	DEALLOCATE( mumps_par%JCN )
	DEALLOCATE( mumps_par%A )
	DEALLOCATE( mumps_par%RHS )			
ENDIF		


ENDDO
!----------------------------------------------------------------------------------

! output singular_value
IF ( mumps_par%MYID .eq. 0 ) THEN
    OPEN(10,file="result_ARPACK.txt")
    DO I = 1, N_omega
    	DO J = 1, nev
    	    write(10,*) singular_value(I,J)
    	ENDDO
    ENDDO
END IF


! Destroy the instance (deallocate internal data structures)
mumps_par%JOB = -2
CALL ZMUMPS(mumps_par)
IF (mumps_par%INFOG(1) <0 ) THEN
    WRITE (0,*) 'error6 in Mumps Solver with INFO(1)=',mumps_par%INFOG(1),' and INFO(2)=',mumps_par%INFOG(2)
    ERROR STOP 200
ENDIF    
CALL MPI_FINALIZE(IERR)

!free the storage allocated inside SuperLU
IF ( mumps_par%MYID .eq. 0 )THEN
	iopt = 3
	CALL c_fortran_zgssv( iopt, NMF, NNZMF, nrhs, values, rowind, colptr, &
                     &    b, ldb, factors, infoSuperLu )
ENDIF     

END 
!=========================================================

!*********************************************************
SUBROUTINE TRANSFORM(x2,nsize,x1)

IMPLICIT NONE

INTEGER::nsize
COMPLEX*16 ::x1(nsize)    !< input
COMPLEX*16 ::x2(nsize)    !< solution vector 

INTEGER::i

DO i =1,nsize
   x2(i) = x1(i)
ENDDO
   

END SUBROUTINE TRANSFORM
!*********************************************************



!*********************************************************
SUBROUTINE MFV(x,NMF,NNZMF,MFIRN,MFJCN,MFA,u)

IMPLICIT NONE

INTEGER::NMF,NNZMF
COMPLEX*16 ::u(NMF)    !< input
COMPLEX*16 ::x(NMF)    !< solution vector 
INTEGER::MFIRN(NNZMF)
INTEGER::MFJCN(NNZMF)
REAL*16::MFA(NNZMF)

INTEGER::i

DO i =1,NMF
   x(i) = (0.0,0.0)
ENDDO

DO i =1, NNZMF
   x(MFIRN(i)) = x(MFIRN(i)) + MFA(i) * u(MFJCN(i))
ENDDO


END SUBROUTINE MFV
!*********************************************************



!*********************************************************
SUBROUTINE MQV(x,NMQ,NNZMQ,MQIRN,MQJCN,MQA,u)

IMPLICIT NONE

INTEGER::NMQ,NNZMQ
COMPLEX*16 ::u(NMQ)    !< input
COMPLEX*16 ::x(NMQ)    !< solution vector 
INTEGER::MQIRN(NNZMQ)
INTEGER::MQJCN(NNZMQ)
REAL*16::MQA(NNZMQ)

INTEGER::i

DO i =1,NMQ
   x(i) = (0.0,0.0)
ENDDO

DO i =1, NNZMQ
   x(MQIRN(i)) = x(MQIRN(i)) + MQA(i) * u(MQJCN(i))
ENDDO


END SUBROUTINE MQV
!*********************************************************


!*********************************************************
SUBROUTINE PV(x2,NMF,NMQ,NNZP,PIRN,PJCN,PA,x1)

IMPLICIT NONE

INTEGER::NMF,NMQ,NNZP
COMPLEX*16 ::x1(NMF)    !< input
COMPLEX*16 ::x2(NMQ)    !< solution vector 
INTEGER::PIRN(NNZP)
INTEGER::PJCN(NNZP)
REAL*16::PA(NNZP)

INTEGER::i

DO i =1,NMQ
   x2(i) = (0.0,0.0)
ENDDO

DO i =1, NNZP
   x2(PIRN(i)) = x2(PIRN(i)) + PA(i) * x1(PJCN(i))
ENDDO


END SUBROUTINE PV
!*********************************************************



!*********************************************************
SUBROUTINE PHV(x2,NMF,NMQ,NNZP,PHIRN,PHJCN,PHA,x1)

IMPLICIT NONE

INTEGER::NMF,NMQ,NNZP
COMPLEX*16 ::x1(NMQ)    !< input
COMPLEX*16 ::x2(NMF)    !< solution vector 
INTEGER::PHIRN(NNZP)
INTEGER::PHJCN(NNZP)
REAL*16::PHA(NNZP)

INTEGER::i

DO i =1,NMF
   x2(i) = (0.0,0.0)
ENDDO

DO i =1, NNZP
   x2(PHIRN(i)) = x2(PHIRN(i)) + PHA(i) * x1(PHJCN(i))
ENDDO


END SUBROUTINE PHV
!*********************************************************





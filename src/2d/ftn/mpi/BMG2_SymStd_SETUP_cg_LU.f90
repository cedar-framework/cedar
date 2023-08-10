      SUBROUTINE BMG2_SymStd_SETUP_cg_LU(&
     &                SO, II, JJ, NStncl, iGs, jGs, NGx, NGy,&
     &                ABD, NABD1, NABD2,&
     &                WS, NMSGr, NProcI, NProcJ, NProc, MyProc,&
     &                ProcGrid, ProcCoord, LocArrSize, MPICOMM)&
     &                BIND(C, NAME='MPI_BMG2_SymStd_SETUP_cg_LU')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_SETUP_cg_LU sets up the matrix on the coarsest grid,
!     and using the LAPACK routine DPBTRF, it forms the LU decomposition
!     of the matrix.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     SO        Refer to BMG2_SymStd_SOLVE_boxmg
!     II        Number of grid points in x direction, including
!               two fictitious points.
!     JJ        Number of grid points in y direction, including
!               two fictitious points.
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!     ABD       Refer to BMG2_SymStd_SOLVE_boxmg
!
! ======================================================================
!  --------------------
!   LOCAL:
!  --------------------
!
!
! ======================================================================
      USE ModInterface
      IMPLICIT NONE

! -----------------------------
!     Includes
!
      INCLUDE 'mpif.h'

      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_parameters_f90.h'
      INCLUDE 'BMG_constants_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), VALUE :: II, JJ, iGs, jGs, NABD1, NABD2,&
           nGx, nGy, NMSGr
      integer(c_int), VALUE :: NProcI, NProcJ, NProc, MyProc,&
           NStncl
      integer, VALUE :: MPICOMM
      integer(c_int) ::  ProcGrid(NProcI, NProcJ), ProcCoord(2,NProc)
      integer(len_t) :: LocArrSize(3,*)
      real(real_t) :: ABD(NABD1,NABD2), SO(II+1,JJ+1,NStncl), WS(NMSGr)

! ----------------------------
!     Local Declarations
!
      INTEGER  IBC, INFO, N, IERR, proc
      integer(len_t) :: I, I1, I2, J, J1, INT_TEMP, IIG, JJG, p_WS
      integer(len_t) :: KK, K, KKMAX, LARGESTNODES
      integer(len_t) :: INT_TEMP1, INT_TEMP2, P1, P2, P1SUM, P2SUM, PXP1, PYP2
      integer cg_comm_type

! ======================================================================

      IBC  = BMG_BCs_definite
      cg_comm_type = BMG_CG_ALLGATHER

! ------------------------------------------------
!     Calculate loop bounds
! ------------------------------------------------

      I1=II-1
      J1=JJ-1

      I2=I1-1
      N=I2*(J1-1)

! ------------------------------------------------
!     Calculate the global number of points
! ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      DO i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1)) - 2
      ENDDO

      IF (IIG.NE.NGx) THEN
         WRITE(*,*) 'What ...IIG = ', IIG, ' NotEqual NGx = ', NGx
      ENDIF
      !
      ! Global number in y
      !
      JJG=2
      DO j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j)) - 2
      ENDDO

      IF (JJG.NE.NGy) THEN
         WRITE(*,*) 'What ...JJG = ', JJG, ' NotEqual NGy = ', NGy
      ENDIF


! ------------------------------------------------
!     Find the largest local array
! ------------------------------------------------

      !
      ! Loop over local x dimensions
      !
      INT_TEMP = LocArrSize(1,1) - 2
      INT_TEMP1 = INT_TEMP
      DO proc=2, NProcI*NProcJ
         INT_TEMP = LocArrSize(1,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP1) THEN
            INT_TEMP1 = INT_TEMP
         END IF
      END DO

      !
      ! Loop over local y dimensions
      !
      INT_TEMP = LocArrSize(2,1) - 2
      INT_TEMP2 = INT_TEMP
      DO proc=2, NProcI*NProcJ
         INT_TEMP = LocArrSize(2,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP2) THEN
            INT_TEMP2 = INT_TEMP
         END IF
      END DO

      !
      ! Conservative: take largest from each
      !
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 * 5 + 3

! ------------------------------------------------
!     Copy all information into the buffer
! ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      INT_TEMP = 4

      IF ( NStncl.EQ.5 ) THEN
         !
         IF ( IBC.EQ.BMG_BCs_definite &
     &        .OR. IBC.EQ.BMG_BCs_indef_nonper ) THEN
            !
            ! Nonperiodic, but possibly indefinite
            !
            DO J = 2, JJ - 1
               DO I = 2, II - 1
                  WS(INT_TEMP)   =  SO(I,J,ko)
                  WS(INT_TEMP+1) = -SO(I,J,kw)
                  WS(INT_TEMP+2) = -SO(I+1,J,knw)
                  WS(INT_TEMP+3) = -SO(I,J,ks)
                  WS(INT_TEMP+4) = -SO(I,J,ksw)
                  INT_TEMP = INT_TEMP + 5
               END DO
            END DO
            !
            ! Indefinite ...
            !
            IF ( IBC.EQ.BMG_BCs_indef_nonper &
     &           .AND. iGs+I1.EQ.NGx .AND. jGs+J1.EQ.NGy ) THEN
               INT_TEMP = INT_TEMP - 5
               WS(INT_TEMP) = WS(INT_TEMP) + SO(I1,J1,KO)
            ENDIF
            !
         ELSE
            !
            IF ( MyProc.EQ.1 ) THEN
               WRITE(*,*) 'BMG2_SymStd_SETUP_cg_LU: Unsupport BC', IBC
            ENDIF
            CALL MPI_Finalize(MPICOMM)
            !
         ENDIF
      ELSE IF ( NStncl.EQ.3 ) THEN
         !
         IF ( IBC.EQ.BMG_BCs_definite &
     &        .OR. IBC.EQ.BMG_BCs_indef_nonper ) THEN
            !
            ! Nonperiodic, but possibly indefinite
            !
            DO J = 2, JJ - 1
               DO I = 2, II - 1
                  WS(INT_TEMP)   =  SO(I,J,ko)
                  WS(INT_TEMP+1) = -SO(I,J,kw)
                  WS(INT_TEMP+2) =  0.0D0
                  WS(INT_TEMP+3) = -SO(I,J,ks)
                  WS(INT_TEMP+4) =  0.0D0
                  INT_TEMP = INT_TEMP + 5
               END DO
            END DO
            !
            IF ( IBC.LT.0 ) THEN
               INT_TEMP = INT_TEMP - 5
               WS(INT_TEMP) = WS(INT_TEMP) + SO(I1,J1,KO)
            ENDIF
            !
         ELSE
            !
            IF ( MyProc.EQ.1 ) THEN
               WRITE(*,*) 'BMG2_SymStd_SETUP_cg_LU: Unsupport BC', IBC
            ENDIF
            CALL MPI_Finalize(MPICOMM)
            !
         ENDIF
      ELSE
         IF (MyProc.EQ.1) THEN
            WRITE(*,*)
            WRITE(*,*) 'NEED: NStncl = 3 or 5 '
            WRITE(*,*) 'HAVE: NStncl = ', NStncl
         END IF

         CALL print_error(C_CHAR_"CG solve failed"//C_NULL_CHAR)
         RETURN

      ENDIF

! ------------------------------------------------
!     Send/Receive information to/from everybody
! ------------------------------------------------

      IF (cg_comm_type .eq. BMG_CG_ALLGATHER) THEN

         CALL MPI_Allgather ( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        MPICOMM, IERR&
     &        )

      ELSE

         CALL MPI_Gather ( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        iZERO, MPICOMM, IERR&
     &        )

      END IF

! ------------------------------------------------
!     Assemble the global coarse grid matrix
! ------------------------------------------------


      IF ( cg_comm_type .eq. BMG_CG_ALLGATHER .OR.&
     &     MyProc .eq. iONE ) THEN


         KKMAX=0

         DO proc=1, NProcI*NProcJ

            INT_TEMP = proc * LARGESTNODES + 1
            INT_TEMP = WS(INT_TEMP)

            P1 = ProcCoord(1,INT_TEMP)
            P2 = ProcCoord(2,INT_TEMP)

            P1SUM = 0
            DO i=1, P1-1
               P1SUM = P1SUM + LocArrSize(1,ProcGrid(i,1)) - 2
            END DO

            P2SUM = 0
            DO j=1, P2-1
               P2SUM = P2SUM + LocArrSize(2,ProcGrid(1,j)) - 2
            END DO

            PXP1 = LocArrSize(1,INT_TEMP) - 2
            PYP2 = LocArrSize(2,INT_TEMP) - 2

            DO j=1, PYP2
               DO i=1, PXP1

                  KK = P1SUM + (i-1) + (IIG-2) * (P2SUM + j - 1) + 1
                  p_WS =  proc*LARGESTNODES+5*( (i-1)+(j-1)*PXP1) + 4

                  ABD(IIG,KK)   = WS( p_WS )
                  ABD(IIG-1,KK) = WS( p_WS + 1 )

                  ABD(3,KK) =  WS( p_WS + 2 )
                  ABD(2,KK) =  WS( p_WS + 3 )
                  ABD(1,KK) =  WS( p_WS + 4 )

                  IF (KK.gt.KKMAX) KKMAX = KK

               END DO
            END DO

         END DO



! ------------------------------------------------
!     DPBTRF is a LAPACK routine
! ------------------------------------------------

         CALL DPBTRF('U', KKMAX, IIG-1, ABD, NABD1, INFO)

         IF (INFO .NE. 0) THEN

            IF (MyProc.EQ.1) THEN
               call print_error(C_CHAR_"Coarse grid Cholesky decomposition failed!"//C_NULL_CHAR)
               WRITE(*,*) 'INFO = ', INFO
            END IF

            RETURN

         ENDIF

      END IF

! ======================================================================


! =============================

      RETURN
      END

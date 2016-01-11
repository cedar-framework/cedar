      SUBROUTINE BMG3_SymStd_SOLVE_cg_LU(&
     &                Q, QF, II, JJ, KK, ABD, BBD, NABD1, NABD2, NOGm,&
     &                NProcI, NProcJ, NProcK, NProc, MyProc,&
     &                ProcGrid, ProcCoord, LocArrSize,&
     &                WS, NMSGr, MPICOMM&
     &                ) BIND(C, NAME='BMG3_SymStd_SOLVE_cg_LU')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SOLVE_cg_LU does a direct solve on the coarsest grid.
!     It uses the LAPACK routine DPBTRS.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     QF        Refer to BMG3_SymStd_SOLVE_boxmg
!
!     II        Number of grid points in x direction, including
!               two fictitious points.
!     JJ        Number of grid points in y direction, including
!               two fictitious points.
!
!     ABD       Refer to BMG3_SymStd_SOLVE_boxmg
!
!     NABD1     Refer to BMG3_SymStd_SOLVE_boxmg
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!     BBD       Workspace, refer to BMG3_SymStd_SOLVE_boxmg
!
! ======================================================================
!  --------------------
!   OUTPUT:
!  --------------------
!
!     Q         Refer to BMG3_SymStd_SOLVE_boxmg
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

      INCLUDE 'mpif.h'
      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_parameters_f90.h'
      INCLUDE 'BMG_constants_f90.h'


! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: ii, jj, kk, nabd1, nabd2, NMSGr
      integer(c_int), value :: NOGm, NProc, NProcI, NProcJ, NProcK,&
           myproc, mpicomm
      integer(c_int) :: ProcGrid(NProcI, NProcJ, NProcK), ProcCoord(3,NProc)
      integer(len_t) :: LocArrSize(3,*)
      real(real_t) :: abd(nabd1,nabd2), bbd(nabd2), q(ii,jj,kk), qf(ii,jj,kk),&
           WS(NMSGr)

! ----------------------------
!     Local Declarations
!
      integer i, ibw, j, k, info
      INTEGER IIG, JJG, KKG
      INTEGER INT_TEMP, INT_TEMP1, INT_TEMP2, INT_TEMP3
      INTEGER KKMAX, iBEG, jBEG, kBEG, iEND, jEND, kEND
      INTEGER proc, P1, P2, P3, P1SUM, P2SUM, P3SUM
      INTEGER IERR, TMP_RANK, IIL, JJL, KKL, KK_t
      INTEGER p_WS, LARGESTNODES, MPI_IERROR
      INTEGER iMyProc, ictr
      INTEGER cg_comm_type

! ======================================================================

      ! hardcoded for now
      cg_comm_type = BMG_CG_ALLGATHER

! ------------------------------------------------
!     Calculate the global number of points
! ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      do i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1,1)) - 2
      enddo

      !
      ! Global number in y
      !
      JJG=2
      do j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j,1)) - 2
      enddo

      !
      ! Global number in z
      !
      KKG=2
      do k=1,NProcK
         KKG = KKG + LocArrSize(3,ProcGrid(1,1,k)) - 2
      enddo

! ------------------------------------------------
!     Find the largest local array
! ------------------------------------------------

      !
      ! Loop over local x dimensions
      !
      INT_TEMP = LocArrSize(1,1) - 2
      INT_TEMP1 = INT_TEMP
      DO proc=2, NProcI*NProcJ*NProcK
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
      DO proc=2, NProcI*NProcJ*NProcK
         INT_TEMP = LocArrSize(2,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP2) THEN
            INT_TEMP2 = INT_TEMP
         END IF
      END DO

      !
      ! Loop over local y dimensions
      !
      INT_TEMP = LocArrSize(3,1) - 2
      INT_TEMP3 = INT_TEMP
      DO proc=2, NProcI*NProcJ*NProcK
         INT_TEMP = LocArrSize(3,proc) - 2
         IF (INT_TEMP.gt.INT_TEMP3) THEN
            INT_TEMP3 = INT_TEMP
         END IF
      END DO


      !
      ! Conservative: take largest from each
      !
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 * INT_TEMP3  + 4

! ------------------------------------------------
!     Copy all information into the buffer
! ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      WS(4) = KK
      INT_TEMP = 5

      DO K = 2, KK - 1
         DO J = 2, JJ - 1
            DO I = 2, II - 1
               WS(INT_TEMP) = QF(I,J,K)
               INT_TEMP = INT_TEMP + 1
            END DO
         END DO
      END DO

! ------------------------------------------------
!     Send/Receive information to/from everybody
! ------------------------------------------------


      IF (cg_comm_type .eq. BMG_CG_ALLGATHER) THEN

         CALL MPI_ALLGATHER( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        MPICOMM, IERR&
     &        )

      ELSE

         CALL MPI_GATHER( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        iZERO, MPICOMM, IERR&
     &        )

      ENDIF

! ------------------------------------------------
!     Assemble the global right hand side
! ------------------------------------------------

      IF ( cg_comm_type .eq. BMG_CG_ALLGATHER .OR. &
     &     MyProc .eq. iONE ) THEN

         KKMAX=0

         DO proc = 1, NProcI * NProcJ * NProcK

            TMP_RANK = WS(proc*LARGESTNODES+1)

            P1 = ProcCoord(1,TMP_RANK)
            P2 = ProcCoord(2,TMP_RANK)
            P3 = ProcCoord(3,TMP_RANK)

            P1SUM = 0
 6          P2SUM = 0
            P3SUM = 0

            DO i=1, P1-1
               P1SUM = P1SUM + LocArrSize(1,ProcGrid(i,1,1)) - 2
            END DO

            DO j=1, P2-1
               P2SUM = P2SUM + LocArrSize(2,ProcGrid(1,j,1)) - 2
            END DO

            DO k=1, P3-1
               P3SUM = P3SUM + LocArrSize(3,ProcGrid(1,1,k)) - 2
            END DO

            IIL = LocArrSize(1,TMP_RANK) - 2
            JJL = LocArrSize(2,TMP_RANK) - 2
            KKL = LocArrSize(3,TMP_RANK) - 2

            DO k = 1, KKL
               DO j = 1, JJL
                  DO i = 1, IIL

                     KK_t = (P1SUM + i-1) + (IIG-2) * (P2SUM + j - 1) &
     &                    + (IIG-2)*(JJG-2)*(P3SUM + k - 1) + 1

                     p_WS =  proc*LARGESTNODES + 5 &
     &                    + ((i-1) + (j-1)*IIL + (k-1)*IIL*JJL)

                     BBD(KK_t) =  WS(p_WS)

                     IF (KK_t.gt.KKMAX) KKMAX = KK_t

                  END DO
               END DO
            END DO


         END DO

!     -------------------------------------------
!     solve the linear system
!     -------------------------------------------

         ibw=(IIG-2)*(JJG-1)+1

         CALL DPBTRS ('U', KKMAX, IBW, 1, ABD, NABD1, BBD, KKMAX, INFO)
         IF (INFO .NE. 0) THEN

            call print_error("Coarse grid solve failed!"//C_NULL_CHAR)
            WRITE(*,510) 'INFO = ', INFO

            RETURN

         ENDIF

      ENDIF

! ------------------------------------------------
!     Extract the local data
! ------------------------------------------------

      IF (cg_comm_type .eq. BMG_CG_GATHER_SCATTER) THEN

         !
         ! Broadcast the solution vector
         !

         CALL MPI_Bcast(BBD(1), NABD2, MPI_DOUBLE_PRECISION, &
     &        iZERO, MPICOMM, IERR)

      END IF


      KK_t=0

      IIL = LocArrSize(1,MyProc) - 2
      JJL = LocArrSize(2,MyProc) - 2
      KKL = LocArrSize(3,MyProc) - 2

      P1 = ProcCoord(1,MyProc)
      P2 = ProcCoord(2,MyProc)
      P3 = ProcCoord(3,MyProc)


      P1SUM = 0
      DO i=1, P1-1
         P1SUM = P1SUM + LocArrSize(1,ProcGrid(i,1,1)) - 2
      END DO

      P2SUM = 0
      DO j=1, P2-1
         P2SUM = P2SUM + LocArrSize(2,ProcGrid(1,j,1)) - 2
      END DO

      P3SUM = 0
      DO k=1, P3-1
         P3SUM = P3SUM + LocArrSize(3,ProcGrid(1,1,k)) - 2
      END DO



      !
      !  Setup loop boundaries in x
      !
      IF ( P1.EQ.1 ) THEN
         iBEG = 1
      ELSE
         iBEG = 0
      END IF

      IF (P1.EQ.NProcI) THEN
         iEND = IIL
      ELSE
         iEND = IIL+1
      END IF

      !
      !  Setup loop boundaries in y
      !
      IF ( P2.EQ.1 ) THEN
         jBEG = 1
      ELSE
         jBEG = 0
      END IF

      IF ( P2.EQ.NProcJ) THEN
         jEND = JJL
      ELSE
         jEND = JJL + 1
      ENDIF

      !
      !  Setup loop boundaries in z
      !
      IF ( P3.EQ.1 ) THEN
         kBEG = 1
      ELSE
         kBEG = 0
      END IF

      IF ( P3.EQ.NProcK) THEN
         kEND = KKL
      ELSE
         kEND = KKL + 1
      ENDIF



      DO k=kBEG, kEND
         DO j=jBEG, jEND
            DO i=iBEG, iEND

               KK_t = (P1SUM + i-1) + (IIG-2) * (P2SUM + j - 1) +&
     &              (IIG-2)*(JJG-2) * (P3SUM + k - 1) + 1

               Q(i+1,j+1,k+1)=BBD(KK_t)

            END DO
         END DO
      END DO

! ======================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SOLVE_cg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

! ===========================================

      RETURN
      END

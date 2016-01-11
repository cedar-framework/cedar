      SUBROUTINE BMG3_SymStd_SETUP_cg_LU(&
     &                SO, II, JJ, KK, NStncl, ABD, NABD1, NABD2,&
     &                WS, NMSGr, NProcI, NProcJ, NProcK, NProc, MyProc,&
     &                ProcGrid, ProcCoord, LocArrSize, MPICOMM&
     &                ) BIND(C, NAME='MPI_BMG3_SymStd_SETUP_cg_LU')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SETUP_cg_LU sets up the matrix on the coarsest grid,
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
!     KK        Number of grid points in z direction, including
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
!     ABD       Refer to BMG3_SymStd_SOLVE_boxmg
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

! ----------------------------
!     Includes
!
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

      INCLUDE 'mpif.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: ii, jj, kk, nabd1, nabd2, NMSGr
      integer(c_int), value :: NStncl, NProcI, NProcj, NProcK, NProc,&
           MyProc, MPICOMM
      integer(c_int) :: ProcGrid(NProcI,NProcJ,NProcK), ProcCoord(3, NProc)
      integer(len_t) :: LocArrSize(3,*)
      real(real_t) :: abd(nabd1,nabd2), so(ii+1,jj+1,kk+1,NStncl)
      real(real_t) :: WS(NMSGr)

! ----------------------------
!     Local Declarations
!
      INTEGER i, i1, i2, ibw, j, j1, k, k1, kl, KK_t, n
      INTEGER IIG, JJG, KKG, proc, KKMAX, P1, P2, P3
      INTEGER P1SUM, P2SUM, P3SUM, PXP1, PYP2, PZP3, p_WS
      INTEGER INT_TEMP, INT_TEMP1, INT_TEMP2, INT_TEMP3
      INTEGER LARGESTNODES, IERR, MPI_IERROR
      INTEGER cg_comm_type, info

! ======================================================================

      ! hardcoded for now
      cg_comm_type = BMG_CG_ALLGATHER
!
!     set up and factor matrix for direct solve on coarsest grid.
!
      info = 0

      i1=ii-1
      j1=jj-1
      k1=kk-1

      i2=i1-1
!      n=i2*(j1-1)*(kk-2)

! ------------------------------------------------
!     Calculate the global number of points
! ------------------------------------------------

      !
      ! Global number in x
      !
      IIG=2
      DO i=1,NProcI
         IIG = IIG + LocArrSize(1,ProcGrid(i,1,1)) - 2
      ENDDO

      !
      ! Global number in y
      !
      JJG=2
      DO j=1,NProcJ
         JJG = JJG + LocArrSize(2,ProcGrid(1,j,1)) - 2
      ENDDO


      !
      ! Global number in z
      !
      KKG=2
      DO k=1,NProcK
         KKG = KKG + LocArrSize(3,ProcGrid(1,1,k)) - 2
      ENDDO

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
      ! Loop over local z dimensions
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
      LARGESTNODES = INT_TEMP1 * INT_TEMP2 * INT_TEMP3 * 14 + 4

! ------------------------------------------------
!     Copy all information into the buffer
! ------------------------------------------------

      WS(1) = MyProc
      WS(2) = II
      WS(3) = JJ
      WS(4) = KK
      INT_TEMP = 5

      IF ( NStncl.EQ.14 ) THEN
         !
         DO K = 2, KK - 1
            DO J = 2, JJ - 1
               DO I = 2, II - 1
                  WS(INT_TEMP)    =  SO(I  ,J  ,k  ,kp)
                  WS(INT_TEMP+ 1) = -SO(I  ,J  ,k  ,kpw)
                  WS(INT_TEMP+ 2) = -SO(I+1,J  ,k  ,kpnw)
                  WS(INT_TEMP+ 3) = -SO(I  ,J  ,k  ,kps)
                  WS(INT_TEMP+ 4) = -SO(I  ,J  ,k  ,kpsw)
                  WS(INT_TEMP+ 5) = -SO(I+1,J+1,k  ,kbne)
                  WS(INT_TEMP+ 6) = -SO(I  ,J+1,k  ,kbn)
                  WS(INT_TEMP+ 7) = -SO(I  ,J+1,k  ,kbnw)
                  WS(INT_TEMP+ 8) = -SO(I+1,J  ,k  ,kbe)
                  WS(INT_TEMP+ 9) = -SO(I  ,J  ,k  ,kb)
                  WS(INT_TEMP+10) = -SO(I  ,J  ,k  ,kbw)
                  WS(INT_TEMP+11) = -SO(I+1,J  ,k  ,kbse)
                  WS(INT_TEMP+12) = -SO(I  ,J  ,k  ,kbs)
                  WS(INT_TEMP+13) = -SO(I  ,J  ,k  ,kbsw)
                  INT_TEMP = INT_TEMP + 14
               END DO
            END DO
         END DO
         !
      ELSE IF ( NStncl.EQ.4 ) THEN
         !
         DO K = 2, KK - 1
            DO J = 2, JJ - 1
               DO I = 2, II - 1
                  WS(INT_TEMP)    =  SO(I  ,J  ,k  ,kp)
                  WS(INT_TEMP+ 1) = -SO(I  ,J  ,k  ,kpw)
                  WS(INT_TEMP+ 2) = rZERO
                  WS(INT_TEMP+ 3) = -SO(I  ,J  ,k  ,kps)
                  WS(INT_TEMP+ 4) = rZERO
                  WS(INT_TEMP+ 5) = rZERO
                  WS(INT_TEMP+ 6) = rZERO
                  WS(INT_TEMP+ 7) = rZERO
                  WS(INT_TEMP+ 8) = rZERO
                  WS(INT_TEMP+ 9) = -SO(I  ,J  ,k  ,kb)
                  WS(INT_TEMP+10) = rZERO
                  WS(INT_TEMP+11) = rZERO
                  WS(INT_TEMP+12) = rZERO
                  WS(INT_TEMP+13) = rZERO
                  INT_TEMP = INT_TEMP + 14
               END DO
            END DO
         END DO
         !
      ELSE
         call print_error(C_CHAR_"Cholesky failed! (incorrect NStncl)"//C_NULL_CHAR)
         RETURN
      ENDIF


! ------------------------------------------------
!     Send/Receive information to/from everybody
! ------------------------------------------------

      IF (cg_comm_type .eq. BMG_CG_ALLGATHER) THEN

         CALL MPI_ALLGATHER ( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        MPICOMM, IERR&
     &        )

      ELSE

         CALL MPI_GATHER ( &
     &        WS(1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        WS(LARGESTNODES+1), LARGESTNODES, MPI_DOUBLE_PRECISION,&
     &        iZERO, MPICOMM, IERR&
     &        )

      ENDIF

!     ------------------------------------------------
!     Assemble the global coarse grid matrix
!     ------------------------------------------------


      IF ( cg_comm_type .eq. BMG_CG_ALLGATHER .OR.&
     &     MyProc .eq. iONE ) THEN

         KKMAX=0

         DO proc=1, NProcI*NProcJ*NProcK

            INT_TEMP = proc * LARGESTNODES + 1
            INT_TEMP = WS(INT_TEMP)

            P1 = ProcCoord(1,INT_TEMP)
            P2 = ProcCoord(2,INT_TEMP)
            P3 = ProcCoord(3,INT_TEMP)

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

            PXP1 = LocArrSize(1,INT_TEMP) - 2
            PYP2 = LocArrSize(2,INT_TEMP) - 2
            PZP3 = LocArrSize(3,INT_TEMP) - 2

            IBW = (IIG-2)*(JJG-1)+1

            DO k=1, PZP3
               DO j=1, PYP2
                  DO i=1, PXP1

                     KK_t = (P1SUM + i - 1) + (IIG-2) * (P2SUM + j - 1) &
     &                    + (IIG-2)*(JJG-2)*(P3SUM + k - 1) + 1

                     p_WS = proc*LARGESTNODES + 5 &
     &                    + 14*((i-1) + (j-1)*PXP1 + (k-1)*PXP1*PYP2)

                     ABD(IBW+1,KK_t)                 = WS( p_WS ) ! kp
                     ABD(IBW  ,KK_t)                 = WS( p_WS + 1 ) !

                     ABD(IBW-IIG+4,KK_t)             = WS( p_WS + 2 ) !
                     ABD(IBW-IIG+3,KK_t)             = WS( p_WS + 3 ) !
                     ABD(IBW-IIG+2,KK_t)             = WS( p_WS + 4 ) !

                     ABD(IBW-(JJG-3)*(IIG-2)+2,KK_t) = WS( p_WS + 5 ) !
                     ABD(IBW-(JJG-3)*(IIG-2)+1,KK_t) = WS( p_WS + 6 ) !
                     ABD(IBW-(JJG-3)*(IIG-2)  ,KK_t) = WS( p_WS + 7 ) !

                     ABD(IBW-(JJG-2)*(IIG-2)+2,kk_t) = WS( p_WS + 8 ) !
                     ABD(IBW-(JJG-2)*(IIG-2)+1,kk_t) = WS( p_WS + 9 ) !
                     ABD(IBW-(JJG-2)*(IIG-2)  ,kk_t) = WS( p_WS + 10 ) !

                     ABD(3,KK_t)                     = WS( p_WS + 11 ) !
                     ABD(2,KK_t)                     = WS( p_WS + 12 ) !
                     ABD(1,KK_t)                     = WS( p_WS + 13 ) !

                     IF (KK_t.gt.KKMAX) KKMAX = KK_t

                  END DO
               END DO
            END DO

         END DO

         CALL DPBTRF('U', KKMAX, IBW, ABD, NABD1, INFO)

         IF (INFO .NE. 0) THEN

            call print_error("Coarse grid Cholesky decomp failed!"//C_NULL_CHAR)

            if (MyProc .eq. 1) then
               WRITE(*,*) 'INFO = ', INFO
            endif

         END IF

      END IF

! ======================================================================

 150  FORMAT(I4,4X,11(E10.4,4X))

! ===========================================

      RETURN
      END

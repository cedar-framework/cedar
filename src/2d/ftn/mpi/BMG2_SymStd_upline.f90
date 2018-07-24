      SUBROUTINE BMG2_SymStd_upline(SOR,Q,II,JJ,&
           JBEG, RWORK, gwork, ngwork, gptr,&
           iface, niface, NPts, NLines,&
           K, NOLX, XCOMM, NSOR, TDG, &
           NMSGr, NOG, TDSX_SOR_PTRS,&
           SIZE, factorize, mp)

      use iso_c_binding, only: c_bool
      use message_passing
      use ModInterface
      IMPLICIT NONE

      INCLUDE 'mpif.h'

      INTEGER II, JJ, JBEG, NPts, NLines

      INTEGER K, NOLX, NOG, NSOR, NMSGr
      integer(c_int), value :: ngwork, niface

      INTEGER XCOMM(2,2:NOLX)
      integer(c_int) :: gptr(2,2:NOLX)

      INTEGER TDSX_SOR_PTRS(NOLX)

      REAL*8 SOR(II,JJ,2), RWORK(NMSGr)
      real(real_t) :: iface(niface), gwork(ngwork)

      REAL*8 Q(II,JJ), TDG(NSOR)

      INTEGER SIZE, N, My_L1_Rank, My_L2_Rank, IERR, kl

      INTEGER I, J, AT, FT, MULT, grank

      logical(c_bool) :: factorize
      type(c_ptr) :: mp

      integer :: icolor

! =============================================================
! Loop over all lower levels
! =============================================================

      icolor = mod(jbeg, 2) + 1
      CALL MPI_COMM_RANK(MPI_COMM_WORLD, grank, IERR)

      !
      ! Get my coarse level rank
      !
      IF (XCOMM(2,2) .NE. MPI_COMM_NULL) THEN
         CALL MPI_COMM_RANK(XCOMM(2,2), My_L1_Rank, IERR)
      END IF

      DO kl = 2, NOLX-1

         IF (XCOMM(2,kl+1) .NE. MPI_COMM_NULL) THEN

            !
            ! Get the size of the fine level process group
            !

            CALL MPI_COMM_SIZE(XCOMM(2,kl+1), SIZE, IERR)

            !
            ! Scatter Solution to processes in this subgroup
            !
            IF (XCOMM(2,kl) .NE. MPI_COMM_NULL) THEN
               call cedar_scatter(mp, gwork(gptr(icolor, kl) + 1),&
                    NLINES*8, MPI_DOUBLE_PRECISION,&
                    iface, NLINES*8, MPI_DOUBLE_PRECISION,&
                    0, XCOMM(2,kl),IERR)

            !
            ! Send interface solution from this level to
            ! LineSolve_C to compute solution.
            !
               CALL C_Wrapper(&
     &              TDG(TDSX_SOR_PTRS(kl)),&
     &              iface, rwork, SIZE,&
     &              NLines, NMSGr,&
     &              My_L1_Rank, factorize)

            END IF

            !
            ! Get my fine level rank
            !
            CALL MPI_COMM_RANK(XCOMM(2,kl+1), My_L2_Rank, IERR)

            !
            ! Copy Solution from TDG(:,:,4) into RWORK
            ! for use on the next level.
            !
            IF (XCOMM(2,kl) .NE. MPI_COMM_NULL) THEN
               CALL BMG2_SymStd_TDG_2_RWork(&
     &              TDG(TDSX_SOR_PTRS(kl)),&
     &              gwork(gptr(icolor,kl+1) + 1), NLines, SIZE, &
     &              My_L1_Rank, My_L2_Rank)
            END IF

            My_L1_Rank = My_L2_Rank

         END IF

      END DO

!     CALL MPI_FINALIZE(IERR)
!     STOP

      !
      ! Get ranks and size for fine process level
      !
      CALL MPI_COMM_RANK(XCOMM(2,NOLX), My_L2_Rank, IERR)
      CALL MPI_COMM_SIZE(XCOMM(2,NOLX), SIZE, IERR)

      !
      ! Scatter solution to next-to-finest level interface
      ! system
      !

      call cedar_scatter(mp, gwork(gptr(icolor,NOLX)+1),NLINES*8,MPI_DOUBLE_PRECISION,&
           iface, NLINES*8, MPI_DOUBLE_PRECISION,&
           0, XCOMM(2,NOLX), IERR)

      !
      ! Pointers into RWORK
      !
      MULT = 0
      if (factorize) then
         DO J=JBEG,JJ-1,2
            CALL BMG2_SymStd_LineSolve_C_ml(SOR(2,J,1),&
                 &        SOR(2,J,2), Q(1,J),&
                 &        iface(MULT*8 + 1),&
                 &        Npts, SIZE,&
                 &        My_L1_Rank)

            MULT = MULT+1
         END DO
      else
         DO J=JBEG,JJ-1,2
            CALL BMG2_SymStd_LineSolve_C_ml_eff(&
                 SOR(2,J,1), SOR(3,J,2), SOR(2,J,2),&
                 Q(1,J), iface(MULT*8 + 1),&
                 RWORK,&
                 Npts, SIZE,&
                 My_L1_Rank)

            MULT = MULT+1
         END DO
      endif

      RETURN
      END

! =================================================================

      SUBROUTINE C_Wrapper(TDG, iface, RWORK, SIZE, &
     &                     NLines, NMSGr, MyRank, factorize)
      use iso_c_binding, only: c_bool
      IMPLICIT NONE

      INTEGER JBEG, JJ, SIZE, NLines, NMSGr, MyRank

      REAL*8 TDG(2*SIZE+2,NLines,4)
      REAL*8 RWORK(NMSGr)
      real*8 iface(nlines*8)
      logical(c_bool) :: factorize

      INTEGER MULT, J, I

      MULT = 0
      if (factorize) then
         DO J=1,NLines
            CALL BMG2_SymStd_LineSolve_C_ml(TDG(2,J,1),&
                 &                    TDG(2,J,2), TDG(1,J,4), &
                 &                    iface(MULT*8 + 1), &
                 &                    2*SIZE, SIZE,&
                 &                    MyRank)

            MULT = MULT + 1
         END DO
      else
         DO J=1,NLines
            CALL BMG2_SymStd_LineSolve_C_ml_eff(&
                 TDG(2,J,1), TDG(3,J,2), TDG(2,J,2), TDG(1,J,4), &
                 iface(MULT*8 + 1), &
                 RWORK, &
                 2*SIZE, SIZE, MyRank)
            MULT = MULT + 1
         END DO
      endif

      RETURN
      END

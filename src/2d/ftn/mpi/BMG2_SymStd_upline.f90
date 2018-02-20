      SUBROUTINE BMG2_SymStd_upline(SOR,Q,II,JJ,&
     &                JBEG, RWORK, NPts, NLines, &
     &                K, NOLX, XCOMM, NSOR, TDG, &
     &                NMSGr, NOG, TDSX_SOR_PTRS,&
     &                SIZE, shm_enabled, shm_buff,&
     &                shm_len, shm_win)

      use iso_c_binding, only: c_bool, c_int
      IMPLICIT NONE

      INCLUDE 'mpif.h'

      INTEGER II, JJ, JBEG, NPts, NLines

      INTEGER K, NOLX, NOG, NSOR, NMSGr

      INTEGER XCOMM(2,2:NOLX)

      INTEGER TDSX_SOR_PTRS(NOLX)

      REAL*8 SOR(II,JJ,2), RWORK(NMSGr)

      REAL*8 Q(II,JJ), TDG(NSOR)

      INTEGER SIZE, N, My_L1_Rank, My_L2_Rank, IERR, kl

      INTEGER I, J, AT, FT, MULT, grank, idx
      integer shm_win
      integer(c_int) :: shm_len
      real*8 shm_buff(shm_len)
      logical(c_bool) :: shm_enabled

! =============================================================
! Loop over all lower levels
! =============================================================

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
               if (My_L1_Rank .eq. 0) then
                  CALL MPI_SCATTER(RWORK(My_L1_Rank*NLINES*8+1),&
                       &              NLINES*8,&
                       &              MPI_DOUBLE_PRECISION,MPI_IN_PLACE, NLINES*8, &
                       &              MPI_DOUBLE_PRECISION,0,XCOMM(2,kl),IERR)
               else
                  CALL MPI_SCATTER(0.0d0,&
                       &              NLINES*8,&
                       &              MPI_DOUBLE_PRECISION,RWORK, NLINES*8, &
                       &              MPI_DOUBLE_PRECISION,0,XCOMM(2,kl),IERR)
               endif

            !
            ! Send interface solution from this level to
            ! LineSolve_C to compute solution.
            !
               CALL C_Wrapper(&
     &              TDG(TDSX_SOR_PTRS(kl)),&
     &              RWORK, SIZE,&
     &              NLines, NMSGr,&
     &              My_L1_Rank)

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
     &              RWORK, NLines, SIZE, &
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

      if (shm_enabled) then
         if (My_L2_Rank .eq. 0) then
            call MPI_Win_lock_all(MPI_MODE_NOCHECK, shm_win, ierr)
            do k=0, SIZE-1
               do i=1, nlines
                  do j=1, 8
                     idx = k * (nlines*8) + (i*8 + j) + 1
                     shm_buff(idx) = rwork(idx)
                  enddo
               enddo
            enddo
            call MPI_Win_sync(shm_win, ierr)
         endif
         call MPI_Barrier(XCOMM(2,NOLX), ierr)
      else
         if (My_L2_Rank .eq. 0) then
            CALL MPI_SCATTER(RWORK(My_L2_Rank*NLINES*8+1),NLINES*8,&
                 &     MPI_DOUBLE_PRECISION,MPI_IN_PLACE, NLINES*8, &
                 &     MPI_DOUBLE_PRECISION,0,XCOMM(2,NOLX),IERR)
         else
            CALL MPI_SCATTER(0.0d0,NLINES*8,&
                 &     MPI_DOUBLE_PRECISION,RWORK, NLINES*8, &
                 &     MPI_DOUBLE_PRECISION,0,XCOMM(2,NOLX),IERR)
         endif
      endif

      ! call ml_relax_shm_up(ml_obj, rwork)

      !
      ! Pointers into RWORK
      !

      if (shm_enabled) then
         MULT = 0
         DO J=JBEG,JJ-1,2
            !MB      DO J=JBEG,JBEG  ! previous line was commented out
            idx = My_L2_Rank * nlines * 8 + MULT * 8 + 1
            CALL BMG2_SymStd_LineSolve_C_ml(SOR(2,J,1),&
                 &        SOR(2,J,2), Q(1,J),&
                 &        shm_buff(idx),&
                 &        Npts, SIZE,&
                 &        My_L1_Rank)

            MULT = MULT+1

         END DO
      else
         MULT = 0
         DO J=JBEG,JJ-1,2
            !MB      DO J=JBEG,JBEG  ! previous line was commented out

            CALL BMG2_SymStd_LineSolve_C_ml(SOR(2,J,1),&
                 &        SOR(2,J,2), Q(1,J),&
                 &        RWORK(MULT*8 + 1),&
                 &        Npts, SIZE,&
                 &        My_L1_Rank)

            MULT = MULT+1

         END DO
      endif

      RETURN
      END

! =================================================================

      SUBROUTINE C_Wrapper(TDG, RWORK, SIZE, &
     &                     NLines, NMSGr, MyRank)

      IMPLICIT NONE

      INTEGER JBEG, JJ, SIZE, NLines, NMSGr, MyRank

      REAL*8 TDG(2*SIZE+2,NLines,4)
      REAL*8 RWORK(NMSGr)

      INTEGER MULT, J, I

      MULT = 0
      DO J=1,NLines
!MB      DO J=1,1  ! previous line was commented out
         CALL BMG2_SymStd_LineSolve_C_ml(TDG(2,J,1),&
     &                    TDG(2,J,2), TDG(1,J,4), &
     &                    RWORK(MULT*8 + 1), &
     &                    2*SIZE, SIZE,&
     &                    MyRank)

         MULT = MULT + 1

      END DO

      RETURN
      END

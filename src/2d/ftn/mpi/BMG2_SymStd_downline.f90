      SUBROUTINE BMG2_SymStd_downline(SOR,Q,II,JJ,&
     &                JBEG, RWORK, NPts, NLines, &
     &                K, NOLX, XCOMM, NSOR, TDG, &
     &                NMSGr, NOG, TDSX_SOR_PTRS,&
     &                pgSIZE, fact_flags)

      use iso_c_binding, only: c_bool
      use ModInterface
      IMPLICIT NONE

      INCLUDE 'mpif.h'

      INTEGER II, JJ, JBEG, NPts, NLines

      INTEGER K, NOLX, NOG, NSOR, NMSGr

      INTEGER XCOMM(2,2:NOLX)

      INTEGER TDSX_SOR_PTRS(NOLX)

      REAL*8 SOR(II,JJ,2), RWORK(NMSGr)

      REAL*8 Q(II,JJ), TDG(NSOR)

      INTEGER pgSIZE, N, MyRank, IERR, kl

      INTEGER I, J, MULT, grank

      REAL*8 D, OD
      integer :: flag_stride
      logical(c_bool) :: fact_flags(2 * NOG)
      ! Flag for the interface system
      logical(c_bool) :: inter_flag

      ! always factorize interface systems
      inter_flag = .true.
      if (jbeg .eq. 2) then
         flag_stride = 0
      else
         flag_stride = NOLX
      endif

      CALL MPI_COMM_RANK(XCOMM(1,NOLX), grank, IERR)

! =======================================================
! Initialize IDs and pgSIZE
! =======================================================

      CALL MPI_COMM_SIZE(XCOMM(2,NOLX), pgSIZE, IERR)
      CALL MPI_COMM_RANK(XCOMM(2,NOLX), MyRank, IERR)

! =======================================================
! Do first step in downward pass
! =======================================================

      !
      ! Compute local interface system
      !
      MULT = 0
      NLines = 0
      DO J=JBEG, JJ-1, 2

         CALL BMG2_SymStd_LineSolve_A_ml(SOR(2,J,1),&
     &             SOR(2,J,2), Q(2,J),&
     &             RWORK(MULT*8+1),&
     &             NPts, fact_flags(K+flag_stride))

         MULT = MULT + 1

      END DO



      NLines = MULT

      !
      ! Gather interface equations to head node
      !
      IF (pgSIZE .GT. 1) THEN
         call ftimer_begin("comm-inter");
         if (myrank .eq. 0) then
            CALL MPI_GATHER(MPI_IN_PLACE, NLines*8, &
                 &        MPI_DOUBLE_PRECISION,&
                 &        RWORK, NLines*8, &
                 &        MPI_DOUBLE_PRECISION,&
                 &        0, XCOMM(2,NOLX), IERR)
         else
            CALL MPI_GATHER(RWORK, NLines*8, &
                 &        MPI_DOUBLE_PRECISION,&
                 &        0.0d0, NLines*8, &
                 &        MPI_DOUBLE_PRECISION,&
                 &        0, XCOMM(2,NOLX), IERR)
         endif
         call ftimer_end("comm-inter");
      END IF

!     CALL DUMP_RWORK(RWORK, 0, pgSIZE, NLines)

!     CALL MPI_FINALIZE(IERR)
!     STOP

! =======================================================
! Loop over the lower levels and coarsen
! =======================================================

      DO kl = NOLX-1, 2, -1

         IF (MyRank .EQ. 0) THEN

            !
            ! If I'm the head node, copy inteface
            ! system into TDG data structure.
            !

            CALL BMG2_SymStd_RWork_2_TDG(&
     &           TDG(TDSX_SOR_PTRS(kl)),&
     &           RWORK, pgSIZE, NLines, N,&
     &           K, kl)

!           write(*,*) ''
!           CALL dump_TDG(
!    &           TDG(TDSX_SOR_PTRS(kl,K)),
!    &           pgSIZE, grank, NLines)

            ! Form the local interface system
            ! for this level
            !
            CALL A_Wrapper(&
     &           TDG(TDSX_SOR_PTRS(kl)),&
     &           RWORK, N, pgSIZE, &
     &           JBEG, JJ, NLines, inter_flag)


            !
            ! Get the pgSIZE and MyRank for this level
            !
            CALL MPI_COMM_SIZE(XCOMM(2,kl), pgSIZE, IERR)
            CALL MPI_COMM_RANK(XCOMM(2,kl), MyRank, IERR)

            !
            ! If there is more than one process in
            ! this group then gather everything to
            ! the group head node
            !
            IF (pgSIZE .GT. 1) THEN
               call ftimer_begin("comm-inter");
               if (myrank .eq. 0) then
                  CALL MPI_GATHER(MPI_IN_PLACE, NLines*8, &
                       &              MPI_DOUBLE_PRECISION,&
                       &              RWORK, NLines*8, &
                       &              MPI_DOUBLE_PRECISION,&
                       &              0, XCOMM(2,kl), IERR)
               else
                  CALL MPI_GATHER(RWORK, NLines*8, &
                       &              MPI_DOUBLE_PRECISION,&
                       &              0.0d0, NLines*8, &
                       &              MPI_DOUBLE_PRECISION,&
                       &              0, XCOMM(2,kl), IERR)
               endif
               call ftimer_end("comm-inter");
            END IF

!           write(*,*) ''
!           CALL DUMP_RWORK(RWORK, 0, pgSIZE, NLines)

         END IF

!        CALL MPI_FINALIZE(IERR)
!        STOP

      END DO

      fact_flags(k+flag_stride) = .false.


      RETURN
      END

! ============================================================

      SUBROUTINE A_Wrapper(TDG, RWORK, Npts, &
     &                     pgSIZE, JBEG, JJ, NLines, FLG)

         use iso_c_binding, only: c_bool
         IMPLICIT NONE

         INTEGER JBEG, JJ, Npts, NLines, pgSIZE

         REAL*8 TDG(2*pgSIZE+2,NLines,4)
         REAL*8 RWORK(8*NLines)
         logical(c_bool) :: FLG

         INTEGER MULT, J

         MULT = 0
         DO J=1, NLines
            CALL BMG2_SymStd_LineSolve_A_ml(TDG(2,J,1),&
     &                TDG(2,J,2), TDG(2,J,4),&
     &                RWORK(MULT*8+1),&
     &                NPts, FLG)

            MULT = MULT + 1
         END DO

      RETURN
      END

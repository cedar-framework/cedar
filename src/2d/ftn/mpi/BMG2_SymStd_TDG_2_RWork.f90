      SUBROUTINE BMG2_SymStd_TDG_2_RWork(&
     &                TDG, RWORK, NLines,&
     &                SIZE, My_L1_Rank, &
     &                My_L2_Rank)

      IMPLICIT NONE

      INCLUDE 'mpif.h'


      INTEGER NLines, NPts, SIZE
      INTEGER My_L1_Rank, My_L2_Rank

      REAL*8 TDG(2*SIZE+2, NLines, 4)
      REAL*8 RWORK(8*NLines*SIZE)

      INTEGER LN, I, J, LOS, POS, PINCR, LINCR

      INTEGER grank, IERR

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,grank,IERR)

      PINCR = 8*NLines
      LINCR = 8

      LOS = 0

      DO J=1, NLines
         POS = 0
         DO I=0,SIZE-1

            RWORK(LOS+POS+1) = TDG(2*I + 1,J,4)
            RWORK(LOS+POS+2) = TDG(2*I + 2,J,4)
            RWORK(LOS+POS+3) = TDG(2*I + 3,J,4)
            RWORK(LOS+POS+4) = TDG(2*I + 4,J,4)

            POS = POS + PINCR

         END DO

         LOS = LOS + LINCR

      END DO

      RETURN
      END

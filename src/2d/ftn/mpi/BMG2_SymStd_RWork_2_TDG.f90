      SUBROUTINE BMG2_SymStd_RWork_2_TDG(&
     &                      LSOR, BUFF, SIZE,&
     &                      NLines, N, K, NOL)

      IMPLICIT NONE

      include 'mpif.h'

      INTEGER SIZE, NLines, K

      REAL*8 LSOR(2*SIZE+2,NLines,4)
      REAL*8 BUFF(SIZE*NLines*8)

      INTEGER N, POS, LOS, PINCR, LINCR

      INTEGER ii, jj

      INTEGER IERR, NOL, grank

      CALL MPI_COMM_RANK(MPI_COMM_WORLD,grank,IERR)

      PINCR = 8*NLines
      LINCR = 8

      LOS = 0

      DO jj=1,NLines
         POS = 0
         N = 1
         DO ii=1,SIZE

            N = N + 1

            LSOR(N,jj,1) = BUFF(LOS+POS+2)
            LSOR(N,jj,2) = BUFF(LOS+POS+1)
            LSOR(N,jj,3) = BUFF(LOS+POS+3)
            LSOR(N,jj,4) = BUFF(LOS+POS+4)

            IF (BUFF(LOS+POS+6) .GT. 0) THEN

               N = N + 1

               LSOR(N,jj,1) = BUFF(LOS+POS+6)
               LSOR(N,jj,2) = BUFF(LOS+POS+5)
               LSOR(N,jj,3) = BUFF(LOS+POS+7)
               LSOR(N,jj,4) = BUFF(LOS+POS+8)

            END IF

            POS = POS + PINCR

         END DO

         LOS = LOS + LINCR

         LSOR(N+1,jj,2) = LSOR(N,jj,3)

      END DO

      N = N - 1

      RETURN
      END

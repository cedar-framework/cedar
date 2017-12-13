      SUBROUTINE BMG2_SymStd_SETUP_ADD_SOR_PTRS(&
     &                             NProc, NLines,&
     &                             Min_GSz, NOL, &
     &                             pOFFSET,&
     &                             pSOR_LEV)&
     BIND(C, NAME='BMG2_SymStd_SETUP_ADD_SOR_PTRS')

         USE ModInterface
         IMPLICIT NONE

! ======================================================
! Intent In Variables
! ======================================================

         integer(c_int), value :: NProc, Min_Gsz
         integer(len_t), value :: NLines
         integer(c_int), value :: NOL
         integer(c_int) :: pOFFSET

! ======================================================
! Intent Out Variables
! ======================================================

         integer(c_int) :: pSOR_LEV(NOL)

! ======================================================
! Local Variables
! ======================================================

         INTEGER INSIZE, &
     &           NGPS,&
     &           REM,&
     &           CLUMP,&
     &           GSIZE,&
     &           RUN

         INTEGER kl

! ======================================================
! Initialize Counters
! ======================================================

         RUN = 0
         INSIZE = NProc

! ======================================================
! Step through levels and record pointers
! ======================================================


         DO kl = NOL-1, 2, -1

            IF ( kl .eq. NOL-1) THEN
               pSOR_LEV(kl) = 0
            ELSE
               pSOR_LEV(kl) = pSOR_LEV(kl+1) + &
     &                           (2*GSIZE+2)*NLines*4
            ENDIF
            !
            ! The number of groups to divide
            ! processes into
            !
            NGPS = INSIZE / Min_GSz

            !
            ! The remaining number of processes
            ! after diving line into groups of
            ! size Min_GSz
            !
            REM = mod(INSIZE, Min_GSz)

            !
            ! This just computes the ceiling of
            ! REM/NGPS.  CLUMP is the increment
            ! in which the remaining processes
            ! will be distributed
            !
            IF (real(REM)/NGPS .GT. REM/NGPS) THEN
               CLUMP = REM/NGPS + 1
            ELSE
               CLUMP = REM/NGPS
            END IF

            !
            ! Compute the size of the group
            !
            GSIZE = Min_GSz + CLUMP

            !
            ! Reset INSIZE
            !
            INSIZE = NGPS

            !
            ! Increment the running total
            !
            RUN = RUN + (2*GSIZE+2)*NLines*4

         END DO

         pOFFSET = RUN

      RETURN
      END

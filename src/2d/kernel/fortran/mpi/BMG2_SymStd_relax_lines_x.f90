      SUBROUTINE BMG2_SymStd_relax_lines_x( &
     &                K, SO, QF, Q, SOR, B,&
     &                II, JJ, iGs, jGs,&
     &                NOG, NStncl, IRELAX_SYM, UPDOWN,&
     &                DATADIST, iWork, NMSGi, pMSG, RWORK, NMSGr,&
     &                MPICOMM, XLINECOMM, YLINECOMM&
     &                ) BIND(C,NAME='MPI_BMG2_SymStd_relax_lines_x')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Perform zebra-line relaxation in x.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
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
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: II, JJ, NMSGi, NMSGr
      integer(c_int), value :: nog, nstncl
      integer(len_t), value :: iGs, jGs
      integer(c_int), value :: IRELAX_SYM, K, UPDOWN
      integer(len_t) :: iWork(NMSGi)
      integer(c_int) :: pMSG(NBMG_pMSG,NOG)
      integer, value :: MPICOMM
      integer(c_int), value :: XLINECOMM, YLINECOMM
      integer(len_t) :: DATADIST(2,*)
      real(real_t) :: B(II,JJ), Q(II,JJ), QF(II,JJ), SO(II+1,JJ+1,NStncl), &
           SOR(II,JJ,2), RWORK(NMSGr)

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, J, J1, JBEG, JEND
      INTEGER JBEG_START, JBEG_END, JBEG_STRIDE
      INTEGER ierror, myid, size

      INTEGER NP, IERR, OFFSET, NLINES
      INTEGER STATUS(MPI_STATUS_SIZE), MULT, TAG
      INTEGER ptrn
      INTEGER myidy, sizey
      INTEGER line_solve_comm_type

! ======================================================================

      line_solve_comm_type = BMG_LINE_SOLVE_COMM_TRADITIONAL
      TAG = 3

      J1=JJ-1
      I1=II-1
!
!     Local Declarations
!     on the way down we relax red lines and then black lines
!     on the way up we relax in the opposite order
!

      call MPI_Comm_Rank(XLINECOMM,myid,ierror)
      call MPI_Comm_Size(XLINECOMM,size,ierror)

      call MPI_Comm_Rank(YLINECOMM,myidy,ierror)
      call MPI_Comm_Size(YLINECOMM,sizey,ierror)
!
      IF ( UPDOWN.EQ.BMG_DOWN .OR. IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) THEN
         !
         !  Red, black => (start,end,stride) = (3,2,-1)
         !  Black, red => (start,end,stride) = (2,3,1)
         !
         JBEG_START  = 2 + MOD(jGs,2)
         JBEG_END    = 2 + MOD(jGs+1,2)
         JBEG_STRIDE = MOD(jGs+1,2) - MOD(jGs,2)
         !
      ELSEIF ( IRELAX_SYM.EQ.BMG_RELAX_SYM ) THEN
         !
         JBEG_START  = 2 + MOD(jGs+1,2)
         JBEG_END    = 2 + MOD(jGs,2)
         JBEG_STRIDE = MOD(jGs,2) - MOD(jGs+1,2)
         !
      ENDIF

      DO JBEG=JBEG_START, JBEG_END, JBEG_STRIDE

         IF ( NStncl.EQ.5 ) THEN
            !
            DO J=JBEG,J1,2
               DO I=2,I1
                  !
                  Q(I,J) = QF(I,J)&
     &                   + SO(I,J,KS)*Q(I,J-1)&
     &                   + SO(I,J+1,KS)*Q(I,J+1)&
     &                   + SO(I,J,KSW)*Q(I-1,J-1)&
     &                   + SO(I+1,J,KNW)*Q(I+1,J-1)&
     &                   + SO(I,J+1,KNW)*Q(I-1,J+1)&
     &                   + SO(I+1,J+1,KSW)*Q(I+1,J+1)
                  !
               ENDDO
            ENDDO
            !
         ELSE
            !
            DO J=JBEG,J1,2
               DO I=2,I1
                  Q(I,J) = QF(I,J)&
     &                   + SO(I,J,KS)*Q(I,J-1)&
     &                   + SO(I,J+1,KS)*Q(I,J+1)
               ENDDO
            ENDDO
            !
         ENDIF

         ! store the index of the last line
         JEND = J-2

         ! ====================================================

         NP = II-2
         CALL BMG2_SymStd_LineSolve_A(SOR,Q,II,JJ,&
     &        JBEG, RWORK, NP, NLINES)

         ! =====================================================

         IF (SIZE.gt.1) THEN    ! if there is more than one processor in
            if (myid .eq. 0) then
               CALL MPI_GATHER(MPI_IN_PLACE,NLINES*8,&
                    &           MPI_DOUBLE_PRECISION,&
                    &           RWORK,NLINES*8,MPI_DOUBLE_PRECISION,&
                    &           0, XLINECOMM, IERR)
            else
               CALL MPI_GATHER(RWORK,NLINES*8,&
                    &           MPI_DOUBLE_PRECISION,&
                    &           RWORK,NLINES*8,MPI_DOUBLE_PRECISION,&
                    &           0, XLINECOMM, IERR)
            endif
         ENDIF

         ! ======================================================

         MULT = 0
         DO J=JBEG,J1,2

            CALL BMG2_SymStd_LineSolve_B (SOR(2,J,1), SOR(3,J,2), &
     &           SOR(2,J,2), Q(2,J), &
     &           RWORK(SIZE*NLINES*8 + 1),&
     &           RWORK(MULT*8+1), &
     &           RWORK(MULT*8+1), &
     &           NP, size, DATADIST, myid, 8*NLINES, 8*NLINES)

            MULT = MULT+1

         ENDDO

         ! ====================================================

         IF (SIZE.gt.1) THEN
            if (myid .eq. 0) then
               CALL MPI_SCATTER(RWORK(MYID*NLINES*8+1),NLINES*8,&
                    &           MPI_DOUBLE_PRECISION,MPI_IN_PLACE, NLINES*8, &
                    &           MPI_DOUBLE_PRECISION,0,XLINECOMM,IERR)
            else
               CALL MPI_SCATTER(RWORK(MYID*NLINES*8+1),NLINES*8,&
                    &           MPI_DOUBLE_PRECISION,RWORK, NLINES*8, &
                    &           MPI_DOUBLE_PRECISION,0,XLINECOMM,IERR)
            endif
         ENDIF

         ! ====================================================

         MULT = 0
         DO J=JBEG,J1,2

            CALL BMG2_SymStd_LineSolve_C (SOR(2,J,1), SOR(3,J,2), &
     &           SOR(2,J,2), Q(1,J),  &
     &           RWORK(MULT*8 + 1),&
     &           RWORK(SIZE*NLINES*8 + 1),&
     &           NP, size, DATADIST,&
     &           myid)

            MULT = MULT+1

         ENDDO

         ! ====================================================

         IF (line_solve_comm_type .EQ.&
     &        BMG_LINE_SOLVE_COMM_TUNED)  THEN

            IF (SIZEY .GT. 1) THEN

               IF (MYIDY .EQ. 0) THEN
                  ! send to the top
                  IF (JEND .EQ. JJ-1) THEN
                     CALL MPI_Send(Q(1,JJ-1),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY+1,0,YLINECOMM,IERR)
                  ENDIF
               ENDIF

               IF (MYIDY .EQ. SIZEY-1) THEN
                  ! receive from the bottom
                  IF (JBEG .EQ. 3) THEN
                     CALL MPI_Recv(Q(1,1),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY-1,0,YLINECOMM,STATUS,IERR)
                  ENDIF
               ENDIF

               IF (MYIDY .GT. 0 .AND. MYIDY .LT. SIZEY-1) THEN
                  ! send to the top and receive from the bottom
                  IF (JBEG .EQ. 3 .AND. JEND .EQ. JJ-1) THEN
                     CALL MPI_Sendrecv(Q(1,JJ-1),II,&
     &                    MPI_DOUBLE_PRECISION,&
     &                    MYIDY+1,0,Q(1,1),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY-1,0,YLINECOMM,STATUS,IERR)

                  ELSE IF (JBEG .NE. 3 .AND. JEND .EQ. JJ-1) THEN
                     ! send to the top
                     CALL MPI_Send(Q(1,JJ-1),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY+1,0,YLINECOMM,IERR)

                  ELSE IF (JBEG .EQ. 3 .AND. JEND .NE. JJ-1) THEN
                     ! receive from the bottom
                     CALL MPI_Recv(Q(1,1),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY-1,0,YLINECOMM,STATUS,IERR)
                  ENDIF

               ENDIF


               IF (MYIDY .EQ. SIZEY-1) THEN
                  ! send to the bottom
                  IF (JBEG .EQ. 2) THEN
                     CALL MPI_Send(Q(1,2),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY-1,0,YLINECOMM,IERR)
                  ENDIF
               ENDIF

               IF (MYIDY .EQ. 0) THEN
                  ! receive from the top
                  IF (JEND .EQ. JJ-2) THEN
                     CALL MPI_Recv(Q(1,JJ),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY+1,0,YLINECOMM,STATUS,IERR)
                  ENDIF
               ENDIF

               IF (MYIDY .GT. 0 .AND. MYIDY .LT. SIZEY-1) THEN
                  ! send to the bottom and receive from the top
                  IF (JBEG .EQ. 2 .AND. JEND .EQ. JJ-2) THEN
                     CALL MPI_Sendrecv(Q(1,2),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY-1,0,Q(1,JJ),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY+1,0,YLINECOMM,STATUS,IERR)

                  ELSE IF (JBEG .EQ. 2 .AND. JEND .NE. JJ-2) THEN
                     ! send to the bottom
                     CALL MPI_Send(Q(1,2),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY-1,0,YLINECOMM,IERR)

                  ELSE IF (JBEG .NE. 2 .AND. JEND .EQ. JJ-2) THEN
                     ! receive from the top
                     CALL MPI_Recv(Q(1,JJ),II,MPI_DOUBLE_PRECISION,&
     &                    MYIDY+1,0,YLINECOMM,STATUS,IERR)
                  ENDIF

               ENDIF

            ENDIF

         ELSE  IF (line_solve_comm_type .EQ.&
     &           BMG_LINE_SOLVE_COMM_TRADITIONAL)  THEN

            ptrn = 1

            call MSG_tbdx_send(Q, rwork, &
     &           iWork(pMSG(ipL_MSG_NumAdjProc,K)),&
     &           iWork(pMSG(ipL_MSG_Proc,K)),&
     &           iWork(pMSG(ipL_MSG_Ipr,K)),&
     &           iWork(pMSG(ipL_MSG_Index,K)),&
     &           ptrn, ierror)

            call MSG_tbdx_receive(Q, rwork,&
     &           iWork(pMSG(ipL_MSG_NumAdjProc,K)),&
     &           iWork(pMSG(ipL_MSG_Proc,K)),&
     &           iWork(pMSG(ipL_MSG_Ipr,K)),&
     &           iWork(pMSG(ipL_MSG_Index,K)),&
     &           ptrn, ierror)

         ELSE

            WRITE(*,*) 'ERROR: invalid value for parameter'
            WRITE(*,*) '      BMG(id_BMG2_LINE_SOLVE_COMM_TYPE)=',&
     &           line_solve_comm_type
            STOP

         ENDIF


      ENDDO

      IF (line_solve_comm_type .EQ.&
     &        BMG_LINE_SOLVE_COMM_TRADITIONAL)  THEN
         !
         ! traditional communications scheme
         !
         ptrn = 1

         call MSG_tbdx_close(Q, rwork,&
     &        iWork(pMSG(ipL_MSG_NumAdjProc,K)),&
     &        iWork(pMSG(ipL_MSG_Proc,K)),&
     &        iWork(pMSG(ipL_MSG_Ipr,K)),&
     &        iWork(pMSG(ipL_MSG_Index,K)),&
     &        ptrn, ierror)

      ENDIF

! ======================================================================

      RETURN
      END

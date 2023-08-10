      SUBROUTINE BMG2_SymStd_relax_lines_y( &
     &                K, SO, QF, Q, SOR, B,&
     &                II, JJ, iGs, jGs,&
     &                NOG, NStncl, IRELAX_SYM, UPDOWN,&
     &                DATADIST, RWORK, NMSGr,&
     &                MPICOMM, XLINECOMM, YLINECOMM, halof&
     &                ) BIND(C,NAME='MPI_BMG2_SymStd_relax_lines_y')

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
      include 'mpif.h'
      include 'MSG_f90.h'

      include 'BMG_workspace_f90.h'


      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: II, JJ, NMSGr
      integer(c_int), value :: NOG, NStncl
      integer(len_t), value :: iGs, jGs
      integer(c_int), value :: IRELAX_SYM, K, UPDOWN
      integer, value :: MPICOMM
      integer(c_int), value :: XLINECOMM, YLINECOMM
      integer(c_int) :: DATADIST(2,*)
      real(real_t) ::   B(JJ,II), Q(II,JJ), QF(II,JJ), SO(II+1,JJ+1,NStncl),&
           SOR(JJ,II,2), RWORK(NMSGr)
      type(c_ptr) :: halof

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, IBEG, J, J1, IEND
      INTEGER IBEG_START, IBEG_END, IBEG_STRIDE
      INTEGER ierror, size, myid, sizex, myidx

      INTEGER NP, IERR, OFFSET, NLINES
      INTEGER STATUS(MPI_STATUS_SIZE), MULT, TAG
      INTEGER ptrn
      INTEGER line_solve_comm_type

! ======================================================================

      TAG = 4
      line_solve_comm_type = BMG_LINE_SOLVE_COMM_TRADITIONAL

      J1=JJ-1
      I1=II-1

      call MPI_Comm_rank(YLINECOMM,myid,ierror)
      call MPI_Comm_size(YLINECOMM,size,ierror)

      call MPI_Comm_rank(XLINECOMM,myidx,ierror)
      call MPI_Comm_size(XLINECOMM,sizex,ierror)

      IF ( UPDOWN.EQ.BMG_DOWN .OR. IRELAX_SYM.EQ.BMG_RELAX_NONSYM ) THEN
         !
         !  Red, black => (start,end,stride) = (3,2,-1)
         !  Black, red => (start,end,stride) = (2,3,1)
         !
         IBEG_START  = 2 + MOD(iGs,2)
         IBEG_END    = 2 + MOD(iGs+1,2)
         IBEG_STRIDE = MOD(iGs+1,2) - MOD(iGs,2)
         !
      ELSEIF ( IRELAX_SYM.EQ.BMG_RELAX_SYM ) THEN
         !
         IBEG_START  = 2 + MOD(iGs+1,2)
         IBEG_END    = 2 + MOD(iGs,2)
         IBEG_STRIDE = MOD(iGs,2) - MOD(iGs+1,2)
         !
      ENDIF

      DO IBEG=IBEG_START,IBEG_END,IBEG_STRIDE

         IF ( NStncl.EQ.5 ) THEN
            !
            DO  J=2,J1
               DO  I=IBEG,I1,2
                  !
                  B(J,I) = QF(I,J)&
     &                   + SO(I,J,KW)*Q(I-1,J)&
     &                   + SO(I+1,J,KW)*Q(I+1,J)&
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
            DO J=2,J1
               DO I=IBEG,I1,2
                  !
                  B(J,I) = QF(I,J)&
     &                   + SO(I,J,KW)*Q(I-1,J)&
     &                   + SO(I+1,J,KW)*Q(I+1,J)
                  !
               ENDDO
            ENDDO
            !
         ENDIF

         ! store the last line index
         IEND = I - 2

         NP = JJ-2
         CALL BMG2_SymStd_LineSolve_A(SOR,B,JJ,II,&
     &        IBEG, RWORK, NP, NLINES)

         ! =====================================================

         IF (SIZE .GT. 1) THEN
            if (myid .eq. 0) then
               CALL MPI_Gather(MPI_IN_PLACE,NLINES*8,&
                    &           MPI_DOUBLE_PRECISION,&
                    &           RWORK,NLINES*8,MPI_DOUBLE_PRECISION,&
                    &           0, YLINECOMM, IERR)
            else
               CALL MPI_Gather(RWORK,NLINES*8,&
                    &           MPI_DOUBLE_PRECISION,&
                    &           0.0d0,NLINES*8,MPI_DOUBLE_PRECISION,&
                    &           0, YLINECOMM, IERR)
            endif
         ENDIF

         ! ======================================================

         MULT = 0
         DO I=IBEG,I1,2

            NP = J1-1

            CALL BMG2_SymStd_LineSolve_B (SOR(2,I,1), SOR(3,I,2), &
     &           SOR(2,I,2), B(2,I),  &
     &           RWORK(SIZE*NLINES*8 + 1),&
     &           RWORK(MULT*8+1), &
     &           RWORK(MULT*8+1), &
     &           NP, size, DATADIST, myid, 8*NLINES, 8*NLINES)

            MULT = MULT+1

         ENDDO

         ! ====================================================

         IF (SIZE .GT. 1) THEN
            if (myid .eq. 0) then
               CALL MPI_Scatter(RWORK(MYID*NLINES*8+1),NLINES*8,&
                    &           MPI_DOUBLE_PRECISION,MPI_IN_PLACE,NLINES*8, &
                    &           MPI_DOUBLE_PRECISION,0,YLINECOMM,IERR)
            else
               CALL MPI_Scatter(0.0d0,NLINES*8,&
                    &           MPI_DOUBLE_PRECISION,RWORK, NLINES*8, &
                    &           MPI_DOUBLE_PRECISION,0,YLINECOMM,IERR)
            endif
         ENDIF

         ! ====================================================

         MULT = 0
         DO I=IBEG,I1,2

            NP = J1-1

            CALL BMG2_SymStd_LineSolve_C (SOR(2,I,1), SOR(3,I,2), &
     &           SOR(2,I,2), B(1,I), &
     &           RWORK(MULT*8 + 1),&
     &           RWORK(SIZE*NLINES*8 + 1),&
     &           NP, size, DATADIST,&
     &           myid)

            MULT = MULT+1

         ENDDO

         DO I=IBEG,I1,2
            DO J=1,JJ
               Q(I,J) = B(J,I)
            ENDDO
         ENDDO


         ! ====================================================

         IF (line_solve_comm_type .EQ.&
     &        BMG_LINE_SOLVE_COMM_TUNED)  THEN

            IF (SIZEX .GT. 1) THEN

               IF (MYIDX .EQ. 0) THEN
                  ! send to the right
                  IF (IEND .EQ. II-1) THEN
                     CALL MPI_Send(B(1,II-1),JJ,MPI_DOUBLE_PRECISION,&
     &                             MYIDX+1,0,XLINECOMM,IERR)
                  ENDIF
               ENDIF

               IF (MYIDX .EQ. SIZEX-1) THEN
                  ! receive from the left
                  IF (IBEG .EQ. 3) THEN
                     CALL MPI_Recv(B(1,1),JJ,MPI_DOUBLE_PRECISION,&
     &                             MYIDX-1,0,XLINECOMM,STATUS,IERR)
                     DO J=1,JJ
                        Q(1,J) = B(J,1)
                     ENDDO
                  ENDIF
               ENDIF


               IF (MYIDX .GT. 0 .AND. MYIDX .LT. SIZEX-1) THEN
                  ! send to the right and receive from the left
                  IF (IBEG .EQ. 3 .AND. IEND .EQ. II-1) THEN
                     CALL MPI_Sendrecv(B(1,II-1),JJ,&
     &                        MPI_DOUBLE_PRECISION,&
     &                        MYIDX+1,0,B(1,1),JJ,MPI_DOUBLE_PRECISION,&
     &                        MYIDX-1,0,XLINECOMM,STATUS,IERR)

                     DO J=1,JJ
                        Q(1,J) = B(J,1)
                     ENDDO

                  ELSE IF (IBEG .NE. 3 .AND. IEND .EQ. II-1) THEN
                     CALL MPI_Send(B(1,II-1),JJ,MPI_DOUBLE_PRECISION,&
     &                        MYIDX+1,0,XLINECOMM,IERR)

                  ELSE IF (IBEG .EQ. 3 .AND. IEND .NE. II-1) THEN
                     CALL MPI_Recv(B(1,1),JJ,MPI_DOUBLE_PRECISION,&
     &                        MYIDX-1,0,XLINECOMM,STATUS,IERR)

                     DO J=1,JJ
                        Q(1,J) = B(J,1)
                     ENDDO
                  ENDIF

               ENDIF

               ! then send to the left and receive from the right

               IF (MYIDX .EQ. SIZEX-1) THEN
                  ! send to the left
                  IF (IBEG .EQ. 2) THEN
                     CALL MPI_Send(B(1,2),JJ,MPI_DOUBLE_PRECISION,&
     &                    MYIDX-1,0,XLINECOMM,IERR)
                  ENDIF
               ENDIF

               IF (MYIDX .EQ. 0) THEN
                  ! receive from the right
                  IF (IEND .EQ. II-2) THEN
                     CALL MPI_Recv(B(1,II),JJ,MPI_DOUBLE_PRECISION,&
     &                    MYIDX+1,0,XLINECOMM,STATUS,IERR)

                     DO J=1,JJ
                        Q(II,J) = B(J,II)
                     ENDDO
                  ENDIF
               ENDIF

               IF (MYIDX .GT. 0 .AND. MYIDX .LT. SIZEX-1) THEN
                  IF (IBEG .EQ. 2 .AND. IEND .EQ. II-2) THEN
                      ! send to the left and receive from the right
                     CALL MPI_Sendrecv(B(1,2),JJ,MPI_DOUBLE_PRECISION,&
     &                    MYIDX-1,0,B(1,II),JJ,MPI_DOUBLE_PRECISION,&
     &                    MYIDX+1,0,XLINECOMM,STATUS,IERR)

                     DO J=1,JJ
                        Q(II,J) = B(J,II)
                     ENDDO

                  ELSE IF (IBEG .EQ. 2 .AND. IEND .NE. II-2) THEN
                     ! send to the left
                     CALL MPI_Send(B(1,2),JJ,MPI_DOUBLE_PRECISION,&
     &                    MYIDX-1,0,XLINECOMM,IERR)

                  ELSE IF (IBEG .NE. 2 .AND. IEND .EQ. II-2) THEN
                     ! receive from the right
                     IF (IEND .EQ. II-2) THEN
                        CALL MPI_Recv(B(1,II),JJ,MPI_DOUBLE_PRECISION,&
     &                       MYIDX+1,0,XLINECOMM,STATUS,IERR)

                        DO J=1,JJ
                           Q(II,J) = B(J,II)
                        ENDDO

                     ENDIF

                  ENDIF

               ENDIF

            ENDIF

         ELSE IF (line_solve_comm_type .EQ.&
     &           BMG_LINE_SOLVE_COMM_TRADITIONAL)  THEN

            call halo_exchange(K, Q, halof)
         ELSE

            WRITE(*,*) 'ERROR: invalid value for parameter'
            WRITE(*,*) '      BMG(id_BMG2_LINE_SOLVE_COMM_TYPE)=',&
     &           line_solve_comm_type
            STOP

         ENDIF


         ENDDO

! ======================================================================

      RETURN
      END

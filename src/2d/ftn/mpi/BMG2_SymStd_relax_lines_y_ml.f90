      SUBROUTINE BMG2_SymStd_relax_lines_y( &
           K, SO, QF, Q, SOR, B,&
           II, JJ, iGs, jGs,&
           NOG, NStncl, IRELAX_SYM, UPDOWN,&
           DATADIST, RWORK, NMSGr,&
           iface, niface, iface_ptrs, gwork, ngwork, gptr,&
           MPICOMM, &
           YCOMM, NOLY,&
           TDSY_SOR_PTRS,&
           NSOR, TDG, fact_flags, factorize, halof, mp)&
     BIND(C,NAME='MPI_BMG2_SymStd_relax_lines_y_ml')

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

      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value ::  II, JJ, NMSGr, NSOR
      integer(c_int), value :: ngwork, niface
      integer(c_int), value :: NOG, NStncl
      integer(len_t), value :: iGs, jGs
      integer(c_int), value :: K, IRELAX_SYM, UPDOWN
      integer, value :: MPICOMM
      integer(c_int) :: DATADIST(2,*)
      integer(c_int), value :: NOLY

      real(real_t) :: B(JJ,II), Q(II,JJ), QF(II,JJ), SO(II+1,JJ+1,NStncl)
      real(real_t) :: SOR(JJ,II,2), RWORK(NMSGr), TDG(NSOR)
      real(real_t) :: iface(niface), gwork(ngwork)

      integer :: YCOMM(2,2:NOLY)
      integer(c_int) :: gptr(2, 2:NOLY), iface_ptrs(2)
      integer(c_int) :: TDSY_SOR_PTRS(NOLY)

      logical(c_bool) :: fact_flags(2 * NOG)
      logical(c_bool), value :: factorize
      type(c_ptr) :: halof, mp

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, IBEG, J, J1, IEND
      INTEGER IBEG_START, IBEG_END, IBEG_STRIDE
      INTEGER ierror, size, myid, sizex, myidx

      INTEGER NPts, IERR, OFFSET, NLINES
      INTEGER STATUS(MPI_STATUS_SIZE), MULT, TAG
      INTEGER ptrn

      INTEGER CP
      INTEGER AT, BT, CT, FT, XT
      integer :: icolor, iface_len, iface_ptr

! ======================================================================

      TAG = 4

      J1=JJ-1
      I1=II-1


      call MPI_Comm_rank(YCOMM(1,NOLY),myid,ierror)
      call MPI_Comm_size(YCOMM(1,NOLY),size,ierror)

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

         NPts = JJ-2

         icolor = mod(ibeg, 2) + 1
         iface_ptr = iface_ptrs(icolor) + 1
         iface_len = ((II-2) / 2 + mod(ibeg,2) * mod(II-2,2)) * 8

         CALL BMG2_SymStd_downline(SOR,B,JJ,II,&
              IBEG, RWORK, gwork, ngwork, gptr,&
              iface(iface_ptr), iface_len, NPts, NLINES,&
              K, NOLY, YCOMM, NSOR, TDG, &
              NMSGr, NOG, TDSY_SOR_PTRS,&
              CP, fact_flags, factorize, mp)

         !
         ! Pointers into RWORK
         !
         AT = 1
         BT = AT + 2*CP + 2
         CT = BT + 2*CP + 2
         FT = CT + 2*CP + 2
         XT = FT + 2*CP + 2

         MULT = 0
         DO I=IBEG,I1,2

            NPts = J1-1

            CALL BMG2_SymStd_LineSolve_B_ml(gwork(gptr(icolor,2) + MULT*8 + 1),&
                 gwork(gptr(icolor,2) + MULT*8+1),&
                 RWORK(AT),RWORK(BT),RWORK(CT),RWORK(FT),RWORK(XT),&
                 2*CP, CP, DATADIST, myid, 8*NLINES, 8*NLINES)

            MULT = MULT+1

         ENDDO

         call BMG2_SymStd_upline(SOR, B, JJ, II,&
              IBEG, RWORK, gwork, ngwork, gptr,&
              iface(iface_ptr), iface_len, NPts, NLINES, &
              K, NOLY, YCOMM, NSOR, TDG, &
              NMSGr, NOG, TDSY_SOR_PTRS, CP, factorize, mp)

         DO I=IBEG,I1,2
            DO J=1,JJ
               Q(I,J) = B(J,I)
            ENDDO
         ENDDO

         call halo_exchange(K, Q, halof)


         ENDDO

! ======================================================================

      RETURN
      END

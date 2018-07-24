      SUBROUTINE BMG2_SymStd_relax_lines_x_ml( &
           K, SO, QF, Q, SOR, B,&
           II, JJ, iGs, jGs,&
           NOG, NStncl, IRELAX_SYM, UPDOWN,&
           DATADIST, RWORK, NMSGr,&
           iface, niface, iface_ptrs, gwork, ngwork, gptr,&
           MPICOMM, &
           XCOMM, NOLX,&
           TDSX_SOR_PTRS, &
           NSOR, TDG, fact_flags, factorize, halof, mp)&
     BIND(C, NAME='MPI_BMG2_SymStd_relax_lines_x_ml')

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

      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: II, JJ, NMSGr, NSOR
      integer(c_int), value :: NOG, NStncl, IRELAX_SYM, K, UPDOWN

      integer(len_t), value :: iGs, jGs

      integer, value :: MPICOMM
      integer :: XCOMM(2, 2:NOLX)
      integer(c_int) :: gptr(2, 2:NOLX), iface_ptrs(2)
      integer(c_int) :: TDSX_SOR_PTRS(NOLX)
      integer(len_t) :: DATADIST(2,*)
      integer(c_int), value :: NOLX
      integer(c_int), value :: ngwork, niface

      real(real_t) :: B(II,JJ), Q(II,JJ), QF(II,JJ), SO(II+1,JJ+1,NStncl)
      real(real_t) :: SOR(II,JJ,2), RWORK(NMSGr), TDG(NSOR)
      real(real_t) :: iface(niface), gwork(ngwork)

      logical(c_bool) :: fact_flags(2 * NOG)
      logical(c_bool), value :: factorize
      type(c_ptr) :: halof, mp

! ----------------------------
!     Local Declarations
!
      INTEGER I, I1, J, J1, JBEG, JEND
      INTEGER JBEG_START, JBEG_END, JBEG_STRIDE
      INTEGER ierror, myid, size

      INTEGER Npts, IERR, OFFSET, NLINES
      INTEGER STATUS(MPI_STATUS_SIZE), MULT, TAG
      INTEGER ptrn

      INTEGER CP
      INTEGER AT, BT, CT, FT, XT
      integer :: icolor, iface_len, iface_ptr

! ======================================================================

      TAG = 3

      J1=JJ-1
      I1=II-1
!
!     Local Declarations
!     on the way down we relax red lines and then black lines
!     on the way up we relax in the opposite order
!

      call MPI_Comm_Rank(XCOMM(1,NOLX),myid,ierror)
      call MPI_Comm_Size(XCOMM(1,NOLX),size,ierror)

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

! ======================================================================

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
                  !
               ENDDO
            ENDDO
            !
         ENDIF

! ======================================================================

         ! store the index of the last line
         JEND = J-2

         ! ====================================================

         Npts = II-2

         icolor = mod(jbeg, 2) + 1
         iface_ptr = iface_ptrs(icolor) + 1
         iface_len = ((JJ-2) / 2 + mod(jbeg,2) * mod(JJ-2,2)) * 8

         ! ====================================================
         ! Multilevel LineSolve
         ! ====================================================

         !
         ! Coarsen until we have system of interface equations
         ! on process 0
         !

         CALL BMG2_SymStd_downline(SOR,Q,II,JJ, &
              JBEG, RWORK, gwork, ngwork, gptr,&
              iface(iface_ptr), iface_len, Npts, NLines, &
              K, NOLX, XCOMM, NSOR, TDG, &
              NMSGr, NOG, TDSX_SOR_PTRS,&
              CP, fact_flags, factorize, mp)


         ! Pointers into RWORK
         !

         AT = 1
         BT = AT + 2*CP + 2
         CT = BT + 2*CP + 2
         FT = CT + 2*CP + 2
         XT = FT + 2*CP + 2


         !
         ! Solve systems of interface equations on process 0
         !

         MULT = 0
         DO J=JBEG,J1,2

            CALL BMG2_SymStd_LineSolve_B_ml(gwork(gptr(icolor,2) + MULT*8+1),&
                 gwork(gptr(icolor,2) + MULT*8+1), &
                 rwork(AT),rwork(BT),rwork(CT),rwork(FT),rwork(XT),&
                 2*CP, CP, DATADIST, myid, 8*NLINES, 8*NLINES)

            MULT = MULT + 1

         END DO


!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
!                     More Multivelel Stuff
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

         CALL BMG2_SymStd_upline(SOR,Q,II,JJ, &
              JBEG, RWORK, gwork, ngwork, gptr,&
              iface(iface_ptr), iface_len, Npts, NLines,&
              K, NOLX, XCOMM, NSOR, TDG, &
              NMSGr, NOG, TDSX_SOR_PTRS,&
              CP, factorize, mp)



         ! ====================================================

         call halo_exchange(K, Q, halof)

      ENDDO

! ======================================================================

      RETURN
      END

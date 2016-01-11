      SUBROUTINE BMG3_SymStd_SOLVE_cg( &
           q, qf, ii, jj, kk, abd, bbd, nabd1, nabd2,&
           IBC) BIND(C, NAME='BMG3_SymStd_SOLVE_cg')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SOLVE_cg does a direct solve on the coarsest grid. it
!     uses the LAPACK routine DPBTRS.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     QF        Refer to BMG3_SymStd_SOLVE_boxmg
!
!     II        Number of grid points in x direction, including
!               two fictitious points.
!     JJ        Number of grid points in y direction, including
!               two fictitious points.
!     KK        Number of grid points in z direction, including
!               two fictitious points.
!
!     ABD       Refer to BMG3_SymStd_SOLVE_boxmg
!
!     NABD1     Refer to BMG3_SymStd_SOLVE_boxmg
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!     BBD       Workspace, refer to BMG3_SymStd_SOLVE_boxmg
!
!  --------------------
!   OUTPUT:
!  --------------------
!
!     Q         Refer to BMG3_SymStd_SOLVE_boxmg
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
      INCLUDE 'BMG_parameters_f90.h'
      INCLUDE 'BMG_constants_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(c_int), value :: IBC
      integer(len_t), value :: ii, jj, kk, nabd1, nabd2
      real(real_t) :: abd(nabd1,nabd2), bbd(nabd2), q(ii,jj,kk), qf(ii,jj,kk)

! ----------------------------
!     Local Declarations
!
      integer i, i1, i2, ibw, j, j1, k, kt, k1, n, info
      integer IPN
      INTEGER PER_x, PER_y, PER_xy, PER_z, PER_xz, PER_yz, PER_xyz
      INTEGER PER_sum_x, Per_sum_y, Per_sum_z
      REAL*8  C, CINT, QINT
! ======================================================================
      IPN = IABS(IBC)
      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)
      PER_z = IABS(BMG_BCs_def_per_z)
      PER_xz = IABS(BMG_BCs_def_per_xz)
      PER_yz = IABS(BMG_BCs_def_per_yz)
      PER_xyz = IABS(BMG_BCs_def_per_xyz)
      Per_sum_x = 0
      IF(IPN.EQ.PER_x .OR. IPN.EQ.PER_xy .OR. IPN .EQ.PER_xz&
     &     .OR. IPN.EQ.PER_xyz) THEN
         PER_sum_x = 1
      ENDIF
      PER_sum_y = 0
      IF(IPN.EQ.PER_y .OR. IPN.EQ.PER_xy .OR. IPN.EQ.PER_yz&
     &     .OR. IPN.EQ.PER_xyz) THEN
         PER_sum_y = 1
      ENDIF
      PER_sum_z = 0
      IF(IPN.EQ.PER_z .OR. IPN.EQ.PER_xz .OR. IPN.EQ.PER_yz&
     &     .OR. IPN.EQ.PER_xyz) THEN
         PER_sum_z = 1
      ENDIF
!
!     direct solve on coarsest grid
!

      i1=ii-1
      j1=jj-1
      k1=kk-1

      i2=i1-1
      n=i2*(j1-1)*(kk-2)
      ibw=i2*j1+1

      kt=0
      do k=2,k1
         do j=2,j1
            do i=2,i1
               kt=kt+1
               bbd(kt)=qf(i,j,k)
            enddo
         enddo
      enddo

! -------------------------------------------------------
!     Solve using the LAPACK routine DPBTRS OR DPOTRS
! -------------------------------------------------------

      IF(Per_sum_x.eq.1 .or. Per_sum_y.eq.1 .or. Per_sum_z.eq.1) THEN
         CALL DPOTRS('U',KT,1,ABD,NABD1,BBD,NABD2,INFO)
      ELSE
         CALL DPBTRS ('U', KT, IBW, 1, ABD, NABD1, BBD, NABD2, INFO)
      ENDIF

      IF (INFO .NE. 0) THEN
         WRITE(*,510) 'INFO = ', INFO
         call print_error(C_CHAR_"Coarse grid solve failed!"//C_NULL_CHAR)
         RETURN
      ENDIF


      kt=0
      do k=2,k1
         do j=2,j1
            do i=2,i1
               kt=kt+1
               q(i,j,k)=bbd(kt)
            enddo
         enddo
      enddo

      IF ( IBC.NE.0 ) THEN

      CINT=rZERO
      QINT=rZERO
      DO k=2,k1
         DO j=2,j1
            DO i=2,i1
               QINT=QINT+q(i,j,k)
               CINT=CINT+1
            ENDDO
         ENDDO
      ENDDO
      C=-QINT/CINT
      DO k = 2,k1
         DO j=2,j1
            DO i=2,i1
               Q(i,j,k)=Q(i,j,k)+C
            ENDDO
         ENDDO
      ENDDO


      ENDIF

      IF(PER_sum_x .EQ. 1) THEN
         DO k = 1,kk
            DO j = 1,jj
               Q(1,j,k) = Q(i1,j,k)
               Q(ii,j,k) = Q(2,j,k)
            ENDDO
         ENDDO
      ENDIF

      IF(PER_sum_y .EQ. 1) THEN
         DO k = 1,kk
            DO i = 1,ii
               Q(i,1,k) = Q(i,j1,k)
               Q(i,jj,k) = Q(i,2,k)
            ENDDO
         ENDDO
      ENDIF

      IF(PER_sum_z .EQ. 1) THEN
         DO j = 1,jj
            DO i = 1,ii
               Q(i,j,1) = Q(i,j,k1)
               Q(i,j,kk) = Q(i,j,2)
            ENDDO
         ENDDO
      ENDIF

! ======================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SOLVE_cg.f',/,5X,A)
 510  FORMAT (5X,A,1X,I3)

! ===========================================

      RETURN
      END

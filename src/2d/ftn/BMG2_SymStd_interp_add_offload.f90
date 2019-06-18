      SUBROUTINE BMG2_SymStd_interp_add(&
     &                       Q ,QC, RES, SO, CI,&
     &                       IIC, JJC, IIF, JJF, NStncl, JPN &
     &                       ) BIND(C, NAME='BMG2_SymStd_interp_add_offload')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_PerSymStd_interp_add interpolates Q from the coarse mesh KC
!     to the fine mesh KF and adds result to Q on fine mesh.
!
! ======================================================================
! $license_flag$
! ======================================================================
!  --------------------
!   INPUT:
!  --------------------
!
!     IIC       Number of grid points in x direction on coarse grid,
!               including two fictitious points.
!     JJC       Number of grid points in y direction on coarse grid,
!               including two fictitious points.
!
!     IIF       Number of grid points in x direction on fine grid,
!               including two fictitious points.
!     JJF       Number of grid points in y direction on fine grid,
!               including two fictitious points.
!
!     QC        Q for coarse grid.
!
!     SOR       Refer to BMG2_SymStd_SOLVE_boxmg.
!     CI        Refer to BMG2_SymStd_SOLVE_boxmg.
!
!     IPN       Refer to BMG2_SymStd_SOLVE_boxmg.
!
! ======================================================================
!  --------------------
!   INPUT/OUTPUT:
!  --------------------
!
!     Q         Refer to BMG2_SymStd_SOLVE_boxmg.
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
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER(len_t), VALUE :: IIC, IIF, JJC, JJF
      INTEGER(C_INT), VALUE :: JPN, NStncl
      REAL(real_t) :: CI(IIC,JJC,8), Q(IIF,JJF), QC(IIC,JJC),&
     &                SO(IIF,JJF,NStncl), RES(IIF,JJF)

! ----------------------------
!     Local Declarations
!
      INTEGER  IC, I, IICF, IICF1, IIC1, IIF1, JC, J, JJCF, JJCF1, &
     &         JJC1, JJF1, IPN, PER_x, PER_y, PER_xy
      REAL*8   A, AQ

! ======================================================================

!
!   interpolate answers from coarse to fine mesh and add
!   to answers on fine mesh.
!
      JJF1=JJF-1
      IIF1=IIF-1
      IIC1=IIC-1
      JJC1=JJC-1
      IICF=(IIF-2)/2+3
      JJCF=(JJF-2)/2+3
      IICF1=IICF-1
      JJCF1=JJCF-1
! --------------------------------------------------
!     NB: division is only possible in the interior
! --------------------------------------------------

      !$omp target teams distribute parallel do simd collapse(2)
      DO J = 2, JJF1
         DO I = 2, IIF1
            RES(I,J)=RES(I,J)/SO(I,J,KO)
         ENDDO
      ENDDO
      !$omp end target teams distribute parallel do

! --------------------------------------------------
!   interpolate answers from coarse to fine mesh
!   and add to answers on fine mesh.
! --------------------------------------------------
      J=2
      I=2
      Q(2,J)=Q(2,J)+QC(2,2)
      !$omp target teams
      !$omp distribute parallel do simd
      DO IC=3,IICF1
         I=(IC-1)*2
         J=2
         Q(I,J)=Q(I,J)+QC(IC,2)
         A=CI(IC,2,LR)*QC(IC,2)+CI(IC,2,LL)*QC(IC-1,2)
         Q(I-1,J)=Q(I-1,J)+A+RES(I-1,J)
      ENDDO
      !$omp end distribute parallel do
      !$omp distribute parallel do simd
      DO JC=3,JJCF1
         J=(JC-1)*2
         I=2
         Q(2,J)=Q(2,J)+QC(2,JC)
         AQ=CI(2,JC,LA)*QC(2,JC)+CI(2,JC,LB)*QC(2,JC-1)
         Q(2,J-1)=Q(2,J-1)+AQ+RES(2,J-1)
      enddo
      !$omp end distribute parallel do
      !$omp end target teams

      !$omp target teams distribute parallel do simd collapse(2)
      DO JC=3,JJCF1
         DO  IC=3,IICF1
            I=(IC-1)*2
            J=(JC-1)*2
            Q(I,J)=Q(I,J)+QC(IC,JC)
            A=CI(IC,JC,LR)*QC(IC,JC)+CI(IC,JC,LL)*QC(IC-1,JC)
            Q(I-1,J)=Q(I-1,J)+A+RES(I-1,J)
            AQ=CI(IC,JC,LA)*QC(IC,JC)+CI(IC,JC,LB)*QC(IC,JC-1)
            Q(I,J-1)=Q(I,J-1)+AQ+RES(I,J-1)
            A=CI(IC,JC,LSW)*QC(IC-1,JC-1)+CI(IC,JC,LNW)*QC(IC-1,JC)&
     &           +CI(IC,JC,LNE)*QC(IC,JC)+CI(IC,JC,LSE)*QC(IC,JC-1)
            Q(I-1,J-1)=Q(I-1,J-1)+A+RES(I-1,J-1)
         ENDDO
      ENDDO
      !$omp end target teams distribute parallel do
      IPN=IABS(JPN)
      IF(IPN.EQ.BMG_BCs_definite .OR. IPN.EQ.BMG_BCs_indef_nonper)&
     & RETURN
      PER_x = IABS(BMG_BCs_def_per_x)
      PER_y = IABS(BMG_BCs_def_per_y)
      PER_xy = IABS(BMG_BCs_def_per_xy)
      IF(IPN.EQ.PER_y.OR.IPN.EQ.PER_xy) THEN
         !$omp target teams distribute parallel do simd
         DO I=1,IIF
            Q(I,1)=Q(I,JJF1)
            Q(I,JJF)=Q(I,2)
         ENDDO
         !$omp end target teams distribute parallel do
      ENDIF
      IF(IPN.EQ.PER_x.OR.IPN.EQ.PER_xy) THEN
         !$omp target teams distribute parallel do simd
         DO J=1,JJF
            Q(1,J)=Q(IIF1,J)
            Q(IIF,J)=Q(2,J)
         ENDDO
         !$omp end target teams distribute parallel do
      ENDIF
! ======================================================================

      RETURN
      END

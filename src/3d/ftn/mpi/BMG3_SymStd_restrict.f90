      SUBROUTINE BMG3_SymStd_restrict(&
     &                KFG, KCG, &
     &                Q, QC, CI, Nx, Ny, Nz, Nxc, Nyc, Nzc,&
     &                iGs, jGs, kGs&
     &                ) BIND(C, NAME='MPI_BMG3_SymStd_restrict')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_restrict computes the restriction of a vector on the
!     fine grid, Q, to a vector on the coarse grid, QC.  The weights
!     involve the transpose of the interpolation operator from the
!     coarse grid to the fine grid.
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
      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: Nx, Nxc, Ny, Nyc, Nz, Nzc,&
           iGs, jGs, kGs
      integer(c_int), value :: KCG, KFg
      real(real_t) :: CI(Nxc,Nyc,Nzc,26), Q(Nx,Ny,Nz), QC(Nxc,Nyc,Nzc)

! ----------------------------
!     Local Declarations
!
      INTEGER  i, ic, j, jc, k, kc
      INTEGER  ISTART, JSTART, KSTART

! ======================================================================

! ---------------------------------------------
!     Restrict the vector Q -> QC:
! ---------------------------------------------

      IF ( MOD(iGs,2).EQ.1 ) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF

      IF ( MOD(jGs,2).EQ.1 ) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF

      IF ( MOD(kGs,2).EQ.1 ) THEN
         KSTART = 0
      ELSE
         KSTART = 1
      ENDIF



      k=KSTART
      DO kc=2, Nzc-1
         k=k+2
         j=JSTART
         DO jc=2, Nyc-1
            j=j+2
            i=ISTART
            DO ic=2, Nxc-1
               i=i+2
               QC(ic,jc,kc) = CI(ic,jc,kc,lxyne)*Q(i-1,j-1,k)&
     &                      + CI(ic,jc,kc,lxya)*Q(i,j-1,k)&
     &                      + CI(ic+1,jc,kc,lxynw)*Q(i+1,j-1,k)&
     &                      + CI(ic,jc,kc,lxyr)*Q(i-1,j,k)&
     &                      + Q(i,j,k)&
     &                      + CI(ic+1,jc,kc,lxyl)*Q(i+1,j,k)&
     &                      + CI(ic,jc+1,kc,lxyse)*Q(i-1,j+1,k)&
     &                      + CI(ic,jc+1,kc,lxyb)*Q(i,j+1,k)&
     &                      + CI(ic+1,jc+1,kc,lxysw)*Q(i+1,j+1,k)&
     &                      + CI(ic,jc,kc,ltne)*Q(i-1,j-1,k-1)&
     &                      + CI(ic,jc,kc,lyznw)*Q(i,j-1,k-1)&
     &                      + CI(ic+1,jc,kc,ltnw)*Q(i+1,j-1,k-1)&
     &                      + CI(ic,jc,kc,lxzne)*Q(i-1,j,k-1)&
     &                      + CI(ic,jc,kc,lxza)*Q(i,j,k-1)&
     &                      + CI(ic+1,jc,kc,lxznw)*Q(i+1,j,k-1)&
     &                      + CI(ic,jc+1,kc,ltse)*Q(i-1,j+1,k-1)&
     &                      + CI(ic,jc+1,kc,lyzne)*Q(i,j+1,k-1)&
     &                      + CI(ic+1,jc+1,kc,ltsw)*Q(i+1,j+1,k-1)&
     &                      + CI(ic,jc,kc+1,lbne)*Q(i-1,j-1,k+1)&
     &                      + CI(ic,jc,kc+1,lyzsw)*Q(i,j-1,k+1)&
     &                      + CI(ic+1,jc,kc+1,lbnw)*Q(i+1,j-1,k+1)&
     &                      + CI(ic,jc,kc+1,lxzse)*Q(i-1,j,k+1)&
     &                      + CI(ic,jc,kc+1,lxzb)*Q(i,j,k+1)&
     &                      + CI(ic+1,jc,kc+1,lxzsw)*Q(i+1,j,k+1)&
     &                      + CI(ic,jc+1,kc+1,lbse)*Q(i-1,j+1,k+1)&
     &                      + CI(ic,jc+1,kc+1,lyzse)*Q(i,j+1,k+1)&
     &                      + CI(ic+1,jc+1,kc+1,lbsw)*Q(i+1,j+1,k+1)
               enddo

            enddo

         enddo

! ======================================================================

      return
      end

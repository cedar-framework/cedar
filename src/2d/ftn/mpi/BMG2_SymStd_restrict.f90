      SUBROUTINE BMG2_SymStd_restrict( &
     &                KF, KC, NOG,&
     &                Q, QC, CI,&
     &                Nx, Ny, Nxc, Nyc, iGs, jGs&
     &                ) BIND(C, NAME='MPI_BMG2_SymStd_restrict')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG2_SymStd_restrict computes the restriction of a vector on the
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
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ----------------------------
!     Argument Declarations
!
      INTEGER(len_t), VALUE :: Nx, Nxc,&
           Ny, Nyc, iGs, jGs
      INTEGER(C_INT), VALUE :: NOG, KC, KF
      REAL(real_t) :: CI(Nxc,Nyc,8), Q(Nx,Ny), QC(Nxc,Nyc)

! ----------------------------
!     Local Declarations
!
      INTEGER  i, ic, j, jc, ISTART, JSTART

! ======================================================================

! ---------------------------------------------
!     Restrict the vector Q -> QC:
! ---------------------------------------------


      IF (mod(iGs,2).eq.1) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF

      IF (mod(jGs,2).eq.1) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF


      j=JSTART
      DO jc=2, Nyc-1
         j=j+2
         i=ISTART
         DO ic=2, Nxc-1
            i=i+2
            QC(ic,jc) = Ci(ic,jc,LNE)*Q(i-1,j-1)&
     &                + Ci(ic,jc,LA)*Q(i,j-1)&
     &                + Ci(ic+1,jc,LNW)*Q(i+1,j-1)&
     &                + Ci(ic,jc,LR)*Q(i-1,j)&
     &                + Q(i,j)&
     &                + Ci(ic+1,jc,LL)*Q(i+1,j)&
     &                + Ci(ic,jc+1,LSE)*Q(i-1,j+1)&
     &                + Ci(ic,jc+1,LB)*Q(i,j+1)&
     &                + Ci(ic+1,jc+1,LSW)*Q(i+1,j+1)
          ENDDO
       ENDDO

       ! Note: no update of the ghost bdry of QC is necessary

!       ptrn = 1

!       call MSG_tbdx_send(QC, MSG_Buffer,
!     &      iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Proc,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Ipr,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Index,KC)),
!     &      ptrn, ierror)

!       call MSG_tbdx_receive(QC, MSG_Buffer,
!     &      iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Proc,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Ipr,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Index,KC)),
!     &      ptrn, ierror)

!       call MSG_tbdx_close(QC, MSG_Buffer,
!     &      iWorkMSG(pMSG(ipL_MSG_NumAdjProc,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Proc,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Ipr,KC)),
!     &      iWorkMSG(pMSG(ipL_MSG_Index,KC)),
!     &      ptrn, ierror)


! ======================================================================

       RETURN
       END

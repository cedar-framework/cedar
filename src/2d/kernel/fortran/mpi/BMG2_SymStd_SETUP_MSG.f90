      SUBROUTINE BMG2_SymStd_SETUP_MSG(&
     &                pMSG, pMSGSO, iMSG_Geom, NMSGi, pSI_MSG,&
     &                IBC, IGRD, NOG, NOGm, NProc, MyProc,  &
     &                DimX, DimY, DimXfine, DimYfine,&
     &                ProcGrid, NProcI, NProcJ, MPICOMM&
     &                ) BIND(C, NAME='BMG2_SymStd_SETUP_MSG')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Create pointers into the integer work space describing the storage
!     of MSG communication setup.  Note that the pointer shift, pSI_MSG,
!     is changed in every call.
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

      INCLUDE 'BMG_parameters_f90.h'
      INCLUDE 'BMG_workspace_f90.h'

! ---------------------------
!    Argument Declarations:
!
      INTEGER(C_INT),  VALUE :: IBC, NOG, NOGm, NProcI, NProcj
      INTEGER(len_t), VALUE :: NMSGi
      INTEGER, VALUE :: MPICOMM
      INTEGER(len_t) :: IGRD(NOGm,NBMG_pIGRD), iMSG_Geom(NMSGi)
      INTEGER(C_INT) :: pMSG(NBMG_pMSG,NOG), pMSGSO(NBMG_pMSG,NOG)
      INTEGER(C_INT) :: DimX(NProcI,NOGm), DimY(NProcJ,NOGm)
      INTEGER(C_INT) :: DimXfine(NProcI), DimYfine(NProcJ)
      INTEGER(C_INT) ::  ProcGrid(NProcI,NProcJ)
      INTEGER(C_INT) :: pSI_MSG
      INTEGER(C_INT), VALUE :: NProc, MyProc

! --------------------------
!     Local Declarations:
!
      INTEGER I, J, iGs, jGs
      Integer ierror, N

! ======================================================================

      call MSG_set_comm_parent(MPICOMM)
      call MSG_enable(MyProc, NProc)
!     Note that the pointer shift, pSI_MSG, is changed in every call

      CALL BMG2_SymStd_SETUP_PtrMSG(&
     &     IGRD(NOG,idL_BMG_NLx), IGRD(NOG,idL_BMG_NLy), pSI_MSG,&
     &     NProc, NOG, pMSG)

      CALL BMG2_SymStd_SETUP_PtrMSGSO(&
     &     IGRD(NOG,idL_BMG_NLx), IGRD(NOG,idL_BMG_NLy), pSI_MSG,&
     &     NProc, NOG, pMSGSO)


! -------------------------------------------
!     Compute the dimensions for all local
!     grids on all processors
! -------------------------------------------

      DO I=1, NProcI
         DimX(I,NOG) = DimXfine(I)
      END DO

      DO J=1, NProcJ
         DimY(J,NOG) = DimYfine(J)
      END DO

      DO N=NOG-1, 1, -1

         iGs = 1

         DO I=1, NProcI

            if (mod(iGs,2).eq.1) then
               DimX(I,N) = (DimX(I,N+1)+1)/2
            else
               IF (mod(DimX(I,N+1),2).eq.1) THEN
                  DimX(I,N) = (DimX(I,N+1)-1)/2
               ELSE
                  DimX(I,N) = (DimX(I,N+1)+1)/2
               ENDIF
            endif

            iGs = iGs + DimX(I,N+1)

         END DO

         jGs = 1

         DO J=1, NProcJ

            if (mod(jGs,2).eq.1) then
               DimY(J,N) = (DimY(J,N+1)+1)/2
            else
               IF (mod(DimY(J,N+1),2).eq.1) THEN
                  DimY(J,N) = (DimY(J,N+1)-1)/2
               ELSE
                  DimY(J,N) = (DimY(J,N+1)+1)/2
               ENDIF
            endif

            jGs = jGs + DimY(J,N+1)

         END DO

      END DO


! -------------------------------------------
!     Create MSG grid information for U
! -------------------------------------------

      DO N=NOG, 1, -1
         !
         CALL BMG2_SymStd_SETUP_MSGGrid(&
     &        IGRD(N,idL_BMG_NGx), IGRD(N,idL_BMG_NGy), IBC,&
     &        iMSG_Geom(pMSG(ipL_MSG_LocalArraySize,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordLocalData,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordActData,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_ActDataStart,N)),&
     &        DimX, DimY, ProcGrid, MyProc, &
     &        NProc, NProcI, NProcJ, NOGm, N,&
     &        MPICOMM)
         !
      ENDDO

! --------------------------------------------
!     Create MSG grid information for SO
! --------------------------------------------

      DO N=NOG, 1, -1
         !
         CALL BMG2_SymStd_SETUP_MSGGridSO(&
     &        IGRD(N,idL_BMG_NGx), IGRD(N,idL_BMG_NGy),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_LocalArraySize,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordLocalData,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordActData,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_ActDataStart,N)),&
     &        DimX, DimY, ProcGrid, NProc, NProcI, NProcJ, NOGm, N,&
     &        MPICOMM )
         !
      ENDDO


!
!     Setup MSG communications for all grids
!

      DO N=NOG, 1, -1

!         write(*,*) N, NProc, MyProc

         iMSG_Geom(pMSG(ipL_MSG_NumAdjProc, N)) = 0

!         write(*,*)MyProc,':', (iWork(pMSG(ipL_MSG_ActDataStart,N)+i),
!     &        i=0,3*NProc-1)

         CALL MSG_tp_setup (&
     &        iMSG_Geom(pMSG(ipL_MSG_LocalArraySize,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_ActDataStart,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordLocalData,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordActData,N)),&
     &        NProc, MyProc,&
     &        iMSG_Geom(pMSG(ipL_MSG_NumAdjProc,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_Proc,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_Ipr,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_Index,N)),&
     &        1,1,ierror)

!         write(*,*) 'after MSG_tp_setup'


!         write(*,*) N,IGRD(N,idL_BMG_ICOORD),IGRD(N,idL_BMG_JCOORD)


         iMSG_Geom(pMSGSO(ipL_MSG_NumAdjProc, N)) = 0

         CALL MSG_tp_setup (&
     &        iMSG_Geom(pMSGSO(ipL_MSG_LocalArraySize,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_ActDataStart,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordLocalData,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordActData,N)),&
     &        NProc, MyProc,&
     &        iMSG_Geom(pMSGSO(ipL_MSG_NumAdjProc,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_Proc,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_Ipr,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_Index,N)),&
     &        1,1,ierror)


      ENDDO

      return
      end

      SUBROUTINE BMG3_SymStd_SETUP_MSG( &
     &                pMSG, pMSGSO, iMSG_Geom, NMSGi, pSI_MSG, IBC,&
     &                IGRD, NOG, NOGm, NProc, MyProc, &
     &                DimX, DimY, DimZ, DimXfine, DimYfine, DimZfine,&
     &                ProcGrid, NProcI, NProcJ, NProcK,&
     &                MPICOMM&
     &                ) BIND(C, NAME='BMG3_SymStd_SETUP_MSG')

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

! ------------------------------------------------
!     Includes
!
      include 'mpif.h'
      include 'MSG_f90.h'

      include 'BMG_constants_f90.h'
      include 'BMG_workspace_f90.h'

! ------------------------------------------------
!     Argument Declarations
!
      integer(c_int), value :: NOG, NOGm, NProcI, NProcJ, NProck, IBC
      integer, value :: MPICOMM
      integer(len_t), value :: NMSGi
      integer(c_int), value :: NProc, MyProc
      integer(c_int) :: pSI_MSG
      integer(len_t) :: IGRD(NOGm,NBMG_pIGRD), iMSG_Geom(NMSGi)
      integer(c_int) :: pMSG(NBMG_pMSG,NOG), pMSGSO(NBMG_pMSG,NOG)
      integer(len_t) :: DimX(NProcI,NOGm), DimY(NProcJ,NOGm), DimZ(NProcK,NOGm)
      integer(len_t) :: DimXfine(NProcI), DimYfine(NProcJ), DimZfine(NProcK)
      integer(c_int) :: ProcGrid(NProcI,NProcJ,NProcK)

! ------------------------------------------------
!     Local Declarations
!
      INTEGER N, NOG_cg, NOG_fg, ierror
      integer(len_t) :: iGs, jGs, kGs, I, J, K

! ======================================================================

      call MSG_set_comm_parent(MPICOMM)
      call MSG_enable(MyProc, NProc)
!
!     Create the pointers into the integer work space
!     describing the storage of MSG communication setup
!

!     Note that the pointer shift, pSI_MSG, is changed in every call

      NOG_fg = NOG   ! Index of the finest grid
      NOG_cg = 1     ! Index of the coarsest grid

      CALL BMG3_SymStd_SETUP_PtrMSG(&
     &                 IGRD(NOG,idL_BMG_NLx), IGRD(NOG,idL_BMG_NLy), &
     &                 IGRD(NOG,idL_BMG_NLz), iZERO, iZERO, iZERO, &
     &                 pSI_MSG, NProc, NOG, NOG_fg, NOG_cg, pMSG&
     &                 )

      CALL BMG3_SymStd_SETUP_PtrMSG(&
     &                 IGRD(NOG,idL_BMG_NLx), IGRD(NOG,idL_BMG_NLy),&
     &                 IGRD(NOG,idL_BMG_NLz), iONE, iONE, iONE, &
     &                 pSI_MSG, NProc, NOG, NOG_fg, NOG_cg, pMSGSO&
     &                 )


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

      DO K=1, NProcK
         DimZ(K,NOG) = DimZfine(K)
      END DO

      DO N=NOG-1, 1, -1

         iGs = 1
         !
         DO I=1, NProcI
            !
            IF ( MOD(iGs,2).EQ.1 ) THEN
               DimX(I,N) = (DimX(I,N+1)+1)/2
            ELSE
               IF ( MOD(DimX(I,N+1),2).EQ.1 ) THEN
                  DimX(I,N) = (DimX(I,N+1)-1)/2
               ELSE
                  DimX(I,N) = (DimX(I,N+1)+1)/2
               ENDIF
            ENDIF

            iGs = iGs + DimX(I,N+1)

         END DO

         jGs = 1
         !
         DO J=1, NProcJ
            !
            IF ( MOD(jGs,2).EQ.1 ) THEN
               DimY(J,N) = (DimY(J,N+1)+1)/2
            ELSE
               IF ( MOD(DimY(J,N+1),2).EQ.1 ) THEN
                  DimY(J,N) = (DimY(J,N+1)-1)/2
               ELSE
                  DimY(J,N) = (DimY(J,N+1)+1)/2
               ENDIF
            ENDIF
            !
            jGs = jGs + DimY(J,N+1)
            !
         END DO


         kGs = 1
         !
         DO K=1, NProcK
            !
            IF ( MOD(kGs,2).EQ.1 ) THEN
               DimZ(K,N) = (DimZ(K,N+1)+1)/2
            ELSE
               IF ( MOD(DimZ(K,N+1),2).EQ.1 ) THEN
                  DimZ(K,N) = (DimZ(K,N+1)-1)/2
               ELSE
                  DimZ(K,N) = (DimZ(K,N+1)+1)/2
               ENDIF
            ENDIF
            !
            kGs = kGs + DimZ(K,N+1)
            !
         END DO
         !
      END DO

!
!     create MSG grid information for SO
!


      DO N=NOG_fg, NOG_cg, -1

         CALL BMG3_SymStd_SETUP_MSGGrid(&
     &        IGRD(N,idL_BMG_NGx), IGRD(N,idL_BMG_NGy),&
     &        IGRD(N,idL_BMG_NGz),IBC,&
     &        iZERO, iZERO, iZERO, &
     &        iMSG_Geom(pMSG(ipL_MSG_LocalArraySize,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordLocalData,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_GlobalCoordActData,N)),&
     &        iMSG_Geom(pMSG(ipL_MSG_ActDataStart,N)),&
     &        DimX, DimY, DimZ, ProcGrid, &
     &        NProc, NProcI, NProcJ, NProcK, NOGm, N,&
     &        MPICOMM )

      ENDDO

      DO N=NOG_fg, NOG_cg, -1

         CALL BMG3_SymStd_SETUP_MSGGrid(&
     &        IGRD(N,idL_BMG_NGx), IGRD(N,idL_BMG_NGy),&
     &        IGRD(N,idL_BMG_NGz),IBC, &
     &        iONE, iONE, iONE, &
     &        iMSG_Geom(pMSGSO(ipL_MSG_LocalArraySize,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordLocalData,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_GlobalCoordActData,N)),&
     &        iMSG_Geom(pMSGSO(ipL_MSG_ActDataStart,N)),&
     &        DimX, DimY, DimZ, ProcGrid, &
     &        NProc, NProcI, NProcJ, NProcK, NOGm, N,&
     &        MPICOMM )

      ENDDO

!
!     Setup MSG communications for all grids
!

      DO N=NOG,1,-1


         iMSG_Geom(pMSG(ipL_MSG_NumAdjProc, N)) = 0

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

! ======================================================================

      RETURN
      END

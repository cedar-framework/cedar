subroutine BMG3_SymStd_UTILS_update_stencil_ghosts(&
     K, SO, II, JJ, KK, iWork, NMSGi, pMSGSO,&
     buffer, NMSGr, NOG&
     ) BIND(C, NAME='BMG3_SymStd_UTILS_update_stencil_ghosts')

  USE ModInterface
  IMPLICIT NONE

  INCLUDE 'mpif.h'
  INCLUDE 'MSG_f90.h'

  INCLUDE 'BMG_constants_f90.h'
  INCLUDE 'BMG_workspace_f90.h'



  integer(len_t), value :: II, JJ, KK, NMSGi, NMSGr
  integer(c_int), value :: NOG, K
  integer(len_t) :: iWork(NMSGi)
  integer(c_int) :: pMSGSO(NBMG_pMSG,NOG)
  real(real_t) :: SO(II+1,JJ+1,KK+1,14), BUFFER(NMSGr)

  integer :: kpz, ptrn, ierror


  DO kpz=1,14

     ptrn = 1
     call MSG_tbdx_send(SO(1,1,1,kpz), buffer, &
          &           iWork(pMSGSO(ipL_MSG_NumAdjProc,K)),&
          &           iWork(pMSGSO(ipL_MSG_Proc,K)),&
          &           iWork(pMSGSO(ipL_MSG_Ipr,K)),&
          &           iWork(pMSGSO(ipL_MSG_Index,K)),&
          &           ptrn, ierror)

     call MSG_tbdx_receive(SO(1,1,1,kpz), buffer,&
          &           iWork(pMSGSO(ipL_MSG_NumAdjProc,K)),&
          &           iWork(pMSGSO(ipL_MSG_Proc,K)),&
          &           iWork(pMSGSO(ipL_MSG_Ipr,K)),&
          &           iWork(pMSGSO(ipL_MSG_Index,K)),&
          &           ptrn, ierror)

     call MSG_tbdx_close(SO(1,1,1,kpz), buffer,&
          &           iWork(pMSGSO(ipL_MSG_NumAdjProc,K)),&
          &           iWork(pMSGSO(ipL_MSG_Proc,K)),&
          &           iWork(pMSGSO(ipL_MSG_Ipr,K)), &
          &           iWork(pMSGSO(ipL_MSG_Index,K)),&
          &           ptrn, ierror)

  ENDDO

endsubroutine BMG3_SymStd_UTILS_update_stencil_ghosts

      SUBROUTINE BMG3_SymStd_SETUP_ITLI07_ex(&
     &                KGF, KGC, SO, SOC, CI, IIF, JJF, KKF, IIC, JJC, &
     &                KKC, iGs, jGs, kGS, NOG, NOGm, IGRD,&
     &                iWORK, NMSGi, PMSGSO,&
     &                BUFFER, NMSGr, NProc, ProcGrid, MyProcI, MyProcJ,&
     &                MyProcK, NProcI, NProcJ, NProcK, DimX, DimY, DimZ,&
     &                MPICOMM&
     &                ) BIND(C, NAME='MPI_BMG3_SymStd_SETUP_ITLI07_ex')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     Construct the variational coarse-grid operator for the case of
!     a 7-point fine-grid discretization using the "block local"
!     technique.
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

! ----------------------------
!     Includes
!
      INCLUDE 'mpif.h'
      INCLUDE 'MSG_f90.h'

      INCLUDE 'BMG_constants_f90.h'
      INCLUDE 'BMG_stencils_f90.h'
      INCLUDE 'BMG_workspace_f90.h'
      INCLUDE 'BMG_parameters_f90.h'

! ----------------------------
!     Argument Declarations
!
      integer(len_t), value :: iGs, jGs, kGs, IIC, JJC, KKC,&
           IIF, JJF, KKF, NMSGi, NMSGr
      integer(c_int), value :: NOG, NOGm, NProc, NProcI, NProcJ,&
           NProcK, MyProcI, MyProcJ, MyProcK, MPICOMM,&
           KGF, KGC
      integer(len_t) :: iWork(NMSGi), IGRD(NOGm,29)
      integer(c_int) :: pMSGSO(NBMG_pMSG,NOGm), ProcGrid(NProcI,NProcJ, NProcK)
      real(real_t) :: CI(IIC,JJC,KKC,26), SO(IIF+1,JJF+1,KKF+1,4),&
           SOC(IIC+1,JJC+1,KKC+1,14), BUFFER(NMSGr)
      integer(c_int) :: DimX(NprocI,NOGm), DimY(NprocJ,NOGm), DimZ(Nprock,NOGm)

! ----------------------------
!     Local Declarations:
!
      INTEGER I, IIC1, IIC2, IC, J, JJC1, JJC2, JC, K, KC, KKC1, kpz,&
     &        ptrn, ierror
      REAL*8  CB, CBE, CBN, CBNW, CBNE, CBS, CBSE, CBSW, CBW,    &
     &        CE, CN, CNW, CNE, CO, CS, CSE, CSW, CW,&
     &        CTE, CTN, CTNW, CTNE, CT, CTS, CTSE, CTSW, CTW

      INTEGER LXGP, RXGP, LYGP, RYGP, LZGP, RZGP
      INTEGER MyProc, ISTART, JSTART, KSTART, MPI_IERROR

! ======================================================================

      MyProc = MSG_MyProc(MPICOMM)

      IIC1 = IIC-1
      IIC2 = IIC-2
      JJC1 = JJC-1
      JJC2 = JJC-2
      KKC1 = KKC-1

      DO KPZ=1,14
         DO k=1,KKC+1
            DO j=1,JJC+1
               DO i=1,IIC+1
                  soc(i,j,k,KPZ)= rZERO
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      IF (MOD(iGs,2).EQ.1) THEN
         ISTART = 0
      ELSE
         ISTART = 1
      ENDIF

      IF (MOD(jGs,2).EQ.1) THEN
         JSTART = 0
      ELSE
         JSTART = 1
      ENDIF

      IF(MOD(kGs,2).EQ.1) THEN
         KSTART = 0
      ELSE
         KSTART = 1
      ENDIF

      K = KSTART

      DO KC = 2,KKC1

         K  = K+2
         J = JSTART

         DO JC = 2,JJC1
            J = J+2
            I = ISTART

            DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KPW)*CI(IC,JC,KC,LXYL)
               CW = - SO(I-1,J,K,KP)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K,KPW)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZNW)
               CNW = -SO(I-1,J+1,K,KP)*CI(IC,JC+1,KC,LXYSW)&
     &              + SO(I-1,J+1,K,KPW)*CI(IC-1,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LTSW)
               CTNW = -SO(I-1,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J+1,K+1,KPW)*CI(IC-1,JC+1,KC+1,LYZSE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYSW)
               CTW = - SO(I-1,J,K+1,KP)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K+1,KPW)*CI(IC-1,JC,KC+1,LXZB)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LBSW)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBNW)
               CTSW = -SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LXZSW)&
     &              + SO(I-1,J-1,K+1,KPW)*CI(IC-1,JC,KC+1,LYZSW)
               CSW = - SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J-1,K,KPW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTNW)
               CBNW = -SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I-1,J+1,K-1,KPW)*CI(IC-1,JC+1,KC,LYZNE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LXYSW)
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J,K-1,KPW)*CI(IC-1,JC,KC,LXZA)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LTSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTNW)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZNW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZNW)
               CT = SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZSW)
               CTN = SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LBSW)
               CTS = SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBNW)
               CN = SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYSW)
               CS = SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYNW)
               CB = SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZNW)
               CBN = SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LTSW)
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTNW)
               SOC(IC,JC,KC,KPW) = CO + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW&
     &              + CI(IC,JC+1,KC+1,LBSE)*CTNW&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW &
     &              + CI(IC,JC,KC+1,LBNE)*CTSW&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC+1,LXZB)*CT&
     &              + CI(IC,JC+1,KC+1,LYZSE)*CTN&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS&
     &              + CI(IC,JC+1,KC,LXYB)*CN&
     &              + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LYZNW)*CBS
              ENDDO
            ENDDO
         ENDDO

         K= KSTART

         DO KC = 2,KKC1
            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CW = SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYSW)
               CS = SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYSW)
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J-1,K,KPW)*CI(IC-1,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K,KPS)*CI(IC,JC-1,KC,LXYL)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBSW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTSW)
               CTW = SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBSW)
               CTS = SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBSW)
               CTSW = SO(I-1,J-1,K+1,KPW)*CI(IC-1,JC,KC+1,LYZSE)&
     &              + SO(I-1,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZSW)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYSW)&
     &              - SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBSW)
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTSW)
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTSW)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTSW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYSW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZNE)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZNW)
               SOC(IC,JC,KC,KPSW) = CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC,JC,KC,LTNE)*CBSW
             ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBSW)
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBSW)
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBSW)
               CBSW = SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZSE)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBSW)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZSW)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYSW)
               SOC(IC,JC,KC,KBSW) = CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = - SO(I,J,K,KP)&
     &              + SO(I,J,K,KPW)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K,KPS)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J,K,KPW)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I,J,K,KPS)*CI(IC,JC,KC,LXYA)&
     &              + SO(I,J,K,KB)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J,K+1,KB)*CI(IC,JC,KC+1,LXZB)
               CW = - SO(I-1,J,K,KP)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J,K,KPW)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZNE)
               CNW = -SO(I-1,J+1,K,KP)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I-1,J+1,K,KPS)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LBSE)
               CN = - SO(I,J+1,K,KP)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I,J+1,K,KPS)&
     &              + SO(I,J+1,K,KPW)*CI(IC,JC+1,KC,LXYSE)&
     &              + SO(I+1,J+1,K,KPW)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J+1,K+1,KB)*CI(IC,JC+1,KC+1,LYZSE)
               CNE = - SO(I+1,J+1,K,KP)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J+1,K,KPW)*CI(IC,JC+1,KC,LXYB)&
     &              + SO(I+1,J+1,K,KPS)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J+1,K+1,KB)*CI(IC+1,JC+1,KC+1,LBSW)
               CE = -SO(I+1,J,K,KP)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J,K,KPW)&
     &              + SO(I+1,J+1,K,KPS)*CI(IC+1,JC+1,KC,LXYSW)&
     &              + SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J,K+1,KB)*CI(IC+1,JC,KC+1,LXZSW)
               CSE = - SO(I+1,J-1,K,KP)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYL)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC+1,LBNW)
               CS = -SO(I,J-1,K,KP)*CI(IC,JC,KC,LXYA)&
     &              + SO(I,J,K,KPS)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC+1,JC,KC,LXYNW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC+1,LYZSW)
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYR)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYA)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTNE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBNE)
               CB = -SO(I,J,K-1,KP)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J,K,KB)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZNW)
               CBW = SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTNE)&
     &              + SO(I-1,J,K,KB)*CI(IC,JC,KC,LXYR)&
     &              - SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZNE)
               CBNW = -SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LXYSE)
               CBN = -SO(I,J+1,K-1,KP)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LTSE)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC,KC,LXZA)&
     &              + SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LXYB)
               CBNE = -SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZNE)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LXYSW)
               CBE = -SO(I+1,J,K-1,KP)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LTSW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXYL)
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZNW)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LXYNW)
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTNE)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LXZA)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LTNW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LXYA)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTNE)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZNE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYNE)
               CT = -SO(I,J,K+1,KP)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J,K+1,KB)&
     &              + SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J,K+1,KPW)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LYZSW)
               CTW = -SO(I-1,J,K+1,KP)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J,K+1,KPW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I-1,J,K+1,KB)*CI(IC,JC,KC,LXYR)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBNE)
               CTNW = - SO(I-1,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I-1,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I-1,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYSE)
               CTN = -SO(I,J+1,K+1,KP)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LBSE)&
     &              + SO(I+1,J+1,K+1,KPW)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I,J+1,K+1,KPS)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I,J+1,K+1,KB)*CI(IC,JC+1,KC,LXYB)
               CTNE = -SO(I+1,J+1,K+1,KP)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J+1,K+1,KPW)*CI(IC,JC+1,KC+1,LYZSE)&
     &              + SO(I+1,J+1,K+1,KPS)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J+1,K+1,KB)*CI(IC+1,JC+1,KC,LXYSW)
               CTE = -SO(I+1,J,K+1,KP)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J,K+1,KPW)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J+1,K+1,KPS)*CI(IC+1,JC+1,KC+1,LBSW)&
     &              + SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I+1,J,K+1,KB)*CI(IC+1,JC,KC,LXYL)
               CTSE = -SO(I+1,J-1,K+1,KP)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC,LXYNW)
               CTS = -SO(I,J-1,K+1,KP)*CI(IC,JC,KC+1,LYZSW)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LXZB)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC+1,JC,KC+1,LBNW)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC,LXYA)
               CTSW = -SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBNE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYNE)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LXZSE)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSW)
               SOC(IC,JC,KC,KP) = - CO - CI(IC,JC,KC,LXYR)*CW&
     &              - CI(IC,JC+1,KC,LXYSE)*CNW&
     &              - CI(IC,JC+1,KC,LXYB)*CN&
     &              - CI(IC+1,JC+1,KC,LXYSW)*CNE&
     &              - CI(IC+1,JC,KC,LXYL)*CE&
     &              - CI(IC+1,JC,KC,LXYNW)*CSE&
     &              - CI(IC,JC,KC,LXYA)*CS&
     &              - CI(IC,JC,KC,LXYNE)*CSW&
     &              - CI(IC,JC,KC,LXZA)*CB&
     &              - CI(IC,JC,KC,LXZNE)*CBW&
     &              - CI(IC,JC+1,KC,LTSE)*CBNW&
     &              - CI(IC,JC+1,KC,LYZNE)*CBN&
     &              - CI(IC+1,JC+1,KC,LTSW)*CBNE&
     &              - CI(IC+1,JC,KC,LXZNW)*CBE&
     &              - CI(IC+1,JC,KC,LTNW)*CBSE&
     &              - CI(IC,JC,KC,LYZNW)*CBS&
     &              - CI(IC,JC,KC,LTNE)*CBSW&
     &              - CI(IC,JC,KC+1,LXZB)*CT&
     &              - CI(IC,JC,KC+1,LXZSE)*CTW&
     &              - CI(IC,JC+1,KC+1,LBSE)*CTNW&
     &              - CI(IC,JC+1,KC+1,LYZSE)*CTN&
     &              - CI(IC+1,JC+1,KC+1,LBSW)*CTNE&
     &              - CI(IC+1,JC,KC+1,LXZSW)*CTE&
     &              - CI(IC+1,JC,KC+1,LBNW)*CTSE&
     &              - CI(IC,JC,KC+1,LYZSW)*CTS&
     &              - CI(IC,JC,KC+1,LBNE)*CTSW
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CBW = SO(I-1,J-1,K-1,KPS)*CI(IC,JC,KC,LTNW)
               CW = SO(I-1,J-1,K,KPS)*CI(IC,JC,KC,LXYNW)
               CBNW = SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZNW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZNW)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYNW)
               CBN = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTNW)
               CN = SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYNW)
               CNW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYNW)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTNW)&
     &              + SO(I-1,J-1,K,KPW)*CI(IC-1,JC,KC,LXYA)&
     &              + SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYL)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBNW)
               CTW = SO(I-1,J-1,K+1,KPS)*CI(IC,JC,KC+1,LBNW)
               CTNW = SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYNW)&
     &              - SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBNW)&
     &              + SO(I-1,J-1,K+1,KPW)*CI(IC-1,JC,KC+1,LYZSW)&
     &              + SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LXZSW)
               CTN = SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBNW)
               SOC(IC,JC,KC,KPNW) = CI(IC,JC-1,KC,LXZNE)*CBW&
     &              + CI(IC,JC-1,KC,LXYR)*CW&
     &              + CI(IC,JC,KC,LTSE)*CBNW&
     &              + CI(IC,JC,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LXYB)*CN&
     &              + CI(IC,JC,KC,LXYSE)*CNW&
     &              + CI(IC,JC-1,KC+1,LXZSE)*CTW&
     &              + CI(IC,JC,KC+1,LBSE)*CTNW&
     &              + CI(IC,JC,KC+1,LYZSE)*CTN
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CW = SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZSW)
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNW)
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J,K-1,KPW)*CI(IC-1,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J,K-1,KB)*CI(IC,JC,KC-1,LXYL)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBSW)
               CBSW = SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSW)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNW)
               CNW = SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBSW)
               CBNW = SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J+1,K-1,KPW)*CI(IC-1,JC+1,KC,LYZSE)&
     &              - SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBSW)&
     &              + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYSW)
               CB = SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZSW)
               CBN = SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBSW)
               CBS = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNW)
               SOC(IC,JC,KC,KBW) = CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LYZNW)*CBS
               ENDDO
             ENDDO
          ENDDO

          K = KSTART

          DO KC = 2,KKC1

             K = K+2
             J = JSTART

             DO JC = 2,JJC1

                J = J+2
                I = ISTART

                DO IC = 2,IIC1
                I = I+2
                CE = SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZSE)
                CSE = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)
                CBE = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZSE)&
     &               + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZB)&
     &               + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBNE)&
     &               + SO(I-1,J,K-1,KB)*CI(IC,JC,KC-1,LXYR)&
     &               + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBSE)
                CBSE = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &               + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)&
     &               + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)&
     &               - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)
                CNE = SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBSE)
                CBN = SO(I-1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBSE)
                CB = SO(I-1,J,K-1,KPW)*CI(IC,JC,KC,LXZSE)
                CBS = SO(I-1,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)
                CBNE = SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &               + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSE)&
     &               - SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBSE)&
     &               + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYSE)
                SOC(IC,JC,KC,KBE) = CI(IC,JC,KC,LXYL)*CE&
     &               + CI(IC,JC,KC,LXYNW)*CSE&
     &               + CI(IC,JC,KC,LXZNW)*CBE&
     &               + CI(IC,JC,KC,LTNW)*CBSE&
     &               + CI(IC,JC+1,KC,LXYSW)*CNE&
     &               + CI(IC-1,JC+1,KC,LYZNE)*CBN&
     &               + CI(IC-1,JC,KC,LXZA)*CB&
     &               + CI(IC-1,JC,KC,LYZNW)*CBS&
     &               + CI(IC,JC+1,KC,LTSW)*CBNE
                ENDDO
             ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CBSW = SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZSE)&
     &              - SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYSE)
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBSE)
               CS = SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZSE)
               CSE = SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBSW)
               CBSE = SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSE)&
     &              - SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBSW)&
     &              + SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC-1,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYSW)
               CBS = SO(I,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZB)&
     &              - SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYB)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBSE)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBSW)
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBSE)
               CB = SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZSE)
               CBE = SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LBSW)
               SOC(IC,JC,KC,KBS) = CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J= J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CBE = SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC,KC,LBNW)
               CB = SO(I,J-1,K-1,KPS)*CI(IC,JC,KC,LYZSW)
               CBW = SO(I-1,J-1,K-1,KPS)*CI(IC,JC,KC,LBNE)
               CNW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)
               CN = SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZSW)
               CNE = SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBNW)
               CBNE = SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYNW)&
     &              - SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)
               CBN = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)&
     &              - SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LXZB)
               CBNW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)
               SOC(IC,JC,KC,KBN) = CI(IC+1,JC-1,KC,LXZNW)*CBE&
     &              + CI(IC,JC-1,KC,LXZA)*CB&
     &              + CI(IC,JC-1,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LXYSE)*CNW&
     &              + CI(IC,JC,KC,LXYB)*CN&
     &              + CI(IC+1,JC,KC,LXYSW)*CNE&
     &              + CI(IC+1,JC,KC,LTSW)*CBNE&
     &              + CI(IC,JC,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LTSE)*CBNW
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CNE = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)
               CBN = SO(I-1,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)
               CBNE = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)
               CBE = SO(I-1,J-1,K-1,KPS)*CI(IC,JC,KC,LBNE)
               SOC(IC,JC,KC,KBNE) = CI(IC,JC,KC,LXYSW)*CNE&
     &              + CI(IC-1,JC,KC,LYZNE)*CBN&
     &              + CI(IC,JC,KC,LTSW)*CBNE&
     &              + CI(IC,JC-1,KC,LXZNW)*CBE
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CNW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNW)
               CBW = SO(I-1,J-1,K-1,KPS)*CI(IC,JC,KC,LBNW)
               CBNW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNW)&
     &              + SO(I-1,J-1,K-1,KPW)*CI(IC-1,JC,KC,LYZSW)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSW)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNW)
               CBN = SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNW)
               SOC(IC,JC,KC,KBNW) = CI(IC,JC,KC,LXYSE)*CNW&
     &              + CI(IC,JC-1,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LTSE)*CBNW&
     &              + CI(IC,JC,KC,LYZNE)*CBN
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CSE = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBSE)
               CBE = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBSE)
               CBSE = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBSE)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYSE)
               CBS = SO(I-1,J-1,K-1,KPW)*CI(IC,JC,KC,LBSE)
               SOC(IC,JC,KC,KBSE) = CI(IC,JC,KC,LXYNW)*CSE&
     &              + CI(IC,JC,KC,LXZNW)*CBE&
     &              + CI(IC,JC,KC,LTNW)*CBSE&
     &              + CI(IC-1,JC,KC,LYZNW)*CBS
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1

               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KB)*CI(IC,JC,KC,LXZB)
               CB = -SO(I,J,K-1,KP)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J,K-1,KB)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZSW)
               CBW = -SO(I-1,J,K-1,KP)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J,K-1,KPW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I-1,J,K-1,KB)*CI(IC,JC,KC-1,LXYR)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LBNE)
               CBNW =-SO(I-1,J+1,K-1,KP)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I-1,J+1,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I-1,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYSE)
               CBN = -SO(I,J+1,K-1,KP)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I,J+1,K-1,KPW)*CI(IC,JC+1,KC,LBSE)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I,J+1,K-1,KPS)*CI(IC,JC,KC,LXZB)&
     &              + SO(I,J+1,K-1,KB)*CI(IC,JC+1,KC-1,LXYB)
               CBNE = -SO(I+1,J+1,K-1,KP)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J+1,K-1,KB)*CI(IC+1,JC+1,KC-1,LXYSW)&
     &              + SO(I+1,J+1,K-1,KPW)*CI(IC,JC+1,KC,LYZSE)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC,KC,LXZSW)
               CBE = -SO(I+1,J,K-1,KP)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J,K-1,KPW)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J+1,K-1,KPS)*CI(IC+1,JC+1,KC,LBSW)&
     &              + SO(I+1,J,K-1,KB)*CI(IC+1,JC,KC-1,LXYL)
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LBNW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LXZSW)&
     &              + SO(I+1,J-1,K-1,KB)*CI(IC+1,JC,KC-1,LXYNW)
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZSW)&
     &              + SO(I,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYA)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LBNE)&
     &              + SO(I,J,K-1,KPS)*CI(IC,JC,KC,LXZB)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LBNW)
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LBNE)&
     &              + SO(I-1,J-1,K-1,KB)*CI(IC,JC,KC-1,LXYNE)&
     &              + SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LXZSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZSW)
               CW = SO(I-1,J,K,KB)*CI(IC,JC,KC,LXZSE)
               CNW = SO(I-1,J+1,K,KB)*CI(IC,JC+1,KC,LBSE)
               CN = SO(I,J+1,K,KB)*CI(IC,JC+1,KC,LYZSE)
               CNE = SO(I+1,J+1,K,KB)*CI(IC+1,JC+1,KC,LBSW)
               CE = SO(I+1,J,K,KB)*CI(IC+1,JC,KC,LXZSW)
               CSE = SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LBNW)
               CS = SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZSW)
               CSW = SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LBNE)
               SOC(IC,JC,KC,KB) = CO + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC+1,KC,LTSE)*CBNW&
     &              + CI(IC,JC+1,KC,LYZNE)*CBN&
     &              + CI(IC+1,JC+1,KC,LTSW)*CBNE&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC+1,KC,LXYSE)*CNW&
     &              + CI(IC,JC+1,KC,LXYB)*CN&
     &              + CI(IC+1,JC+1,KC,LXYSW)*CNE&
     &              + CI(IC+1,JC,KC,LXYL)*CE&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE&
     &              + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               ENDDO
            ENDDO
         ENDDO

         K = KSTART

         DO KC = 2,KKC1

            K = K+2
            J = JSTART

            DO JC = 2,JJC1
               J = J+2
               I = ISTART

               DO IC = 2,IIC1
               I = I+2
               CO = SO(I,J,K,KPS)*CI(IC,JC,KC,LXYB)
               SOC(IC,JC,KC,KPS) = CO
               CS = -SO(I,J-1,K,KP)*CI(IC,JC,KC,LXYB)&
     &              + SO(I,J-1,K,KPS)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC+1,LYZSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXYA)*CS
               CSW = -SO(I-1,J-1,K,KP)*CI(IC,JC,KC,LXYSE)&
     &              + SO(I,J-1,K,KPW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I-1,J-1,K,KPS)*CI(IC,JC-1,KC,LXYR)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LTSE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC+1,LBSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXYNE)*CSW
               CTSW = - SO(I-1,J-1,K+1,KP)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I-1,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZSE)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I-1,J-1,K+1,KB)*CI(IC,JC,KC,LXYSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW
               CTS = -SO(I,J-1,K+1,KP)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I,J-1,K+1,KB)*CI(IC,JC,KC,LXYB)&
     &              + SO(I,J-1,K+1,KPS)*CI(IC,JC-1,KC+1,LXZB)&
     &              + SO(I,J-1,K+1,KPW)*CI(IC,JC,KC+1,LBSE)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC+1,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS
               CTSE = -SO(I+1,J-1,K+1,KP)*CI(IC+1,JC,KC+1,LBSW)&
     &              + SO(I+1,J-1,K+1,KPW)*CI(IC,JC,KC+1,LYZSE)&
     &              + SO(I+1,J-1,K+1,KPS)*CI(IC+1,JC-1,KC+1,LXZSW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC+1,LBNW)*CTSE
               CBSW = -SO(I-1,J-1,K-1,KP)*CI(IC,JC,KC,LTSE)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I-1,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZNE)&
     &              + SO(I-1,J-1,K,KB)*CI(IC,JC,KC,LXYSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LTNE)*CBSW
               CBS = -SO(I,J-1,K-1,KP)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I,J-1,K-1,KPS)*CI(IC,JC-1,KC,LXZA)&
     &              + SO(I,J-1,K-1,KPW)*CI(IC,JC,KC,LTSE)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I,J-1,K,KB)*CI(IC,JC,KC,LXYB)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LYZNW)*CBS
               CBSE = -SO(I+1,J-1,K-1,KP)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J-1,K-1,KPW)*CI(IC,JC,KC,LYZNE)&
     &              + SO(I+1,J-1,K-1,KPS)*CI(IC+1,JC-1,KC,LXZNW)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE
               CW = SO(I-1,J,K,KPS)*CI(IC,JC,KC,LXYSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXYR)*CW
               CTW = SO(I-1,J,K+1,KPS)*CI(IC,JC,KC+1,LBSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW
               CT = SO(I,J,K+1,KPS)*CI(IC,JC,KC+1,LYZSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC+1,LXZB)*CT
               CTE = SO(I+1,J,K+1,KPS)*CI(IC+1,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC+1,LXZSW)*CTE
               CBW = SO(I-1,J,K-1,KPS)*CI(IC,JC,KC,LTSE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXZNE)*CBW
               CB = SO(I,J,K-1,KPS)*CI(IC,JC,KC,LYZNE)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC,JC,KC,LXZA)*CB
               CBE = SO(I+1,J,K-1,KPS)*CI(IC+1,JC,KC,LTSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE
               CE = SO(I+1,J,K,KPS)*CI(IC+1,JC,KC,LXYSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC,LXYL)*CE
               CSE = -SO(I+1,J-1,K,KP)*CI(IC+1,JC,KC,LXYSW)&
     &              + SO(I+1,J-1,K,KPW)*CI(IC,JC,KC,LXYB)&
     &              + SO(I+1,J-1,K,KPS)*CI(IC+1,JC-1,KC,LXYL)&
     &              + SO(I+1,J-1,K,KB)*CI(IC+1,JC,KC,LTSW)&
     &              + SO(I+1,J-1,K+1,KB)*CI(IC+1,JC,KC+1,LBSW)
               SOC(IC,JC,KC,KPS) = SOC(IC,JC,KC,KPS)&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE
               SOC(IC,JC,KC,KPS) = CO + CI(IC,JC,KC,LXYA)*CS&
     &              + CI(IC,JC,KC,LXYNE)*CSW&
     &              + CI(IC,JC,KC+1,LBNE)*CTSW&
     &              + CI(IC,JC,KC+1,LYZSW)*CTS&
     &              + CI(IC+1,JC,KC+1,LBNW)*CTSE&
     &              + CI(IC,JC,KC,LTNE)*CBSW&
     &              + CI(IC,JC,KC,LYZNW)*CBS&
     &              + CI(IC+1,JC,KC,LTNW)*CBSE&
     &              + CI(IC,JC,KC,LXYR)*CW&
     &              + CI(IC,JC,KC+1,LXZSE)*CTW&
     &              + CI(IC,JC,KC+1,LXZB)*CT&
     &              + CI(IC+1,JC,KC+1,LXZSW)*CTE&
     &              + CI(IC,JC,KC,LXZNE)*CBW&
     &              + CI(IC,JC,KC,LXZA)*CB&
     &              + CI(IC+1,JC,KC,LXZNW)*CBE&
     &              + CI(IC+1,JC,KC,LXYL)*CE&
     &              + CI(IC+1,JC,KC,LXYNW)*CSE
              ENDDO
            ENDDO
         ENDDO

      DO kpz=1,14


            ptrn = 1
            call MSG_tbdx_send(SOC(1,1,1,kpz), buffer, &
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_receive(SOC(1,1,1,kpz), buffer,&
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)),&
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)

            call MSG_tbdx_close(SOC(1,1,1,kpz), buffer,&
     &           iWork(pMSGSO(ipL_MSG_NumAdjProc,KGC)),&
     &           iWork(pMSGSO(ipL_MSG_Proc,KGC)),&
     &           iWork(pMSGSO(ipL_MSG_Ipr,KGC)), &
     &           iWork(pMSGSO(ipL_MSG_Index,KGC)),&
     &           ptrn, ierror)



      ENDDO

! ======================================================================

      RETURN
      END

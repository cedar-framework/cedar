SUBROUTINE BMG2_SymStd_SETUP_TAUSCH(&
     &                NOG, &
     &                DimX, DimY, DimXfine, DimYfine,&
     &                NProcI, NProcJ &
     &                ) BIND(C, NAME='BMG2_SymStd_SETUP_Tausch')

  use ModInterface
  IMPLICIT NONE

  integer(C_INT), VALUE :: NOG, NProcI, NProcJ
  INTEGER(len_t) :: DimX(NProcI,NOG), DimY(NProcJ,NOG)
  INTEGER(len_t) :: DimXfine(NProcI), DimYfine(NProcJ)

  integer(len_t) :: I, J, iGs, jGs
  integer N

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

end SUBROUTINE BMG2_SymStd_SETUP_TAUSCH

SUBROUTINE BMG3_SymStd_SETUP_TAUSCH(&
     NOG,&
     DimX, DimY, Dimz, DimXfine, DimYfine, DimZfine,&
     NProcI, NProcJ, NProcK &
     ) BIND(C, NAME='BMG3_SymStd_SETUP_Tausch')

  use ModInterface
  IMPLICIT NONE

  integer(C_INT), VALUE :: NOG, NProcI, NProcJ, NProcK
  INTEGER(len_t) :: DimX(NProcI,NOG), DimY(NProcJ,NOG), DimZ(NProcK,NOG)
  INTEGER(len_t) :: DimXfine(NProcI), DimYfine(NProcJ), DimZfine(NProcK)

  integer(len_t) :: I, J, K, iGs, jGs, kGs
  integer N

  DO I=1, NProcI
     DimX(I,NOG) = DimXfine(I)
  END DO

  DO J=1, NProcJ
     DimY(J,NOG) = DimYfine(J)
  END DO

  do k=1,nprock
     dimy(k,nog) = dimzfine(k)
  enddo

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

     kgs = 1

     do k=1,nprock

        if (mod(kgs,2) .eq. 1) then
           dimz(k,n) = (dimz(k,n+1)+1) / 2
        else
           if (mod(dimz(k,n+1),2) .eq. 1) then
              dimz(k,n) = (dimz(k,n+1)-1)/2
           else
              dimz(k,n) = (dimz(k,n+1)+1)/2
           endif
        endif

        kgs = kgs + dimz(k,n+1)
     enddo

  END DO

end SUBROUTINE BMG3_SymStd_SETUP_TAUSCH

      SUBROUTINE BMG3_SymStd_SETUP_nog( &
     &                NLx, NLy, NLz, NLXYZc, NOG, &
     &                iGs, jGs, kGs, MPICOMM&
     &                ) BIND(C, NAME='BMG3_SymStd_SETUP_nog')

! ======================================================================
!  --------------------
!   DESCRIPTION:
!  --------------------
!
!     BMG3_SymStd_SETUP_nog computes the number of grids (NOG) based on
!     a local processor coarsening limit NLXYZc.
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

! ---------------------------
!     Includes:
!
      INCLUDE 'mpif.h'

! ---------------------------
!     Argument Declarations:
!
      integer(len_t), value :: iGs, jGs, kGs,&
           NLx, NLy, NLz, NLXYZc
      integer(c_int), value :: MPICOMM
      integer(c_int) :: NOG

! ---------------------------
!     Local Declarations:
!
      INTEGER  iGs_c, jGs_c, kGs_c, kg, &
     &         NLx_c, NLy_c, NLz_c, ierror

! ======================================================================

! ======================================================================
! -----------------------------------
!     Compute the number of grids:
! -----------------------------------

      !
      ! Initilize coarse grid counter (local)
      !
      kg = 1
      !
      ! Initialize coarse grid dimensions (local)
      !
      iGs_c = iGs
      jGs_c = jGs
      kGs_c = kGs
      NLx_c = NLx
      NLy_c = NLy
      NLz_c = NLz

 100  CONTINUE     ! >>>>>>>>>>>> LOOP BOUNDARY: computing NOG

         IF ( MOD(iGs_c,2).EQ.1 ) THEN
            iGs_c = (iGs_c+1)/2
            NLx_c = (NLx_c+1)/2
         ELSE
            iGs_c = iGs_c/2+1
            IF ( MOD(NLx_c,2).EQ.1 ) THEN
               NLx_c = (NLx_c-1)/2
            ELSE
               NLx_c = (NLx_c+1)/2
            ENDIF
         ENDIF

         IF ( MOD(jGs_c,2).EQ.1 ) THEN
            jGs_c = (jGs_c+1)/2
            NLy_c = (NLy_c+1)/2
         ELSE
            jGs_c = jGs_c/2+1
            IF ( MOD(NLy_c,2).EQ.1 ) THEN
               NLy_c = (NLy_c-1)/2
            ELSE
               NLy_c = (NLy_c+1)/2
            ENDIF
         ENDIF

         IF ( MOD(kGs_c,2).EQ.1 ) THEN
            kGs_c = (kGs_c+1)/2
            NLz_c = (NLz_c+1)/2
         ELSE
            kGs_c = kGs_c/2+1
            IF ( MOD(NLz_c,2).EQ.1 ) THEN
               NLz_c = (NLz_c-1)/2
            ELSE
               NLz_c = (NLz_c+1)/2
            ENDIF
         ENDIF

         IF ( MIN(NLx_c,NLy_c,NLz_c).GE.NLXYZc ) then
            kg=kg+1
            GOTO 100
         endif

 110  CONTINUE     ! >>>>>>>>>>>> LOOP BOUNDARY: computing NOG

      !
      !  Global minimum is stored in NOG
      !
      CALL MPI_Allreduce( kg, NOG, 1,&
     &                    MPI_INTEGER, MPI_MIN, MPICOMM, ierror &
     &                   )

! ======================================================================

 500  FORMAT (/,'FATAL ERROR: BMG3_SymStd_SETUP_nog',/,5X,A)
 510  FORMAT (5X,A,I7)
 520  FORMAT (5X,A,I7,/)
 530  FORMAT (/,2X,I1,1X,A,/)

! ===========================================

      RETURN
      END

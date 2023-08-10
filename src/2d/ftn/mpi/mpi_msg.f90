!________________________________________________________________
!
!
!   Message-Passing Toolkit for Structured Grid Communications
!
!                     Communications Routines
!
!---------------------------------------------------------------
!
!   Contents:
!-----------
!
!
!
!   MSG_enable - gives the configuration information
!
!   MSG_disable - disables the message-passing environment
!
!   MSG_nproc - gives the number of processors involved
!
!   MSG_myproc - gives my processor number (from 1 to MAX_PROCS)
!
!   MSG_set_comm_parent - sets the parent communicator
!
!   MSG_tbdx_send - boundary information exchange routine which makes use of
!                   repeated communication patterns (a "half channel" option)
!                   for a tensor product grid;
!                   it sends the boundary data out
!
!   MSG_tbdx_receive - receives the boundary data
!
!   MSG_tbdx_close - closes the communication pattern (channel)
!
!   MSG_tbdx_gather - gathers boundary data in a single buffer
!
!   MSG_tbdx_scatter - places the boundary data into the local array
!
!---------------------------------------------------------------
!
!   written by A. Malevsky, last modified:  Aug 8, 2023
!
!   record of updates to this version:
!
!       MSG_COMM_PARENT_FLAG added to mpi_param_f90.h
!       to signal a modification of the default parent communicator
!       in odrer to fix the problem with some compilers whcih do not
!       take multiple initialisations with DATA in common blocks
!       May 14, 1997
!
!       MSG_tbdx_send, MSG_tbdx_receive, and MSG_tbdx_close
!       exit immediately if the number of adjacent processors is 0
!       May 16, 1997
!
!       Source updated for freeform Fortran 90/95
!       Aug 8, 2023
!
!________________________________________________________________
!
!

subroutine MSG_enable(MyProc, NumProc)
  !-----------------------------------------------------------------
  !
  !   MSG_enable initiates a separate message-passing environment
  !   for the MSG data transfer functions and obtains information
  !   on the number of processors in the system and the current
  !   the current processor's number
  !
  !____________________________ OUTPUT _____________________________
  !
  !   MyProc - my processor number (from 1 to NumProc)
  !   NumProc - number of available processors
  !
  !-----------------------------------------------------------------
  implicit none
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  include 'MSG_f90.h'
  integer ierror, MyProc, NumProc
  logical flag
  MSG_VERSION = 2
  !
  !     initialize the MPI if it has not been already initialized
  !

  MSG_BLOCKING = 1

  call MPI_initialized(flag, ierror)
  if(.not.flag) call MPI_init(ierror)
  !
  !     create a new communicator to separate the MSG communications
  !
  if(MSG_COMM_PARENT_FLAG.ne.MSG_COMM_PARENT_MODIFIED) then
     !        no parent communicator has been specified
     MSG_COMM_PARENT = MPI_COMM_WORLD
     call MPI_COMM_DUP(MPI_COMM_WORLD, MSG_COMM, ierror)
  else
     !        a parent communicator has been specified
     call MPI_COMM_DUP(MSG_COMM_PARENT, MSG_COMM, ierror)
  endif
  MSG_COMM_PARENT_FLAG = 0
  call MPI_COMM_RANK(MSG_COMM, MyProc, ierror)
  MyProc = MyProc + 1
  call MPI_COMM_SIZE(MSG_COMM, NumProc, ierror)

  return
end subroutine MSG_enable


subroutine MSG_set_comm_parent(comm)
  implicit none
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  integer comm
  MSG_COMM_PARENT_FLAG = MSG_COMM_PARENT_MODIFIED
  MSG_COMM_PARENT = comm
  return
end subroutine MSG_set_comm_parent


integer function MSG_myproc()
  implicit none
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  integer ierror, MyProc
  call MPI_COMM_RANK(MSG_COMM, MyProc, ierror)
  MSG_myproc = MyProc + 1
  return
end function MSG_myproc


integer function MSG_nproc()
  implicit none
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  integer ierror, NumProc
  call MPI_COMM_SIZE(MSG_COMM, NumProc, ierror)
  MSG_nproc = NumProc
  return
end function MSG_nproc



subroutine MSG_comm_type (T)
  implicit none
  INTEGER T

  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'

  IF (T .LT. 0  .OR. T .GT. 1) THEN
     print *,'ERROR in MSG_comm_type(T)'
     print *,'  T = ',T
     print *,'  T must be either 0 or 1!'
     STOP
  ELSE
     MSG_BLOCKING = T
  ENDIF

  RETURN
END subroutine MSG_comm_type







subroutine MSG_disable(ierror)
  implicit none
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  integer ierror
  call MPI_BARRIER(MSG_COMM, ierror)
  call MPI_COMM_FREE(MSG_COMM, ierror)
  return
end subroutine MSG_disable





subroutine MSG_tbdx_send(x,y,nproc,proc,ipr,index,ptrn,ierr)
  implicit none
  integer nproc,ptrn,ierr
  integer proc(nproc), ipr(*), index(*)
  REAL(real_t) x(*), y(*)
  !-----------------------------------------------------------------------
  !     interface information exchange routine for a repeated
  !     communication pattern on a tensor product grid
  !     MSG version 2.0
  !
  !-----------------------------------------------------------------------
  !
  !
  ! arguments:
  !------------
  !
  !     x     = input array (a multidimensional prism), its boundaries
  !             will be updated with the data from the other processors
  !     y     = work array of the size of at least two times
  !             the maximal boundary segment
  !     nproc = number of adjacent processors (input)
  !     proc  = array of size nproc containing the numbers (IDs)
  !             of neighboring processors (input).
  !     ipr   = array containing pointers to the beginnings
  !             of each segment in the array of indices index (input)
  !     index = array of indices of boundary elements (input)
  !             a negative index of the first element of a segment
  !             indicates that the segment is contiguous and will
  !             be processed in place
  !     ptrn  = indicates a pattern to use;
  !             ptrn must be between 1 and MAX_PATTERNS,
  !             Only a limited number of patterns can be allocated.
  !             A pattern must be explicitly deallocated if it is not
  !             needed anymore. Another pattern can be later opened
  !             with the same ptrn.
  !
  ! return code:
  !-------------
  !
  !       ierr =     0 --- boundary information has been sent,
  !                        pattern stays open
  !                 -1 --- error opening the channels
  !                 -2 --- the number of pattern specified
  !                        is zero or larger than the maximum
  !                        allowed - increase MAX_PATTERNS and
  !                        recompile
  !                 -3 --- the specified number of adjacent processors is
  !                        wrong, either it is larger than the allowed
  !                        maximum (increase MAX_PROCS and recompile)
  !                        or is less than 0
  !                 >0 --- error in MPI functions (see MPI error codes)
  !
  !-----------------------------------------------------------------------
  !     local variables
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  integer iproc, InSegmentStart, InSegmentSize, myproc, &
  &        OutSegmentSize, OutSegmentStart
  ! external MSG_myproc
  ! logical first_call
  ! data first_call/.true./
  ! save first_call
  !
  ierr = 0

  if (MSG_BLOCKING .eq. 1) THEN  ! use blocking communication

     if(nproc.eq.0) return
     myproc = MSG_Myproc()

     if(proc(1).lt.0) then
        proc(1) = -proc(1)
     endif

     do iproc = 1, nproc
        OutSegmentSize = ipr(2*iproc) - ipr(2*iproc-1)
        OutSegmentStart = ipr(2*iproc-1)

        InSegmentSize = ipr(2*iproc+1) - ipr(2*iproc)
        InSegmentStart = ipr(2*iproc)


        if (myproc .lt. proc(iproc)) then

           if (index(OutSegmentStart).ge.0) then

              call MSG_tbdx_gather(x, y, iproc, ipr, index)
              call MPI_Send(y, OutSegmentSize, &
              &                 8, &
              &                 proc(iproc)-1,ptrn,MSG_COMM, &
              &                 ierr)
              if(ierr.ne.MPI_SUCCESS) return

           else

              call MPI_Send(x(-index(OutSegmentStart)), &
              &                 OutSegmentSize, &
              &                 8, &
              &                 proc(iproc)-1, ptrn, MSG_COMM, &
              &                 ierr )

           endif

           if(index(InSegmentStart).ge.0) then
              call MPI_Recv(y, InSegmentSize, &
              &                 8, &
              &                 proc(iproc)-1, ptrn, MSG_COMM, &
              &                 RecvStatus,ierr)
              if(ierr.ne.MPI_SUCCESS) return
              call MSG_tbdx_scatter(x, y, iproc, ipr, index)

           else

              call MPI_Recv(x(-index(InSegmentStart)), &
              &                 InSegmentSize, &
              &                 8, &
              &                 proc(iproc)-1, ptrn, MSG_COMM, &
              &                 RecvStatus, ierr)
              if(ierr.ne.MPI_SUCCESS) return

           endif

        else

           if(index(InSegmentStart).ge.0) then
              call MPI_Recv(y, InSegmentSize, &
              &                 8, &
              &                 proc(iproc)-1, ptrn, MSG_COMM, &
              &                 RecvStatus,ierr)
              if(ierr.ne.MPI_SUCCESS) return

              call MSG_tbdx_scatter(x, y, iproc, ipr, index)

           else

              call MPI_Recv(x(-index(InSegmentStart)), &
              &                 InSegmentSize, &
              &                 8, &
              &                 proc(iproc)-1, ptrn, MSG_COMM, &
              &                 RecvStatus, ierr)
              if(ierr.ne.MPI_SUCCESS) return

           endif

           if (index(OutSegmentStart).ge.0) then

              call MSG_tbdx_gather(x, y, iproc, ipr, index)
              call MPI_Send(y, OutSegmentSize, &
              &                 8, &
              &                 proc(iproc)-1,ptrn,MSG_COMM, &
              &                 ierr)
              if(ierr.ne.MPI_SUCCESS) return

           else

              call MPI_Send(x(-index(OutSegmentStart)), &
              &                 OutSegmentSize, &
              &                 8, &
              &                 proc(iproc)-1, ptrn, MSG_COMM, &
              &                 ierr )

           endif

        end if
     END DO

  ELSE ! use non-blocking communication

     if(nproc.eq.0) return
     if( ptrn.gt.MAX_PATTERNS.or.ptrn.le.0 ) then
        ierr = -2
        return
     endif
     if( nproc.gt.MAX_PROCS.or.nproc.lt.0) then
        ierr = -3
        return
     endif
     ! if(first_call) then
        do iproc=1,MAX_PATTERNS
           MSG_sendid(1,iproc)=0
           MSG_recvid(1,iproc)=0
        enddo
        ! first_call = .false.
     !endif
     if(MSG_sendid(1,ptrn).eq.0.and.MSG_recvid(1,ptrn).eq.0) then
        !
        !     open the communication channels
        !
        !     set the type of data transfer for this pattern
        !
        if(proc(1).lt.0) then
           MSG_TRANSFER_TYPE(ptrn) = 1
           proc(1) = -proc(1)
        else
           MSG_TRANSFER_TYPE(ptrn) = 0
        endif
        !
        !     find the maximal size of outgoing segment in order to
        !     find a "safe" place within the array y to put the
        !     buffer for the incoming data
        !
        MSGSegment(ptrn) = 0
        do iproc = 1, nproc
           OutSegmentSize = ipr(2*iproc) - ipr(2*iproc-1)
           if(MSGSegment(ptrn).lt.OutSegmentSize) &
                &              MSGSegment(ptrn) = OutSegmentSize
        enddo
        MSGSegment(ptrn) = MSGSegment(ptrn) + 1
        !
        !     open up channels
        !
        do iproc = 1, nproc
           OutSegmentSize = ipr(2*iproc) - ipr(2*iproc-1)
           OutSegmentStart = ipr(2*iproc-1)
           if(OutSegmentSize.gt.0) then
              if(index(OutSegmentStart).ge.0) then
                 !     noncontiguous memory segment: give the buffer's address
                 call MPI_send_init(y, OutSegmentSize, &
                      &                    REAL(real_t)_PRECISION, &
                      &                    proc(iproc)-1, ptrn, MSG_COMM, &
                      &                    MSG_sendid(iproc,ptrn), ierr)
                 if(ierr.ne.MPI_SUCCESS) return
              else
                 !     contiguous memory segment: give the data address
                 call MPI_send_init(x(-index(OutSegmentStart)), &
                      &                    OutSegmentSize, &
                      &                    REAL(real_t)_PRECISION, &
                      &                    proc(iproc)-1, ptrn, MSG_COMM, &
                      &                    MSG_sendid(iproc,ptrn), ierr)
                 if(ierr.ne.MPI_SUCCESS) return
              endif
           endif
           InSegmentSize = ipr(2*iproc+1) - ipr(2*iproc)
           InSegmentStart = ipr(2*iproc)
           if(InSegmentSize.gt.0) then
              if(index(InSegmentStart).ge.0) then
                 !     noncontiguous memory segment: give the buffer's address
                 call MPI_recv_init(y(MSGSegment(ptrn)), &
                      &                    InSegmentSize, &
                      &                    REAL(real_t)_PRECISION, &
                      &                    proc(iproc)-1, ptrn, MSG_COMM, &
                      &                    MSG_recvid(iproc,ptrn), ierr)
                 if(ierr.ne.MPI_SUCCESS) return
              else
                 !     contiguous memory segment: give the data address
                 call MPI_recv_init(x(-index(InSegmentStart)), &
                      &                    InSegmentSize, &
                      &                    REAL(real_t)_PRECISION, &
                      &                    proc(iproc)-1, ptrn, MSG_COMM, &
                      &                    MSG_recvid(iproc,ptrn), ierr)
                 if(ierr.ne.MPI_SUCCESS) return
              endif
           endif
        enddo
     endif
     !
     if(MSG_TRANSFER_TYPE(ptrn).eq.1) then
        !
        !     exchange data through the channels using all to all
        !
        !     send all the messages out one by one
        do iproc = 1, nproc
           OutSegmentSize = ipr(2*iproc) - ipr(2*iproc-1)
           if(OutSegmentSize.gt.0) then
              !     gather the outgoing data in the outgoing buffer
              call MSG_tbdx_gather(x, y, iproc, ipr, index)
              !     start sending the outgoing data to iproc
              call MPI_start(MSG_sendid(iproc,ptrn), ierr)
              if(ierr.ne.MPI_SUCCESS) return
              call MPI_wait(MSG_sendid(iproc,ptrn),SendStatus,ierr)
              if(ierr.ne.MPI_SUCCESS) return
           endif
        enddo
     else
        !
        !     exchange data through the channels using shifts
        !
        if(nproc.gt.1) then
           do iproc = 1, nproc - 1
              OutSegmentSize = ipr(2*iproc) - ipr(2*iproc-1)
              InSegmentSize = ipr(2*iproc+1) - ipr(2*iproc)
              !     start receiving the incoming data from iproc
              if(InSegmentSize.gt.0) then
                 call MPI_start(MSG_recvid(iproc,ptrn), ierr)
                 if(ierr.ne.MPI_SUCCESS) return
              endif
              if(OutSegmentSize.gt.0) then
                 !     gather the outgoing data in the outgoing buffer
                 call MSG_tbdx_gather(x, y, iproc, ipr, index)
                 !     start sending the outgoing data to iproc
                 call MPI_start(MSG_sendid(iproc,ptrn), ierr)
                 if(ierr.ne.MPI_SUCCESS) return
              endif
              if(InSegmentSize.gt.0) then
                 !     wait for the incoming data to be received
                 call MPI_wait(MSG_recvid(iproc,ptrn), &
                 &                    RecvStatus, ierr)
                 if(ierr.ne.MPI_SUCCESS) return
                 !     place the incoming data into the local array
                 call MSG_tbdx_scatter(x, y(MSGSegment(ptrn)), &
                 &                    iproc, ipr, index)
              endif
              !     wait for the outgoing data to be sent
              if(OutSegmentSize.gt.0) then
                 call MPI_wait(MSG_sendid(iproc,ptrn), &
                 &                    SendStatus, ierr)
                 if(ierr.ne.MPI_SUCCESS) return
              endif
           enddo
        endif
        iproc = nproc
        OutSegmentSize = ipr(2*iproc) - ipr(2*iproc-1)
        InSegmentSize = ipr(2*iproc+1) - ipr(2*iproc)
        !     start receiving the incoming data from nproc
        if(InSegmentSize.gt.0) then
           call MPI_start(MSG_recvid(iproc,ptrn), ierr)
           if(ierr.ne.MPI_SUCCESS) return
        endif
        if(OutSegmentSize.gt.0) then
           !     gather the outgoing data in the outgoing buffer
           call MSG_tbdx_gather(x, y, iproc, ipr, index)
           !     start sending the outgoing data to nproc
           call MPI_start(MSG_sendid(iproc,ptrn), ierr)
           if(ierr.ne.MPI_SUCCESS) return
        endif
        !
     endif

  ENDIF

  return
end subroutine MSG_tbdx_send



subroutine MSG_tbdx_close(x,y,nproc,proc,ipr,index,ptrn,ierr)
  implicit none
  integer nproc,ptrn,ierr
  integer proc(nproc),ipr(*),index(*)
  REAL(real_t) x(*),y(*)
  !-----------------------------------------------------------------------
  !
  !   deallocates a communication pattern (closes a channel)
  !   for a tensor product grid
  !   for the details, see MSG_tdbx_send
  !   MPI version 2.0
  !
  !-----------------------------------------------------------------------
  !
  !     local variables
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  integer ii




  ierr = 0
  IF ( MSG_BLOCKING .eq. 0 ) THEN

     if(nproc.eq.0) return

     do ii=1,nproc
        if(MSG_sendid(ii,ptrn).ne.0) then
           call MPI_request_free(MSG_sendid(ii,ptrn), ierr)
           if(ierr.ne.MPI_SUCCESS) return
        endif
        if(MSG_recvid(ii,ptrn).ne.0) then
           call MPI_request_free(MSG_recvid(ii,ptrn), ierr)
           if(ierr.ne.MPI_SUCCESS) return
        endif
     enddo

     MSG_sendid(1,ptrn) = 0
     MSG_recvid(1,ptrn) = 0
     ierr = 0

  ENDIF

  return
end subroutine MSG_tbdx_close


subroutine MSG_tbdx_receive(x,y,nproc,proc,ipr,index,ptrn,ierr)
  implicit none
  integer nproc,ptrn,ierr
  integer proc(nproc), ipr(*), index(*)
  REAL(real_t) x(*), y(*)
  !-----------------------------------------------------------------------
  !
  !     interface information exchange routine for a repeated
  !     communication pattern and a tensor product grid
  !     for the details, see MSG_tdbx_send
  !     MSG version 2.0
  !
  ! arguments:
  !------------
  !
  !     x     = output vector (a multidimensional prism), its boundaries
  !             will be updated with the data from the other processors
  !     y     = input array containing the boundary information
  !     nproc = number of adjacent processors (input)
  !     proc  = array of size nproc containing the numbers (IDs)
  !             of neighboring processors (input).
  !     ipr   = array containing pointers to the beginnings
  !             of each segment in the buffer y (input)
  !     index = array of indices of boundary elements(input)
  !     ptrn  = indicates a pattern to use;
  !             ptrn must be between 1 and MAX_PATTERNS,
  !             Only a limited number of patterns can be allocated.
  !             A pattern must be explicitly deallocated if it is not
  !             needed anymore. Another pattern can be later opened
  !             with the same ptrn.
  !     ierr  = error code produced by the MPI routines (0 if success)
  !             (output)
  !
  !-----------------------------------------------------------------------
  !
  !     local variables
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  integer OutSegmentSize, InSegmentSize, iproc
  !


  ierr = 0

  IF ( MSG_BLOCKING .eq. 0 ) THEN

     if(nproc.eq.0) return
     if(MSG_TRANSFER_TYPE(ptrn).eq.1) then
        !
        !     use "all to all" mode
        !
        do iproc = 1, nproc
           InSegmentSize = ipr(2*iproc+1) - ipr(2*iproc)
           if(InSegmentSize.gt.0) then

              call MPI_start(MSG_recvid(iproc,ptrn), ierr)
              if(ierr.ne.MPI_SUCCESS) return
              call MPI_wait(MSG_recvid(iproc,ptrn),RecvStatus,ierr)
              if(ierr.ne.MPI_SUCCESS) return

              !     scatter the incoming data
              call MSG_tbdx_scatter(x, y(MSGSegment(ptrn)), &
              &                 iproc, ipr, index)
           endif
        enddo
     else
        !
        !     use the "series of shifts" mode
        !
        OutSegmentSize = ipr(2*nproc) - ipr(2*nproc-1)
        InSegmentSize = ipr(2*nproc+1) - ipr(2*nproc)
        !
        !     wait until the exchange with the last processor on the
        !     neighbors list is completed and
        !     reset the channels for the next transaction
        !
        if(InSegmentSize.gt.0) then

           call MPI_wait(MSG_recvid(nproc,ptrn), RecvStatus, ierr)
           if(ierr.ne.MPI_SUCCESS) return

           call MSG_tbdx_scatter(x, y(MSGSegment(ptrn)), &
                &                nproc, ipr, index)
        endif

        if(OutSegmentSize.gt.0) then
           call MPI_wait(MSG_sendid(nproc,ptrn), SendStatus, ierr)
           if(ierr.ne.MPI_SUCCESS) return
        endif
     endif

  ENDIF

  return
end subroutine MSG_tbdx_receive


subroutine MSG_tbdx_gather(x, y, iproc, ipr, index)
  implicit none
  integer iproc, ipr(*), index(*)
  REAL(real_t) x(*), y(*)
  !-----------------------------------------------------------------------
  !
  !     gathers the outgoing boundary information into a single buffer
  !
  !-----------------------------------------------------------------------
  !
  !
  ! arguments:
  !------------
  !
  !     x     = input array (local data - a multidimensional prism)
  !     y     = buffer to be filled
  !     iproc = index of the adjacent procesor within the neighbors list
  !             (input)
  !     ipr   = array containing pointers to the beginnings
  !             of each segment in the buffer y (input)
  !     index = array of indices of boundary elements(input)
  !             if the first index of the segment is negative
  !             then the data is processed in place
  !
  !-----------------------------------------------------------------------
  !
  !     local variables
  integer OutSegmentStart, OutSegmentEnd, OutSegmentSize, j, k
  OutSegmentStart = ipr(2*iproc-1)
  if(index(OutSegmentStart).lt.0) return
  OutSegmentSize = ipr(2*iproc) - ipr(2*iproc-1)
  OutSegmentEnd = OutSegmentStart + OutSegmentSize - 1
  k=1
  do j = OutSegmentStart, OutSegmentEnd
     y(k)=x(index(j))
     k = k + 1
  enddo
  return
end subroutine MSG_tbdx_gather


subroutine MSG_tbdx_scatter(x, y, iproc, ipr, index)
  implicit none
  integer iproc, ipr(*), index(*)
  REAL(real_t) x(*), y(*)
  !-----------------------------------------------------------------------
  !
  !     places boundary information into the local array
  !
  !
  ! arguments:
  !------------
  !
  !     x     = output array (local data - a multidimensional prism)
  !     y     = address of the buffer for incoming data
  !     iproc = index of the adjacent procesor within the neighbors list
  !             (input)
  !     ipr   = array containing pointers to the beginnings
  !             of each segment in the buffer y
  !     index = array of indices of boundary elements(input)
  !             if the first index of a segment is negative than
  !             the segment is processed in place
  !
  !-----------------------------------------------------------------------
  !
  !     local variables
  integer InSegmentSize, InSegmentStart, InSegmentEnd, j, k
  InSegmentStart = ipr(2*iproc)
  if(index(InSegmentStart).lt.0) return
  InSegmentSize = ipr(2*iproc+1) - ipr(2*iproc)
  InSegmentEnd = InSegmentStart + InSegmentSize - 1
  k=1
  do j = InSegmentStart, InSegmentEnd
     x(index(j))=y(k)
     k = k + 1
  enddo
  return
end subroutine MSG_tbdx_scatter


logical function MSG_shift_in_trouble(nproc, ipr, gc_eid, gc_ld, &
     &                                FirstGlobalIndex, LastGlobalIndex, MyPeriodic)
  implicit none
  integer nproc, ipr(*), gc_eid(6), gc_ld(6), &
       & FirstGlobalIndex(3), LastGlobalIndex(3), MyPeriodic(3)
  ! local data
  !
  include 'geom_param_fort_f90.h'
  include 'mpi_param_fort_f90.h'
  integer MyEID(6), MyLDSize, MyEIDSize, MyHaloSize, ierr, &
       &        MyInSize, iproc, LocType, GlobalType
  !
  MSG_shift_in_trouble = .false.
  LocType = 0
  !
  ! find the size of the halo=eid-ld segment
  !
  MyLDSize = (gc_ld(2) - gc_ld(1) + 1)* &
       &           (gc_ld(4) - gc_ld(3) + 1)* &
       &           (gc_ld(6) - gc_ld(5) + 1)
  call MSG_to_contiguous(MyEID, gc_eid, &
       &     FirstGlobalIndex, LastGlobalIndex, MyPeriodic)
  MyEIDSize = (MyEID(2) - MyEID(1) + 1)* &
       &            (MyEID(4) - MyEID(3) + 1)* &
       &            (MyEID(6) - MyEID(5) + 1)
  MyHaloSize = MyEIDSize - MyLDSize
  !
  ! find the size of my incoming segment
  !
  MyInSize = 0
  do iproc=1,nproc
     MyInSize = MyInSize + ipr(2*iproc+1) - ipr(2*iproc)
  enddo
  if(MyHaloSize.gt.MyInSize) LocType = 1
  call MPI_Allreduce(LocType, GlobalType, 1, MPI_INTEGER, MPI_MAX, MSG_COMM, ierr)
  MSG_shift_in_trouble = GlobalType.eq.1
  return
end function MSG_shift_in_trouble

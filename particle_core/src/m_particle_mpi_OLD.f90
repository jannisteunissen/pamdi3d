module m_particle_mpi
  use m_particle_core

  implicit none
  public

  integer :: PM_myrank, PM_ntasks, PM_comm

contains

    !> MPI: Share the particles evenly over the tasks, with no domain decomposition
  subroutine shareParticles(myrank, ntasks, commm)
    use mpi
    use m_utils
    integer, intent(IN) :: myrank, ntasks, comm

    integer :: iHasMore, iHasLess, nPartSend, nRecvs, nSends
    integer :: ierr, meanNP, iMin, n
    integer :: sendReqs(ntasks), recvReqs(ntasks)
    integer :: nPartPerTask(0:ntasks-1)

    ! Get the number of particles each task has
    call MPI_ALLGATHER(PM_nPart, 1, MPI_integer, nPartPerTask, 1, MPI_integer, comm, ierr)
    
    meanNP = (sum(nPartPerTask) + ntasks - 1) / ntasks
    nSends = 0
    nRecvs = 0
    iHasLess = 0
    iHasMore = 0

    do while (iHasMore < ntasks)

       if (nPartPerTask(iHasMore) > meanNP) then
          ! Find first task with less than mean particles
          iHasLess = UT_ixFirstTrue(nPartPerTask < meanNP)

          ! Determine amount of particles to send
          nPartSend = min(nPartPerTask(iHasMore) - meanNP, meanNP - nPartPerTask(iHasLess))

          if (myrank == iHasMore) then ! Send particles
             iMin     = nPartPerTask(iHasMore) - nPartSend + 1
             nSends   = nSends + 1
             call MPI_ISEND( PM_pList(iMin), nPartSend, PM_pTypeMPI, iHasLess, iHasMore, &
                  & comm, sendReqs(nSends), ierr)
          else if (myrank == iHasLess) then ! Receive particles
             iMin     = nPartPerTask(iHasLess) + 1
             call checkNumParticles(iMin + nPartSend - 1)
             nRecvs   = nRecvs + 1
             call MPI_IRECV( PM_pList(iMin), nPartSend, PM_pTypeMPI, iHasMore, iHasMore, &
                  comm, recvReqs(nRecvs), ierr)
          end if

          ! Update the number of particles for each task
          nPartPerTask(iHasLess) = nPartPerTask(iHasLess) + nPartSend
          nPartPerTask(iHasMore) = nPartPerTask(iHasMore) - nPartSend
       end if

       if (nPartPerTask(iHasMore) <= meanNP) iHasMore = iHasMore + 1
    end do

    if (nRecvs > 0) call MPI_WAITALL(nRecvs, recvReqs, MPI_STATUSES_IGNORE, ierr)
    call MPI_BARRIER(comm, ierr)
    PM_nPart = nPartPerTask(myrank)

  end subroutine shareParticles


  !> Divide the particles over the tasks based on their cell index
  subroutine PM_divideParticles(myrank, ntasks)
    use mpi
    integer, intent(in) :: myrank, ntasks

    integer :: ll, ix, rank, ierr, iMin, iMax, nPartBehind, nPartInInterval, nPartSend
    integer :: tag, partCnt
    integer :: nSends, nRecvs, sender, recver, meanNP, nptemp, nSteps
    integer :: sendReqs(ntasks), recvReqs(ntasks)
    integer :: nPartPerTask(0:ntasks-1), nPartPerIntervalPerTask(0:ntasks-1, 0:ntasks-1)
    integer :: splitIx(-1:ntasks), minPos, maxPos, middlePos, maxCellIx
    integer, allocatable :: nPartPerCellIx(:)
    real(dp), allocatable :: cellIxs(:)

    !       print *, "Before divide ", myrank, " has ", PM_nPart, " particles"

    ! Get the number of particles each task has
    call MPI_ALLGATHER(PM_nPart, 1, MPI_integer, nPartPerTask, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    meanNP      = (sum(nPartPerTask) + ntasks - 1) / ntasks
    maxCellIx   = EFA_getMaxCellIx()

    allocate( nPartPerCellIx(maxCellIx) )
    allocate( cellIxs(PM_nPart) )

    nPartPerCellIx = 0

    ! Set the cell indexes for the particles
    do ll = 1, PM_nPart
       ix                   = EFA_getCellIndexAt(PM_pList(ll)%x)
       cellIxs(ll)          = ix
       nPartPerCellIx(ix)   = nPartPerCellIx(ix) + 1
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, nPartPerCellIx, maxCellIx, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    call heapsortParticleByDbles(PM_pList(1:PM_nPart),cellIxs)

    ! Find the splitIx(:) s.t. task n will hold particles with ix > splitIx(n-1) and ix <= splitIx(n)
    splitIx(-1)                   = 0
    splitIx(ntasks)               = maxCellIx
    rank                          = 0
    nPartBehind                   = 0
    nPartPerIntervalPerTask(:,:)  = 0
    nPartInInterval               = 0

    do ix = 1, maxCellIx
       nPartInInterval = nPartInInterval + nPartPerCellIx(ix)

       if (nPartInInterval >= meanNP .or. ix == maxCellIx) then
          splitIx(rank)     = ix
          nPartInInterval   = 0
          partCnt           = 0

          do ll = nPartBehind+1, PM_nPart
             !                print *, ll, cellIxs(ll)
             if (cellIxs(ll) <= ix) then
                partCnt = partCnt + 1
             else
                exit
             end if
          end do

          nPartPerIntervalPerTask(rank, myrank) = partCnt
          nPartBehind       = nPartBehind + partCnt
          rank              = rank + 1

          if (rank == ntasks) exit
       end if

    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, nPartPerIntervalPerTask, ntasks*ntasks, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Send particles to task that holds a region
    nSends = 0
    nRecvs = 0

    do recver = 0, ntasks - 1

       ! Receive the particles in tasks' region from all other tasks
       if (myrank == recver) then

          do sender = 0, ntasks - 1
             if (sender == myrank) cycle ! Don't have to send to ourselves

             nPartSend = nPartPerIntervalPerTask(recver, sender)

             if (nPartSend > 0) then
                nRecvs = nRecvs + 1
                iMin   = PM_nPart + 1
                PM_nPart = PM_nPart + nPartSend
                call checkNumParticles(PM_nPart)
                tag    = sender

                !                   print *, myrank, ": receives", nPartSend, "from", sender, "iMin", iMin
                call MPI_IRECV( PM_pList(iMin), nPartSend, PM_pTypeMPI, sender, tag, &
                     & MPI_COMM_WORLD, recvReqs(nRecvs), ierr)
             end if
          end do


       else ! Send particles to recver

          nPartSend = nPartPerIntervalPerTask(recver, myrank)

          if (nPartSend > 0) then
             nSends      = nSends + 1
             tag         = myrank

             ! Find index of first particle to send
             iMin        = sum( nPartPerIntervalPerTask(0:recver-1, myrank) ) + 1

             !                print *, myrank, ": sends", nPartSend, "to", recver, "iMin", iMin
             call MPI_ISEND( PM_pList(iMin), nPartSend, PM_pTypeMPI, recver, tag, &
                  & MPI_COMM_WORLD, sendReqs(nSends), ierr)

             ! Important, wait until the particles are send over before deactivating them
             call MPI_WAIT(sendReqs(nSends), MPI_STATUS_IGNORE, ierr)

             ! Mark the particles that we have send away as inactive
             do ll = iMin, iMin + nPartSend - 1
                call PM_removePart(ll)
             end do

          end if
       end if

    end do

    if (nRecvs > 0) call MPI_WAITALL(nRecvs, recvReqs, MPI_STATUSES_IGNORE, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !       print *, "After divide ", myrank, " has ", PM_nPart, "PM_nPart (some dead)"
    call removeDeadParticles()
    !       print *, "After divide ", myrank, " has ", PM_nPart, " particles"

    deallocate( nPartPerCellIx )
    deallocate( cellIxs )

  end subroutine PM_divideParticles

  subroutine createPartTypeMPI(PM_pTypeMPI)
    integer, intent(INOUT) :: PM_pTypeMPI
    integer, parameter :: nf = 3
    integer :: blocklen(nf), types(nf), ierr, type_size, mpi_type_size
    integer(KIND=MPI_ADDRESS_KIND) :: displ(nf), lb, extent
    type(PM_pType) :: electron_example

    blocklen = (/11, 1, 1/) ! 11 doubles, 1 integer, 1 logical

    displ(1) = 0
    call MPI_TYPE_GET_EXTENT(MPI_DOUBLE_PRECISION, lb, extent, ierr)
    displ(2) = displ(1) + blocklen(1) *  extent
    call MPI_TYPE_GET_EXTENT(MPI_integer, lb, extent, ierr)
    displ(3) = displ(2) + blocklen(2) * extent

    types = (/MPI_DOUBLE_PRECISION, MPI_integer, MPI_logical/)

    call MPI_TYPE_CREATE_STRUCT(nf, blocklen, displ, types, PM_pTypeMPI, ierr)
    call MPI_TYPE_COMMIT(PM_pTypeMPI, ierr)

    type_size = storage_size(electron_example)
    mpi_type_size = blocklen(1) * storage_size(0.0d0) + blocklen(2) * storage_size(1) + blocklen(3) * storage_size(.true.)

    if (type_size /= mpi_type_size .or. ierr /= 0) call ERR_show("createPartTypeMPI size error!")
  end subroutine createPartTypeMPI

  
  !> Return the number of real particles
  real(dp) function PM_getNumRealPart()
    use mpi
    integer  :: ierr
    real(dp) :: part_sum
    part_sum = sum(PM_pList(1:PM_nPart)%weight)
    call MPI_ALLREDUCE(part_sum, PM_getNumRealPart, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
  end function PM_getNumRealPart

  !> Return the number of simulation particles
  integer function PM_getNumSimPart()
    use mpi
    integer :: ierr, 
    call MPI_ALLREDUCE(PM_nPart, PM_getNumSimPart, 1, MPI_integer, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
  end function PM_getNumSimPart


end module m_particle_mpi

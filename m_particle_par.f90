module m_particle_par
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: PP_get_num_real_part
  public :: PP_get_num_sim_part
  public :: PP_share_mpi
  public :: PP_reorder_by_bins_mpi

contains

  integer function PP_get_num_sim_part(pc)
    use mpi
    use m_particle_core
    type(PC_t), intent(in) :: pc
    integer                :: ierr, n_part_sum

    n_part_sum = pc%get_num_sim_part()
    call MPI_ALLREDUCE(n_part_sum, PP_get_num_sim_part, 1, MPI_integer, &
         MPI_SUM, MPI_COMM_WORLD, ierr)
  end function PP_get_num_sim_part

  real(dp) function PP_get_num_real_part(pc)
    use mpi
    use m_particle_core
    type(PC_t), intent(in) :: pc
    integer                :: ierr
    real(dp)               :: n_part_sum

    n_part_sum = pc%get_num_real_part()
    call MPI_ALLREDUCE(n_part_sum, PP_get_num_real_part, 1, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
  end function PP_get_num_real_part

  subroutine PP_share_mpi(pc, myrank, ntasks)
    use mpi
    type(PC_t), intent(inout) :: pc
    integer, intent(IN) :: myrank, ntasks

    integer :: i_more, i_less, n_send
    integer :: ierr, mean_n_part, i_min, i_max, tag
    integer, dimension(0:ntasks-1) :: n_part_task

    ! print *, myrank, "Before sharing: pc%n_part", pc%n_part
    ! Get the number of particles each task has
    call MPI_ALLGATHER(pc%n_part, 1, MPI_integer, n_part_task, 1, &
         MPI_integer, MPI_COMM_WORLD, ierr)
    mean_n_part = (sum(n_part_task) + ntasks - 1) / ntasks

    i_less = 0
    i_more = 0

    do while (i_more < ntasks)
       if (n_part_task(i_more) > mean_n_part) then

          ! Find first task with less than mean particles
          do while (n_part_task(i_less) >= mean_n_part)
             i_less = i_less + 1
          end do

          ! Determine amount of particles to send
          n_send = min(n_part_task(i_more) - mean_n_part, &
               mean_n_part - n_part_task(i_less))

          if (myrank == i_more) then ! Send
             i_min = n_part_task(i_more) - n_send + 1
             i_max = i_min + n_send - 1
             tag   = i_more
             call send_parts_mpi(pc%particles(i_min:i_max), i_less, tag)
          else if (myrank == i_less) then ! Receive
             i_min = n_part_task(i_less) + 1
             i_max = i_min + n_send - 1
             call pc%check_space(i_max)
             tag   = i_more
             call recv_parts_mpi(pc%particles(i_min:i_max), i_more, tag)
          end if

          ! Update the number of particles for each task
          n_part_task(i_less) = n_part_task(i_less) + n_send
          n_part_task(i_more) = n_part_task(i_more) - n_send
       end if

       if (n_part_task(i_more) <= mean_n_part) i_more = i_more + 1
    end do

    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
    pc%n_part = n_part_task(myrank)
    ! print *, myrank, "After sharing: pc%n_part", pc%n_part
  end subroutine PP_share_mpi

  ! Divide the particles over the tasks based on their bin index
  subroutine PP_reorder_by_bins_mpi(pc, binner, myrank, ntasks)
    use mpi
    use m_lookup_table
    use m_mrgrnk
    type(PC_t), intent(inout)   :: pc
    class(PC_bin_t), intent(in) :: binner
    integer, intent(in)         :: myrank, ntasks
    integer                     :: ll, ix, dst, ierr, n_behind
    integer                     :: n_part_interval, n_send
    integer                     :: i_min, i_max
    integer                     :: tag, part_count
    integer                     :: n_mean_left, n_part_left
    integer                     :: n_sends, n_recvs, sender, recver
    integer                     :: n_part_task(0:ntasks-1)
    integer                     :: n_part_dst_src(0:ntasks-1, 0:ntasks-1)
    integer, allocatable        :: n_part_per_bin(:), ix_list(:)
    real(dp), allocatable       :: bin_ixs(:)

    print *, "Before divide ", myrank, " has ", pc%n_part, " particles"

    ! Get the number of particles each task has
    call MPI_ALLGATHER(pc%n_part, 1, MPI_integer, n_part_task, 1, &
         MPI_integer, MPI_COMM_WORLD, ierr)

    allocate(n_part_per_bin(binner%n_bins))
    allocate(bin_ixs(pc%n_part))
    allocate(ix_list(pc%n_part))
    n_part_per_bin = 0

    ! Set the bin indexes for the particles
    do ll = 1, pc%n_part
       ix                   = binner%bin_func(pc%particles(ll))
       bin_ixs(ll)          = ix
       n_part_per_bin(ix)   = n_part_per_bin(ix) + 1
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, n_part_per_bin, binner%n_bins, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Sort by bin_ixs
    call mrgrnk(bin_ixs, ix_list)
    bin_ixs                   = bin_ixs(ix_list)
    pc%particles(1:pc%n_part) = pc%particles(ix_list)

    dst                 = 0
    n_behind            = 0
    n_part_dst_src(:,:) = 0
    n_part_interval     = 0
    n_part_left         = sum(n_part_task)
    n_mean_left         = ceiling(sum(n_part_task) / real(ntasks, dp))

    do ix = 1, binner%n_bins
       n_part_interval = n_part_interval + n_part_per_bin(ix)

       if (n_part_interval >= n_mean_left .or. ix == binner%n_bins) then
          part_count = 0
          do ll = n_behind+1, pc%n_part
             if (bin_ixs(ll) > ix) exit
             part_count = part_count + 1
          end do

          n_part_dst_src(dst, myrank) = part_count
          n_behind                    = n_behind + part_count
          dst                         = dst + 1
          n_part_left                 = n_part_left - n_part_interval
          n_part_interval             = 0

          if (dst == ntasks) exit ! dst runs from from 0 to ntasks-1
          n_mean_left = ceiling(n_part_left / real(ntasks-dst, dp))
       end if

    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, n_part_dst_src, ntasks*ntasks, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Send particles to task that holds a region
    n_sends = 0
    n_recvs = 0

    do recver = 0, ntasks - 1
       if (myrank == recver) then ! Receive
          do sender = 0, ntasks - 1
             if (sender == myrank) cycle ! Skip ourselves
             n_send = n_part_dst_src(recver, sender)

             if (n_send > 0) then
                n_recvs   = n_recvs + 1
                tag       = sender
                i_min     = pc%n_part + 1
                i_max     = pc%n_part + n_send
                pc%n_part = i_max
                call pc%check_space(i_max)
                ! print *, myrank, ": receives", n_send, "from", sender, "i_max", i_max
                call recv_parts_mpi(pc%particles(i_min:i_max), sender, tag)
             end if
          end do
       else ! Send
          n_send = n_part_dst_src(recver, myrank)

          if (n_send > 0) then
             n_sends = n_sends + 1
             tag     = myrank

             ! Find index of first particle to send
             i_min    = sum(n_part_dst_src(0:recver-1, myrank)) + 1
             i_max = i_min + n_send - 1

             ! print *, myrank, ": sends", n_send, "to", recver, "i_min", i_min
             call send_parts_mpi(pc%particles(i_min:i_max), recver, tag)

             ! Mark the particles that we have send away as inactive
             do ll = i_min, i_max
                call pc%remove_part(ll)
             end do

          end if
       end if
    end do

    ! print *, "After divide ", myrank, " has ", pc%n_part, "pc%n_part (some dead)"
    call pc%clean_up()
    print *, "After divide ", myrank, " has ", pc%n_part, " particles"
  end subroutine PP_reorder_by_bins_mpi

  subroutine send_parts_mpi(parts, recver, tag)
    use mpi
    type(PC_part_t), intent(in) :: parts(:)
    integer, intent(in) :: recver, tag

    integer :: n_send, req_size, rel_size, ierr
    real(dp), allocatable :: r_buf(:)

    ! We use a buffer because we cannot send the part type *neatly* with mpi
    rel_size = storage_size(parts(1)) / storage_size(1.0_dp)
    n_send   = size(parts)
    req_size = n_send * rel_size
    allocate(r_buf(req_size))

    r_buf = transfer(parts, r_buf, req_size)
    call MPI_SEND(r_buf, req_size, MPI_DOUBLE_PRECISION, recver, tag, &
         & MPI_COMM_WORLD, ierr)
  end subroutine send_parts_mpi

  subroutine recv_parts_mpi(parts, sender, tag)
    use mpi
    type(PC_part_t), intent(inout) :: parts(:)
    integer, intent(in) :: sender, tag

    integer :: n_send, req_size, rel_size, ierr
    integer :: mpi_status(MPI_STATUS_SIZE)
    real(dp), allocatable :: r_buf(:)

    ! We use a buffer because we cannot send the part type *neatly* with mpi
    n_send   = size(parts)
    rel_size = storage_size(parts(1)) / storage_size(1.0_dp)
    req_size = n_send * rel_size
    allocate(r_buf(req_size))
    call MPI_RECV(r_buf, req_size, MPI_DOUBLE_PRECISION, sender, tag, &
         & MPI_COMM_WORLD, mpi_status, ierr)
    parts = transfer(r_buf, parts, n_send)
  end subroutine recv_parts_mpi

end module m_particle_par
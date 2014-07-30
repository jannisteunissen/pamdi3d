module m_particle
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: PM_particles_to_density
  public :: PM_fld_error
  public :: PM_get_num_part_mpi
  public :: PM_divide_particles
  public :: PM_get_max_dt

contains

  subroutine PM_fld_error(pc, rng, n_samples, fld_err, only_store)
    use m_efield_amr
    use m_random
    type(PC_t), intent(in) :: pc
    type(RNG_t), intent(inout) :: rng
    integer, intent(in) :: n_samples
    real(dp), intent(out) :: fld_err
    logical, intent(in) :: only_store

    integer :: i, n
    real(dp) :: diff, total, fld(3)
    real(dp), allocatable, save :: pos_samples(:, :)
    real(dp), allocatable, save :: fld_samples(:, :)

    if (pc%n_part < 1) then
       print *, "PM_fld_error: not enough particles"
       stop
    end if

    if (only_store) then
       fld_err = 0
       if (allocated(fld_samples)) deallocate(fld_samples)
       if (allocated(pos_samples)) deallocate(pos_samples)
       allocate(fld_samples(3, n_samples))
       allocate(pos_samples(3, n_samples))

       do i = 1, n_samples
          ! Randomly select a particle
          n = int(rng%uni_01() * pc%n_part) + 1
          pos_samples(:,i) = pc%particles(n)%x
          fld_samples(:,i) = E_get_field(pos_samples(:,i))
       end do
    else
       total = 0
       diff = 0

       do i = 1, n_samples
          fld = E_get_field(pos_samples(:,i))
          total = total + norm2(fld_samples(:,i) )
          diff  = diff + norm2(fld - fld_samples(:,i) )
       end do

       fld_err = diff / (total + epsilon(1.0_dp))
    end if
  end subroutine PM_fld_error

  subroutine PM_particles_to_density(pc)
    use m_efield_amr
    type(PC_t), intent(in) :: pc
    integer                :: ll

    call E_set_vars((/E_i_elec/), (/0.0_dp/))

    do ll = 1, pc%n_part
       call E_add_to_var(E_i_elec, pc%particles(ll)%x, pc%particles(ll)%w)
    end do
  end subroutine PM_particles_to_density

  integer function PM_get_num_part_mpi(pc)
    use mpi
    use m_particle_core
    type(PC_t), intent(in) :: pc
    integer :: ierr, n_part_sum

    n_part_sum = pc%n_part
    call MPI_ALLREDUCE(n_part_sum, PM_get_num_part_mpi, 1, MPI_integer, MPI_SUM, &
         MPI_COMM_WORLD, ierr)
  end function PM_get_num_part_mpi

  function PM_get_max_dt(pc, myrank, root) result(dt_max)
    type(PC_t), intent(in) :: pc
    integer, intent(in) :: myrank, root
    real(dp) :: dt_max
    ! Do something with cfl..
    print *, "TODO CFL"
    dt_max = 1.0e-12_dp
  end function PM_get_max_dt

  ! Divide the particles over the tasks based on their cell index
  subroutine PM_divide_particles(pc, myrank, ntasks)
    include 'mpif.h'
    type(PC_t), intent(inout) :: pc
    integer, intent(in) :: myrank, ntasks

    integer :: ll, ix, rank, ierr, iMin, n_behind, n_part_interval, n_send
    integer :: tag, part_count
    integer :: n_sends, n_recvs, sender, recver, n_part_mean
    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: n_part_task(0:ntasks-1)
    integer :: n_part_interval_task(0:ntasks-1, 0:ntasks-1)
    integer :: splitIx(-1:ntasks)
    integer :: rel_size_part_t, max_buf_size, req_buf_size
    integer, parameter :: n_cells = 10000
    integer, allocatable :: nPartPerCellIx(:)
    real(dp), allocatable :: cellIxs(:)
    real(dp), allocatable :: r_buf(:)
    type(PC_part_t) :: part_temp(1)
    real(dp) :: x_min, x_max, dx, inv_dx

    ! print *, "Before divide ", myrank, " has ", pc%n_part, " particles"

    ! Get the number of particles each task has
    call MPI_ALLGATHER(pc%n_part, 1, MPI_integer, n_part_task, 1, MPI_integer, MPI_COMM_WORLD, ierr)

    x_min = minval(pc%particles(1:pc%n_part)%x(1))
    x_max = maxval(pc%particles(1:pc%n_part)%x(1))
    call MPI_ALLREDUCE(MPI_IN_PLACE, x_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE(MPI_IN_PLACE, x_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    ! Add small extra space to avoid indexing problems at boundaries
    x_min = x_min * (1-1e-6_dp)
    x_max = x_max * (1+1e-6_dp)
    dx = (x_max - x_min) / n_cells
    inv_dx = 1 / dx

    n_part_mean      = (sum(n_part_task) + ntasks - 1) / ntasks
    allocate( nPartPerCellIx(n_cells) )
    allocate( cellIxs(pc%n_part) )

    nPartPerCellIx = 0

    ! Set the cell indexes for the particles
    do ll = 1, pc%n_part
       ix                   = int((pc%particles(ll)%x(1) - x_min) * inv_dx) + 1
       cellIxs(ll)          = ix
       nPartPerCellIx(ix)   = nPartPerCellIx(ix) + 1
    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, nPartPerCellIx, n_cells, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Sort by cell_ixs
    call heapsortParticleByDbles(pc%particles(1:pc%n_part), cellIxs)

    ! Find the splitIx(:) s.t. task n will hold particles with ix > splitIx(n-1) and ix <= splitIx(n)
    splitIx(-1)                   = 0
    splitIx(ntasks)               = n_cells
    rank                          = 0
    n_behind                   = 0
    n_part_interval_task(:,:)  = 0
    n_part_interval               = 0

    do ix = 1, n_cells
       n_part_interval = n_part_interval + nPartPerCellIx(ix)

       if (n_part_interval >= n_part_mean .or. ix == n_cells) then
          splitIx(rank)     = ix
          n_part_interval   = 0
          part_count           = 0

          do ll = n_behind+1, pc%n_part
             !                print *, ll, cellIxs(ll)
             if (cellIxs(ll) <= ix) then
                part_count = part_count + 1
             else
                exit
             end if
          end do

          n_part_interval_task(rank, myrank) = part_count
          n_behind       = n_behind + part_count
          rank              = rank + 1

          if (rank == ntasks) exit
       end if

    end do

    call MPI_ALLREDUCE(MPI_IN_PLACE, n_part_interval_task, ntasks*ntasks, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Send particles to task that holds a region
    n_sends = 0
    n_recvs = 0

    ! We use a buffer because we cannot send the type directly with Fortran MPI
    rel_size_part_t = storage_size(part_temp(1)) / storage_size(1.0_dp)
    if (rel_size_part_t * storage_size(1.0_dp) /= storage_size(part_temp(1))) then
       print *, "PC_part_t has size not divisible by size of double prec."
       stop
    end if

    max_buf_size = rel_size_part_t * max(&
         maxval(n_part_interval_task(myrank, :)), &
         maxval(n_part_interval_task(:, myrank)))
    allocate(r_buf(max_buf_size))
    print *, "Max buf size", max_buf_size

    do recver = 0, ntasks - 1

       ! Receive the particles in tasks' region from all other tasks
       if (myrank == recver) then

          do sender = 0, ntasks - 1
             if (sender == myrank) cycle ! Don't have to send to ourselves

             n_send = n_part_interval_task(recver, sender)
             req_buf_size = n_send * rel_size_part_t

             if (n_send > 0) then
                n_recvs = n_recvs + 1
                tag    = sender

                print *, myrank, ": receives", n_send, "from", sender, "iMin", iMin
                call MPI_RECV(r_buf, req_buf_size, MPI_DOUBLE_PRECISION, sender, tag, &
                     & MPI_COMM_WORLD, mpi_status, ierr)
                iMin   = pc%n_part + 1
                pc%n_part = pc%n_part + n_send
                call checkNumParticles(pc%n_part)
                pc%particles(iMin:iMin+n_send-1) = &
                     transfer(r_buf(1:req_buf_size), part_temp, n_send)
             end if
          end do


       else ! Send particles to recver
          n_send = n_part_interval_task(recver, myrank)
          req_buf_size = n_send * rel_size_part_t

          if (n_send > 0) then
             n_sends      = n_sends + 1
             tag         = myrank

             ! Find index of first particle to send
             iMin        = sum( n_part_interval_task(0:recver-1, myrank) ) + 1

             print *, myrank, ": sends", n_send, "to", recver, "iMin", iMin
             r_buf(1:req_buf_size) = transfer(pc%particles(iMin:iMin+n_send-1), &
                  r_buf, req_buf_size)
             call MPI_SEND(r_buf, req_buf_size, MPI_DOUBLE_PRECISION, recver, tag, &
                  & MPI_COMM_WORLD, ierr)

             ! Mark the particles that we have send away as inactive
             do ll = iMin, iMin + n_send - 1
                call LL_add(pc%clean_list, ll)
             end do

          end if
       end if

    end do

    ! if (n_recvs > 0) call MPI_WAITALL(n_recvs, recvReqs, MPI_STATUSES_IGNORE, ierr)
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !       print *, "After divide ", myrank, " has ", pc%n_part, "pc%n_part (some dead)"
    call pc%clean_up()
    ! print *, "After divide ", myrank, " has ", pc%n_part, " particles"

  end subroutine PM_divide_particles

end module m_particle

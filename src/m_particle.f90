module m_particle
  use m_particle_core

  implicit none
  private

  integer, parameter :: dp = kind(0.0d0)

  ! Public methods
  public :: PM_particles_to_density
  public :: PM_fld_error
  public :: PM_get_num_part_mpi
  public :: PM_get_max_dt
  public :: PM_bin_func

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

  function PM_bin_func(my_part, real_args) result(i_bin)
    type(PC_part_t), intent(in) :: my_part
    real(dp), intent(in) :: real_args(:)
    integer :: i_bin
    real(dp) :: x, x_min, inv_dx

    x      = my_part%x(1)
    x_min  = real_args(1)
    inv_dx = real_args(2)
    i_bin  = 1 + int((x - x_min) * inv_dx)
  end function PM_bin_func
end module m_particle

!> Module which contains all particle_core modules, so that a user does not have to
!> include them separately.
module m_pc_all
  use m_units_constants
  use m_cross_sec
  use m_gas
  use m_particle_core
  use m_random ! external from lib rng_fortran

  implicit none
  public

end module m_pc_all

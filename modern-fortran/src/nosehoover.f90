module NoseHoover

  use iso_fortran_env, only: int64, real64

  implicit none

  private
  ! parameters
  integer n
  real dt
  real Q
  real T

  ! state
  real eta
  real K

  ! derived parameters
  real dt48
  real dof
  real invQ

contains

  subroutine initialize(n_particles, delta_t, Q, target_T)

    implicit none

    integer, intent(in) :: n_particles
    real, intent(in) :: delta_t
    real, intent(in) :: Q
    real, intent(in) :: target_T

    ! parameters
    n = n_particles
    dt = delta_t
    T = target_T

    ! derived parameters
    dt48 = 48.0*delta_t
    invQ = 1.0/Q
    dof = 3.0*n - 6.0

    ! initialize state
    eta = 0.0
  end subroutine initialize

  subroutine onestep(rx, ry, rz, vx, vy, vz, fx, fy, fz)
    implicit none

    real, allocatable, intent(in out) :: rx(:)
    real, allocatable, intent(in out) :: ry(:)
    real, allocatable, intent(in out) :: rz(:)
    real, allocatable, intent(in out) :: vx(:)
    real, allocatable, intent(in out) :: vy(:)
    real, allocatable, intent(in out) :: vz(:)
    real, allocatable, intent(in) :: fx(:)
    real, allocatable, intent(in) :: fy(:)
    real, allocatable, intent(in) :: fz(:)

    real :: eta1, eta2, ssum

    K = 0.5*dt*eta
    eta1 = 1.0 - K
    eta2 = 1.0/(1.0+K)

    ! update velocities using array operations
    vx = (vx*eta1 + dt48*fx) * eta2
    vy = (vy*eta1 + dt48*fy) * eta2
    vz = (vz*eta1 + dt48*fz) * eta2

    ! update positions
    rx = rx + vx*dt
    ry = ry + vy*dt
    rz = rz + vz*dt

    ! kinetic energy sum (elemental arrays)
    ssum = sum(vx**2 + vy**2 + vz**2)

    eta = eta + (ssum - dof*T) * invQ*dt
  end subroutine onestep

end module NoseHoover

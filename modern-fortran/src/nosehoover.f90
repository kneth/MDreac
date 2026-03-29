module NoseHoover

  use iso_fortran_env, only: int64, real64

  implicit none

  ! parameters
  integer(int64) n
  real(real64) dt
  real(real64) Q
  real(real64) T

  ! state
  real(real64) eta
  real(real64) K

  ! derived parameters
  real(real64) dt48
  real(real64) dof
  real(real64) invQ

contains

  subroutine initalize(n_particles, delta_t, Q, target_T)
    implicit none

    integer(int64), intent(in) n_particles
    real(real64), intent(in) delta_t
    real(real64), intent(in) Q
    real(real64), intent(in) target_T


    ! parameters
    n = n_particles
    dt = delta_t

    ! derived parameters
    dt48 = 48.0*delta_t
    invQ = 1.0/Q
    dof = 3.0*n - 6.0

    ! initialize state
    eta = 0.0
  end subroutine initalize

  subroutine onestep(rx, ry, rz, vx, vy, vz, fx, fy, fz)
    implicit none

    real(real64), intent(in, out) rx(:)
    real(real64), intent(in, out) ry(:)
    real(real64), intent(in, out) rz(:)
    real(real64), intent(in, out) vx(:)
    real(real64), intent(in, out) vy(:)
    real(real64), intent(in, out) vz(:)
    real(real64), intent(in) fx(:)
    real(real64), intent(in) fy(:)
    real(real64), intent(in) fz(:)

    integer(int64) i
    real(real64) sum

    K = 0.5*dt*eta
    eta1 = 1.0 - K
    eta2 = 1.0/(1.0+K)
    sum = 0.0

    do i=1, n
       vx(i) = (vx(i)*eta1 + dt48*fx(i)) * eta2
       vy(i) = (vy(i)*eta1 + dt48*fy(i)) * eta2
       vz(i) = (vz(i)*eta1 + dt48*fz(i)) * eta2
    end do

    do i=1, n
       rx(i) = rx(i) + vx(i)*dt
       ry(i) = ry(i) + vy(i)*dt
       rz(i) = rz(i) + vz(i)*dt
    end do

    do i=1, n
       sum = sum + vx(i)**2 + vy(i)**2 + vz(i)**2
    end do

    eta = eta + (sum-dof*T) * invQ*dt
  end subroutine onestep

end module NoseHoover

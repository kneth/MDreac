program mmd

  use iso_fortran_env, only: int64, real64
  use NoseHoover, only: initialize, onestep
  use RandGen, only: gauss

  implicit none

  ! configuration
  integer :: N
  real :: T, rho

  ! simulation parameters
  integer :: steps
  real :: dt, Q

  ! system
  real, allocatable :: rx(:), ry(:), rz(:)
  real, allocatable :: vx(:), vy(:), vz(:)
  real, allocatable :: fx(:), fy(:), fz(:)

  ! derived parameters
  real :: L

  call ReadConfiguration(N, T, rho, steps, dt, Q)
  call InitializeSystem(N, T, rho, L)
  call NoseHoover%initialize(N, dt, Q, T)

  do i=1, steps
    if (mod(i, 10) == 0) then
      call Links()
      call UpdateLinkCell()
    end if
    call ComputeForces()
    call NoseHoover%onestep(rx, ry, rz, vx, vy, vz, fx, fy, fz)
  end do

end program mmd

  subroutine ReadConfiguration(N, T, rho, steps, dt, Q)
    implicit none
    integer, intent(out) :: N, steps
    real, intent(out) :: T, rho, dt, Q

    ! Read configuration from file or user input
    ! For now, we will use dummy values
    N = 100
    T = 3.0
    rho = 0.8
    steps = 1000
    dt = 0.001
    Q = 1.0
  end subroutine ReadConfiguration

  subroutine InitializeSystem(N, T, rho, L)
    implicit none

    integer, intent(in) :: N
    real, intent(in) :: T, rho
    real, intent(out) :: L

    real :: sumv, factor

    ! Allocate memory for system arrays
    allocate(rx(N), ry(N), rz(N))
    allocate(vx(N), vy(N), vz(N))
    allocate(fx(N), fy(N), fz(N))

    L = (real(N) / rho)**(1.0/3.0)

    ! Initialize positions (even spacing along x for now)
    rx = real([(i-1, i=1,N)]) * (L / N)
    ry = 0.0
    rz = 0.0

    ! Fill velocities with Gaussian numbers and scale to target temperature
    vx = [(gauss(0.0, 1.0), i=1,N)]
    vy = [(gauss(0.0, 1.0), i=1,N)]
    vz = [(gauss(0.0, 1.0), i=1,N)]

    sumv = sum(vx**2 + vy**2 + vz**2)
    factor = sqrt(3.0 * (real(N) - 1.0) * T / sumv)
    vx = vx * factor
    vy = vy * factor
    vz = vz * factor

    ! initialize forces; set them to zero
    fx = 0.0
    fy = 0.0
    fz = 0.0

  end subroutine InitializeSystem
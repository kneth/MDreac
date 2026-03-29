program mmd

  use iso_fortran_env, only: int64, real64

  implicit none

  real(real64) pi
  real(real64) r
  real(real64) a

  pi = 3.1415926535
  r = 2.5
  a = pi * r**2

  print *, "Hello, World", a

end program mmd

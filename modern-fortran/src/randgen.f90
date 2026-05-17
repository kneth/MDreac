! random number generator

module RandGen
  implicit none

contains
    function gauss(mean, stddev) result(r)
        real, intent(in) :: mean, stddev
        real :: r

        ! Generate a Gaussian random variable using Box-Muller transform
        real :: u1, u2
        real, parameter :: pi = 3.14159265358979323846

        call random_number(u1)
        call random_number(u2)
        r = mean + stddev * sqrt(-2.0 * log(u1)) * cos(2.0 * pi * u2)
    end function gauss

end module RandGen

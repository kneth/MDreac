module LinkCells
  implicit none

  private
    integer, allocatable :: head(:)

contains

    subroutine Links
        implicit none
        integer :: i, M

        do i = 1, M
            head(i) = 0
        end do


    end subroutine Links

end module LinkCells
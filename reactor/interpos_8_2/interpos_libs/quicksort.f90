! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module qsort_c_module

implicit none
public :: QsortC
private :: Partition

contains

recursive subroutine QsortC(A)
  real, dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call Partition(A, iq)
     call QsortC(A(:iq-1))
     call QsortC(A(iq:))
  endif
end subroutine QsortC

subroutine Partition(A, marker)
  real, dimension(:) :: A
  integer :: marker
  integer :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine Partition

end module qsort_c_module

program sortdriver
  use qsort_c_module
  implicit none
  integer, parameter :: r = 10
  real, dimension(1:r) :: myarray = &        ! (1:r)
     (/0, 50, 20, 25, 90, 10, 5, 1, 99, 75/)
  print *, "myarray is ", myarray
  call QsortC(myarray)
  print *, "sorted array is ", myarray
end program sortdriver

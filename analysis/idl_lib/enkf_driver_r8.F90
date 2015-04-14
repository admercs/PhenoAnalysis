subroutine enkf_driver_r8(ndim, nrens, nrobs, A, HA, obsdat, obserr) 

  use enkf_module_r8

  implicit none

  ! arguments
  integer, intent(in) :: ndim, nrens, nrobs
  real(kind=8), intent(inout) :: A(ndim,nrens)
  real(kind=8), intent(inout) :: HA(nrobs,nrens)
  real(kind=8), intent(in) :: obsdat(nrobs)
  real(kind=8), intent(in) :: obserr(nrobs)

  ! local variables
  integer :: i, s, clock
  integer, allocatable :: seed(:)

  ! parameters
  logical, parameter :: use_sqrt = .true.
  
  ! init random number generator
  call random_seed(size = s)
  allocate(seed(s))
  call system_clock(count=clock)
  i=1
  seed = clock + 37 * (/ (i - 1, i = 1, s) /)

  ! uncomment line below to generate 1:1 reproducible
  ! assimilation experiments (it always initializes 
  ! the random number generator with the same numbers)
  seed(:) = 200

  call random_seed(put = seed)
  deallocate(seed)
 
  write(*,'(I8,I8,I8)') ndim, nrens, nrobs
  
  if (use_sqrt) then
     allocate(ROT(nrens,nrens))
     call calc_rot(nrens)
     call enkf(HA, obsdat, obserr, 23, A=A)
     deallocate(ROT)
  else
     call enkf(HA, obsdat, obserr, 13, A=A)
  endif
 
end subroutine enkf_driver_r8

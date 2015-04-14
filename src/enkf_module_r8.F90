! Copyleft 2006-2014 NASA and Blue Marble Research (http://www.bluemarble.ch)
! Author(s): Reto Stockli
! 
! This file is part of <phenoanalysis> and was started within the
! NASA Energy and Water Cycle Study (NEWS) grant No. NNG06CG42G. It
! now is a open source software project with code and documentation
! found on: http://phenoanalysis.sourceforge.net
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! IMPORTANT NOTE: Intellectual Property Right and Author statements in the
! some of the code may supersede the above Copyleft / GNU GPL license.

module enkf_module_r8

  ! 2013/03/12 Reto Stockli (Blue Marble Research) 

  ! based on Geir Evensen's 2003 paper and his public code from 2007
  ! This modules is self-consistent and includes the full ensemble kalman filter 
 
  ! DEPENDENCIES: 
  ! lapack/BLAS compatible linear algebra routines. 
  ! - Lapack is available here: http://www.netlib.org/lapack/
  ! - BLAS is available here: http://www.netlib.org/blas/
  ! - Both Lapack and BLAS are part of OSX by use of link flag -framework veclib.
  ! - Some Lapack/BLAS routines can be substituted by IBM's optimized ESSL package on XLF compilers

  ! Ensemble rotation matrix used during SQRT EnKF
  real(kind=8), allocatable :: ROT(:,:)

  contains

subroutine calc_rot(nrens)
  ! update random rotation only for first process when using 
  ! distributed global or several local analyses since all grid
  ! points need to use the same rotation.

  ! only needed for SQRT EnKF
  ! Make sure to allocate ROT(nrens,nrens) externally before calling this routine!

  implicit none

  ! Arguments
  integer, intent(in) :: nrens

  ! Local variables
  integer :: i

  ROT=0.d0
  do i=1,nrens
     ROT(i,i)=1.d0
  enddo
  ! call randrot(ROT,nrens)  ! old non mean preserving random rotation
  call mean_preserving_rotation(ROT,nrens)

end subroutine calc_rot

subroutine enkf(HA, obsdat, obserr, mode, X, A)
  ! Computes the analysed ensemble for A using the EnKF or square root schemes.
  ! This code is mostly based on the one by G. Evensen published on his website in June 2007

  ! Please allocate ROT(nrens,nrens) and call calc_rot() before using this subroutine with mode=?3 (SQRT EnKF)

  implicit none

  ! arguments
  real(kind=8), intent(inout) :: HA(:,:)  ! model states and parameters in observation space
  real(kind=8), intent(in) :: obsdat(:)   ! observation vector
  real(kind=8), intent(in) :: obserr(:)   ! observation variance vector
  integer, intent(in) :: mode             
  ! first integer means:
  !    1=EnKF
  !    2=SQRT
  ! Second integer is pseudo inversion:
  !    1=eigen value pseudo inversion of SS'+(N-1)R
  !    2=SVD subspace pseudo inversion of SS'+(N-1)R
  !    3=SVD subspace pseudo inversion of SS'+EE'
  real(kind=8), intent(out), optional :: X(:,:)     ! ensemble covariance matrix for updating A externally
  real(kind=8), intent(inout), optional :: A(:,:)   ! model state and parameter matrix for internal update

  ! local variables
  integer :: ndim              ! dimension of model state
  integer :: nrens             ! number of ensemble members
  integer :: nrobs             ! number of observations

  integer :: i,j,k,m
  integer :: nrmin
  integer :: iblkmax
  logical :: lreps
  logical :: use_sqrt          ! use Square-Root EnKF instead of regular EnKF

  real(kind=8) :: ensemble_mean
  real(kind=8) :: observation_variance
  real(kind=8), allocatable :: HAmean(:)
  real(kind=8), allocatable :: observation_mean(:)

  ! internal EnKF arrays
  real(kind=8), allocatable :: eig(:)
  real(kind=8), allocatable :: W(:,:)
  real(kind=8), allocatable :: X2(:,:)
  real(kind=8), allocatable :: X3(:,:)
  real(kind=8), allocatable :: Reps(:,:)
  real(kind=8), allocatable :: X5(:,:)      ! Ensemble covariance matrix (coefficients)
  real(kind=8), allocatable :: E(:,:)       ! observation error covariance matrix
  real(kind=8), allocatable :: IM(:,:)      ! innovation matrix
  real(kind=8), allocatable :: D(:,:)       ! perturbed observation matrix
  real(kind=8), allocatable :: R(:,:)       ! matrix holding error covariance matrix for observations R 
  real(kind=8), allocatable :: innov(:)     ! vector holding d-H*mean(A)

  ! parameters
  logical, parameter :: verbose = .false.           ! Printing some diagnostic output
  real(kind=8), parameter :: truncation = 0.99      ! The ratio of variance retained in pseudo inversion (0.99)  

  ! check mode
  select case (mode)
  case (11,12,13)
     if (verbose) write(*,'(A)') 'analysis: Regular EnKF'
     use_sqrt = .false.
  case (21,22,23)
     if (verbose) write(*,'(A)') 'analysis: SQRT EnKF'
     use_sqrt = .true.
  case default
     write(*,'(A,I4)') 'analysis: unknown flag for mode: ',mode
     return
  end select

  ! get and check dimensions
  nrobs = size(HA,1)
  nrens = size(HA,2)

  if (present(A)) then
     if (verbose) write(*,'(A)') 'analysis: updating States'
     ndim = size(A,1)
  endif

  if (size(obsdat).ne.nrobs) then
     write(*,'(A,I8,I8)') 'analysis: nrobs of obsdat does not match HA : ',size(obsdat),nrobs
     return
  endif
  if (size(obserr).ne.nrobs) then
     write(*,'(A,I8,I8)') 'analysis: nrobs of obserr does not match HA : ',size(obserr),nrobs
     return
  endif

  ! allocate arrays
  if ((.not.use_sqrt).or.(nrobs.eq.1)) then
     allocate(D(nrobs,nrens))
     allocate(IM(nrobs,nrens))
  endif

  allocate(E(nrobs,nrens))

  ! prepare EnKF analysis

  ! Construct observation perturbations E from 
  ! random fields and observation uncertainty vector
  call random(E,nrobs*nrens)
  
  allocate(observation_mean(nrobs))
  observation_mean=0.d0
  observation_variance=0.d0
  do j=1,nrens
     do m=1,nrobs
        observation_mean(m)=observation_mean(m)+E(m,j)
     enddo
  enddo
  observation_mean=observation_mean/float(nrens)
  
  do j=1,nrens
     do m=1,nrobs
        E(m,j)=E(m,j)-observation_mean(m)
        observation_variance = observation_variance + E(m,j)**2.
        E(m,j)=sqrt(obserr(m))*E(m,j)
     enddo
  enddo
  deallocate(observation_mean)

  observation_variance=sqrt(observation_variance/dble(nrens*nrobs))
  E=E/observation_variance

  ! compute covariance matrix of observation errors
  if ((mode.eq.11).or.(mode.eq.12).or.(mode.eq.21).or.(mode.eq.22).or.(nrobs.eq.1)) then
     allocate(R(nrobs,nrobs)) 
     R=matmul(E,transpose(E))/dble(nrens-1)
  endif

  ! full observation matrix is not needed any more in modes ?3 of the analysis
  ! SQRT EnKF saves memory for heavy observation counts
  if ((.not.use_sqrt).or.(nrobs.eq.1)) then

     ! Construct ensemble of measurements D=d+E
     do j=1,nrens
        do m=1,nrobs
           D(m,j)=obsdat(m)+E(m,j)
        enddo
     enddo
     
     ! Compute innovation IM=D-HA
     IM=D-HA
  endif

  ! Compute mean(HA) 
  allocate(HAmean(nrobs))
  HAmean=0.d0
  do j=1,nrens
     do m=1,nrobs
        HAmean(m)=HAmean(m)+HA(m,j)
     enddo
  enddo
  HAmean=(1.d0/dble(nrens))*HAmean

  ! compute I = D - H*mean(A)
  if (use_sqrt) then
     allocate(innov(nrobs))
     do k = 1,nrobs
        call ensmean(HA(k,:),ensemble_mean)
        innov(k) = obsdat(k) - ensemble_mean
     end do
  endif

  ! Compute HA'=HA-mean(HA)
  do j=1,nrens
     HA(:,j)=HA(:,j)-HAmean(:)
  enddo
  deallocate(HAmean)

  lreps=.false.
  if (verbose) write(*,'(A)') 'analysis: verbose is on'

  ! Pseudo inversion of C=SS' +(N-1)*R
  if (verbose) write(*,'(A)') 'analysis: Inversion of C'
  if (nrobs == 1) then
     nrmin=1
     allocate(W(1,1))
     allocate(eig(1))
     eig(1)=dot_product(HA(1,:),HA(1,:))+dble(nrens-1)*R(1,1)
     eig(1)=1.d0/eig(1)
     W(1,1)=1.d0

  else
     select case (mode)
     case(11,21)
        nrmin=nrobs
        ! Evaluate R= S*S` + (nrens-1)*R

        call dgemm('n','t',nrobs,nrobs,nrens, &
             1.d0, HA, nrobs, &
             HA, nrobs, &
             dble(nrens-1), R, nrobs)
        
        ! Compute eigenvalue decomposition of R -> W*eig*W` 
        allocate(W(nrobs,nrobs))
        allocate(eig(nrobs))
        call eigC(R,nrobs,W,eig)
        call eigsign(eig,nrobs,truncation)

     case(12,22)
        nrmin=min(nrobs,nrens)
        allocate(W(nrobs,nrmin))
        allocate(eig(nrmin))

        call lowrankCinv(HA,R,nrobs,nrens,nrmin,W,eig,truncation)
       
     case(13,23)
        nrmin=min(nrobs,nrens)
        allocate(W(nrobs,nrmin))
        allocate(eig(nrmin))
        call lowrankE(HA,E,nrobs,nrens,nrmin,W,eig,truncation)
       
     case default
        write(*,'(A,I4)') 'analysis: Unknown mode: ',mode
     end select
  endif

  if ((mode.eq.11).or.(mode.eq.12).or.(mode.eq.21).or.(mode.eq.22).or.(nrobs.eq.1)) then
     deallocate(R)
  endif
  deallocate(E)

  ! Generation of X5 (or representers in EnKF case with few measurements)
  if (verbose) write(*,'(A)') 'analysis: Generation of X5'
  allocate(X5(nrens,nrens))

  select case (mode)
  case(11,12,13)
     allocate(X3(nrobs,nrens))
     if (nrobs > 1) then
        call genX3(nrens,nrobs,nrmin,eig,W,IM,X3)
     else
        X3=IM*eig(1)
     endif

     deallocate(eig)
     deallocate(W)

!     if (2_8*ndim*nrobs < 1_8*nrens*(nrobs+ndim)) then
!        !        Code for few observations ( m<nN/(2n-N) )
!        if (verbose) write(*,'(A)') 'analysis: Representer approach is used'
!        lreps=.true.
!        allocate (Reps(ndim,nrobs))
        !        Reps=matmul(A,transpose(S))
!        call dgemm('n','t',ndim,nrobs,nrens,1.d0,A,ndim,HA,nrobs,0.d0,Reps,ndim)
!     else
        if (verbose) write(*,'(A)') 'analysis: X5 approach is used'
        !        X5=matmul(transpose(HA),X3)
        call dgemm('t','n',nrens,nrens,nrobs,1.d0,HA,nrobs,X3,nrobs,0.d0,X5,nrens)
        deallocate(X3)
        do i=1,nrens
           X5(i,i)=X5(i,i)+1.d0
        enddo
!     endif
     
  case(21,22,23)
     ! Mean part of X5
     call meanX5(nrens,nrobs,nrmin,HA,W,eig,innov,X5)
  
     deallocate(innov)
     
     ! Generating X2
     allocate(X2(nrmin,nrens))

     call genX2(nrens,nrobs,nrmin,HA,W,eig,X2)

     deallocate(eig)
     deallocate(W)
     
     ! Generating X5 matrix
     call X5sqrt(X2,nrobs,nrens,nrmin,X5,mode)

     deallocate(X2)     
  case default
     write(*,'(A,I4)') 'analysis: Unknown flag for mode: ',mode
  end select

  ! Final ensemble update
  if (present(A)) then
     if (verbose) write(*,'(A)') 'analysis: Final ensemble update'
     if (lreps) then
        A=A+matmul(Reps,X3)
        call dgemm('n','n',ndim,nrens,nrobs,1.d0,Reps,ndim,X3,nrobs,1.d0,A,ndim)
        deallocate(X3)
        deallocate(Reps)
     else
        iblkmax=min(ndim,200)
        call multa(A, X5, ndim, nrens, iblkmax )
     endif
  endif

  ! Return ensemble matrix for updating states externally
  if (present(X)) then
     X = X5
  endif
  deallocate(X5)

  ! clean up
  if ((.not.use_sqrt).or.(nrobs.eq.1)) then
     deallocate(D)
     deallocate(IM)
  endif

end subroutine enkf

subroutine kernel_smooth(A)
  ! updates vector A with a kernel smoothing algorithm by West 1993

  implicit none

  ! arguments
  real(kind=8), intent(inout) :: A(:)

  ! local variables
  real(kind=8) :: ensemble_mean, ensemble_variance
  real(kind=8), allocatable :: random_ensemble(:)
  integer :: nrens
  real(kind=8), parameter :: alpha = 0.995 ! shrinkage factor

  nrens = size(A,1)

  allocate(random_ensemble(nrens))

  call ensmean(A,ensemble_mean)
  call ensmean(alpha*A + (1.-alpha)*ensemble_mean,ensemble_mean)
  call ensvar(A,ensemble_mean,ensemble_variance)
  call random(random_ensemble,nrens)
  A= ensemble_mean  + (1.- alpha**2)*random_ensemble*ensemble_variance

  deallocate(random_ensemble)

end subroutine kernel_smooth

subroutine kernel_inflate(A,var,alpha)
  ! inflates variance of vector A to the minimum variance var*alpha
  ! (nothing changes if variance is above var*alpha)

  implicit none

  ! arguments
  real(kind=8), intent(inout) :: A(:)
  real(kind=8), intent(in) :: var
  real(kind=8), intent(in) :: alpha ! minimum variance factor

  ! local variables
  real(kind=8) :: ensemble_mean, ensemble_variance
  real(kind=8) :: varscale
  integer :: nrens

  nrens = size(A,1)

  call ensmean(A,ensemble_mean)
  call ensvar(A,ensemble_mean,ensemble_variance)

  varscale = max(var/max(ensemble_variance,1.d-8)*alpha,1.d0)

  A = (A - ensemble_mean) * sqrt(varscale) + ensemble_mean

end subroutine kernel_inflate

subroutine kernel_deflate(A,var,alpha)
  ! deflates variance of vector A to the maximum of variance var*alpha
  ! (nothing changes if variance is below var*alpha)

  implicit none

  ! arguments
  real(kind=8), intent(inout) :: A(:)
  real(kind=8), intent(in) :: var ! prior ensemble variance
  real(kind=8), intent(in) :: alpha ! minimum variance factor

  ! local variables
  real(kind=8) :: ensemble_mean, ensemble_variance
  real(kind=8) :: varscale
  integer :: nrens

  nrens = size(A,1)

  call ensmean(A,ensemble_mean)
  call ensvar(A,ensemble_mean,ensemble_variance)

  varscale = min(var/max(ensemble_variance,1.d-8)*alpha,1.d0)

  A = (A - ensemble_mean) * sqrt(varscale) + ensemble_mean

end subroutine kernel_deflate

subroutine kernel_rescale(A,varscale)
  ! rescales variance of vector A with scale factor varscale

  implicit none

  ! arguments
  real(kind=8), intent(inout) :: A(:)
  real(kind=8), intent(in) :: varscale ! variance scale factor

  ! local variables
  real(kind=8) :: ensemble_mean
  integer :: nrens

  nrens = size(A,1)

  call ensmean(A,ensemble_mean)

  A = (A - ensemble_mean) * sqrt(varscale) + ensemble_mean

end subroutine kernel_rescale

subroutine kernel_limit(A,A_min,A_max)
  ! limit the mean of vector by lower and upper bounds. If limits are crossed, then the 
  ! mean of the vector is shifted back to the limit. Individual members of the vector can 
  ! therefore still be beyond limits. This guarantees that the variance of Vector A is maintained.
  
  implicit none

  ! arguments
  real(kind=8), intent(inout) :: A(:)
  real(kind=8), intent(in) :: A_min
  real(kind=8), intent(in) :: A_max

  ! local variables
  real(kind=8) :: A_mean

  call ensmean(A,A_mean)

  if (A_mean.lt.A_min) A = A + (A_min - A_mean)
  if (A_mean.gt.A_max) A = A + (A_max - A_mean)

end subroutine kernel_limit

subroutine ensmean(A,ave)
  ! calculates the ensemble mean of a vector A

  implicit none

  ! arguments
  real(kind=8), intent(in)  :: A(:)
  real(kind=8), intent(out) :: ave

  ! local variables
  integer :: nrens
  integer :: j
  
  nrens = size(A,1)

  ave=A(1)
  do j=2,nrens
     ave=ave+A(j)
  enddo
  ave=(1.d0/dble(nrens))*ave
  
end subroutine ensmean

subroutine ensvar(A,ave,var)
  ! calculates the ensemble variance of a vector A

  implicit none

  ! arguments
  real(kind=8), intent(in)  :: A(:)
  real(kind=8), intent(in)  :: ave
  real(kind=8), intent(out) :: var

  ! local variables
  integer :: nrens
  integer :: j
  
  nrens = size(A,1)

  var=0.d0
  do j=1,nrens
     var=var+(A(j)-ave)*(A(j)-ave)
  enddo
  var=(1.d0/dble(nrens-1))*var
  
end subroutine ensvar

subroutine multa(A, X, ndim, nrens, iblkmax)

  implicit none

  ! arguments
  integer, intent(in) :: ndim
  integer, intent(in) :: nrens
  integer, intent(in) :: iblkmax
  real(kind=8), intent(in)    :: X(nrens,nrens)
  real(kind=8), intent(inout) :: A(ndim,nrens)

  ! local variables
  real(kind=8) :: v(iblkmax,nrens)  ! Automatic work array 
  integer :: ia,ib

  do ia = 1,ndim,iblkmax
     ib = min(ia+iblkmax-1,ndim)
     v(1:ib-ia+1,1:nrens) = A(ia:ib,1:nrens)
     call dgemm('n','n', ib-ia+1, nrens, nrens, &
          1.d0, v(1,1), iblkmax, &
          X(1,1), nrens, &
          0.d0, A(ia,1), ndim)
  enddo

end subroutine multa

subroutine random(work1,n)
  !  Returns a vector of random values N(variance=1,mean=0)

  implicit none

  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: work1(n)

  ! local variables
  real(kind=8), allocatable :: work2(:)
  real(kind=8), parameter   ::  pi=3.141592653589d0
  
  allocate (work2(n))
  
  call random_number(work1)
  call random_number(work2)

  work1 = sqrt(-2.d0*log(work1+tiny(0.d0)))*cos(2.d0*pi*work2)
  
  deallocate(work2)

end subroutine random

subroutine randrot(Q,nrens)

  implicit none

  ! arguments
  integer, intent(in)  :: nrens
  real(kind=8),    intent(out) :: Q(nrens,nrens)
  
  ! local variables
  real(kind=8) ::  A(nrens,nrens)
  real(kind=8) :: B(nrens,nrens)
  real(kind=8) :: sigma(nrens)
  real(kind=8) :: work(10*nrens)
  real(kind=8), parameter :: pi=3.14159253589
  integer :: ierr
  
  call random_number(B)
  call random_number(A)
  Q = sqrt(-2.d0*log(A+tiny(0.d0))) * cos(2.d0*pi*B)

!$OMP END CRITICAL
! QR factorization
  call dgeqrf(nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
  if (ierr /= 0) write(*,'(A,I6)')  'randrot: dgeqrf ierr=',ierr
  
! Construction of Q
  call dorgqr(nrens, nrens, nrens, Q, nrens, sigma, work, 10*nrens, ierr )
  if (ierr /= 0) write(*,'(A,I6)')  'randrot: dorgqr ierr=',ierr
!$OMP END CRITICAL

end subroutine randrot

subroutine lowrankE(S,E,nrobs,nrens,nrmin,W,eig,truncation)

  implicit none

  ! arguments
  integer, intent(in)  :: nrobs
  integer, intent(in)  :: nrens
  integer, intent(in)  :: nrmin
  real(kind=8),    intent(in)  :: S(nrobs,nrens)
  real(kind=8),    intent(in)  :: E(nrobs,nrens)
  real(kind=8),    intent(out) :: W(nrobs,nrmin)
  real(kind=8),    intent(out) :: eig(nrmin)
  real(kind=8),    intent(in)  :: truncation
  
  ! local variables
  real(kind=8) :: U0(nrobs,nrmin)
  real(kind=8) :: sig0(nrmin)
  real(kind=8) :: X0(nrmin,nrens)
  real(kind=8) :: U1(nrmin,nrmin)
  real(kind=8) :: VT1(1,1)
  integer :: i,j  
  integer :: lwork
  integer :: ierr
  real(kind=8), allocatable :: work(:)

  ! Compute SVD of S=HA`  ->  U0, sig0

  call  svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)

  ! Compute X0=sig0^{*T} U0^T E 

  ! X0= U0^T R
  call dgemm('t','n',nrmin,nrens,nrobs, 1.d0,U0,nrobs, E,nrobs, 0.d0,X0,nrmin)
    
  do j=1,nrens
     do i=1,nrmin
        X0(i,j)=sig0(i)*X0(i,j)
     enddo
  enddo
  
  ! Compute singular value decomposition  of X0(nrmin,nrens)
  lwork=2*max(3*nrens+nrobs,5*nrens)
  allocate(work(lwork))
  eig=0.d0
  
  call dgesvd('S', 'N', nrmin, nrens, X0, nrmin, eig, U1, nrmin, VT1, 1, work, lwork, ierr)
  deallocate(work)
  if (ierr /= 0) then
     write(*,'(A,I6)') '(lowrankE): ierr from call dgesvd 1= ',ierr
     stop
  endif
  
  do i=1,nrmin
     eig(i)=1.d0/(1.d0+eig(i)**2)
  enddo
  
  ! W = U0 * sig0^{-1} * U1
  do j=1,nrmin
     do i=1,nrmin
        U1(i,j)=sig0(i)*U1(i,j)
     enddo
  enddo
  
  call dgemm('n','n',nrobs,nrmin,nrmin, 1.d0,U0,nrobs, U1,nrmin, 0.d0,W,nrobs)

end subroutine

subroutine eigC(R,nrobs,Z,eig)
  ! Compute eigenvalue decomposition of R -> Z*eig*Z` 

  implicit none

  ! arguments
  integer, intent(in) :: nrobs
  real(kind=8), intent(in)    :: R(nrobs,nrobs)
  real(kind=8), intent(out) :: Z(nrobs,nrobs)
  real(kind=8), intent(out)   :: eig(nrobs)

  ! local variables
  real(kind=8) ::  RR(nrobs,nrobs)
  real(kind=8) :: fwork(8*nrobs)  
  integer :: iwork(5*nrobs)
  integer :: ifail(nrobs)
  real(kind=8) :: abstol,ddum
  integer :: idum,neig,ierr
  real(kind=8), external :: DLAMCH
  
!#ifdef XLF
!  real(kind=8), allocatable :: ap(:)
!  integer k
!#endif
  
  idum=1
  
!#ifdef XLF
! Upper packed storage as in ESSL manual
!  allocate (ap(nrobs*(nrobs+1)/2) )
!  k=0
!  do j=1,nrobs
!     do i=1,j
!        k=k+1
!        ap(k)=R(i,j)
!     enddo
!  enddo
!  call dspev(21,ap,eig,Z,nrobs,nrobs,fwork,2*nrobs)
!  deallocate(ap)
!#else
  abstol=2.d0*DLAMCH('S')
  RR=R
  call dsyevx('V', 'A', 'U', nrobs, RR, nrobs, ddum, ddum, idum, idum, abstol, &
       neig, eig, Z, nrobs, fwork, 8*nrobs, iwork, ifail, ierr )
!#endif
  
end subroutine

subroutine eigsign(eig,nrobs,truncation)
  ! Returns the inverse of the truncated eigenvalue spectrum

  implicit none

  ! arguments
  integer, intent(in)    :: nrobs
  real(kind=8),    intent(inout) :: eig(nrobs)
  real(kind=8),    intent(in)    :: truncation
  
  ! local variables
  integer :: i,nrsigma
  real(kind=8) :: sigsum,sigsum1
  
  ! Significant eigenvalues
  sigsum=sum( eig(1:nrobs) )
  sigsum1=0.d0
  nrsigma=0
  do i=nrobs,1,-1
     if (sigsum1/sigsum < truncation) then
        nrsigma=nrsigma+1
        sigsum1=sigsum1+eig(i)
        eig(i) = 1.d0/eig(i)
     else
        eig(1:i)=0.d0
        exit
     endif
  enddo

end subroutine

subroutine genX2(nrens,nrobs,i_dim,S,W,eig,X2)
  ! Generate X2= (I+eig)^{-0.5} * W^T * S

  implicit none
  integer, intent(in) :: nrens
  integer, intent(in) :: nrobs
  integer, intent(in) :: i_dim ! i_dim=nrobs for A4 and nrmin for A5
  real(kind=8), intent(in)    :: W(i_dim,nrens)
  real(kind=8), intent(in)    :: S(nrobs,nrens)
  real(kind=8), intent(in)    :: eig(i_dim)
  real(kind=8), intent(out)   :: X2(i_dim,nrens)
  integer :: i,j
  
  call dgemm('t','n',i_dim,nrens,nrobs,1.d0,W,nrobs, S,nrobs, 0.d0,X2,i_dim)
  
  do j=1,nrens
     do i=1,i_dim
        X2(i,j)=sqrt(eig(i))*X2(i,j)
     enddo
  enddo
  
end subroutine genX2

subroutine genX3(nrens,nrobs,nrmin,eig,W,D,X3)

  implicit none

  ! arguments
  integer, intent(in) :: nrens
  integer, intent(in) :: nrobs
  integer, intent(in) :: nrmin
  real(kind=8),    intent(in) :: eig(nrmin)
  real(kind=8),    intent(in) :: W(nrobs,nrmin)
  real(kind=8),    intent(in) :: D(nrobs,nrens)
  real(kind=8),    intent(out) :: X3(nrobs,nrmin)
  
  ! local variables
  real(kind=8) :: X1(nrmin,nrobs)
  real(kind=8) :: X2(nrmin,nrens)
  integer :: i,j
  
  do i=1,nrmin
     do j=1,nrobs
        X1(i,j)=eig(i)*W(j,i)
     enddo
  enddo
  
!     X2=matmul(X1,D)
  call dgemm('n','n',nrmin,nrens,nrobs,1.d0,X1,nrmin,D ,nrobs,0.d0,X2,nrmin)

!     X3=matmul(W,X2)
  call dgemm('n','n',nrobs,nrens,nrmin,1.d0,W ,nrobs,X2,nrmin,0.d0,X3,nrobs)

end subroutine genX3

subroutine meanX5(nrens,nrobs,nrmin,S,W,eig,innov1,X5)

  implicit none

  ! arguments
  integer, intent(in) :: nrens
  integer, intent(in) :: nrobs
  integer, intent(in) :: nrmin
  real(kind=8), intent(in)    :: W(nrmin,nrmin)
  real(kind=8), intent(in)    :: S(nrobs,nrens)
  real(kind=8), intent(in)    :: eig(nrmin)
  real(kind=8), intent(in)    :: innov1(nrobs)
  real(kind=8), intent(out)   :: X5(nrens,nrens)
  
  ! local variables
  real(kind=8) :: y1(nrmin)
  real(kind=8) :: y2(nrmin)
  real(kind=8) :: y3(nrobs)
  real(kind=8) :: y4(nrens) 
  integer :: i
  
  if (nrobs==1) then
     y1(1)=W(1,1)*innov1(1)
     y2(1)=eig(1)*y1(1)
     y3(1)=W(1,1)*y2(1)
     y4(:)=y3(1)*S(1,:)
  else
     call dgemv('t',nrobs,nrmin,1.d0,W,nrobs,innov1,1,0.d0,y1 ,1)
     y2=eig*y1  
     call dgemv('n',nrobs,nrmin,1.d0,W ,nrobs,y2,1,0.d0,y3 ,1)
     call dgemv('t',nrobs,nrens,1.d0,S ,nrobs,y3,1,0.d0,y4 ,1)
  endif
  
  do i=1,nrens
     X5(:,i)=y4(:)
  enddo

! X5=enN + (I - enN) X5  = enN + X5
  X5=1.d0/dble(nrens) + X5

end subroutine meanX5

subroutine X5sqrt(X2,nrobs,nrens,nrmin,X5,mode)

  implicit none

  ! arguments
  integer, intent(in) :: nrobs
  integer, intent(in) :: nrens
  integer, intent(inout) :: nrmin ! note that nrmin=nrobs in a4
  real(kind=8), intent(in)    :: X2(nrmin,nrens)
  real(kind=8), intent(inout) :: X5(nrens,nrens)
  integer, intent(in) :: mode
  
  ! local variables
  real(kind=8) :: X3(nrens,nrens)
  real(kind=8) :: X33(nrens,nrens)
  real(kind=8) :: X4(nrens,nrens)
  real(kind=8) :: IenN(nrens,nrens)    
  real(kind=8) :: U(nrmin,1),sig(nrmin),VT(nrens,nrens)
  integer :: i,j,lwork,ierr
  real(kind=8), allocatable, dimension(:)   :: work,isigma
  
  ! SVD of X2
  lwork=2*max(3*nrens+nrens,5*nrens)
 
  allocate(work(lwork))
  sig=0.d0
  call dgesvd('N', 'A', nrmin, nrens, X2, nrmin, sig, U, nrmin, VT, nrens, work, lwork, ierr)
  deallocate(work)
  if (ierr /= 0) then
     write(*,'(A,I6)') 'X5sqrt: ierr from call dgesvd = ',ierr
     stop
  endif
  
  if (mode == 21) nrmin=min(nrens,nrobs)
  allocate(isigma(nrmin))
  isigma=1.d0
  do i=1,nrmin
     if ( sig(i) > 1.d0 ) write(*,'(A,I6,G14.6)') 'X5sqrt: WARNING (m_X5sqrt): sig > 1',i,sig(i)
     isigma(i)=sqrt( max(1.d0-sig(i)**2,0.d0) )
  enddo
  
  do j=1,nrens
     X3(:,j)=VT(j,:)
  enddo
  
  do j=1,nrmin
     X3(:,j)=X3(:,j)*isigma(j)
  enddo

  ! Multiply  X3* V' = (V*sqrt(I-sigma*sigma) * V' to ensure symmetric sqrt and 
  ! mean preserving rotation.   Sakov paper eq 13
  call dgemm('n','n',nrens,nrens,nrens,1.d0,X3,nrens,VT,nrens,0.d0,X33,nrens)
  ! Check mean preservation X33*1_N = a* 1_N (Sakov paper eq 15)
  !   do i=1,nrens
  !      write(*,'(A,I6,G14.6)') 'sum(X33)= ',i,sum(X33(i,:))
  !   enddo
   
  !!X33=X3

  call dgemm('n','n',nrens,nrens,nrens,1.d0,X33,nrens,ROT,nrens,0.d0,X4,nrens)
  
  IenN=-1.d0/dble(nrens)
  do i=1,nrens
     IenN(i,i)=  IenN(i,i) + 1.d0
  enddo
  
  call dgemm('n','n',nrens,nrens,nrens,1.d0,IenN,nrens,X4,nrens,1.d0,X5,nrens)
  
  deallocate(isigma)
  
end subroutine X5sqrt

subroutine lowrankCinv(S,R,nrobs,nrens,nrmin,W,eig,truncation)

  implicit none

  ! arguments
  integer, intent(in)  :: nrobs
  integer, intent(in)  :: nrens
  integer, intent(in)  :: nrmin
  real(kind=8),    intent(in)  :: S(nrobs,nrens)
  real(kind=8),    intent(in)  :: R(nrobs,nrobs)
  real(kind=8),    intent(out) :: W(nrobs,nrmin)
  real(kind=8),    intent(out) :: eig(nrmin)
  real(kind=8),    intent(in)  :: truncation
  
  ! local variables
  real(kind=8) :: U0(nrobs,nrmin)
  real(kind=8) :: sig0(nrmin)
  real(kind=8) :: B(nrmin,nrmin)
  real(kind=8) :: Z(nrmin,nrmin)
  integer :: i,j

  ! Compute SVD of S=HA`  ->  U0, sig0
  call  svdS(S,nrobs,nrens,nrmin,U0,sig0,truncation)

  ! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}
  call lowrankCee(B,nrmin,nrobs,nrens,R,U0,sig0)
    
  ! Compute eigenvalue decomposition  of B(nrmin,nrmin)
  call eigC(B,nrmin,Z,eig)

  ! Compute inverse diagonal of (I+Lamda)
  do i=1,nrmin
     eig(i)=1.d0/(1.d0+eig(i))
  enddo

  ! W = U0 * sig0^{-1} * Z
  do j=1,nrmin
     do i=1,nrmin
        Z(i,j)=sig0(i)*Z(i,j)
     enddo
  enddo
  
  call dgemm('n','n',nrobs,nrmin,nrmin, 1.d0,U0,nrobs, Z,nrmin, 0.d0,W,nrobs)
  
end subroutine lowrankCinv

subroutine lowrankCee(B,nrmin,nrobs,nrens,R,U0,sig0)

  implicit none

  ! arguments
  integer, intent(in) :: nrmin
  integer, intent(in) :: nrobs
  integer, intent(in) :: nrens
  real(kind=8), intent(inout) :: B(nrmin,nrmin)
  real(kind=8), intent(in)    :: R(nrobs,nrobs)
  real(kind=8), intent(in)    :: U0(nrobs,nrmin)
  real(kind=8), intent(in)    :: sig0(nrmin)

  ! local variables
  real(kind=8) :: X0(nrmin,nrobs)
  integer ::  i,j
  
  ! Compute B=sig0^{-1} U0^T R U0 sig0^{-1}

  ! X0= U0^T R
  call dgemm('t','n',nrmin,nrobs,nrobs, 1.d0,U0,nrobs, R,nrobs, 0.d0,X0,nrmin)
  
  ! B= X0 U0
  call dgemm('n','n',nrmin,nrmin,nrobs, 1.d0,X0,nrmin, U0,nrobs, 0.d0,B,nrmin)
  
  do j=1,nrmin
     do i=1,nrmin
        B(i,j)=sig0(i)*B(i,j)
     enddo
  enddo
  
  do j=1,nrmin
     do i=1,nrmin
        B(i,j)=sig0(j)*B(i,j)
     enddo
  enddo
  
  B=dble(nrens-1)*B
  
end subroutine lowrankCee

subroutine svdS(S, nrobs,nrens,nrmin,U0,sig0,truncation)

  implicit none
  
  ! arguments
  integer, intent(in)  :: nrobs
  integer, intent(in)  :: nrens
  integer, intent(in)  :: nrmin
  real(kind=8),    intent(in)  :: S(nrobs,nrens)
  real(kind=8),    intent(out) :: sig0(nrmin)
  real(kind=8),    intent(inout)  :: U0(nrobs,nrmin)
  real(kind=8),    intent(in)  :: truncation

  ! local variables
  real(kind=8) :: S0(nrobs,nrens)
  real(kind=8) :: VT0(1,1)
  integer :: ierr
  integer :: lwork
  integer :: nrsigma,i
  real(kind=8) :: sigsum,sigsum1
  real(kind=8), allocatable, dimension(:) :: work

  ! Compute SVD of S=HA`  ->  U0, sig0
  lwork=2*max(3*nrens+nrobs,5*nrens)
  allocate(work(lwork))

  ! important (need a local and variable array for dgesv)
  S0=S

  sig0=1.d0
  call dgesvd('S', 'N', nrobs, nrens, S0, nrobs, sig0, U0, nrobs, VT0, nrens, work, lwork, ierr)
  deallocate(work)
  if (ierr /= 0) then
     write(*,'(A,I6)') 'svdS: ierr from call dgesvd 0= ',ierr
     stop
  endif

  sigsum=0.d0
  do i=1,nrmin
     sigsum=sigsum+sig0(i)**2
  enddo
  
  sigsum1=0.d0
  ! Significant eigenvalues.
  nrsigma=0
  do i=1,nrmin                       
     if (sigsum1/sigsum < truncation) then
        nrsigma=nrsigma+1
        sigsum1=sigsum1+sig0(i)**2
     else
        sig0(i:nrmin)=0.d0
        exit
     endif
  enddo

!  write(*,'(a,i5,g13.5)') ' dominant sing. values and share ',nrsigma,sigsum1/sigsum
!  write(*,'(5g11.3)')sig0
  
  do i=1,nrsigma
     sig0(i) = 1.d0/sig0(i)
  enddo

end subroutine svdS

subroutine mean_preserving_rotation(Up,nrens)
  ! Generates the mean preserving random rotation for the EnKF SQRT algorithm
  ! using the algorithm from Sakov 2006-07.  I.e, generate rotation Up such that
  ! Up*Up^T=I and Up*1=1 (all rows have sum = 1)  see eq 17.
  ! From eq 18,    Up=B * Upb * B^T 
  ! B is a random orthonormal basis with the elements in the first column equals 1/sqrt(nrens)
  ! Upb = | 1  0 |
  !       | 0  U |
  ! where U is an arbitrary orthonormal matrix of dim nrens-1 x nrens-1  (eq. 19)
  
  implicit none
  
  ! arguments
  integer, intent(in)    :: nrens
  real(kind=8), intent(out) :: Up(nrens,nrens)

  ! local variables
  real(kind=8) :: B(nrens,nrens)
  real(kind=8) :: Q(nrens,nrens)
  real(kind=8) :: R(nrens,nrens)
  real(kind=8) :: U(nrens-1,nrens-1)
  real(kind=8) :: Upb(nrens,nrens)
  integer :: j,k
 
  ! Generating the B matrix
  ! Starting with a random matrix with the correct 1st column
  call random_number(B)
  B(:,1)=1.d0/sqrt(dble(nrens))

  ! with overwriting of B
  do k=1,nrens
     R(k,k)=sqrt(dot_product(B(:,k),B(:,k)))
     B(:,k)=B(:,k)/R(k,k)
     do j=k+1,nrens
        R(k,j)=dot_product(B(:,k),B(:,j))
        B(:,j)=B(:,j)- B(:,k)*R(k,j)
     enddo
  enddo

  ! Creating the orthonormal nrens-1 x nrens-1 U matrix
  call randrot(U,nrens-1)

  ! Creating the orthonormal nrens x nrens Upb matrix
  Upb(2:nrens,2:nrens)=U(1:nrens-1,1:nrens-1)
  Upb(1,1)=1.d0
  Upb(2:nrens,1)=0.d0
  Upb(1,2:nrens)=0.d0
  
  
  ! Creating the random orthonormal mean preserving nrens x nrens Upb matrix: Up=B^T Upb B
  call dgemm('n','n',nrens,nrens,nrens,1.d0,B,nrens,Upb,nrens,0.d0,Q,nrens)
  call dgemm('n','t',nrens,nrens,nrens,1.d0,Q,nrens,B,nrens,0.d0,Up,nrens)
  
end subroutine mean_preserving_rotation

end module enkf_module_r8

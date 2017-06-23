! *************************************************
! schro 1D
! *************************************************
! solve the stationnary Schrödinger equation
! for 1 dof system using finite differences
! *************************************************
! All units are atomic units
! *************************************************
! Copyright (C) 2017 Fabien Brieuc under the terms
! of the GNU General Public License.
! *************************************************

! module containing parameters
module param
integer, parameter:: dp = kind(0.d0) ! double precision floating numbers
! constants
real(dp),parameter ::TWOPI = 2.d0 * acos(-1.d0) ! 2*Pi
real(dp),parameter ::kB = 8.617d-5 / 27.211d0   ! Boltzmann constant (au)
real(dp),parameter ::hbar = 1.0d0               ! reduced Planck constant (au)
! calculation parameters
real(dp) :: m = 1.d0     ! mass (au)
real(dp) :: x0 = 1.d0    ! position of the bottom of the wells (bohr)
real(dp) :: V0 = 1.d0    ! barrier height for DW (eV)
real(dp) :: A = 1.d0     ! for quartic potential V(x)=A*x^4 (au)
real(dp) :: f0 = 10.d0   ! Frequency of the HO (THz)
real(dp) :: w0           ! Pulsation of the HO
real(dp) :: T = 300.d0   ! temperature in K
real(dp) :: Dm = 1.d0    ! Dissociation energy of the Morse
                         ! potential (in eV)
real(dp) :: alpha = 1.d0 ! width of the Morse pot. (1/bohr)
real(dp) :: xi = 2.d0    ! limit of the interval on x (bohr)
integer  :: n = 401      ! nb of pts in x = nb of eigenvalues (energies)
                         ! and eigenvectors (wave functions)
integer  :: k = 400      ! nb of states actually taken into account for average
character(len=20)  :: pot = 'harmonic' ! defines which external potential is used
                           ! pot = 'harmonic', 'double-well', 'quartic', 'morse'
namelist / parameters / n, k, xi
namelist / calculation / m, x0, V0, A, f0, T, Dm, alpha, pot
end module param

!Principal program
program Schrodinger
use param
implicit none

!Declarations
!Variables
integer  :: info,i,j
real(dp) :: dx,x,xm
real(dp),dimension(:), allocatable   :: d,e   ! for Matrix Diagonalisation (Lapack)
                                              ! d: diagonal, e: sub-diagonal
real(dp),dimension(:,:), allocatable :: vp    ! eigenvectors
real(dp),dimension(:), allocatable   :: work  ! work array for diagonalisation
real(dp),dimension(:), allocatable   :: proba ! proba density of position (1/bohr)
real(dp) :: Z=0.
!Functions
real(dp) :: V ! potential well

!read input file
open(10,file='input.dat')
read(10, calculation)
V0 = V0 / 27.211d0
f0 = f0 / 4.1341d4
w0 = TWOPI * f0
Dm = Dm / 27.211d0
read(10, parameters)
close(10)

! allocation
allocate(d(n))
allocate(e(n))
allocate(work(2*n+2))
allocate(proba(n))
allocate(vp(n,n))

! initialization
proba = 0.d0
work = 0.d0
dx = 2 * xi / dfloat(n-1)

do i = 1, n
   x = -xi + (i - 1) * dx
   d(i) =  hbar**2/(m*dx**2) + V(x)
   e(i) = -hbar**2/(2.d0*m*dx**2)
enddo
e(1) = 0.d0

do i = 1, n
   do j = 1, n
      if(i == j) then
         vp(i,j) = 1.d0
      else
         vp(i,j) = 0.d0
      endif
   enddo
enddo

! diagonalisation
call dsteqr('V',n,d,e,vp,n,work,info)

open(10,file='energies.res')
write(10,*) '# quantum nb (n), associated energy (eV), associated Boltzman coeff'
open(11,file='wavefunctions.res')
write(11,*) '# x (bohr), unormalized wave function psi(x)'
open(12,file='proba-density.res')
write(12,*) '# x (bohr), unormalized proba. density |psi(x)|^2'
open(13,file='position-proba-density.res')
write(13,*) '# x (bohr), normalized proba. density of position at temperature T'
open(14,file='potential.res')
write(14,*) '# x (bohr), potential well V(x) (eV)'

! write results in output files
do i = 1, n
   x = -xi + (i - 1) * dx
   write(10,*) i, d(i) * 27.211d0, exp(-d(i)/(kB*T)) ! energies
   write(11,*) x, vp(i,:) ! states (Wavefunctions)
   write(12,*) x, vp(i,:)**2 ! proba density of states
   do j = 1, k
      proba(i) = proba(i) + (vp(i,j)**2) * dexp(-d(j)/(kB*T))
   enddo
enddo

! normalised proba density
Z  = sum(proba) * dx
xm = 0.D0

do i = 1, n
    x = -xi + (i - 1) * dx
   write(13,*) x, proba(i) / Z
   write(14,*) x, V(x)
   ! average position
   xm = xm + x * proba(i) / Z * dx
enddo

print*,'Average position:', xm, 'bohr'

close(10)
close(11)
close(13)
close(14)
close(15)

end program Schrodinger

real(dp) function V(x)
use param
implicit none
real(dp), intent(in) :: x

if (trim(pot) == 'harmonic') then
   V=0.5d0*m*(w0**2)*x**2 !Harmonic Oscillator
else if (trim(pot) == 'double-well') then
   V = V0*((x/x0)**2-1.d0)**2 !Double well
else if (trim(pot) == 'quartic') then
   V=A*x**4  !Quartic potential
else if (trim(pot) == 'morse') then
   V=Dm*(Dexp(-alpha*x)-1.d0)**2 !Morse potenital
else
   V=0
   print*, 'error in potential function'
   stop
endif

end function V

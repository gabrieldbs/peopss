!     The routine kpsol is the preconditioner solve routine. It must have
!     that specific name be used in order that the c code can find and link
!     to it.  The argument list must also be as illustrated below:

subroutine fkpsol(udata, uscale, fdata, fscale, vv, ftem, ier)
use mkinsol
implicit none

integer neq, i
double precision udata(*), uscale(*), fdata(*), fscale(*)
double precision vv(*), ftem(*)
integer *4 ier ! Kinsol error flag

common /psize/ neq

neq = 4

do  i = 1, neq
   vv(i) = vv(i) * pp(i)
enddo
ier = 0

return
end

!* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!* * *
!c     The routine kpreco is the preconditioner setup routine. It must have
!c     that specific name be used in order that the c code can find and link
!c     to it.  The argument list must also be as illustrated below:

subroutine fkpset(udata, uscale, fdata, fscale,vtemp1,vtemp2, ier)
use mkinsol
use system
implicit none
integer *4 ier ! Kinsol error flag
integer *4 neq, i
double precision udata(*), uscale(*), fdata(*), fscale(*)
double precision vtemp1(*), vtemp2(*)

common /psize/ neq

do i = 1, neq
   pp(i) = 0.1 / (1.0+exp(-udata(i)))
enddo
   ier = 0
return
end

subroutine call_kinsol(x1_old, xg1_old, ier)
implicit none
integer *4 ier ! Kinsol error flag
integer i
real*8 x1(4), xg1(4)
real*8 x1_old(4), xg1_old(4)
integer*8 iout(15) ! Kinsol additional output information
real*8 rout(2) ! Kinsol additional out information
integer*8 msbpre
real*8 fnormtol, scsteptol
real*8 scale(4)
real*8 constr(4)
integer*4  globalstrat, maxl, maxlrst
integer neq ! Kinsol number of equations
integer*4 max_niter
common /psize/ neq ! Kinsol
integer ierr

neq=4

! INICIA KINSOL

msbpre  = 10 ! maximum number of iterations without prec. setup (?)
fnormtol = 1.0d-3 !Function-norm stopping tolerance
scsteptol = 1.0d-8 ! Function-norm stopping tolerance

maxl = 10000 ! maximum Krylov subspace dimesion (?!?!?!) ! Esto se usa para el preconditioner
maxlrst = 5 ! maximum number of restarts
max_niter = 100000
globalstrat = 0

call fnvinits(3, neq, ier) ! fnvinits inits NVECTOR module
if (ier .ne. 0) then       ! 3 for Kinsol, neq ecuantion number, ier error flag (0 is OK)
  print*, 'call_kinsol: SUNDIALS_ERROR: FNVINITS returned IER = ', ier
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
endif

call fkinmalloc(iout, rout, ier)    ! Allocates memory and output additional information
if (ier .ne. 0) then
   print*, 'call_kinsol: SUNDIALS_ERROR: FKINMALLOC returned IER = ', ier
   call MPI_FINALIZE(ierr) ! finaliza MPI
   stop
 endif

call fkinsetiin('MAX_SETUPS', msbpre, ier)  ! Additional input information
call fkinsetrin('FNORM_TOL', fnormtol, ier)
call fkinsetrin('SSTEP_TOL', scsteptol, ier)
call fkinsetiin('MAX_NITER', max_niter, ier)

do i = 1, neq  !constraint vector
   constr(i) = 1.0
enddo

call fkinsetvin('CONSTR_VEC', constr, ier) ! constraint vector
! CALL FKINSPTFQMR (MAXL, IER)
call fkinspgmr(maxl, maxlrst, ier) !  Scale Preconditioned GMRES solution of linear system (???)
!call fkinspbcg(maxl, ier) !  Scale Preconditioned BCG
!call fkindense(neq,ier)

if (ier .ne. 0) then
  print*, 'call_kinsol: SUNDIALS_ERROR: FKINSPGMR returned IER = ', ier
  call fkinfree ! libera memoria
  call MPI_FINALIZE(ierr) ! finaliza MPI
  stop
endif
call fkinspilssetprec(1, ier) ! preconditiones

do i = 1, neq ! scaling vector
  scale(i) = 1.0
enddo

do i = 1, neq ! Initial guess
      x1(i) = x1_old(i)
      xg1(i) = x1(i)  
enddo

print*,'x1',x1
call fkinsol(x1, globalstrat, scale, scale, ier)         ! Llama a kinsol

print*, 'IER:', ier

if (ier .lt. 0) then
      print*, 'call_kinsol: SUNDIALS_ERROR: FKINSOL returned IER = ', ier
      print*, 'call_kinsol: Linear Solver returned IER = ', iout(9)
!      call fkinfree
!      stop
endif

do i = 1, neq ! output
  x1_old(i) = x1(i)
  xg1_old(i) = x1(i)
enddo

call fkinfree
return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subrutina que llama a kinsol


subroutine call_fkfun(x1_old)
use MPI

integer i

real*8 x1_old(4)
real*8 x1(4)
real*8 f(4)

! MPI

integer tag
parameter(tag = 0)
integer err

x1 = 0.0
do i = 1,neq
  x1(i) = x1_old(i)
enddo

CALL MPI_BCAST(x1, neq, MPI_DOUBLE_PRECISION,0, MPI_COMM_WORLD,err)

call fkfun(x1,f, ier) ! todavia no hay solucion => fkfun 
end


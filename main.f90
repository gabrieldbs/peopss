!  ###############################################################################!     
  !     PEO/Pss Molecular Theory Program 
  !    
  !###############################################################################
!use pks
use system
use const
use solver
use results
  implicit none
  integer i, flagcrash,iii

real*8 csalcoef
  !print*, 'GIT Version: ', _VERSION
ntot=4
  call readinput
  call allocation
  vab=1.
  vpol=vpolcero/vsol

  vneg=4./3.*pi*rsal**3/vsol !volume of anion in units of vsol
  vpos=4./3.*pi*rsal**3/vsol !volume of cation in units of vsol
  yes=0 ! es para  chequear si encuentra o no xalpha, xbeta
  flagcrash=1

  Kd=10**(-pKD)
  KA=10**(-pKA)
  KEO=10**(-pKEO)

  do i = 1, ncsal ! loop in csal

  csalcoef = csalini + (csalfin-csalini)/float(ncsal-1)*float(i-1)  !Na
!  xNaalpha= 10**(-csalcoef)

  xmNaalpha=csalcoef*Na/(1.0d24)
  K0A = (KA)*(vsol*Na/1.0d24)! 	CHECKEAR
  K0EO = (KEO)/(vsol*Na/1.0d24)!   ''
  K0D = (KD)*(vsol*Na/1.0d24)**2 !!!!!!!!!!!!!!!!!!!!!!!
iter=0
  call solve(flagcrash)

!if (flaggg==1)then
!yes=yes+1
!!PAra guardar data
!arraymNa(1,yes)=xmNaalpha/Na*1.d24 
!arraymNa(2,yes)=xmNabeta/Na*1.d24 

!arraymCl(1,yes)=xmClalpha/Na*1.d24 
!arraymCl(2,yes)=xmClbeta/Na*1.d24 

!arraymcsal(1,yes)=arraymNa(1,yes)+arraymCl(1,yes)
!arraymcsal(2,yes)=arraymNa(2,yes)+arraymCl(2,yes)


!arraymA(1,yes)=xmAalpha*1.e24/Na  ! concentracion molar de monomeros, el "*2" da cuenta de que es conc. polA+polB  
!arraymA(2,yes)=xmAbeta*1.e24/Na

!arraymEO(1,yes)=xmEOalpha*1.e24/Na  ! concentracion molar de monomeros, el "*2" da cuenta de que es conc. polA+polB  
!arraymEO(2,yes)=xmEObeta*1.e24/Na

!arraympoltot(1,yes)=xmpoltotalalpha
!arraympoltot(2,yes)=xmpoltotalbeta

!arrayratioEOA(1,yes)=ratioEOAalpha 
!arrayratioEOA(2,yes)=ratioEOAbeta 

!endif


  enddo ! i



open (unit=3,file='csal_poltot_mol_alpha.txt',status='replace')

do iii=1,yes
   write (3,*) arraympoltot(1,iii), arraymcsal(1,iii)
end do

open (unit=4,file='csal_poltot_mol_beta.txt',status='replace')

do iii=1,yes
   write (4,*) arraympoltot(2,iii), arraymcsal(2,iii)
end do

open (unit=40,file='cpoltot_ratioEOA_mol_alpha.txt',status='replace')

do iii=1,yes
   write (40,*) arrayratioEOA(1,iii), arraympoltot(1,iii)
end do

open (unit=30,file='cpoltot_ratioEOA_mol_beta.txt',status='replace')

do iii=1,yes
   write (30,*) arrayratioEOA(2,iii), arraympoltot(2,iii)
end do

open (unit=400,file='csal_ratioEOA_mol_alpha.txt',status='replace')

do iii=1,yes
   write (400,*) arrayratioEOA(1,iii), arraymcsal(1,iii)
end do

open (unit=300,file='csal_ratioEOA_mol_beta.txt',status='replace')

do iii=1,yes
   write (300,*) arrayratioEOA(2,iii), arraymcsal(2,iii)
end do


!open (unit=300,file='cratiocsalalpha.txt',status='replace')

!do iii=1,yes
 !  write (300,*)  arrayratioEOA(1,iii), arraycsal(1,iii)
!end do

!open (unit=400,file='cratiocsalbeta.txt',status='replace')

!do iii=1,yes
!   write (400,*)  arrayratioEOA(2,iii), arraycsal(2,iii)
!end do


 
  call endall     ! clean up and terminate


end 
subroutine endall
stop
end subroutine

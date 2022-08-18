subroutine solve(flagcrash)

!use pks 
use system
use results
use solver
use const

implicit none
real*8 tolerancia,criterio,check_Ka_alpha,check_ka_beta,checkb_Kai,KK0check,KKaAcheckplus,kkaBcheckmin
integer nmax
real*8 x1(4)
real*8 x1g(4)
real*8 checkresults
integer ier, i,newt, j, ii, jj
integer flagcrash
integer flag
integer ngrid,iii
real*8 xiteri

   nmax=npasos     ! npasos por consola i
   criterio=1E-6!criterio pra la norma
   tolerancia=1E-6!criterio pra la diferencia de concentraciones relativa


linearsolver = 2
flaggg=0

   do i=0, npasosrhotot
 flaggg=0 
      xiteri = (rhototmax-rhototmin)*float(i)/float(npasosrhotot) + rhototmin
      xmpoltotalalpha= xiteri  !Segunda variable que fijamos  xmpoltotalalpha

      x1(1)=ratioEOAalphainitial! Ratio inicial alpha 1:1
      x1g(1)=x1(1)  !
      x1(2)=xmNabetainitial     !xNabeta  inicial
      x1g(2)=x1(2)
      x1(3)=xmpoltotalbetainitial!- xmpoltotalalpha  !xxmpoltotalbeta initicial
      x1g(3)=x1(3)
      x1(4)=ratioEOAbetainitial           !xratioeobeta inicial 
      x1g(4)=x1(4)
!      print*, 10**(-x1(1)), 10**(-x1(2))
      checkresults=0.
      call call_kinsol(x1, x1g, ier)
print*,'x1',xmNaalpha,xmpoltotalalpha,x1

      checkresults=abs(x1(1)- x1(4)/(x1(1)+x1(4))) ! ratioalpha - ratiobeta
      checkresults=checkresults+abs(x1(2)- xmNaalpha/(x1(2)+xmNaalpha)) ! rhoNabeta- rhoNaalpha
      checkresults=checkresults+abs(x1(3)- xmpoltotalalpha/(x1(3)+xmpoltotalalpha))!rhopoltotalalpha -rhopoltotalbeta
   !   if ((norma.lt.criterio).and.(checkresults.gt.tolerancia))exit
   !enddo

      

       if ((norma.lt.criterio).and.( checkresults .gt.tolerancia)) then ! encuentra solucion
         flaggg=1
         print*,'Grid Point OK',yes

      ratioEOAalphainitial=x1(1)! Ratio inicial alpha 1:1
      xmNabetainitial=x1(2)     !xNabeta  inicial
      xmpoltotalbetainitial=x1(3)!- xmpoltotalalpha  !xxmpoltotalbeta initicial
      ratioEOAbetainitial=x1(4)           !xratioeobeta inicial 
!      
     !    if(xmpoltotalalpha.gt.x1(3)) then !! definicion de alpha y beta para guardarla data ordenada
            xmpolalpha = xmpoltotalalpha ! volume fraction
            xmpolbeta = x1(3)
            xmsalalpha = xmNaalpha+xmClalpha ! number density
            xmsalbeta = xmNabeta +xmClbeta
            xratioEOAalpha=ratioEOAalpha
            xratioEOAbeta=ratioEOAbeta            
     !    else
     !       xmpolalpha = x1(3) ! volume fraction
     !       xmpolbeta = xmpoltotalalpha
     !       xmsalalpha = xmNabeta+xmClbeta ! number density
     !       xmsalbeta = xmNaalpha +xmClalpha 
        !    xratioEOAalpha=ratioEOAalpha
        !    xratioEOAbeta=ratioEOAbeta
    
      !endif                     


if (flaggg==1)then
yes=yes+1
!!PAra guardar data
arraymNa(1,yes)=xmNaalpha/Na*1.d24
arraymNa(2,yes)=xmNabeta/Na*1.d24

arraymCl(1,yes)=xmClalpha/Na*1.d24
arraymCl(2,yes)=xmClbeta/Na*1.d24

arraymcsal(1,yes)=arraymNa(1,yes)+arraymCl(1,yes)!+faltta la parte asoc.
arraymcsal(2,yes)=arraymNa(2,yes)+arraymCl(2,yes)!falta parte asoc


arraymA(1,yes)=xmAalpha/Na*1.e24  ! concentracion molar de monomeros, el "*2" da cuenta de que es conc. polA+polB  
arraymA(2,yes)=xmAbeta/Na*1.e24

arraymEO(1,yes)=xmEOalpha/Na*1.e24  ! concentracion molar de monomeros, el "*2" da cuenta de que es conc. polA+polB  
arraymEO(2,yes)=xmEObeta/Na*1.e24

!arraympoltot(1,yes)=xmpoltotalbeta
!arraympoltot(2,yes)=xmpoltotalbeta

!arrayratioEOA(1,yes)=ratioEOAalpha
!arrayratioEOA(2,yes)=ratioEOAbeta

arraympoltot(1,yes)=arraymA(1,yes)+arraymEO(1,yes)
arraympoltot(2,yes)=arraymA(2,yes)+arraymEO(2,yes)

arrayratioEOA(1,yes)=arraymEO(1,yes)/arraymA(1,yes)
arrayratioEOA(2,yes)=arraymEO(2,yes)/arraymA(2,yes)


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
!  call endall     ! clean up and terminate
close(3)
close(4)
close(30)
close(40)
close(300)
close(400)
 endif

 endif
enddo
return
end subroutine


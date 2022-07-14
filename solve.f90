subroutine solve(flagcrash)

!use pks 
use system
use results
use solver
use const

implicit none
real*16 tolerancia,criterio,check_Ka_alpha,check_ka_beta,checkb_Kai,KK0check,KKaAcheckplus,kkaBcheckmin
integer nmax
real*16 xiteri
real*16 xiterii
real*16 xiterj
real*16 xiterjj
real*8 x1(4)
real*8 x1g(4)
real*16 checkresults
integer ier, i,newt, j, ii, jj
integer flagcrash
real*16 xatest(2),xbtest(2),segsolalpha(2),segsolbeta(2),yy(2),yybet(2)
real*16 frac(8),fracbet(8)
real*16 phi0,xphisales,xsolvit
integer flag
integer ngrid
real*16 xiter

   nmax=npasos     ! npasos por consola i
   criterio=1E-8!criterio pra la norma
   tolerancia=0.01!criterio pra la diferencia de concentraciones relativa

! #### GRID SEARCH

linearsolver = 2
flaggg=0
! ida

   do i=0, npasosrhotot
  
      xiteri = (log10(rhototmax)-log10(rhototmin))*float(i)/float(npasosrhotot) + log10(rhototmin)
      xmpoltotalalpha=10**(-xiteri)  !Segunda variable que fijamos  xmpoltotalalpha

      x1(1)=1.0! Ratio inicial alpha 1:1
      x1g(1)=x1(1)  !
      x1(2)=0.01     !xNabeta  inicial
      x1g(2)=x1(2)
      x1(3)=0.01!- xmpoltotalalpha  !xxmpoltotalbeta initicial
      x1g(3)=x1(3)
      x1(4)=0.01           !xratioeobeta inicial 
      x1g(4)=x1(4)
!      print*, 10**(-x1(1)), 10**(-x1(2))
      checkresults=0.
      call call_kinsol(x1, x1g, ier)
      checkresults=abs(x1(1)- x1(4)/(x1(1)+x1(4))) ! ratioalpha - ratiobeta
      checkresults=checkresults+abs(x1(2)- xmNaalpha/(x1(2)+xmNaalpha)) ! rhoNabeta- rhoNaalpha
      checkresults=checkresults+abs(x1(3)- xmpoltotalalpha/(x1(3)+xmpoltotalalpha))!rhopoltotalalpha -rhopoltotalbeta
      if ((norma.lt.criterio).and.(checkresults.gt.tolerancia))exit
   enddo

      

       if ((norma.lt.criterio).and.( checkresults .gt.tolerancia)) then ! encuentra solucion

         flaggg=1
         print*,'Grid Point OK',yes

       
         if(xmpoltotalalpha.gt.x1(3)) then !! definicion de alpha y beta para guardarla data ordenada
            xmpolalpha = xmpoltotalalpha ! volume fraction
            xmpolbeta = x1(3)
            xmsalalpha = xmNaalpha+xmClalpha ! number density
            xmsalbeta = xmNabeta +xmClbeta
            xratioEOAalpha=ratioEOAalpha
            xratioEOAbeta=ratioEOAbeta            
         else
            xmpolalpha = x1(3) ! volume fraction
            xmpolbeta = xmpoltotalalpha
            xmsalalpha = xmNabeta+xmClbeta ! number density
            xmsalbeta = xmNaalpha +xmClalpha 
            xratioEOAalpha=ratioEOAalpha
            xratioEOAbeta=ratioEOAbeta
    
         endif                     


       endif

return
end subroutine


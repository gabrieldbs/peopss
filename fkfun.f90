subroutine fkfun(x,f,ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User provided routine for kinsol
! x is the input vector
! f is the output vector, kinsol will change x in order to get f = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use solver
use system
!use molecules
!use bulk
use const
use results
implicit none
 
real*16 xmsalt
real*16 Penality

integer*4 ier
real*16 vectfalpha(3),vectfbeta(3),vectfrac(6)
real*16 x(4),f(4)
real*16 potA,potEo,potNa,free_ener
real*16 muAalpha,muAbeta,muEOalpha,muEObeta,fealpha,febeta
real*16 potquimA,elib,potquimEO,potquimNa,muNaalpha,muNabeta
integer i

print*,'x1',x
if (iter==0)then
x(1)=0.1
x(2)=0.1
x(3)=0.01
x(4)=1
endif
! xmnaalpha , xmpoltotalpha  son inputs
ratioEOAalpha=x(1)    ! phipeo/phiA en alpha

xmNAbeta=x(2)        !rhoNa en beta
xNabeta=xmNabeta*(vpos*vsol)
xmpoltotalbeta=x(3)    ! phitot= phi peo+ phi A en beta
ratioEOAbeta=x(4)   !phipeo/phiA en beta

xmAalpha=xmpoltotalalpha/(1.+ratioEOAalpha)
xmEOalpha=xmAalpha*ratioEOAalpha
xmAbeta=xmpoltotalbeta/(1.+ratioEOAbeta)
xmEObeta=xmAbeta *ratioEOAbeta

print*,'alpha',xmNaalpha,xmpoltotalalpha,ratioEOAalpha
print*,'beta',xmNabeta,xmpoltotalbeta,ratioEOAbeta


vectfalpha(1)=xmAalpha
vectfalpha(2)=xmEOalpha
vectfalpha(3)=xmNaalpha
vectfbeta(1)=xmAbeta
vectfbeta(2)=xmEObeta
vectfbeta(3)=xmNabeta

call fractions(vectfalpha,vectfrac)
!vectfractions-> (1) fA_aspol,(2)fEO_aspol,(3) fA_asion (4)fEO_Asion,(5)fA_unas,f(6)fEO_unas
fA_aspol_alpha=vectfrac(1)
fEO_aspol_alpha=vectfrac(2)
fA_asion_alpha=vectfrac(3)
fEO_asion_alpha=vectfrac(4)
fA_unas_alpha=vectfrac(5)
fEO_unas_alpha=vectfrac(6)

!print*,'alpha',vectfrac

call fractions(vectfbeta,vectfrac)
!vectfractions-> (1) fA_aspol,(2)fEO_aspol,(3) fA_asion (4)fEO_Asion,(5)fA_unas,f(6)fEO_unas
fA_aspol_beta=vectfrac(1)
fEO_aspol_beta=vectfrac(2)
fA_asion_beta=vectfrac(3)
fEO_asion_beta=vectfrac(4)
fA_unas_beta=vectfrac(5)
fEO_unas_beta=vectfrac(6)
!print*,'beta',vectfrac
!!!!

xmClalpha=-fA_unas_alpha*Ma* xmAalpha  + fEO_asion_alpha*Meo*xmEOalpha + xmNaalpha  ! rhoCl en alpha
xmClbeta=-fA_unas_beta*Ma*xmAbeta + fEO_asion_beta*Meo*xmEObeta + xmNabeta  ! rhoCl en alpha
print*,'cl',xmclalpha,xmclbeta
xSolventalpha=1.-Ma*vpol*vsol*xmAalpha -Meo*vpol*vsol*xmEOalpha&
-vpos*vsol*((fA_asion_alpha+fA_aspol_alpha)*Ma*xmAalpha+fEO_asion_alpha*Meo*xmEOalpha)&
-xmNaalpha *vpos*vsol -xmClalpha*vneg*vsol

xSolventbeta=1.- Ma*vpol*vsol*xmAbeta -Meo*vpol*vsol*xmEObeta&
-vpos*vsol*((fA_asion_beta+fA_aspol_beta)*Ma*xmAbeta+fEO_asion_beta*Meo*xmEObeta)&
-xmNabeta*vpos*vsol -xmClbeta*vneg*vsol

xmSolventalpha=xSolventalpha/vsol
xmSolventbeta =xSolventbeta/vsol
print*,'solv',xmsolventalpha,xmsolventbeta

packconst=(1./vsol)*(log(xSolventalpha)-log(xSolventbeta) ) ! para no poner pi 
neutralconst=log(xmClbeta*vsol)+vneg*vsol*packconst-log(xmClalpha*vsol) ! phi
print*,'const',packconst,neutralconst

potquimNa=0.
call muNa(potquimNa)
potNa=potquimNa
print*,'pot',potNA

potquimA=0.
call muA(potquimA)
potA = potquimA


potquimEO=0.
call muEO(potquimEO)
potEO = potquimEO

elib=0.
call fe(elib)
free_ener=elib

Penality=((xmNabeta-xmNaalpha)**2+(xmAalpha-xmAbeta)**2+(xmEOalpha-xmEObeta)**2)

 f(1)=0*free_ener/Penality

 f(2)=0*potNa/Penality

 f(3)=0*potA/Penality

 f(4)=0*potEO/Penality

iter = iter + 1
norma = 0.0

do i = 1, 4
norma = norma +(f(i))**2    
enddo
print*, xNaalpha,xmpoltotalalpha,x(1)
print*, x(2), x(3), x(4)
print*, norma,penality,potA

stop 
ier = 0.0
return
end subroutine

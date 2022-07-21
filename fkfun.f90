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
 
real*8 xmsalt
real*8 Penality,testKa,testkeo,testkd

integer*4 ier
real*8 vectfalpha(3),vectfbeta(3),vectfrac(6)
real*8 x(4),f(4)
real*8 potA,potEo,potNa,free_ener
real*8 muAalpha,muAbeta,muEOalpha,muEObeta,fealpha,febeta
real*8 potquimA,elib,potquimEO,potquimNa,muNaalpha,muNabeta
integer i

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


vectfalpha(1)=xmAalpha
vectfalpha(2)=xmEOalpha
vectfalpha(3)=xmNaalpha
vectfbeta(1)=xmAbeta
vectfbeta(2)=xmEObeta
vectfbeta(3)=xmNabeta

call fractions(vectfalpha,vectfrac)
!vectfractions-> (1) fA_aspol,(2)fEO_aspol,(3) fA_asion (4)fEO_Asion,(5)fA_unas,f(6)fEO_unas
fEO_aspol_alpha=vectfrac(1)
fA_aspol_alpha=vectfrac(2)
fEO_asion_alpha=vectfrac(3)
fA_asion_alpha=vectfrac(4)
fEO_unas_alpha=vectfrac(5)
fA_unas_alpha=vectfrac(6)

!testeo fracciones
!testKa=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fA_asion_alpha-fA_aspol_alpha)/fA_asion_alpha )-pKa
!testkeo=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fEO_asion_alpha-fEO_aspol_alpha)/fEO_asion_alpha )-pKEO
!testkd=-log10((Na/1.0d24)*xmnaalpha*vsol*fA_unas_alpha*Ma*xmAalpha*vab*(1.-fEo_asion_alpha-fEO_aspol_alpha)/fEO_aspol_alpha) -pKd
!print*,'testalpha',testKa,testkeo,testkd
!stop
!print*,'alpha',vectfrac
call fractions(vectfbeta,vectfrac)
!vectfractions-> (1) fA_aspol,(2)fEO_aspol,(3) fA_asion (4)fEO_Asion,(5)fA_unas,f(6)fEO_unas
fEO_aspol_beta=vectfrac(1)
fA_aspol_beta=vectfrac(2)
fEO_asion_beta=vectfrac(3)
fA_asion_beta=vectfrac(4)
fEo_unas_beta=vectfrac(5)
fA_unas_beta=vectfrac(6)
!testeo fracciones
!testKa=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fA_asion_alpha-fA_aspol_alpha)/fA_asion_alpha )-pKa
!testkeo=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-feO_asion_alpha-fEo_aspol_alpha)/feO_asion_alpha )-pKEO
!testkd=-log10((Na/1.0d24)*xmNabeta*vsol*fA_unas_beta*Ma*xmAbeta*vab*(1-fEo_asion_beta-fEO_aspol_beta)/fEO_aspol_beta) -pKd
!print*,'testbeta',testKa,testkeo,testkd
!stop
!print*,'beta',vectfrac
!!!!

xmClalpha=-fA_unas_alpha*Ma*xmAalpha + fEO_asion_alpha*Meo*xmEOalpha + xmNaalpha  ! rhoCl en alpha
xmClbeta=-fA_unas_beta*Ma*xmAbeta + fEO_asion_beta*Meo*xmEObeta + xmNabeta  ! rhoCl en alpha
if (xmClbeta.le.0)then
xmClbeta=1E-10
print*,'Cl beta problem'

endif
if (xmClalpha.le.0)then
xmClalpha=1E-10
print*,'Cl alpha problem'
endif

!print*,'cl',xmclalpha,xmclbeta
xSolventalpha=1. -Ma*vpol*vsol*xmAalpha -Meo*vpol*vsol*xmEOalpha&
-vpos*vsol*((fA_asion_alpha+fA_aspol_alpha)*Ma*xmAalpha+fEO_asion_alpha*Meo*xmEOalpha)&
-xmNaalpha*vpos*vsol -xmClalpha*vneg*vsol

xSolventbeta=1. - Ma*vpol*vsol*xmAbeta -Meo*vpol*vsol*xmEObeta&
-vpos*vsol*((fA_asion_beta+fA_aspol_beta)*Ma*xmAbeta+fEO_asion_beta*Meo*xmEObeta)&
-xmNabeta*vpos*vsol -xmClbeta*vneg*vsol

if (xSolventbeta.le.0)then
xsolventbeta=1E-10
print*,'Solvent beta problem'
endif
if (xsolventalpha.le.0)then
xsolventalpha=1E-10
print*,'Solvent alpha problem'
endif

xmSolventalpha=xSolventalpha/vsol
xmSolventbeta =xSolventbeta/vsol
!print*,'solv',xmsolventalpha,xmsolventbeta

packconst=(1./vsol)*(log(xSolventalpha)-log(xSolventbeta) ) ! para no poner pi 
neutralconst=log(xmClbeta*vsol)+vneg*vsol*packconst-log(xmClalpha*vsol) ! phi
!print*,'const',packconst,neutralconst

!stop
potquimNa=0.
call muNa(potquimNa)
potNa=potquimNa

call muA(potquimA)
potA = potquimA


potquimEO=0.
call muEO(potquimEO)
potEO = potquimEO

elib=0.
call fe(elib)
free_ener=elib

penality=1
Penality=(xmNabeta-xmNaalpha)**2+(xmpoltotalalpha-xmpoltotalbeta)**2
Penality=Penality+(ratioEOAalpha-ratioEOAbeta)**2
penality=penality/(xmNaalpha**2+xmpoltotalalpha**2+ratioeoaalpha**2)
!penality=penality**0.5
 f(1)=free_ener/Penality

 f(2)=potNa/Penality

 f(3)=potA/Penality

 f(4)=potEO/Penality

iter = iter + 1
norma = 0.0

do i = 1, 4
norma = norma +(f(i))**2    
enddo
print*,'Norma , penality', norma,penality
!if (Norma.lt.1E-5)then
!print*,'xmNaalpha,xpoltotalAlpha',xmNaalpha,xmpoltotalAlpha
!print*,'x',x
!endif
ier = 0.0
return
end subroutine

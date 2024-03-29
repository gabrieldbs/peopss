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
real*8 testneuta,testneutb,testpcka,testpckb
integer*4 ier
real*8 vectfalpha(3),vectfbeta(3),vectfrac(6)
real*8 x(4),f(4)
real*8 potA,potEo,potNa,free_ener
real*8 muAalpha,muAbeta,muEOalpha,muEObeta,fealpha,febeta
real*8 potquimA,elib,potquimEO,potquimNa,muNaalpha,muNabeta
real*8 penalityA,penalityNa,penalityEO
integer i

! xmnaalpha , xmpoltotalpha  son inputs
!print*,'x',x
if(x(1).lt.0)then
x(1)=1E-15
print*,'error en x(1)!!!'
!stop
endif

if(x(2).lt.0)then
x(2)=1E-15
print*,'error en x(2)!!!'
!stop
endif

if(x(3).lt.0)then
x(3)=1E-15
print*,'error en x(3)!!!!'
!stop
endif

ratioEOAalpha=x(1)    !rhopeo/rhoA en alpha

xmNAbeta=x(2)        !rhoNa en beta
xNabeta=xmNabeta*(vpos*vsol)
xmpoltotalbeta=x(3)    ! phitot= phi peo+ phi A en beta
ratioEOAbeta=x(4)   !rhopeo/rhoA en beta

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
!print*,'vectalpha',vectfalpha,ratioEoAalpha,xmpoltotalalpha
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
if (xmClalpha.lt.0.0 .or. xmClbeta.lt.0.0)then
print*,'test  Cl fail',xmclalpha,xmclbeta
stop
endif
!print*,'cl',xmclalpha,xmclbeta,xmNabeta,xmAbeta,xmEObeta
xSolventalpha=1. -Ma*vpol*vsol*xmAalpha -Meo*vpol*vsol*xmEOalpha&
-vpos*vsol*((fA_asion_alpha+fA_aspol_alpha)*Ma*xmAalpha+fEO_asion_alpha*Meo*xmEOalpha)&
-xmNaalpha*vpos*vsol -xmClalpha*vneg*vsol

xSolventbeta=1. - Ma*vpol*vsol*xmAbeta -Meo*vpol*vsol*xmEObeta&
-vpos*vsol*((fA_asion_beta+fA_aspol_beta)*Ma*xmAbeta+fEO_asion_beta*Meo*xmEObeta)&
-xmNabeta*vpos*vsol -xmClbeta*vneg*vsol


xmSolventalpha=xSolventalpha/vsol
xmSolventbeta =xSolventbeta/vsol
!print*,'solv',xmsolventalpha,xmsolventbeta
if (xmsolventalpha.lt.0.0 .or. xmsolventbeta.lt.0.0)then
print*,'test  Solvent fail',xmsolventalpha,xmsolventbeta
stop
endif
packconst=(1./vsol)*(log(xSolventalpha)-log(xSolventbeta) ) ! betapi 
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
Penality=abs(xmNabeta-xmNaalpha)/(xmNabeta*0.5+xmNaalpha*0.5)
Penality=Penality+abs(xmpoltotalalpha-xmpoltotalbeta)/(xmpoltotalalpha*0.5+xmpoltotalbeta*0.5)
Penality=Penality+abs(ratioEOAalpha-ratioEOAbeta)/(ratioeoaalpha*0.5+ratioeoabeta*0.5)
!penality=penality/(xmNaalpha +xmpoltotalalpha+ratioeoAalpha)
 f(1)=free_ener/Penality

 f(2)=potNa/Penality

 f(3)=potA/Penality

 f(4)=potEO/Penality

iter = iter + 1
norma = 0.0

do i = 1, 4
!print*,'f',f(i)
norma = norma +(f(i))**2    
enddo

!testneuta=xmNaalpha -xmClalpha-Ma*xmAalpha*fA_unas_alpha+fEO_asion_alpha*Meo*xmEOalpha
!testneutb=xmNabeta -xmClbeta-Ma*xmAbeta*fA_unas_beta+fEO_asion_beta*Meo*xmEObeta

!testpcka=1.-Ma*xmAalpha*vpol*vsol-xmSolventalpha*vsol -Meo*xmEOalpha*vpol*vsol-xmNaalpha*vpos*vsol-xmClalpha*vneg*vsol
!testpcka=testpcka -vneg*vsol*(Ma*xmAalpha*(fA_asion_alpha+fA_aspol_alpha)+Meo*xmEOalpha*fEO_asion_alpha)
 
!testpckb=1.-Ma*xmAbeta*vpol*vsol-xmSolventbeta*vsol -Meo*xmEObeta*vpol*vsol-xmNabeta*vpos*vsol-xmClbeta*vneg*vsol
!testpckb=testpckb -vneg*vsol*(Ma*xmAbeta*(fA_asion_beta+fA_aspol_beta)+Meo*xmEObeta*fEO_asion_beta)
!print*,'testneutra',testneuta,testneutb,testpcka,testpckb


print*,'Norma , penality', norma,penality !, f(1), f(2), f(3) , f(4) !if (Norma.lt.1E-5)then
!print*,'xmNaalpha,xpoltotalAlpha',xmNaalpha,xmpoltotalAlpha
!print*,'x',x
!endif
ier = 0.0
return
end subroutine

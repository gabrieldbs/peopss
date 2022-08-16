subroutine fractions(vectf,vectfrac)
use system
use const
 implicit none
integer i
real*8 vectf(3),vectfrac(6)
real*8 xphiA,xphiEO
real*8 xmphiA,xmphiEO,xmphiNa
real*8 aa, bb, cc,auxA,auxB

xmphiA=vectf(1)
xmphiEo=vectf(2)
xmphiNa=vectf(3) 

vectfrac=0.;


  auxA = Meo*xmphiEO/(Ma*xmphiA)
  auxB = xmphiNa*vsol*Ma*xmphiA*vab/K0D
  auxB = auxB*K0A/(K0A+xmphiNa*vsol)
  auxB = AuxB*K0EO/(K0EO+xmphiNa*vsol)

  aa=auxA
  bb=-1.-auxA-1./auxB
  cc=1.
!print*,'aa',auxb,xmphiNA,xmphieo,xmphia 

  vectfrac(1) = (-bb - SQRT(bb**2 - 4.0*aa*cc))/(2.0*auxA) ! fEO_aspol
  vectfrac(2) = vectfrac(1)*auxA                        ! fA_aspol

if (abs(vectfrac(1)).gt.1. .or. abs(vectfrac(2)).gt.1.)then
print*,'frac',vectfrac(1),vectfrac(2),xmphiA,xmphiEo,xmphiNa
stop
endif

  vectfrac(3) = (1.0 - vectfrac(1))*(xmphiNA*vsol/(K0EO+xmphiNa*vsol)) ! fEO_asion
  vectfrac(4) = (1.0 - vectfrac(2))*(xmphiNA*vsol/(K0A +xmphiNa*vsol)) ! fa_asion

if (abs(vectfrac(3)).gt.1. .or. abs(vectfrac(4)).gt.1.)then
print*,'frac',vectfrac(3),vectfrac(4)
stop
endif

  vectfrac(5) = (1.0 - vectfrac(1)-vectfrac(3)) ! fEO_unas
  vectfrac(6) = (1.0 - vectfrac(2)-vectfrac(4)) ! fA_unas
if (abs(vectfrac(5)).gt.1. .or. abs(vectfrac(6)).gt.1.)then
print*,'frac',vectfrac(5),vectfrac(6)
stop
endif

end subroutine

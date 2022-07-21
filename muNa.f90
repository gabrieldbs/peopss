
subroutine muNa(potquimNa)
use const
use system
use results
implicit none
real*8 potquimNa
potquimNA=0
!xmNaalpha=vectNa(1)/(vpos*vsol)
!xmNabeta=vectNa(2)/(vpos*vsol)
!packconst=vectNa(3)
!neutralconst=vectNa(4)

potquimNa= log(xmNaalpha*vsol) +packconst +neutralconst -log(xmNabeta*vsol)

end  subroutine 

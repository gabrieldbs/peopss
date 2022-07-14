
subroutine muEO(potquimEO)
use const
use system
use results
implicit none
real*16 potquimEO

!xmAalpha=vectpotEO(1)
!xmEOalpha=vectpotEO(2)
!xmAbeta=vectpotEO(3) 
!xmEObeta=vectpotEO(4)
!fEo_unas_alpha=vectpotEO(5)
!fEO_unas_beta=vectpotEO(6)
!packconst=vectpotEO(7)

potquimEO=log(xmEOalpha*vsol)-log(xmEObeta*vsol)-chi*MeO*(Ma*(xmAalpha-xmAbeta)&
+Meo*(xmEOalpha-xmEObeta))+Meo*(log(fEO_unas_alpha)-log(fEO_unas_beta))-packconst*Meo*vpol*vsol

end subroutine 

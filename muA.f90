
subroutine muA(potquimA)
use const
use system
use results
implicit none
real*8 potquimA

!xmAalpha=vectpotA(1)
!xmEOalpha=vectpotA(2)
!xmAbeta=vectpotA(3) 
!xmEObeta=vectpotA(4) 
!fA_unas_alpha=vectpotA(5)
!fA_unas_beta=vectpotA(6)
!packconst=vectpotA(7)
!neutralconst=vectpoA(8)
if (fa_unas_alpha.lt.1E-15)then
fa_unas_alpha=1E-20
endif
if (fa_unas_beta.lt.1E-15)then
fa_unas_beta=1E-15
endif


potquimA= log(xmAalpha*vsol)-log(xmAbeta*vsol)- chi*MA*(MA*(xmAalpha-xmAbeta)&
+MEO*(xmEOalpha-xmEObeta))+MA*(log(fA_unas_alpha)-log(fA_unas_beta))&
-packconst*Ma*vpol*vsol+Ma*neutralconst
!print* ,potquimA
!stop
end  subroutine 

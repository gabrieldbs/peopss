subroutine fe(elib)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use results
use system
use const
implicit none

real*16 Free_energy,elib 

! Calculation of xsolvent
!xmAalpha=vectfe(1)
!xmEOalpha=vectfe(2) 
!xmSolventalpha=vectfe(3)
!xmNaalpha=vectfe(4) 
!xmClalpha=vectfe(5) 
!xmAbeta=vectfe(6) 
!xmEobeta=vectfe(7) 
!xmSolventbeta=vectfe(8) 
!xmNabeta=vectfe(9) 
!xmClbeta=vectfe(10)

!fA_aspol_alpha=vectfe(11)
!fA_aspol_beta=vectfe(12)

!packconst=vectfe(13)

Free_Energy = 0.0

! 1. solvent entropy

  Free_Energy=Free_Energy-xmAalpha-xmEOalpha-xmSolventalpha-xmNaalpha-xmClalpha

  Free_Energy=Free_Energy + 0.5*chi*(xmAalpha*Ma+xmEOalpha*MEO)**2+xmAalpha*Ma*fA_aspol_alpha

  Free_Energy=Free_Energy+xmAbeta+xmEObeta+xmSolventbeta+xmNabeta+xmClbeta

  Free_Energy=Free_Energy+packconst-0.5*chi*(Ma*xmAbeta+Meo*xmEObeta)**2-xmAbeta*Ma*fA_aspol_beta

elib=Free_Energy

return 
end subroutine


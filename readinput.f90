subroutine readinput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  This routine reads variables from fort.8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use pks
use system
use solver
!use kai

implicit none
integer i
character basura

! read starts here, not that read is performed sequentially! 

read(8,*), basura
read(8,*), Ma    ! Ma  polA

read(8, *), basura
read(8, *), MEO  !   MEO  polEO

read(8, *), basura! vp
read(8, *), vpolcero !

read(8, *), basura! rsal
read(8, *), rsal !

read(8, *), basura
read(8, *), rhototmin, rhototmax  !  rhototal maximo y minimo diluida

read(8, *), basura
read(8, *), npasosrhotot ! paso en el ntotal diluida

read(8, *), basura
read(8, *), csalini,csalfin,ncsal  ! salt concentration in bulk (Molar) diluida
!
read(8, *), basura
read(8, *), pKD     ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura
read(8, *), pKA     ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura
read(8, *), pKEO     ! polymer-polymer attraction strenght in kBT

read(8, *), basura
read(8, *), chi  ! cutoff for porr sv interaction in lattice sites

read(8, *), basura
read(8, *), ratioEOAalphainitial !diluida 

read(8, *), basura
read(8, *), xmNabetainitial,xmpoltotalbetainitial,ratioEOAbetainitial ! concentrado

!read(8, *), basura
!read(8, *), infile ! read input from file:?

end subroutine


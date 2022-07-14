!  ###############################################################################!     
  !     PEO/Pss Molecular Theory Program 
  !    
  !###############################################################################
!use pks
use system
use const
use solver
use results
  implicit none
  integer i, flagcrash,iii

  real*16 pKw,vsal,KAinput,Kbind,KBinput,KANainput,KBClinput,kte
  real*16 aa, bb, cc, xNaClbulk
  real*16 xposbulk,xnegbulk,xHplusbulk,xOHminbulk,Kw,xsalt
  real*16 xsolvit,phialphamol,phibetamol
  real*16 xsolvalpha,xsolvbeta,saleihbetamol
real*16 csalcoef
  !print*, 'GIT Version: ', _VERSION
iter=0
  call readinput
  vab=1.
  vpol=vpol/vsol

  vneg=4./3.*pi*rsal**3/vsol !volume of anion in units of vsol
  vpos=4./3.*pi*rsal**3/vsol !volume of cation in units of vsol
  yes=0 ! es para  chequear si encuentra o no xalpha, xbeta
  flagcrash=1

  Kd=10**(-pKD)
  KA=10**(-pKA)
  KEO=10**(-pKEO)

  do i = 1, ncsal ! loop in csal

  csalcoef = csalini + (csalfin-csalini)/float(ncsal-1)*float(i-1)  !Na

  xNaalpha= 10**(-csalcoef)
  xmNaalpha=xNaalpha/(vpos*vsol)
print*,'xmnaalpha',xmnaalpha
  K0A = (KA/vsol)/(Na/1.0d24)! 	CHECKEAR
  K0EO = (KEO/vsol)/(Na/1.0d24)!   ''
  K0D = (KD)*(1.0d24/Na) !!!!!!!!!!!!!!!!!!!!!!!
  call solve(flagcrash)
  
  flaggg=0
  

if (flaggg==1)then
yes=yes+1
!!PAra guardar data
!arrayNa(1,yes)=xmNaalpha/Na*1.d24 
!arrayNa(2,yes)=xmNabeta/Na*1.d24 

!arrayCl(1,yes)=xmClalpha/Na*1.d24 
!arrayCl(2,yes)=xmClbeta/Na*1.d24 

!arraycsal(1,yes)=arrayNa(1,yes)+arrayCl(1,yes)
!arraycsal(2,yes)=arrayNa(2,yes)+arrayCl(2,yes)


!arrayA(1,yes)=xAalpha*1.e24/Na/(MA*vpol)  ! concentracion molar de monomeros, el "*2" da cuenta de que es conc. polA+polB  
!arrayA(2,yes)=xAbeta*1.e24/Na/(Ma*vpol)

!arrayEO(1,yes)=xEOalpha*1.e24/Na/(Meo*vpol)  ! concentracion molar de monomeros, el "*2" da cuenta de que es conc. polA+polB  
!arrayEO(2,yes)=xEObeta*1.e24/Na/(MEO*vpol)

!arrayrhotot(1,yes)=arrayA(1,yes)+arrayEO(1,yes)
!arrayrhotot(2,yes)=arrayA(2,yes)+arrayEO(2,yes)

!arrayratioEOA(1,yes)=arrayEO(1,yes)arrayA(1,yes)
!arrayratioEOA(2,yes)=arrayEO(2,yes)/arrayA(2,yes)

!arrayfA_asion(1,yes)=fA_asion_alpha
!arrayfA_asion(2,yes)=fA_asion_beta

!arrayfEO_asion(1,yes)=fEO_asion_alpha
!arrayfEO_asion(2,yes)=fEO_asion_beta

!arrayfA_aspol(1,yes)=fA_aspol_alpha
!arrayfA_aspol(2,yes)=fA_aspol_beta

!arrayfEO_aspol(1,yes)=fEO_aspol_alpha
!arrayfEO_aspol(2,yes)=fEO_aspol_beta

!arrayfA_unas(1,yes)=fA_unas_alpha
!arrayfA_unas(2,yes)=fA_unas_beta

!arrayfEO_unas(1,yes)=fEO_unas_alpha
!arrayfEO_unas(2,yes)=fEO_unas_beta

endif


  enddo ! i


!open (unit=40,file='crhototratiobeta.txt',status='replace')

!do iii=1,yes
!   write (40,*) arrayrhotot(2,iii), arrayratioEOA(2,iii)
!end do

!open (unit=300,file='cratiocsalalpha.txt',status='replace')

!do iii=1,yes
 !  write (300,*)  arrayratioEOA(1,iii), arraycsal(1,iii)
!end do

!open (unit=400,file='cratiocsalbeta.txt',status='replace')

!do iii=1,yes
!   write (400,*)  arrayratioEOA(2,iii), arraycsal(2,iii)
!end do


 
  call endall     ! clean up and terminate


end 
subroutine endall
stop
end subroutine

module system
!real*16 delta   ! delta is the discretization lenght in z direction
!real*16 sigmaA
integer Ma,Meo
real*8 rsal
integer ntot
real*8 vpolcero
real*8 ratioEOaalphainitial
real*8 xmNabetainitial,xmpoltotalbetainitial,ratioEOAbetainitial
real*8 vpol,vpos,vneg,vab
integer npasos
integer npasosgrid
real*8 phimin, phimax
real*8 csalt
real*8 chi
integer yes
real*8 pKA,pKEo,pkD
real*8 K0A, K0EO,K0D,Ka,Keo,Kd
real*8 rhototmin,rhototmax
integer npasosrhotot
real*8 csalini,csalfin
integer ncsal
real*8 cNaplus,cClmin
real*8 xNaalpha,xNabeta,xmNaalpha,xmNabeta
real*8 xClalpha,xClbeta,xmClalpha,xmClbeta
real*8 xsolventalpha,xsolventbeta
real*8 xmsolventalpha,xmsolventbeta
real*8 xAalpha,xAbeta,xEoalpha,xEobeta
real*8 xmAalpha,xmAbeta,xmEoalpha,xmEobeta
real*8 xmpoltotalalpha,xmpoltotalbeta
real*8 ratioEOAalpha,ratioEOAbeta 
real*8 xratioEOAalpha,xratioEOAbeta,xmpolalpha,xmpolbeta
real*8 xmsalalpha,xmsalbeta
endmodule

module results
integer conteo, flaggg
real*8 fA_aspol_alpha,fA_asion_alpha,fA_unas_alpha
real*8 fEO_aspol_alpha,fEO_asion_alpha,fEo_unas_alpha
real*8 fA_aspol_beta,fA_asion_beta,fA_unas_beta
real*8 fEO_aspol_beta,fEO_asion_beta,fEo_unas_beta
real*8 packconst,neutralconst
real*8 arraympoltot(2,100000)
real*8 arraymcsal(2,100000)
real*8 arraymNa(2,100000)
real*8 arraymCl(2,100000)
real*8 arraymA(2,100000)
real*8 arraymEO(2,100000)
real*8 arrayratioEOA(2,100000)

integer cont
endmodule

module const
real*8, parameter :: pi = 3.14159 ! pi 
real*8, parameter :: Na = 6.02d23 ! Avogadro's number
real*8, parameter :: vsol = 0.03  ! bjerrum lenght in water in nm
!real*16 constq
!real*16 pKw
endmodule


module solver
real*8 norma
integer iter
integer linearsolver
endmodule


module mkinsol
double precision, allocatable :: pp(:)
endmodule

module system
!real*16 delta   ! delta is the discretization lenght in z direction
!real*16 sigmaA
real*16 Ma,Meo
real*16 rsal
real*16 vpol,vpos,vneg,vab
integer npasos
integer npasosgrid
real*16 phimin, phimax
real*16 csalt
real*16 chi
integer yes
real*16 pKA,pKEo,pkD
real*16 K0A, K0EO,K0D,Ka,Keo,Kd
real*16 rhototmin,rhototmax
integer npasosrhotot
real*16 csalini,csalfin
integer ncsal
real*16 cNaplus,cClmin
real*16 xNaalpha,xNabeta,xmNaalpha,xmNabeta
real*16 xClalpha,xClbeta,xmClalpha,xmClbeta
real*16 xsolventalpha,xsolventbeta
real*16 xmsolventalpha,xmsolventbeta
real*16 xAalpha,xAbeta,xEoalpha,xEobeta
real*16  xmAalpha,xmAbeta,xmEoalpha,xmEobeta
real*16 xmpoltotalalpha,xmpoltotalbeta
real*16 ratioEOAalpha,ratioEOAbeta 
real*16 xratioEOAalpha,xratioEOAbeta,xmpolalpha,xmpolbeta
real*16 xmsalalpha,xmsalbeta
real*16 expmupos, expmuneg 
endmodule

module results
integer conteo, flaggg
real*16 fA_aspol_alpha,fA_asion_alpha,fA_unas_alpha
real*16 fEO_aspol_alpha,fEO_asion_alpha,fEo_unas_alpha
real*16 fA_aspol_beta,fA_asion_beta,fA_unas_beta
real*16 fEO_aspol_beta,fEO_asion_beta,fEo_unas_beta
real*16 packconst,neutralconst

real*16 arrayNa(2,400000)
real*16 arrayCl(2,400000)
real*16 arraycsal(2,400000)
real*16 arrayA(2,400000)
real*16 arrayEO(2,400000)
real*16 arrayrhotot(2,400000)
real*16 arrayratioEOA(2,400000)
real*16 arraysolv(2,400000)

real*16 arrayfA_aspol(2,400000)
real*16 arrayfA_asion(2,400000)
real*16 arrayfA_unas(2,400000)
real*16 arrayfEO_aspol(2,400000)
real*16 arrayfEO_asion(2,400000)
real*16 arrayfEO_unas(2,400000)

integer cont
endmodule

module const
real*16, parameter :: pi = 3.14159 ! pi 
real*16, parameter :: Na = 6.02d23 ! Avogadro's number
real*16, parameter :: vsol = 0.03  ! bjerrum lenght in water in nm
!real*16 constq
!real*16 pKw
endmodule


module solver
real*16 norma
integer iter
integer linearsolver
endmodule


module mkinsol
double precision, allocatable :: pp(:)
endmodule

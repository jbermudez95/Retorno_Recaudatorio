/*
nombre: multiplicadores.do
Descripción: realiza estimaciones para el retorno recaudatorio del gasto público
			 en Honduras utilizando un SVAR estructural y proyecciones locales.
Autores: Jose Carlo Bermúdez y Roldan Enamorado.
Fecha: Septiembre 2020
*/

clear all

global path ""	//Inserte la dirección de su computadora donde se encuentra localizada la base de datos.
global out " "  //Inserte la dirección de su computadora donde desee guardar las estimaciones.
cd "$path"

use "data_multi.dta", replace

keep gasto_inversion_real_sa ing_trib_real_sa pib_real_sa transfer_corr_real_sa

rename gasto_inversion_real_sa exp_inv
rename ing_trib_real_sa		   ing_trib 
rename pib_real_sa			   pib 
rename transfer_corr_real_sa   transfer 

label var exp_inv	"Gasto en inversión real y desestacionalizada"
label var ing_trib	"Recaudación real y desestacionalizada"
label var pib 		"PIB real desestacionalizado"
label var transfer	"Transferencias corrientes reales y desestacionalizadas"

tsmktim date, start(2000q1) 			//Variables de tiempo y dummies estacionales
	tsset date
	generate q=quarter(dofq(date))
	tabulate q, generate(dum)
	gen t  = _n
	gen tq = (t)^2
	gen d1 = 1
	
tsfilter hp pib_hp = pib , smooth(1600) trend(trendhp_pib)
foreach var of varlist ing_trib exp_inv transfer pib  {
	g `var'_ = `var' / trendhp_pib
}
	
keep if t < 82										//Muestra hasta el primer trimestre 2020 

forval q = 1/5 {
	g du`q' = 0
		replace du`q' = 1 if t == (35 + `q')		//Dummies por la crisis financiera
}

forval j = 6/9 {
	g du`j' = 0
		replace du`j' = 1 if t == (51 + `j')		//Dummies por la reforma tributaria
}

local z "exp_inv ing_trib pib transfer"
	foreach k of local z {
        gen l`k' =ln(`k')
}

/* 
EL SVAR A ESTIMAR PARA EL GASTO CORRIENTE Y EN INVERSIÓN, RESPECTIVAMENTE

G(t) = a1 + a2*G(t-p)  + a3*T(t-p)  + a4*Y(t-p)  + e_s(t)
T(t) = a5 + a6*G(t-p)  + a7*T(t-p)  + a8*Y(t-p)	 + e_y(t)
Y(t) = a9 + a10*G(t-p) + a11*T(t-p) + a12*Y(t-p) + e_t(t)

Lueo, los choques estructurales son capturados a la Blanchard y Perotti (2002) a partir de:
| 1   0  -b1 | | e_s(t) |   | 1  b2 0 | | u_s(t) |
| 0   1  -a1 | | e_y(t) | = | a2  1 0 | | u_t(t) |
| -c2 -c1  1 | | e_t(t) |   | 0   0 1 | | u_y(t) |

En este caso, se asume b1=b2=0; c1, c2 son estimados mediante two stage least squares; a2 se imputa conforme a la boyanza del MMFMP; a2 se estima dentro del SVAR mediante factorización de Cholesky
*/

gen	   ratio_inv_tax = exp_inv / ing_trib
gen    ratio_tra_tax = transfer/ ing_trib
egen   inv_tax 		 = mean(ratio_inv_tax)
egen   tra_tax 		 = mean(ratio_tra_tax)
scalar a1 	   		 = 1.1						//Boyanza de impuestos a PIB declarada dentro del MMFMP

**********************************************************************************
*****************CHOQUE DE GASTO EN INVERSIÓN -VAR1-******************************
**********************************************************************************

qui var lexp_inv ling_trib lpib, lags(2) exog(t tq d1 du1 du2 du3 du4 du5 du6 du7 du8 du9 dum1 dum2 dum3) 
*varsoc				//Information criteria sugiere 2 lags 				
predict es_var1, residuals equation(#1)
predict et_var1, residuals equation(#2)
predict ey_var1, residuals equation(#3)

gen     et_i = et_var1 - (a1 * ey_var1)					//Ajuste cíclico a la recaudación conforme a B y P (2002)

qui ivregress 2sls ey_var1 et_var1 (es_var1 = et_i)		//Estimación de c1 y c2 mediante variables instrumentales
scalar c1_i = _b[et_var1]
scalar c2_i = _b[es_var1]

matrix A_i = (1,0,0\0,1,-a1\-c2_i,-c1_i,1)				//Matriz de choques con identificación a la B y P (2002)
matrix B_i = (1,0,0\.,1,0\0,0,1)

qui svar lexp_inv ling_trib lpib, lags(1/2) aeq(A_i) beq(B_i) exog(t tq d1 du1 du2 du3 du4 du5 du6 du7 du8 du9)
	irf create svar1, step(8) set(var1) replace
	irf table fevd,  irf(svar1) impulse(lexp_inv) response(ling_trib)
	irf table coirf, irf(svar1) impulse(lexp_inv) response(ling_trib) std
	irf table coirf, irf(svar1) impulse(lexp_inv) response(lexp_inv) std
	irf graph coirf, irf(svar1) impulse(lexp_inv) response(ling_trib) yline(0,lcolor(black)) ///
	xlabel(0(2)8) byopts(graphregion(color(white))) ylabel(, nogrid) byopts(legend(off)) byopts(yrescale) ///
	subtitle("") byopts(note("")) xtitle("") ciopts1(recast(rline) lpattern(dash) lcolor(red) lwidth(thick)) 		 

**********************************************************************************
*******************CHOQUE DE GASTO CORRIENTE -VAR2-*******************************
**********************************************************************************

qui var ltransfer ling_trib lpib, lags(3) exog(t tq d1 du1 du2 du3 du4 du5 du6 du7 du8 du9 dum1 dum2 dum3) 
*varsoc				//Information criteria sugiere 3 lags 
predict es_var2, residuals equation(#1)
predict et_var2, residuals equation(#2)
predict ey_var2, residuals equation(#3)

g et_t = et_var2 - (a1 * ey_var2)					      //Ajuste cíclico a la recaudación conforme a B y P (2002)

qui ivregress 2sls ey_var2 et_var2 (es_var2 = et_t)		 //Estimación de c1 y c2 mediante variables instrumentales
scalar c1_t = _b[et_var2]
scalar c2_t = _b[es_var2]

matrix A_t = (1,0,0\0,1,-a1\-c1_t,-c2_t,1)				 //Matriz de choques con identificación a la B y P (2002)
matrix B_t = (1,0,0\.,1,0\0,0,1)

qui svar ltransfer ling_trib lpib, lags(1/3) aeq(A_t) beq(B_t) exog(t tq d1 du1 du2 du3 du4 du5 du6 du7 du8 du9) 
	irf create svar2, step(8) set(var2) replace
	irf table fevd,  irf(svar2) impulse(ltransfer) response(ling_trib)
	irf table coirf, irf(svar2) impulse(ltransfer) response(ling_trib) std
	irf table coirf, irf(svar2) impulse(ltransfer) response(ltransfer) std
	irf graph coirf, irf(svar2) impulse(ltransfer) response(ling_trib) yline(0,lcolor(black)) ///
	xlabel(0(2)8) byopts(graphregion(color(white))) ylabel(, nogrid) byopts(legend(off)) ///
	subtitle("") byopts(note("")) xtitle("") ciopts1(recast(rline) lpattern(dash) lcolor(red) lwidth(thick))
	
/* 
ESTIMACIÓN DE MODELO MEDIANTE LOCAL PROJECTIONS
*/

g h = t-1

forvalues i = 0/8 {									//Variable cumulativa para estimar proyecciones locales
	gen cum`i'_ing_trib = f`i'.ing_trib - l.ing_trib
}

foreach var in exp_inv_ ing_trib_ transfer_ {
  gen cum_`var' = 0
	gen mult_`var' = 0
		gen semult_`var' = 0
  forvalues i = 0/8 {
	gen cum`i'_`var' = F`i'.`var' + cum_`var'
		replace cum_`var' = cum`i'_`var'
  }	
}

**********************************************************************************
*********************CHOQUE DE GASTO EN INVERSIÓN*********************************
**********************************************************************************
qui reg lexp_inv L(1/2).lexp_inv L(1/2).ling_trib  t 
predict inv_shock, resid

gen mult_inv   = 0
gen semult_inv = 0

forvalues i = 0/8 {
	qui ivreg2 cum`i'_ing_trib_ (cum`i'_exp_inv_ = inv_shock) L(1/2).inv_shock L(1/2).ing_trib_ L(1/2).exp_inv_ t, robust bw(auto)							 
		quietly replace mult_inv = _b[cum`i'_exp_inv_] if h==`i'
			quietly replace semult_inv = _se[cum`i'_exp_inv_] if h==`i'
}

qui gen 	up95_inv  = .
qui gen 	lo95_inv  = .
qui replace up95_inv  = mult_inv + (1.96*semult_inv)
qui replace lo95_inv  = mult_inv - (1.96*semult_inv)
tw (rarea up95_inv lo95_inv h, bcolor(gs14) clw(medthin medthin))(scatter mult_inv h, c(l ) clp(l ) ms(i ) clc(blue) mc(black) clw(medthick)) if h<=8 ///
, graphregion(color(white)) legend(off) xtitle("") ylabel(, nogrid) yline(0, lcolor(black))

**********************************************************************************
*********************CHOQUE DE GASTO EN TRANSFERENCIAS****************************
**********************************************************************************
qui reg ltransfer L(1/2).ltransfer L(1/2).ling_trib t 
predict tra_shock, resid

gen mult_tra   = 0
gen semult_tra = 0

forvalues i = 0/8 {
	qui ivreg2 cum`i'_ing_trib_ (cum`i'_transfer_ = tra_shock) L(1/2).tra_shock L(1/2).ing_trib_ L(1/2).transfer_ t, robust bw(auto)							 
		quietly replace mult_tra = _b[cum`i'_transfer_] if h==`i'
			quietly replace semult_tra = _se[cum`i'_transfer_] if h==`i'
}

qui gen 	up95_tra  = .
qui gen 	lo95_tra  = .
qui replace up95_tra  = mult_tra + (1.96*semult_tra)
qui replace lo95_tra  = mult_tra - (1.96*semult_tra)
tw (rarea up95_tra lo95_tra h, bcolor(gs14) clw(medthin medthin))(scatter mult_tra h, c(l ) clp(l ) ms(i ) clc(blue) mc(black) clw(medthick)) if h<=8 ///
, graphregion(color(white)) legend(off) xtitle("") ylabel(, nogrid) yline(0, lcolor(black))


***FIN***

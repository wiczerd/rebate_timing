clear all
set more off



set seed 56
use individual

gen holdmonth=month
replace month=12+month if year==2009
replace month=24+month if year==2010

capture drop *rebate* *rebdum*

gen rebate2=0
gen rebdum2=0
gen rebmonth2=0
gen lagrebate2=0
gen lagrebdum2=0


for num 1/3 13/15: gen amtX=0
for num 1/3 13/15: gen rebX=0

foreach num of numlist 5(1)15 {
	local num2=`num'-1
	replace rebate2=amt`num'+amt`num2' if month==`num'
	replace rebdum2=reb`num'+reb`num2' if month==`num'
	local num2=`num'-3
	local num3=`num'-2
	
	replace lagrebate2=amt`num2'+amt`num3' if month==`num'  
	replace lagrebdum2=reb`num2'+reb`num3' if month==`num' 
	
}

foreach var of varlist *rebate* *rebdum* {
	drop if `var'==.
}

sort ssuid epppn year month

egen minmonth=min(month), by(ssuid epppnum)
egen pp=group(ssuid epppnum)
replace famsize=6 if famsize>6


replace tpearn=0 if tpearn<0
xtset pp month


keep if swave<3
bys pp: keep if _N==8
egen time=group(famsize month spouse)

gen positiveearn=(tpearn>0)
gen unpaidabs=(eawop==1)
gen ill=(eabre==3 | eabre==4)
gen vacation=(eabre==7)
gen layoff=(eabre==1)
gen slack=(eabre==2)
gen childcare=(eabre==6)


tab time, gen(hh) nof
for var rebate2 lagrebate2: replace X=X/1000

est clear
local nnn=0

**monthly labor force status
foreach num of numlist 1/5 {
	gen dep=(rmesr==`num')
	xtivreg2 dep hh* (rebate2 lagrebate2 = rebdum2 lagrebdum2) [aw=whf], i(pp) fe cluster(pp) fwl(hh*)
	est store est`nnn', title(rmesr`nnn')
	local nnn=`nnn'+1
	drop dep
}

*postive earnings
xtivreg2 positive hh* (rebate2 lagrebate2 = rebdum2 lagrebdum2) [aw=whf], i(pp) fe cluster(pp) fwl(hh*)
est store est`nnn', title(positive)
local nnn=`nnn'+1




foreach var of varlist rebate2 lagrebate2 rebdum2 lagrebdum2 {
	egen temp=max(`var'), by(pp swave)
	replace `var'=temp
	drop temp
}


capture drop minmonth
egen minmonth=min(month), by(pp swave)


***do analysis by wave
bys pp swave: keep if _n==1
drop hh* time
egen time=group(famsize minmonth spouse swave)
tab time, gen(hh) nof



foreach var of varlist unpaidabs ill vacation layoff slack childcare {
	xtivreg2 `var' hh* (rebate2 lagrebate2 = rebdum2 lagrebdum2) [aw=whf], i(pp) fe cluster(pp) fwl(hh*)
	est store est`nnn', title(`var'`nnn')
	local nnn=`nnn'+1

}


estout est* using Table6.txt, order(rebate2 lagrebate2) keep(rebate2 lagrebate2) replace cells(b(star fmt(%9.3f)) se(par(`"="("'`")""'))) stats(N, fmt(%9.0fc)) starl(* 0.1 ** 0.05 *** 0.01) style(fixed) label




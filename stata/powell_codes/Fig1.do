***main CDF estimates

clear all
set more off


set seed 56

use household, clear

gen holdmonth=month
replace month=12+month if year==2009
replace month=24+month if year==2010


gen rebate=0
gen rebdum=0

gen lagrebate=0
gen lagrebdum=0

gen leadrebate=0
gen leadrebdum=0

for num 1/3 13/20: gen amtX=0
for num 1/3 13/20: gen rebX=0

foreach num of numlist 5(1)15 {
	local num2=`num'-1
	replace rebate=amt`num'+amt`num2' if month==`num'
	replace rebdum=reb`num'+reb`num2' if month==`num'
	local num2=`num'-3
	local num3=`num'-2
	
	replace lagrebate=amt`num2'+amt`num3' if month==`num'  
	replace lagrebdum=reb`num2'+reb`num3' if month==`num' 

	local num2=`num'+2
	local num3=`num'+1

	replace leadrebate=amt`num2'+amt`num3' if month==`num'  
	replace leadrebdum=reb`num2'+reb`num3' if month==`num' 
	
}



sort ssuid epppn year month

egen minmonth=min(month), by(ssuid epppnum)
egen pp=group(ssuid epppnum)
replace famsize=6 if famsize>6



replace tpearn= tfearn
replace tpearn=0 if tpearn<0


gen rebate_month=4 if reb4==1
for num 5/12: replace rebate_month=X if rebX==1


keep if swave<3
bys pp: keep if _N==8
egen time=group(famsize month spouse minmonth)

tab time, gen(hh) nof

for var *rebate: replace X=X/1000
gen beta=.
gen se=.
gen beta2=.
gen se2=.
gen beta3=.
gen se3=.

gen num=_n


foreach num of numlist 1(1)120 {
	gen dp=(tpearn<`num'*100)

	xtivreg2 dp hh* (rebate lagrebate = rebdum lagrebdum) , i(pp) fe cluster(pp) fwl(hh*)
	replace beta=_b[rebate] if num==`num'
	replace se=_se[rebate] if num==`num'
	replace beta2=_b[lagrebate] if num==`num'
	replace se2=_se[lagrebate] if num==`num'

	drop dp
}

keep beta* se* num
drop *3
keep if num<121
gen iter=0
save CDFmain.dta, replace


use CDFmain, clear

local bw=1.96

gen quan=num*100
gen high=beta+se*`bw'
gen low=beta-se*`bw'

for num 2: gen highX=betaX+seX*`bw'
for num 2: gen lowX=betaX-seX*`bw'



for var high*: replace X=.012 if X>.012
for var low*: replace X=-.005 if X<-.005


twoway (connected beta quan, sort) (rcap high low quan, sort), legend(label(1 "Coefficient") label(2 "95% Confidence Interval") rows(1)) ylabel(-.005(.005).012) xlabel(0(2000)12000, angle(45)) yline(0) xtitle("Monthly Household Earnings ($)") ytitle("Effect of $1,000 Rebate on" "Probability of Household Earnings Below x", axis(1))  
gr save CDFmain.gph, replace
gr export CDFmain.eps, replace

twoway (connected beta2 quan, sort) (rcap high2 low2 quan, sort), legend(label(1 "Coefficient") label(2 "95% Confidence Interval") rows(1)) ylabel(-.005(.005).012) xlabel(0(2000)12000, angle(45)) yline(0) xtitle("Monthly Household Earnings ($)") ytitle("Effect of Lagged $1,000 Rebate on" "Probability of Household Earnings Below x", axis(1))  
gr save CDFmain_lag.gph, replace
gr export CDFmain_lag.eps, replace



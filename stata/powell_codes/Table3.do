/***This file gets subsamples
	while using main.dta (from Fig3.do)
	for main estimates
***/

clear all
set more off


set seed 56

use household
gen holdmonth=month
replace month=12+month if year==2009
replace month=24+month if year==2010


gen rebate=0
gen rebdum=0

gen lagrebate=0
gen lagrebdum=0




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
	
}

foreach var of varlist *rebate* *rebdum* {
	drop if `var'==.
}

sort ssuid epppn year month

egen minmonth=min(month), by(ssuid epppnum)
egen pp=group(ssuid epppnum)
replace famsize=6 if famsize>6



replace tpearn= tfearn
replace tpearn=0 if tpearn<0



keep if swave<3
bys pp: keep if _N==8
egen time=group(famsize month spouse minmonth)


rename whf www

keep tpearn rebate lagrebate rebdum lagrebdum pp time spouse month www famsize

for num 1/2: gen betaX=.
for num 1/2: gen seX=.
for num 1/2: gen cnsX=.
gen samplesize=.

gen mmark=0
foreach var of varlist rebdum lagrebdum {
	replace mmark=1 if `var'==1
}


bys pp: gen ff=1 if _n==1


foreach hhh of numlist 1/1000 {
	capture drop rktmp
	gen temp=uniform() if ff==1
	egen rktemp=rank(temp)
	replace rktemp=. if ff~=1
	egen rktmp=max(rktemp), by(pp)
	drop temp-rktemp
	local firr=1

foreach tau of numlist 0.06(.01).99 {
	preserve
	*subsample 3000 households
	keep if rktmp<=3000
	local tau2=`tau'*100

	sum time
	local mxx=r(max)
	qui foreach nmm of numlist 1/`mxx' {
		centile tpearn if time==`nmm' & mmark==0, centile(`tau2')
		if (r(c_1)<100 | r(c_1)>12000) {
			drop if time==`nmm'
		}

	}
	count
	replace samplesize=r(N)

	egen time2=group(time)

	gen quan=round(`tau2',1)
	qregpd tpearn rebate lagrebate  [aw=www], instr(rebdum lagrebdum) quantile(`tau') id(pp) fix(time2)  grid1(-.6(.01).3) optimize(grid) grid2(-.6(.01).3)


	matrix A=e(b)
	matrix B=e(V)
	scalar a1=A[1,1]
	scalar a2=A[1,2]
	replace beta1=a1 
	replace beta2=a2
	scalar b1=B[1,1] 
	scalar b2=B[2,2]
	replace se1=b1^.5
	replace se2=b2^.5

	gen temp=tpearn-a1*rebate-a2*lagrebate
	_pctile temp if spouse==0 [aw=www], percentiles(`tau2')
	replace cns1=r(r1)
	_pctile temp if spouse==1 [aw=www], percentiles(`tau2')
	replace cns2=r(r1)
	keep if _n==1
	keep beta1 beta2 se1 se2 quan cns* samplesize

	if (`firr'==1) {
		save sub`hhh'.dta, replace
		local firr=2
	}
	else {
		append using sub`hhh'
		save sub`hhh', replace
	}
	restore
}
}



****generate Table 3 estimates


clear all
use main
keep if quan>5
egen temp=sum(beta1)
replace temp=temp+beta1 if quan==99
replace temp=temp/100

rename temp meanest
sum meanest

egen temp=sum(beta2)
replace temp=temp+beta2 if quan==99
replace temp=temp/100

rename temp meanest_lag 


gen totaleffect=(meanest+meanest_l)*2*96000000000
***MEAN AND TOTAL ESTIMATES
sum meanest* totaleffect if quan==99

use sub1, clear
gen mark=1
local num2=1

qui foreach num of numlist 2/1000 {
	append using sub`num'
	replace mark=`num' if mark==.
	local num2=`num2'+1
}

keep if quan>5
egen temp=sum(beta1), by(mark)
replace temp=temp+beta1 if quan==99
replace temp=temp/100

rename temp meanest

egen temp=sum(beta2), by(mark)
replace temp=temp+beta2 if quan==99
replace temp=temp/100

rename temp meanest_lag 


gen totaleffect=(meanest+meanest_l)*2*96000000000
***STANDARD ERRORS
foreach var of varlist meane* totalef {
	sum `var' if quan==99
	display r(sd)*(3000/22998)^.5
}


***use samplesize variable - assume no effect for censored cells
use main, clear
keep if quan>5
egen temp=sum(beta1*samplesize)
replace temp=temp+beta1*samplesize if quan==99
egen temp2=sum(samplesize)
replace temp2=temp2+samplesize+5*22998*8 if quan==99

replace temp2=100*22998*8
replace temp=temp/temp2

rename temp meanest
sum meanest


egen temp=sum(beta2*samplesize)
replace temp=temp+beta2*samplesize if quan==99
replace temp=temp/temp2

rename temp meanest_lag 


gen totaleffect=(meanest+meanest_l)*2*96000000000
**PANEL B ESTIMATES
sum meanest* totaleffect if quan==99

use sub1, clear
gen mark=1
local num2=1

qui foreach num of numlist 2/196 {
	append using sub`num'
	replace mark=`num' if mark==.
	local num2=`num2'+1
}



keep if quan>5
egen temp=sum(beta1*samplesize), by(mark)
replace temp=temp+beta1*samplesize if quan==99
egen temp2=sum(samplesize), by(mark)
replace temp2=temp2+samplesize+5*3000*8 if quan==99
replace temp2=100*3000*8
replace temp=temp/temp2

rename temp meanest
sum meanest


egen temp=sum(beta2*samplesize), by(mark)
replace temp=temp+beta2*samplesize if quan==99
replace temp=temp/temp2

rename temp meanest_lag 


gen totaleffect=(meanest+meanest_l)*2*96000000000
foreach var of varlist meane* totalef {
	sum `var' if quan==99
	***PANEL B STANDARD ERRORS
	display r(sd)*(3000/22998)^.5
}




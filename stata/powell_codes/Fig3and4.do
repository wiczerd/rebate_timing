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

local firr=1
foreach tau of numlist 0.06(.01).99 {
	preserve
	local tau2=`tau'*100

	centile tpearn if spouse==1, centile(`tau2')
	local spouse=r(c_1)
	centile tpearn if spouse==0, centile(`tau2')
	local single=r(c_1)
	centile tpearn, centile(`tau2')
	local all=r(c_1)

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
		save main.dta, replace
	}
	else {
		append using main
		save main, replace
	}
	restore

	preserve
	gen temp=tpearn-a1*rebate-a2*lagrebate
	centile temp if spouse==1, centile(`tau2')
	local spouse_un=r(c_1)
	centile temp if spouse==0, centile(`tau2')
	local single_un=r(c_1)
	centile temp , centile(`tau2')
	local all_un=r(c_1)
	gen spouse_treat=`spouse'
	gen spouse_untreat=`spouse_un'
	gen single_treat=`single'
	gen single_untreat=`single_un'	
	gen all_treat=`all'
	gen all_untreat=`all_un'	
	gen quan=round(`tau2',1)
	keep spouse_* single_* quan all_*
	keep if _n==1
	if (`firr'==1) {
		save main_counter.dta, replace
		local firr=2
	}
	else {
		append using main_counter
		save main_counter, replace
	}
	restore	
}



use main_counter, clear
keep if quan<96 & quan>6
gen diff=all_treat-all_untreat

twoway (connected diff quan, msize(vtiny) sort) , legend(label(1 "Change in Distribution") rows(1)) xlabel( 10(5)95, angle(45)) ylabel(-100(25)0) xtitle("Quantile") ytitle("Change in Household Earnings ($)" "Due to Rebates", axis(1))  yline(0)
gr save counterfactual_diff.gph, replace
gr export counterfactual_diff.eps, replace

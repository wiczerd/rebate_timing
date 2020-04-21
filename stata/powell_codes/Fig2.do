clear all
set more off


set seed 56
use household

gen holdmonth=month
replace month=12+month if year==2009
replace month=24+month if year==2010


gen rebate=0
gen rebdum=0

for num 1/2: gen lagrebateX=0
for num 1/2: gen lagrebdumX=0

for num 1/3 13/22: gen amtX=0
for num 1/3 13/22: gen rebX=0

foreach num of numlist 5(1)15 {


	local num2=`num'-1
	replace rebate=amt`num'+amt`num2' if month==`num'
	replace rebdum=reb`num'+reb`num2' if month==`num'
	local num2=`num'-3
	local num3=`num'-2
	
	replace lagrebate1=amt`num2'+amt`num3' if month==`num'  
	replace lagrebdum1=reb`num2'+reb`num3' if month==`num' 

	if (`num'>5) {
		local num2=`num'-5
		local num3=`num'-4
		replace lagrebate2=amt`num2'+amt`num3' if month==`num'  
		replace lagrebdum2=reb`num2'+reb`num3' if month==`num' 
	}
	else {
		local num3=`num'-4
		replace lagrebate2=amt`num3' if month==`num'  
		replace lagrebdum2=reb`num3' if month==`num' 
	}


	local num2=`num'-4 
	
	if `num'> 4 {
		local nnn=1
		while `nnn'<=`num2' {
			replace lagrebate2=lagrebate2+amt`nnn'  if month==`num'  
			replace lagrebdum2=lagrebdum2+reb`nnn' if month==`num' 
			local nnn=`nnn'+1
		}
	}
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


tab time, gen(hh) nof


gen beta=.
gen se=.
gen beta2=.
gen se2=.
gen beta3=.
gen se3=.
gen beta4=.
gen se4=.
gen beta5=.
gen se5=.
for num 6/12: gen betaX=.
for num 6/12: gen seX=.

gen num=_n

for var rebate lagrebate*: replace X=X/1000

qui foreach num of numlist 1(1)120 {
	gen dp=(tpearn<`num'*100)

	xtivreg2 dp hh* (rebate lagrebate* = rebdum lagrebdum*) [aw=whf], i(pp) fe cluster(pp) fwl(hh*) 

	replace beta3=_b[rebate] if num==`num'
	replace se3=_se[rebate] if num==`num'
	replace beta4=_b[lagrebate1] if num==`num'
	replace se4=_se[lagrebate1] if num==`num'
	replace beta5=_b[lagrebate2] if num==`num'
	replace se5=_se[lagrebate2] if num==`num'

	drop dp
}
keep beta* se* num
keep if num<121


local bw=1.96



gen quan=num*100


for num 3/6: gen highX=betaX+seX*`bw'
for num 3/6: gen lowX=betaX-seX*`bw'

for var high*: replace X=.016 if X>.016
for var low*: replace X=-.01 if X<-.01

twoway (connected beta3 quan, sort) (rcap high3 low3 quan, sort), legend(label(1 "Coefficient") label(2 "95% Confidence Interval") rows(1)) xlabel(0(2000)12000, angle(45)) yline(0) ylabel(-.01(0.002)0.016, axis(1) angle(0)) xtitle("Monthly Household Earnings ($)") ytitle("Effect of $1,000 Rebate on" "Probability of Household Earnings Below x", axis(1))  
gr save CDFlags2_rebate.gph, replace
gr export CDFlags2_rebate.eps, replace

twoway (connected beta4 quan, sort) (rcap high4 low4 quan, sort), legend(label(1 "Coefficient") label(2 "95% Confidence Interval") rows(1)) xlabel(0(2000)12000, angle(45)) yline(0) ylabel(-.01(0.002)0.016, axis(1) angle(0)) xtitle("Monthly Household Earnings ($)") ytitle("Effect of $1,000 Lagged Rebate on" "Probability of Household Earnings Below x", axis(1))  
gr save CDFlags2_lagged1rebate.gph, replace
gr export CDFlags2_lagged1rebate.eps, replace

twoway (connected beta5 quan, sort) (rcap high5 low5 quan, sort), legend(label(1 "Coefficient") label(2 "95% Confidence Interval") rows(1)) xlabel(0(2000)12000, angle(45)) yline(0) ylabel(-.01(0.002)0.016, axis(1) angle(0)) xtitle("Monthly Household Earnings ($)") ytitle("Effect of $1,000 Lagged (4+ months) Rebate on" "Probability of Household Earnings Below x", axis(1))  
gr save CDFlags2_lagged2rebate.gph, replace
gr export CDFlags2_lagged2rebate.eps, replace




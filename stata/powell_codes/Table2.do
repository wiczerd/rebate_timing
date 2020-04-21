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

gen everreb=reb4
for num 5/12: replace everreb=everreb+rebX

keep if swave<3
bys pp: keep if _N==8
egen time=group(famsize month spouse minmonth)

rename whf www

keep tpearn rebate lagrebate rebdum lagrebdum pp time spouse month www famsize everr reb* amt*

sort pp month
gen change_spouse=0 if pp==pp[_n-1]
replace change_sp=1 if pp==pp[_n-1] & spouse~=spouse[_n-1]


sort pp month
by pp: gen first=1 if _n==1
egen total=sum(first), by(spouse)
bys spouse: sum total
drop total
count


*present summary statistics by marital status in first observation
gen temp=spouse if first==1
egen temp2=max(temp), by(pp)
egen temp3=min(temp), by(pp)
count if first==1
count if first==1 & temp2==temp3
replace spouse=temp2
drop temp*

gen noearn=(tpearn==0)
bys spouse: sum everr [aw=www] if first==1
bys spouse: sum rebate [aw=www] if rebate>0 & first==1
bys spouse: sum famsize if first==1 [aw=www]
bys spouse: tab noearn [aw=www]
bys spouse: sum tpearn [aw=www], d


keep if everr==1
tab spouse if first==1
bys spouse: sum everr [aw=www]
bys spouse: sum rebate [aw=www] if rebate>0
bys spouse: sum famsize [aw=www] if first==1
bys spouse: tab noearn [aw=www]
bys spouse: tab noearn
bys spouse: sum tpearn [aw=www], d

preserve
keep if reb4+reb5==1
tab spouse if first==1
egen total=sum(first), by(spouse)
bys spouse: sum total
drop total
bys spouse: sum everr [aw=www]
bys spouse: sum rebate [aw=www] if rebate>0
bys spouse: sum famsize [aw=www] if first==1
bys spouse: tab noearn [aw=www]
bys spouse: sum tpearn [aw=www], d
restore



keep if reb8+reb9+reb10+reb11+reb12==1
tab spouse if first==1
egen total=sum(first), by(spouse)
bys spouse: sum total
drop total
bys spouse: sum everr [aw=www]
bys spouse: sum rebate [aw=www] if rebate>0
bys spouse: sum famsize [aw=www] if first==1
bys spouse: tab noearn [aw=www]
bys spouse: sum tpearn [aw=www], d






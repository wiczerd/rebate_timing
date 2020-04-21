/* cross sectional analysis */
/* Simple transition matrices */
use ${sipp_rebate}data/rebatedata_cleaned.dta, clear
keep if family_num1==1 				// keep household who have a single subfamily (85% of observations anyway) in the first wave/month


keep if inlist(wave,1,2) 			// he only uses waves 1,2
count


* keep only households whos head is between 25 and 60 in the first month
keep if inrange(hhhead_age,25,60)
count

egen flag = tag(hhid)			// so we only use one observation per month per hh
sum  hh_rebate_ever if hhhead_errp==1 & inlist(ms,1,2) 
sum  hh_rebate_ever if hhhead_errp==1 & inlist(ms,3,6) & flag==1

// calculate sum of rebate over the month


// Mean regressions
keep if inlist(hhhead_errp,1)


local lhs_vars foodstamp welfare foodass pubhousing energy ui
local nlhs : word count `lhs_vars'
*matrix transitionrts = J()
foreach var in `lhs_vars' {
	tab hh_`var'_ind [aw=whfnwgt] if L.hh_`var'_ind==1 & hh_rebate_ever==1 & L.hh_rebate_ind == 1
	tab hh_`var'_ind [aw=whfnwgt] if L.hh_`var'_ind==1 & hh_rebate_ever==1 & L.hh_rebate_ind == 0
	tab hh_`var'_ind [aw=whfnwgt] if L.hh_`var'_ind==0 & hh_rebate_ever==1 & L.hh_rebate_ind == 1
	tab hh_`var'_ind [aw=whfnwgt] if L.hh_`var'_ind==0 & hh_rebate_ever==1 & L.hh_rebate_ind == 0
	
}






*** Davide: Maybe we should use hh_rebate_ind instead hh_rebate_ever?
* and then loop over time.. but I might be wrong, I'm trying to understand the structure of the data.
xtset hhid date

/*
local lhs_vars foodstamp welfare foodass pubhousing energy ui
local nlhs : word count `lhs_vars'
foreach var in `lhs_vars' {
xttrans hh_`var'_ind if hh_rebate_ever == 1, freq
xttrans hh_`var'_ind if hh_rebate_ind == 1 & hh_rebate_ever == 1, freq
xttrans hh_`var'_ind if hh_rebate_ind == 0 & hh_rebate_ever == 1, freq
*xttrans hh_`var'_ind [aw=whfnwgt] if hh_rebate_ever == 1
*xttrans hh_`var'_ind [aw=whfnwgt] if hh_rebate_ind == 1 & hh_rebate_ever == 1
*xttrans hh_`var'_ind [aw=whfnwgt] if hh_rebate_ind == 0 & hh_rebate_ever == 1
}
*/

tab hh_foodstamp_ind [aw=whfnwgt] if L.hh_foodstamp_ind == 0 & hh_rebate_ind == 1 & datedum2 == 1
tab hh_foodstamp_ind [aw=whfnwgt] if L.hh_foodstamp_ind == 0 & hh_rebate_ind == 0 & datedum2 == 1

reg hh_foodstamp_ind L.hh_foodstamp_ind [aw=whfnwgt] if hh_rebate_ind == 1
reg hh_foodstamp_ind L.hh_foodstamp_ind [aw=whfnwgt] if hh_rebate_ind == 0

bys hhid: gen lagfs = L.hh_foodstamp_ind
gen lagfsrebY = lagfs*hh_rebate_ind

reg hh_foodstamp_ind lagfs lagfsrebY [aw=whfnwgt]
xtabond hh_foodstamp_ind lagfsrebY  // weights not allowed
reg hh_foodstamp_ind hh_rebate_ind [aw=whfnwgt] if lagfs == 0



/* david's dropbox 
gl sipp_rebate ~/Dropbox/rebate_timing/
*/

// laura's dropbox 
gl sipp_rebate /Users/gregorjarosch/Dropbox/rebate_timing/
use ${sipp_rebate}data/sippsets08_MPC.dta, clear

//  copy the rebate month to every observation for that person
replace erbamth=. if erbamth==-1
bys uid: egen rebate_month = max(erbamth)

// indicator for month of rebate receipt for the individual
gen p_rebate_ind =(rebate_month == month)
label var p_rebate_ind "indicator for individual rebate reciept this month"

// replace missing or no rebate with 0
replace erbatamt = 0 if (rebate_month!=month)
replace erbatamt = 0 if erbatamt==.

*gen hhhead_errp1 = (inlist(errp,1,2) & epppnum=="0101" & wave==1 & srefmon==1) 			// does this id head-of household? David's guess: but counts don't line up later
gen hhhead_errp1 = (inlist(errp,1,2) & wave==1 & srefmon==1)			// Laura's new definition: just the reference person
bys uid: egen hhhead_errp = max(hhhead_errp1)

// Create sub/primary family ID
egen sfid = group(ssuid shhadid rfid2 rsid)

gen temp = age if hhhead_errp==1 & wave==1 & srefmon==1	 			// will be missing for everyone except the hhhead in the first appearance
bys hhid: egen hhhead_age = min(temp)			// age of your household head first time they in sample
label var hhhead_age "age of the household head at survey entry"

drop temp

recode erebate (2=0), gen(p_rebate_ever)			// now its a dummy=1 if you got the rebate ever

recode ms (1/2=1) (nonmiss=0), gen(married)

// did anyone in the household receive a rebate this month:
bys hhid date: egen hh_rebate_ind  = max(p_rebate_ind)
label var hh_rebate_ind "indicator for household rebate receipt this month"

gen hh_foodstamp_ind =(thfdstp>0)
label var hh_foodstamp_ind "indicator for household food stamp receipt this month"
sort uid date
tsset uid date
gen hh_foodstamp_new_ind=(L.hh_foodstamp_ind==0 & hh_foodstamp_ind==1)
label var hh_foodstamp_new_ind "indicator for household started food stamps this month"
gen hh_foodstamp_off_ind=(L.hh_foodstamp_ind==1 & hh_foodstamp_ind==0)
label var hh_foodstamp_off_ind "indicator for household went off food stamps this month"


gen hh_foodass_ind =(efoodtp1==1 | efoodtp2==2 | efoodtp3==1 | efoodtp4==1)
label var hh_foodass_ind "indicator for household food assistance receipt this month"

gen hh_foodass_new_ind=(L.hh_foodass_ind==0 & hh_foodass_ind==1)
label var hh_foodass_new_ind "indicator for household started food assistance this month"
gen hh_foodass_off_ind=(L.hh_foodass_ind==1 & hh_foodass_ind==0)
label var hh_foodass_off_ind "indicator for household went off food assistance this month"



recode eafood1 (1/2=0) (3/4=1), gen(hh_food_insuf)
label var hh_food_insuf "indicator for household food insufficiency this month (AW)"

bys hhid: egen hh_rebate_ever  = max(hh_rebate_ind)
label var hh_rebate_ever "indicator for household rebate receipt ever"

bys hhid date: egen hh_erbatamt = total(erbatamt)		// household-level rebate values
replace hh_erbatamt = 0 if (hh_rebate_ind!=1)
label var hh_erbatamt "household rebate value"

bys hhid date: egen hh_size = count(errp)		// household size in each date
recode hh_size (1=1) (2=2) (3=3) (4=4) (5=5) (nonmiss=6), gen(hh_size2)


gen cohabiting = (errp==10)
bys hhid date: egen cohab2=max(cohabiting)

save ${sipp_rebate}data/rebatedata_cleaned.dta, replace

*************************************************************
* let's try to replicate Powell's stuff 
*************************************************************
local cohab_drop 0

use ${sipp_rebate}data/rebatedata_cleaned.dta, clear
keep if inlist(wave,1,2) 			// he only uses waves 1,2
count

* drop if missing household labor earnings in any of the 8 months you are present *
bys hhid uid: egen p_count_nonmiss = count(thearn) // non-missing household earnings count
gen head_nonmiss_thearn = p_count_nonmiss if hhhead_errp==1
bys hhid: egen hh_count_nonmiss = max(head_nonmiss_thearn) // non-missing household earnings count
keep if hh_count_nonmiss==8
count
drop hh_count_nonmiss

* exclude household with cohabiting non-married couples 
if `cohab_drop'==1 {
	gen temp = (inrange(ms,3,6) & cohab2==1 & wave==1 & srefmon==1 & hhhead_errp==1)
	bys hhid: egen cohab3 = max(temp)
	drop if cohab3==1
	count
	}

* keep only households whos head is between 25 and 60
keep if inrange(hhhead_age,25,60)
count




egen flag = tag(hhid)			// so we only use one observation per month per hh
sum  hh_rebate_ever if hhhead_errp==1 & inlist(ms,1,2) 
sum  hh_rebate_ever if hhhead_errp==1 & inlist(ms,3,6) & flag==1

// calculate sum of rebate over the month
tsset uid date
sort uid date
gen hhRit = hh_erbatamt+L.hh_erbatamt
gen hhLRit = L3.hh_erbatamt+L4.hh_erbatamt
gen hhRit_ind= (hhRit>0)
gen hhLRit_ind = (hhLRit>0)

gen Rit = erbatamt+L.erbatamt
gen LRit = L3.erbatamt+L4.erbatamt
gen Rit_ind = (Rit>0)
gen LRit_in = (LRit>0)

// Table 2 stats //
sum hh_size2 if inrange(ms,3,6) & wave==1 & srefmon==1 & hhhead_errp==1
sum hh_size2 if inrange(ms,1,2) & wave==1 & srefmon==1 & hhhead_errp==1
tab hh_rebate_ever if inrange(ms,3,6) & wave==1 & srefmon==1 & hhhead_errp==1
tab hh_rebate_ever if inlist(ms,1,2) & wave==1 & srefmon==1 & hhhead_errp==1
sum erbatamt hh_erbatamt if erbatamt>0 & inrange(ms,3,6) & hhhead_errp==1
sum erbatamt hh_erbatamt if erbatamt>0 & inrange(ms,1,2) & hhhead_errp==1

gen have_earning = (thearn>0)
tab have_earning if inrange(ms,3,6)  & hhhead_errp==1
tab have_earning if inrange(ms,1,2)  & hhhead_errp==1

sum thearn if inrange(ms,3,6)  & hhhead_errp==1, det
sum thearn if inrange(ms,1,2)  & hhhead_errp==1, det
count if inrange(ms,3,6) & wave==1 & srefmon==1 & hhhead_errp==1 
count if inlist(ms,1,2) & wave==1 & srefmon==1 & hhhead_errp==1 

// Mean regressions
keep if inlist(hhhead_errp,1)

// hh rebate variables
qui areg thearn hhRit hhLRit i.date#i.srefmon#i.married#i.hh_size2 [aw=whfnwgt], absorb(hhid)
est table, star(.05 .01 .001) keep(hhRit hhLRit) 

qui areg thearn hhRit hhLRit i.date#i.srefmon#i.hh_size2  if married==0 [aw=whfnwgt], absorb(hhid)
est table, star(.05 .01 .001) keep(hhRit hhLRit) 

// hh head rebate variables
qui areg thearn Rit LRit i.date#i.srefmon#i.married#i.hh_size2  [aw=whfnwgt], absorb(hhid)
est table, star(.05 .01 .001) keep(Rit LRit) 

qui areg thearn Rit LRit i.date#i.srefmon#i.hh_size2  if married==0 [aw=whfnwgt], absorb(hhid)
est table, star(.05 .01 .001) keep(Rit LRit) 


qui areg thearn Rit LRit i.date#i.srefmon#i.married#i.hh_size2 [aw=whfnwgt], absorb(hhid)
est table, star(.05 .01 .001) keep(Rit LRit) 


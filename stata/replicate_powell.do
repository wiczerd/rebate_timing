*************************************************************
* let's try to replicate Powell's stuff 
*************************************************************
local cohab_drop 1

use ${sipp_rebate}data/rebatedata_cleaned.dta, clear

*gen temp = family_num2 if wave==1 & srefmon==1
*bys hhid: egen family_num2a = min(temp)
*drop temp
keep if family_num2==1 // keep household who have a single subfamily (85% of observations anyway) in the first wave/month


keep if inlist(wave,1,2) 			// he only uses waves 1,2
count

* drop if missing household labor earnings in any of the 8 months you are present *
bys uid: egen p_count_nonmiss = count(tfearn) // non-missing household earnings count
gen head_nonmiss_tfearn = p_count_nonmiss if hhhead_errp==1
bys hhid: egen hh_count_nonmiss = max(head_nonmiss_tfearn) // non-missing household earnings count
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

* keep only households whos head is between 25 and 60 in the first month
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

gen have_earning = (tfearn>0)
tab have_earning if inrange(ms,3,6)  & hhhead_errp==1
tab have_earning if inrange(ms,1,2)  & hhhead_errp==1

sum tfearn if inrange(ms,3,6)  & hhhead_errp==1, det
sum tfearn if inrange(ms,1,2)  & hhhead_errp==1, det
count if inrange(ms,3,6) & wave==1 & srefmon==1 & hhhead_errp==1 
count if inlist(ms,1,2) & wave==1 & srefmon==1 & hhhead_errp==1 

// Mean regressions
keep if inlist(hhhead_errp,1)

// hh rebate variables
qui areg tfearn Rit LRit i.date#i.srefmon#i.married#i.hh_size2 [aw=whfnwgt], absorb(hhid) vce(cl hhid)
est table, keep(Rit LRit) se stats(N)


qui areg tfearn Rit LRit i.date#i.srefmon#i.hh_size2 if married==0  [aw=whfnwgt], absorb(hhid)  vce(cl hhid)
est table, keep(Rit LRit) se stats(N)


// laura's dropbox 
gl sipp_rebate /Users/gregorjarosch/Dropbox/rebate_timing/


use ${sipp_rebate}data/sippsets08_MPC.dta, clear
drop if date==. // need to find out from david whether the months are ok


bys uid: egen rebate_month =  erebatemnth
gen got_rebate = 1 if rebate_month==month
replace got_rebate = 0 if got_rebate!=1

*isid uid date // these uniquely identify observations

*************************************************************
* let's try to replicate Powell's stuff 
*************************************************************
keep if inlist(wave,1,2) 			// he only uses waves 1,2
egen hhid = group(ssuid shhadid) 	// identify unique households 

destring epppnum, gen(pnum)			// this is the person number
bys hhid: egen hh_head = min(pnum)	// this identifies the pnum 0101 in each household
keep if hh_head==pnum				// I am assuming this is the head?

* drop if missing household labor earnings in any of the 8 months you are present *
bys hhid: egen count_nonmiss = count(thearn)
keep if count_nonmiss==8
drop count_nonmiss

* exclude household with cohabiting non-married couples
drop if (inrange(ms,3,6) & cohab==1) 

* keep only households whos head is between 25 and 60
keep if inrange(age,25,60)

drop hhid
egen hhid = group(ssuid shhadid) // identify unique households 
sum hhid // this is so I can see how many unique hhhids there are

replace erebate = 0 if erebate==.

bys hhid: egen rebate_ever = max(erebate)
tab rebate_ever if flag==1 & inlist(ms,3,6)


gen erbatamt2 = erbatamt
replace erbatamt2 = 0 if erbatamt2==.

tsset hhid date
sort hhid date
gen Rit = erbatamt2+L.erbatamt2
gen LRit = L3.erbatamt2+L4.erbatamt2


gen rebind = 1 if Rit>0
replace rebind=0 if Rit==0

gen Lrebind = 1 if LRit>0
replace Lrebind=0 if LRit==0


egen flag = tag(hhid)
count if inlist(ms,1,2) & wave==1 & srefmon==1 & flag==1 // number of married households (see Table 2 of powell)
count if inrange(ms,3,6) & wave==1 & srefmon==1 & flag==1 // number of married households (see Table 2 of powell)



qui areg thearn Rit LRit i.month#i.ms#i.srefmon  [aw=whfnwgt], absorb(hhid)
est table, star(.05 .01 .001) keep(Rit LRit) 


use ${sipp_rebate}data/rebatedata_cleaned.dta, clear
keep if family_num1==1 				// keep household who have a single subfamily (85% of observations anyway) in the first wave/month

keep if inlist(wave,6,9)

* keep only households whos head is between 25 and 60 in the first month
keep if inrange(hhhead_age,25,60)
count

egen flag = tag(hhid)			// so we only use one observation per month per hh
sum  hh_rebate_ever if hhhead_errp==1 & inlist(ms,1,2) 
sum  hh_rebate_ever if hhhead_errp==1 & inlist(ms,3,6) & flag==1


// Mean regressions
keep if inlist(hhhead_errp,1)
keep if srefmon==4 // LHS variables are from TM module, only in refmon 4

histogram time_since_lastjobR if erbatamt>0, bin(5)



local longrunvars eabmeet eabrent eabevct eabgas eabdent eafbaln eafchld eafskip eafless
foreach var in `longrunvars' {
	eststo base1: reg  `var' i.date i.married i.hh_size2 if hh_rebate_ever==1 [aw=whfnwgt] 
	eststo base2: reg  `var' i.date i.married i.hh_size2 i.time_since_lastjobR if hh_rebate_ever==1 [aw=whfnwgt] 
	
	esttab base1 base2 using ${sipp_rebate}output/`var'_LRregressions.tex, replace  ///
	nomtitles compress ///
	coeflabels(_cons "intercept" ///
	i0.employed_ind "unemployed upon receipt"  ///
	i.time_since_lastjobR "unemployed for X months")  ///
	drop(*.married *.date *.hh_size2) ///
	substitute(# x _ \_ \sym{***} \onepc \sym{**} \fivepc \sym{*} \tenpc htbp h!) se r2
	}

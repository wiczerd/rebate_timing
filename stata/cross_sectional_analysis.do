/* cross sectional analysis */
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



foreach var in `lhs_vars' {

**** Exit rates:
	# delimit  ;
	eststo base1a: reg hh_`var'_ind [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==1  ;
	# delimit cr ;
	
	# delimit  ;
	eststo base1: areg hh_`var'_ind 
	i.srefmon [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==1 , absorb(date) cl(hhid) ;
	# delimit cr ;

	# delimit  ;
		eststo base2: areg hh_`var'_ind  i1.hhLRit_ind
	i.srefmon [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==1 , absorb(date) cl(hhid) ;
	# delimit cr ;
	

**** Entry rates:	
	# delimit  ;
	eststo base1a: reg hh_`var'_ind [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==0  ;
	# delimit cr ;
	
	# delimit  ;
	eststo base1: areg hh_`var'_ind 
	i.srefmon [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==0 , absorb(date) cl(hhid) ;
	# delimit cr ;

	# delimit  ;
		eststo base2: areg hh_`var'_ind  i1.hhLRit_ind
	i.srefmon [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==0 , absorb(date) cl(hhid) ;
	# delimit cr ;
	
	
*	eststo base1: areg hh_`var'_ind i.date 
*	i.srefmon i.married i.hh_size2 [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==1 , a(hhid) cl(hhid) ;

*	eststo base2: areg hh_`var'_ind i1.hhRit_ind i.date 
*	i.srefmon i.married i.hh_size2 [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==1, a(hhid) cl(hhid) ;
	*# delimit cr ;	
	
/*
	# delimit  ;
	eststo base3: areg hh_`var'_ind i1.hhRit_ind  i.date 
	i.srefmon i.married i.hh_size2 [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==1, a(hhid) cl(hhid) ;
	# delimit cr ;	
	
	# delimit  ;
	eststo base1a_2: areg hh_`var'_ind [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==0  ;
	# delimit cr ;

	# delimit  ;
	eststo base1_2: areg hh_`var'_ind i.date 
	i.srefmon i.married i.hh_size2 [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==0 , a(hhid) cl(hhid) ;
	# delimit cr ;

	# delimit  ;
	eststo base2_2: areg hh_`var'_ind i1.hhRit_ind i.date 
	i.srefmon i.married i.hh_size2 [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==0, a(hhid) cl(hhid) ;
	# delimit cr ;	
	
	# delimit  ;
	eststo base3_2: areg hh_`var'_ind i1.hhRit_ind  i.date 
	i.srefmon i.married i.hh_size2 [aw=whfnwgt] if hh_rebate_ever==1 & L.hh_`var'_ind==0, a(hhid) cl(hhid) ;
	# delimit cr ;	
	
*/		
		
/*		
	*esttab base1a base1 base2 base3 using ${sipp_rebate}output/`var'_regressions_onlast.tex, replace  ///
	esttab base1 base2 using ${sipp_rebate}output/`var'_regressions_onlast.tex, replace  ///
	nomtitles compress ///
	coeflabels(_cons "intercept" ///
	i.srefmon "reference month"  ///
	i0.L.hh_`var'_ind "recepient last month"  ///
	i0.L.hh_`var'_ind#i0.employed_ind "recepient last month and unemployed") ///
	drop(*.date) ///
	substitute(# x _ \_ \sym{***} \onepc \sym{**} \fivepc \sym{*} \tenpc htbp h!) se r2
	
	esttab base1a_1 base1_1 base2_1 base3_1 using ${sipp_rebate}output/`var'_regressions_offlast.tex, replace  ///
	nomtitles compress ///
	coeflabels(_cons "intercept" ///
	i0.L.hh_`var'_ind "recepient last month"  ///
	i0.L.hh_`var'_ind#i0.employed_ind "recepient last month and unemployed") ///
	drop(*.married *.srefmon *.date *.hh_size2) ///
	substitute(# x _ \_ \sym{***} \onepc \sym{**} \fivepc \sym{*} \tenpc htbp h!) se r2
	*/

}


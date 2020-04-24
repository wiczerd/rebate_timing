/* david's dropbox 
gl sipp_rebate ~/Dropbox/rebate_timing/
*/

/* davide's dropbox 
gl sipp_rebate /Users/rcedxm19/Dropbox/rebate_timing/
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


gen temp = age if hhhead_errp==1 & wave==1 & srefmon==1	 			// will be missing for everyone except the hhhead in the first appearance
bys hhid: egen hhhead_age = min(temp)			// age of your household head first time they in sample
label var hhhead_age "age of the household head at survey entry"

drop temp

recode erebate (2=0), gen(p_rebate_ever)			// now its a dummy=1 if you got the rebate ever

recode ms (1/2=1) (nonmiss=0), gen(married)

// did anyone in the household receive a rebate this month:
bys hhid date: egen hh_rebate_ind  = max(p_rebate_ind)
label var hh_rebate_ind "indicator for household rebate receipt this month"


// indicators for welfare programs
gen hh_foodstamp_ind =(thfdstp>0)
label var hh_foodstamp_ind "indicator for household food stamp receipt this month"

gen hh_housing_ind = (hsg_r==1)
label var hh_housing_ind "indicator for household housing assistance receipt this month"

gen hh_foodass_ind =(efoodtp1==1 | efoodtp2==2 | efoodtp3==1 | efoodtp4==1)
label var hh_foodass_ind "indicator for household food assistance receipt this month"

gen hh_welfare_ind =(welfare_r==1)
label var hh_welfare_ind "indicator for household welfare receipt this month"

gen hh_ui_ind =(ui_r==1)
label var hh_ui_ind "indicator for received ui this month"

gen hh_energy_ind =(enrgy_r==1)
label var hh_energy_ind "indicator for household Energy Assistance receipt this month"

gen hh_pubhousing_ind =(pubhou_r==1)
label var hh_pubhousing_ind "indicator for public housing receipt this month"


gen hh_rentsubsidy_ind =(lowrent_r==1)
label var hh_rentsubsidy_ind "indicator for rent subsidy this month"


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

// number of unique rfid2 values within hhhid and date
egen family_num1 = nvals(rfid), by(hhid date) 
egen family_num2 = nvals(rfid2), by(hhid date) 


tsset uid date
sort uid date

gen hhRit_ind= (hh_erbatamt>0)
gen hhLRit_ind = (L.hh_erbatamt>0)

gen Rit_ind = (erbatamt>0)

recode esr (1/5=1) (nonmiss=0), gen(employed_ind)


// generate time since layoff
egen time_since_lastjob = rowtotal(ntime utime), missing

gen temp=time_since_lastjob if (erbatamt>0 & employed_ind==0)
replace temp=-1 if (erbatamt>0 & employed_ind==1)

bys uid: egen time_since_lastjobR = mean(temp)
label var time_since_lastjobR "months since last job upon rebate receipt"

drop temp

gen temp=employed_ind if (erbatamt>0)
bys uid: egen employed_indR = mean(temp)
label var employed_indR "employment status upon rebate receipt"


// attach displacement indicators to separations
foreach dvar of varlist displaced displaced_layoff displaced_slackbiz displaced_empclosed {
	gen `dvar'_stint = .
	replace `dvar'_stint = `dvar' if (EU==1 | EN==1) & srefmon>1
	replace `dvar'_stint = max(l.`dvar',`dvar') if (EU==1 | EN==1) & srefmon==1

	bys uid uspellid: egen displaced_utmp = max(`dvar'_stint) if uspellid<.
	bys uid nspellid: egen displaced_ntmp = max(`dvar'_stint) if nspellid<.

	replace `dvar'_stint = displaced_utmp  if uspellid<.
	replace `dvar'_stint = displaced_ntmp  if nspellid<.

	drop displaced_utmp displaced_ntmp 
	xtset uid date
}

// create wave-level versions of stuff:

// Wave aggregate these things:
by  uid:      gen erebate_1tmp = sum(erebate)
bys uid wav: egen erebate_wave = max(erebate_1tmp==1)
drop erebate_1tmp

bys uid wave: egen nmo_wv = count(srefmon)
foreach vi in foodstamp housing foodass welfare ui energy pubhousing rentsubsidy {
	by uid wave: egen hh_`vi'_wvmax  = max(hh_`vi'_ind)
	by uid wave: egen hh_`vi'_wvct   = total(hh_`vi'_ind)
	gen hh_`vi'_wvfrac = hh_`vi'_wvct/ nmo_wv
	gen hh_`vi'_wvall  = hh_`vi'_wvct==nmo_wv
}
// pull across the adult-welfare variables
foreach vi of varlist eabmeet-aafday{
	if(substr("`vi'",1,1) == "a"){
		drop `vi'
	}
	else{
		gen `vi'6 = `vi' if wave==6
		by uid : egen `vi'wv6 = max(`vi'6)
		gen `vi'9 = `vi' if wave==9
		by uid : egen `vi'wv9 = max(`vi'9)
		drop `vi'6 `vi'9
	}
}

by uid: egen estimuse_max = max(estimuse)
by uid: egen enotable_max = max(enotable)

save ${sipp_rebate}data/rebatedata_cleaned.dta, replace



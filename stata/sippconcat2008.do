// concatenate data for the 2008 panel, inc rebate, MPC and basic info

global coredat_dir /media/david/DATAext4/data/SIPP
global TMdat_dir   ~/Dropbox/MPCwhose/
global outdir     ~/Dropbox/rebate_timing/data/
global sippset sippsets08	
	
use "${coredat_dir}/${sippset}/set_d", clear
di "------------------------------------------------------------"
di "${sippset}"
di "------------------------------------------------------------"

di "Dropping blank IDs:"
drop if id == ""
// previously, person_time was group(id srefmon wave age)
egen person_time = group(id srefmon wave)

di "Dropping if duplicate (id, srefmon, wave, job) combo"

/*if("${sippset}"=="sippsets90"|"${sippset}"=="sippsets91"|"${sippset}"=="sippsets92"|"${sippset}"=="sippsets93"){
	rename job job_orig
	rename jobid_revised job
	replace job = job_orig if job==. & job_orig<.
}*/
bysort person_time jobno: gen copy_num = _n
keep if copy_num == 1
drop copy_num

bysort person_time: egen hourstie = sd(hours_m)
replace hourstie = (hourstie <= 0) 

bysort person_time: egen num_jobs_t = count(job)

di "Dropping if a row has missing job but a job was given elsewhere for that person and time period"
drop if job == . & num_jobs_t > 0

bysort id job: egen prim_count = sum(primjob)
bysort person_time: egen max_prim_count_t = max(prim_count)

di "Dropping non-primary job if more than one job given in time period and not an hours tie"
drop if num_jobs_t > 1 & primjob != 1 & hourstie != 1

di "Dropping non-primary job if more than one job given in time period and not an hours tie"
drop if num_jobs_t > 1 & primjob == 1 & hourstie == 1 & prim_count < max_prim_count

di "Dropping non-primary job if one has lower earnings than the other"
duplicates tag person_time , gen(duptag)
bysort person_time: drop if duptag>0 & primjob!=1
	
drop person_time num_jobs_t prim_count max_prim_count_t hourstie duptag

saveold "${coredat_dir}/${sippset}/${sippset}ABDFG", replace 
foreach X in a b f g{
	use "${coredat_dir}/${sippset}/set_`X'", clear
	di "Dropping blank IDs:"
	drop if id == ""
	egen unique_id = group(id srefmon wave)
	di "Dropping duplcate (id, srefmon, wave) combinations:"
	bysort unique_id: keep if _n == 1
	drop unique_id
	merge 1:1 id srefmon wave age using "${coredat_dir}/${sippset}/${sippset}ABDFG", nogen force
	saveold "${coredat_dir}/${sippset}/${sippset}ABDFG", replace version(13)
}

// merge in the topical modules:
//need to merge in sippp08putm1.dta  sippp08putm2.dta  sippp08putm4.dta  sippp08putm6.dta
use ${coredat_dir}/${sippset}/${sippset}ABDFG, clear

egen uid = group(id)
label var uid "integer-type unique (within panel) person id"
gen date = ym(year,month)
format date %tm

merge m:1 id wave using ${TMdat_dir}/sippp08putm1.dta , gen(_merge_tm1)
merge m:1 id wave using ${TMdat_dir}/sippp08putm2.dta , gen(_merge_tm2) update
merge m:1 id wave using ${TMdat_dir}/sippp08putm4.dta , gen(_merge_tm4) update
merge m:1 id wave using ${TMdat_dir}/sippp08putm6.dta , gen(_merge_tm6) update
merge m:1 id wave using ${TMdat_dir}/sippp08putm9.dta , gen(_merge_tm9) update

drop if date>=.
// duplicates are all from _merge_tm2 updates
xtset uid date


// labor force status ----------------------------
//ESR 1 -- With job entire month, worked all weeks.
//ESR 2 -- With job entire month, missed one or more weeks but not because of a layoff.
//ESR 3 -- With job entire month, missed 1 or more weeks because of layoff.
//ESR 4 -- With job part of month, but not because of layoff or looking for work.
//ESR 5 -- With job part of month, some time spent on layoff or looking for work.
//ESR 6 -- No job in month, spent entire month on layoff or looking for work.
//ESR 7 -- No job in month, spent part of month on layoff or looking for work.
//ESR 8 -- No job in month, no time spent on layoff or looking for work
gen     lfstat = 1 if esr<=5
replace lfstat = 2 if esr>=6 & esr <=7
replace lfstat = 3 if esr>7 & esr< .
label var lfstat "E or U or N"

gen EU = l.lfstat==1 & lfstat==2 if l.lfstat==1
gen EN = l.lfstat==1 & lfstat==3 if l.lfstat==1

gen UE = l.lfstat==2 & lfstat==1 if l.lfstat==2
gen NE = l.lfstat==3 & lfstat==1 if l.lfstat==3

label var EU "Employment to unemployment transition last period"
label var EN "Employment to non-participation transition last period"

label var UE "Unemployment to employment transition last period"
label var NE "Non-participation to employment to transition last period"

gen Eend = year==eyear & month == emonth
gen Estart=year==syear & month == smonth
gen JCstart = job != l.job if lfstat==1 & l.lfstat==1
gen JCend   = job != f.job if lfstat==1 & f.lfstat==1
bysort uid wave: egen JCend_any = max(JCend)
bysort uid wave: egen JCstart_any = max(JCstart)
xtset uid date
gen fEE = (Eend==1 & JCend_any==1) | (Estart==1 & JCstart_any==1) if lfstat==1 & f.lfstat==1
gen EE = l.fEE==1 if lfstat==1 & l.lfstat==1

label var EE "J2J to transition last period, using start-end dates and job numbers"
drop Eend Estart JCstart JCend

// defined following Flaaen, Shapiro Sorkin (AEJ Macro 2019) https://www.aeaweb.org/articles?id=10.1257/mac.20170162
gen displaced_layoff = ersend==1
gen displaced_empclosed = ersend==9 | ersend==10
gen displaced_slackbiz = ersend==13
gen displaced = (displaced_layoff | displaced_empclosed | displaced_slackbiz) // <- we might need to take max over the wave for this... I don't recall if its monthly report or wave-based

compress
label data "MPC SIPP08 dataset" 
saveold "${outdir}${sippset}_MPC", replace version(13)

//quick check
gen delpy = log(tptotinc)- log(l.tptotinc)
gen lerbatamt = log(erbatamt)

reg delpy l.lerbatamt  if erebate==1

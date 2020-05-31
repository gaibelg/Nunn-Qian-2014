qui { // Generic preperations

	version 13.0
	clear all
	capture log close
	set more off
	set matsize 7000
	set scheme s1mono
	
}
*

* Main folder must include the data file "FAid_Final.dta" available at the AER website and the placebo file "placebos.dta"
cd "D:\replication\Replication Files"

log using "FAid_GG.log", replace


************************************
*** Some basic coding & cleaning ***
************************************

use "FAid_Final.dta", clear

tsset obs year // Defining the dataset as a time series panel. 'obs' is a variable of coded countries

qui { // Converting wheat aid measure & production to thousands of tonnes - coefficients are easier to read

	replace wheat_aid=wheat_aid/1000 	
	replace US_wheat_production=US_wheat_production/1000
	replace recipient_wheat_prod=recipient_wheat_prod/1000
	replace recipient_cereals_prod=recipient_cereals_prod/1000
	replace real_usmilaid=real_usmilaid/1000
	replace real_us_nonfoodaid_ecaid=real_us_nonfoodaid_ecaid/1000
	replace non_us_oda_net=non_us_oda_net/1000
	replace non_us_oda_net2=non_us_oda_net2/1000
	replace world_wheat_aid=world_wheat_aid/1000
	replace world_cereals_aid=world_cereals_aid/1000 
	replace non_US_wheat_aid=non_US_wheat_aid/1000 
	replace non_US_cereals_aid=non_US_cereals_aid/1000
	
}
*

gen l_US_wheat_production=l.US_wheat_production // Lagged US wheat production


foreach x of varlist all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec { // Interacting weather variables with the propensity of receiving aid
	drop `x'_faavg
	gen `x'_faavg=`x'*fadum_avg
}
*


/* Creating instruments */

gen instrument=l.US_wheat_production*fadum_avg // Main instrument: interaction of lagged wheat production in the US and the propensity of recieving aid
la var instrument "Baseline interaction instrument: US wheat prod (t-1) x avg food aid prob (1971-2006)"

gen instrument2=l.US_wheat_production // Secondary instrument: the more basic version - lagged US wheat production only.
la var instrument2 "Alternative unteracted instrument: US wheat production (t-1)"


/* Restricting the sample to relevant years */

keep if inrange(year, 1971, 2006)

replace year=year-1970 // This won't harm the regressions, but will allow us to use quadratic time trends.

/* Coding regions */

gen cont=.
replace cont=1 if wb_region=="East Asia & Pacific"
replace cont=2 if wb_region=="Europe & Central Asia"
replace cont=3 if wb_region=="Latin America & Caribbean"
replace cont=4 if wb_region=="Middle East & North Africa"
replace cont=5 if wb_region=="South Asia"
replace cont=6 if wb_region=="Sub-Saharan Africa"


/* Generating required controls and interations */

gen USA_ln_income = ln(USA_rgdpch) // USA log GDP

gen US_income_fadum_avg=USA_ln_income*fadum_avg // Interacting USA log GDP and the propensity of recieving aid

tab year, gen(ydum) // Year dummies
tab cont, gen(contdum) // Region dummies
tab risocode, gen(cdum) // Country dummies

* Intreactions of region and years
forval x=1/36{
	gen cont1_y`x'=contdum1*ydum`x'
}
forval x=1/36{
	gen cont2_y`x'=contdum2*ydum`x'
}
forval x=1/36{
	gen cont3_y`x'=contdum3*ydum`x'
}
forval x=1/36{
	gen cont4_y`x'=contdum4*ydum`x'
}
forval x=1/36{
	gen cont5_y`x'=contdum5*ydum`x'
}
forval x=1/36{
	gen cont6_y`x'=contdum6*ydum`x'
}
*


forval x=1/36{ // Interacting average recipient per capita cereals production in recepient countries and years
	gen rcereal_y`x'=recipient_pc_cereals_prod_avg*ydum`x'
}
*

forval x=1/36{ // Interacting average total cereal import quantity in recepient countries and years
	gen rimport_y`x'=cereal_pc_import_quantity_avg*ydum`x'
}
*

forval x=1/36{ // Interacting average per capita real U.S. non-food aid ($) and years
	gen usec_y`x'=real_us_nonfoodaid_ecaid_avg*ydum`x'
}
*

forval x=1/36{ // Interacting average per capita real U.S. military aid and years
	gen usmil_y`x'=real_usmilaid_avg*ydum`x'
}
*

bysort risocode: egen ln_rgdpch_avg=mean(ln_rgdpch) // Avereage log GDP of countries across time

forval x=1/36{ // Interacting Avereage log GDP of countries with years

	gen gdp_y`x'=ln_rgdpch_avg*ydum`x'
	
}
*

gen oil_fadum_avg=oil_price_2011_USD*fadum_avg // Interacting oil prices (base year - 2011) and the propensity of recieving aid

gen US_democ_pres_fadum_avg=US_president_democ*fadum_avg  // Interacting dummy variable for a USA democratic president and the propensity of recieving aid


sort risocode year

save "sample_0.dta", replace


/* Defining appropriate macro of controls for convinient and compact representations of the regressions below. 
The local is defined according to the researchers' specifications. */

local baseline_controls "oil_fadum_avg US_income_fadum_avg US_democ_pres_fadum_avg gdp_y2-gdp_y36 usmil_y2-usmil_y36 usec_y2-usec_y36 rcereal_y2-rcereal_y36 rimport_y2-rimport_y36 all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec all_Precip_jan_faavg-all_Precip_dec_faavg all_Temp_jan_faavg-all_Temp_dec_faavg"


**************************
*** Summary Statistics ***
**************************

/* Generating in-sample indicator so that all specifications have the same number of observations. 
That is, the researchers chose to restrict all their data by the most restrictive specification */

qui: xi: ivreg2 intra_state (wheat_aid=instrument) `baseline_controls' i.risocode i.year*i.wb_region, cluster(risocode)
gen in_sample=1 if e(sample)==1

*** Table 1 ***
***************
* Replication of table 1 in the article 
estpost summarize any_war wheat_aid fadum_avg instrument2 if in_sample==1
esttab using table1.rtf, cells("count mean sd") noobs title("Descriptive statistics") alignment(c) ///
sfmt(%9.0g) replace 

preserve
	
	replace year=year+1970 
	
	*** Figure 1 ***
	****************
	* Amount of aid as a function of lagged US production over time

	collapse (mean) wheat_aid (mean) instrument2 (mean) any_war (mean) US_wheat_production, by(year)
	
	gr two (scatter wheat_aid instrument2, mlabel(year)) (lfit wheat_aid instrument2, lcolor(red)), legend(off) ///
	title("") xtitle("Lagged US wheat production (mil MT)") ///
	ytitle("Amount of wheat aid (mil MT)") scale(0.9)
	
	graph export fig1.png, replace
	
	*** Figure 2 ***
	****************
	* Wheat production time trend

	gr two (scatter US_wheat_production year) (qfit US_wheat_production year, lcolor(red)), legend(off) ///
	title("") xtitle("year") ///
	ytitle("US wheat production (mil MT)") scale(0.9)
	
	graph export fig2.png, replace
	
	* Avarage conflict incidence time trend
	
	gr two (scatter any_war year) (qfit any_war year, lcolor(red)), legend(off) ///
	title("") xtitle("year") ///
	ytitle("Avarage conflict incidence") scale(0.9)
	
	graph export fig3.png, replace
		
restore


*** Table 2 ***
***************
**************************************************************
*** Replication of collomn 5 in table 2: Baseline OLS & IV ***
**************************************************************
 
cap eststo clear

*** OLS Estimates ***

eststo M1: reg any_war wheat_aid `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs)

*** Reduced Form ***
preserve

	replace any_war=any_war*1000 // Just for convinience

	eststo M2: reg any_war instrument `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs)


restore

*** Second Stage of IV ***

eststo M3: ivreg2 any_war (wheat_aid=instrument) `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs) ffirst

*** First stage of IV *****

eststo M4: reg wheat_aid instrument `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs)

esttab using table2.rtf, b(%12.5f) se(%12.5f) long nobase noomitted stats(N, labels(Observations)) star(* 0.1 ** 0.05 *** 0.01) ///
keep(wheat_aid instrument) mtitles("OLS" "Reduced Form" "Second Stage" "First Stage") title("Baseline specification - replication of colomn 5 of table 2") ///
nonumbers addnotes("Standard errors are clustered at the country level") alignment(c) replace

eststo clear

*** Table 2 ***
***************
*****************************************************
*** Replication of collomn 5 in table 3: OLS & IV ***
*****************************************************

/* Generating analogous control variables variables; linear time trends interactions */
gen x1=year*real_usmilaid_avg
gen x2=year*real_us_nonfoodaid_ecaid_avg
gen x3=year*ln_rgdpch_avg
gen x4=year*recipient_pc_cereals_prod_avg
gen x5=year*cereal_pc_import_quantity_avg

local controls "oil_price_2011_USD USA_ln_income US_president_democ all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec x1 x2 x3 x4 x5"

/* Generating in-sample indicator so that all specifications have the same number of observations */
qui: ivreg2 any_war (wheat_aid=instrument2) `controls' i.obs c.year#i.cont, cluster(obs)
gen in_sample2=1 if e(sample)==1


*** OLS Estimates ***

eststo M1: reg any_war wheat_aid `controls' i.obs c.year#i.cont if in_sample2==1, cluster(obs)

*** Reduced Form ***
preserve

	replace any_war=any_war*1000 // Just for convinience

	eststo M2: reg any_war instrument2 `controls' i.obs c.year#i.cont if in_sample2==1, cluster(obs)

restore

*** Second Stage of IV ***

eststo M3: ivreg2 any_war (wheat_aid=instrument2) `controls' i.obs c.year#i.cont if in_sample2==1, cluster(obs) ffirst

*** First stage of IV *****

eststo M4: reg wheat_aid instrument2 `controls' i.obs c.year#i.cont if in_sample2==1, cluster(obs)

esttab using table3.rtf, b(%12.5f) se(%12.5f) long nobase noomitted stats(N, labels(Observations)) star(* 0.1 ** 0.05 *** 0.01) ///
keep(wheat_aid instrument2) mtitles("OLS" "Reduced Form" "Second Stage" "First Stage") ///
title("Uninteraucted specification - replication of colomn 5 of table 3") nonumbers  ///
addnotes("Standard errors are clustered at the country level") alignment(c) replace

eststo clear







*****************
*** Extension ***
*****************






*** Figure 3 ***
****************
// Plotting the problem of combining the interaction of lagged US wheat production and D with time-region FE //

	
bysort cont year: egen fe_ir=mean(any_war) // Calculating the region-time fixed effect values

gen deviations=any_war -fe_ir // Deviations from the region-time fixed effect

bysort fadum_avg year: egen mean_deviations=mean(deviations) // For every level of propensity to recieve aid (D) we calculate the average deviation


preserve
	
	replace year=year+1970
	
	collapse (mean) mean_deviations, by(fadum_avg year)

	xtile quintile_fadum_avg = fadum_avg, nq(5) // Deviding propensity to recieve aid to five level groups

	bysort quintile_fadum_avg year: egen mean_deviations_quintile_fadum=mean(mean_deviations) // In each group obtaining the mean by years
	
	gr two (scatter mean_deviations_quintile_fadum year, msize(small) mcolor(blue)) ///
	(qfit mean_deviations_quintile_fadum year, lcolor(red)), by(quintile_fadum_avg) ///
	legend(off) title("") xtitle("year") ytitle("Mean conflict net of region-time FE") scale(0.9)
	
	graph export fig4.png, replace
	
restore


*** Table 4 ***
***************
// Presenting in each region the disitribution of D (propensity to recieve aid) //

xtile levels_of_fadum = fadum_avg, nq(5) // Deviding propensity to recieve aid to five level groups

preserve // First, we filter just country included in the regressions. Then, calculating the share of each D catagory amond included countries for each region.

	local baseline_controls "oil_fadum_avg US_income_fadum_avg US_democ_pres_fadum_avg gdp_y2-gdp_y36 usmil_y2-usmil_y36 usec_y2-usec_y36 rcereal_y2-rcereal_y36 rimport_y2-rimport_y36 all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec all_Precip_jan_faavg-all_Precip_dec_faavg all_Temp_jan_faavg-all_Temp_dec_faavg"

	qui reg any_war wheat_aid `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs)

	keep if e(sample)==1

	collapse (lastnm) levels_of_fadum (mean) cont, by(obs)
	
	bysort cont: tab levels_of_fadum

restore


*** Table 3 ***
***************
// Estimating the models with different groups or with simple quadratic time trends, with and without controls //

local baseline_controls "oil_fadum_avg US_income_fadum_avg US_democ_pres_fadum_avg gdp_y2-gdp_y36 usmil_y2-usmil_y36 usec_y2-usec_y36 rcereal_y2-rcereal_y36 rimport_y2-rimport_y36 all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec all_Precip_jan_faavg-all_Precip_dec_faavg all_Temp_jan_faavg-all_Temp_dec_faavg"
local controls "oil_price_2011_USD USA_ln_income US_president_democ all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec x1 x2 x3 x4 x5"

local groups levels_of_fadum

{ // Baseline specification:
	

	*** OLS Estimates ***

	eststo M1: reg any_war wheat_aid `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs)

	*** Reduced Form ***
	preserve

		replace any_war=any_war*1000 // Just for convinience

		eststo M2: reg any_war instrument `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs)


	restore

	*** Second Stage of IV ***

	eststo M3: ivreg2 any_war (wheat_aid=instrument) `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs) ffirst

	*** First stage of IV *****

	eststo M4: reg wheat_aid instrument `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs)

	esttab using table4.rtf, b(%12.5f) se(%12.5f) long nobase noomitted stats(N, labels(Observations)) star(* 0.1 ** 0.05 *** 0.01) ///
	keep(wheat_aid instrument) mtitles("OLS" "Reduced Form" "Second Stage" "First Stage") title("Baseline specification - groups by D") ///
	nonumbers addnotes("Standard errors are clustered at the country level") alignment(c) replace

eststo clear	

}
*

{ // Baseline specification without controls:
	

	*** OLS Estimates ***

	eststo M1: reg any_war wheat_aid i.obs i.year#i.`groups' if in_sample==1, cluster(obs)

	*** Reduced Form ***
	preserve

		replace any_war=any_war*1000 // Just for convinience

		eststo M2: reg any_war instrument i.obs i.year#i.`groups' if in_sample==1, cluster(obs)


	restore

	*** Second Stage of IV ***

	eststo M3: ivreg2 any_war (wheat_aid=instrument) i.obs i.year#i.`groups' if in_sample==1, cluster(obs) ffirst

	*** First stage of IV *****

	eststo M4: reg wheat_aid instrument i.obs i.year#i.`groups' if in_sample==1, cluster(obs)

	esttab using table5.rtf, b(%12.5f) se(%12.5f) long nobase noomitted stats(N, labels(Observations)) star(* 0.1 ** 0.05 *** 0.01) ///
	keep(wheat_aid instrument) mtitles("OLS" "Reduced Form" "Second Stage" "First Stage") title("Baseline specification - groups by D without controls") ///
	nonumbers addnotes("Standard errors are clustered at the country level") alignment(c) replace

eststo clear	
}
*

{ // Alternative specification:

	*** OLS Estimates ***

	eststo M1: reg any_war wheat_aid `controls' i.obs c.year##c.year if in_sample2==1, cluster(obs)

	*** Reduced Form ***
	preserve

		replace any_war=any_war*1000 // Just for convinience

		eststo M2: reg any_war instrument2 `controls' i.obs c.year##c.year if in_sample2==1, cluster(obs)

	restore

	*** Second Stage of IV ***

	eststo M3: ivreg2 any_war (wheat_aid=instrument2) `controls' i.obs c.year##c.year if in_sample2==1, cluster(obs) ffirst

	*** First stage of IV *****

	eststo M4: reg wheat_aid instrument2 `controls' i.obs c.year##c.year if in_sample2==1, cluster(obs)

	esttab using table6.rtf, b(%12.5f) se(%12.5f) long nobase noomitted stats(N, labels(Observations)) star(* 0.1 ** 0.05 *** 0.01) ///
	keep(wheat_aid instrument2) mtitles("OLS" "Reduced Form" "Second Stage" "First Stage") ///
	title("Uninteraucted specification - global quadratic time trends") nonumbers  ///
	addnotes("Standard errors are clustered at the country level") alignment(c) replace

eststo clear
}
*	

{ // Alternative specification without controls: 

	*** OLS Estimates ***

	eststo M1: reg any_war wheat_aid i.obs c.year##c.year if in_sample2==1, cluster(obs)

	*** Reduced Form ***
	preserve

		replace any_war=any_war*1000 // Just for convinience

		eststo M2: reg any_war instrument2 i.obs c.year##c.year if in_sample2==1, cluster(obs)

	restore

	*** Second Stage of IV ***

	eststo M3: ivreg2 any_war (wheat_aid=instrument2) i.obs c.year##c.year if in_sample2==1, cluster(obs) ffirst

	*** First stage of IV *****

	eststo M4: reg wheat_aid instrument2 i.obs c.year##c.year if in_sample2==1, cluster(obs)

	esttab using table7.rtf, b(%12.5f) se(%12.5f) long nobase noomitted stats(N, labels(Observations)) star(* 0.1 ** 0.05 *** 0.01) ///
	keep(wheat_aid instrument2) mtitles("OLS" "Reduced Form" "Second Stage" "First Stage") ///
	title("Uninteraucted specification - global quadratic time trends without controls") nonumbers  ///
	addnotes("Standard errors are clustered at the country level") alignment(c) replace

eststo clear
}
*

*** Table A2 ***
****************
// Placebo tests //

replace year=year+1970

merge m:1 year using "placebos.dta", nogen

replace year=year-1970

replace star_wars=star_wars*10^7 // Scaling the varibles in order to ease the reading
replace economics=economics*10^5
replace econometrics=econometrics*10^6

foreach var in star_wars economics econometrics { // Baseline specification placebo checks:
 
	qui{
		
		gen placebo=`var'*fadum_avg
		
		local baseline_controls "oil_fadum_avg US_income_fadum_avg US_democ_pres_fadum_avg gdp_y2-gdp_y36 usmil_y2-usmil_y36 usec_y2-usec_y36 rcereal_y2-rcereal_y36 rimport_y2-rimport_y36 all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec all_Precip_jan_faavg-all_Precip_dec_faavg all_Temp_jan_faavg-all_Temp_dec_faavg"

		*** Reduced Form ***
		preserve

			replace any_war=any_war*1000 // Just for convinience

			eststo M1: reg any_war placebo `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs)


		restore

		*** Second Stage of IV ***

		eststo M2: ivreg2 any_war (wheat_aid=placebo) `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs) ffirst
		
			scalar F_baseline_placebo_`var'= e(widstat)
			
		*** First stage of IV *****

		eststo M3: reg wheat_aid placebo `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs)

		esttab using table8_`var'.rtf, b(%12.5f) se(%12.5f) long nobase noomitted stats(N, labels(Observations)) star(* 0.1 ** 0.05 *** 0.01) ///
		keep(wheat_aid placebo) mtitles("Reduced Form" "Second Stage" "First Stage") title("Baseline specification - `var'") ///
		nonumbers addnotes("Standard errors are clustered at the country level") alignment(c) replace

		eststo clear
		
		drop placebo
	}
	*
	
	dis "F statistic `var' (KP): " F_baseline_placebo_`var'

}
*

*** Table 5 ***
***************
// Monte Carlo - instrumenting with random walks - baseline spcification//

set seed 100

cap snapshot erase _all
snapshot save, label("sample_0")
	

cap postclose Monte_Carlo
postfile Monte_Carlo rep beta_RF_a se_RF_a F_FS_a beta_FS_a se_FS_a beta_IV_a se_IV_a beta_RF_b se_RF_b F_FS_b beta_FS_b se_FS_b beta_IV_b se_IV_b using "MC_Exports", replace

global R = 1000
global rep = 1

while $rep <= $R {
	
	qui { 
		
		bysort year: gen epsilon = rnormal(0,1) // The randomized component
		
		* The DGP is defined by the equation: RW_t = RW_(t-1) + epsilon
			
		gen random_walk_${rep}=0
		replace random_walk_${rep} = epsilon if year==1
		replace random_walk_${rep} = random_walk_${rep}[_n-1] + epsilon if year>1
			
		gen instrument_MC=random_walk_${rep}*fadum_avg // The placebo
			
// Baseline specification: model a //
		
		local baseline_controls "oil_fadum_avg US_income_fadum_avg US_democ_pres_fadum_avg gdp_y2-gdp_y36 usmil_y2-usmil_y36 usec_y2-usec_y36 rcereal_y2-rcereal_y36 rimport_y2-rimport_y36 all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec all_Precip_jan_faavg-all_Precip_dec_faavg all_Temp_jan_faavg-all_Temp_dec_faavg"
			
		// Reduced form //
		reg any_war instrument_MC `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs)
		
			scalar beta_RF_a=_b[instrument_MC]
			scalar se_RF_a=_se[instrument_MC]
		
		// First stage //
		reg wheat_aid instrument_MC `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs)
		
			scalar beta_FS_a=_b[instrument_MC]
			scalar se_FS_a=_se[instrument_MC]
		
		// Second stage // 
		ivreg2 any_war (wheat_aid=instrument_MC) `baseline_controls' i.obs i.year#i.cont if in_sample==1, cluster(obs) ffirst
		
			scalar F_FS_a= e(widstat)
			scalar beta_IV_a=_b[wheat_aid]
			scalar se_IV_a=_se[wheat_aid]
			
// Baseline specification with interacted D groups: model b //

		local baseline_controls "oil_fadum_avg US_income_fadum_avg US_democ_pres_fadum_avg gdp_y2-gdp_y36 usmil_y2-usmil_y36 usec_y2-usec_y36 rcereal_y2-rcereal_y36 rimport_y2-rimport_y36 all_Precip_jan-all_Precip_dec all_Temp_jan-all_Temp_dec all_Precip_jan_faavg-all_Precip_dec_faavg all_Temp_jan_faavg-all_Temp_dec_faavg"
		local groups levels_of_fadum

		// Reduced form //
		reg any_war instrument_MC `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs)
		
			scalar beta_RF_b=_b[instrument_MC]
			scalar se_RF_b=_se[instrument_MC]
		
		// First stage //
		reg wheat_aid instrument_MC `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs)
		
			scalar beta_FS_b=_b[instrument_MC]
			scalar se_FS_b=_se[instrument_MC]
		
		// Second stage // 
		ivreg2 any_war (wheat_aid=instrument_MC) `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs) ffirst
		
			scalar F_FS_b= e(widstat)
			scalar beta_IV_b=_b[wheat_aid]
			scalar se_IV_b=_se[wheat_aid]
			
		post Monte_Carlo (${rep}) (beta_RF_a) (se_RF_a) (F_FS_a) (beta_FS_a) (se_FS_a) (beta_IV_a) (se_IV_a) (beta_RF_b) (se_RF_b) (F_FS_b) (beta_FS_b) (se_FS_b) (beta_IV_b) (se_IV_b)
			
		drop epsilon instrument_MC
			
		global rep=${rep}+1
			
	}
	*
	
	dis "Iterations left: " ${R}-${rep}
	
}
*

postclose Monte_Carlo

snapshot save, label("sample_0_with_random_walks")

use MC_Exports, clear

foreach i in a b { // Generating regressions statistics

	gen T_RF_`i' = beta_RF_`i'/se_RF_`i'
	gen p_val_RF_`i'=round((1-normal(abs(T_RF_`i')))*2,0.001)

	gen T_FS_`i' = beta_FS_`i'/se_FS_`i'
	gen p_val_FS_`i'=round((1-normal(abs(T_FS_`i')))*2,0.001)

	gen T_IV_`i' = beta_IV_`i'/se_IV_`i'
	gen p_val_IV_`i'=round((1-normal(abs(T_IV_`i')))*2,0.001)

}
*

* Results: calculating the proportions of the relevant groups

count if p_val_RF_a<=0.05 
count if p_val_RF_a<=0.01 


count if p_val_IV_a<=0.05 & F_FS_a>=10
count if p_val_IV_a<=0.01 & F_FS_a>=10
count if F_FS_a>=10

count if p_val_RF_b<=0.05 
count if p_val_RF_b<=0.01


count if p_val_IV_b<=0.05 & F_FS_b>=10
count if p_val_IV_b<=0.01 & F_FS_b>=10
count if F_FS_b>=10
 
save "MC_Exports.dta", replace
snapshot save, label("Exports")


*** Table A1 ***
****************
// Robustness check of the specification with interacted D groups //

snapshot restore 1

* Generating new segments variables with different cardinality

xtile levels_of_fadum3 = fadum_avg, nq(3)
xtile levels_of_fadum4 = fadum_avg, nq(4)
xtile levels_of_fadum5 = fadum_avg, nq(5)
xtile levels_of_fadum6 = fadum_avg, nq(6)
xtile levels_of_fadum7 = fadum_avg, nq(7)

forval i=3(1)7 { // Estimations for every number of levels
	
	local groups levels_of_fadum`i'
	
	*** Reduced Form ***
	preserve

		replace any_war=any_war*1000 // Just for convinience

		eststo M2: reg any_war instrument `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs)


	restore

	*** Second Stage of IV ***

	eststo M3: ivreg2 any_war (wheat_aid=instrument) `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs) ffirst

	*** First stage of IV *****

	eststo M4: reg wheat_aid instrument `baseline_controls' i.obs i.year#i.`groups' if in_sample==1, cluster(obs)

	esttab using table9_`i'.rtf, b(%12.5f) se(%12.5f) long nobase noomitted stats(N, labels(Observations)) star(* 0.1 ** 0.05 *** 0.01) ///
	keep(wheat_aid instrument) mtitles("Reduced Form" "Second Stage" "First Stage") title("Baseline specification - groups by D - different cardinilities") ///
	nonumbers addnotes("Standard errors are clustered at the country level") alignment(c) replace

eststo clear	

}
*

log close

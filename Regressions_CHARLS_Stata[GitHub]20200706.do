* =========================================================================
* Status in Quo of Social Participation and Associated Sociodemographic
* 	Characteristics
* 
* Author: Chenxin Tan
*	Master Student in Applied Quantitative Research, New York University
* 	Personal URL: https://github.com/ChenxinTan
*	Email: chenxin.tan@nyu.edu
*
* Data: China Health and Retirement Longitudinal Study (CHARLS)
*	Website: http://charls.pku.edu.cn/index/en.html
*
* Stata Version: SE 15.1
* =========================================================================

* === < Programmming Setting and Load Data > ===

clear all
clear matrix
capture log off
set more off

set seed 2048

global path "/Users/chenxintan/Desktop/New York University/AQA I/finalpaper"
global datadir "${path}/CHARLS2015/"

cd "${path}"

* <! NOTE !> The data used in this .do file was derived from a R script and is for
* 			the purpose of doing regresssions with sampling design conviniently.

use "${datadir}complete_charls2015.dta", clear

compress

* === < Logit Regressions > ===

replace participation = participation - 1
lab define participation 0 "0 No" 1 "1 Yes", modify

replace male = male - 1
lab define male 0 "0 female" 1 "1 male", modify

recode hukou2015 (1 = 0 "0 rural")(2 = 1 "1 urban"), gen(urban)

replace married = married - 1
lab define married 0 "0 unmarried (widowed)" 1 "1 married", modify

global basicvar     male age agesq urban married hh_member
global physicalvar  male age agesq urban married hh_member i.srhealth phyfunction ADLscore IADLscore MMSEscore CESDscore i.lifesatisfy
global sociodemovar male age agesq urban married hh_member i.srhealth phyfunction ADLscore IADLscore MMSEscore CESDscore i.lifesatisfy i.edulv loghh_perincome

logit participation $basicvar [pw = INDV_weight_ad2], vce(cluster communityID)
est store logit_basic

logit participation $physicalvar [pw = INDV_weight_ad2], vce(cluster communityID)
est store logit_physical

logit participation $sociodemovar [pw = INDV_weight_ad2], vce(cluster communityID)
est store logit_sociodemo

est tab logit_basic logit_physical logit_sociodemo, star(0.05 0.01 0.001) b(%8.3f)

outreg2 [logit_basic logit_physical logit_sociodemo] ///
	using "${path}/table_and_figure/LogitRegressions.xls", replace excel ///
	stats(coef se) cti(raw) dec(3) alpha(0.001, 0.01, 0.05)
outreg2 [logit_basic logit_physical logit_sociodemo] ///
	using "${path}/table_and_figure/LogitRegressions.xls", replace excel ///
	stats(coef se) eform cti(odds ratio) dec(3) alpha(0.001, 0.01, 0.05)

* === < Multinomial Logit Regressions > ====

mlogit spclass $basicvar [pw = INDV_weight_ad2], vce(cluster communityID)
est store mlogit_basic

mlogit spclass $physicalvar [pw = INDV_weight_ad2], vce(cluster communityID)
est store mlogit_physical

mlogit spclass $sociodemovar [pw = INDV_weight_ad2], vce(cluster communityID)
est store mlogit_sociodemo

est tab mlogit_basic mlogit_physical mlogit_sociodemo, star(0.05 0.01 0.001) b(%8.3f)

outreg2 [mlogit_basic mlogit_physical mlogit_sociodemo] ///
	using "${path}/table_and_figure/MultinomialRegressions.xls", replace excel ///
	stats(coef se) eform cti(odds ratio) dec(3) alpha(0.001, 0.01, 0.05)



capture log close

* =========================================================================
* Author: Chenxin Tan
* Last Modified: July 06, 2020
* =========================================================================

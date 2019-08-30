*! Tim P Morris | 05sep2019
version 14
* Runs simulation study to produce 
* 1. estimates data
* 2. states data at start of each rep and end of final rep
* Note: requires survsim package to be installed (ssc install survsim)

quietly {
set seed 1
local nsim 2
local nobs 100 500
local gamma1 1    // DGM is exponential (and Weibull)
local gamma2 1.5  // DGM is Weibull (not exponential)
tempname estimates states
postfile `estimates' int(rep_id) int(n_obs) byte(gamma method) float(theta se) using PostEstimates.dta, replace
postfile `states' int(rep_id) int(n_obs) str2000 s1 str2000 s2 str1100 s3 using PostStates.dta, replace

noi _dots 0, title("Simulation running...")
forvalues i = 1/`nsim' {
    * store initial rngstate

    foreach n of local nobs {
        post `states' (`i') (`n') (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))
        clear
        * declare your sample size
        set obs `n'
        * generate a binary treatment group (0/1), with Prob(0.5) of being in each arm 
        gen byte trt = rbinomial(1,0.5)

        * DGM 
        forvalues j=1/2 {
            * Simulate survival times from Weibull, under proportional hazards, with administrative censoring at 5 years
            survsim stime`j' event`j', dist(weibull) lambda(0.1) gamma(`gamma`j'') cov(trt -0.5) maxt(5)
            * Declare the data to be survival data
            stset stime`j', failure(event`j'=1)
            capture streg trt, dist(exp) nohr
                post `estimates' (`i') (`n') (`j') (1) (_b[trt]) (_se[trt])
            capture streg trt, dist(weibull) nohr
                post `estimates' (`i') (`n') (`j') (2) (_b[trt]) (_se[trt])
            capture stcox trt, estimate
                post `estimates' (`i') (`n') (`j') (3) (_b[trt]) (_se[trt])
        }
    }
    noi _dots `i' 0
}
post `states' (`=`nsim'+1') (.) (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))
postclose `estimates'
postclose `states'
}

* Label estimates data and re-save
use PostEstimates, clear
    label variable rep_id "Rep num"
    label variable n_obs "n{sub:obs}"
    label variable gamma "True γ"
    label variable method "Method"
    label variable theta "θᵢ"
    label variable se "SE(θᵢ)"
    label define gammalab 1 "γ=1" 2 "γ=1.5"
        label values gamma gammalab
    label define methodlab 1 "Exponential" 2 "Weibull" 3 "Cox"
        label values method methodlab
    sort rep_id n_obs gamma method
save PostEstimatesClean, replace
list, noobs sepby(rep_id)

* to load your dataset of random number states
use PostStates, replace
* to extract the rngstate and set it for repetition 2 with n_obs==500 (will simulate data for each gamma)
local i 2 // identify rep_id of interest
local n_obs 500 // identify n_obs you want to replicate
local gamma1 1    // DGM is exponential (and Weibull)
local gamma2 1.5  // DGM is Weibull (not exponential)

keep if n_obs == `n_obs' & rep_id == `i'
local statei = s1[1]+s2[1]+s3[1]
clear
set rngstate `statei'
set obs `n_obs'
gen byte trt = rbinomial(1,0.5)
forvalues j=1/2 {
    survsim stime`j' event`j', dist(weibull) lambda(0.1) gamma(`gamma`j'') cov(trt -0.5) maxt(5)
    stset stime`j', failure(event`j'=1)
    streg trt, dist(exp) nohr
}
* Observe that above estimates match identically those produced when the whole lot were run
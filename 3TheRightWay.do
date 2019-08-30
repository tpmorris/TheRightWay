*! Tim P Morris | 05sep2019
version 14
* Runs simulation study to produce 
* 1. estimates data
* 2. states data at start of each rep and end of final rep
* Note: requires survsim package to be installed (ssc install survsim)

set seed 1
local nsim 2
local nobs 100 500
local gamma1 1    // DGM is exponential (and Weibull)
local gamma2 1.5  // DGM is Weibull (not exponential)

* Define program to run one repetition
capture program drop toysim
program define toysim, rclass
    version 14
    syntax [ , nobs(integer 100) gamma(real 1) theta(real -.5) ]

    clear
    set obs `nobs'
    gen byte trt = rbinomial(1,0.5)
    survsim stime event, dist(weibull) lambda(0.1) gamma(`gamma') cov(trt `theta') maxt(5)
    stset stime, failure(event=1)
    streg trt, dist(exponential) nohr
        return scalar loghr_exp = _b[trt]
        return scalar se_exp = _se[trt]
    streg trt, dist(weibull) nohr
        return scalar loghr_wei = _b[trt]
        return scalar se_wei = _se[trt]
    stcox trt, nohr
        return scalar loghr_cox = _b[trt]
        return scalar se_cox = _se[trt]
end

tempname estimates states
postfile `estimates' int(rep_id) int(n_obs) float(gamma) byte(method) float(theta se) using RightEstimates.dta, replace
postfile `states' int(rep_id) int(n_obs) float(gamma) str2000 s1 str2000 s2 str1100 s3 using RightStates.dta, replace

quietly {
noi _dots 0, title("Simulation running...")
forval i = 1/`nsim' {
    foreach n of numlist 100 500 {
        foreach gamma of numlist 1 1.5 {
            post `states' (`i') (`n') (`gamma') (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))
            toysim, nobs(`n') gamma(`gamma')
                post `estimates' (`i') (`n') (`gamma') (1) (r(loghr_exp)) (r(se_exp))
                post `estimates' (`i') (`n') (`gamma') (2) (r(loghr_wei)) (r(se_wei))
                post `estimates' (`i') (`n') (`gamma') (3) (r(loghr_cox)) (r(se_cox))
        }
    }
    noi _dots `i' 0
}
post `states' (`=`nsim'+1') (.) (.) (substr(c(rngstate),1,2000)) (substr(c(rngstate),2001,2000)) (substr(c(rngstate),4001,.))
postclose `estimates'
postclose `states'
}

* Label estimates data and re-save
use RightEstimates, clear
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
save RightEstimatesClean, replace
list, noobs sepby(rep_id)

* to load your dataset of random number states
use RightStates, replace
* to extract the rngstate and set it for repetition 2 with n_obs==500 (will simulate data for each gamma)
local i 2 // identify rep_id of interest
local n_obs 500 // identify n_obs you want to replicate
local gamma 1    // DGM is exponential (and Weibull)

keep if n_obs == `n_obs' & rep_id == `i' & gamma == `gamma'
local statei = s1[1]+s2[1]+s3[1]
clear
set rngstate `statei'
set obs `n_obs'
gen byte trt = rbinomial(1,0.5)
survsim stime`j' event`j', dist(weibull) lambda(0.1) gamma(`gamma') cov(trt -0.5) maxt(5)
stset stime`j', failure(event`j'=1)
streg trt, dist(exp) nohr
* Observe that above estimates match identically those produced when the whole lot were run
*! Tim P Morris | 05sep2019
version 14
* Runs simulation study to produce 
* 1. estimates data
* Note: requires survsim package to be installed (ssc install survsim)

capture program drop toysim
program define toysim, rclass
    version 14
    syntax [ , nobs(integer 100) gamma(real 1) theta(real -.5) ]

    * I want to be able to return strings as well as scalars!
    //return scalar rngstream = c(rngstream)
    //return local state = c(rngstate)
    return scalar nobs = `nobs'
    return scalar gamma = `gamma'

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

foreach n of numlist 100 500 {
    foreach gamma of numlist 1 1.5 {
        if `gamma' == 1.5 local ests SimulateEstimates_nobs_`n'_gamma1_5 // *facepalm* cannot name dataset "...gamma_1.5" because "."
        else local ests SimulateEstimates_nobs_`n'_gamma`gamma'
        simulate n_obs = r(nobs) gamma = r(gamma) ///
            theta_exp = r(loghr_exp) se_exp = r(se_exp)    ///
            theta_wei = r(loghr_wei) se_wei = r(se_wei) ///
            theta_cox = r(loghr_cox) se_cox = r(se_cox) ///
            , reps(2) saving(`ests', replace) seed(1) : toysim, nobs(`n') gamma(`gamma')
    }
}

* Now append, label etc.
use SimulateEstimates_nobs_100_gamma1, clear
    append using SimulateEstimates_nobs_500_gamma1
    append using SimulateEstimates_nobs_100_gamma1_5
    append using SimulateEstimates_nobs_500_gamma1_5

bysort n_obs gamma: generate int rep_id = _n // This is an artificial repetition no because simulate will not post it!
rename *_exp *1
rename *_wei *2
rename *_cox *3
reshape long theta se, i(rep_id n_obs gamma) j(method)

label define methodlab 1 "Exponential" 2 "Weibull" 3 "Cox"
    label values method methodlab
sort rep_id n_obs gamma method


save SimulateEstimatesAll, replace
list
list, noobs sepby(rep_id)

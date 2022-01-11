####################################################################################################
# Example of code written by Jia Wei (Nuffield Department of Medicine, University of Oxford)       #
# for analyses of antibody waning following second vaccination.                                    #
# Accompanying paper/preprint: SARS-CoV-2 anti-spike IgG antibody responses after second dose of   #
# ChAdOx1 or BNT162b2 and correlates of protection in the UK general population.                   #
####################################################################################################

library(tidyverse)
library(brms)


###ChAdOx1

priors=c(
  set_prior("normal(7.5,1)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,1)",coef="age",class="b"),
  set_prior("normal(0,1)",coef="sex1",class="b"),
  set_prior("normal(0,1)",coef="ethnicity1",class="b"),
  set_prior("normal(0,1)",coef="lthc1",class="b"),
  set_prior("normal(0,1)",coef="hcw1",class="b"),
  set_prior("normal(0,1)",coef="IMD",class="b"),
  set_prior("normal(0,1)",coef="dur",class="b"),
  set_prior("normal(0,1)",coef="prior1",class="b"),
  
  set_prior("normal(0,0.1)",coef="time:age",class="b"),
  set_prior("normal(0,0.1)",coef="time:sex1",class="b"),
  set_prior("normal(0,0.1)",coef="time:ethnicity1",class="b"),
  set_prior("normal(0,0.1)",coef="time:lthc1",class="b"),
  set_prior("normal(0,0.1)",coef="time:hcw1",class="b"),
  set_prior("normal(0,0.01)",coef="time:IMD",class="b"),
  set_prior("normal(0,0.1)",coef="time:dur",class="b"),
  set_prior("normal(0,0.1)",coef="time:prior1",class="b"),
  
  set_prior("cauchy(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",coef="Intercept",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",class="sigma"),
  set_prior("lkj_corr_cholesky(2)",class="L")
)

set_inits=function(seed=1){
  set.seed(seed)
  list(
    Intercept=rnorm(1,7.5,1),
    b=c(rnorm(1,-0.01,0.01),#time
        rnorm(1,0,0.1),#age
        rnorm(1,0,0.1),#sex
        rnorm(1,0,0.1),#ethnicity
        rnorm(1,0,0.1),#lthc
        rnorm(1,0,0.1),#hcw
        rnorm(1,0,0.1),#IMD
        rnorm(1,0,0.1),#dur
        rnorm(1,0,0.1),#prior
        
        rnorm(1,0,0.01),#age
        rnorm(1,0,0.01),#sex
        rnorm(1,0,0.01),#ethnicity
        rnorm(1,0,0.01),#lthc
        rnorm(1,0,0.01),#hcw
        rnorm(1,0,0.01),#IMD
        rnorm(1,0,0.01),#dur
        rnorm(1,0,0.01)#prior
    ),
    sigma=runif(1,0,1),
    sd_1=c(runif(1,0,1),runif(1,0,0.02)),
    z_1=matrix(rep(c(7.5,-0.01),100639),2,100639)
  )
}

inits_list=list(
  set_inits(1),
  set_inits(2),
  set_inits(3),
  set_inits(4)
)

fit_az <- brm(formula=log2(assay)|cens(cen1)~
                 1+time*age+time*sex+time*ethnicity+time*lthc+time*hcw+
                 time*IMD+time*dur+time*prior+(1+time|ID),
               data=dfaz,cores=4,family=gaussian(),
               prior = priors,inits=inits_list,
               chains = 4, iter=4000, warmup = 2000, seed=42,
               control = list(adapt_delta=0.95))


####BNT162b2

priors=c(
  set_prior("normal(9.5,1)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,1)",coef="age",class="b"),
  set_prior("normal(0,1)",coef="sex1",class="b"),
  set_prior("normal(0,1)",coef="ethnicity1",class="b"),
  set_prior("normal(0,1)",coef="lthc1",class="b"),
  set_prior("normal(0,1)",coef="hcw1",class="b"),
  set_prior("normal(0,1)",coef="IMD",class="b"),
  set_prior("normal(0,1)",coef="dur",class="b"),
  set_prior("normal(0,1)",coef="week31",class="b"),
  set_prior("normal(0,1)",coef="prior1",class="b"),
  
  set_prior("normal(0,0.1)",coef="time:age",class="b"),
  set_prior("normal(0,0.1)",coef="time:sex1",class="b"),
  set_prior("normal(0,0.1)",coef="time:ethnicity1",class="b"),
  set_prior("normal(0,0.1)",coef="time:lthc1",class="b"),
  set_prior("normal(0,0.1)",coef="time:hcw1",class="b"),
  set_prior("normal(0,0.01)",coef="time:IMD",class="b"),
  set_prior("normal(0,0.1)",coef="time:dur",class="b"),
  set_prior("normal(0,0.1)",coef="time:week31",class="b"),
  set_prior("normal(0,0.1)",coef="time:prior1",class="b"),
  
  set_prior("cauchy(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",coef="Intercept",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",class="sigma"),
  set_prior("lkj_corr_cholesky(2)",class="L")
)


set_inits=function(seed=1){
  set.seed(seed)
  list(
    Intercept=rnorm(1,9.5,0.2),
    b=c(rnorm(1,-0.01,0.01),#time
        rnorm(1,0,0.1),#age
        rnorm(1,0,0.1),#sex
        rnorm(1,0,0.1),#ethnicity
        rnorm(1,0,0.1),#lthc
        rnorm(1,0,0.1),#hcw
        rnorm(1,0,0.1),#IMD
        rnorm(1,0,0.1),#dosing
        rnorm(1,0,0.1),#3week
        rnorm(1,0,0.1),#prior
        
        rnorm(1,0,0.01),#age
        rnorm(1,0,0.01),#sex
        rnorm(1,0,0.01),#ethnicity
        rnorm(1,0,0.01),#lthc
        rnorm(1,0,0.01),#hcw
        rnorm(1,0,0.01),#IMD
        rnorm(1,0,0.01),#dosing
        rnorm(1,0,0.01),#3week
        rnorm(1,0,0.01)#prior
    ),
    sigma=runif(1,0,1),
    sd_1=c(runif(1,0,1),runif(1,0,0.02)),
    z_1=matrix(rep(c(9.5,-0.01),55053),2,55053)
  )
}

inits_list=list(
  set_inits(1),
  set_inits(2),
  set_inits(3),
  set_inits(4)
)

fit_pf <- brm(formula=log2(assay)|cens(cen1)~
                 1+time*age+time*sex+time*ethnicity+time*lthc+time*hcw+time*IMD+
                 time*week3+time*dur+time*prior+(1+time|ID),
               data=dfpf,cores=4,family=gaussian(),
               prior = priors,
               chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits_list,
               control = list(adapt_delta=0.95))



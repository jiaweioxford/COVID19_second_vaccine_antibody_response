####################################################################################################
# Example of code written by Koen Pouwels (Health Economics Research Centre, University of Oxford) #
# for analyses of correlates of protection against SARS-Cov-2 infection.                           #
# Accompanying paper/preprint: SARS-CoV-2 anti-spike IgG antibody responses after second dose of   #
# ChAdOx1 or BNT162b2 and correlates of protection in the UK general population.                   #
####################################################################################################

library(mgcv)
library(dplyr)

# Not vaccinated vs reference level:
gam_model_lag21_59_rounded_vactype <- bam(result_mk ~ s(assay_21_59,k=160, by=vaccine_type) +
                                          
                                            te(study_day, age_at_visit, k=c(20,20), by=region_model) +
                                          
                                            vaccine_type +
                                            ethnicity_wo + sex + 
                                            region_model + 
                                            imd_samp + multigen + hhsizegroup + rural_urban_class +
                                            ever_care_home_worker + patient_facing_clean_ever + 
                                            ever_personfacing_socialcare + ever_lthc +
                                            visit_freq + smoke_now 
                                            contact_hospital + contact_carehome
                                          , data = data_all, 
                                          family = binomial,
                                          method="fREML",
                                          discrete=TRUE,
                                          nthreads=8
)


gam_model_lag21_59_rounded_vactype$aic

newdata = data_all[1,] # if first row at reference level for everything, otherwise change this as needed
newdata <- rbind(newdata,newdata)
newdata$vaccine_type <- "No Vaccine"
newdata$assay_21_59_round <- c(1,8)
newdata <- newdata%>%dplyr::select(assay_21_59, 
                                   study_day, age_at_visit, vaccine_type, 
                                   ethnicity_wo , sex , region_model , imd_samp , multigen , 
                                   hhsizegroup , rural_urban_class , ever_care_home_worker , 
                                   patient_facing_clean_ever , ever_personfacing_socialcare , 
                                   ever_lthc , visit_freq , smoke_now , contact_hospital , contact_carehome)

Xp <- predict(gam_model_lag21_59_rounded_vactype, newdata = newdata, type="lpmatrix")

getor_vaccine_type <- function(newdata = newdata, contrast = c(2,8), model=model) {
  newdata$assay_21_59 <- c(contrast[1],contrast[2])
  Xp <- predict(model, newdata = newdata, type="lpmatrix")
  diff <- (Xp[2,]-Xp[1,])
  dly <- t(diff)%*%coef(model)
  se.dly <- sqrt(t(diff)%*%vcov(model)%*%diff)
  c(exp(dly), exp(dly -2*se.dly), exp(dly + 2*se.dly))
}

dat_ors_unvac <- data.frame(assay_round = seq(3,800,1))
dat_ors_unvac$or <- NA
dat_ors_unvac$ll <- NA
dat_ors_unvac$ul <- NA

getor_vaccine_type(newdata = newdata,  contrast = c(2, 3),
                   model=gam_model_lag21_59_rounded_vactype)


for (i in 1:nrow(dat_ors_unvac)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(2,dat_ors_unvac$assay_round[i]),
                             model = gam_model_lag21_59_rounded_vactype)
  dat_ors_unvac$or[i] <- temp[1]
  dat_ors_unvac$ll[i] <- temp[2]
  dat_ors_unvac$ul[i] <- temp[3]
  
  
}

dat_ors_unvac

###########################################################################################################
###############################################################################################
# AZ no prior inf vs not vaccinated: 

summary(data_all$vaccine_type)
#### get odds ratios out:
# data at which predictions to be compared (don't bother about covariates, 
# only vaccination status really matters)
newdata = data_all[1,]
newdata <- rbind(newdata,newdata)
newdata$vaccine_type <- c("No Vaccine","vaccine_AZ_noinf")
newdata$assay_21_59 <- c(1,8)
newdata <- newdata%>%dplyr::select(assay_21_59, 
                                   study_day, age_at_visit, vaccine_type, 
                                   ethnicity_wo , sex , region_model , imd_samp , multigen , 
                                   hhsizegroup , rural_urban_class , ever_care_home_worker , 
                                   patient_facing_clean_ever , ever_personfacing_socialcare , 
                                   ever_lthc , visit_freq , smoke_now , contact_hospital , contact_carehome)


Xp <- predict(gam_model_lag21_59_rounded_vactype, newdata = newdata, type="lpmatrix")


dat_ors_AZnoinfvsunvac <- data.frame(assay_round = seq(3,800,1))
dat_ors_AZnoinfvsunvac$or <- NA
dat_ors_AZnoinfvsunvac$ll <- NA
dat_ors_AZnoinfvsunvac$ul <- NA

getor_vaccine_type(newdata = newdata,  contrast = c(2, 800),
                   model=gam_model_lag21_59_rounded_vactype)

for (i in 1:nrow(dat_ors_AZnoinfvsunvac)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(2,dat_ors_AZnoinfvsunvac$assay_round[i]),
                             model = gam_model_lag21_59_rounded_vactype)
  dat_ors_AZnoinfvsunvac$or[i] <- temp[1]
  dat_ors_AZnoinfvsunvac$ll[i] <- temp[2]
  dat_ors_AZnoinfvsunvac$ul[i] <- temp[3]
  
}

dat_ors_AZnoinfvsunvac

#####################################################################################################
###################################################################################################################
# Pfizer no prior inf vs unvaccinated: 

#### get odds ratios out:
# data at which predictions to be compared (don't bother about covariates, 
# only vaccination status really matters)
newdata = data_all[1,]
newdata <- rbind(newdata,newdata)
newdata$vaccine_type <- c("No Vaccine","vaccine_PF_noinf")
newdata$assay_21_59 <- c(1,8)
newdata <- newdata%>%dplyr::select(assay_21_59, 
                                   study_day, age_at_visit, vaccine_type, 
                                   ethnicity_wo , sex , region_model , imd_samp , multigen , 
                                   hhsizegroup , rural_urban_class , ever_care_home_worker , 
                                   patient_facing_clean_ever , ever_personfacing_socialcare , 
                                   ever_lthc , visit_freq , smoke_now , contact_hospital , contact_carehome)

summary(gam_model_lag21_59_rounded_vactype)
Xp <- predict(gam_model_lag21_59_rounded_vactype, newdata = newdata, type="lpmatrix")



dat_ors_PFnoinfvsunvac <- data.frame(assay_round = seq(3,800,1))
dat_ors_PFnoinfvsunvac$or <- NA
dat_ors_PFnoinfvsunvac$ll <- NA
dat_ors_PFnoinfvsunvac$ul <- NA


for (i in 1:nrow(dat_ors_PFnoinfvsunvac)) {
  temp <- getor_vaccine_type(newdata=newdata,  contrast=c(output.small$novacc_2.5perc,dat_ors_PFnoinfvsunvac$assay_round[i]),
                             model = gam_model_lag21_59_rounded_vactype)
  dat_ors_PFnoinfvsunvac$or[i] <- temp[1]
  dat_ors_PFnoinfvsunvac$ll[i] <- temp[2]
  dat_ors_PFnoinfvsunvac$ul[i] <- temp[3]
  
}

dat_ors_PFnoinfvsunvac
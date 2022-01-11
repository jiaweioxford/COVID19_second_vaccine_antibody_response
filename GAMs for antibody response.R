####################################################################################################
# Example of code written by Jia Wei (Nuffield Department of Medicine, University of Oxford)       #
# for analyses of antibody response following first and second vaccination.                        #
# Accompanying paper/preprint: SARS-CoV-2 anti-spike IgG antibody responses after second dose of   #
# ChAdOx1 or BNT162b2 and correlates of protection in the UK general population.                   #
####################################################################################################

library(tidyverse)
library(mgcv)

##ChAdOx1, no prior infection
xrange <- range(df_az$time)
yrange <- range(as.numeric(df_az$age))
zrange <- range(df_az$dur)

xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                  seq(from = yrange[1], to = yrange[2],1),
                  seq(from = zrange[1], to = zrange[2],1))  
colnames(xz) <- c("time", "age","dur")
plot_data <- data.frame(xz) 

plot_data<-plot_data %>% 
  filter(time2>= -14-dur)

m_az <- bam(log10(assay)~te(time,age,dur,bs="bs",k=c(20,10,10)),
                method = "fREML", data = df_az, discrete = T, nthreads = 12)

# predicing: 
pred <- predict(m_az, newdata=plot_data, se.fit=TRUE, discrete=F,
                newdata.guaranteed=TRUE)
plot_data$prediction <- pred$fit
plot_data$se <- pred$se.fit
plot_data$ll <- (pred$fit-1.96*pred$se.fit)
plot_data$ul <- (pred$fit+1.96*pred$se.fit)


##BNT162b2, no prior infection

xrange <- range(df_pf$time)
yrange <- range(as.numeric(df_pf$age))
zrange <- range(as.numeric(df_pf$dur))

xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                  seq(from = yrange[1], to = yrange[2],1),
                  seq(from = zrange[1], to = zrange[2],1))  
colnames(xz) <- c("time", "age","dur")
plot_data <- data.frame(xz) 

plot_data<-plot_data %>% 
  filter(time>= -14-dur)

m_pf <- bam(log10(assay)~te(time,age,dur,bs="bs",k=c(25,10,10)),method = "fREML", 
                data = df_pf,discrete = T, nthreads = 12)

# predicing: 
pred <- predict(m_pf, newdata=plot_data, se.fit=TRUE, discrete=F,
                newdata.guaranteed=TRUE)
plot_data$prediction <- pred$fit
plot_data$se <- pred$se.fit
plot_data$ll <- (pred$fit-1.96*pred$se.fit)
plot_data$ul <- (pred$fit+1.96*pred$se.fit)

#BNT162b2, 3-week, no prior infction

xrange <- range(df_pf2$time)
yrange <- range(as.numeric(df_pf2$age))

xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                  seq(from = yrange[1], to = yrange[2],1))  
colnames(xz) <- c("time", "age")
plot_data2 <- data.frame(xz) 

m_pf2 <- bam(log10(assay)~te(time,age,bs="bs",k=c(30,10)),method = "fREML", 
                 data = df_pf2,discrete = T, nthreads = 12)

# predicing: 
pred <- predict(m_pf2, newdata=plot_data2, se.fit=TRUE, discrete=F,
                newdata.guaranteed=TRUE)
plot_data2$prediction <- pred$fit
plot_data2$se <- pred$se.fit
plot_data2$ll <- (pred$fit-1.96*pred$se.fit)
plot_data2$ul <- (pred$fit+1.96*pred$se.fit)

plot_data2$dur=21

plot_data=rbind(plot_data,plot_data2)



#####ChAdOx1, prior infection
xrange <- range(df_az_p$time)
yrange <- range(as.numeric(df_az_p$age))
zrange <- range(df_az_p$dur)

xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                  seq(from = yrange[1], to = yrange[2],1),
                  seq(from = zrange[1], to = zrange[2],1))  
colnames(xz) <- c("time", "age","dur")
plot_data <- data.frame(xz) 

plot_data<-plot_data %>% 
  filter(time>= -90-dur)

m_az_p <- bam(log10(assay)~te(time,age,dur,bs="bs",k=c(30,10,10)),
                method = "fREML", data = df_az_p, discrete = T, nthreads = 12)

# predicing: 
pred <- predict(m_az_p, newdata=plot_data, se.fit=TRUE, discrete=F,
                newdata.guaranteed=TRUE)
plot_data$prediction <- pred$fit
plot_data$se <- pred$se.fit
plot_data$ll <- (pred$fit-1.96*pred$se.fit)
plot_data$ul <- (pred$fit+1.96*pred$se.fit)


#####BNT162b2, prior infection
xrange <- range(df_pf_p$time)
yrange <- range(as.numeric(df_pf_p$age))
zrange <- range(as.numeric(df_pf_p$dur))

xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                  seq(from = yrange[1], to = yrange[2],1),
                  seq(from = zrange[1], to = zrange[2],1))  
colnames(xz) <- c("time", "age","dur")
plot_data <- data.frame(xz) 

plot_data<-plot_data %>% 
  filter(time>= -90-dur)

m_pf_p <- bam(log10(assay)~te(time,age,dur,bs="bs",k=c(30,10,10)),method = "fREML", 
                data = df_pf_p,discrete = T, nthreads = 12)

# predicing: 
pred <- predict(m_pf_p, newdata=plot_data, se.fit=TRUE, discrete=F,
                newdata.guaranteed=TRUE)
plot_data$prediction <- pred$fit
plot_data$se <- pred$se.fit
plot_data$ll <- (pred$fit-1.96*pred$se.fit)
plot_data$ul <- (pred$fit+1.96*pred$se.fit)


#BNT162b2, 3-week, prior infction

xrange <- range(df_pf_p2$time)
yrange <- range(as.numeric(df_pf_p2$age))

xz <- expand.grid(seq(from = xrange[1], to = xrange[2],1), 
                  seq(from = yrange[1], to = yrange[2],1))  
colnames(xz) <- c("time", "age")
plot_data2 <- data.frame(xz) 

m_pf_p2 <- bam(log10(assay)~te(time,age,bs="bs",k=c(30,10)),method = "fREML", 
                 data = df_pf_p2,discrete = T, nthreads = 12)

# predicing: 
pred <- predict(m_pf_p2, newdata=plot_data2, se.fit=TRUE, discrete=F,
                newdata.guaranteed=TRUE)
plot_data2$prediction <- pred$fit
plot_data2$se <- pred$se.fit
plot_data2$ll <- (pred$fit-1.96*pred$se.fit)
plot_data2$ul <- (pred$fit+1.96*pred$se.fit)

plot_data2$dur=21

plot_data=rbind(plot_data,plot_data2)


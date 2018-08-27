http://www.metafor-project.org/doku.php/tips:multiple_factors_interactions
http://www.metafor-project.org/doku.php/analyses:konstantopoulos2011
http://www.wvbauer.com/lib/exe/fetch.php/talks:2016_viechtbauer_goettingen_multilevel_multivariate_network_ma.pdf
https://stats.stackexchange.com/questions/117497/comparing-between-random-effects-structures-in-a-linear-mixed-effects-model

setwd('C:/Users/Vernon/Dropbox/Work/SEEC/Research/Review_seasonality_grasses/')

library(readODS)
library(ggplot2)
library(tidyr)

dat = read_ods('Meta-analysis/Meta-analysis_Q_season.ods', sheet=1)

#Get additional covariates:
library(raster)
warmtocold = raster('C:/Users/Vernon/Dropbox/GIS/Earth_observations/Climate/CHELSA/Seasonality/warmToCold.tif')
dat$warmtocold = raster::extract(warmtocold, cbind(dat$lon,dat$lat))
continents = shapefile('C:/Users/Vernon/Dropbox/GIS/Mapping_schemes/TDWG_political_maps/Level1/level1.shp')
coords = dat[,c('lon','lat')]
coordinates(coords) =~lon+lat
dat$continent = over(coords, continents)$ CONTINENT_

#Calculate response (Hedge's d) for control and treatment biomass:
library(metafor)
d = escalc(measure="SMD", m1i=Biomass_treatment, m2i=Biomass_control, sd1i=SD_treatment, sd2i=SD_control, n1i=N_treatment,
           n2i=N_control, data=dat, var.names=c("HedgeD","d_var"), digits=4)
# #Calculate response (Hedge's d) for control and treatment biomass:
# dat$d = c(dat$Biomass_treatment - dat$Biomass_control) / 
#   sqrt(c(c(c(dat$N_treatment-1)*dat$SD_treatment^2 + c(dat$N_control-1)*dat$SD_control^2) / 
#            c(dat$N_treatment + dat$N_control - 2))) *
#   c(1 - c(3/c(4*c(dat$N_treatment + dat$N_control -2)-1)))

#Exclude Niu et al. 2016
d = d[d$Ref!='Niu et al. 2016',]
# d = d[d$Ref!='Sullivan et al. 2016',]

#Summary statistics:
length(levels(as.factor(d$Study)))
length(levels(as.factor(d$TreatmentID)))


#Standardise study length
sl = tapply(d$Study_length, INDEX=d$TreatmentID, function(x) length(unique(x)))
sldat = data.frame(TreatmentID=names(sl), study_length=sl)
ggplot(data=sldat, aes(sldat$study_length)) + 
  geom_histogram() +
  scale_x_continuous(breaks=seq(0,10,1)) +
  labs(x="Study length", y="Count")
summary(as.factor(sldat$study_length))
d$Study_length_std = scale(d$Study_length)

hist(d$HedgeD)

#Check correlation among predictors 
library(GGally)
ggpairs(d[,c('Exp_perc_summer','Exp_perc_autumn','Exp_perc_winter','Exp_perc_spring','Exp_perc_year','warmtocold')])


names(d)

#Intercept-only model = overall mean of effect sizes
m1.1 = rma.mv(yi=HedgeD, V=d_var, random= ~1|Study/TreatmentID/Study_length_std, struct="CS", method="REML", digits=4, data=d)
summary(m1.1)

m1.2 = rma.mv(yi=HedgeD, V=d_var, random = list(~1| Study_length_std, ~1| TreatmentID, ~1| Study), struct="CS", method="REML", digits=4, data=d)
summary(m1.2)

m1.2 = rma.mv(yi=HedgeD, V=d_var, mods = ~ factor(Study_length), random = ~Study_length|Study/TreatmentID, struct="HAR", method="REML", digits=4, data=d)
summary(m1.2)
rma.mv(yi, V, mods = ~ factor(time) - 1, data = dat.long,
       random = ~ time | study, struct = "HAR")


mod = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_summer + Exp_perc_autumn + Exp_perc_winter + Exp_perc_spring,
              random= ~1|Study/TreatmentID/Study_length_std,struct="CS",method="REML",digits=4, data=d)
summary(mod)


profile(m1.1, sigma2=1)
profile(m1.1, sigma2=2)
profile(m1.1, sigma2=3)

profile(m1.2, sigma2=1)
profile(m1.2, sigma2=2)
profile(m1.2, sigma2=3)

# Build a 3-level model without within-study variance.
modelnovar2 <- rma.mv(yi=HedgeD, V=d_var, random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d,
                      sigma2=c(0,NA,NA), tdist=TRUE)
anova(m1.1,modelnovar2)

# Build a 3-level model without within-treatment variance.
modelnovar3 <- rma.mv(yi=HedgeD, V=d_var, random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d,
                      sigma2=c(NA,0,NA), tdist=TRUE)
anova(m1.1,modelnovar3)

# Build a 3-level model without within-year variance.
modelnovar4 <- rma.mv(yi=HedgeD, V=d_var, random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d,
                      sigma2=c(NA,NA,0), tdist=TRUE)
anova(m1.1,modelnovar4)

# Determining how the total variance is distributed over the
# three levels of the meta-analytic model;
# Print the results in percentages on screen.
n <- nrow(d)
list.inverse.variances <- 1 / (d$d_var)
sum.inverse.variances <- sum(list.inverse.variances)
squared.sum.inverse.variances <- (sum.inverse.variances) ^ 2
list.inverse.variances.square <- 1 / (d$d_var^2)
sum.inverse.variances.square <-
  sum(list.inverse.variances.square)
numerator <- (n - 1) * sum.inverse.variances
denominator <- squared.sum.inverse.variances -
  sum.inverse.variances.square
estimated.sampling.variance <- numerator / denominator
I2_1 <- (estimated.sampling.variance) / (m1.1$sigma2[1]
                                         + m1.1$sigma2[2] + estimated.sampling.variance)
I2_2 <- (m1.1$sigma2[1]) / (m1.1$sigma2[1]
                            + m1.1$sigma2[2] + estimated.sampling.variance)
I2_3 <- (m1.1$sigma2[2]) / (m1.1$sigma2[1]
                            + m1.1$sigma2[2] + estimated.sampling.variance)
amountvariancelevel1 <- I2_1 * 100
amountvariancelevel2 <- I2_2 * 100
amountvariancelevel3 <- I2_3 * 100

amountvariancelevel1
amountvariancelevel2
amountvariancelevel3

#------------------------------------------------------------------------------------
#Seasonal models:
#Summer
mSum = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_summer,
              random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d[d$Exp_perc_summer!=100,])
summary(mSum)
preds = predict(mSum, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
ggplot(preds, aes(x=X.Exp_perc_summer,y=pred)) + 
  geom_point() + #position=position_jitter(width=0.5)
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
#Autumn
mAut = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_autumn,
              random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d[d$Exp_perc_autumn!=100,])
summary(mAut)
#Winter
mWin = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_winter,
              random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d[d$Exp_perc_winter!=100,])
summary(mWin)
#Spring
mSpr = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_spring,
              random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d[d$Exp_perc_spring!=100,])
summary(mSpr)
#All
mSeasons = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_summer + Exp_perc_winter + Exp_perc_autumn + Exp_perc_spring,
                  random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(mSeasons)
forest(mSeasons$beta[2:5], sei=mSeasons$se[2:5], slab=c('Summer change','Autumn change','Winter change','Spring change'),
       measure='Biomass change effect size')
preds = predict(mSeasons, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
ggplot(preds, aes(x=X.Exp_perc_summer,y=pred)) + 
  geom_point() + #position=position_jitter(width=0.5)
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
ggplot(preds, aes(x=X.Exp_perc_autumn,y=pred)) + 
  geom_point() + #position=position_jitter(width=0.5)
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
ggplot(preds, aes(x=X.Exp_perc_winter,y=pred)) + 
  geom_point() + #position=position_jitter(width=0.5)
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
ggplot(preds, aes(x=X.Exp_perc_spring,y=pred)) + 
  geom_point() + #position=position_jitter(width=0.5)
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")

seasRes = preds %>% gather(season, exp_perc, c(X.Exp_perc_summer,X.Exp_perc_autumn,X.Exp_perc_winter,X.Exp_perc_spring))
seasRes$season = factor(seasRes$season, levels=c('X.Exp_perc_summer','X.Exp_perc_autumn','X.Exp_perc_winter','X.Exp_perc_spring'))
labels <- c(X.Exp_perc_summer='Summer % change',X.Exp_perc_autumn='Autumn % change',X.Exp_perc_winter='Winter % change',
            X.Exp_perc_spring="Spring % change")
seasRes$RainWintOrSum = d$RainWintOrSum
jpeg('Meta-analysis/Figures_results/PrecipChange_vs_EffectSize_predicted.jpg', width=12, height=8, res=300, units='cm')
ggplot(seasRes, aes(exp_perc, pred, colour=RainWintOrSum)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  ylim(-2,2) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  scale_color_manual(name='Rain in winter\n or summer', 
                     values=c("#d7191c", "#2b83ba", "#fdae61"),
                     labels=c('Summer','Winter','All year')) +
  labs(x='% change in precipitation', y="Predicted effect size") +
  theme_bw() +
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="black")
dev.off()



par(mar=c(4,4,1,2))
forest(sav$pred, sei=sav$se, slab=slabs, transf=transf.ztor, xlab="Correlation", xlim=c(-.4,.7))
text(-.4, 8, "Structure/Type",       pos=4, font=2)
text(.7,  8, "Correlation [95% CI]", pos=2, font=2)

#------------------------------------------------------------------------------------
#Seasonal models with RainWintorSum:
#Summer
mSumRWS = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_summer * RainWintOrSum - 1,
                 random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d) #[d$Exp_perc_summer!=100,]
summary(mSumRWS)
preds = predict(mSumRWS, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
sumRWS = preds %>% gather(RWS, RWSdummy, c(X.RainWintOrSums,X.RainWintOrSumw,X.RainWintOrSumy))
sumRWS = sumRWS[sumRWS$RWSdummy==1,]
jpeg('Meta-analysis/Figures_results/ChangeSummer_vs_BiomassChange_RWS_predicted.jpg', width=16, height=12, res=300, units='cm')
ggplot(sumRWS, aes(X.Exp_perc_summer, pred, colour=RWS)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  labs(x='Summer % change', y="Predicted effect size") +
  scale_colour_discrete(name="Rain season",
                        breaks=c('X.RainWintOrSums','X.RainWintOrSumw','X.RainWintOrSumy'),
                        labels=c("Summer", "Winter", "All year")) +
  theme_bw()
dev.off()

qqnorm(residuals(mSumRWS,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mSumRWS,type="pearson"),col="red")

#Autumn
mAutRWS = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_autumn * RainWintOrSum - 1,
                 random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d) #[d$Exp_perc_summer!=100,]
summary(mAutRWS)
preds = predict(mAutRWS, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
autRWS = preds %>% gather(RWS, RWSdummy, c(X.RainWintOrSums,X.RainWintOrSumw,X.RainWintOrSumy))
autRWS = autRWS[autRWS$RWSdummy==1,]
# jpeg('Meta-analysis/Figures_results/ChangeAutmer_vs_BiomassChange_RWS_predicted.jpg', width=16, height=12, res=300, units='cm')
ggplot(autRWS, aes(X.Exp_perc_autumn, pred, colour=RWS)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  labs(x='Autumn % change', y="Predicted effect size") +
  scale_colour_discrete(name="Rain season",
                        breaks=c('X.RainWintOrSums','X.RainWintOrSumw','X.RainWintOrSumy'),
                        labels=c("Summer", "Winter", "All year")) +
  theme_bw()
# dev.off()

#Winter
mWinRWS = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_winter * RainWintOrSum - 1,
                 random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d) #[d$Exp_perc_summer!=100,]
summary(mWinRWS)
preds = predict(mWinRWS, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
winRWS = preds %>% gather(RWS, RWSdummy, c(X.RainWintOrSums,X.RainWintOrSumw,X.RainWintOrSumy))
winRWS = winRWS[winRWS$RWSdummy==1,]
jpeg('Meta-analysis/Figures_results/ChangeWinter_vs_BiomassChange_RWS_predicted.jpg', width=16, height=12, res=300, units='cm')
ggplot(winRWS, aes(X.Exp_perc_winter, pred, colour=RWS)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  labs(x='Winter % change', y="Predicted effect size") +
  scale_colour_discrete(name="Rain season",
                        breaks=c('X.RainWintOrSums','X.RainWintOrSumw','X.RainWintOrSumy'),
                        labels=c("Summer", "Winter", "All year")) +
  theme_bw()
dev.off()



#Spring
mSprRWS = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_spring * RainWintOrSum - 1,
                 random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d) #[d$Exp_perc_summer!=100,]
summary(mSprRWS)
preds = predict(mSprRWS, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
sprRWS = preds %>% gather(RWS, RWSdummy, c(X.RainWintOrSums,X.RainWintOrSumw,X.RainWintOrSumy))
sprRWS = sprRWS[sprRWS$RWSdummy==1,]
jpeg('Meta-analysis/Figures_results/ChangeSpring_vs_BiomassChange_RWS_predicted.jpg', width=16, height=12, res=300, units='cm')
ggplot(sprRWS, aes(X.Exp_perc_spring, pred, colour=RWS)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  labs(x='Spring % change', y="Predicted effect size") +
  scale_colour_discrete(name="Rain season",
                        breaks=c('X.RainWintOrSums','X.RainWintOrSumw','X.RainWintOrSumy'),
                        labels=c("Summer", "Winter", "All year")) +
  theme_bw()
dev.off()


#------------------------------------------------------------------------------------
#Seasonal models with photosynthetic type:
#Summer
mSumPtype = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_summer * C3_C4 - 1,
                   random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d) #[d$Exp_perc_summer!=100,]
summary(mSumPtype)
preds = predict(mSumPtype, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
sumPtype = preds %>% gather(ptype, ptypedummy, c(X.C3_C4Both,X.C3_C4C3,X.C3_C4C4))
sumPtype = sumPtype[sumPtype$ptypedummy==1,]
labels <- c(X.Exp_perc_summer='Summer % change',X.Exp_perc_autumn='Autumn % change',X.Exp_perc_winter='Winter % change',
            X.Exp_perc_spring="Spring % change")
# jpeg('Meta-analysis/Figures_results/PrecipChange_vs_EffectSize_predicted.jpg', width=20, height=16, res=300, units='cm')
ggplot(sumPtype, aes(X.Exp_perc_summer, pred, colour=ptype)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  labs(x='Summer % change', y="Predicted effect size") +
  theme_bw()
# dev.off()


#Autumn
mAut = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_autumn,
              random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d[d$Exp_perc_autumn!=100,])
summary(mAut)

#Winter
mWinPtype = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_winter * C3_C4 - 1,
                   random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(mWinPtype)
preds = predict(mWinPtype, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
winPtype = preds %>% gather(ptype, ptypedummy, c(X.C3_C4Both,X.C3_C4C3,X.C3_C4C4))
winPtype = winPtype[winPtype$ptypedummy==1,]
# jpeg('Meta-analysis/Figures_results/PrecipChange_vs_EffectSize_predicted.jpg', width=20, height=16, res=300, units='cm')
ggplot(winPtype, aes(X.Exp_perc_winter, pred, colour=ptype)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  labs(x='Winter % change', y="Predicted effect size") +
  theme_bw()
# dev.off()

#Spring
mSprPtype = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_spring * C3_C4 - 1,
                   random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(mSprPtype)

#All
mSeasPtype = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_spring + Exp_perc_summer + Exp_perc_winter + Exp_perc_autumn + C3_C4 - 1,
                    random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(mSeasPtype)
preds = predict(mSeasPtype, addx=TRUE, digits=2)
preds = do.call(cbind.data.frame, preds)
ggplot(preds, aes(x=X.Exp_perc_winter,y=pred,colour=as.factor(X.C3_C4C3))) + 
  geom_point() + #position=position_jitter(width=0.5)
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
ggplot(preds, aes(x=X.Exp_perc_summer,y=pred,colour=as.factor(X.C3_C4C3))) + 
  geom_point() + #position=position_jitter(width=0.5)
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")


#---------------------------------------------------------------------------------------------------------------
#Plotting

#Plot response of biomass to annual precip:
jpeg('Meta-analysis/Figures_results/Precip_vs_biomass_control.jpg', width=16, height=12, res=300, units='cm')
ggplot(d, aes(Control_rain, Biomass_control, colour=continent)) + #factor(RainWintOrSum, labels=c('Summer','Winter','All year'))
  geom_point() +
  geom_text(aes(label=TreatmentID),hjust=0, vjust=0, size=2) +
  labs(x='Control rainfall', y=expression(Control~biomass~(g.m^-2))) +
  theme_bw() #+  stat_smooth(method="lm",formula=y~x,fullrange=T,se=F,size=1,colour="red")
dev.off()

jpeg('Meta-analysis/Figures_results/Precip_vs_biomass_treatment.jpg', width=16, height=12, res=300, units='cm')
ggplot(d, aes(Treatment_rain, Biomass_treatment, colour=continent)) + #factor(RainWintOrSum, labels=c('Summer','Winter','All year'))
  geom_point() +
  geom_text(aes(label=TreatmentID),hjust=0, vjust=0, size=2) +
  labs(x='Treatment rainfall', y=expression(Treatment~biomass~(g.m^-2))) +
  theme_bw() #+ stat_smooth(method="lm",formula=y~x,fullrange=T,se=F,size=1,colour="red")
dev.off()

#Plot biomass difference vs precip difference:
plotdat = d %>% gather(season, exp_perc, c(Exp_perc_summer,Exp_perc_autumn,Exp_perc_winter,Exp_perc_spring,Exp_perc_year))
plotdat$season = factor(plotdat$season, levels=c('Exp_perc_summer','Exp_perc_autumn','Exp_perc_winter','Exp_perc_spring','Exp_perc_year'))
labels <- c(Exp_perc_summer='Summer % change',Exp_perc_autumn='Autumn % change',Exp_perc_winter='Winter % change',Exp_perc_spring="Spring % change",
            Exp_perc_year='Year % change')
jpeg('Meta-analysis/Figures_results/PrecipChange_vs_EffectSize.jpg', width=20, height=16, res=300, units='cm')
ggplot(plotdat, aes(exp_perc, HedgeD)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
dev.off()
jpeg('Meta-analysis/Figures_results/PrecipChange_vs_EffectSize_RainWintorSum.jpg', width=20, height=16, res=300, units='cm')
ggplot(plotdat, aes(exp_perc, HedgeD, colour=RainWintOrSum)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
dev.off()
jpeg('Meta-analysis/Figures_results/PrecipChange_vs_EffectSize_TreatmentID.jpg', width=20, height=16, res=300, units='cm')
ggplot(plotdat, aes(exp_perc, HedgeD, colour=TreatmentID)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  geom_text(aes(label=TreatmentID, colour=TreatmentID), hjust=0, vjust=0, size=2) +
  scale_colour_discrete(guide = 'none') +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
dev.off()
jpeg('Meta-analysis/Figures_results/PrecipSpringChange_vs_EffectSize_TreatmentID.jpg', width=12, height=8, res=300, units='cm')
ggplot(d, aes(Exp_perc_spring, HedgeD, colour=TreatmentID)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  geom_text(aes(label=TreatmentID, colour=TreatmentID), hjust=0, vjust=0, size=2) +
  scale_colour_discrete(guide = 'none') +
  labs(x='Spring % change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
dev.off()
#Plot frequency distributions of variables:
#Response (change in biomass - Hedge's D)
ggplot(data=d, aes(d$HedgeD)) + 
  geom_histogram() +
  labs(x="Hedge's D", y="Count")
#Summer % change
ggplot(data=d, aes(d$Exp_perc_summer)) + 
  geom_histogram() +
  labs(x="Summer % change", y="Count")
#Autumn % change
ggplot(data=d, aes(d$Exp_perc_autumn)) + 
  geom_histogram() +
  labs(x="Autumn % change", y="Count")
#Winter % change
ggplot(data=d, aes(d$Exp_perc_winter)) + 
  geom_histogram() +
  labs(x="Winter % change", y="Count")
#Spring % change
ggplot(data=d, aes(d$Exp_perc_spring)) + 
  geom_histogram() +
  labs(x="Spring % change", y="Count")


ggplot(plotdat, aes(exp_perc, HedgeD, colour=Rain_data_full_year)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
boxplot(d$HedgeD~d$Rain_data_full_year)

ggplot(plotdat, aes(exp_perc, HedgeD, colour=Change_seasonal_shift_or_seasonal_change_or_annual)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
boxplot(d$HedgeD~d$Change_seasonal_shift_or_seasonal_change_or_annual)

ggplot(plotdat, aes(exp_perc, HedgeD, colour=Study_length)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
ggplot(d, aes(Study_length, HedgeD)) +
  geom_point() +
  geom_line(aes(group=TreatmentID))  +
  geom_text(aes(label=TreatmentID, colour=TreatmentID),hjust=0, vjust=0, size=2, show.legend = FALSE) +
  labs(x='Study_length', y="Hedge's D (biomass change)") +
  theme_bw()

ggplot(plotdat, aes(exp_perc, HedgeD, colour=as.factor(Grazing_code))) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
boxplot(d$HedgeD~d$Grazing_code)

ggplot(plotdat, aes(exp_perc, HedgeD, colour=as.factor(C3_C4))) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
boxplot(d$HedgeD~d$C3_C4)

ggplot(plotdat, aes(exp_perc, HedgeD, colour=as.factor(native_alien))) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels)) +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()
ggplot(plotdat, aes(native_alien, HedgeD)) +
  geom_boxplot() +
  # scale_colour_gradientn(colours=rev(c('#2c7bb6','#abd9e9','#ffffbf','#fee090','#fdae61','#d73027'))) +
  facet_wrap(~season, ncol=2, labeller=labeller(season = labels))# +
# labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
# theme_bw()



ggplot(plotdat[plotdat$season=='Exp_perc_summer' & plotdat$exp_perc<=95 | plotdat$season=='Exp_perc_summer' & plotdat$exp_perc>=105,], 
       aes(exp_perc, HedgeD, colour=RainWintOrSum)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()

ggplot(plotdat[plotdat$season=='Exp_perc_spring' & plotdat$exp_perc<=95 | plotdat$season=='Exp_perc_spring' & plotdat$exp_perc>=105,], 
       aes(exp_perc, HedgeD, colour=RainWintOrSum)) +
  geom_vline(xintercept=100) + geom_hline(yintercept=0) +
  geom_point() +
  labs(x='% change in precipitation', y="Hedge's D (Change in biomass)") +
  theme_bw()


#---------------------------------------------------------------------------------------------------------------
#Random

m1 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_winter + native_alien - 1,
            random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(m1)
anova(m1, btt=2:3)
anova(m1, L=c(0,0,0,0,1,-1))

m1 = rma.mv(yi=HedgeD, V=d_var, mods= ~ native_alien + Exp_perc_summer* C3_C4 + Exp_perc_winter* C3_C4 + Exp_perc_spring* C3_C4 + Exp_perc_autumn * C3_C4-1, 
            random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(m1)
preds = predict(m1, addx=TRUE)
preds = do.call(cbind.data.frame, preds)
library(ggplot2)
ggplot(preds, aes(x=X.Exp_perc_winter,y=pred)) + 
  geom_point(position=position_jitter(width=0.5)) +
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
ggplot(preds, aes(x=X.Exp_perc_winter,y=pred,colour=as.factor(X.native_alienboth))) + 
  geom_point() +
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
ggplot(preds, aes(x=X.Exp_perc_spring,y=pred)) + 
  geom_point(position=position_jitter(width=0.5)) +
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
ggplot(preds, aes(x=X.Exp_perc_autumn,y=pred)) + 
  geom_point(position=position_jitter(width=0.5)) +
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")
ggplot(preds, aes(x=X.warmtocold,y=pred)) + 
  geom_point(position=position_jitter(width=0.5)) +
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")


m1ml = rma.mv(yi=HedgeD, V=d_var, mods= ~ Treatment_rain + Exp_perc_summer + Exp_perc_winter + Exp_perc_spring + Exp_perc_autumn + C3_C4 +
                continent, 
              random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="ML",digits=4, data=d)
m2ml = rma.mv(yi=HedgeD, V=d_var, mods= ~ Treatment_rain + Exp_perc_summer + Exp_perc_winter + Exp_perc_spring + Exp_perc_autumn + continent, 
              random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="ML",digits=4, data=d)
anova(m1ml,m2ml)
anova(m1, btt=1)
anova(m1, L=c(0,0,0,0,1,-1))

m1 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_winter, 
            random= ~1|Ref/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(m1)

m1 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_summer * C3_C4, 
            random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(m1)
anova(m1, btt=3:4)
anova(m1, L=c(0,0,1,-1,0))

m1 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_winter:RainWintOrSum - 1, 
            random= ~1|Study/TreatmentID/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(m1)
anova(m1, L=c(0,1,-1))
anova(m1, L=c(1,-1,0))
anova(m1, L=c(1,0,-1))

m1 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Exp_perc_spring * Exp_perc_summer, 
            random= ~1|Ref/Obs_year,struct="CS",method="REML",digits=4, data=d)
summary(m1)
anova(m1, L=c(0,1,-1))
anova(m1, L=c(1,-1,0))
anova(m1, L=c(1,0,-1))

plot(coef(m1)[1:3], type="o", pch=19, bty="l")
axis(side=1, at=1:3, labels=c("none","some","high"))
lines(coef(res.i2)[4:6], type="o", pch=15, lty="dotted")
legend("topright", legend=c("Test Administrator Blind", "Test Administrator Aware"), lty=c("solid", "dotted"), pch=c(19,15), inset=0.01)
title("Estimated Average Effects based on the Interaction Model")


library(multcomp)
#For pairwise comparison of RainWintOrSum:
summary(glht(m1, linfct=rbind(c(0,0,1,0,0,0), c(0,0,0,1,0,0), c(0,0,1,-1,0,0))), test=adjusted("none"))


m1 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Control_rain + 
              RainWintOrSum +
              Exp_perc_summer + 
              Exp_perc_winter + 
              Exp_perc_spring + 
              Exp_perc_autumn, 
            random= ~1|Ref/Obs_year,struct="CS",method="REML",digits=4, data=d)

m2 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Control_rain + 
              Exp_perc_summer + 
              Exp_perc_winter + 
              Exp_perc_spring + 
              Exp_perc_autumn, 
            random= ~1|Ref/Obs_year,struct="CS",method="REML",digits=4, data=d)

m2 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Control_rain + 
              RainWintOrSum +
              Exp_perc_summer + #Exp_perc_summer:RainWintOrSum +
              Exp_perc_winter + Exp_perc_winter:RainWintOrSum +
              Exp_perc_spring + Exp_perc_spring:RainWintOrSum + 
              Exp_perc_autumn + Exp_perc_autumn:RainWintOrSum-1, 
            random= ~1|Ref/Obs_year,struct="CS",method="ML",digits=4, data=d)


summary(m1)

library(multcomp)
#For pairwise comparison of RainWintOrSum:
summary(glht(m1, linfct=rbind(c(0,0,1,0,0,0,0,0), c(0,0,0,1,0,0,0,0), c(0,0,1,-1,0,0,0,0))), test=adjusted("none"))
#For pairwise comparison of RainWintOrSum:
summary(glht(m1, linfct=rbind('' c(0,0,0,0,1,-1,0,0), c(0,0,0,0,1,0,-1,0), c(0,0,0,0,1,0,0,-1), c(0,0,0,0,0,1,-1,0),
                              c(0,0,0,0,0,1,0,-1), c(0,0,0,0,0,0,1,-1)
)), test=adjusted("none"))

m1 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Control_rain + 
              RainWintOrSum +
              Exp_perc_summer + 
              Exp_perc_winter + 
              Exp_perc_spring + 
              Exp_perc_autumn, 
            random= ~1|Ref/Obs_year,struct="CS",method="ML",digits=4, data=d)
m2 = rma.mv(yi=HedgeD, V=d_var, mods= ~ Control_rain + 
              RainWintOrSum +
              Exp_perc_summer + 
              Exp_perc_winter + 
              Exp_perc_spring + 
              Exp_perc_autumn, 
            random= ~1|Ref/Obs_year,struct="CS",method="ML",digits=4, data=d)
anova(m1,m2)
anova(m1, btt=6)

qqnorm(residuals(m1,type="pearson"),main="QQ plot: residuals")
qqline(residuals(m1,type="pearson"),col="red")

### Use model to make predictions

preds = predict(m1, levels=0, addx=TRUE)
preds = do.call(cbind.data.frame, preds)
# names(preds)[9] = 'Exp_perc_summer'

library(ggplot2)
ggplot(preds, aes(x=X.Exp_perc_winter,y=pred)) + 
  stat_smooth(method="lm",formula=y~x,fullrange=T,se=FALSE,size=1,colour="red")+
  geom_point(data=d,aes(x=Exp_perc_winter,y=HedgeD),position=position_jitter(width=0.2),pch=19,colour="blue",size=3)+
  # labs(x="Rain",y ="Effect size") +
  theme(axis.title.x=element_text(colour="black",face="bold",size=15,vjust=-0.5),
        axis.title.y=element_text(colour="black",face="bold",size=15,vjust=1),
        axis.text.y=element_text(colour="black",face="bold",size=12),
        axis.text.x=element_text(colour="black",face="bold",size=12),
        panel.background =element_rect(fill="transparent",colour="black"),
        panel.border=element_rect(fill=NA,colour="black"))

preds$RainWintOrSum = d$RainWintOrSum
ggplot(preds, aes(x=X.Exp_perc_winter,y=pred, colour=RainWintOrSum)) + 
  geom_point(position=position_jitter(width=0.5)) +
  geom_line(aes(y=predict(m1, newmods=cbind(0,0,0,0,0,preds$X.Exp_perc_winter,0,0,1,0,0,0,0))$pred)) +
  geom_line(aes(y=predict(m1, newmods=cbind(0,0,0,0,0,preds$X.Exp_perc_winter,0,0,0,1,0,0,0))$pred))



newdat <- expand.grid(RainWintOrSum=unique(d$RainWintOrSum),
                      X.Exp_perc_winter=c(min(d$Exp_perc_winter),
                                          max(d$Exp_perc_winter)))
newdat <- model.matrix(~ RainWintOrSum + X.Exp_perc_winter, data=newdat)

d$X.Exp_perc_winter = d$Exp_perc_winter
ggplot(preds, aes(x=X.Exp_perc_winter, y=HedgeD, colour=RainWintOrSum)) +
  geom_point(size=3) +
  geom_line(data = preds,aes(x = X.Exp_perc_winter,y = preds$pred, size = 1))

ggplot(preds) +
  geom_line(aes(x = X.Exp_perc_winter,y = pred, size = 1)) +
  geom_line(aes(y=predict(m1)$pred, group=RainWintOrSum)) +
  geom_line(data=newdat, aes(y=predict(m1, level=0, newdata=newdat)$pred))


predict(m1, newmods=cbind(0,0,0,0,1:200,0,0,0,1,0,0,0,0))$pred

newdat = data.frame(X.Control_rain=mean(d$Control_rain),
                    X.RainWintOrSumw=1,
                    X.RainWintOrSumy=0,
                    X.Exp_perc_summer=mean(d$Exp_perc_summer),
                    X.Exp_perc_winter=d$X.Exp_perc_winter,
                    X.Exp_perc_spring=mean(d$Exp_perc_spring),
                    X.Exp_perc_autumn=mean(d$Exp_perc_autumn),
                    X.RainWintOrSumw.Exp_perc_winter=1,
                    X.RainWintOrSumy.Exp_perc_winter=0,
                    X.RainWintOrSumw.Exp_perc_spring=0,
                    X.RainWintOrSumy.Exp_perc_spring=0,
                    X.RainWintOrSumw.Exp_perc_autumn=0)
predict(m1, newdata=newdat)
ggplot(preds, aes(x=X.Exp_perc_winter,y=pred, colour=RainWintOrSum)) + 
  geom_point(position=position_jitter(width=0.5)) +
  geom_line(aes(y=predict(m1, newdata=newdat)$pred))


addmargins(table(dat$weeks, dat$tester))



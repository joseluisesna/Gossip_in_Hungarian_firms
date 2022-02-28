########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## Multilevel analysis (4)
## R script written by Jose Luis Estevez (Masaryk University)
## Date: Feb 27h 2022
########################################################################################################################

# R PACKAGES REQUIRED
library(lme4)
# DATA LOADING
rm(list=ls())
load('modellingdata.RData')

########################################################################################################################

# 1) PREPARATION OF THE DATA

# Creation of single large dataset for triadic relations (all units together)
triad_data <- do.call('rbind',triad_data)

# Creation of variables
triad_data$neg_gossip <- 1*(triad_data$gossip == -1) # negative gossip
triad_data$pos_gossip <- 1*(triad_data$gossip == 1) # positive gossip

triad_data$SR_pos <- 1*(triad_data$SR == 1) # positive tie (SR)
triad_data$ST_pos <- 1*(triad_data$ST == 1) # positive tie (ST)
triad_data$RT_pos <- 1*(triad_data$RT == 1) # positive tie (RT)
triad_data$SR_neg <- 1*(triad_data$SR == -1) # negative tie (SR)
triad_data$ST_neg <- 1*(triad_data$ST == -1) # negative tie (ST)
triad_data$RT_neg <- 1*(triad_data$RT == -1) # negative tie (RT)

# re-level (nonbroker as the reference level)
triad_data$sender_role <- factor(triad_data$sender_role,levels=c('nonbroker','broker'))
triad_data$receiver_role <- factor(triad_data$receiver_role,levels=c('nonbroker','broker'))
triad_data$target_role <- factor(triad_data$target_role,levels=c('nonbroker','broker'))

# Addition of gender and hierarchical position
triad_data <- merge(x=triad_data,y=attributes[,c('responder','woman')],by.x='sender',by.y='responder',all.x=TRUE)
triad_data <- merge(x=triad_data,y=attributes[,c('responder','woman')],by.x='receiver',by.y='responder',all.x=TRUE)
triad_data <- merge(x=triad_data,y=attributes[,c('responder','woman')],by.x='target',by.y='responder',all.x=TRUE)
triad_data <- merge(x=triad_data,y=attributes[,c('responder','hr_leader')],by.x='sender',by.y='responder',all.x=TRUE)
triad_data <- merge(x=triad_data,y=attributes[,c('responder','hr_leader')],by.x='receiver',by.y='responder',all.x=TRUE)
triad_data <- merge(x=triad_data,y=attributes[,c('responder','hr_leader')],by.x='target',by.y='responder',all.x=TRUE)
names(triad_data) <- c(names(triad_data)[1:22],
                       'sender_woman','receiver_woman','target_woman','sender_boss','receiver_boss','target_boss')

# Create an effect for the isolates (to complement the group)
for(i in seq_along(networks_mtx)){
  networks_mtx[[i]]$isolates <- rowSums(networks_mtx[[i]]$positive,na.rm=TRUE)+
    colSums(networks_mtx[[i]]$positive,na.rm=TRUE)
  networks_mtx[[i]]$isolates <- networks_mtx[[i]]$isolates[networks_mtx[[i]]$isolates == 0]
}

isolates <- c(networks_mtx[[1]]$isolates,networks_mtx[[2]]$isolates,networks_mtx[[3]]$isolates,
              networks_mtx[[4]]$isolates,networks_mtx[[5]]$isolates,networks_mtx[[6]]$isolates)
triad_data$sender_iso <- 1*triad_data$sender %in% names(isolates)
triad_data$receiver_iso <- 1*triad_data$receiver %in% names(isolates)
triad_data$target_iso <- 1*triad_data$target %in% names(isolates)

########################################################################################################################

# 2) MODELLING 

# 2.0) Null model (only random effects)
results_pos0 <- glmer(data=triad_data,pos_gossip ~ 1 +
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_pos0)
results_neg0 <- glmer(data=triad_data,neg_gossip ~ 1 +
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_neg0)

# 2.1) Model 1 (dyadic effects and individuals-level controls)
results_pos1 <- glmer(data=triad_data,pos_gossip ~
                        SR_pos + ST_pos + RT_pos +
                        SR_neg + ST_neg + RT_neg +
                        receiver_boss + sender_boss + target_boss +
                        receiver_woman + sender_woman + target_woman +
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_pos1)
results_neg1 <- glmer(data=triad_data,neg_gossip ~
                        SR_pos + ST_pos + RT_pos +
                        SR_neg + ST_neg + RT_neg +
                        receiver_boss + sender_boss + target_boss +
                        receiver_woman + sender_woman + target_woman +
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_neg1)

# Model comparison
anova(results_pos0,results_pos1);anova(results_neg0,results_neg1)

# 2.2) Model 2 (group effects)
results_pos2 <- glmer(data=triad_data,pos_gossip ~
                        SR_pos + ST_pos + RT_pos +
                        SR_neg + ST_neg + RT_neg +
                        receiver_boss + sender_boss + target_boss +
                        receiver_woman + sender_woman + target_woman +
                        samegroup_SR + samegroup_ST + samegroup_RT + 
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_pos2)
results_neg2 <- glmer(data=triad_data,neg_gossip ~
                        SR_pos + ST_pos + RT_pos +
                        SR_neg + ST_neg + RT_neg +
                        receiver_boss + sender_boss + target_boss +
                        receiver_woman + sender_woman + target_woman +
                        samegroup_SR + samegroup_ST + samegroup_RT +
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_neg2)

# Model comparison
anova(results_pos1,results_pos2);anova(results_neg1,results_neg2)

# 2.3) Model 3 (brokering effects)
results_pos3 <- glmer(data=triad_data,pos_gossip ~
                        SR_pos + ST_pos + RT_pos +
                        SR_neg + ST_neg + RT_neg +
                        receiver_boss + sender_boss + target_boss +
                        receiver_woman + sender_woman + target_woman +   
                        samegroup_SR + samegroup_ST + samegroup_RT +
                        receiver_role + sender_role + target_role +
                        receiver_iso + sender_iso + target_iso +
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_pos3)
results_neg3 <- glmer(data=triad_data,neg_gossip ~
                        SR_pos + ST_pos + RT_pos +
                        SR_neg + ST_neg + RT_neg +
                        receiver_boss + sender_boss + target_boss +
                        receiver_woman + sender_woman + target_woman +   
                        samegroup_SR + samegroup_ST + samegroup_RT +
                        receiver_role + sender_role + target_role +
                        receiver_iso + sender_iso + target_iso + 
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_neg3)

# Model comparison
anova(results_pos2,results_pos3);anova(results_neg2,results_neg3)

########################################################################################################################

# Save image
save.image('ML_results.RData')

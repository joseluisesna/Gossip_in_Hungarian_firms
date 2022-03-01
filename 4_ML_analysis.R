########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## Multilevel analysis (4)
## R script written by Jose Luis Estevez (Masaryk University)
## Date: Mar 1st 2022
########################################################################################################################

# R PACKAGES REQUIRED
library(lme4);library(effectsize);library(ggplot2)
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
triad_data$sender_role <- factor(triad_data$sender_role,levels=c('non-broker','broker'))
triad_data$receiver_role <- factor(triad_data$receiver_role,levels=c('non-broker','broker'))
triad_data$target_role <- factor(triad_data$target_role,levels=c('non-broker','broker'))

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

# Gossip by broker status
gossip_broker <- as.data.frame(matrix(NA,nrow=6,ncol=8))
rownames(gossip_broker) <- c('Unit A','Unit B','Unit C','Unit D','Unit E','Unit F')
colnames(gossip_broker) <- c('pos_sen','pos_rec','pos_tar','pos','neg_sen','neg_rec','neg_tar','neg')

gossip_broker[,1] <- table(triad_data[triad_data$gossip == 1,]$unit,triad_data[triad_data$gossip == 1,]$sender_role)[,'broker']
gossip_broker[,2] <- table(triad_data[triad_data$gossip == 1,]$unit,triad_data[triad_data$gossip == 1,]$receiver_role)[,'broker']
gossip_broker[,3] <- table(triad_data[triad_data$gossip == 1,]$unit,triad_data[triad_data$gossip == 1,]$target_role)[,'broker']
gossip_broker[,4] <- table(triad_data[triad_data$gossip == 1,]$unit)

gossip_broker[,5] <- table(triad_data[triad_data$gossip == -1,]$unit,triad_data[triad_data$gossip == -1,]$sender_role)[,'broker']
gossip_broker[,6] <- table(triad_data[triad_data$gossip == -1,]$unit,triad_data[triad_data$gossip == -1,]$receiver_role)[,'broker']
gossip_broker[,7] <- table(triad_data[triad_data$gossip == -1,]$unit,triad_data[triad_data$gossip == -1,]$target_role)[,'broker']
gossip_broker[,8] <- table(triad_data[triad_data$gossip == -1,]$unit)

gossip_broker
colSums(gossip_broker)

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
                        samegroup_SR + samegroup_ST + samegroup_RT + samegroup_SR:samegroup_ST +
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_pos2)
results_neg2 <- glmer(data=triad_data,neg_gossip ~
                        SR_pos + ST_pos + RT_pos +
                        SR_neg + ST_neg + RT_neg +
                        receiver_boss + sender_boss + target_boss +
                        receiver_woman + sender_woman + target_woman +
                        samegroup_SR + samegroup_ST + samegroup_RT + samegroup_SR:samegroup_ST +
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
                        samegroup_SR + samegroup_ST + samegroup_RT + samegroup_SR:samegroup_ST +
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
                        samegroup_SR + samegroup_ST + samegroup_RT + samegroup_SR:samegroup_ST +
                        receiver_role + sender_role + target_role +
                        receiver_iso + sender_iso + target_iso + 
                        (1|unit) + (1|unit:receiver) + (1|unit:sender) + (1|unit:target),
                      family=binomial(link='logit'))
summary(results_neg3)

# Model comparison
anova(results_pos2,results_pos3);anova(results_neg2,results_neg3)

########################################################################################################################

# 3) STANDARDISED COEFFICIENTS
# Positive gossip
effsize_pos1 <- effectsize(results_pos1)
effsize_pos2 <- effectsize(results_pos2)
effsize_pos3 <- effectsize(results_pos3)
# Negative gossip
effsize_neg1 <- effectsize(results_neg1)
effsize_neg2 <- effectsize(results_neg2)
effsize_neg3 <- effectsize(results_neg3)

effsize_pos3
effsize_neg3

# Visualisation
effsize_pos3$gossiptype <- 'Positive gossip'
effsize_neg3$gossiptype <- 'Negative gossip'
effsize_plot <- rbind(effsize_pos3,effsize_neg3)

effsize_plot$Parameter <- factor(effsize_plot$Parameter,
                                 levels=c('target_iso','receiver_iso','sender_iso',
                                          'target_rolebroker','receiver_rolebroker','sender_rolebroker',
                                          'samegroup_SR:samegroup_ST','samegroup_RT','samegroup_ST','samegroup_SR',
                                          'target_boss','receiver_boss','sender_boss',
                                          'target_woman','receiver_woman','sender_woman',
                                          'RT_neg','ST_neg','SR_neg',
                                          'RT_pos','ST_pos','SR_pos','(Intercept)'),
                                 labels=c('Isolate (target)','Isolate (receiver)','Isolate (sender)',
                                          'Broker (target)','Broker (receiver)','Broker (sender)',
                                          'Same group (sender-receiver-target)','Same group (receiver-target)',
                                          'Same group (sender-target)','Same group (sender-receiver)',
                                          'Manager (target)','Manager (receiver)','Manager (sender)',
                                          'Woman (target)','Woman (receiver)','Woman (sender)',
                                          'Negative tie (receiver-target)','Negative tie (sender-target)','Negative tie (sender-receiver)',
                                          'Positive tie (receiver-target)','Positive tie (sender-target)','Positive tie (sender-receiver)',
                                          'Intercept'))
effsize_plot$gossiptype <- factor(effsize_plot$gossiptype,levels=c('Positive gossip','Negative gossip'))

grid.background <- theme_bw()+
  theme(plot.background=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+
  theme(axis.line=element_line(color='black'))+
  theme(strip.text.x=element_text(colour='white',face='bold'))+
  theme(strip.background=element_rect(fill='black'))

jpeg(filename='Results.jpeg',width=10,height=8,units='in',res=500)
ggplot(data=effsize_plot[effsize_plot$Parameter != 'Intercept',], aes(x=Parameter, y=Std_Coefficient, ymin=CI_low, ymax=CI_high)) +
  geom_hline(yintercept= 0, lty=2,colour='red') +
  geom_pointrange(position=position_dodge(width = 0.5)) +
  facet_wrap(~gossiptype) +
  coord_flip() +
  xlab("") + ylab("Standardised coefficient (95% CI)") +
  theme_bw() +
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12)) +
  grid.background
dev.off()

########################################################################################################################

# Save image
save.image('ML_results.RData')

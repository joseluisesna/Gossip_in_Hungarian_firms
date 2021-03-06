########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## Descriptive analysis (3)
## R script written by Jose Luis Estevez (Linkoping University)
## Date: April 9th 2021
########################################################################################################################

# R PACKAGES REQUIRED
library(sna);library(igraph)

# DATA LOADING
rm(list=ls())
load('tidieddata2.RData')

########################################################################################################################

# 0) BASIC DESCRIPTIVE STATS OF THE SAMPLE

# Number of nodes per unit
summary(attributes$group)

# Women per unit
for(i in organisation_ID){
  print(sum(attributes[attributes$group == i,]$woman == 1,na.rm=TRUE))
}

# Managers per unit
for(i in organisation_ID){
  print(sum(attributes[attributes$group == i,]$hr_leader == 1,na.rm=TRUE))
}

# Managers per unit
for(i in organisation_ID){
  print(sum(attributes[attributes$group == i & attributes$woman == 1,]$hr_leader == 1,na.rm=TRUE))
}

# Age
for(i in organisation_ID){
  print(mean(2018 - attributes[attributes$group == i,]$birth_year,na.rm=TRUE))
  print(min(2018 - attributes[attributes$group == i,]$birth_year,na.rm=TRUE))
  print(max(2018 - attributes[attributes$group == i,]$birth_year,na.rm=TRUE))
}

# Tenure
attributes$tenure_year <- as.numeric(substr(attributes$hr_work_start,1,4))
attributes$tenure_month <- as.numeric(substr(attributes$hr_work_start,6,7))
attributes$tenure_day <- as.numeric(substr(attributes$hr_work_start,9,10))
attributes$tenure <- round(2018 - attributes$tenure_year - attributes$tenure_month/12 - attributes$tenure_day/30,1)

for(i in organisation_ID){
  print(mean(attributes[attributes$group == i,]$tenure,na.rm=TRUE))
  print(min(attributes[attributes$group == i,]$tenure,na.rm=TRUE))
  print(max(attributes[attributes$group == i,]$tenure,na.rm=TRUE))
}

########################################################################################################################

# 1) DESCRIPTIVE STATS (POSITIVE & NEGATIVE TIES)

# Object to store summary statistics
pos_mtx_sum <- neg_mtx_sum <- matrix(NA,9,length(organisation_ID),
                                     dimnames=list(c('nodes','missing-tie','density','recip','closure','isolates',
                                                     'avg. degree','SD out','SD in'),organisation_ID))

for(i in organisation_ID){
  # positive ties
  pos_mtx_sum[1,i] <- nrow(networks_mtx[[i]]$positive)
  pos_mtx_sum[2,i] <- 1-sum(!is.na(networks_mtx[[i]]$positive))/(nrow(networks_mtx[[i]]$positive)*(nrow(networks_mtx[[i]]$positive)-1))
  pos_mtx_sum[3,i] <- edge_density(graph_from_adjacency_matrix(networks_mtx[[i]]$positive,mode='directed',diag=FALSE))
  pos_mtx_sum[4,i] <- reciprocity(graph_from_adjacency_matrix(networks_mtx[[i]]$positive,mode='directed',diag=FALSE))
  pos_mtx_sum[5,i] <- transitivity(graph_from_adjacency_matrix(networks_mtx[[i]]$positive,mode='directed',diag=FALSE))
  pos_mtx_sum[6,i] <- sum(sna::degree(networks_mtx[[i]]$positive,cmode='freeman') == 0)
  pos_mtx_sum[7,i] <- mean(rowSums(networks_mtx[[i]]$positive,na.rm=TRUE))
  pos_mtx_sum[8,i] <- sd(rowSums(networks_mtx[[i]]$positive,na.rm=TRUE))
  pos_mtx_sum[9,i] <- sd(colSums(networks_mtx[[i]]$positive,na.rm=TRUE))
  # negative ties
  neg_mtx_sum[1,i] <- nrow(networks_mtx[[i]]$negative)
  neg_mtx_sum[2,i] <- 1-sum(!is.na(networks_mtx[[i]]$negative))/(nrow(networks_mtx[[i]]$negative)*(nrow(networks_mtx[[i]]$negative)-1))
  neg_mtx_sum[3,i] <- edge_density(graph_from_adjacency_matrix(networks_mtx[[i]]$negative,mode='directed',diag=FALSE))
  neg_mtx_sum[4,i] <- reciprocity(graph_from_adjacency_matrix(networks_mtx[[i]]$negative,mode='directed',diag=FALSE))
  neg_mtx_sum[5,i] <- transitivity(graph_from_adjacency_matrix(networks_mtx[[i]]$negative,mode='directed',diag=FALSE))
  neg_mtx_sum[6,i] <- sum(sna::degree(networks_mtx[[i]]$negative,cmode='freeman') == 0)
  neg_mtx_sum[7,i] <- mean(rowSums(networks_mtx[[i]]$negative,na.rm=TRUE))
  neg_mtx_sum[8,i] <- sd(rowSums(networks_mtx[[i]]$negative,na.rm=TRUE)) 
  neg_mtx_sum[9,i] <- sd(colSums(networks_mtx[[i]]$negative,na.rm=TRUE)) 
}
write.table(round(pos_mtx_sum,3),'pos_ties.csv',row.names=TRUE,sep=';')
write.table(round(neg_mtx_sum,3),'neg_ties.csv',row.names=TRUE,sep=';')

# Out- and in-degree distributions (positive ties)
par(mfrow=c(3,2))

for(i in seq_along(networks_mtx)){
  hist(rowSums(networks_mtx[[i]]$positive,na.rm=TRUE),main=paste('Out-degrees in unit ',i,sep=''))
}
for(i in seq_along(networks_mtx)){
  hist(colSums(networks_mtx[[i]]$positive,na.rm=TRUE),main=paste('In-degrees in unit ',i,sep=''))
}

########################################################################################################################

# 2) DESCRIPTIVE STATS (POSITIVE & NEGATIVE GOSSIP)

# Projections of the positive and negative gossip cubes (sender-receiver, sender-target, and receiver-target)
gos_pos_sr <- gos_pos_st <- gos_pos_rt <- gos_pos
gos_neg_sr <- gos_neg_st <- gos_neg_rt <- gos_neg

for(i in seq_along(gos_pos)){
  gos_pos_sr[[i]] <- apply(gos_pos[[i]],c(1,2),max,na.rm=TRUE)
  gos_pos_st[[i]] <- apply(gos_pos[[i]],c(1,3),max,na.rm=TRUE)
  gos_pos_rt[[i]] <- apply(gos_pos[[i]],c(2,3),max,na.rm=TRUE)
  gos_neg_sr[[i]] <- apply(gos_neg[[i]],c(1,2),max,na.rm=TRUE)
  gos_neg_st[[i]] <- apply(gos_neg[[i]],c(1,3),max,na.rm=TRUE)
  gos_neg_rt[[i]] <- apply(gos_neg[[i]],c(2,3),max,na.rm=TRUE)
}

gossip_ntw <- list(pos_sr=gos_pos_sr,pos_st=gos_pos_st,pos_rt=gos_pos_rt,
                   neg_sr=gos_neg_sr,neg_st=gos_neg_st,neg_rt=gos_neg_rt)

# NA allocation
for(i in seq_along(gossip_ntw)){
  for(j in seq_along(gossip_ntw[[i]])){
    gossip_ntw[[i]][[j]][is.infinite(gossip_ntw[[i]][[j]])] <- NA
  }
}

# SUmmary table
gossip_sum_ind <- matrix(NA,12,6,dimnames=list(c('pos_s','pos_r','pos_t','neg_s','neg_r','neg_t',
                                                'pos_sr','neg_sr','pos_st','neg_st','pos_rt','neg_rt'),organisation_ID))

for(i in 1:ncol(gossip_sum_ind)){
  # Number of positive gossip senders, receivers, and targets per unit
  gossip_sum_ind[1,i] <- sum(rowSums(gossip_ntw$pos_sr[[i]],na.rm = TRUE) != 0)
  gossip_sum_ind[2,i] <- sum(rowSums(gossip_ntw$pos_rt[[i]],na.rm = TRUE) != 0)
  gossip_sum_ind[3,i] <- sum(colSums(gossip_ntw$pos_rt[[i]],na.rm = TRUE) != 0)
  # Number of negative gossip senders, receivers, and targets per unit
  gossip_sum_ind[4,i] <- sum(rowSums(gossip_ntw$neg_sr[[i]],na.rm = TRUE) != 0)
  gossip_sum_ind[5,i] <- sum(rowSums(gossip_ntw$neg_rt[[i]],na.rm = TRUE) != 0)
  gossip_sum_ind[6,i] <- sum(colSums(gossip_ntw$neg_rt[[i]],na.rm = TRUE) != 0)
  # Number of gossip ties per gossip projection (sr, st, rt)
  gossip_sum_ind[7,i] <- sum(gossip_ntw$pos_sr[[i]],na.rm = TRUE)
  gossip_sum_ind[8,i] <- sum(gossip_ntw$neg_sr[[i]],na.rm = TRUE)
  gossip_sum_ind[9,i] <- sum(gossip_ntw$pos_st[[i]],na.rm = TRUE)
  gossip_sum_ind[10,i] <- sum(gossip_ntw$neg_st[[i]],na.rm = TRUE)
  gossip_sum_ind[11,i] <- sum(gossip_ntw$pos_rt[[i]],na.rm = TRUE)
  gossip_sum_ind[12,i] <- sum(gossip_ntw$neg_rt[[i]],na.rm = TRUE)
}

gossip_sum_ind

########################################################################################################################

# 3) BIVARIATE STATS (NETWORK LEVEL: JACCARD INDICES FOR MATRIX OVERLAP)

# Dyadic: Jaccard indices (matrix overlap)
jacc_ind <- array(NA,dim=c(6,2,6),
                  dimnames=list(c('pos_sr','pos_st','pos_rt','neg_sr','neg_st','neg_rt'),
                                c('positive','negative'),
                                organisation_ID))

for(i in seq_along(organisation_ID)){
  jacc_ind[1,1,i] <- round(Jaccard(gossip_ntw$pos_sr[[i]],networks_mtx[[i]]$positive),3)*100
  jacc_ind[2,1,i] <- round(Jaccard(gossip_ntw$pos_st[[i]],networks_mtx[[i]]$positive),3)*100
  jacc_ind[3,1,i] <- round(Jaccard(gossip_ntw$pos_rt[[i]],networks_mtx[[i]]$positive),3)*100
  jacc_ind[4,1,i] <- round(Jaccard(gossip_ntw$neg_sr[[i]],networks_mtx[[i]]$positive),3)*100
  jacc_ind[5,1,i] <- round(Jaccard(gossip_ntw$neg_st[[i]],networks_mtx[[i]]$positive),3)*100
  jacc_ind[6,1,i] <- round(Jaccard(gossip_ntw$neg_rt[[i]],networks_mtx[[i]]$positive),3)*100
  jacc_ind[1,2,i] <- round(Jaccard(gossip_ntw$pos_sr[[i]],networks_mtx[[i]]$negative),3)*100
  jacc_ind[2,2,i] <- round(Jaccard(gossip_ntw$pos_st[[i]],networks_mtx[[i]]$negative),3)*100
  jacc_ind[3,2,i] <- round(Jaccard(gossip_ntw$pos_rt[[i]],networks_mtx[[i]]$negative),3)*100
  jacc_ind[4,2,i] <- round(Jaccard(gossip_ntw$neg_sr[[i]],networks_mtx[[i]]$negative),3)*100
  jacc_ind[5,2,i] <- round(Jaccard(gossip_ntw$neg_st[[i]],networks_mtx[[i]]$negative),3)*100
  jacc_ind[6,2,i] <- round(Jaccard(gossip_ntw$neg_rt[[i]],networks_mtx[[i]]$negative),3)*100
}

jacc_ind
round(apply(jacc_ind,c(1,2),mean),1); apply(jacc_ind,c(1,2),min); apply(jacc_ind,c(1,2),max)

########################################################################################################################

# 4) CLUSTER DETECTION (Louvain method)

for(i in seq_along(networks_mtx)){
  networks_mtx[[i]]$clust_louv <- cluster_louvain(graph_from_adjacency_matrix(networks_mtx[[i]]$positive,
                                                                              mode='undirected',diag=FALSE))
  networks_mtx[[i]]$clusters <- 1*outer(networks_mtx[[i]]$clust_louv$membership,networks_mtx[[i]]$clust_louv$membership,
                                        FUN='==')
  rownames(networks_mtx[[i]]$clusters) <- colnames(networks_mtx[[i]]$clusters) <- rownames(networks_mtx[[i]]$positive)
}

# Creation of triadic data to check how many gossip is intra- vs. inter-group
# Creation of large data set N(N-1)(N-2)
triad_data <- vector('list',length=length(org_subject))
names(triad_data) <- names(org_subject)

for(i in seq_along(triad_data)){
  triad_data[[i]] <- data.frame(matrix(NA,nrow=length(org_subject[[i]])^3,ncol=5))
  colnames(triad_data[[i]]) <- c('unit','sender','receiver','target','gossip')
  triad_data[[i]]$unit <- organisation_ID[[i]]
  triad_data[[i]]$sender <- rep(org_subject[[i]],each=(length(org_subject[[i]]))^2)
  triad_data[[i]]$receiver <- rep(org_subject[[i]],times=length(org_subject[[i]]),each=length(org_subject[[i]]))
  triad_data[[i]]$target <- rep(org_subject[[i]],times=(length(org_subject[[i]])^2))
  # Exclusion of diagonals (if sender = receiver, sender = target, or receiver = target)
  triad_data[[i]] <- triad_data[[i]][triad_data[[i]]$sender != triad_data[[i]]$receiver & 
                                       triad_data[[i]]$sender != triad_data[[i]]$target & 
                                       triad_data[[i]]$receiver != triad_data[[i]]$target,]
}

# Extraction of gossip triads from the gossip cube
for(i in seq_along(gos_pos)){
  for(s in rownames(gos_pos[[i]])){
    for(r in rownames(gos_pos[[i]])){
      for(t in rownames(gos_pos[[i]])){
        if(s != r & s != t & r != t){
          triad_data[[i]][triad_data[[i]]$sender == s & triad_data[[i]]$receiver == r & triad_data[[i]]$target == t,
          ]$gossip <- gos_pos[[i]][s,r,t] - gos_neg[[i]][s,r,t]
        }
      }
    }
  }
}

# Extraction of positive and negative ties, and whether respondents belong to the same group or not
for(i in seq_along(triad_data)){
  for(j in 1:nrow(triad_data[[i]])){
    # Positive and negative ties
    triad_data[[i]]$SR[j] <- networks_mtx[[i]]$positive[triad_data[[i]]$sender[j],triad_data[[i]]$receiver[j]] - 
      networks_mtx[[i]]$negative[triad_data[[i]]$sender[j],triad_data[[i]]$receiver[j]]
    triad_data[[i]]$ST[j] <- networks_mtx[[i]]$positive[triad_data[[i]]$sender[j],triad_data[[i]]$target[j]] - 
      networks_mtx[[i]]$negative[triad_data[[i]]$sender[j],triad_data[[i]]$target[j]]
    triad_data[[i]]$RT[j] <- networks_mtx[[i]]$positive[triad_data[[i]]$receiver[j],triad_data[[i]]$target[j]] - 
      networks_mtx[[i]]$negative[triad_data[[i]]$receiver[j],triad_data[[i]]$target[j]]
    # Group membership
    triad_data[[i]]$samegroup_SR[j] <- networks_mtx[[i]]$clusters[triad_data[[i]]$sender[j],triad_data[[i]]$receiver[j]]
    triad_data[[i]]$samegroup_ST[j] <- networks_mtx[[i]]$clusters[triad_data[[i]]$sender[j],triad_data[[i]]$target[j]]
    triad_data[[i]]$samegroup_RT[j] <- networks_mtx[[i]]$clusters[triad_data[[i]]$receiver[j],triad_data[[i]]$target[j]]
  }
}

########################################################################################################################

# 5) BROKERS' DETECTION

# Creation of matrices with two roles: broker (betweenness-central actors) and non-broker
for(i in seq_along(networks_mtx)){
  # Edge betweenness
  networks_mtx[[i]]$sqrt_btw <- betweenness(graph_from_adjacency_matrix(networks_mtx[[i]]$positive,
                                                                        mode='directed',diag=FALSE),directed=TRUE) 
  # Distance object
  networks_mtx[[i]]$sqrt_btw <- as.dist(abs(outer(networks_mtx[[i]]$sqrt_btw,networks_mtx[[i]]$sqrt_btw,'-')))
  # Hierarchical clustering (Ward D method)
  networks_mtx[[i]]$sqrt_btw <- hclust(networks_mtx[[i]]$sqrt_btw,method='ward.D')
  plot(networks_mtx[[i]]$sqrt_btw)
  rect.hclust(networks_mtx[[i]]$sqrt_btw,2,border='blue')
  # Extraction of roles
  networks_mtx[[i]]$roles <- cutree(networks_mtx[[i]]$sqrt_btw,k=2)
}

# Labelling the roles
networks_mtx[[1]]$roles <- factor(networks_mtx[[1]]$roles,levels=c(1,2),labels=c('broker','nonbroker'))
networks_mtx[[2]]$roles <- factor(networks_mtx[[2]]$roles,levels=c(1,2),labels=c('broker','nonbroker'))
networks_mtx[[3]]$roles <- factor(networks_mtx[[3]]$roles,levels=c(2,1),labels=c('broker','nonbroker'))
networks_mtx[[4]]$roles <- factor(networks_mtx[[4]]$roles,levels=c(1,2),labels=c('broker','nonbroker'))
networks_mtx[[5]]$roles <- factor(networks_mtx[[5]]$roles,levels=c(2,1),labels=c('broker','nonbroker'))
networks_mtx[[6]]$roles <- factor(networks_mtx[[6]]$roles,levels=c(2,1),labels=c('broker','nonbroker'))

# Extraction of roles for the triadic data
for(i in seq_along(triad_data)){
  triad_data[[i]]$target_role <- triad_data[[i]]$receiver_role <- triad_data[[i]]$sender_role <- NA
}

for(i in seq_along(networks_mtx)){
  for(j in names(networks_mtx[[i]]$roles)){
    triad_data[[i]]$sender_role[triad_data[[i]]$sender == j] <- as.character(networks_mtx[[i]]$roles[j])
    triad_data[[i]]$receiver_role[triad_data[[i]]$receiver == j] <- as.character(networks_mtx[[i]]$roles[j])
    triad_data[[i]]$target_role[triad_data[[i]]$target == j] <- as.character(networks_mtx[[i]]$roles[j])
  }
}

# Creation of matrices: broker-broker, broker-non-broker, etc.
for(x in seq_along(networks_mtx)){
  networks_mtx[[x]]$cc <- networks_mtx[[x]]$cp <- networks_mtx[[x]]$pc <- networks_mtx[[x]]$pp <- networks_mtx[[x]]$clusters*0
  for(i in seq_along(networks_mtx[[x]]$roles)){
    for(j in seq_along(networks_mtx[[x]]$roles)){
      if(networks_mtx[[x]]$roles[i] == 'broker' & networks_mtx[[x]]$roles[j] == 'broker'){
        networks_mtx[[x]]$cc[i,j] <- 1
      }
      if(networks_mtx[[x]]$roles[i] == 'broker' & networks_mtx[[x]]$roles[j] == 'nonbroker'){
        networks_mtx[[x]]$cp[i,j] <- 1
      }
      if(networks_mtx[[x]]$roles[i] == 'nonbroker' & networks_mtx[[x]]$roles[j] == 'broker'){
        networks_mtx[[x]]$pc[i,j] <- 1
      }
      if(networks_mtx[[x]]$roles[i] == 'nonbroker' & networks_mtx[[x]]$roles[j] == 'nonbroker'){
        networks_mtx[[x]]$pp[i,j] <- 1
      }
    }
  }
}

########################################################################################################################

# Visualisation of the networks with communities and core-peripheral nodes detected 
ntw_plot <- networks_mtx

for(i in seq_along(ntw_plot)){
  # Clustering using Newman's method based on positive ties only
  ntw_plot[[i]]$group <- ntw_plot[[i]]$positive 
  ntw_plot[[i]]$group <- graph_from_adjacency_matrix(ntw_plot[[i]]$group,mode='undirected',diag=FALSE)
  ntw_plot[[i]]$group <- igraph::cluster_louvain(ntw_plot[[i]]$group) # clustering (Louvain method)
  # Show both positive and negative ties
  ntw_plot[[i]]$vis <- graph_from_adjacency_matrix(ntw_plot[[i]]$positive,mode='directed',diag=FALSE)
  # Layout based only on positive ties
  set.seed(0708)
  ntw_plot[[i]]$layout <- layout_with_fr(graph_from_adjacency_matrix(ntw_plot[[i]]$positive,mode='directed'))
  # Customisation of nodes and ties
  V(ntw_plot[[i]]$vis)$shape <- ifelse(attributes[attributes$group == organisation_ID[[i]],]$hr_leader == 1,'square','circle')
  V(ntw_plot[[i]]$vis)$color <- ifelse(networks_mtx[[i]]$roles == 'broker','red','yellow')
  V(ntw_plot[[i]]$vis)$size <- ifelse(networks_mtx[[i]]$roles == 'broker',10,6)
  E(ntw_plot[[i]]$vis)$color <- 'darkgreen'
  E(ntw_plot[[i]]$vis)$width <- 1
  # Visualisation
  jpeg(filename=paste('Unit',i,'.jpeg',sep=''),width=4,height=4,units='in',res=1000)
  plot(ntw_plot[[i]]$vis,mark.groups=ntw_plot[[i]]$group,
       edge.arrow.size=.2,vertex.label=NA,layout=ntw_plot[[i]]$layout)
  dev.off()
}

########################################################################################################################

# 6) DESCRIPTIVE STATISTICS

# 6.1) Positive gossip
# Keep only the gossip triads
CP_descrip <- triad_data

for(i in seq_along(CP_descrip)){
  CP_descrip[[i]] <- CP_descrip[[i]][!is.na(CP_descrip[[i]]$gossip) & CP_descrip[[i]]$gossip == 1,]
}

com_descrip <- CP_descrip

# In-group vs. out-group: descriptive (SR, ST, RT)
for(i in seq_along(com_descrip)){
  com_descrip[[i]] <- as.vector(table(com_descrip[[i]]$samegroup_SR,
                                      com_descrip[[i]]$samegroup_ST,
                                      com_descrip[[i]]$samegroup_RT))
}
com_descrip <- do.call('rbind',com_descrip)
com_descrip <- com_descrip[,c(8,4,6,2,7,3,5,1)]
colnames(com_descrip) <- c('III','II0','IOI','IOO','OII','OIO','OOI','OOO')
com_descrip <- com_descrip[,c('III','IOO','OIO','OOI','OOO')]
com_descrip
round(colSums(com_descrip)/sum(com_descrip)*100,2)

# Centre-periphery: descriptive (sender, receiver, target)
for(i in seq_along(CP_descrip)){
  CP_descrip[[i]] <- as.vector(table(CP_descrip[[i]]$sender_role,
                                     CP_descrip[[i]]$receiver_role,
                                     CP_descrip[[i]]$target_role))
}
CP_descrip <- do.call('rbind',CP_descrip)
CP_descrip <- CP_descrip[,c(1,5,3,7,2,6,4,8)]
colnames(CP_descrip) <- c('BBB','BBN','BNB','BNN','NBB','NBN','NNB','NNN')
CP_descrip
round(colSums(CP_descrip)/sum(CP_descrip)*100,2)

# 6.2) Negative gossip
# Keep only the gossip triads
CP_descrip <- triad_data

for(i in seq_along(CP_descrip)){
  CP_descrip[[i]] <- CP_descrip[[i]][!is.na(CP_descrip[[i]]$gossip) & CP_descrip[[i]]$gossip == -1,]
}

com_descrip <- CP_descrip

# In-group vs. out-group: descriptive (SR, ST, RT)
for(i in seq_along(com_descrip)){
  com_descrip[[i]] <- as.vector(table(com_descrip[[i]]$samegroup_SR,
                                      com_descrip[[i]]$samegroup_ST,
                                      com_descrip[[i]]$samegroup_RT))
}
com_descrip <- do.call('rbind',com_descrip)
com_descrip <- com_descrip[,c(8,4,6,2,7,3,5,1)]
colnames(com_descrip) <- c('III','II0','IOI','IOO','OII','OIO','OOI','OOO')
com_descrip <- com_descrip[,c('III','IOO','OIO','OOI','OOO')]
com_descrip
round(colSums(com_descrip)/sum(com_descrip)*100,2)

# Centre-periphery: descriptive (sender, receiver, target)
for(i in seq_along(CP_descrip)){
  CP_descrip[[i]] <- as.vector(table(CP_descrip[[i]]$sender_role,
                                     CP_descrip[[i]]$receiver_role,
                                     CP_descrip[[i]]$target_role))
}
CP_descrip <- do.call('rbind',CP_descrip)
CP_descrip <- CP_descrip[,c(1,5,3,7,2,6,4,8)]
colnames(CP_descrip) <- c('BBB','BBN','BNB','BNN','NBB','NBN','NNB','NNN')
CP_descrip
round(colSums(CP_descrip)/sum(CP_descrip)*100,2)

########################################################################################################################

# 7) BROKERS' PROPERTIES

brokers <- pos_odeg <- pos_ideg <- neg_odeg <- neg_ideg <- networks_mtx

for(i in seq_along(brokers)){
  brokers[[i]] <- brokers[[i]]$roles
  pos_odeg[[i]] <- rowSums(pos_odeg[[i]]$positive,na.rm=TRUE)
  pos_ideg[[i]] <- colSums(pos_ideg[[i]]$positive,na.rm=TRUE)
  neg_odeg[[i]] <- rowSums(neg_odeg[[i]]$negative,na.rm=TRUE)
  neg_ideg[[i]] <- colSums(neg_ideg[[i]]$negative,na.rm=TRUE)
}

brokers <- c(brokers[[1]],brokers[[2]],brokers[[3]],brokers[[4]],brokers[[5]],brokers[[6]])
pos_odeg <- c(pos_odeg[[1]],pos_odeg[[2]],pos_odeg[[3]],pos_odeg[[4]],pos_odeg[[5]],pos_odeg[[6]])
pos_ideg <- c(pos_ideg[[1]],pos_ideg[[2]],pos_ideg[[3]],pos_ideg[[4]],pos_ideg[[5]],pos_ideg[[6]])
neg_odeg <- c(neg_odeg[[1]],neg_odeg[[2]],neg_odeg[[3]],neg_odeg[[4]],neg_odeg[[5]],neg_odeg[[6]])
neg_ideg <- c(neg_ideg[[1]],neg_ideg[[2]],neg_ideg[[3]],neg_ideg[[4]],neg_ideg[[5]],neg_ideg[[6]])


brokers <- factor(brokers,levels=c(1,2),labels=c('broker','nonbroker'))
brokers <- as.data.frame(brokers)
brokers$pos_odeg <- pos_odeg
brokers$pos_ideg <- pos_ideg
brokers$neg_odeg <- neg_odeg
brokers$neg_ideg <- neg_ideg
brokers$ID <- rownames(brokers)

attributes <- merge(x=attributes,y=brokers,by.x='responder',by.y='ID',all.x=TRUE)

# Descriptive stats
summary(attributes$brokers)

# Female brokers
sum(attributes[attributes$brokers == 'broker',]$woman)
sum(attributes[attributes$brokers == 'nonbroker',]$woman)
prop.test(x=c(8,33),n=c(28,100),alternative='two.sided')

# Brokers in management
sum(attributes[attributes$brokers == 'broker',]$hr_leader)
sum(attributes[attributes$brokers == 'nonbroker',]$hr_leader)
prop.test(x=c(16,18),n=c(28,100),alternative='two.sided')

# Age of brokers
t.test(x=(2018 - attributes[attributes$brokers == 'broker',]$birth_year),
       y=(2018 - attributes[attributes$brokers == 'nonbroker',]$birth_year),alternative='two.sided')
range(2018 - attributes[attributes$brokers == 'broker',]$birth_year)
range(2018 - attributes[attributes$brokers == 'nonbroker',]$birth_year)

# Tenure of brokers
t.test(x=(attributes[attributes$brokers == 'broker',]$tenure),
       y=(attributes[attributes$brokers == 'nonbroker',]$tenure),alternative='two.sided')
range(attributes[attributes$brokers == 'broker',]$tenure,na.rm=TRUE)
range(attributes[attributes$brokers == 'nonbroker',]$tenure,na.rm=TRUE)

# Degree differences
# Positive outdegree
t.test(x=(attributes[attributes$brokers == 'broker',]$pos_odeg),
       y=(attributes[attributes$brokers == 'nonbroker',]$pos_odeg),alternative='two.sided')
range(attributes[attributes$brokers == 'broker',]$pos_odeg,na.rm=TRUE)
range(attributes[attributes$brokers == 'nonbroker',]$pos_odeg,na.rm=TRUE)
# Positive indegree
t.test(x=(attributes[attributes$brokers == 'broker',]$pos_ideg),
       y=(attributes[attributes$brokers == 'nonbroker',]$pos_ideg),alternative='two.sided')
range(attributes[attributes$brokers == 'broker',]$pos_ideg,na.rm=TRUE)
range(attributes[attributes$brokers == 'nonbroker',]$pos_ideg,na.rm=TRUE)
# Negative outdegree
t.test(x=(attributes[attributes$brokers == 'broker',]$neg_odeg),
       y=(attributes[attributes$brokers == 'nonbroker',]$neg_odeg),alternative='two.sided')
range(attributes[attributes$brokers == 'broker',]$neg_odeg,na.rm=TRUE)
range(attributes[attributes$brokers == 'nonbroker',]$neg_odeg,na.rm=TRUE)
# Negative indegree
t.test(x=(attributes[attributes$brokers == 'broker',]$neg_ideg),
       y=(attributes[attributes$brokers == 'nonbroker',]$neg_ideg),alternative='two.sided')
range(attributes[attributes$brokers == 'broker',]$neg_ideg,na.rm=TRUE)
range(attributes[attributes$brokers == 'nonbroker',]$neg_ideg,na.rm=TRUE)

########################################################################################################################

# Removal of unnecessary objects
rm(pos_mtx_sum);rm(neg_mtx_sum);rm(gos_pos_sr);rm(gos_pos_st);rm(gos_pos_rt);rm(gos_neg_sr);rm(gos_neg_st);rm(gos_neg_rt)
rm(gossip_sum_ind);rm(jacc_ind);rm(i);rm(j);rm(ntw_plot);rm(s);rm(r);rm(t);rm(x);rm(com_descrip);rm(CP_descrip)
rm(brokers);rm(pos_odeg);rm(pos_ideg);rm(neg_odeg);rm(neg_ideg)

# Save image
save.image('modellingdata.RData')
########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## Data tidying (1)
## R script written by Jose Luis Estevez (Linkoping University)
## Date: December 22nd, 2020
########################################################################################################################

# R PACKAGES REQUIRED
library(ggplot2);library(sna)

# DATA LOADING AND DATA TIDYING
load('data.RData')
rm(list=ls()[-c(1,13,20)]) # keep only attributes, networks, and gossip
attributes <- attributes[,1:75] # # exclusion of network variables from attributes

########################################################################################################################

# Some general information
organisation_ID <- unique(na.omit(attributes$group)) # 9 units or groups
employee_ID <- attributes$responder # N=225, but observe there are 2 duplicated employees: 1F6032 and 1F6033
networks_available <- unique(networks$network_questionID_EN) # 43 networks available in the data set

# Fixing duplicated employees in all three objects: attributes, networks, and gossip
'%!in%' <- function(x,y)!('%in%'(x,y))
attributes <- attributes[attributes$responder %!in% c('1F6032_F106b','1F6033_F106b'),] # N=223
attributes$responder[122] <- '1F6032'
attributes$responder[123] <- '1F6033'

duplicated_actors <- c('1F6032_F106a','1F6032_F106b','1F6033_F106a','1F6033_F106b')

for(i in 1:nrow(networks)){
  if(networks$sender[i] %in% duplicated_actors){
    networks$sender[i] <- substr(networks$sender[i],1,6)
  }
  if(networks$receiver[i] %in% duplicated_actors){
    networks$receiver[i] <- substr(networks$receiver[i],1,6)
  }
}

for(i in 1:nrow(gossip)){
  if(gossip$sender[i] %in% duplicated_actors){
    gossip$sender[i] <- substr(gossip$sender[i],1,6)
  }
  if(gossip$receiver[i] %in% duplicated_actors){
    gossip$receiver[i] <- substr(gossip$receiver[i],1,6)
  }
  if(gossip$target[i] %in% duplicated_actors){
    gossip$target[i] <- substr(gossip$target[i],1,6)
  }
}

# Respondents missing
apply(attributes[,],2,is.na)*1 -> missing_data
attributes$respondent_missing <- as.vector(apply(missing_data,1,sum)) # items non-responded
attributes$respondent_missing <- attributes$respondent_missing > 50
(missing_respondents <- attributes[attributes$respondent_missing == TRUE,'responder']) # 18 non-respondents

########################################################################################################################

# CREATION OF GOSSIP-CUBES

# Subjects' IDs by organisation: for the dimensions of matrices and cubes
org_subject <- vector('list',length=length(organisation_ID))
names(org_subject) <- organisation_ID

for(i in seq_along(org_subject)){
  org_subject[[i]] <- attributes[attributes$group == names(org_subject)[i],]$responder
}

# Exclusion of gossip triples when not all three parties (sender, receiver, target) are part of the same network
valid_triplets <- vector(length=nrow(gossip))

for(i in 1:nrow(gossip)){
  if(gossip$sender[i] %in% org_subject[[gossip$group[[i]]]] &
     gossip$receiver[i] %in% org_subject[[gossip$group[[i]]]] &
     gossip$target[i] %in% org_subject[[gossip$group[[i]]]])
  valid_triplets[i] <- TRUE
}
gossip <- gossip[valid_triplets,]

# Any duplicated?
any(duplicated(gossip[,c('sender','receiver','target')],)) 
which(duplicated(gossip[,c('sender','receiver','target')],))
gossip[c(598,599,676,677,770,771),] # neutral gossips to positive or negative
gossip <- gossip[-c(599,677,770),]

# Creation of gossip cubes (sender, receiver, target)
gossip_cube <- org_subject
for(i in seq_along(gossip_cube)){
  gossip_cube[[i]] <- array(NA,dim=c(rep(length(org_subject[[i]]),3)),
                            dimnames=list(org_subject[[i]],org_subject[[i]],org_subject[[i]]))
}

# Allocation of the gossip triplets in the cubes
for(i in 1:nrow(gossip)){
  gossip_cube[[gossip$group[i]]][gossip$sender[i],gossip$receiver[i],gossip$target[i]] <- gossip$info[i]
}

########################################################################################################################

# CREATION OF RELATIONAL MATRICES

# Exclusion of ties between different networks
valid_ties <- vector(length=nrow(networks))

for(i in 1:nrow(networks)){
  if(networks$sender[i] %in% org_subject[[networks$group[[i]]]] &
     networks$receiver[i] %in% org_subject[[networks$group[[i]]]])
    valid_ties[i] <- TRUE
}
networks <- networks[valid_ties,]

# Creation of matrices (relational data)
networks_mtx <- org_subject

for(i in seq_along(networks_mtx)){
  networks_mtx[[i]] <- vector('list',length=length(networks_available))
  names(networks_mtx[[i]]) <- networks_available
  for(j in seq_along(networks_mtx[[i]])){
    networks_mtx[[i]][[j]] <- array(0,dim=c(rep(length(org_subject[[i]]),2)),
                                    dimnames=list(org_subject[[i]],org_subject[[i]]))
    # Allocation of missing data
    diag(networks_mtx[[i]][[j]]) <- NA # diagonal
    for(k in 1:nrow(networks_mtx[[i]][[j]])){
      if(rownames(networks_mtx[[i]][[j]])[k] %in% missing_respondents){
        networks_mtx[[i]][[j]][k,] <- NA # respondents missing
      }
    }
  }
}

# Allocation of the ties into the matrices
networks$tie <- 1
for(i in 1:nrow(networks)){
  gro <- networks$group[i]
  mtx <- networks$network_questionID_EN[i]
  row <- networks$sender[i]
  col <- networks$receiver[i]
  networks_mtx[[gro]][[mtx]][row,col] <- networks$tie[i]
}

rm(duplicated_actors);rm(missing_data);rm(valid_triplets);rm(valid_ties);rm(employee_ID);rm(networks)
rm(i);rm(j);rm(gro);rm(mtx);rm(row);rm(col)

########################################################################################################################

# CHECKING THE GOSSIP CUBE

gos <- gos_pos <- gos_neg <- gos_neu <- gossip_cube

for(x in seq_along(gossip_cube)){
  for(i in 1:nrow(gossip_cube[[x]])){
    for(j in 1:nrow(gossip_cube[[x]])){
      for(k in 1:nrow(gossip_cube[[x]])){
        if(i != j & i != k & j != k){
          # the gossip cube (no matter the tone)
          if(!is.na(gossip_cube[[x]][i,j,k])){
            gos[[x]][i,j,k] <- 1
          }else{
            gos[[x]][i,j,k] <- 0
          }
          # the positive gossip cube
          if(!is.na(gossip_cube[[x]][i,j,k]) & gossip_cube[[x]][i,j,k] == 1){
            gos_pos[[x]][i,j,k] <- 1
          }else{
            gos_pos[[x]][i,j,k] <- 0
          }
          # the negative gossip cube
          if(!is.na(gossip_cube[[x]][i,j,k]) & gossip_cube[[x]][i,j,k] == -1){
            gos_neg[[x]][i,j,k] <- 1
          }else{
            gos_neg[[x]][i,j,k] <- 0
          }
          # the neutral gossip cube
          if(!is.na(gossip_cube[[x]][i,j,k]) & gossip_cube[[x]][i,j,k] == 0){
            gos_neu[[x]][i,j,k] <- 1
          }else{
            gos_neu[[x]][i,j,k] <- 0
          }
        }
      }
    }
  }
}

# How many triad in the data set? How many participants reported gossip?
rec_sen <- gos
rec_sen_pos <- gos_pos
rec_sen_neg <- gos_neg
rec_sen_neu <- gos_neu

for(x in seq_along(rec_sen)){
  rec_sen[[x]] <- apply(rec_sen[[x]],c(2,1),max,na.rm=TRUE) # receiver-sender matrix
  rec_sen_pos[[x]] <- apply(rec_sen_pos[[x]],c(2,1),max,na.rm=TRUE)
  rec_sen_neg[[x]] <- apply(rec_sen_neg[[x]],c(2,1),max,na.rm=TRUE)
  rec_sen_neu[[x]] <- apply(rec_sen_neu[[x]],c(2,1),max,na.rm=TRUE)
  diag(rec_sen[[x]]) <- 0
  diag(rec_sen_pos[[x]]) <- 0
  diag(rec_sen_neg[[x]]) <- 0
  diag(rec_sen_neu[[x]]) <- 0
  rec_sen[[x]] <- rowSums(rec_sen[[x]])
  rec_sen_pos[[x]] <- rowSums(rec_sen_pos[[x]])
  rec_sen_neg[[x]] <- rowSums(rec_sen_neg[[x]])
  rec_sen_neu[[x]] <- rowSums(rec_sen_neu[[x]])
}

gossip_sum <- as.data.frame(matrix(NA,nrow=length(rec_sen),ncol=10,
                                   dimnames=list(names(rec_sen),c('triads','triads_pos','triads_neg','triads_neu',
                                                                  'triads_total','reporters','reporters_pos',
                                                                  'reporters_neg','reporters_neu','reporters_total'))))
for(x in seq_along(rec_sen)){
  # summary of triads: with gossip, positive or negative, and total possible
  gossip_sum[x,'triads'] <- sum(gos[[x]],na.rm=TRUE)
  gossip_sum[x,'triads_pos'] <- sum(gos_pos[[x]],na.rm=TRUE)
  gossip_sum[x,'triads_neg'] <- sum(gos_neg[[x]],na.rm=TRUE)
  gossip_sum[x,'triads_neu'] <- sum(gos_neu[[x]],na.rm=TRUE)
  gossip_sum[x,'triads_total'] <- nrow(gos[[x]]) * (nrow(gos[[x]])-1) * (nrow(gos[[x]])-2)
  # number of respondents reporting those triads
  gossip_sum[x,'reporters'] <- sum(rec_sen[[x]] != 0)
  gossip_sum[x,'reporters_pos'] <- sum(rec_sen_pos[[x]] != 0)
  gossip_sum[x,'reporters_neg'] <- sum(rec_sen_neg[[x]] != 0)
  gossip_sum[x,'reporters_neu'] <- sum(rec_sen_neu[[x]] != 0)
  gossip_sum[x,'reporters_total'] <- length(rec_sen[[x]])
}

View(gossip_sum)

########################################################################################################################

# SOME DESCRIPTIVES OF THE GOSSIP DATA

desc_receiver <- vector('list',length=length(gossip_cube))
for(i in seq_along(gossip_cube)){
  desc_receiver[[i]] <- as.data.frame(matrix(NA,nrow=nrow(gossip_cube[[i]]),ncol=5))
  rownames(desc_receiver[[i]]) <- rownames(gossip_cube[[i]])
  colnames(desc_receiver[[i]]) <- c('network','dyads','dyads_pos','dyads_neg','dyads_neu')
  desc_receiver[[i]]$network <- names(gossip_cube)[i]
}

for(x in seq_along(desc_receiver)){
  for(i in 1:nrow(desc_receiver[[x]])){
    # Receiver-specific 
    desc_receiver[[x]]$dyads[i] <- sum(gos[[x]][,i,],na.rm=TRUE)
    desc_receiver[[x]]$dyads_pos[i] <- sum(gos_pos[[x]][,i,],na.rm=TRUE)
    desc_receiver[[x]]$dyads_neg[i] <- sum(gos_neg[[x]][,i,],na.rm=TRUE)
    desc_receiver[[x]]$dyads_neu[i] <- sum(gos_neu[[x]][,i,],na.rm=TRUE)
  }
}

desc_receiver <- do.call('rbind',desc_receiver)
desc_receiver$subject <- rownames(desc_receiver)

# Visualisations
grid.background <- theme_bw()+
  theme(plot.background=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+
  theme(axis.line=element_line(color='black'))+
  theme(strip.text.x=element_text(colour='white',face='bold'))+
  theme(strip.background=element_rect(fill='black'))+
  theme(axis.text.x = element_blank())

jpeg(filename='Receiver-specific gossip.jpeg',width=11,height=6,units='in',res=1000)
ggplot(data=desc_receiver)+
  geom_point(aes(x=subject,y=dyads_neu),colour='royalblue',size=2)+
  geom_point(aes(x=subject,y=dyads_neu),colour='white',size=1,shape='N')+
  geom_point(aes(x=subject,y=dyads_neg),colour='red',size=2)+
  geom_point(aes(x=subject,y=dyads_neg),size=2,shape='-')+
  geom_point(aes(x=subject,y=dyads_pos),colour='chartreuse',size=2)+
  geom_point(aes(x=subject,y=dyads_pos),size=2,shape='+')+
  facet_wrap(~network,nrow=3,scales='free')+
  xlab('Respondent')+ylab('Sender-target dyads reported')+
  grid.background
dev.off()

########################################################################################################################

# Need to exclude networks F106b-c-d (gossip reported by less than 1/2 of the sample in these three units)
excl <- c('F106b','F106c','F106d')
attributes <- attributes[attributes$group %!in% excl,]
missing_respondents <- missing_respondents[missing_respondents %in% attributes$responder]
gossip <- gossip[gossip$group %!in% excl,]
gossip_cube <- gossip_cube[organisation_ID %!in% excl]
gos <- gos[organisation_ID %!in% excl]
gos_pos <- gos_pos[organisation_ID %!in% excl]
gos_neg <- gos_neg[organisation_ID %!in% excl]
gossip_sum <- gossip_sum[rownames(gossip_sum) %!in% excl,]
networks_mtx <- networks_mtx[organisation_ID %!in% excl]
org_subject <- org_subject[organisation_ID %!in% excl]
organisation_ID <- organisation_ID[organisation_ID %!in% excl]

########################################################################################################################

# INCONSISTENCIES IN GOSSIP?

# Senders conveying both positive and negative sentiments about the same target?
gos_pos <- as.data.frame(gossip[gossip$info == 1,])
gos_neg <- as.data.frame(gossip[gossip$info == -1,])

for(i in 1:nrow(gos_pos)){
  gos_pos$st_dyad[i] <- paste(gos_pos$sender[i],gos_pos$target[i],sep=';') 
}
for(i in 1:nrow(gos_neg)){
  gos_neg$st_dyad[i] <- paste(gos_neg$sender[i],gos_neg$target[i],sep=';') 
}

intersect_posneg <- intersect(gos_pos$st_dyad,gos_neg$st_dyad) # 48 inconsistencies
inconst_gos <- rbind(gos_pos[gos_pos$st_dyad %in% intersect_posneg,],gos_neg[gos_neg$st_dyad %in% intersect_posneg,])
inconst_gos$neg <-  1*inconst_gos$info == -1 # negative

# Inconsistencies for a total of 130 cases regarding...
length(unique(inconst_gos$sender)) # 29 unique senders
length(unique(inconst_gos$target)) # 36 unique targets

########################################################################################################################

# Removal of unnecessary objects
rm(rec_sen);rm(rec_sen_pos);rm(rec_sen_neg);rm(rec_sen_neu);rm(desc_receiver);rm(gos);rm(excl);rm(i);rm(j);rm(k);rm(x)
rm(grid.background);rm(gos_pos);rm(gos_neg);rm(gos_neu);rm(intersect_posneg);rm(inconst_gos)

# Save image
save.image('tidieddata.RData')

########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## Data tidying (1)
## R script written by Jose Luis Estevez (University of Linkoping)
## Date: December 22nd, 2020
########################################################################################################################

# R PACKAGES REQUIRED
library(ggplot2)

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

gos <- gos_pos <- gos_neg <- gossip_cube

for(x in seq_along(gossip_cube)){
  for(i in 1:nrow(gossip_cube[[x]])){
    for(j in 1:nrow(gossip_cube[[x]])){
      for(k in 1:nrow(gossip_cube[[x]])){
        if(i != j & i != k & j != k){
          # the gossip cube (no matter if positive or negative)
          if(!is.na(gossip_cube[[x]][i,j,k]) & gossip_cube[[x]][i,j,k] %in% c(-1,1)){
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
        }
      }
    }
  }
}

# How many triad in the data set? How many participants reported gossip?
rec_sen <- gos
rec_sen_pos <- gos_pos
rec_sen_neg <- gos_neg

for(x in seq_along(rec_sen)){
  rec_sen[[x]] <- apply(rec_sen[[x]],c(2,1),max,na.rm=TRUE) # receiver-sender matrix
  rec_sen_pos[[x]] <- apply(rec_sen_pos[[x]],c(2,1),max,na.rm=TRUE)
  rec_sen_neg[[x]] <- apply(rec_sen_neg[[x]],c(2,1),max,na.rm=TRUE)
  
  diag(rec_sen[[x]]) <- 0
  diag(rec_sen_pos[[x]]) <- 0
  diag(rec_sen_neg[[x]]) <- 0
  
  rec_sen[[x]] <- rowSums(rec_sen[[x]])
  rec_sen_pos[[x]] <- rowSums(rec_sen_pos[[x]])
  rec_sen_neg[[x]] <- rowSums(rec_sen_neg[[x]])
}

gossip_sum <- as.data.frame(matrix(NA,nrow=length(rec_sen),ncol=8,
                                   dimnames=list(names(rec_sen),c('triads','triads_pos','triads_neg','triads_total',
                                                                  'reporters','reporters_pos','reporters_neg',
                                                                  'reporters_total'))))

for(x in seq_along(rec_sen)){
  # summary of triads: with gossip, positive or negative, and total possible
  gossip_sum[x,'triads'] <- sum(gos[[x]],na.rm=TRUE)
  gossip_sum[x,'triads_pos'] <- sum(gos_pos[[x]],na.rm=TRUE)
  gossip_sum[x,'triads_neg'] <- sum(gos_neg[[x]],na.rm=TRUE)
  gossip_sum[x,'triads_total'] <- nrow(gos[[x]]) * (nrow(gos[[x]])-1) * (nrow(gos[[x]])-2)
  # number of respondents reporting those triads
  gossip_sum[x,'reporters'] <- sum(rec_sen[[x]] != 0)
  gossip_sum[x,'reporters_pos'] <- sum(rec_sen_pos[[x]] != 0)
  gossip_sum[x,'reporters_neg'] <- sum(rec_sen_neg[[x]] != 0)
  gossip_sum[x,'reporters_total'] <- length(rec_sen[[x]])
}

gossip_sum

# Need to exclude networks F106b-c-d (gossip reported by less than 1/3 of the sample in these three units)
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

# RECIPROCITY OF THE GOSSIP NETWORKS (TARGET-SPECIFIC)

desc_target <- vector('list',length=length(gossip_cube))
for(i in seq_along(gossip_cube)){
  desc_target[[i]] <- as.data.frame(matrix(NA,nrow=nrow(gossip_cube[[i]]),ncol=7))
  rownames(desc_target[[i]]) <- rownames(gossip_cube[[i]])
  colnames(desc_target[[i]]) <- c('network','dyads','gosip','dyads_pos','gosip_pos','dyads_neg','gosip_neg')
  desc_target[[i]]$network <- names(gossip_cube)[i]
}

for(x in seq_along(desc_target)){
  for(i in 1:nrow(desc_target[[x]])){
    desc_target[[x]]$dyads[i] <- sum(gos[[x]][,,i],na.rm=TRUE)
    desc_target[[x]]$gosip[i] <- sna::grecip(gos[[x]][,,i],measure='edgewise')
    desc_target[[x]]$dyads_pos[i] <- sum(gos_pos[[x]][,,i],na.rm=TRUE)
    desc_target[[x]]$gosip_pos[i] <- sna::grecip(gos_pos[[x]][,,i],measure='edgewise')
    desc_target[[x]]$dyads_neg[i] <- sum(gos_neg[[x]][,,i],na.rm=TRUE)
    desc_target[[x]]$gosip_neg[i] <- sna::grecip(gos_neg[[x]][,,i],measure='edgewise') 
  }
}

desc_target <- do.call('rbind',desc_target)

# Visualisation
grid.background <- theme_bw()+
  theme(plot.background=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+
  theme(axis.line=element_line(color='black'))+
  theme(strip.text.x=element_text(colour='white',face='bold'))+
  theme(strip.background=element_rect(fill='black'))

jpeg(filename='Target-specific gossip.jpeg',width=9,height=6,units='in',res=1000)
ggplot(data=desc_target)+
  geom_point(aes(x=rownames(desc_target),y=dyads_pos,colour=gosip_pos),size=5,shape='+')+
  geom_point(aes(x=rownames(desc_target),y=dyads_neg,colour=gosip_neg),size=5,shape='-')+
  facet_wrap(~network,nrow=3,scales='free')+
  scale_colour_gradient(name='Mutual',low='cornflowerblue',high='firebrick3',n.break=10)+
  ggtitle('Target-specific gossip by sentiment expressed')+xlab('Target')+ylab('Sender-goseiver dyads')+
  grid.background+theme(axis.text.x=element_blank())
dev.off()

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
rm(rec_sen);rm(rec_sen_pos);rm(rec_sen_neg);rm(desc_target);rm(gos);rm(grid.background);rm(excl)
rm(gos_pos);rm(gos_neg);rm(intersect_posneg);rm(inconst_gos);rm(i);rm(j);rm(k);rm(x)

# Save image
save.image('tidieddata.RData')

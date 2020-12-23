########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## Creation of composite networks (2)
## R script written by Jose Luis Estevez (Linkoping University)
## Date: December 23rd 2020
########################################################################################################################

# R PACKAGES REQUIRED
library(lpSolve);library(irr);library(ape);library(lattice);library(ggplot2);library(sna);library(viridis)

# DATA LOADING AND DATA TIDYING
load('tidieddata.RData')

########################################################################################################################

# COMPOSITVE NETWORK CREATION (FOLLOWING VOROS & SNIJDERS, 2017)

# 1) Exclusion of networks concerning knowledge about others' earnings (networks 15:23)
for(i in seq_along(networks_mtx)){
  networks_mtx[[i]] <- networks_mtx[[i]][c(1:14,24:43)]
}

# 2) When the matrices capture grades of the same dimension, binarization
for(i in seq_along(networks_mtx)){
  # regular_personal_conversation = 1 if at least several times per week, otherwise 0
  networks_mtx[[i]]$regular_personal_conversation <- networks_mtx[[i]]$personal_conversation_eight_a_week + 
    networks_mtx[[i]]$personal_conversation_several_times_a_week
  # regular_work_conversation = 1 if at least several times per week, otherwise 0
  networks_mtx[[i]]$regular_work_conversation <- networks_mtx[[i]]$work_conversation_eight_a_week + 
    networks_mtx[[i]]$work_conversation_several_times_a_week
}

for(i in seq_along(networks_mtx)){
  networks_mtx[[i]] <- networks_mtx[[i]][c(1:4,10:27,33:36)]
}

# 3) Matrix overlap (Jaccard indices) between each pair of matrices
Jaccard <- function(matrix1,matrix2){
  shared_ties <- matrix1*matrix2
  diff_ties <- 1*((matrix1+matrix2)==1)
  denominator <- sum(shared_ties,na.rm=TRUE)+sum(diff_ties,na.rm=TRUE)
  outcome <- ifelse(denominator==0,0,sum(shared_ties,na.rm=TRUE)/denominator)
  return(outcome)
}

mtx_overlap <- networks_mtx # Object where to allocate the matrix overlap
for(i in seq_along(mtx_overlap)){
  mtx_overlap[[i]] <- array(NA,dim=c(rep(length(mtx_overlap[[i]]),2)),
                            dimnames=list(names(mtx_overlap[[i]]),names(mtx_overlap[[i]])))
}

for(mtx in seq_along(mtx_overlap)){
  for(i in 1:nrow(mtx_overlap[[mtx]])){
    for(j in 1:ncol(mtx_overlap[[mtx]])){
      mtx_overlap[[mtx]][i,j] <- Jaccard(networks_mtx[[mtx]][[i]],networks_mtx[[mtx]][[j]])
    }
  }
}

# 4) Consistency of similarities across groups: Kendall's W
matrix_selection <- function(list_mtx_overlap){
  items_to_cluster <- rownames(list_mtx_overlap[[1]])
  clust_ranks <- vector('list',length=length(items_to_cluster))
  names(clust_ranks) <- items_to_cluster
    
  for(i in seq_along(clust_ranks)){
    clust_ranks[[i]] <- matrix(NA,nrow=length(items_to_cluster),ncol=length(organisation_ID),
                               dimnames=list(names(clust_ranks),organisation_ID))
  }
    
  for(i in seq_along(list_mtx_overlap)){
    for(j in 1:nrow(list_mtx_overlap[[i]])){
      clust_ranks[[j]][,i] <- rank(list_mtx_overlap[[i]][j,]) # ranks by row of the matrix
    }
  }
    
  for(i in seq_along(clust_ranks)){
    clust_ranks[[i]] <- clust_ranks[[i]][-i,] # remove the row that is the item itself
    clust_ranks[[i]] <- kendall(clust_ranks[[i]],TRUE)$value # Kendall W for concordance
  }
    
  min_k <- min(unlist(clust_ranks))
  item_to_exclude <- which.min(unlist(clust_ranks))
  item_name <- names(clust_ranks)[which.min(unlist(clust_ranks))]
    
  print(paste('Exclude item ',item_to_exclude,' (',item_name,'). Kendall W = ', round(min_k,3),sep=''))
}

# Application of the function
matrix_selection(mtx_overlap)
for(i in seq_along(mtx_overlap)){
  mtx_overlap[[i]] <- mtx_overlap[[i]][-24,-24]
}

matrix_selection(mtx_overlap)
for(i in seq_along(mtx_overlap)){
  mtx_overlap[[i]] <- mtx_overlap[[i]][-8,-8]
}

matrix_selection(mtx_overlap)

# 5) Hierarchical clustering of the mean Jaccard matrix
# Mean Jaccard matrix
mean_jaccard <- array(NA,dim=c(rep(nrow(mtx_overlap[[1]]),2),length(mtx_overlap)),
                      dimnames=list(rownames(mtx_overlap[[1]]),rownames(mtx_overlap[[1]]),names(mtx_overlap)))

for(i in seq_along(mtx_overlap)){
  mean_jaccard[,,i] <- mtx_overlap[[i]] # Allocation of the distance matrices into a 3D array
}
mean_jaccard <- apply(mean_jaccard,c(1,2),mean) # Mean values across all the groups

rownames(mean_jaccard) <- c('Deserve wage increase','Cooperate well','Cooperate with job duties',
                            'Would not cooperate','Deserve wage cut','I listen to them','Are a friend',
                            'Do their job well','Other colleagues listen to them','Other colleagues ask them for help',
                            'Other colleagues despise them','Other colleagues appreciate them','Are neutral to me',
                            'Are trustworthy','I appreciate them','Share negative information about me',
                            'Are not a friend','Does not belong to the team','Are not suitable for the job',
                            'Are popular','I turn for help','Receive personal information',
                            'Regular personal conversations','Regular work conversations')
colnames(mean_jaccard) <- rownames(mean_jaccard)

# Hierarchical clustering with Ward's method
cluster_model <- hclust(as.dist(1-mean_jaccard),method='ward.D')
ape_model <- as.phylo(cluster_model)
(cutrees <- cutree(cluster_model,k=4))

# Visualisation of the tree
colours <- c('darkgreen','red','royalblue','darkgrey')
jpeg(filename='Cluster tree.jpeg',width=9,height=6,units='in',res=1000)
plot(ape_model,cex=1,label.offset=.1,tip.color=colours[cutrees])
dev.off()

# 6) Check the fit with a heatmap (using the order of the hierarchical clustering above)
clust_order <- c('Receive personal information','Are popular','I turn for help',
                 'Other colleagues ask them for help','I appreciate them','Other colleagues appreciate them',
                 'Other colleagues listen to them','I listen to them','Are trustworthy','Do their job well',
                 'Are a friend','Cooperate well','Cooperate with job duties','Deserve wage increase',
                 'Regular work conversations','Regular personal conversations','Share negative information about me',
                 'Deserve wage cut','Does not belong to the team','Other colleagues despise them',
                 'Are not suitable for the job','Are not a friend','Would not cooperate','Are neutral to me')
mean_jaccard2 <- matrix(NA,nrow=nrow(mean_jaccard),ncol=ncol(mean_jaccard),dimnames=list(clust_order,clust_order))

for(i in rownames(mean_jaccard2)){
  for(j in colnames(mean_jaccard2)){
    mean_jaccard2[i,j] <- mean_jaccard[i,j]
  }
}

jpeg(filename='Heatmap.jpeg',width=8,height=8,units='in',res=1000)
levelplot(mean_jaccard2,cuts=7,at=seq(0,.7,.1),col.regions=heat.colors(100),
          xlab='',ylab='',scales=list(x=list(rot=90)),contour=FALSE,
          panel = function(...){
            panel.levelplot(...)
            panel.abline(h = 8.5)
            panel.abline(v = 8.5)
            panel.abline(h = 16.5)
            panel.abline(v = 16.5)
            panel.abline(h = 23.5)
            panel.abline(v = 23.5)
          })
dev.off()

# 7) Composite networks
# Three latent dimensions: affection, respect, and negative
affective <- c('wage_increasing','cooperate_well','cooperate_job_duties','friend','does_job_well','trustworthy',
               'regular_personal_conversation','regular_work_conversation')
respect <- c('popular','received_personal_info','listen_to_her','colleagues_listen_to_her','colleagues_ask_for_her_help',
             'colleagues_appreciate','appreciation','turn_for_her_help')
negative <- c('would_not_cooperate','wage_reduction','colleagues_despise','shares_negative_info_about_me','not_friend',
              'not_suitable_for_job','belong_to_team')

comp_networks <- vector('list',length=length(networks_mtx))
names(comp_networks) <- names(networks_mtx)

for(i in seq_along(comp_networks)){
  comp_networks[[i]] <- vector('list',length=3)
  names(comp_networks[[i]]) <- c('affective','respect','negative')
  for(j in seq_along(comp_networks[[i]])){
    if(j == 1){
      comp_networks[[i]][[j]] <- array(NA,dim=c(rep(nrow(networks_mtx[[i]][[1]]),2),8),
                                       dimnames=list(rownames(networks_mtx[[i]][[1]]),rownames(networks_mtx[[i]][[1]]),
                                                     affective))   
    }else if(j == 2){
      comp_networks[[i]][[j]] <- array(NA,dim=c(rep(nrow(networks_mtx[[i]][[1]]),2),8),
                                       dimnames=list(rownames(networks_mtx[[i]][[1]]),rownames(networks_mtx[[i]][[1]]),
                                                     respect)) 
    }else{
      comp_networks[[i]][[j]] <- array(NA,dim=c(rep(nrow(networks_mtx[[i]][[1]]),2),7),
                                       dimnames=list(rownames(networks_mtx[[i]][[1]]),rownames(networks_mtx[[i]][[1]]),
                                                     negative)) 
    }
  }
}

for(i in seq_along(comp_networks)){
  for(j in seq_along(comp_networks[[i]])){
    for(k in 1:dim(comp_networks[[i]][[j]])[3]){
      items <- colnames(comp_networks[[i]][[j]][1,,])
      comp_networks[[i]][[j]][,,k] <- networks_mtx[[i]][items][[k]]
    }
    comp_networks[[i]][[j]] <- apply(comp_networks[[i]][[j]],c(1,2),sum)
  }
}

# 8) Threshold selection
# Out-degree by cutoff point
cuts <- vector('list',length=8)
for(i in seq_along(cuts)){
  cuts[[i]] <- comp_networks
}

for(i in seq_along(cuts)){
  for(j in seq_along(cuts[[i]])){
    for(k in seq_along(cuts[[i]][[j]])){
      cuts[[i]][[j]][[k]] <- mean(sna::degree(1*(cuts[[i]][[j]][[k]] >= i),cmode='outdegree'))
    }
    cuts[[i]][[j]] <- unlist(cuts[[i]][[j]])
  }
  cuts[[i]] <- do.call('rbind',cuts[[i]])
}

degree_sum <- as.data.frame(do.call('rbind',cuts))
degree_sum$Unit <- rep(organisation_ID,8)
degree_sum$cutoff <- rep(1:8,each=6)
degree_sum <- tidyr::gather(degree_sum,key="tie",value='degree',c('affective','respect','negative'))

# Visualisation
grid.background <- theme_bw()+
  theme(plot.background=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+
  theme(axis.line=element_line(color='black'))+
  theme(strip.text.x=element_text(colour='white',face='bold'))+
  theme(strip.background=element_rect(fill='black'))

jpeg(filename='Composite network thresholds.jpeg',width=9,height=6,units='in',res=1000)
ggplot(data=degree_sum)+
  geom_line(aes(x=cutoff,y=degree,group=Unit,colour=Unit),linetype='solid',size=1.5)+
  geom_point(aes(x=cutoff,y=degree,),colour='black',size=4)+
  geom_point(aes(x=cutoff,y=degree,colour=Unit),size=3)+
  geom_vline(data=data.frame(xint=c(7,2,7),tie=c('affective','negative','respect')),
             aes(xintercept=xint),linetype='solid',colour='red',size=3,alpha=.5)+ 
  facet_wrap(~tie,nrow=1,scales='free_x')+
  xlab('Number of network items')+ylab('Average out-degree in the composite network')+
  scale_x_continuous(breaks=1:8)+
  scale_colour_manual(values=viridis(6))+
  grid.background
dev.off()

# 9) Final ties
networks_mtx <- comp_networks 

for(i in seq_along(networks_mtx)){
  for(j in seq_along(networks_mtx[[i]])){
    if(j %in% 1:2){
      networks_mtx[[i]][[j]] <- 1*(networks_mtx[[i]][[j]] >= 7)
    }else{
      networks_mtx[[i]][[j]] <- 1*(networks_mtx[[i]][[j]] >= 2)
    }
  }
}

# Inspection the overlap between these three type of ties
mtx_overlap <- networks_mtx 
for(i in seq_along(mtx_overlap)){
  mtx_overlap[[i]] <- array(NA,dim=c(rep(length(mtx_overlap[[i]]),2)),
                            dimnames=list(names(mtx_overlap[[i]]),names(mtx_overlap[[i]])))
}

for(mtx in seq_along(mtx_overlap)){
  for(i in 1:nrow(mtx_overlap[[mtx]])){
    for(j in 1:ncol(mtx_overlap[[mtx]])){
      mtx_overlap[[mtx]][i,j] <- Jaccard(networks_mtx[[mtx]][[i]],networks_mtx[[mtx]][[j]])
    }
  }
  diag(mtx_overlap[[mtx]]) <- NA
}

mtx_overlap

########################################################################################################################

# Removal of unnecessary objects
rm(cluster_model);rm(ape_model);rm(mean_jaccard);rm(mean_jaccard2);rm(networks_available);rm(clust_order);rm(Jaccard)
rm(matrix_selection);rm(mtx);rm(i);rm(j);rm(k);rm(colours);rm(cutrees);rm(cuts);rm(affective);rm(respect);rm(negative)
rm(items);rm(degree_sum);rm(grid.background);rm(mtx_overlap)

# Save image
save.image('tidieddata2.RData')

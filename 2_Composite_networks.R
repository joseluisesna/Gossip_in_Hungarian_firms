########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## Creation of composite networks (2)
## R script written by Jose Luis Estevez (Linkoping University)
## Date: March 23rd 2021
########################################################################################################################

# R PACKAGES REQUIRED
library(lpSolve);library(irr);library(ape);library(lattice);library(ggplot2);library(sna);library(viridis);library(igraph)

# DATA LOADING AND DATA TIDYING
rm(list=ls())
load('tidieddata.RData')

########################################################################################################################

# COMPOSITVE NETWORK CREATION (FOLLOWING VOROS & SNIJDERS, 2017)

# 1) Exclusion of networks concerning knowledge about others' earnings (networks 15:23) 
for(i in seq_along(networks_mtx)){
  networks_mtx[[i]] <- networks_mtx[[i]][c(1:14,24:42)]
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
  networks_mtx[[i]] <- networks_mtx[[i]][c(1:4,10:27,33:35)]
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
  mtx_overlap[[i]] <- mtx_overlap[[i]][-23,-23]
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
                            'Are popular','I turn for help','Regular personal conversations','Regular work conversations')
colnames(mean_jaccard) <- rownames(mean_jaccard)

# Hierarchical clustering with Ward's method
cluster_model <- hclust(as.dist(1-mean_jaccard),method='ward.D')
ape_model <- as.phylo(cluster_model)
(cutrees <- cutree(cluster_model,k=2))

# Visualisation of the tree
colours <- c('darkgreen','red')
jpeg(filename='Cluster tree.jpeg',width=9,height=6,units='in',res=1000)
plot(ape_model,cex=1,label.offset=.1,tip.color=colours[cutrees])
dev.off()

# 6) Check the fit with a heatmap (using the order of the hierarchical clustering above)
clust_order <- c('I turn for help','Other colleagues ask them for help','I appreciate them',
                 'Other colleagues appreciate them','Other colleagues listen to them','I listen to them','Are popular',
                 'Are trustworthy','Do their job well','Are a friend','Cooperate well','Cooperate with job duties',
                 'Deserve wage increase','Regular work conversations','Regular personal conversations',
                 'Share negative information about me','Deserve wage cut','Does not belong to the team',
                 'Other colleagues despise them','Are not suitable for the job','Are not a friend','Would not cooperate',
                 'Are neutral to me')
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
            panel.abline(h = 15.5)
            panel.abline(v = 15.5)
          })
dev.off()

# 7) Composite networks
# Two latent dimensions: positive and negative
positive <- c('wage_increasing','cooperate_well','cooperate_job_duties','friend','does_job_well','trustworthy',
               'regular_personal_conversation','regular_work_conversation','popular','listen_to_her',
              'colleagues_listen_to_her','colleagues_ask_for_her_help','colleagues_appreciate','appreciation','turn_for_her_help')
negative <- c('would_not_cooperate','wage_reduction','colleagues_despise','shares_negative_info_about_me','not_friend',
              'not_suitable_for_job','belong_to_team','is_neutral')

comp_networks <- vector('list',length=length(networks_mtx))
names(comp_networks) <- names(networks_mtx)

for(i in seq_along(comp_networks)){
  comp_networks[[i]] <- vector('list',length=2)
  names(comp_networks[[i]]) <- c('positive','negative')
  for(j in seq_along(comp_networks[[i]])){
    if(j == 1){
      comp_networks[[i]][[j]] <- array(NA,dim=c(rep(nrow(networks_mtx[[i]][[1]]),2),15),
                                       dimnames=list(rownames(networks_mtx[[i]][[1]]),rownames(networks_mtx[[i]][[1]]),
                                                     positive))  
    }else{
      comp_networks[[i]][[j]] <- array(NA,dim=c(rep(nrow(networks_mtx[[i]][[1]]),2),8),
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
cuts <- vector('list',length=max(c(length(positive),length(negative))))
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
degree_sum$Unit <- rep(organisation_ID,length(cuts))
degree_sum$cutoff <- rep(1:length(cuts),each=6)
degree_sum <- tidyr::gather(degree_sum,key="tie",value='degree',c('positive','negative'))
degree_sum$tie <- factor(degree_sum$tie,levels=c('positive','negative'),labels=c('Positive','Negative'))

# Visualisation
grid.background <- theme_bw()+
  theme(plot.background=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+
  theme(axis.line=element_line(color='black'))+
  theme(strip.text.x=element_text(colour='white',face='bold'))+
  theme(strip.background=element_rect(fill='black'))

degree_sum$Unit <- ifelse(degree_sum$Unit == 'A104','Unit A',
                          ifelse(degree_sum == 'F101','Unit B',
                                 ifelse(degree_sum == 'F103','Unit C',
                                        ifelse(degree_sum == 'F105','Unit D',
                                               ifelse(degree_sum == 'F106a','Unit E','Unit F')))))

jpeg(filename='Composite network thresholds.jpeg',width=9,height=6,units='in',res=1000)
ggplot(data=degree_sum)+
  geom_line(aes(x=cutoff,y=degree,group=Unit,colour=Unit),linetype='solid',size=1.5)+
  geom_point(aes(x=cutoff,y=degree,),colour='black',size=4)+
  geom_point(aes(x=cutoff,y=degree,colour=Unit),size=3)+
  geom_vline(data=data.frame(xint=c(NA,3),tie=factor(c('Positive','Negative'),levels=c('Positive','Negative'))),
             aes(xintercept=xint),linetype='solid',colour='red',size=6,alpha=.33)+
  geom_vline(data=data.frame(xint=c(10,NA),tie=factor(c('Positive','Negative'),levels=c('Positive','Negative'))),
             aes(xintercept=xint),linetype='solid',colour='red',size=6,alpha=.33)+
  facet_wrap(~tie,nrow=1)+
  xlab('Number of network items')+ylab('Average out-degree in the composite network')+
  scale_x_continuous(breaks=1:length(cuts))+
  scale_colour_manual(values=viridis(6))+
  grid.background
dev.off()

# 9) Final ties
networks_mtx <- comp_networks 

for(i in seq_along(networks_mtx)){
  for(j in seq_along(networks_mtx[[i]])){
    if(j == 1){
      networks_mtx[[i]][[j]] <- 1*(networks_mtx[[i]][[j]] >= 10) # positive (10/15 positive ties needed)
    }else{
      networks_mtx[[i]][[j]] <- 1*(networks_mtx[[i]][[j]] >= 3) # negative (3/8 negative tie needed)
    }
  }
}

# 10) Massive out-degree? (data correction)
# Individuals with massive out-degree in the positive network (more than 10 ties sent)
which(rowSums(networks_mtx[[1]]$positive,na.rm=TRUE)>=10)
which(rowSums(networks_mtx[[2]]$positive,na.rm=TRUE)>=10) -> massive_F101 # 1 case
which(rowSums(networks_mtx[[3]]$positive,na.rm=TRUE)>=10) -> massive_F103 # 7 cases
which(rowSums(networks_mtx[[4]]$positive,na.rm=TRUE)>=10)
which(rowSums(networks_mtx[[5]]$positive,na.rm=TRUE)>=10)
which(rowSums(networks_mtx[[6]]$positive,na.rm=TRUE)>=10)

# Use only reciprocated ties instead
rec_F101 <- 1*(networks_mtx[[2]]$positive[massive_F101,] + t(networks_mtx[[2]]$positive[,massive_F101]) == 2)
rec_F103 <- 1*(networks_mtx[[3]]$positive[massive_F103,] + t(networks_mtx[[3]]$positive[,massive_F103]) == 2)
networks_mtx[[2]]$positive[massive_F101,] <- rec_F101
networks_mtx[[3]]$positive[massive_F103,] <- rec_F103

# Inspection the overlap between positive and negative ties (Jaccard)
mtx_overlap <- networks_mtx
for(i in seq_along(mtx_overlap)){
  mtx_overlap[[i]] <- Jaccard(networks_mtx[[i]]$positive,networks_mtx[[i]]$negative)
}
unlist(mtx_overlap)

# If overlap between positive and negative, no tie exists
for(i in seq_along(networks_mtx)){
  networks_mtx[[i]]$positive[!is.na(networks_mtx[[i]]$positive + networks_mtx[[i]]$negative) & 
                          networks_mtx[[i]]$positive + networks_mtx[[i]]$negative == 2] <- 0
  networks_mtx[[i]]$negative[!is.na(networks_mtx[[i]]$positive + networks_mtx[[i]]$negative) & 
                          networks_mtx[[i]]$positive + networks_mtx[[i]]$negative == 2] <- 0
}

########################################################################################################################

# Visualisation of the networks
ntw_plot <- networks_mtx

for(i in seq_along(ntw_plot)){
  # Show both positive and negative ties
  ntw_plot[[i]]$vis <- graph_from_adjacency_matrix(ntw_plot[[i]]$positive - ntw_plot[[i]]$negative,
                                                   mode='directed',diag=FALSE,weighted=TRUE)
  # Layout based only on positive ties
  set.seed(0708)
  ntw_plot[[i]]$layout <- layout_with_fr(graph_from_adjacency_matrix(ntw_plot[[i]]$positive,mode='directed'))
  # Customisation of nodes and ties
  V(ntw_plot[[i]]$vis)$color <- ifelse(attributes[attributes$group == organisation_ID[[i]],]$woman == 1,'magenta','skyblue')
  V(ntw_plot[[i]]$vis)$shape <- ifelse(attributes[attributes$group == organisation_ID[[i]],]$hr_leader == 1,'square','circle')
  V(ntw_plot[[i]]$vis)$size <- ifelse(attributes[attributes$group == organisation_ID[[i]],]$hr_leader == 1,10,7)
  
  E(ntw_plot[[i]]$vis)$color <- as.character(factor(E(ntw_plot[[i]]$vis)$weight,levels=c(-1,1),labels=c('red','darkgreen')))
  E(ntw_plot[[i]]$vis)$width <- 1
  # Visualisation
  jpeg(filename=paste('Unit',i,'.jpeg',sep=''),width=4,height=4,units='in',res=1000)
  plot(ntw_plot[[i]]$vis,edge.arrow.size=.2,vertex.label=NA,layout=ntw_plot[[i]]$layout)
  dev.off()
}

########################################################################################################################

# Removal of unnecessary objects
rm(cluster_model);rm(ape_model);rm(mean_jaccard);rm(mean_jaccard2);rm(networks_available);rm(clust_order)
rm(matrix_selection);rm(mtx);rm(i);rm(j);rm(k);rm(colours);rm(cutrees);rm(cuts);rm(positive);rm(negative)
rm(items);rm(degree_sum);rm(grid.background);rm(mtx_overlap);rm(ntw_plot);rm(comp_networks);
rm(rec_F101);rm(rec_F103);rm(massive_F101);rm(massive_F103)

# Save image
save.image('tidieddata2.RData')
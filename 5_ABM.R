########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## ABM
## R script written by Jose Luis Estevez (Linkoping University)
## Date: February 24th 2020
########################################################################################################################

rm(list=ls())

# R PACKAGES REQUIRED
library(ggplot2);library(igraph)

########################################################################################################################

# gossip-based network evolution algorithm
ntw_evolution <- function(matrix,
                          reps=5000, # iterations
                          pos_gossip=0.5, # proportion of positive gossip
                          receiver_selectivity=.1, # how selective actors are when choosing receivers [0: very selective; 1: non-selective]
                          target_selectivity=.1, # how selective actors are when choosing targets [0: very selective; 1: non-selective]
                          force_sr=.05, # effect on the sender-receiver tie
                          force_rt=.05){ # effect on the receiver-target tie
  
  # creation of the output object: a list of matrices, as many as iterations has the ABM
  obj <- vector('list',length=reps)
  obj[[1]] <- matrix # the initiation matrix
  
  # the gossip procedure
  for(i in 2:reps){
    obj[[i]] <- obj[[i-1]]
    
    # type of gossip sent: positive or negative (by default 50% each)
    type <- rbinom(1,1,pos_gossip) 
    # sender selection (at random)
    sender <- as.numeric(sample(rownames(obj[[i-1]]),1)) 
    
    # negative gossip 
    if(type == 0){ 
      
      # target selection
      # gossipers can gossip about n_targets (this number ranging from 1 to all other nodes in the network (n-1))
      n_targets <- 1 + rbinom(1,size=(nrow(obj[[i-1]])-1),prob=target_selectivity) 
      # choose one of the n_targets who are least liked by the sender
      potential_targets <- order(obj[[i-1]][sender,-sender],decreasing=FALSE)[1:n_targets]  
      target <- as.numeric(sample(potential_targets,1))
      
      # receiver selection
      # gossipers can gossip to n_receivers (this number ranging from 1 to all other nodes in the network (n-2))
      n_receivers <- 1 + rbinom(1,size=(nrow(obj[[i-1]])-2),prob=receiver_selectivity) 
      # choose one of the n_receivers who like the target the least 
      potential_receivers <- order(obj[[i-1]][target,-c(sender,target)],decreasing=FALSE)[1:n_receivers] 
      receiver <- as.numeric(sample(potential_receivers,1))
      
      # Network evolution part
      # SR: attraction force
      sr_effect <- 1/(1+exp(-(log(obj[[i]][sender,receiver]/(1-obj[[i]][sender,receiver]))+force_sr))) 
      obj[[i]][sender,receiver] <- obj[[i]][receiver,sender] <- sr_effect
      # RT: repulsion force
      rt_effect <- 1/(1+exp(-(log(obj[[i]][receiver,target]/(1-obj[[i]][receiver,target]))-force_rt))) 
      obj[[i]][receiver,target] <- obj[[i]][target,receiver] <- rt_effect
      
      # positive gossip 
    }else{ 
      
      # target selection
      # gossipers can gossip about n_targets (this number ranging from 1 to all other nodes in the network (n-1))
      n_targets <- 1 + rbinom(1,size=(nrow(obj[[i-1]])-1),prob=target_selectivity)
      # choose one of the n_targets who are most liked by the sender 
      potential_targets <- order(obj[[i-1]][sender,-sender],decreasing=TRUE)[1:n_targets] 
      target <- as.numeric(sample(potential_targets,1))
      
      # receiver selection
      # gossipers can gossip to n_receivers (this number ranging from 1 to all other nodes in the network (n-2))
      n_receivers <- 1 + rbinom(1,size=(nrow(obj[[i-1]])-2),prob=receiver_selectivity)
      # choose one of the n_receivers who like the target the most 
      potential_receivers <- order(obj[[i-1]][target,-c(sender,target)],decreasing=TRUE)[1:n_receivers] 
      receiver <- as.numeric(sample(potential_receivers,1))
      
      # Network evolution part
      # SR: attraction force
      sr_effect <- 1/(1+exp(-(log(obj[[i]][sender,receiver]/(1-obj[[i]][sender,receiver]))+force_sr))) 
      obj[[i]][sender,receiver] <- obj[[i]][receiver,sender] <- sr_effect
      # RT: attraction force
      rt_effect <- 1/(1+exp(-(log(obj[[i]][receiver,target]/(1-obj[[i]][receiver,target]))+force_rt))) 
      obj[[i]][receiver,target] <- obj[[i]][target,receiver] <- rt_effect
    }
  }
  return(obj)
}

# Function used (visual)
force_visual <- data.frame(x=seq(0,1,by=.001))
force_visual$y <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x)))))
force_visual$attraction_.05 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))+.05)))
force_visual$repulsion_.05 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))-.05)))
force_visual$attraction_.1 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))+.1)))
force_visual$repulsion_.1 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))-.1)))
force_visual$attraction_.5 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))+.5)))
force_visual$repulsion_.5 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))-.5)))

ggplot(data=force_visual)+
  geom_line(aes(x,y),colour='black',size=1)+
  geom_line(aes(x,attraction_.05),colour='green',size=1)+
  geom_line(aes(x,attraction_.1),colour='green2',size=1)+
  geom_line(aes(x,attraction_.5),colour='green3',size=1)+
  geom_line(aes(x,repulsion_.05),colour='red',size=1)+
  geom_line(aes(x,repulsion_.1),colour='red2',size=1)+
  geom_line(aes(x,repulsion_.5),colour='red3',size=1)

########################################################################################################################

# Initiation-matrix (rnorm with mean .5 and sd .05)
initiation_mtx <- matrix(data=NA,nrow=20,ncol=20)

set.seed(290691)
for(i in 1:nrow(initiation_mtx)){
  for(j in 1:ncol(initiation_mtx)){
    if(i > j){
      initiation_mtx[i,j] <- rnorm(1,0.5,.05)
    }
  }
}

initiation_mtx <- sna::symmetrize(initiation_mtx,'lower') # symmetrisation 
rownames(initiation_mtx) <- colnames(initiation_mtx) <- 1:nrow(initiation_mtx)
hist(initiation_mtx[lower.tri(initiation_mtx,diag=FALSE)],breaks=25,xlim=c(0,1),main="Time 1",xlab='Tie strength') 

# Network evolution (how many network evolutions: now set to 10)
mtx <- vector('list',length=10)
for(i in seq_along(mtx)){
  mtx[[i]] <- initiation_mtx
}

########################################################################################################################

# ABMs
set.seed(290691)
mtx_50posgos <- lapply(mtx,ntw_evolution)

########################################################################################################################

# Results (Keep only results at rates 500, 1000, 2000, 3000, 4000, 5000)
for(i in seq_along(mtx_50posgos)){
  mtx_50posgos[[i]] <- mtx_50posgos[[i]][c(500,1000,2000,3000,4000,5000)]
  names(mtx_50posgos[[i]]) <- c(500,1000,2000,3000,4000,5000)
}

# Visualistions (histograms)
visual_ntws <- hist_ties <- vector('list')

par(mfrow=c(2,3))
for(i in seq_along(mtx_50posgos)){
  for(j in seq_along(mtx_50posgos[[i]])){
    hist(mtx_50posgos[[i]][[j]][lower.tri((mtx_50posgos[[i]][[j]]),diag=FALSE)],
         breaks=25,xlim=c(0,1),
         main=paste('Time ',names(mtx_50posgos[[i]])[[j]],sep=''),xlab='Tie strength')
  }
}

# Visualistions (networks) (only at rates 3000, 4000, and 5000)
par(mfrow=c(1,3))
for(i in seq_along(mtx_50posgos)){
  for(j in c(4,5,6)){
    x <- graph_from_adjacency_matrix(1*(mtx_50posgos[[i]][[j]] > .95),mode='undirected',diag=FALSE)
    layo <- layout_with_fr(x)
    groups <- cluster_edge_betweenness(x)
    
    
    y <- graph_from_adjacency_matrix(1*(mtx_50posgos[[i]][[j]] > .95) - 1*(mtx_50posgos[[i]][[j]] < .05),
                                     mode='undirected',weighted=TRUE,diag=FALSE)
    E(y)$color <- as.character(factor(E(y)$weight,levels=c(-1,1),labels=c('red','navyblue')))
    plot(y,layout=layo,mark.groups=groups,main=paste('Rate ',names(mtx_50posgos[[i]])[[j]],sep=''))
  }
}
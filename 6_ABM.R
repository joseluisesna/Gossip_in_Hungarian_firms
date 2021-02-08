########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## ABM
## R script written by Jose Luis Estevez (Linkoping University)
## Date: February 3rd 2020
########################################################################################################################

rm(list=ls())

# R PACKAGES REQUIRED
library(igraph)

########################################################################################################################

# gossip-based network evolution algorithm
ntw_evolution <- function(matrix,reps=2500,force=.1,pos_gossip=.5,n_receivers=5,n_targets=5){
  obj <- vector('list',length=reps)
  obj[[1]] <- matrix
  # the gossip procedure
  for(i in 2:reps){
    obj[[i]] <- obj[[i-1]]
    # type of gossip sent: positive or negative (by default 50% each)
    type <- rbinom(1,1,pos_gossip) 
    # sender selection
    sender <- as.numeric(sample(rownames(obj[[i-1]]),1)) 
    if(type == 0){ # negative goosip 
      # target selection: the n_targets who are least liked by the sender 
      potential_targets <- order(obj[[i-1]][sender,-sender],decreasing=FALSE)[1:n_targets] 
      target <- as.numeric(sample(potential_targets,1))
      # receiver selection: the n_receivers who least like the target 
      potential_receivers <- order(obj[[i-1]][target,-c(sender,target)],decreasing=FALSE)[1:n_receivers] 
      receiver <- as.numeric(sample(potential_receivers,1))
      
      # SR: attraction force
      sr_effect <- 1/(1+exp(-(log(obj[[i]][sender,receiver]/(1-obj[[i]][sender,receiver]))+force))) 
      obj[[i]][sender,receiver] <- obj[[i]][receiver,sender] <- sr_effect
      # ST: repulsion force
      st_effect <- 1/(1+exp(-(log(obj[[i]][sender,target]/(1-obj[[i]][sender,target]))-force))) 
      obj[[i]][sender,target] <- obj[[i]][target,sender] <- st_effect
      # RT: repulsion force
      rt_effect <- 1/(1+exp(-(log(obj[[i]][receiver,target]/(1-obj[[i]][receiver,target]))-force))) 
      obj[[i]][receiver,target] <- obj[[i]][target,receiver] <- rt_effect
      
    }else{ # positive gossip 
      # target selection: the n_targets who are most liked by the sender 
      potential_targets <- order(obj[[i-1]][sender,-sender],decreasing=TRUE)[1:n_targets] 
      target <- as.numeric(sample(potential_targets,1))
      # receiver selection: the n_receivers who most like the target 
      potential_receivers <- order(obj[[i-1]][target,-c(sender,target)],decreasing=TRUE)[1:n_receivers] 
      receiver <- as.numeric(sample(potential_receivers,1))
      # SR: attraction force
      sr_effect <- 1/(1+exp(-(log(obj[[i]][sender,receiver]/(1-obj[[i]][sender,receiver]))+force))) 
      obj[[i]][sender,receiver] <- obj[[i]][receiver,sender] <- sr_effect
      # ST: attraction force
      st_effect <- 1/(1+exp(-(log(obj[[i]][sender,target]/(1-obj[[i]][sender,target]))+force))) 
      obj[[i]][sender,target] <- obj[[i]][target,sender] <- st_effect
      # RT: attraction force
      rt_effect <- 1/(1+exp(-(log(obj[[i]][receiver,target]/(1-obj[[i]][receiver,target]))+force))) 
      obj[[i]][receiver,target] <- obj[[i]][target,receiver] <- rt_effect
    }
  }
  return(obj)
}

# Function used (visual)
force_visual <- data.frame(x=seq(0,1,by=.001))
force_visual$y <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x)))))
force_visual$attraction_.1 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))+.1)))
force_visual$repulsion_.1 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))-.1)))
force_visual$attraction_.5 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))+.5)))
force_visual$repulsion_.5 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))-.5)))
force_visual$attraction_1 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))+1)))
force_visual$repulsion_1 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))-1)))
force_visual$attraction_1.5 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))+1.5)))
force_visual$repulsion_1.5 <- 1/(1+exp(-(log(force_visual$x/(1-force_visual$x))-1.5)))

ggplot(data=force_visual)+
  geom_line(aes(x,y),colour='black',size=1)+
  geom_line(aes(x,attraction_.1),colour='green4',size=1)+
  geom_line(aes(x,attraction_.5),colour='green3',size=1)+
  geom_line(aes(x,attraction_1),colour='green2',size=1)+
  geom_line(aes(x,attraction_1.5),colour='green1',size=1)+
  geom_line(aes(x,repulsion_.1),colour='red4',size=1)+
  geom_line(aes(x,repulsion_.5),colour='red3',size=1)+
  geom_line(aes(x,repulsion_1),colour='red2',size=1)+
  geom_line(aes(x,repulsion_1.5),colour='red1',size=1)

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

# Network evolution (how many network evolutions: now set to 5)
mtx <- vector('list',length=5)
for(i in seq_along(mtx)){
  mtx[[i]] <- initiation_mtx
}

########################################################################################################################

# ABMs
set.seed(290691)
mtx_50posgos <- lapply(mtx,ntw_evolution,reps=2500,force=.1,pos_gossip=0.5,n_receivers=3,n_targets=3)

########################################################################################################################

# Results (Keep only results at rates 100, 500, 1000, 1500, 2000, 2500)
for(i in seq_along(mtx_50posgos)){
  mtx_50posgos[[i]] <- mtx_50posgos[[i]][c(100,500,1000,1500,2000,2500)]
  names(mtx_50posgos[[i]]) <- c(100,500,1000,1500,2000,2500)
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

# Visualistions (networks) (only at rates 1000, 1500, 2000, 2500)
par(mfrow=c(2,2))
for(i in seq_along(mtx_50posgos)){
  for(j in 3:6){
    x <- graph_from_adjacency_matrix(1*(mtx_50posgos[[i]][[j]] > .95),mode='undirected',diag=FALSE)
    layo <- layout_with_fr(x)
    groups <- cluster_edge_betweenness(x)
    y <- graph_from_adjacency_matrix(1*(mtx_50posgos[[i]][[j]] > .95) - 1*(mtx_50posgos[[i]][[j]] < .05),
                                     mode='undirected',weighted=TRUE,diag=FALSE)
    E(y)$color <- as.character(factor(E(y)$weight,levels=c(-1,1),labels=c('red','navyblue')))
    plot(y,layout=layo,mark.groups=groups,main=paste('Rate ',names(mtx_50posgos[[i]])[[j]],sep=''))
  }
}
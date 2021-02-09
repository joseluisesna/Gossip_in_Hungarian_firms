########################################################################################################################
## GOSSIP IN HUNGARIAN FIRMS
## Descriptive analysis (3)
## R script written by Jose Luis Estevez (Linkoping University)
## Date: January 15th 2021
########################################################################################################################

# R PACKAGES REQUIRED
library(ggplot2)

# DATA LOADING AND DATA TIDYING
load('tidieddata2.RData')

########################################################################################################################

# TRIADIC DESCRIPTION 

# Creation of large data set N(N-1)(N-2)
triad_data <- vector('list',length=length(org_subject))
names(triad_data) <- names(org_subject)

# Extraction of all the triads
for(i in seq_along(triad_data)){
  triad_data[[i]] <- data.frame(matrix(NA,nrow=length(org_subject[[i]])^3,ncol=4))
  colnames(triad_data[[i]]) <- c('group','sender','receiver','target')
  triad_data[[i]]$group <- organisation_ID[[i]]
  triad_data[[i]]$sender <- rep(org_subject[[i]],each=(length(org_subject[[i]]))^2)
  triad_data[[i]]$receiver <- rep(org_subject[[i]],times=length(org_subject[[i]]),each=length(org_subject[[i]]))
  triad_data[[i]]$target <- rep(org_subject[[i]],times=(length(org_subject[[i]])^2))
  # Exclusion if sender = receiver, sender = target, or receiver = target    
  triad_data[[i]] <- triad_data[[i]][triad_data[[i]]$sender != triad_data[[i]]$receiver & 
                                       triad_data[[i]]$sender != triad_data[[i]]$target & 
                                       triad_data[[i]]$receiver != triad_data[[i]]$target,]
  triad_data[[i]] <- merge(triad_data[[i]],gossip[,-1],by=c('group','sender','receiver','target'),all.x=TRUE)
  # Allocation of missing data
  for(j in 1:nrow(triad_data[[i]])){
    if(triad_data[[i]]$receiver[j] %in% missing_respondents){
      triad_data[[i]]$missing[j] <- TRUE
    }else{
      triad_data[[i]]$missing[j] <- FALSE
    }
  }
}

# Extraction of ties from the composite networks and allocation in triad_data
for(i in seq_along(triad_data)){
  for(j in 1:nrow(triad_data[[i]])){
    # Positive ties are reciprocated
    triad_data[[i]]$sr_pos[j] <- sym_mtx[[i]]$expressive[triad_data[[i]]$sender[j],triad_data[[i]]$receiver[j]] 
    triad_data[[i]]$st_pos[j] <- sym_mtx[[i]]$expressive[triad_data[[i]]$sender[j],triad_data[[i]]$target[j]] 
    triad_data[[i]]$rt_pos[j] <- sym_mtx[[i]]$expressive[triad_data[[i]]$receiver[j],triad_data[[i]]$target[j]] 
    # Negative ties are directed
    triad_data[[i]]$sr_neg[j] <- networks_mtx[[i]]$negative[triad_data[[i]]$sender[j],triad_data[[i]]$receiver[j]] 
    triad_data[[i]]$rs_neg[j] <- networks_mtx[[i]]$negative[triad_data[[i]]$receiver[j],triad_data[[i]]$sender[j]] 
    triad_data[[i]]$st_neg[j] <- networks_mtx[[i]]$negative[triad_data[[i]]$sender[j],triad_data[[i]]$target[j]] 
    triad_data[[i]]$rt_neg[j] <- networks_mtx[[i]]$negative[triad_data[[i]]$receiver[j],triad_data[[i]]$target[j]] 
  }
}

triad_data <- do.call('rbind',triad_data) 

triad_data$sr_neg_any <- 1*((triad_data$sr_neg + triad_data$rs_neg) > 0) # the tie sender-receiver is negative in one direction

# Bivariate description (type of gossip by kind of relational triad)

triad_data$sr <- triad_data$sr_pos - triad_data$sr_neg_any
triad_data$st <- triad_data$st_pos - triad_data$st_neg
triad_data$rt <- triad_data$rt_pos - triad_data$rt_neg

biv_descrip <- data.frame(table(SR=triad_data$sr,ST=triad_data$st,RT=triad_data$rt)) 
gossip_triads <- data.frame(table(SR=triad_data$sr,ST=triad_data$st,RT=triad_data$rt,gossip=triad_data$info)) 

biv_descrip$pos_gos <- gossip_triads[gossip_triads$gossip == 1,'Freq']
biv_descrip$neg_gos <- gossip_triads[gossip_triads$gossip == -1,'Freq']
biv_descrip$neu_gos <- gossip_triads[gossip_triads$gossip == 0,'Freq']

biv_descrip$prop_pos <- round(biv_descrip$pos_gos / biv_descrip$Freq,3)*100
biv_descrip$prop_neg <- round(biv_descrip$neg_gos / biv_descrip$Freq,3)*100
biv_descrip$prop_neu <- round(biv_descrip$neu_gos / biv_descrip$Freq,3)*100

biv_descrip

# Removal of unnecessary objects
rm(gossip_triads);rm(biv_descrip);rm(i);rm(j)

########################################################################################################################

# BIPARTITE APPROACH (sender-receiver dyads vs. gossip targets)

# Find  those cases where both parties (sender and receiver) reported gossip about the same target
gossip$RST <- paste(gossip$receiver,gossip$sender,gossip$target,sep=';')
gossip$SRT <- paste(gossip$sender,gossip$receiver,gossip$target,sep=';')

oneside_reported <- gossip[gossip$RST %!in% gossip$SRT,] # only one party reported (N = 1323)
twosides_reported <- gossip[gossip$RST %in% gossip$SRT,] # both parties did (N = 238/2)

twosides_reported$info2 <- NA
for(i in 1:nrow(twosides_reported)){
  x <- twosides_reported$RST[i]
  twosides_reported$info2[i] <- twosides_reported[twosides_reported$SRT == x,]$info
}

# Split when both parties agreed or disagreed
twosides_agree <- twosides_reported[twosides_reported$info == twosides_reported$info2,] # Agree (N = 126/")
twosides_disagree <- twosides_reported[twosides_reported$info != twosides_reported$info2,] # Disagree (N = 112/")

# If one party said "neutral" and the other "positive" or "negative", keep the valence of the second
twosides_disagree$info <- twosides_disagree$info + twosides_disagree$info2
sum(twosides_disagree$info == 0) # How many disagreements left? 18/2
twosides_disagree[twosides_disagree$info == 0,]$info <- -1 # Consider negative if one positive and the other negative

gossip_dy <- rbind(oneside_reported[,c('group','receiver','sender','target','info')],
                   twosides_agree[,c('group','receiver','sender','target','info')],
                   twosides_disagree[,c('group','receiver','sender','target','info')])

table(gossip$info);table(gossip_dy$info) # Now we have 30 negative gossips and 17 positive gossips more

# Symmetrisation (if either party reported gossip)
gossip_dy2 <- gossip_dy
names(gossip_dy2) <- c('group','sender','receiver','target','info')
gossip_dy <- merge(gossip_dy,gossip_dy2,all=TRUE)

# Creation of large data (gossip-pairs and targets)
bipar_data <- vector('list',length=length(org_subject))
names(bipar_data) <- names(org_subject)

# Extraction of all the triads
for(i in seq_along(bipar_data)){
  bipar_data[[i]] <- data.frame(matrix(NA,nrow=length(org_subject[[i]])^3,ncol=4))
  colnames(bipar_data[[i]]) <- c('group','sender','receiver','target')
  bipar_data[[i]]$group <- organisation_ID[[i]]
  bipar_data[[i]]$sender <- rep(org_subject[[i]],each=(length(org_subject[[i]]))^2)
  bipar_data[[i]]$receiver <- rep(org_subject[[i]],times=length(org_subject[[i]]),each=length(org_subject[[i]]))
  bipar_data[[i]]$target <- rep(org_subject[[i]],times=(length(org_subject[[i]])^2))
  # Exclusion if sender < receiver, sender = target, or receiver = target (half the cube)    
  bipar_data[[i]] <- bipar_data[[i]][bipar_data[[i]]$sender < bipar_data[[i]]$receiver & 
                                       bipar_data[[i]]$sender != bipar_data[[i]]$target & 
                                       bipar_data[[i]]$receiver != bipar_data[[i]]$target,]
  bipar_data[[i]] <- merge(bipar_data[[i]],gossip_dy,by=c('group','sender','receiver','target'),all.x=TRUE)
  # Allocation of missing data (if either sender or receiver is missing)
  for(j in 1:nrow(bipar_data[[i]])){
    if(bipar_data[[i]]$receiver[j] %in% missing_respondents | bipar_data[[i]]$sender[j] %in% missing_respondents){
      bipar_data[[i]]$missing[j] <- TRUE
    }else{
      bipar_data[[i]]$missing[j] <- FALSE
    }
  }
}

bipar_data <- do.call('rbind',bipar_data) 

# Addition of positive and negative ties
bipar_data <- merge(bipar_data,triad_data[,c('group','receiver','sender','target','sr','st','rt')],all.x=TRUE)

# Dyads and dyad-to-target ties
bipar_data$dyads <- paste(bipar_data$sender,bipar_data$receiver,sep=';')
bipar_data$dyad_target <- NA
for(i in 1:nrow(bipar_data)){
  if(!is.na(bipar_data$st[i]) & !is.na(bipar_data$rt[i])){
    if(bipar_data$st[i] == 1 & bipar_data$rt[i] == 1){
      bipar_data$dyad_target[i] <- '+2'
    }else if((bipar_data$st[i] == 1 & bipar_data$rt[i] == 0) | (bipar_data$st[i] == 0 & bipar_data$rt[i] == 1)){
      bipar_data$dyad_target[i] <- '+1'
    }else if((bipar_data$st[i] == -1 & bipar_data$rt[i] == 1) | (bipar_data$st[i] == 1 & bipar_data$rt[i] == -1)){
      bipar_data$dyad_target[i] <- '+/-1'
    }else if((bipar_data$st[i] == -1 & bipar_data$rt[i] == 0) | (bipar_data$st[i] == 0 & bipar_data$rt[i] == -1)){
      bipar_data$dyad_target[i] <- '-1'
    }else if(bipar_data$st[i] == -1 & bipar_data$rt[i] == -1){
      bipar_data$dyad_target[i] <- '-2'
    }else{
      bipar_data$dyad_target[i] <- '0'
    }
  }
}

bipar_data$dyad_target <- factor(bipar_data$dyad_target,levels=c('0','-2','-1','+/-1','+1','+2'))

table_dyads <- data.frame(table(dyad_to_target=bipar_data$dyad_target,SR=bipar_data$sr))
dyads_gos <- data.frame(table(dyad_to_target=bipar_data$dyad_target,SR=bipar_data$sr,gossip=bipar_data$info))

table_dyads$pos_gos <- dyads_gos[dyads_gos$gossip == 1,'Freq']
table_dyads$neg_gos <- dyads_gos[dyads_gos$gossip == -1,'Freq']
table_dyads$neu_gos <- dyads_gos[dyads_gos$gossip == 0,'Freq']

table_dyads$prop_pos <- round(table_dyads$pos_gos / table_dyads$Freq,3)*100
table_dyads$prop_neg <- round(table_dyads$neg_gos / table_dyads$Freq,3)*100
table_dyads$prop_neu <- round(table_dyads$neu_gos / table_dyads$Freq,3)*100

table_dyads

# Visualisation
bipar_data$Gossip <- ifelse(bipar_data$missing==TRUE,NA,
                            ifelse((bipar_data$missing == FALSE & !is.na(bipar_data$info) & bipar_data$info == 0),'Neutral',
                                   ifelse((bipar_data$missing == FALSE & !is.na(bipar_data$info) & bipar_data$info == 1),'Positive',
                                          ifelse((bipar_data$missing == FALSE & !is.na(bipar_data$info) & bipar_data$info == -1),'Negative','No gossip'))))

bipar_data$sr <- factor(bipar_data$sr,levels=c(0,1,-1),labels=c('No tie','Positive tie','Negative tie'))

grid.background <- theme_bw()+
  theme(plot.background=element_blank(),panel.grid.minor=element_blank(),panel.border=element_blank())+
  theme(axis.line=element_line(color='black'))+
  theme(strip.text.x=element_text(colour='white',face='bold'))+
  theme(strip.background=element_rect(fill='black'))

jpeg(filename='Bivariate.jpeg',width=7,height=4,units='in',res=1000)
ggplot(data=bipar_data[bipar_data$missing==FALSE,],aes(x=dyad_target,fill=Gossip,na.rm=TRUE))+
  geom_bar(colour='black',position='fill')+
  facet_wrap(~sr,nrow=1)+
  scale_fill_manual(values=c('red','orange','grey','chartreuse'))+
  xlab('Relations of the dyad with the gossip target')+ylab('Proportion')+
  grid.background
dev.off()

# Removal of unnecessary objects
rm(oneside_reported);rm(twosides_reported);rm(twosides_agree);rm(twosides_disagree);rm(gossip_dy2);rm(dyads_gos)
rm(table_dyads);rm(i);rm(j);rm(x)

save.image('modellingdata.RData') 

########################################################################################################################

# NODE-LEVEL APPROACH (BASED ON DEGREES)

# Gossip targets (how many colleagues engaged in gossip about a certain person x)
gossip_dy <- gossip_dy[,c('receiver','target','info')] # in gossip_dy the gossip is reciprocal (s-r <-> r-s)
gossip_dy <- gossip_dy[!duplicated(gossip_dy),] # hence, exclusion of duplicated

nodes <- data.frame(table(gossip_dy$target,gossip_dy$info)) 
nodes <- tidyr::spread(nodes,key=Var2,value=Freq)
names(nodes) <- c('ID','neg_target','neu_target','pos_target')

# Degree centrality (symmetrizided networks for positive ties)
for(mtx in seq_along(sym_mtx)){
  for(item in seq_along(sym_mtx[[mtx]])){
    sym_mtx[[mtx]][[item]] <- colSums(sym_mtx[[mtx]][[item]],na.rm=TRUE)
  }
  sym_mtx[[mtx]] <- do.call('cbind',sym_mtx[[mtx]])
}
sym_mtx <- do.call('rbind',sym_mtx)

# In-degree centrality (directed networks for negative tie)
for(mtx in seq_along(networks_mtx)){
  for(item in seq_along(networks_mtx[[mtx]])){
    networks_mtx[[mtx]][[item]] <- colSums(networks_mtx[[mtx]][[item]],na.rm=TRUE)
  }
  networks_mtx[[mtx]] <- do.call('cbind',networks_mtx[[mtx]])
}
networks_mtx <- do.call('rbind',networks_mtx)

colnames(sym_mtx) <- c('pos_deg','neg_deg')
colnames(networks_mtx) <- c('pos_indeg','neg_indeg')

nodes <- do.call('cbind',list(nodes,sym_mtx,networks_mtx,attributes[,-1]))

# Visualisation
load('modellingdata.RData') 

# Colours for the nodes based on how often one is a negative target
for(i in 1:nrow(nodes)){
  gradient <- colorRampPalette(c('grey','red'))
  x <- nodes[nodes$group == nodes$group[[i]],]
  gradient <- gradient(max(x$neg_target)+1)
  nodes$color[i] <- gradient[nodes$neg_target[i]+1]
}

# Visualisation of the networks
ntw_plot <- sym_mtx
for(i in seq_along(ntw_plot)){
  # Clustering using Newman's method based on positive ties only
  ntw_plot[[i]]$group <- ntw_plot[[i]]$expressive 
  ntw_plot[[i]]$group <- graph_from_adjacency_matrix(ntw_plot[[i]]$group,mode='undirected',diag=FALSE)
  ntw_plot[[i]]$group <- igraph::cluster_edge_betweenness(ntw_plot[[i]]$group) # clustering
  # Show both expressive and negative ties
  ntw_plot[[i]]$vis <- graph_from_adjacency_matrix(sym_mtx[[i]]$expressive - networks_mtx[[i]]$negative,
                                                   mode='directed',diag=FALSE,weighted=TRUE)
  # Layout based only on expressive ties
  ntw_plot[[i]]$layout <- layout_with_kk(graph_from_adjacency_matrix(ntw_plot[[i]]$expressive,mode='directed'))
  # Customisation of nodes and ties
  V(ntw_plot[[i]]$vis)$color <- nodes[nodes$group == organisation_ID[[i]],]$color
  V(ntw_plot[[i]]$vis)$shape <- ifelse(attributes[attributes$group == organisation_ID[[i]],]$hr_leader == 1,'square','circle')
  V(ntw_plot[[i]]$vis)$size <- ifelse(attributes[attributes$group == organisation_ID[[i]],]$hr_leader == 1,10,7)
  
  E(ntw_plot[[i]]$vis)$color <- as.character(factor(E(ntw_plot[[i]]$vis)$weight,
                                                    levels=c(-1,1),labels=c('red','darkgreen')))
  E(ntw_plot[[i]]$vis)$width <- 1
  # Visualisation
  jpeg(filename=paste('Unit',i,'.jpeg',sep=''),width=4,height=4,units='in',res=1000)
  plot(ntw_plot[[i]]$vis,mark.groups=ntw_plot[[i]]$group,
       vertex.label=NA,edge.arrow.size=.1,layout=ntw_plot[[i]]$layout)
  dev.off()
}

save.image('modellingdata.RData') 
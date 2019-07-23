library(igraph)
library(GraphAlignment)
setwd("~/")

genes <- read.table("new_res_sub_s4",header=F,sep="\t");
g=as.matrix(genes)
g = gsub(" ","",g);
y <- graph_from_edgelist(g[,1:2], directed=FALSE)
E(y)$weight = as.numeric(g[,3,drop=F]);

y <- simplify(y)
bad.vs<-V(y)[degree(y) == 1] 
# remove isolated nodes
y <-delete.vertices(y, bad.vs)

is.connected(y)

diameter(y, directed = FALSE,unconnected = FALSE, weights = NULL);
get.diameter(y, directed = FALSE, unconnected = FALSE, weights = NULL);

sg1 <- spinglass.community(y, weights = E(y)$weight, update.rule="config",start.temp=1, stop.temp=.30, gamma =2)

edgelist <- tapply(seq_along(membership(sg1)), membership(sg1), function(xx) xx)
comList <- tapply(membership(sg1), membership(sg1), names)

length(comList)                         ## number of communities
comsize <- sapply(comList, length)
comsize
bigComIndex <- which(comsize == max(comsize))
bigComIndex
h <- comList[bigComIndex]
h

degree = data.frame(unlist(degree(y)));

hub_genes <- NULL;
for(i in 1:length(comList)){
  tmp = data.frame(matrix(comList[[i]],ncol=1),degree[comList[[i]],]);
  tmp_order = tmp[order(tmp$degree.comList..i.....),];
  hub = matrix(tmp_order[(nrow(tmp_order)-9):nrow(tmp_order),1],ncol=1);
  hub_genes = rbind(hub_genes,hub);
}
write.table(cc, "new_res_sub_s4_spinglass.txt",quote = F,row.names = F,col.names = F);
write.table(hub_genes, "new_res_sub_s4_hubgenes.txt",quote = F,col.names = F,row.names = F)

library(MatrixEQTL)
library(ggplot2)
library(preprocessCore);
setwd("~/Data/");

### Stage data propcessing:
miR <- as.data.frame(read.table("stage4_mirna_matrix.txt", sep="\t", header=TRUE));
rna <- as.data.frame(read.table("stage4_rna_matrix.txt", sep="\t", header=TRUE));

mi <- as.matrix(miR);
r <- as.matrix(rna);
y = gsub(" ","",r);

sampleID = which(colnames(mi) %in% colnames(y));
x = mi[,sampleID];

### Dimension of x1: stage1 757 * 32
### Dimension of x1: stage2 807 * 82
### Dimension of x1: stage3 779 * 59
### Dimension of x1: stage4 727 * 23
x1 <- NULL;
for(i in 1:nrow(x)){
  #cat(i,"\n");
  index <- as.numeric(x[i,(2:ncol(x))]);
  if(sum(index!=0)){x1 <- rbind(x1,x[i,,drop=F]);}
}

### Dimension of x2: stage1 345 * 32
### Dimension of x2: stage2 339 * 82
### Dimension of x2: stage3 337 * 59
### Dimension of x2: stage4 339 * 23
x2 <- NULL;
for(i in 1:nrow(x1)){
  #cat(i,"\n");
  index <- as.numeric(x1[i,(2:ncol(x1))]);
  if(sum(index!=0) > 0.9 * (ncol(x) -1)){ x2 <- rbind(x2,x1[i,]);}
  else{print(i)}
}

### Dimension of y1: stage1 19551 * 32
### Dimension of y1: stage2 19717 * 82
### Dimension of y1: stage3 19684 * 59
### Dimension of y1: stage4 19363 * 23
y1 <- NULL;
for(i in 1:nrow(y)){
  cat(i,"\n");
  index <- as.numeric(y[i,(2:ncol(y))]);
  if(sum(index!=0)){y1 <- rbind(y1,y[i,]);}
}
### Dimension of y2: stage1 15824 * 32
### Dimension of y2: stage2 15717 * 82
### Dimension of y2: stage3 15849 * 59
### Dimension of y2: stage4 15942 * 23
y2 <- NULL;
for(i in 1:nrow(y1)){
  cat(i,',',sep="");
  index <- as.numeric(y1[i,(2:ncol(y1))]);
  if(sum(index!=0) > 0.9 * (ncol(y) -1)){ y2 <- rbind(y2,y1[i,]);}
}

write.table(x2, file="./stage4_mirna.txt",col.names =T , sep="\t",quote = F,row.names = F);
write.table(y2, file="./stage4_rna.txt",col.names=T, sep="\t",quote = F,row.names = F);

sub_x2 = x2[,2:ncol(x2)];
class(sub_x2) <- "numeric";
sub_y2 = y2[,2:ncol(y2)];
class(sub_y2) <- "numeric";
z.mirna <- t(scale(t(sub_x2), center = TRUE, scale = TRUE));
z.rna <- t(scale(t(sub_y2), center = TRUE, scale = TRUE));

mirna_norm = normalize.quantiles(z.mirna,copy=TRUE);
rna_norm = normalize.quantiles(z.rna,copy=TRUE);
colnames(mirna_norm) = colnames(sub_x2);
colnames(rna_norm) = colnames(sub_y2);
rownames(mirna_norm) = x2[,1];
rownames(rna_norm) = y2[,1];
write.table(mirna_norm, file="./stage4_mirna_norm.txt",col.names =T , quote = F, row.names = T,sep="\t");
write.table(rna_norm, file="./stage4_rna_norm.txt",col.names=T, quote = F,row.names = T,sep="\t");


### PCA analysis on miRNA:
{
  data = as.matrix(mirna_norm);
  cvrt0 = rbind(matrix(1,1,ncol(data)));
  cvrt0 = t(qr.Q(qr(t(cvrt0))));
  
  data = data - tcrossprod(data,cvrt0) %*% cvrt0;
  
  e = eigen(crossprod(data))
  pdf("./stage4_mirna_pc.pdf")
  plot(e$values/sum(e$values), main='miRNA Scree Plot', ylab = 'Explained variaton', pch = 19, col = 'blue');
  dev.off()
  for( i in 1:3 ) {
    plot(e$vectors[i,], col = 'blue', pch = 19, main = paste0('CNV PC ',i))
  }
  sum(e$values > 3*median(e$values))
  svcvrt = t(e$vectors[,1:15]);
}

### PCA analysis on gene expression:
{
 #  library(impute)
  data = as.matrix(rna_norm);
  #data <- impute.knn(data ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)$data;
  cvrt0 = rbind(matrix(1,1,ncol(data)));
  cvrt0 = t(qr.Q(qr(t(cvrt0))));
  
  data = data - tcrossprod(data,cvrt0) %*% cvrt0;
  
  e = eigen(crossprod(data))
  pdf("./stage4_rna_PC.pdf");
  plot(e$values/sum(e$values), main='RNA Scree Plot', ylab = 'Explained variaton', pch = 19, col = 'blue');
  dev.off();
  for( i in 1:6 ) {
    plot(e$vectors[i,], col = 'blue', pch = 19, main = paste0('Ribosome PC ',i))
  }
  sum(e$values > 4*median(e$values));
  moleculartraitcvrt = t(e$vectors[,1:20]);
}

cvrtWithPCs = rbind(moleculartraitcvrt, svcvrt);
cvrtWithPCs <- cbind(1:nrow(cvrtWithPCs),cvrtWithPCs);
write.table(cvrtWithPCs,"./cvrtWithPCs_stage4",col.names=T,row.names=F,quote=T,sep='\t');


### Run MatrixeQTL for each stage:
library(MatrixEQTL);
rm(list=ls());
source('~/Dropbox (UNC Charlotte)/eQTL_data_Li/mxeqtl.R');
cvrtWithPCs <- read.table("./cvrtWithPCs_stage4",header=T);
max <- matrix(".",20,15);
for(i in 1:20){
  for(j in 1:15){
    cvrtWithPCs_new <- rbind(cvrtWithPCs[i,],cvrtWithPCs[21:(21+j-1),]);
    write.table(cvrtWithPCs_new,"./cvrtWithPCs_stage4_sub",quote=F,col.names=T,row.names=F,sep='\t');
    me<-mxeqtl('./stage4_mirna_norm.txt','mirna_posi','./stage4_rna_norm.txt','./gene.posi',covariates="./cvrtWithPCs_stage4_sub",
               cis_output_file=paste('./stage4_cisresults_(',i,",",j,')',sep=""),cis_pval=0.05,model="linear", MAF=0, cis_dist=250000000, missing="NA",
               trans_output_file=paste('./stage4_transresults_(',i,",",j,')',sep=""),trans_pval=1e-4);
    a <- me$cis$eqtls[which(me$cis$eqtls$FDR<0.05),];
    num_cis <-nrow(a);
    b <- me$trans$eqtls[which(me$trans$eqtls$FDR<0.05),];
    num_trans <- nrow(b);
    max[i,j] = num_cis + num_trans;
    cvrtWithPCs_new = NULL;
  }
}
write.table(max,"./max_stage4",quote=F,col.names = T,row.names = T,sep="\t")

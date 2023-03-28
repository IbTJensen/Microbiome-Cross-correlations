library(data.table)
OTU <- fread(file = "/home/au681900/GenomeDK/faststorage/CaseB/Microbiome_nitrate.csv")
OTU <- OTU[genotype == "Gifu" & Nitrate_supplied == "no" & compartments == "root"]
M <- as.matrix(OTU[,-(1:5)])
M <- M[,colSums(M) != 0]
template_OTU <- apply(M, 2, function(x) weighted.mean(x, ifelse(x == 0, 0.5, 1) ))
# template_OTU <- colMeans(M, na.rm = T); template_OTU <- template_OTU[template_OTU!=0]
weigted.var <- function(x, w){
  wv <- weighted.mean(x^2, w) - weighted.mean(x, w)^2
  return(wv)
}
# vars <- apply(M, 2, var, na.rm = T)
vars <- apply(M, 2, function(x) weigted.var(x, ifelse(x == 0, 0.5, 1) ))
zp <- jitter(apply(M, 2, function(x) mean(x==0))/2)
zp[zp<0] <- 0

pi <- zp
sigma2 <- log(1 + vars/template_OTU^2)
mu <- log(template_OTU) - sigma2/2

D <- data.table(mu = mu, sigma2 = sigma2, pi = pi)
D <- D[order(mu)][71:1070]
fwrite(D, "/home/au681900/GenomeDK/faststorage/Correlation paper/Simulations/Templates/OTU_template.csv")

Gene <- fread(file = "/home/au681900/GenomeDK/faststorage/CaseC/RNA-seq_rawcounts.csv")
M2 <- as.matrix(Gene[,-1])
M2 <- M2[rowSums(M2) > 100,]
# template_gene <- rowMeans(M2, na.rm = T); template_gene <- template_gene[template_gene!=0]
template_gene <- apply(M2, 1, function(x) weighted.mean(x, ifelse(x == 0, 0.5, 1) ))
# vars <- apply(M2, 1, var, na.rm = T)
vars <- apply(M2, 1, function(x) weigted.var(x, ifelse(x == 0, 0.5, 1) ))
zp <- jitter(apply(M2, 1, function(x) mean(x==0))/2)
zp[zp<0] <- 0

pi <- zp
sigma2 <- log(1 + vars/template_gene^2)
mu <- log(template_gene) - sigma2/2

D <- data.table(mu = mu, sigma2 = sigma2, pi = pi)
fwrite(D, "/home/au681900/GenomeDK/faststorage/Correlation paper/Simulations/Templates/Gene_template.csv")

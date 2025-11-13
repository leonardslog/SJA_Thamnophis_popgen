setwd("~/Documents/sji_thamnophis/")
library(terra)
library(scales)
library(adegenet)
library(plotrix)
library(poppr)

# range of sampling coordinates
# 47.932753,48.697517
# -123.179626,-122.403992

# set map plotting parameters based off above
x <- c(-123.19,-122.42)
y <- c(48.4, 48.72)

# plot North America
north_america <- vect("~/Documents/sji_thamnophis/NorthAmerica_shape/North_America.shp")
plot(north_america,xlim=x, ylim=y,col ="grey")
sja_ext <- ext(-123.25, -122.73, 48.4145278930001, 48.728694916)

sja_shape <- crop(north_america, sja_ext)





#palette colors       
piecolors.ord <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")
piecolors.sir <- c("#e6ab02","#66a61e","#e7298a","#7570b3")
piecolors.ele <- c("#d95f02","#66a61e","#e7298a")

ordinoides.col <- c("black","#d73027")
sirtalis.col <- c("#1f78b4", "#a6cee3","#b2df8a", "#33a02c")
elegans.col <- c("gold","brown")

#Read data from structure
ord.snps <- import2genind("~/Documents/sji_thamnophis/stacks/ref_sirtalis/ordinoides/ordinoides_SJI/unlinked_p1/ord_usnps_p1r80.gen")
ord.popnames <- c("Henry", "Lopez", "San Juan", "Orcas", "Shaw", "Waldron")

sir.snps <- import2genind("~/Documents/sji_thamnophis/stacks/ref_sirtalis/sirtalis/sirtalis_SJI/unlinked_p1/sir_usnps_p1r80.gen")
sir.popnames <- c("Henry", "Orcas", "Lopez", "San Juan", "Shaw", "Waldron")

ele.snps <- import2genind("~/Documents/sji_thamnophis/stacks/ref_sirtalis/elegans/elegans_SJI/unlinked_p1/ele_usnps_p1r80.gen")
ele.popnames <- c("San Juan", "Shaw","Orcas")

# sample info
inddata <- read.csv("sample_localities.csv")
indcoords <- inddata[,c(7,4,5)]

##############
# ordinoides #
##############

#remove minor allele frequency < 5% threshold
ord.maf <- informloci(ord.snps, MAF=0.02)

ord.maf.popnames <- c("Henry", "Lopez", "San Juan", "Orcas", "Shaw", "Waldron")
popNames(ord.maf) <- ord.maf.popnames


# dapc under island assignments
ord.maf.dapc <- dapc(x = ord.maf, pop=ord.maf@pop)
ord.maf.scatter <- scatter.dapc(ord.maf.dapc, scree.da=TRUE, scree.pca=TRUE, solid=1, 
                                posi.pca="bottomright", posi.da = "topright", 
                                col=piecolors.ord)
# dapc under found clusters
ord.maf.clusters <- find.clusters(ord.maf)
ord.maf.dapc <- dapc(x = ord.maf, pop=ord.maf.clusters$grp)
ord.maf.posterior <- ord.maf.dapc$posterior

# match coords to posterior
index <- sapply(row.names(ord.maf.posterior), function(x) {
  which(indcoords[,1] == x)
})

# dimensions of this will change depending on how many clusters you choose
ord.maf.nclust <- ncol(ord.maf.posterior)
ord.maf.data <- cbind(ord.maf.posterior, indcoords[index, -1])

# Plot data
plot(sja_shape, box=TRUE, xlab="Longitude", ylab="Latitude", buffer=FALSE)
apply(ord.maf.data,1, function(z) {
  zz <- data.matrix(z[1:ord.maf.nclust])
  index <- which(zz !=0)
  floating.pie(xpos = z[ord.maf.nclust+2], ypos = z[ord.maf.nclust+1], x = zz[index], radius = .010, col = ordinoides.col[index])})

##############################
# dapc of island populations #
##############################

ord.maf.islands <- dapc(ord.maf, pop=ord.maf@pop)

ord1 <- dapc(ord.maf, n.da=200, n.pca=10) # 7
ord1.tmp <- optim.a.score(ord1)

ord2 <- dapc(ord.maf, n.da=200, n.pca=20) # 6
ord2.tmp <- optim.a.score(ord2)

ord3 <- dapc(ord.maf, n.da=200, n.pca=30) # 7
ord3.tmp <- optim.a.score(ord3)

ord4 <- dapc(ord.maf, n.da=200, n.pca=40) # 7
ord4.tmp <- optim.a.score(ord4)

ord5 <- dapc(ord.maf, n.da=200, n.pca=50) # 7
ord5.tmp <- optim.a.score(ord5)

ord6 <- dapc(ord.maf, n.da=200, n.pca=60) # 7
ord6.tmp <- optim.a.score(ord6)


ord.maf.islands <- dapc(ord.maf, pop=ord.maf@pop)

ord.maf.islands.scatter <- scatter.dapc(ord.maf.islands, scree.pca = TRUE, 
                                        scree.da=TRUE, solid=1, posi.da = "topright", 
                                        col=piecolors.ord)

# match coords to posterior
index <- sapply(row.names(ord.maf.islands.posterior), function(x) {
  which(indcoords[,1] == x)
})

# Grab posterior as data matrix
ord.maf.islands.posterior <- ord.maf.islands$posterior

ord.maf.islands.nclust <- ncol(ord.maf.islands.posterior)
ord.maf.islands.data <- cbind(ord.maf.islands.posterior, indcoords[index, -1])

# Plot data
plot(sja_shape, box=TRUE, xlab="Longitude", ylab="Latitude", buffer=FALSE)
apply(ord.maf.islands.data,1, function(z) {
  zz <- data.matrix(z[1:ord.maf.islands.nclust])
  index <- which(zz !=0)
  floating.pie(xpos = z[ord.maf.islands.nclust+2], ypos = z[ord.maf.islands.nclust+1], x = zz[index], radius = .010, col = piecolors.ord[index])})
inset()


###########
# elegans #
###########

#remove minor allele frequency < 5% threshold
ele.maf <- informloci(ele.snps, MAF=0.02)

ele.maf.popnames <- c("San Juan", "Shaw", "Orcas")
popNames(ele.maf) <- ele.maf.popnames

# dapc under island assignments
piecolors.ele <- c("#e7298a", "#66a61e", "#7570b3")
ele.maf.dapc <- dapc(x = ele.maf, pop=ele.maf@pop)
ele.maf.scatter <- scatter.dapc(ele.maf.dapc, scree.da=TRUE, scree.pca=TRUE, solid=1, 
                                posi.pca="bottomright", posi.da = "topright", 
                                col=piecolors.ele)

# dapc under found clusters
ele.maf.clusters <- find.clusters(ele.maf)
ele.maf.dapc <- dapc(x = ele.maf, pop=ele.maf.clusters$grp)
ele.maf.posterior <- ele.maf.dapc$posterior

# match coords to posterior
index <- sapply(row.names(ele.maf.posterior), function(x) {
  which(indcoords[,1] == x)
})

# dimensions of this will change depending on how many clusters you choose
ele.maf.nclust <- ncol(ele.maf.posterior)
ele.maf.data <- cbind(ele.maf.posterior, indcoords[index, -1])

# Plot data
plot(sja_shape, box=TRUE, xlab="Longitude", ylab="Latitude", buffer=FALSE)
apply(ele.maf.data,1, function(z) {
  zz <- data.matrix(z[1:ele.maf.nclust])
  index <- which(zz !=0)
  floating.pie(xpos = z[ele.maf.nclust+2], ypos = z[ele.maf.nclust+1], x = zz[index], radius = .010, col = elegans.col[index])})

############
# sirtalis #
############

#remove minor allele frequency < 5% threshold
sir.maf <- informloci(sir.snps, MAF=0.02)

sir.maf.popnames <- c("Henry", "Orcas", "Lopez", "San Juan", "Shaw", "Waldron")

popNames(sir.maf) <- sir.maf.popnames

# dapc under island assignments
piecolors.ord <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02")
piecolors.sir <- c("#1b9e77", "#7570b3", "#d95f02","#e7298a", "#66a61e", "#e6ab02")
# piecolors.sir <- c("#e6ab02","#66a61e","#e7298a","#7570b3")
sir.maf.dapc <- dapc(x = sir.maf, pop=sir.maf@pop)
sir.maf.scatter <- scatter.dapc(sir.maf.dapc, scree.da=TRUE, scree.pca=TRUE, solid=1, 
                                posi.pca="bottomleft", posi.da = "topleft", 
                                col=piecolors.sir)

# dapc under found clusters
sir.maf.clusters <- find.clusters(sir.maf)
sir.maf.dapc <- dapc(x = sir.maf, pop=sir.maf.clusters$grp)
sir.maf.posterior <- sir.maf.dapc$posterior

# match coords to posterior
index <- sapply(row.names(sir.maf.posterior), function(x) {
  which(indcoords[,1] == x)
})

# dimensions of this will change depending on how many clusters you choose
sir.maf.nclust <- ncol(sir.maf.posterior)
sir.maf.data <- cbind(sir.maf.posterior, indcoords[index, -1])

# Plot data
plot(sja_shape, box=TRUE, xlab="Longitude", ylab="Latitude", buffer=FALSE)
apply(sir.maf.data,1, function(z) {
  zz <- data.matrix(z[1:sir.maf.nclust])
  index <- which(zz !=0)
  floating.pie(xpos = z[sir.maf.nclust+2], ypos = z[sir.maf.nclust+1], x = zz[index], radius = .010, col = sirtalis.col[index])})

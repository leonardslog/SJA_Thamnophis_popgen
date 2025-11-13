library(adegenet)
library(dartR)
library(poppr)

###############################
# IBD calculation using dartR #
###############################

# locality data
inddata <- read.csv("~/Documents/sji_thamnophis/sample_localities.csv")

# ordinoides 
ord.genind <- import2genind("~/Documents/sji_thamnophis/stacks/ref_sirtalis/ordinoides/ordinoides_SJI/unlinked_p1/ord_usnps_p1r80.gen")

ord.keep <- which(isPoly(ord.genind))
ord.keep <-names(ord.keep)
ord.polymorphic <- ord.genind[loc=ord.keep]

#remove minor allele frequency < 5% threshold
ord.maf <- poppr::informloci(ord.polymorphic, cutoff=0.05)
ord.maf.popnames <- c("Henry", "Lopez", "Orcas","San Juan", "Shaw", "Waldron")
popNames(ord.maf) <- ord.maf.popnames

ord.gl <- gi2gl(ord.maf, parallel = FALSE, verbose = NULL)

# ind coordinate data
ord.data <- subset(inddata, inddata$species=="Thamnophis ordinoides")
ord.coords <- ord.data[,c(7,4,5)]

coords.reorder <- ord.coords[match(ord.gl@ind.names, ord.coords$dapc),]
coords.reorder <- coords.reorder[,2:3]
colnames(coords.reorder) <- c("lat","lon")

ord.ibd <- gl.ibd(x=ord.gl, distance="propShared", coordinates = coords.reorder, 
                  plot.out=TRUE)


# sirtalis
sir.genind <- import2genind("~/Documents/sji_thamnophis/stacks/ref_sirtalis/sirtalis/sirtalis_SJI/unlinked_p1/sir_usnps_p1r80.gen")

sir.keep <- which(isPoly(sir.genind))
sir.keep <-names(sir.keep)
sir.polymorphic <- sir.genind[loc=sir.keep]

#remove minor allele frequency < 5% threshold
sir.maf <- informloci(sir.polymorphic, cutoff=0.05)
sir.maf.popnames <- c("Henry", "Orcas", "Lopez", "San Juan", "Shaw", "Waldron")
popNames(sir.maf) <- sir.maf.popnames

sir.gl <- gi2gl(sir.maf, parallel = FALSE, verbose = NULL)

# ind coordinate data
sir.data <- subset(inddata, inddata$species=="Thamnophis sirtalis")
sir.coords <- sir.data[,c(7,4,5)]

coords.reorder <- sir.coords[match(sir.gl@ind.names, sir.coords$dapc),]
coords.reorder <- coords.reorder[,2:3]
colnames(coords.reorder) <- c("lat","lon")

sir.ibd <- gl.ibd(x=sir.gl, distance="propShared", coordinates = coords.reorder, 
                  plot.out=TRUE)


# elegans #
ele.genind <- import2genind("~/Documents/sji_thamnophis/stacks/ref_sirtalis/elegans/elegans_SJI/unlinked_p1/ele_usnps_p1r80.gen")

ele.keep <- which(isPoly(ele.genind))
ele.keep <-names(ele.keep)
ele.polymorphic <- ele.genind[loc=ele.keep]

#remove minor allele frequency < 5% threshold
ele.maf <- informloci(ele.polymorphic, cutoff=0.05)
ele.maf.popnames <- c("San Juan", "Shaw", "Orcas")
popNames(ele.maf) <- ele.maf.popnames

ele.gl <- gi2gl(ele.maf, parallel = FALSE, verbose = NULL)

# ind coordinate data
ele.data <- subset(inddata, inddata$species=="Thamnophis elegans")
ele.coords <- ele.data[,c(7,4,5)]

coords.reorder <- ele.coords[match(ele.gl@ind.names, ele.coords$dapc),]
coords.reorder <- coords.reorder[,2:3]
colnames(coords.reorder) <- c("lat","lon")

ele.ibd <- gl.ibd(x=ele.gl, distance="propShared", coordinates = coords.reorder, 
                  plot.out=TRUE)

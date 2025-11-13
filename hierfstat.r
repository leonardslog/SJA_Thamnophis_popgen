library(hierfstat)
library(adegenet)
setwd("~/Documents/sji_thamnophis/")
sample_info <- read.csv("sample_localities.csv")

# ordinoides
vcf <- read.VCF("stacks/ref_sirtalis/elegans/elegans_SJI/unlinked_p1/populations.snps.vcf")

ord.data <- adegenet::read.structure("stacks/ref_sirtalis/ordinoides/ordinoides_SJI/unlinked_p1/ord_usnps_p1r80.str", 
                                 n.ind = 58,
                                 n.loc = 3813,
                                 onerowperind = FALSE,
                                 col.lab = 0,
                                 col.pop = 2
)

ord.pops <- data.frame(ord.data@pop)
ord.loci <- data.frame(ord.data@tab)
ord.cbind.data <- cbind(ord.pops,ord.loci)


sir.data <- adegenet::read.structure("stacks/ref_sirtalis/sirtalis/sirtalis_SJI/unlinked_p1/sir_usnps_p1r80.str", 
                                     n.ind = 35,
                                     n.loc = 4414,
                                     onerowperind = FALSE,
                                     col.lab = 0,
                                     col.pop = 2
)

sir.pops <- data.frame(sir.data@pop)
sir.loci <- data.frame(sir.data@tab)
sir.cbind.data <- cbind(sir.pops,sir.loci)


ele.data <- adegenet::read.structure("stacks/ref_sirtalis/elegans/elegans_SJI/unlinked_p1/ele_usnps_p1r80.str", 
                                     n.ind = 15,
                                     n.loc = 2919,
                                     onerowperind = FALSE,
                                     col.lab = 0,
                                     col.pop = 2
)

ele.pops <- data.frame(ele.data@pop)
ele.loci <- data.frame(ele.data@tab)
ele.cbind.data <- cbind(ele.pops,ele.loci)

ord.wc84 <- genet.dist(ord.cbind.data, diploid=TRUE, method="WC84")

# ordinoides #
#         Henry      Lopez    SanJuan      Orcas       Shaw
# Lopez   0.12130006                                            
# SanJuan 0.13505481 0.16987480                                 
# Orcas   0.16396807 0.21799114 0.12381882                      
# Shaw    0.16472527 0.21479362 0.13411228 0.08993209           
# Waldron 0.11604055 0.10694933 0.15720808 0.19759608 0.19658755

# mean = 0.1539968

sir.wc84 <- genet.dist(sir.cbind.data, diploid=TRUE, method="WC84")

# sirtalis #
#         Henry     Orcas     Lopez   SanJuan      Shaw
# Orcas   0.2130876                                        
# Lopez   0.2370408 0.1424143                              
# SanJuan 0.0000000 0.1377069 0.2137424                    
# Shaw    0.2702384 0.1789503 0.1344992 0.3083050          
# Waldron 0.2357994 0.1155894 0.1322017 0.1250459 0.1972688

# mean = 0.176126

ele.wc84 <- genet.dist(ele.cbind.data, diploid=TRUE, method="WC84")

# elegans #
#       SanJuan      Shaw
# Shaw  0.1928937          
# Orcas 0.2266491 0.2208097

# mean = 0.2134508

ord.basic.stats <- basic.stats(ord.cbind.data, diploid = TRUE, digits=4)
# $overall
# Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.3651  0.2164  0.2493  0.0329  0.2559  0.0395  0.1320  0.1543 -0.6867  0.0504 

sir.basic.stats <- basic.stats(sir.cbind.data, diploid = TRUE, digits=4)
# $overall
# Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.3731  0.2223  0.2679  0.0456  0.2770  0.0547  0.1701  0.1974 -0.6787  0.0703 

ele.basic.stats <- basic.stats(ele.cbind.data, diploid = TRUE, digits=4)
# $overall
# Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.4158  0.2608  0.3054  0.0446  0.3277  0.0669  0.1460  0.2041 -0.5941  0.0905 

####################
# 3 pop comparison #
####################

ord.data <- adegenet::read.structure("stacks/ref_sirtalis/ordinoides/ordinoides_3islands/unlinked_p1/populations.str", 
                                     n.ind = 40,
                                     n.loc = 2335,
                                     onerowperind = FALSE,
                                     col.lab = 0,
                                     col.pop = 2
)

ord.pops <- data.frame(ord.data@pop)
ord.loci <- data.frame(ord.data@tab)
ord.3islands.cbind.data <- cbind(ord.pops,ord.loci)


sir.data <- adegenet::read.structure("stacks/ref_sirtalis/sirtalis/sirtalis_3islands/unlinked_p1/populations.str", 
                                     n.ind = 16,
                                     n.loc = 3117,
                                     onerowperind = FALSE,
                                     col.lab = 0,
                                     col.pop = 2
)

sir.pops <- data.frame(sir.data@pop)
sir.loci <- data.frame(sir.data@tab)
sir.3islands.cbind.data <- cbind(sir.pops,sir.loci)


ele.data <- adegenet::read.structure("stacks/ref_sirtalis/elegans/elegans_SJI/unlinked_p1/ele_usnps_p1r80.str", 
                                     n.ind = 15,
                                     n.loc = 2919,
                                     onerowperind = FALSE,
                                     col.lab = 0,
                                     col.pop = 2
)

ele.pops <- data.frame(ele.data@pop)
ele.loci <- data.frame(ele.data@tab)
ele.cbind.data <- cbind(ele.pops,ele.loci)

##########

ord.3islands.wc84 <- genet.dist(ord.3islands.cbind.data, diploid=TRUE, method="WC84")

# ordinoides #
# SanJuan     Orcas
# Orcas 0.1367200          
# Shaw  0.1481722 0.1151375

# mean = 

sir.3islands.wc84 <- genet.dist(sir.3islands.cbind.data, diploid=TRUE, method="WC84")

# sirtalis #
# SanJuan       Shaw
# Shaw  0.27854376           
# Orcas 0.06861188 0.18633204

# mean = 

ele.wc84 <- genet.dist(ele.cbind.data, diploid=TRUE, method="WC84")

# elegans #
#       SanJuan      Shaw
# Shaw  0.1928937          
# Orcas 0.2266491 0.2208097

# mean = 0.2134508

ord.3islands.basic.stats <- basic.stats(ord.3islands.cbind.data, diploid = TRUE, digits=4)
# $overall
# Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.3798  0.2348  0.2590  0.0242  0.2711  0.0363  0.0933  0.1338 -0.6173  0.0474 

sir.3islands.basic.stats <- basic.stats(sir.3islands.cbind.data, diploid = TRUE, digits=4)
# $overall
# Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.4149  0.2582  0.2949  0.0367  0.3133  0.0551  0.1244  0.1757 -0.6067  0.0742 

ele.basic.stats <- basic.stats(ele.cbind.data, diploid = TRUE, digits=4)
# $overall
# Ho      Hs      Ht     Dst     Htp    Dstp     Fst    Fstp     Fis    Dest 
# 0.4158  0.2608  0.3054  0.0446  0.3277  0.0669  0.1460  0.2041 -0.5941  0.0905 

###################################################################################
# R code for 2018 work on Lednia tetonica etc., created by S. Hotaling. 
# Title: "Mountain stoneflies may tolerate warming streams: evidence from organismal physiology and gene expression" 
# Journal: Global Change Biology, 2020
# Authors: Scott Hotaling, Alisha A. Shah, et. al.
# R version 3.6.1 (2019-07-05)
###################################################################################

# Read in the data
setwd('/Users/scotthotaling/Desktop/Stonefly_RNAseq/EdgeR/Test_for_GitHub/')
data <- read.csv("gene_count_matrix_N17.csv",row.names=1) # removed sample #9_GATCGT


#############################
#### Change Column Names ####
#############################

#control
names(data)[1] <- "1_TTAGGC_rep1"
names(data)[2] <- "10_ACTGTG_rep1"
names(data)[3] <- "12_CATGCA_rep1"
names(data)[4] <- "15_ACTCAG_rep1"
names(data)[5] <- "16_GTGGCC_rep1"
names(data)[6] <- "17_GATCCA_rep1"
names(data)[7] <- "18_TACAGC_rep1"
names(data)[8] <- "4_TGTCAG_rep1"

#treatment
names(data)[9] <- "11_GTTCGA_rep1"
names(data)[10] <- "13_AGACTC_rep1"
names(data)[11] <- "14_ATCTGG_rep1"
names(data)[12] <- "2_CTAGCA_rep1"
names(data)[13] <- "3_AACTCG_rep1"
names(data)[14] <- "5_GCCAAT_rep1"
names(data)[15] <- "6_CAGATC_rep1"
names(data)[16] <- "7_CGTACG_rep1"
names(data)[17] <- "8_ATGTCG_rep1"


##############################
#### Package Installation ####
##############################

# install EdgeR and limma, only have to do this the first time unless you update the version of R

# BiocManager::install("biocLite")
# BiocManager::install("edgeR")
# BiocManager::install("limma")
# BiocManager::install("RUVSeq")
# install.packages("openxlsx")
# install.packages("tximport")


#######################
#### Load Packages ####
#######################

require("limma")
require("edgeR")
require("RUVSeq")
require("openxlsx")


#### Filter ####

dim(data)
filter <- apply(data, 1, function(x) length(x[x>5])>=2)
	# Filtered out genes by requiring > 5 reads and at least 2 samples
filtered <- data[filter,]
dim(filtered)
	# 52,954 // 17 samples


#### Group/organize ####

x <- as.factor(c(rep("Control",8),rep("Treatment",9)))
set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))
set


#### Make some plots of the raw data ####

library(RColorBrewer)
colors <- brewer.pal(3, "Set2")
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set, col=colors[x], cex=1.2)


#### Normalize the data and re-plot ####

set2 <- betweenLaneNormalization(set, which="upper")
	# Normalization method that forces the median of the upper quartile in the data to be the same
	# The two groups here are treatment/control (not different lanes)
	
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x])
	# Plotting the Relative Log Expression (RLE) of each library
	
plotPCA(set2, col=colors[x], cex=1.2)


#### Control for unwanted variation with RUVseq ####

design <- model.matrix(~x, data=pData(set))
	# x is the treatment groups described above
y <- DGEList(counts=counts(set), group=x)
	# perform DGE analysis and call it y
y <- calcNormFactors(y, method="upperquartile")
	# normalize libraries
	
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

fit <- glmFit(y, design)
	# Fit a Generalized Linear Model  to the data
	# Good for non-normally distributed response data
	
lrt <- glmLRT(fit, coef=2)
	# Apply a Likelihood-ratio Test (LRT) to the GLM

top <- topTags(lrt, n=nrow(set))$table
empirical <- rownames(set)[which(!(rownames(set) %in% rownames(top)[1:5000]))]

set3 <- RUVg(set, empirical, k=1)

pData(set3)
dim(set3)
	## 52,954 transcripts // 17 samples


#### Plot the data with unwanted variation removed ####
plotRLE(set3, outline=FALSE, ylim=c(-4, 4), col=colors[x])
plotPCA(set3, col=colors[x], cex=1.2)


#### Identify DEGs across groups ####

# Set groups
group <- c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2)
y <- DGEList(counts=counts(set3), group=group)

y <- calcNormFactors(y) # calculates normalization factors to scale raw library sizes, TMM is default method
y$samples ## Shows library sizes and normalization factors


#### Make various groups ####

# Gusher Population
# Gusher = Mt. St. John
gusher <- set3[,c(1,8,12,13,16)]

# Lunch Creek (also L. tumana only)
lunch <- set3[,c(2,4,7,9,11,17)]
lunch_control <- set3[,c(2,4,7)]
lunch_treatment <- set3[,c(9,11,17)]

# Tetonica Pond
pond <- set3[,c(3,5,6,10,14,15)]

# All lednia
lednia <- set3[,1:17]

# Tetonica only
# "Gusher" and "Tetonica Pond"
tetonica <- set3[,c(1,3,5,6,8,10,12,13,14,15,16)]
tetonica_control <- set3[,c(1,3,5,6,8)]
tetonica_treatment <- set3[,c(10,12,13,14,15,16)]


########################################################################################################################
#### ALL LEDNIA #########################################################################################################
########################################################################################################################

#### Create groups and filter ####
# Create grouping factor
group_lednia = c(rep("control",8),rep("treatment",9))

# Check dimensions
dim(lednia)
	# 52954 // 17 samples

total_lednia <- lednia[rowSums(counts(lednia)) > 0,] # trim out any rows with no reads counts
dim(total_lednia)
	# 52954 // 17 samples
	# no zeros because all samples are included


#### Make DGE list ####

# DGE list for lednia only
y_lednia <- DGEList(counts=counts(total_lednia), group=group_lednia)

y_lednia <- calcNormFactors(y_lednia) # calculates normalization factors to scale raw library sizes, TMM is default method

y_lednia$samples # shows library sizes and normalization factors

# pull CPM matrix for all samples
lednia_cpm <- cpm(y_lednia, normalized.lib.sizes=TRUE, log=TRUE)
# write.csv(x=lednia_cpm, file = "Lednia_CPM_Matrix-Oct19.csv")

# Getting raw normalized counts for each sample
# write.csv(x=y_lednia$counts, file = "NormalizedCounts_AllLednia-Oct19.csv")

# create color vector for graph
colors_lednia <- c(rep("blue",8),rep("red",9))

# las=2 makes label text perpendicular to axis
par(las=2)

# generate plot
pdf("LibrarySize_Lednia2.pdf")
library_lednia <- barplot(y_lednia$samples$lib.size, col=colors_lednia, cex.axis=.8, ylab = "Library Size", xlab = "Sample")
legend("topright", cex = 0.82, legend=c("Control", "Treatment"), bty = "n", pt.cex=1.5, fill=c("blue", "red"))
print(library_lednia)
dev.off()


########## Design Matrix ###########

fac_lednia <- c(rep("Control",8), rep("Treatment",9)) # use this to get groups
fac_lednia <- factor(fac_lednia) # renames fac names with numbers

design_lednia <- model.matrix(~0+fac_lednia) # 0 is the intercept
colnames(design_lednia) <- levels(fac_lednia)
design_lednia


####### Estimating Dispersions - CR method ########

# Uses GLM instead of qCML method because testing multiple factors here

# Uses Cox-Reid profile-adjusted likelihood

# estimates common dispersion and tagwise dispersions in one run
y_lednia <- estimateDisp(y_lednia, design_lednia)

# common dispersions
y_lednia$common.dispersion
	# Dispersion = 1.401553

# tagwise (gene-specific dispersions)
summary(y_lednia$tagwise.dispersion)

# estimated prior degrees of freedom
y_lednia$prior.df

# output logCPM for each sample (this is overall for everything, not by sample)
write.csv(x=y_lednia$AveLogCPM, file = "LogCPM_AllLednia.csv")


####### BCV Plot ##########

# BCV is the coefficient of variation with which the (unknown) true abundance of the gene varies between replicate RNA samples

pdf("BCV_Lednia2.pdf")
BCV_lednia <- plotBCV(y_lednia)
print(BCV_lednia)
dev.off()


####### DGE - Quasi-Likelihood F-test #########

# quasi-likelihood F-test better for bulk RNA-seq data because stricter error rate control, accounts for
		# uncertainty in dispersion estimation

group_lednia <- factor(c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2))
# group_lednia_individuals <- factor(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17))

design_lednia <- model.matrix(~0+group_lednia)
# design_lednia_ind <- model.matrix(~0+group_lednia_individuals)

# QL model representing the study design fitted to the data
fit_lednia <- glmQLFit(y_lednia, design_lednia)
# fit_lednia_ind <- glmQLFit(y_lednia, design_lednia_ind)

# tests use FDR < 0.05

# compare Control (1) vs Treatment (2)
# for top 10000 genes
qlf_lednia <- glmQLFTest(fit_lednia, contrast=c(-1,1))
control_v_treat_lednia <- topTags(qlf_lednia, n = 10000, p.value=0.05)
control_v_treat_lednia <- data.frame(control_v_treat_lednia)
#require(xlsx)
write.csv(x=control_v_treat_lednia, file = "LedniaAll_DEGs.csv")

# Might want glf_lednia$fitted.values for normalized read counts...
qlf_lednia$fitted.values


########################################################################################################################
#### TETONICA ONLY #########################################################################################################
########################################################################################################################

#### Create groups and filter ####

# Create grouping factor
group_tetonica = c(rep("Control",5),rep("Treatment",6))

# Check dimensions
dim(tetonica)
	# 52954 // 11 samples

total_tetonica <- tetonica[rowSums(counts(tetonica)) > 0,] # trim out any rows with no reads counts
dim(total_tetonica)
	# 52280 // 11 samples


#### Make DGE list and normalize ####

# DGE list for tetonica only
y_tetonica <- DGEList(counts=counts(total_tetonica), group=group_tetonica)

y_tetonica <- calcNormFactors(y_tetonica) # calculates normalization factors to scale raw library sizes, TMM is default method

y_tetonica$samples # shows library sizes and normalization factors


# create color vector for graph
colors_tetonica <- c(rep("blue",5),rep("red",6))

# las=2 makes label text perpendicular to axis
par(las=2)

# generate plot
pdf("LibrarySize_Tetonica2.pdf")
library_tetonica <- barplot(y_tetonica$samples$lib.size, col=colors_tetonica, cex.axis=.8, ylab = "Library Size", xlab = "Sample")
legend("topright", cex = 0.82, legend=c("Control", "Treatment"), bty = "n", pt.cex=1.5, fill=c("blue", "red"))
print(library_tetonica)
dev.off()


########## Design Matrix ###########

fac_tetonica <- c(rep("Control",5), rep("Treatment",6)) # use this to get groups
fac_tetonica <- factor(fac_tetonica) # renames fac names with numbers

design_tetonica <- model.matrix(~0+fac_tetonica) # 0 is the intercept
colnames(design_tetonica) <- levels(fac_tetonica)
design_tetonica


####### Estimating Dispersions - CR method ########

# Uses GLM instead of qCML method because testing multiple factors here

# Uses Cox-Reid profile-adjusted likelihood

# estimates common dispersion and tagwise dispersions in one run
y_tetonica <- estimateDisp(y_tetonica, design_tetonica)

# common dispersions
y_tetonica$common.dispersion
	# Dispersion = 0.9174679

# tagwise (gene-specific dispersions)
summary(y_tetonica$tagwise.dispersion)

# estimated prior degrees of freedom
y_tetonica$prior.df


####### BCV Plot ##########

# BCV is the coefficient of variation with which the (unknown) true abundance of the gene varies between replicate RNA samples

pdf("BCV_Tetonica2.pdf")
BCV_tetonica <- plotBCV(y_tetonica)
print(BCV_tetonica)
dev.off()


####### DGE - Quasi-Likelihood F-test #########

# quasi-likelihood F-test better for bulk RNA-seq data because stricter error rate control, accounts for
		# uncertainty in dispersion estimation

group_tetonica <- factor(c(1,1,1,1,1,2,2,2,2,2,2))

design_tetonica <- model.matrix(~0+group_tetonica)

# QL model representing the study design fitted to the data
fit_tetonica <- glmQLFit(y_tetonica, design_tetonica)

# tests use FDR < 0.05

# compare Control (1) vs Treatment (2)
qlf_tetonica <- glmQLFTest(fit_tetonica, contrast=c(-1,1))
control_v_treat_tetonica <- topTags(qlf_tetonica, n = 10000, p.value=0.05)
control_v_treat_tetonica <- data.frame(control_v_treat_tetonica)

#require(xlsx)
write.csv(x=control_v_treat_tetonica, file = "Tetonica_DEGs.csv")


########################################################################################################################
#### GUSHER = Mt. St. John ########################################################################################################
########################################################################################################################

#### Create groups and filter ####

# Create grouping factor
group_gusher = c(rep("control",2),rep("treatment",3))

# Check dimensions
dim(gusher)
	# 52954 // 5 samples

total_gusher <- gusher[rowSums(counts(gusher)) > 0,] # trim out any rows with no reads counts
dim(total_gusher)
	# 49423 // 5 samples


#### Make DGE list and normalize ####

# DGE list for Gusher only
y_gusher <- DGEList(counts=counts(total_gusher), group=group_gusher)

y_gusher <- calcNormFactors(y_gusher) # calculates normalization factors to scale raw library sizes, TMM is default method

y_gusher$samples # shows library sizes and normalization factors


# create color vector for graph
colors_gusher <- c(rep("blue",2),rep("red",3))

# las=2 makes label text perpendicular to axis
par(las=2)

# generate plot
pdf("Library_Gusher.pdf")
library_gusher <- barplot(y_gusher$samples$lib.size, col=colors_gusher, cex.axis=.8, ylab = "Library Size", xlab = "Sample")
legend("topright", cex = 0.82, legend=c("Control", "Treatment"), bty = "n", pt.cex=1.5, fill=c("blue", "red"))
print(library_gusher)
dev.off()


########## Design Matrix ###########

fac_gusher <- c(rep("Control",2), rep("Treatment",3)) # use this to get groups

fac_gusher <- factor(fac_gusher) # renames fac names with numbers

design_gusher <- model.matrix(~0+fac_gusher) # 0 is the intercept

colnames(design_gusher) <- levels(fac_gusher)

design_gusher


####### Estimating Dispersions - CR method ########

# Uses GLM instead of qCML method because testing multiple factors here

# Uses Cox-Reid profile-adjusted likelihood

# estimates common dispersion and tagwise dispersions in one run
y_gusher <- estimateDisp(y_gusher, design_gusher)

# common dispersions
y_gusher$common.dispersion
	# Dispersion = 0.5657731

# tagwise (gene-specific dispersions)
summary(y_gusher$tagwise.dispersion)

# estimated prior degrees of freedom
y_gusher$prior.df


####### BCV Plot ##########

# BCV is the coefficient of variation with which the (unknown) true abundance of the gene varies between replicate RNA samples

pdf("BCV_Gusher.pdf")
BCV_gusher <- plotBCV(y_gusher)
print(BCV_gusher)
dev.off()


####### DGE - Quasi-Likelihood F-test #########

# quasi-likelihood F-test better for bulk RNA-seq data because stricter error rate control, accounts for
		# uncertainty in dispersion estimation

group_gusher <- factor(c(1,1,2,2,2))

design_gusher <- model.matrix(~0+group_gusher)

# QL model representing the study design fitted to the data
fit_gusher <- glmQLFit(y_gusher, design_gusher)

# tests use FDR < 0.05

# compare Control (1) vs Treatment (2)
qlf_gusher <- glmQLFTest(fit_gusher, contrast=c(-1,1))
control_v_treat_gusher <- topTags(qlf_gusher, n = 10000, p.value=0.05)
control_v_treat_gusher <- data.frame(control_v_treat_gusher)
#require(xlsx)
write.csv(x=control_v_treat_gusher, file = "Gusher_DEGs.csv")


########################################################################################################################
#### LUNCH = Lunch Creek ########################################################################################################
########################################################################################################################

#### Create groups and filter ####

# Create grouping factor
group_lunch = c(rep("control",3),rep("treatment",3))

# Check dimensions
dim(lunch)
	# 52954 // 6 samples

total_lunch <- lunch[rowSums(counts(lunch)) > 0,] # trim out any rows with no reads counts
dim(total_lunch)
	# 52306 // 6 samples


#### Make DGE list and normalize ####

# DGE list for lunch only
y_lunch <- DGEList(counts=counts(total_lunch), group=group_lunch)

y_lunch <- calcNormFactors(y_lunch) # calculates normalization factors to scale raw library sizes, TMM is default method

y_lunch$samples # shows library sizes and normalization factors


# create color vector for graph
colors_lunch <- c(rep("blue",3),rep("red",3))

# las=2 makes label text perpendicular to axis
par(las=2)

# generate plot
pdf("LibrarySize_Lunch2.pdf")
library_lunch <- barplot(y_lunch$samples$lib.size, col=colors_lunch, cex.axis=.8, ylab = "Library Size", xlab = "Sample")
legend("topright", cex = 0.82, legend=c("Control", "Treatment"), bty = "n", pt.cex=1.5, fill=c("blue", "red"))
print(library_lunch)
dev.off()


########## Design Matrix ###########

fac_lunch <- c(rep("Control",3), rep("Treatment",3)) # use this to get groups

fac_lunch <- factor(fac_lunch) # renames fac names with numbers

design_lunch <- model.matrix(~0+fac_lunch) # 0 is the intercept

colnames(design_lunch) <- levels(fac_lunch)

design_lunch


####### Estimating Dispersions - CR method ########

# Uses GLM instead of qCML method because testing multiple factors here

# Uses Cox-Reid profile-adjusted likelihood

# estimates common dispersion and tagwise dispersions in one run
y_lunch <- estimateDisp(y_lunch, design_lunch)

# common dispersions
y_lunch$common.dispersion
	# Dispersion = 1.442635

# tagwise (gene-specific dispersions)
summary(y_lunch$tagwise.dispersion)

# estimated prior degrees of freedom
y_lunch$prior.df


####### BCV Plot ##########

# BCV is the coefficient of variation with which the (unknown) true abundance of the gene varies between replicate RNA samples

pdf("BCV_Lunch2.pdf")
BCV_lunch <- plotBCV(y_lunch)
print(BCV_lunch)
dev.off()


####### DGE - Quasi-Likelihood F-test #########

# quasi-likelihood F-test better for bulk RNA-seq data because stricter error rate control, accounts for
		# uncertainty in dispersion estimation

group_lunch <- factor(c(1,1,1,2,2,2))

design_lunch <- model.matrix(~0+group_lunch)

# QL model representing the study design fitted to the data
fit_lunch <- glmQLFit(y_lunch, design_lunch)

# tests use FDR < 0.05

# compare Control (1) vs Treatment (2)
qlf_lunch <- glmQLFTest(fit_lunch, contrast=c(-1,1))
control_v_treat_lunch <- topTags(qlf_lunch, n = 10000, p.value=0.05)
control_v_treat_lunch <- data.frame(control_v_treat_lunch)
#require(xlsx)
write.csv(x=control_v_treat_lunch, file = "LunchCreek_DEGs.csv")


########################################################################################################################
#### POND = Tetonica Pond #######################################################################################################
########################################################################################################################

#### Create groups and filter ####

# Create grouping factor
group_pond = c(rep("control",3),rep("treatment",3))

# Check dimensions
dim(pond)
	# 52954 // 6 samples

total_pond <- pond[rowSums(counts(pond)) > 0,] # trim out any rows with no reads counts
dim(total_pond)
	# 51796 // 6 samples


#### Make DGE list and normalize ####

# DGE list for pond only
y_pond <- DGEList(counts=counts(total_pond), group=group_pond)

y_pond <- calcNormFactors(y_pond) # calculates normalization factors to scale raw library sizes, TMM is default method

y_pond$samples # shows library sizes and normalization factors


# create color vector for graph
colors_pond <- c(rep("blue",3),rep("red",3))

# las=2 makes label text perpendicular to axis
par(las=2)

# generate plot
pdf("LibrarySize_Pond2.pdf")
library_pond <- barplot(y_pond$samples$lib.size, col=colors_pond, cex.axis=.8, ylab = "Library Size", xlab = "Sample")
legend("topright", cex = 0.82, legend=c("Control", "Treatment"), bty = "n", pt.cex=1.5, fill=c("blue", "red"))
print(library_pond)
dev.off()


########## Design Matrix ###########

fac_pond <- c(rep("Control",3), rep("Treatment",3)) # use this to get groups
fac_pond <- factor(fac_pond) # renames fac names with numbers
design_pond <- model.matrix(~0+fac_pond) # 0 is the intercept
colnames(design_pond) <- levels(fac_pond)
design_pond


####### Estimating Dispersions - CR method ########

# Uses GLM instead of qCML method because testing multiple factors here

# Uses Cox-Reid profile-adjusted likelihood

# estimates common dispersion and tagwise dispersions in one run
y_pond <- estimateDisp(y_pond, design_pond)

# common dispersions
y_pond$common.dispersion
	# Dispersion = 0.8336669

# tagwise (gene-specific dispersions)
summary(y_pond$tagwise.dispersion)

# estimated prior degrees of freedom
y_pond$prior.df


####### BCV Plot ##########

# BCV is the coefficient of variation with which the (unknown) true abundance of the gene varies between replicate RNA samples

pdf("BCV_Pond2.pdf")
BCV_pond <- plotBCV(y_pond)
print(BCV_pond)
dev.off()


####### DGE - Quasi-Likelihood F-test #########

# quasi-likelihood F-test better for bulk RNA-seq data because stricter error rate control, accounts for
		# uncertainty in dispersion estimation

group_pond <- factor(c(1,1,1,2,2,2))

design_pond <- model.matrix(~0+group_pond)

# QL model representing the study design fitted to the data
fit_pond <- glmQLFit(y_pond, design_pond)

# tests use FDR < 0.05

# compare Control (1) vs Treatment (2)
qlf_pond <- glmQLFTest(fit_pond, contrast=c(-1,1))
control_v_treat_pond <- topTags(qlf_pond, n = 10000, p.value=0.05)
control_v_treat_pond <- data.frame(control_v_treat_pond)
#require(xlsx)
write.csv(x=control_v_treat_pond, file = "TetonicaPond_DEGs.csv")


####################################################################################################################################################
#### Separate significantly upregulated and downregulated genes ####################################################################################
####################################################################################################################################################

# Lednia
lednia_up = rownames(control_v_treat_lednia[control_v_treat_lednia$logFC > 0,])
lednia_down = rownames(control_v_treat_lednia[control_v_treat_lednia$logFC < 0,])

# Tetonica
tetonica_up = rownames(control_v_treat_tetonica[control_v_treat_tetonica$logFC > 0,])
tetonica_down = rownames(control_v_treat_tetonica[control_v_treat_tetonica$logFC < 0,])


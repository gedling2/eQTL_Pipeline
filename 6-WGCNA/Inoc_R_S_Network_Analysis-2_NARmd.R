### Load Libraries ###
library('dplyr')
library('edgeR')
library(WGCNA)
library(cluster)
library("flashClust")
library("data.table")

### WGCNA Options ###
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

### Input Raw Counts from FeatureCounts Subread for inoculated samples ###
WGCNA_NC_data_inoc<-read.csv(url("https://osu.box.com/shared/static/e1ub70or1kzbzct56pblm2lz7mnzuazg.csv"), header = T, row.names = 1)

### Assign i to inoculated samples ###
colnames(WGCNA_NC_data_inoc)<-paste("i", colnames(WGCNA_NC_data_inoc), sep="")

### Importing external information ### 
urlBlup<-"https://osu.box.com/shared/static/wccfgknahatdq1c005gmnpvor4w02dta.csv"
urlNorm<-"https://osu.box.com/shared/static/392pk6dl3m4m75pardd4wkjvpwi24ffd.csv"
lesionValue<-"https://osu.box.com/shared/static/yb0qluwmjpqogzfsgauvdp05ojh4tenf.csv"


### Importing External Inofrmation (Normailzed P. sojae counts) ###
normValues<-read.csv(url(urlNorm, "r"), header=T)
normValues=mutate(normValues, IND=gsub("[(]", "",normValues$individuals))
normValues=mutate(normValues, IND2=gsub("[)]", ".",normValues$IND))
normValues.df<-select(normValues, -individuals, -IND)
colnames(normValues.df)<-tolower(colnames(normValues.df))
normValues.df=mutate(normValues.df, ind=tolower(ind2))
normValues.df=select(normValues.df, -ind2)

### Importing External Inofrmation (BLUPs) ###
blupValues<-read.csv(url(urlBlup, "r"))
blupValues=mutate(blupValues, IND=gsub("[(]", "",blupValues$INDIVIDUALS))
blupValues=mutate(blupValues, IND2=gsub("[)]", ".",blupValues$IND))
blupValues=mutate(blupValues, IND3=sub("([R|S])0","\\1", IND2))
blupValues.df<-select(blupValues, X.Intercept., IND3)
colnames(blupValues.df)<-tolower(c("intercept", "IND3"))
blupValues.df=mutate(blupValues.df, ind=tolower(ind3))
blupValues.df=select(blupValues.df, -ind3)

### Importing External Inofrmation (Lesion Length Values) ###
lesionValue.df<-read.csv(url(lesionValue, "r"))
lesionValue.df=mutate(lesionValue.df, IND=gsub("[(]", "",lesionValue.df$Individuals))
lesionValue.df=mutate(lesionValue.df, IND2=gsub("[)]", ".",lesionValue.df$IND))
lesionValue.df<-select(lesionValue.df,Lesion.Length, IND2)
lesionValue.df=mutate(lesionValue.df, ind=tolower(IND2)) %>%select(-IND2)

### Combine external traits into data frame ###
traits.df.tmp<-merge(blupValues.df, normValues.df)
traits.df<-merge(traits.df.tmp, lesionValue.df)
rownames(traits.df) <- traits.df$ind
traits.df$ind <- NULL

### Repeat Network Analysis 20 times ###
for (i in seq_len(20)){

print(paste(".....Starting Round.... ", i))
if (dir.exists("Ref")){
  unlink("Ref")
}else {
  dir.create("Ref")
}
pdf(paste0("Ref/Network-Analysis", i, ".pdf"),height=8,width=8)

### Normalizing raw counts using EdgeR ###
print(paste(".....EdgeR.... ", i))
counts<-WGCNA_NC_data_inoc
group<-c("R", "S",rep("R", 45), rep("S", 48))
cds<-DGEList(counts, group=group)

### Filtering cpm >1 across 50 samples ###
keep <- rowSums(cpm(cds)>1) >= 50
cds <- cds[keep, , keep.lib.sizes=FALSE]
cds <- calcNormFactors( cds)
plotMDS.DGEList( cds , main = "MDS Plot for Count Data",labels = colnames( cds$counts ) )
logcpm<-cpm(cds, prior.count = 0.25, log=T)

### Caluclating total normalized filtered count per Indivdual resistance and suceptible separatly ###
colsums_inoc<-as.data.frame(colSums(logcpm))
colsums_inoc$Sample<-rownames(colsums_inoc)
colnames(colsums_inoc)<-c("Sum", "Sample")
inoc_R_rownames<-grep("^iR", colsums_inoc$Sample )
inoc_S_rownames<-grep("^iS", colsums_inoc$Sample)
inoc_Resistance<-colsums_inoc[inoc_R_rownames,]
inoc_Susceptible<-colsums_inoc[inoc_S_rownames,]


### Calucaltuate minimum standard deviation of 30 random samples of resistant individuals and sucptible separartly ###
set.seed(1001*i)

#### Calucluate minimun standard deviation cutoff ###
sd_minimizer<-function(data, cutoff=1, rep=2000){
  sample_holder<-data.frame(Rep=as.numeric(), sd=as.numeric())
  for (i in seq_len(rep)) {
    r_sample_no<-sample(nrow(data), 30)
    pick_30_sample<-data[r_sample_no,]
    sds<-sd(pick_30_sample[,-2], na.rm = T)
    
    if (is.na(sds)){ 
      next
    } else{
      sample_holder<-rbind(sample_holder, data.frame(i, sds))
      
      
    }
  }
  
  
  sample_holder
}


##### Select 30 indidviduals with standard deviation equal to or lower to the standard deviation cutoff ####
lowest_sd_sample_picker<-function(data, cutoff=1, rep=2000){
  for (i in seq_len(rep)) {
    r_sample_no<-sample(nrow(data), 30)
    pick_30_sample<-data[r_sample_no,]
    sds<-sd(pick_30_sample[,-2], na.rm = T)
    
    if (is.na(sds)){ 
      next
    } else if (sds<=cutoff){
      sample_holder<-pick_30_sample
      return(sample_holder)
      print("Done")
      
    }
  }
  
}
### Calling standrad minimzer to get least standard devation for 10,000 repititions for susceptible and resistant indviduals separatly ###
inoc_Sus_outlier_rmd<-sd_minimizer(inoc_Susceptible,20000, 10000)
inoc_lowest_Sus_sd<-range(inoc_Sus_outlier_rmd$sds)[1]
inoc_Res_outlier_rmd<-sd_minimizer(inoc_Resistance,20000, 10000)
inoc_lowest_Res_sd<-range(inoc_Res_outlier_rmd$sds)[1]

### Calling individuals less than the standrad deviation cutoff for suceptible and resistant separatly ### 
inoc_lowest_sd_Sus<-lowest_sd_sample_picker(inoc_Susceptible,inoc_lowest_Sus_sd, 100000)
while (is.null(inoc_lowest_sd_Sus)){
  
  inoc_lowest_sd_Sus<-lowest_sd_sample_picker(inoc_Susceptible,inoc_lowest_Sus_sd, 100000)
}

inoc_lowest_sd_Res<-lowest_sd_sample_picker(inoc_Resistance,inoc_lowest_Res_sd, 1000000)

while (is.null(inoc_lowest_sd_Res)){
  
  inoc_lowest_sd_Res<-lowest_sd_sample_picker(inoc_Resistance,inoc_lowest_Res_sd, 1000000)
}

WGCNA_NC_data_inoc_R<-WGCNA_NC_data_inoc[,inoc_lowest_sd_Res$Sample]
WGCNA_NC_data_inoc_S<-WGCNA_NC_data_inoc[,inoc_lowest_sd_Sus$Sample]



### Normailzing raw counts of selected 60 Individuals (susceptible and resitant indidviduals separtaly) using EdgeR ###
fun_logcpm<-function(counts, grp, name){
 as.data.frame(colSums((counts)))
 group<-grp
 cds<-DGEList(counts, group=group)
  
### Filtering normalized counts by cpm > 1 acorss 25 samples ###  
  
  keep <- rowSums(cpm(cds)>1) >= 25
  cds <- cds[keep, , keep.lib.sizes=FALSE]
  dim(cds)
  cds <- calcNormFactors( cds)
  cds$samples
  plotMDS.DGEList( cds , main = "MDS Plot for Count Data",labels = colnames( cds$counts ) )
  logcpm<-cpm(cds, prior.count = .25, log=T)
  return(logcpm)
  
}

### Creating data strucutre of Normailzed filtered counts of 60 selected Individuals (susceptible and resitant indidviduals separtaly) ###
counts_S<-WGCNA_NC_data_inoc_S
grp_S<-rep("S", 30)

counts_R<-WGCNA_NC_data_inoc_R
grp_R<-rep("R", 30)

logcpm_S<-fun_logcpm(counts_S,grp_S, "SUS")
logcpm_R<-fun_logcpm(counts_R,grp_R, "RES")

logcpm.df.S<-as.data.frame(logcpm_S)

### Removed Indidivual 280, due to no genotype data ###
if ( "iS280" %in% names(logcpm.df.S)){
logcpm.df.S=select(logcpm.df.S, -iS280)}
logcpm.df.R<-as.data.frame(logcpm_R)

### Renaming headers ###
clean_headers<-function(logcpm.df){
colnames(logcpm.df)<-gsub("iX.", "", colnames(logcpm.df))
colnames(logcpm.df)<-tolower(colnames(logcpm.df))
colnames(logcpm.df)<-gsub("i", "",colnames(logcpm.df))
colnames(logcpm.df)
return(logcpm.df)
}
logcpm.df.S<-clean_headers(logcpm.df.S)
logcpm.df.R<-clean_headers(logcpm.df.R)

### Counts number of rows in Resistnace samples ###
dim_r_R=dim(logcpm.df.R)[1]
### Counts number of columns in Resistnace samples ###
dim_c_R=dim(logcpm.df.R)[2]

### Counts number of rows in Susceptible samples ###
dim_r_S=dim(logcpm.df.S)[1]
### Counts number of columns in Susceptible samples ###
dim_c_S=dim(logcpm.df.S)[2]

### Merging Resistant and scueptible data structures ###
interS<-intersect(rownames(logcpm.df.S), rownames(logcpm.df.R))
logcpm.df.S.fltrd<-logcpm.df.S[interS,]
logcpm.df.R.fltrd<-logcpm.df.R[interS,]
logcpm.df.combined<-cbind(logcpm.df.S.fltrd, logcpm.df.R.fltrd)


### Selecting external traits to selected 60 Individuals ####
filter_trait_data<-function(logcpm.df){
filter_traits<-traits.df[colnames(logcpm.df),]
datTraits<-filter_traits
datExpr<-as.data.frame(t(logcpm.df))
datTraits<-datTraits[rownames(datExpr),]
return(list(datTraits, datExpr))
}
datTraits<-filter_trait_data(logcpm.df.combined)[[1]]
datExpr<-filter_trait_data(logcpm.df.combined)[[2]]
### Return TRUE if datasets align correctly ###
table(rownames(datTraits)==rownames(datExpr)) 

### WGCNA Network Analysis ###

### Cluster samples by expression ###
clusret_sample_expr<-function(datExpr, datTraits, thresholdZ.k, name){

par(mar= c(6, 8.5, 3, 3))
  ### Calculates the whole network connectivity ###
A = adjacency(t(datExpr),type="signed") 
#### Standardized connectivity ###
k = as.numeric(apply(A,2,sum))-1 
Z.k = scale(k)
outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
### Convert traits to a color representation where red indicates high values ###
traitColors = data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]] = paste(names(datTraits))
datColors = data.frame(outlier = outlierColor,traitColors)

plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample Dendrogram and Trait Heatmap")

remove.samples = Z.k<thresholdZ.k | is.na(Z.k)
datExpr= datExpr[!remove.samples,]
datTraits = datTraits[!remove.samples,]
return(list(datTraits, datExpr))
}

### Default -2.5 ###
thresholdZ.k = -2.5 
datTraits.F<-clusret_sample_expr(datExpr, datTraits, thresholdZ.k, "SUS_RES")[[1]]
datExpr.F<-clusret_sample_expr(datExpr, datTraits, thresholdZ.k, "SUS_RES")[[2]]

### Step-by-dtep network construction and module detection ###
### Calculate Soft-thresholding power ###
NETWORKTYPE<-"signed hybrid"
RSQUARED_CUTOFF<-0.9

step_by_step_network_c<-function(datExpr, NETWORKTYPE, RSQUARED_CUTOFF){
  
### Choose a set of soft-thresholding powers ###
powers = c(c(1:10), seq(from = 1, to=60, by=2))

### Call the network topology analysis function ###
sft=pickSoftThreshold(datExpr,powerVector=powers, networkType=NETWORKTYPE, corFnc= "bicor",corOptions=list(use = 'p',maxPOutliers=0.1))

### Plot the results ###
par(mfrow = c(1,2));
cex1 = 0.8;

### Scale-free topology fit index as a function of the soft-thresholding power ###
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

### Corresponds to using an R^2 cut-off of h ###
abline(h=RSQUARED_CUTOFF,col="red")

### Mean connectivity as a function of the soft-thresholding power ###
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
return(sft)
}

sft<-step_by_step_network_c(datExpr.F, NETWORKTYPE, RSQUARED_CUTOFF)

### Set the Soft Power, Lowest index ###
softPower = sft$powerEstimate;

if(is.na(softPower)){
  
  softPower=14
  
}

### Save Data image and history ###
save.image(file="SoftPower.RData")
savehistory(file="SoftPower.Rhistory")


#3. Co-expression similarity and adjacency
#Following setps will take a place here
#Topology overlap Matrix
#Clustering using TOM

network_analysis<-function(datExpr, datTraits, softPower, NETWORKTYPE, MIN_MODULE_SIZE, DEEP_SPLIT, MEDISSTHRES, name){
#Calculate adjacencies

adjacencyP = adjacency(datExpr, power = softPower, type=NETWORKTYPE);


#Topological Overlap Matrix (TOM)

TOM = TOMsimilarity(adjacencyP, TOMType="signed", verbose=1);
dissTOM = 1-TOM

#Clustering using TOM
# Call the hierarchical clustering function
geneTree = flashClust(as.dist(dissTOM), method = "average");

# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04);

# We like large modules, so we set the minimum module size relatively high:
#MIN_MODULE_SIZE=30
#DEEP_SPLIT=0 #0-3 0- Smaller number of larger modules 3 Larger number smaller modules

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = DEEP_SPLIT, pamRespectsDendro = FALSE,
                            minClusterSize = MIN_MODULE_SIZE);

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)

# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#Merging of modules whose expression profiles are very similar
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
#MEDiss = 1-cor(MEs);
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = flashClust(as.dist(MEDiss), method = "average");
# Plot the result
#sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",xlab = "", sub = "")


# Plot the cut line into the dendrogram
abline(h=MEDISSTHRES, col = "red")
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDISSTHRES, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;


plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;

# Correlate traits

#Define number of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(", 
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))
#display the corelation values with a heatmap plot


labeledHeatmap(Matrix= moduleTraitCor, 
               xLabels= names(datTraits), 
               yLabels= names(MEs), 
               ySymbols= names(MEs), 
               colorLabels= FALSE, 
               colors= blueWhiteRed(50), 
               textMatrix= textMatrix, 
               setStdMargins= FALSE, 
               cex.text= 0.5, 
               zlim= c(-1,1), 
               main= paste("Module-trait relationships"))


# 3.b Gene relationship to trait and important modules: Gene Significance and Module Membership
# We quantify associations of individual genes with our trait of interest (Lesion length) by defining Gene Significance GS as (the absolute value of) 
# the correlation between the gene and the trait. For each module, 
# we also define a quantitative measure of module membership MM as the correlation of the module eigengene and the gene expression profile. 
# This allows us to quantify the similarity of all genes on the array to every module.

# Define variable lesion length containing the lesion length column of datTrait
lesion_length = as.data.frame(datTraits$Lesion.Length);
names(lesion_length) = "lesion"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, lesion_length, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(lesion_length), sep="");
names(GSPvalue) = paste("p.GS.", names(lesion_length), sep="");

#Read the annotation file
annotation<-fread("Gmax_275_Wm82.a2.v1.annotation_merged.txt", header = T)
probes = names(datExpr)
probes2annot = match(probes, annotation$locusName)

print(length(probes))
geneInfo0 = data.frame(Probe = probes, 
                  pacid = annotation$pacId[probes2annot],
                  peptideName = annotation$peptideName[probes2annot],
                  pfam = annotation$Pfam[probes2annot],
                  panther = annotation$Panther[probes2annot],
                  kog = annotation$KOG[probes2annot],
                  kegg_ec = annotation$`KEGG/ec`[probes2annot],
                  go = annotation$GO[probes2annot],
                  v2 = annotation$V2[probes2annot],
                  v3 = annotation$V3[probes2annot],
                  moduleColor = moduleColors,
                  geneTraitSignificance,
                  GSPvalue);

# Order modules by their significance for Lesion length value ###
modOrder = order(-abs(cor(MEs, lesion_length, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.lesion));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = paste0("Ref/geneInfo",i,".csv"))

eign_lesion<-cbind(MEs, datTraits)
write.csv(eign_lesion, file = paste0("Ref/eign_lesion",i,".csv"))

}


MEDISSTHRES_R = 0.1
MEDISSTHRES_S = 0.1
MIN_MODULE_SIZE=30
#0-3 0- Smaller number of larger modules 3 Larger number smaller modules
DEEP_SPLIT=0 
network_analysis(datExpr.F, datTraits.F, softPower, NETWORKTYPE, MIN_MODULE_SIZE, DEEP_SPLIT, MEDISSTHRES_R, "SUS_RES")
dev.off();
print(paste("Finish Round", i))
}



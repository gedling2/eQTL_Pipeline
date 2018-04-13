library('dplyr')
library('edgeR')
library(WGCNA)
library(cluster)
require("flashClust")
library(data.table)

### WGCNA Options ###
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

### Input Raw Counts from FeatureCounts Subread for inoculated samples ###
WGCNA_NC_data_inoc<-read.csv(url("https://osu.box.com/shared/static/e1ub70or1kzbzct56pblm2lz7mnzuazg.csv"), header = T, row.names = 1)

### Assign "i" to inoculated samples ###
colnames(WGCNA_NC_data_inoc)<-paste("i", colnames(WGCNA_NC_data_inoc), sep="")

### Importing external information ###
urlBlup<-"https://osu.box.com/shared/static/e74v8c4nzbr5asj3qydkflcef8mlqh8y.csv"
urlNorm<-"https://osu.box.com/shared/static/govplg315eho67qxxpi8ja6ntswuikph.csv"
lesionValue<-"https://osu.box.com/shared/static/yb0qluwmjpqogzfsgauvdp05ojh4tenf.csv"
if(!file.exists("lesion.csv")){
  download.file(lesionValue, "lesion.csv", method = "wget")
}

### Importing External Inofrmation (Normailzed P. sojae counts) ###
sample_boostraping<-function(i, urlBlup, urlNorm, WGCNA_NC_data_inoc){
normValues<-read.csv(url(urlNorm, "r"), header=T)
normValues=mutate(normValues, IND=gsub("[(]", "",normValues$individuals))
normValues=mutate(normValues, IND2=gsub("[)]", ".",normValues$IND))
normValues.df<-select(normValues, -individuals, -IND)
colnames(normValues.df)<-tolower(colnames(normValues.df))
normValues.df=mutate(normValues.df, ind=tolower(ind2))
normValues.df=select(normValues.df, -ind2)
normValues.df

### Importing External Inofrmation (BLUPs) ###
blupValues<-read.csv(url(urlBlup, "r"))
blupValues=mutate(blupValues, IND=gsub("[(]", "",blupValues$INDIVIDUALS))
blupValues=mutate(blupValues, IND2=gsub("[)]", ".",blupValues$IND))
blupValues=mutate(blupValues, IND3=sub("([R|S])0","\\1", IND2))
blupValues.df<-select(blupValues, X.Intercept., IND3)
colnames(blupValues.df)<-tolower(c("intercept", "IND3"))
blupValues.df=mutate(blupValues.df, ind=tolower(ind3))
blupValues.df=select(blupValues.df, -ind3)
blupValues.df

### Importing External Information (Lesion Length Values) ###
lesionValue.df<-read.csv("lesion.csv")
lesionValue.df=mutate(lesionValue.df, IND=gsub("[(]", "",lesionValue.df$Individuals))
lesionValue.df=mutate(lesionValue.df, IND2=gsub("[)]", ".",lesionValue.df$IND))
lesionValue.df<-select(lesionValue.df,Lesion.Length, IND2)
lesionValue.df=mutate(lesionValue.df, ind=tolower(IND2)) %>%select(-IND2)

### Combine external traits into data frame ###
traits.df.tmp<-merge(blupValues.df, normValues.df)
traits.df<-merge(traits.df.tmp, lesionValue.df)
rownames(traits.df) <- traits.df$ind
traits.df$ind <- NULL


### Network Analysis ###
  
  print(paste(".....Starting Round.... ", i))
 
  
### Normalizing raw counts using EdgeR ###

  print(paste(".....EdgeR.... ", i))
  counts<-WGCNA_NC_data_inoc
  length(colnames(counts))
  print(as.data.frame(colSums((counts))))
  group<-c("R", "S",rep("R", 47), rep("S", 46))
  print(head(counts))
  cds<-DGEList(counts, group=group)
  
  cds$samples
  ### Filtering cpm >1 across 50 samples ### 
  keep <- rowSums(cpm(cds)>1) >= 50
  cds <- cds[keep, , keep.lib.sizes=FALSE]
  cds <- calcNormFactors( cds)
  cds$samples
  
  plotMDS.DGEList( cds , main = "MDS Plot for Count Data",labels = colnames( cds$counts ) )
  
  logcpm<-cpm(cds, prior.count = 0.25, log=T)
  logcpm
  
  ### Caluclating total normalized filtered count per Indivdual resistance and suceptible separatly ###
  colsums_inoc<-as.data.frame(colSums(logcpm))
  colsums_inoc$Sample<-rownames(colsums_inoc)
  colnames(colsums_inoc)<-c("Sum", "Sample")
  inoc_R_rownames<-grep("^iR", colsums_inoc$Sample )
  inoc_S_rownames<-grep("^iS", colsums_inoc$Sample)
  inoc_Resistance<-colsums_inoc[inoc_R_rownames,]
  inoc_Susceptible<-colsums_inoc[inoc_S_rownames,]
  
  #### Picking same Indivduals from i=10, 6, 3, or 9 ###
  if(i == 10) {
    set.seed(1001*10)
  }
  else if (i == 6){
    set.seed(1001*6)
    
  }else if (i == 3){
    set.seed(1001*3)
    
  }else if (i == 9){
    set.seed(1001*9)
  }

  ### Calucaltuate minimum standard deviation of 30 random samples of resistant individuals and sucptible separartly ###
  sd_minimizer<-function(data, cutoff=1, rep=2000){
    sample_holder<-data.frame(Rep=as.numeric(), sd=as.numeric())
    for (i in seq_len(rep)) {
      r_sample_no<-sample(nrow(data), 30)
      pick_15_sample<-data[r_sample_no,]
      sds<-sd(pick_15_sample[,-2], na.rm = T)
      
      if (is.na(sds)){ 
        next
      } else{
        sample_holder<-rbind(sample_holder, data.frame(i, sds))
        
        
      }
    }
    
    
    sample_holder
  }
  
  
  
  lowest_sd_sample_picker<-function(data, cutoff=1, rep=2000){
    for (i in seq_len(rep)) {
      r_sample_no<-sample(nrow(data), 30)
      pick_15_sample<-data[r_sample_no,]
      sds<-sd(pick_15_sample[,-2], na.rm = T)
      
      if (is.na(sds)){ 
        next
      } else if (sds<=cutoff){
        sample_holder<-pick_15_sample
        return(sample_holder)
        print("Done")
        
      }
    }
    
  }
  
  inoc_Sus_outlier_rmd<-sd_minimizer(inoc_Susceptible,20000, 10000)
  inoc_lowest_Sus_sd<-range(inoc_Sus_outlier_rmd$sds)[1]
  inoc_Res_outlier_rmd<-sd_minimizer(inoc_Resistance,20000, 10000)
  inoc_lowest_Res_sd<-range(inoc_Res_outlier_rmd$sds)[1]
  
  
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
  
### Normailzing raw counts of selected 60 Individuals using EdgeR ### 
  fun_logcpm<-function(counts, grp, name){
    length(colnames(counts))
    as.data.frame(colSums((counts)))
    group<-grp
    cds<-DGEList(counts, group=group)
    cds$samples
    
 ### Filtering normalized counts by cpm > 1 acorss 25 samples ###      
    keep <- rowSums(cpm(cds)>1) >= 25
    cds <- cds[keep, , keep.lib.sizes=FALSE]
    cds <- calcNormFactors( cds)
    cds$samples
    plotMDS.DGEList( cds , main = "MDS Plot for Count Data",labels = colnames( cds$counts ) )
    logcpm<-cpm(cds, prior.count = .25, log=T)
    return(logcpm)
    
  }
  
  counts_S<-WGCNA_NC_data_inoc_S
  grp_S<-rep("S", 30)
  counts_R<-WGCNA_NC_data_inoc_R
  grp_R<-rep("R", 30)
  
  logcpm_S<-fun_logcpm(counts_S,grp_S, "SUS")
  logcpm_R<-fun_logcpm(counts_R,grp_R, "RES")
  
  logcpm.df.S<-as.data.frame(logcpm_S)
  if ( "iS280" %in% names(logcpm.df.S)){
  logcpm.df.S=select(logcpm.df.S, -iS280)}
  
  logcpm.df.R<-as.data.frame(logcpm_R)
  
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
  ### Counts number of rows in Susceptible samples ###
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

return(list(datTraits, datExpr) )
}

### WGCNA Network Analysis ###
### Cluster samples by expression ###
NETWORKTYPE<-"signed hybrid"
RSQUARED_CUTOFF<-0.9
thresholdZ.k = -0.8

#### Saranga fix me!!! ####
j=1

#### Selecting postitivley correlated networks 10 and 9  ####
return_l_boostra_i10<-sample_boostraping(i=5, urlBlup, urlNorm, WGCNA_NC_data_inoc)
return_l_boostra_i9<-sample_boostraping(i=3, urlBlup, urlNorm, WGCNA_NC_data_inoc)

datTraits.F_10<-return_l_boostra_i10[[1]]
datExpr.F_10<-return_l_boostra_i10[[2]]


datTraits.F_9<-return_l_boostra_i9[[1]]
datExpr.F_9<-return_l_boostra_i9[[2]]

traitData<-unique(do.call("rbind", list(datTraits.F_10,datTraits.F_9 )))

### Get common genes in both networks ###
common_genes<-intersect( names(datExpr.F_10), names(datExpr.F_9))
datExpr.F_9=datExpr.F_9 %>% select(which(common_genes %in% names(datExpr.F_9) ))
datExpr.F_10=datExpr.F_10 %>% select(which(common_genes %in% names(datExpr.F_10) ))

### Selected two networks ###
nSets = 2;

### For easier labeling of plots, create a vector holding descriptive names of the two networks ###
setLabels = c("i10 samples", "i9 samples")
shortLabels = c("i10", "i9")

### Form multi-set expression data: columns starting from 9 contain actual expression data. ###
multiExpr = vector(mode = "list", length = nSets)
multiExpr[[1]] = list(data = datExpr.F_10);
names(multiExpr[[1]]$data) = names(datExpr.F_10);
rownames(multiExpr[[1]]$data) = rownames(datExpr.F_10);
multiExpr[[2]] = list(data = datExpr.F_9);
names(multiExpr[[2]]$data) = names(datExpr.F_9);
rownames(multiExpr[[2]]$data) = rownames(datExpr.F_9);
exprSize = checkSets(multiExpr)

### Check that all genes and samples have sufficiently low numbers of missing values. ###
gsg = goodSamplesGenesMS(multiExpr, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  ### Print information about the removed genes ###
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(multiExpr[[1]]$data)[!gsg$goodGenes],
                                              # Check that all genes and samples have sufficiently low numbers of missing values.
                                              gsg = goodSamplesGenesMS(multiExpr, verbose = 3),
                                              gsg$allOK,
                                              collapse = ",")))
  for (set in 1:exprSize$nSets)
  {
    if (sum(!gsg$goodSamples[[set]]))
      printFlush(paste("In set", setLabels[set], "removing samples",
                       paste(rownames(multiExpr[[set]]$data)[!gsg$goodSamples[[set]]], collapse = ", ")))
    ### Remove the offending genes and samples ###
    multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes];
  }
  ### Update exprSize ##3
  exprSize=checkSets(multiExpr)
}

### Cluster the samples on their Euclidean distance, separately in each set. ##

sampleTrees = list()
for (set in 1:nSets){
  sampleTrees[[set]] = hclust(dist(multiExpr[[set]]$data), method = "average")
}
if (!dir.exists("Consensus")){
  dir.create("Consensus")
}

pdf(file = paste0("Consensus/SampleClustering",".pdf"), width = 12, height = 12);

par(mfrow=c(2,1))
par(mar = c(0, 4, 2, 0))
for (set in 1:nSets){
  plot(sampleTrees[[set]], main = paste("Sample clustering on all genes in", setLabels[set]),
       xlab="", sub="", cex = 0.7);}

dev.off()

#### Network anaylisis for i10 and i9 ###
#### Form a multi-set structure that will hold the external traits. ###
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets){
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, rownames(traitData));
  Traits[[set]] = list(data = traitData[traitRows,]);
  rownames(Traits[[set]]$data) = rownames(traitData[traitRows,]);
}

par(mar= c(6, 8.5, 3, 3))
for (set in 1:nSets){
  ### Calculate the whole network connectivity ##
  A = adjacency(t(multiExpr[[set]]$data),type="signed") 
  k = as.numeric(apply(A,2,sum))-1 # standardized connectivity
  Z.k = scale(k)
  outlierColor = ifelse(Z.k<thresholdZ.k,"red","black")
  sampleTree = flashClust(as.dist(1-A), method = "average")
  #### Convert traits to a color representation where red indicates high values ###
  traitColors = data.frame(numbers2colors(Traits[[set]]$data,signed=FALSE))
  dimnames(traitColors)[[2]] = paste(names(Traits[[set]]$data))
  datColors = data.frame(outlier = outlierColor,traitColors)
  
  plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                      colors=datColors,main="Sample Dendrogram and Trait Heatmap")
  
  remove.samples = Z.k<thresholdZ.k | is.na(Z.k)
  multiExpr[[set]]$data = multiExpr[[set]]$data[!remove.samples, ]
  
}

collectGarbage();
### Check the size of the leftover data ###
exprSize = checkSets(multiExpr)
exprSize
dev.off()

### Network Analysis for consensus of i10 and i9 ###
#### Form a multi-set structure that will hold the external traits.###
Traits = vector(mode="list", length = nSets);
for (set in 1:nSets){
  setSamples = rownames(multiExpr[[set]]$data);
  traitRows = match(setSamples, rownames(traitData));
  Traits[[set]] = list(data = traitData[traitRows,]);
  rownames(Traits[[set]]$data) = rownames(traitData[traitRows,]);
}
collectGarbage();
### Define data set dimensions ###
nGenes = exprSize$nGenes;
nSamples = exprSize$nSamples

save(multiExpr, Traits, nGenes, nSamples, setLabels, shortLabels, exprSize,file = paste0("Consensus/Consensus-dataInput",".RData"))

### 2.a Step-by-step network construction and module detection ##3
### Choose a set of soft-thresholding powers ###
powers = c(seq(4,10,by=1), seq(12,20, by=2));
### Initialize a list to hold the results of scale-free analysis ##
powerTables = vector(mode = "list", length = nSets);
### Call the network topology analysis function for each set in turn ###
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,verbose = 2, networkType=NETWORKTYPE, corFnc= "bicor",corOptions=list(use = 'p',maxPOutliers=0.1))[[2]]);
collectGarbage();
### Plot the results ###
colors = c("black", "red", "green")
### Plot these columns of the returned scale free analysis tables ###
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
### Get the minima and maxima of the plotted points ##
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
### Plot the quantities in the chosen columns vs. the soft thresholding power ###
sizeGrWindow(8, 6)
pdf(file = paste0("Consensus/SoftThreshold",j,".pdf"), width = 12, height = 12);
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1) {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid(); }
  if (col==1) {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off()

### 2.a.2 Calculation of network adjacencies ###
softPower = 4;

### Initialize an appropriate array to hold the adjacencies ###
adjacencies = array(0, dim = c(nSets, nGenes, nGenes));

### Calculate adjacencies in each individual data set ###
for (set in 1:nSets){
  adjacencies[set, , ] = adjacency(multiExpr[[set]]$data, power = softPower, type=NETWORKTYPE)
  }

### 2.a.3 Calculation of Topological Overlap ###

### Initialize an appropriate array to hold the TOMs ###
TOM = array(0, dim = c(nSets, nGenes, nGenes));

### Calculate TOMs in each individual data set ###
for (set in 1:nSets){
  TOM[set, , ] = TOMsimilarity(adjacencies[set, , ], TOMType="signed", verbose=1);}

# 2.a.4 Scaling of Topological Overlap Matrices to make them comparable across sets
### Define the reference percentile ###
scaleP = 0.95

### Sample sufficiently large number of TOM entries ###
nSamples = as.integer(1/(1-scaleP) * 1000);

### Choose the sampled TOM entries ##
scaleSample = sample(nGenes*(nGenes-1)/2, size = nSamples)
TOMScalingSamples = list();

### These are TOM values at reference percentile ###
scaleQuant = rep(1, nSets)

### Scaling powers to equalize reference TOM values ###
scalePowers = rep(1, nSets)

### Loop over sets ###
for (set in 1:nSets)
{
  ### Select the sampled TOM entries ##
  TOMScalingSamples[[set]] = as.dist(TOM[set, , ])[scaleSample]
  ### Calculate the 95th percentile ###
  scaleQuant[set] = quantile(TOMScalingSamples[[set]],
                             probs = scaleP, type = 8);
  ### Scale the male TOM ###
  if (set>1) {
    scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set]);
    TOM[set, ,] = TOM[set, ,]^scalePowers[set];
  }
}

### For plotting, also scale the sampled TOM entries ###
scaledTOMSamples = list();
for (set in 1:nSets)
  scaledTOMSamples[[set]] = TOMScalingSamples[[set]]^scalePowers[set]

### Open a suitably sized graphics window ###
sizeGrWindow(6,6)
pdf(file = paste0("Consensus/TOMScaling-QQPlot",j,".pdf"), width = 6, height = 6);

### qq plot of the unscaled samples ###
qqUnscaled = qqplot(TOMScalingSamples[[1]], TOMScalingSamples[[2]], plot.it = TRUE, cex = 0.6,
                    xlab = paste("TOM in", setLabels[1]), ylab = paste("TOM in", setLabels[2]),main = "Q-Q plot of TOM", pch = 20)

### qq plot of the scaled samples ###
qqScaled = qqplot(scaledTOMSamples[[1]], scaledTOMSamples[[2]], plot.it = FALSE)
points(qqScaled$x, qqScaled$y, col = "red", cex = 0.6, pch = 20);
abline(a=0, b=1, col = "blue")
legend("topleft", legend = c("Unscaled TOM", "Scaled TOM"), pch = 20, col = c("black", "red"))
dev.off()

###2.a.5 Calculation of consensus Topological Overlap ###
consensusTOM = pmin(TOM[1, , ], TOM[2, , ]);

### 2.a.6 Clustering and module identification ###

### Clustering ###
consTree=flashClust(as.dist(1-consensusTOM), method = "average");

### Trying to obtain large modules, so we set the minimum module size relatively high: ###
minModuleSize = 30;
DEEP_SPLIT=0;

### Module identification using dynamic tree cut ###
unmergedLabels = cutreeDynamic(dendro = consTree, distM = 1-consensusTOM,
                               deepSplit = DEEP_SPLIT, cutHeight = 0.995,
                               minClusterSize = minModuleSize,
                               pamRespectsDendro = FALSE );
unmergedColors = labels2colors(unmergedLabels)
sizeGrWindow(8,6)
pdf(file = paste0("Consensus/TreeConstruction",j,".pdf"), width = 8, height = 6);
plotDendroAndColors(consTree, unmergedColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

###  2.a.7 Merging of modules whose expression profiles are very similar ##
### Calculate module eigengenes ###
unmergedMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = unmergedColors)

### Calculate consensus dissimilarity of consensus module eigengenes ###
consMEDiss = consensusMEDissimilarity(unmergedMEs);

### Cluster consensus modules ###
consMETree = hclust(as.dist(consMEDiss), method = "average");

### Plot the result ###
sizeGrWindow(7,6)
par(mfrow = c(1,1))
plot(consMETree, main = "Consensus clustering of consensus module eigengenes",
     xlab = "", sub = "")
abline(h=0.8, col = "red")
merge = mergeCloseModules(multiExpr, unmergedLabels, cutHeight = 0.8, verbose = 3)

### Numeric module labels ###
moduleLabels = merge$colors;

### Convert labels to colors ##
moduleColors = labels2colors(moduleLabels)

### Eigengenes of the new merged modules ###
consMEs = merge$newMEs
sizeGrWindow(9,6)
plotDendroAndColors(consTree, cbind(unmergedColors, moduleColors),
                    c("Unmerged", "Merged"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
save(consMEs, moduleColors, moduleLabels, consTree, file = paste0("Consensus/Consensus-NetworkConstruction-man",j,".RData"))

### Set up variables to contain the module-trait correlations ###
moduleTraitCor = list();
moduleTraitPvalue = list();

### Calculate the correlations ###
for (set in 1:nSets)
{
  moduleTraitCor[[set]] = cor(consMEs[[set]]$data, Traits[[set]]$data, use = "p");
  moduleTraitPvalue[[set]] = corPvalueFisher(moduleTraitCor[[set]], exprSize$nSamples[set]);
}

### Convert numerical lables to colors for labeling of modules in the plot ##
MEColors = labels2colors(as.numeric(substring(names(consMEs[[1]]$data), 3)));
MEColorNames = paste("ME", MEColors, sep="");

### Open a suitably sized window ##
sizeGrWindow(10,7)

### Plot the module-trait relationship table for Network Analysis 10 ###
pdf(file = paste0("Consensus/ModuleTraitRelationships","-i10.pdf"), wi = 7, he = 10);
set = 1
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
par(mar = c(6, 8.8, 3, 2.2));
labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off();

### Plot the module-trait relationship table for network analysis 9 ###

set = 2
textMatrix = paste(signif(moduleTraitCor[[set]], 2), "\n(",
                   signif(moduleTraitPvalue[[set]], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,7)

pdf(file = paste0("Consensus/ModuleTraitRelationships","i3.pdf"), wi = 7, he = 10);
par(mar = c(6, 8.8, 3, 2.2));

labeledHeatmap(Matrix = moduleTraitCor[[set]],
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module--trait relationships in", setLabels[set]))
dev.off()



### Initialize matrices to hold the consensus correlation and p-value ##
consensusCor = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));
consensusPvalue = matrix(NA, nrow(moduleTraitCor[[1]]), ncol(moduleTraitCor[[1]]));

### Find consensus negative correlations ###
negative = moduleTraitCor[[1]] < 0 & moduleTraitCor[[2]] < 0;
consensusCor[negative] = pmax(moduleTraitCor[[1]][negative], moduleTraitCor[[2]][negative]);
consensusPvalue[negative] = pmax(moduleTraitPvalue[[1]][negative], moduleTraitPvalue[[2]][negative]);

### Find consensus positive correlations ###
positive = moduleTraitCor[[1]] > 0 & moduleTraitCor[[2]] > 0;
consensusCor[positive] = pmin(moduleTraitCor[[1]][positive], moduleTraitCor[[2]][positive]);
consensusPvalue[positive] = pmax(moduleTraitPvalue[[1]][positive], moduleTraitPvalue[[2]][positive]);

textMatrix = paste(signif(consensusCor, 2), "\n(",
                   signif(consensusPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[[set]])
sizeGrWindow(10,15)
pdf(file = paste0("Consensus/ModuleTraitRelationships","-consensus.pdf"), wi = 12, he = 20);
par(mar = c(8, 8.8, 3, 2.2));
labeledHeatmap(Matrix = consensusCor,
               xLabels = names(Traits[[set]]$data),
               yLabels = MEColorNames,
               ySymbols = MEColorNames,
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Consensus module--trait relationships across\n",
                            paste(setLabels, collapse = " and ")))

dev.off()

### 4.a Exporting results of the Consensus network analysis ####
annotation<-fread("Gmax_275_Wm82.a2.v1.annotation_merged.txt", header = T)
probes = names(multiExpr[[1]]$data)
probes2annot = match(probes, annotation$locusName)
consMEs.unord = multiSetMEs(multiExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}

GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2);
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);

GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 6*nTraits)
rownames(GSmat) = probes;
colnames(GSmat) = spaste(
  c("GS.i10.", "GS.i9.", "p.GS.i10.", "p.GS.i9.", "Z.GS.meta.", "p.GS.meta"),
  rep(traitNames, rep(6, nTraits)))

### Same code for kME ##
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
  c("kME.i10.", "kME.i9.", "p.kME.i10.", "p.kME.i9.", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(6, nMEs)))

info = data.frame(Probe = probes, pacid = annotation$pacId[probes2annot],
                  peptideName = annotation$peptideName[probes2annot],
                  pfam = annotation$Pfam[probes2annot],
                  panther = annotation$Panther[probes2annot],
                  kog = annotation$KOG[probes2annot],
                  kegg_ec = annotation$`KEGG/ec`[probes2annot],
                  go = annotation$GO[probes2annot],
                  v2 = annotation$V2[probes2annot],
                  v3 = annotation$V3[probes2annot],
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  GSmat,
                  kMEmat);
write.table(info, file = paste0("Consensus/consensusAnalysis-CombinedNetworkResults",".txt"),row.names = FALSE, quote = FALSE,  sep="\t");




### 5 Comparing eigengene networks Lesion length ###
### Create a variable lesion length values ###
lesion.length = vector(mode = "list", length = nSets);
for (set in 1:nSets)
{
  lesion.length[[set]] = list(data = as.data.frame(Traits[[set]]$data$Lesion.Length));
  names(lesion.length[[set]]$data) = "LesionLength"
  rownames(lesion.length[[set]]$data) = rownames(Traits[[set]]$data)
}

### Recalculate consMEs to give them color names ###
consMEsC = multiSetMEs(multiExpr, universalColors = moduleColors);

### We add the Lesion length value to the eigengenes and order them by consesus hierarchical clustering ###
MET = consensusOrderMEs(addTraitToMEs(consMEsC, lesion.length));

sizeGrWindow(15,10);
pdf(file = paste0("Consensus/EigengeneNetworks",".pdf"), width= 10, height = 12);
par(cex = 0.8)
plotEigengeneNetworks(MET, setLabels, marDendro = c(3,5,5,3), marHeatmap = c(3,5,5,3),
                      zlimPreservation = c(0.5, 1), xLabelsAngle = 0)
dev.off()

eign_lesion<-addTraitToMEs(consMEsC, lesion.length)

### Writting consensus module egigengene values to lesion length comparison ###
write.csv(eign_lesion[[1]]$data, file = paste0("Consensus/Consensus_eign_lesion",".csv"))

#### Plotting Resulting ###
library("ggplot2")
plot_data<-as.data.frame(eign_lesion[[1]]$data)
ggplot(plot_data, aes(x=LesionLength, y= MEthistle1)) + geom_area( fill="pink", alpha=.2)+ geom_line()

ggplot(plot_data, aes(x=LesionLength, y=MEdarkgreen)) + geom_area( fill="green", alpha=.2)+geom_line()


dev.off()


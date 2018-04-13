### Loading Libraries ###
library("data.table")
library('dplyr')

#### Loading Inoculated PPA hotspot file ###
PPA_Hotspots_file<-"PPA_hotspots_transposed.csv"
ppa_table<-read.csv(PPA_Hotspots_file, na.strings = "N/A", header=TRUE)

### Signficant Postivily correlated WGCNA modules ###
## Filter experssion values for genes in selected module ###
datExprCons<-info %>% filter("honeydew1" ==  ModuleColor) %>% select(Probe, GS.i10.normalized.total.counts, GS.i9.normalized.total.counts)
names(datExprCons)<-c("Probe", "Reference2", "Reference1")
datExprModColor<- datExprCons
rownames(ppa_table)<-ppa_table$GI
ppa_table$GI<-NULL
temp<-as.data.frame(t(ppa_table))
temp['GI']<-rownames(temp)
datExprModColor['GI']<-datExprModColor$Probe
datExprModColor_PPA<-datExprModColor %>% filter(GI %in% temp$GI)
ppa_table_matched_ModColor<-temp %>% filter(GI %in% datExprModColor$GI)
rownames(datExprModColor_PPA)<-unlist(datExprModColor_PPA$GI)
datExprModColor_PPA$GI<-NULL
datExprModColor_PPA$Probe<-NULL
rownames(ppa_table_matched_ModColor)<-unlist(ppa_table_matched_ModColor$GI)
ppa_table_matched_ModColor$GI<-NULL
ppa_table_matched_ModColor_t<-as.matrix(t(ppa_table_matched_ModColor))

# This defines the Hotspots based gene significance measure ###
GSSNP=data.frame(matrix(NA, nrow=dim(ppa_table_matched_ModColor_t)[[1]],ncol=dim(datExprModColor_PPA)[[2]] ))
for (i in c(1:dim(ppa_table_matched_ModColor_t)[[1]]) ){GSSNP[i,]= as.numeric(abs(cor(as.numeric(ppa_table_matched_ModColor_t[i,]), datExprModColor_PPA,use="p")))}
row.names(GSSNP)=paste(unlist(colnames(ppa_table_matched_ModColor)))
GSSNP_P_values<-corPvalueStudent(as.matrix(GSSNP),28)
textMatrixSNP= paste(signif(as.matrix(GSSNP), 2), "\n(", 
                  signif(GSSNP_P_values, 1), ")", sep= "")
dim(textMatrixSNP)= dim(GSSNP)
pdf( "Consensus/Heatmaps_Eigengen_Expr_Cor_PPA_honeydew1.pdf", width=6, height=15)
par(mar= c(6, 12,3,3))

# display the corelation values with a heatmap plot
labeledHeatmap(Matrix= GSSNP, 
               xLabels= names(datExprModColor_PPA), 
               yLabels= rownames(GSSNP), 
               ySymbols= rownames(GSSNP), 
               colorLabels= FALSE, 
               colors= greenWhiteRed(50), 
               textMatrix= textMatrixSNP, 
               setStdMargins= FALSE, 
               cex.text= 0.9, 
               zlim= c(-1,1), 
               main= paste("Honeydew1-S-PPA relationships"))
dev.off()


## LOAD LIBRARIES PACKAGES
library(gplots)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(kableExtra)
library(knitr)
library(stringr)
library(gridExtra)
library(grid)
library(lattice)
library(gtable)
library(dplyr)
library(xml2)
library(ggrepel) #geom_text_repel or geom_label_repel
library(xtable)

## FOR DIFFERENTIAL EXPRESSION & ANALYSIS
library(DESeq2)
library(enrichR)
library(UpSetR)
library(biomaRt)

## PACKAGES FOR PATHWAY ANALYSIS
library(PoiClaClu)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pathview)
library(gage)
library(gageData)

## ## POWER CALCULATION
library(PROPER)

## OPTIONS
options(repr.matrix.max.rows = 100)
columns(org.Hs.eg.db)

## SET WORKING DIRECTORY (ASSUMING YOU HAVE EVERYTHING SETUP IN THE FOLDER YOU WORKING IN)
setwd(getwd())

## FOR HTSEQ ANALYSIS - COMMENT BELOW
htseq.counts.dir <- paste0(getwd(),"/results_htseq_counts")
if (file.exists(htseq.counts.dir)){
    in.file <- "nf-rnaSeqCount-scleroderma/htseqCounts/gene_counts_final.txt"
    out.dir <- paste0(htseq.counts.dir,'/')
}else {
    dir.create(file.path(htseq.counts.dir), showWarnings = FALSE)
    in.file <- "nf-rnaSeqCount-scleroderma/htseqCounts/gene_counts_final.txt"
    out.dir <- paste0(htseq.counts.dir,'/')
}

## ## FOR FEATURESEQ ANALYSIS - COMMENT ABOVE
## feature.counts.dir <- paste0(getwd(),"/results_feature_counts")
## if (file.exists(feature.counts.dir)){
##     in.file <- "nf-rnaSeqCount-scleroderma/featureCounts/gene_counts_final.txt"
##     out.dir <- paste0(feature.counts.dir,'/')
## }else {
##     dir.create(file.path(feature.counts.dir), showWarnings = FALSE)
##     in.file <- "nf-rnaSeqCount-scleroderma/featureCounts/gene_counts_final.txt"
##     out.dir <- paste0(feature.counts.dir,'/')
## }


## ------------------------------ LOAD COUNTS MATRIX DATA & ASSIGN NAMES, CONDITION, SITE AND GENDER ETC.
the.data <- read.csv(file=in.file, header=TRUE, row.names=1, sep = '\t')
the.data <- the.data[,order(names(the.data))]
the.data <- the.data[!row.names(the.data) %in% c("__no_feature", "__ambiguous", "__too_low_aQual","__not_aligned","__alignment_not_unique"),]

## ASSIGN SAMPLE NAMES
new_names <-  c('P1', 'B1', 'P2', 'B2', 'P3', 'B3', 'P4', 'B4', 'P5', 'B5', 'P6', 'B6', 'P7', 'B7', 'P8', 'B8', 'P9', 'B9',
                'C1', 'C2', 'R3', 'R5', 'C6', 'C7', 'R8')
colnames(the.data) <- new_names
the.data <- as.matrix(the.data)

## SAMPLE CONRITION
condition <- factor(c(rep('case',18),rep('ctrl',7)))

## SAMPLE SITE
site <- factor(c(rep(c('arm', 'back'),9),
                 'arm', 'arm', 'back', 'back', 'arm', 'arm', 'back'))

## GENDER
gender <- factor(c(rep('F',16), rep('M',2),
                   'F', 'M', 'F', 'F', 'M', 'M', 'F'))

## AGE
age <- c(rep('48',2), rep('54',2), rep('48',2), rep('52',2), rep('49',2), rep('55',2), rep('58',2), rep('46',2), rep('64',2),
               '29', '33', '40', '38', '31', '34', '30')

## DISEASE DURATION
disease.duration <- c(rep('0',2), rep('36',2), rep('36',2), rep('60',2), rep('36',2), rep('36',2), rep('24',2), rep('12',2), rep('12',2),
                      rep('0',7))

## SKIN SCORE (MRSS)
mrss <- c(rep('0',2), rep('34',2), rep('4',2), rep('9',2), rep('16',2), rep('32',2), rep('5',2), rep('25',2), rep('22',2),
          rep('0',7))

## AFFECTED SKIN SCORE
affected.mrss <- c(rep('0',2), '3', '0', '1', '0', '2', '0', '2', '0', '3', '0', '1', '0', '2', '0', '1', '0',
                   rep('0',7))

## DISEASE ACTIVITY SCORE
disease.activity <- c(rep('0',2), rep('2.5',2), rep('3',2), rep('3',2), rep('4',2), rep('6.5',2), rep('4.5',2), rep('4',2), rep('7.7',2),
                            rep('0',7))

## ERYTHROCYTE SEDIMENTATION RATE
esr <- c(rep('0',2), rep('0',2), rep('1',2), rep('1',2), rep('1',2), rep('1',2), rep('0',2), rep('1',2), rep('1',2),
         rep('0',7))

## C-REACTIVE PROTEIN (CRP)
crp <- c(rep('0',2), rep('5',2), rep('8.8',2), rep('12.5',2), rep('25',2), rep('13',2), rep('7',2), rep('21',2), rep('8',2),
         rep('0',7))

## LUNG DISEASE
lung <- factor(c(rep('0',2), rep('1',2), rep('0',2), rep('0',2), rep('1',2), rep('1',2), rep('0',2), rep('1',2), rep('0',2),
                 rep('0',7)))

## MYOSITIS
myositis <- factor(c(rep('0',2), rep('1',2), rep('0',2), rep('1',2), rep('0',2), rep('1',2), rep('0',2), rep('0',2), rep('0',2),
              rep('0',7)))

## DIFFUSING CAPACITY OF THE LUNG FOR CARBON MONOXIDE (DLCO)
dlco <- c(rep('0',2), rep('0.62',2), rep('0.64',2), rep('1.09',2), rep('0.45',2), rep('0.50',2), rep('0.83',2), rep('0.40',2), rep('0.65',2),
          rep('0',7))

## ANTI-SCL 70
anti.scl <- factor(c(rep('0',2), rep('1',2), rep('0',2), rep('0',2), rep('0',2), rep('1',2), rep('0',2), rep('1',2), rep('0',2),
                     rep('0',7)))

## ASSIGN - CREATE THE COLUMN DATA (DESCRIPTIVE/CLINICAL INFORMATION FOR SAMPLES)
the.coldata <- DataFrame(row.names=colnames(the.data), condition, site, gender, age, disease.duration, lung, myositis, anti.scl, mrss, affected.mrss, disease.activity,
                     esr, crp, dlco)


## FILTERING TABLE - POPULATED AT EACH FILTERING STAGE
filtering.table <- data.frame(genes=rep(length(row.names(the.data)),9),
                              row.names=c("ALL","CASE.ARM","CASE.BACK","CASE",
                                          "SEVERE.ARM","SEVERE.BACK",
                                          "MILD.ARM","MILD.BACK","SEVERE.MILD"))

## REMOVE SAMPLES P9, B9 AND C1. THEY ARE WOLVES AMONGST SHEEP
the.data <- the.data[, !colnames(the.data) %in% c('P1', 'B1', 'P9', 'B9', 'C1')]
the.data <- the.data[,order(colnames(the.data))]

the.data <- the.data[which(rowSums(the.data) >= 1),]
filtering.table$step.1 <- rep(length(row.names(the.data)),9)

the.data <- the.data[which(rowSums(the.data >= 5) >= 3),]
filtering.table$step.2 <- rep(length(row.names(the.data)),9)

the.coldata <- the.coldata[!rownames(the.coldata) %in% c('P1', 'B1', 'P9', 'B9', 'C1'),]
the.coldata <- the.coldata[order(rownames(the.coldata)),]

## ASSIGN INDIVIDUALS NUMBERS (1-13 INDIVIDUALS IN THESE ANALYSIS)
the.coldata$ind <- factor(c(seq(1,7,1), 8, 11, 12, seq(1,7,1), 9, 10, 13))

## ASSIGN INDIVIDUALS IN NESTED WITHIN GROUPS (CONDITION, LUNG DISEASE, MYOSITIS, ANTI-SCL) - FOR DESIGN FORMULA
the.coldata$ind.n <- factor(c(seq(1,7,1), seq(1,3,1), seq(1,7,1), seq(4,6,1)))
the.coldata$lung.n <- factor(c(1, 1, 2, 2, 3, 3, 4, rep("NA",3), 1, 1, 2, 2, 3, 3, 4, rep("NA",3)))
the.coldata$myositis.n <- factor(c(1, 1, 2, 2, 3, 3, 4, rep("NA",3), 1, 1, 2, 2, 3, 3, 4, rep("NA",3)))
the.coldata$anti.scl.n <- factor(c(1, 1, 2, 3, 2, 4, 3, rep("NA",3), 1, 1, 2, 3, 2, 4, 3, rep("NA",3)))

## LOAD ENSEMBL DATASETS FOR HUMAN
ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

## LOAD FILE WITH THE FUNCTIONS USED IN THIS SCRIPT
source('scleroderma_functions.R')

## GET ENSEMBLE ANNOTATIONS
ensembl.annotations <- get_ensembl_annotations(row.names(the.data))
ensembl.annotations <- ensembl.annotations[which(!ensembl.annotations$HGNC %in% c("C12orf74","LINC00595","CCL3L3")),]

## REMOVE X & Y CHROMOSOME GENES FROM THE DATA
xy.chrom <- ensembl.annotations[which(ensembl.annotations$chromosome == 'Y'),]$ensembl_id
the.data <- the.data[which(!row.names(the.data) %in% xy.chrom),]
filtering.table$step.3 <- rep(length(row.names(the.data)),9)

## ============================= DIFFERENTIAL EXPRESSION ANALYSIS STARTS HERE! ============================= ##
## CREATE DATA GROUPS
case_arm <- c("P2","P3","P4","P5","P6","P7","P8")
case_back <- c("B2","B3","B4","B5","B6","B7","B8")
ctrl_arm <- c("C2","C6","C7")
ctrl_back <- c("R3","R5","R8")

scl.case_arm <- c("P2","P6","P8")
scl.case_back <- c("B2","B6","B8")
scl.ctrl_arm <- c("P3","P4","P5","P7")
scl.ctrl_back <- c("B3","B4","B5","B7")

## COMPARE SUBSETS - LIST OF LISTS
compare.ALL <- list(c(case_arm, case_back), c(ctrl_arm, ctrl_back))
compare.CASE.ARM <- list(case_arm, c(ctrl_arm, ctrl_back))
compare.CASE.BACK <- list(case_back, c(ctrl_arm, ctrl_back))
compare.CASE <- list(case_arm, case_back)
compare.SEVERE.ARM <- list(scl.case_arm, c(ctrl_arm, ctrl_back))
compare.SEVERE.BACK <- list(scl.case_back, c(ctrl_arm, ctrl_back))
compare.MILD.ARM <- list(scl.ctrl_arm, c(ctrl_arm, ctrl_back))
compare.MILD.BACK <- list(scl.ctrl_back, c(ctrl_arm, ctrl_back))
compare.SEVERE.MILD <- list(c(scl.case_arm, scl.case_back), c(scl.ctrl_arm, scl.ctrl_back))

## PATH FOR TABLES
tables.out.dir <- paste0(out.dir, "tables/")
if (file.exists(tables.out.dir)){
    system(paste0('rm -r ','"',tables.out.dir,'"',"*"))
}else {
    dir.create(file.path(tables.out.dir), showWarnings = TRUE)
}

## PATH FOR FIGURES
figures.out.dir <- paste0(out.dir, "figures/")
if (file.exists(figures.out.dir)){
    system(paste0('rm -r ','"',figures.out.dir,'"',"*"))
}else {
    dir.create(file.path(figures.out.dir), showWarnings = TRUE)
}

## CREATE ENRICHR DIR
enrichr.out.dir <- paste0(out.dir, "enrichr/")
if (file.exists(enrichr.out.dir)){
    system(paste0('rm -r ','"',enrichr.out.dir,'"',"*"))
}else {
    dir.create(file.path(enrichr.out.dir), showWarnings = TRUE)
}

## CREATE OUTPUT DIRECTORY CALLED "PATHWAYS" INSIDE THE FOLDER OF OUTPUT
path.out.dir <- paste0(out.dir, "pathways/")
if (file.exists(path.out.dir)){
    system(paste0('rm -r ','"',path.out.dir,'"',"*"))
}else {
    dir.create(file.path(path.out.dir), showWarnings = TRUE)
}

## KEGG DOWNLOADED PATHWAYS
kegg.out.dir <- paste0(out.dir, "pathways/kegg_downloads/")
if (file.exists(kegg.out.dir)){
    system(paste0('rm -r ','"',kegg.out.dir,'"',"*"))
}else {
    dir.create(file.path(kegg.out.dir), showWarnings = TRUE)
}

## PROPER POWER CALCULATIONS
sim.opts.Cheung <- RNAseq.SimOptions.2grp(ngenes = 25000,
                                          p.DE=0.05, lOD="cheung",
                                          lBaselineExpr="cheung",
                                          lfc=2)

simres = runSims(Nreps = c(3, 4, 6, 7, 14,7),
                 Nreps2 = c(6, 6, 8, 6, 6,7),
                 sim.opts=sim.opts.Cheung,
                 DEmethod="DESeq2", nsims=20)

powers = comparePower(simres, alpha.type="fdr", alpha.nominal=0.05,
                      stratify.by="expr", strata.filtered=1,
                      target.by = 'lfc')
summaryPower(powers)


power.table <- as.data.frame(summaryPower(powers))
power.table$`Actual FDR` <- formatC(power.table$`Actual FDR`, format = "g", digits = 2)
power.table$`Marginal power` <- formatC(power.table$`Marginal power`, format = "g", digits = 2)
power.table$FDC <- formatC(power.table$FDR, format = "g", digits = 2)

power.table

out.table = paste0(tables.out.dir,"table.powers.tex")
table.head <- "\\begin{longtable}{ C{2em} C{2em} C{6.5em} C{6.5em} C{6.5em} C{6.5em} C{6.5em} C{2.2em} }"
table.caption <- paste0("\\caption[Summary of the power calculations performed with \\proper{}]{\\textbf{Summary of the power calculations performed with \\proper{}.}}")
table.label <- paste0("\\label{tab:summary.power}\\\\")

write.table(table.head, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" & ", eol=" \n")
write.table(table.caption, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
write.table(table.label, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")

write.table("\\toprule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
write.table(paste0("\\rowcolor{gray!25} ", paste(names(power.table),collapse=" & ")),
            file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
write.table("\\midrule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
write.table(power.table, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
write.table("\\bottomrule\n\\end{longtable}", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")

pdf(paste0(figures.out.dir,'plot.powers.pdf'), paper='a4r', onefile=FALSE)
plotPower(powers)
dev.off()

## ----------------------------- REMOVE GENDER/SITE BIAS ("SANITIZE")????? USE CONTROL MALE & FEMALES
sanitize <- the.data[, colnames(the.data) %in% c(ctrl_arm,ctrl_back)]
sanitize.coldata <- the.coldata[rownames(the.coldata) %in% c(ctrl_arm,ctrl_back),]

SANITIZE <- DESeqDataSetFromMatrix(countData=the.data, colData=the.coldata, design = ~ site)
SANITIZE$site <- relevel(SANITIZE$site, ref="back")
SANITIZE <- estimateSizeFactors(SANITIZE)
SANITIZE <- SANITIZE[ rowSums(counts(SANITIZE, normalized = TRUE) >= 5) > 3,]

dds.SANITIZE <- DESeq(SANITIZE)

res.SANITIZE <- results(dds.SANITIZE, alpha=0.05, lfcThreshold=1.1, altHypothesis="greaterAbs")
summary(res.SANITIZE)
sum(res.SANITIZE$padj < 0.05, na.rm=TRUE)

sig.SANITIZE <- get_sig(res.SANITIZE, 0.05, 1.1)
sanitize_genes <- rownames(sig.SANITIZE)

discarded_genes <- get_gene_names(sanitize_genes)
discarded_genes

the.data <- the.data[ !rownames(the.data) %in% sanitize_genes, ]

filtering.table$step.4 <- rep(length(row.names(the.data)),9)

## ------------------------------ GET ENSEMBL BIOTYPES BEFORE RUNNING ANY ANALYSES

## GET THE GENE LIST FROM THE DATA LOADED ABOVE
gene_list <- row.names(the.data)

## GET THE BIOTYPES USING THE ENSEMBL GENE IDS
biotypes <- get_gene_names(gene_list)[,c('ensembl_id','biotype')]

## GET THE FREQUENCY TABLE OF BIOTYPES
freq.table <- as.data.frame(table(biotypes$biotype))
freq.table <- freq.table[order(freq.table$Freq, decreasing=TRUE),]
colnames(freq.table) <- c("Gene Biotype","Number of Genes")

## PLOT BIOTYPES
pdf(paste0(figures.out.dir,'plot.biotypes_before.pdf'), paper='a4', onefile=FALSE)
ggplot(data=freq.table, aes(x=reorder(`Gene Biotype`, -`Number of Genes`), y=`Number of Genes`)) +
    geom_bar(stat="identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    xlab("Gene Biotypes") +
    ylab("Counts") +
    ggtitle("Categories of Genes")
dev.off()

## ------------------------------ CREATED THE DESEQ DATASETS FOR EACH COMPARISON

## 1. ALL CASES VS ALL CONTROLS
title.ALL <- 'ALL: Cases VS Controls'
ALL <- DESeqDataSetFromMatrix(countData=the.data, colData=the.coldata, design = ~ site + site:ind.n + site:condition)
ALL$condition <- relevel(ALL$condition, ref="ctrl")
ALL <- estimateSizeFactors(ALL)
ALL <- ALL[ rowSums(counts(ALL,normalized=TRUE) >= 5) >= 3, ]

## 2. case_arm vs all_ctrl
title.CASE.ARM <- 'CASE.ARM: Cases (Arm) VS Controls'
case.arm <- the.data[, colnames(the.data) %in% c(case_arm,ctrl_arm,ctrl_back)]
case.arm.coldata <- the.coldata[rownames(the.coldata) %in% c(case_arm,ctrl_arm,ctrl_back),]

CASE.ARM <- DESeqDataSetFromMatrix(countData=case.arm, colData=case.arm.coldata, design = ~ site + condition)
CASE.ARM$condition <- relevel(CASE.ARM$condition, ref="ctrl")
CASE.ARM <- estimateSizeFactors(CASE.ARM)
CASE.ARM <- CASE.ARM[ rowSums(counts(CASE.ARM,normalized=TRUE) >= 5) >= 3, ]

## 3. case_back vs all_ctrl
title.CASE.BACK <- 'CASE.BACK: Cases (Back) VS Controls'
case.back <- the.data[, colnames(the.data) %in% c(case_back,ctrl_arm,ctrl_back)]
case.back.coldata <- the.coldata[rownames(the.coldata) %in% c(case_back,ctrl_arm,ctrl_back),]

CASE.BACK <- DESeqDataSetFromMatrix(countData=case.back, colData=case.back.coldata, design = ~ site + condition)
CASE.BACK$condition <- relevel(CASE.BACK$condition, ref="ctrl")
CASE.BACK <- estimateSizeFactors(CASE.BACK)
CASE.BACK <- CASE.BACK[ rowSums(counts(CASE.BACK,normalized=TRUE) >= 5) >= 3, ]

## 4. case_arm vs case_back
title.CASE <- 'CASE: Cases (Arm) VS Case (Back)'
case <- the.data[, colnames(the.data) %in% c(case_arm,case_back)]
case.coldata <- the.coldata[rownames(the.coldata) %in% c(case_arm,case_back),]

CASE <- DESeqDataSetFromMatrix(countData=case, colData=case.coldata, design = ~ ind.n + site)
CASE$site <- relevel(CASE$site, ref="back")
CASE <- estimateSizeFactors(CASE)
CASE <- CASE[ rowSums(counts(CASE,normalized=TRUE) >= 5) >= 3, ]

## ## 5. scle_case_arm vs all_ctrl
title.SEVERE.ARM <- 'SEVERE.ARM: Severe (Arm) VS Controls'
severe.arm <- the.data[, colnames(the.data) %in% c(scl.case_arm,ctrl_arm,ctrl_back)]
severe.arm.coldata <- the.coldata[rownames(the.coldata) %in% c(scl.case_arm,ctrl_arm,ctrl_back),]

SEVERE.ARM <- DESeqDataSetFromMatrix(countData=severe.arm, colData=severe.arm.coldata, design = ~ site + condition)
SEVERE.ARM$condition <- relevel(SEVERE.ARM$condition, ref="ctrl")
SEVERE.ARM <- estimateSizeFactors(SEVERE.ARM)
SEVERE.ARM <- SEVERE.ARM[ rowSums(counts(SEVERE.ARM,normalized=TRUE) >= 5) >= 3, ]

## 6. scle_case_back vs all_ctrl
title.SEVERE.BACK <- 'SEVERE.BACK: Severe (Back) VS Controls'
severe.back <- the.data[, colnames(the.data) %in% c(scl.case_back,ctrl_arm,ctrl_back)]
severe.back.coldata <- the.coldata[rownames(the.coldata) %in% c(scl.case_back,ctrl_arm,ctrl_back),]

SEVERE.BACK <- DESeqDataSetFromMatrix(countData=severe.back, colData=severe.back.coldata, design = ~ site + condition)
SEVERE.BACK$condition <- relevel(SEVERE.BACK$condition, ref="ctrl")
SEVERE.BACK <- estimateSizeFactors(SEVERE.BACK)
SEVERE.BACK <- SEVERE.BACK[ rowSums(counts(SEVERE.BACK,normalized=TRUE) >= 5) >= 3, ]

## 7. scle_ctrl_arm vs all_ctrl
title.MILD.ARM <- 'MILD.ARM: Mild (Arm) VS Controls'
mild.arm <- the.data[, colnames(the.data) %in% c(scl.ctrl_arm,ctrl_arm,ctrl_back)]
mild.arm.coldata <- the.coldata[rownames(the.coldata) %in% c(scl.ctrl_arm,ctrl_arm,ctrl_back),]

MILD.ARM <- DESeqDataSetFromMatrix(countData=mild.arm, colData=mild.arm.coldata, design = ~ site + condition)
MILD.ARM$condition <- relevel(MILD.ARM$condition, ref="ctrl")
MILD.ARM <- estimateSizeFactors(MILD.ARM)
MILD.ARM <- MILD.ARM[ rowSums(counts(MILD.ARM,normalized=TRUE) >= 5) >= 3, ]

## 8. scle_ctrl_back vs all_ctrl
title.MILD.BACK <- 'MILD.BACK: Mild (Back) VS Controls'
mild.back <- the.data[, colnames(the.data) %in% c(scl.ctrl_back,ctrl_arm,ctrl_back)]
mild.back.coldata <- the.coldata[rownames(the.coldata) %in% c(scl.ctrl_back,ctrl_arm,ctrl_back),]

MILD.BACK <- DESeqDataSetFromMatrix(countData=mild.back, colData=mild.back.coldata, design = ~ site + condition)
MILD.BACK$condition <- relevel(MILD.BACK$condition, ref="ctrl")
MILD.BACK <- estimateSizeFactors(MILD.BACK)
MILD.BACK <- MILD.BACK[ rowSums(counts(MILD.BACK,normalized=TRUE) >= 5) >= 3, ]

## ## 9. scle_case vs scle_ctrl
title.SEVERE.MILD <- 'SEVERE.MILD: Severe VS Mild'
severe.mild <- the.data[, colnames(the.data) %in% c(scl.case_arm,scl.case_back,scl.ctrl_arm,scl.ctrl_back)]
severe.mild.coldata <- the.coldata[rownames(the.coldata) %in% c(scl.case_arm,scl.case_back,scl.ctrl_arm,scl.ctrl_back),]

SEVERE.MILD <- DESeqDataSetFromMatrix(countData=severe.mild, colData=severe.mild.coldata, design = ~ site + site:anti.scl.n + site:anti.scl)
SEVERE.MILD$anti.scl <- relevel(SEVERE.MILD$anti.scl, ref="0")
SEVERE.MILD <- estimateSizeFactors(SEVERE.MILD)
SEVERE.MILD <- SEVERE.MILD[ rowSums(counts(SEVERE.MILD,normalized=TRUE) >= 5) >= 3, ]

filtering.table$step.5 <- c(length(ALL),
                            length(CASE.ARM),
                            length(CASE.BACK),
                            length(CASE),
                            length(SEVERE.ARM),
                            length(SEVERE.BACK),
                            length(MILD.ARM),
                            length(MILD.BACK),
                            length(SEVERE.MILD))

## ------------------------------ RUN DIFFERENTIAL EXPRESSION

dds.ALL <- DESeq(ALL)
dds.CASE.ARM <- DESeq(CASE.ARM)
dds.CASE.BACK <- DESeq(CASE.BACK)
dds.CASE <- DESeq(CASE)
dds.SEVERE.ARM <- DESeq(SEVERE.ARM)
dds.SEVERE.BACK <- DESeq(SEVERE.BACK)
dds.MILD.ARM <- DESeq(MILD.ARM)
dds.MILD.BACK <- DESeq(MILD.BACK)
dds.SEVERE.MILD <- DESeq(SEVERE.MILD)

## EVALUATE THE BEST LFC THRESHOLD TO USE BY CREATING A TABLE FOR LFC AT 0.01 FDR FOR ALL COMPARISONS
check_lfc <- function(dds) {
    the_list <- vector()
    for (lfc in c(seq(0,2,0.1))) {
        res <- results(dds, alpha=0.01, lfcThreshold=lfc, altHypothesis="greaterAbs")
        ## summary(res)
        sig <- sum(res$padj < 0.01, na.rm=TRUE)
        the_list <- append(the_list, sig)
    }
    return(the_list)
}

lfc.table <- data.frame(row.names=(seq(0,2,0.1)))

lfc.table$ALL <- check_lfc(dds.ALL)
lfc.table$CASE.ARM <- check_lfc(dds.CASE.ARM)
lfc.table$CASE.BACK <- check_lfc(dds.CASE.BACK)
lfc.table$CASE <- check_lfc(dds.CASE)
lfc.table$SEVERE.ARM <- check_lfc(dds.SEVERE.ARM)
lfc.table$SEVERE.BACK <- check_lfc(dds.SEVERE.BACK)
lfc.table$MILD.ARM <- check_lfc(dds.MILD.ARM)
lfc.table$MILD.BACK <- check_lfc(dds.MILD.BACK)
lfc.table$SEVERE.MILD <- check_lfc(dds.SEVERE.MILD)

lfc.table$`FC Cutoff` <- format(round(as.numeric(row.names(lfc.table)), 2), nsmall = 2)
lfc.table <- lfc.table[,c(10,1:9)]

colnames(lfc.table) <- c("FC Cutoff","ALL", "CASE ARM", "CASE BACK", "CASE", "SEVERE ARM","SEVERE BACK", "MILD ARM", "MILD BACK", "SEVERE MILD")
names(lfc.table) <- paste0("\\textbf{",colnames(lfc.table),"}")

colnames(lfc.table)[1] <- paste0("\\multirow{-3}{*}{\\parbox{3.5em}{\\centering{",colnames(lfc.table[1]),"}}}")
colnames(lfc.table)[2] <- paste0("\\multirow{-2}{*}{\\parbox{3em}{\\centering{",colnames(lfc.table[2]),"}}}")
colnames(lfc.table)[c(3:5,8:9)] <- paste0("\\multirow{-2}{*}{\\parbox{3.5em}{\\centering{",colnames(lfc.table[c(3:5,8:9)]),"}}}")
colnames(lfc.table)[c(6:7,10)] <- paste0("\\multirow{-2}{*}{\\parbox{4.75em}{\\centering{",colnames(lfc.table[c(6:7,10)]),"}}}")

latex.columns <- "\\begin{tabular}{ L{3.5em} R{3em} R{3.5em} R{3.5em} R{3.5em} R{4.75em} R{4.75em} R{3.5em} R{3.5em} R{4.75em} }"

## \rowcolor{gray!25} & \multicolumn{9}{ c }{\textbf{# of singnificant genes in each comparison}} \\\hhline{>{\arrayrulecolor{gray!25}}->{\arrayrulecolor{black}}---------}
below.toprule <- c("\\rowcolor{gray!25} & \\multicolumn{9}{ c }{\\textbf{\\# of singnificant genes in each comparison}} \\\\\\hhline{>{\\arrayrulecolor{gray!25}}->{\\arrayrulecolor{black}}---------}","\\rowcolor{gray!25} &&&&&&&&&\\\\")
write.latex.table(lfc.table, latex.columns, below.toprule, paste0(tables.out.dir,"table.fc.cutoff",".tex"))

## GET THE DIFFERENTIALLY EXPRESSED GENES FROM THE DDS CREATED (USE LFC AND FDR DETERMINED ABOVE)
res.ALL <- results(dds.ALL, alpha=0.01, lfcThreshold=2, altHypothesis="greaterAbs")
summary(res.ALL)
sum(res.ALL$padj < 0.01, na.rm=TRUE)

res.CASE.ARM <- results(dds.CASE.ARM, alpha=0.01, lfcThreshold=2, altHypothesis="greaterAbs")
summary(res.CASE.ARM)
sum(res.CASE.ARM$padj < 0.01, na.rm=TRUE)

res.CASE.BACK <- results(dds.CASE.BACK, alpha=0.01, lfcThreshold=2, altHypothesis="greaterAbs")
summary(res.CASE.BACK)
sum(res.CASE.BACK$padj < 0.01, na.rm=TRUE)

res.CASE <- results(dds.CASE, alpha=0.01, lfcThreshold=1.1, altHypothesis="greaterAbs")
summary(res.CASE)
sum(res.CASE$padj < 0.01, na.rm=TRUE)

res.SEVERE.ARM <- results(dds.SEVERE.ARM, alpha=0.01, lfcThreshold=2, altHypothesis="greaterAbs")
summary(res.SEVERE.ARM)
sum(res.SEVERE.ARM$padj < 0.01, na.rm=TRUE)

res.SEVERE.BACK <- results(dds.SEVERE.BACK, alpha=0.01, lfcThreshold=2, altHypothesis="greaterAbs")
summary(res.SEVERE.BACK)
sum(res.SEVERE.BACK$padj < 0.01, na.rm=TRUE)

res.MILD.ARM <- results(dds.MILD.ARM, alpha=0.01, lfcThreshold=2, altHypothesis="greaterAbs")
summary(res.MILD.ARM)
sum(res.MILD.ARM$padj < 0.01, na.rm=TRUE)

res.MILD.BACK <- results(dds.MILD.BACK, alpha=0.01, lfcThreshold=2, altHypothesis="greaterAbs")
summary(res.MILD.BACK)
sum(res.MILD.BACK$padj < 0.01, na.rm=TRUE)

res.SEVERE.MILD <- results(dds.SEVERE.MILD, alpha=0.01, lfcThreshold=1.1, altHypothesis="greaterAbs")
summary(res.SEVERE.MILD)
sum(res.SEVERE.MILD$padj < 0.01, na.rm=TRUE)

filtering.table$step.6 <- c(sum(res.ALL$padj < 0.01, na.rm=TRUE),
                            sum(res.CASE.ARM$padj < 0.01, na.rm=TRUE),
                            sum(res.CASE.BACK$padj < 0.01, na.rm=TRUE),
                            sum(res.CASE$padj < 0.01, na.rm=TRUE),
                            sum(res.SEVERE.ARM$padj < 0.01, na.rm=TRUE),
                            sum(res.SEVERE.BACK$padj < 0.01, na.rm=TRUE),
                            sum(res.MILD.ARM$padj < 0.01, na.rm=TRUE),
                            sum(res.MILD.BACK$padj < 0.01, na.rm=TRUE),
                            sum(res.SEVERE.MILD$padj < 0.01, na.rm=TRUE))

filtering.table$Set <- row.names(filtering.table)
filtering.table <- filtering.table[,c(8,1:7)]

colnames(filtering.table) <- c("Set", "Genes", "Step 1\\tnote{*}", "Step 2\\tnote{$\\dagger$}", "Step 3\\tnote{$\\ddagger$}", "Step 4\\tnote{$\\mathsection$}", "Step 5\\tnote{$\\mathparagraph$}", "Step 6\\tnote{$\\ddagger$}")

names(filtering.table) <- paste0("\\textbf{",colnames(filtering.table),"}")

names(filtering.table)[1] <- paste0("\\multirow{-2}{*}{\\parbox{8.5em}{\\centering{",colnames(filtering.table[1]),"}}}") 

latex.columns <- "\\begin{tabular}{ L{8.5em} R{3.75em} R{3.75em} R{3.75em} R{3.75em} R{3.75em} R{3.75em} R{3.75em} }"


below.toprule <- "\\rowcolor{gray!25} & \\multicolumn{7}{ c }{\\textbf{\\# of genes at each filtering step}} \\\\\\hhline{>{\\arrayrulecolor{gray!25}}->{\\arrayrulecolor{black}}-------}"
write.latex.table(filtering.table, latex.columns, below.toprule, paste0(tables.out.dir,"table.filtering",".tex"))

## ------------------------------ CREATE LOG TRANSFORMATIONS AND VARIANCE STABILISING TRANSFORMATION

rld.ALL <- rlog(ALL, blind=FALSE)

rld.CASE.ARM <- rlog(CASE.ARM, blind=FALSE)
rld.CASE.BACK <- rlog(CASE.BACK, blind=FALSE)
rld.CASE <- rlog(CASE, blind=FALSE)
rld.SEVERE.ARM <- rlog(SEVERE.ARM, blind=FALSE)
rld.SEVERE.BACK <- rlog(SEVERE.BACK, blind=FALSE)
rld.MILD.ARM <- rlog(MILD.ARM, blind=FALSE)
rld.MILD.BACK <- rlog(MILD.BACK, blind=FALSE)
rld.SEVERE.MILD <- rlog(SEVERE.MILD, blind=FALSE)

## ------------------------------ GET DIFFERENTIALLY EXPRESSED RESULTS (FDR=0.05 or 5%, LOGFOLDCHANGE DEPENDS)

sig.ALL <- get_sig(res.ALL, 0.01, 2)
sig.CASE.ARM <- get_sig(res.CASE.ARM, 0.01, 2)
sig.CASE.BACK <- get_sig(res.CASE.BACK, 0.01, 2)
sig.CASE <- get_sig(res.CASE, 0.01, 1.1)
sig.SEVERE.ARM <- get_sig(res.SEVERE.ARM, 0.01, 2)
sig.SEVERE.BACK <- get_sig(res.SEVERE.BACK, 0.01, 2)
sig.MILD.ARM <- get_sig(res.MILD.ARM, 0.01, 2)
sig.MILD.BACK <- get_sig(res.MILD.BACK, 0.01, 2)
sig.SEVERE.MILD <- get_sig(res.SEVERE.MILD, 0.01, 1.1)

## CONFIRM - SUMMARISE FOR EACH COMPARISON
summary(sig.ALL)
summary(sig.CASE.ARM)
summary(sig.CASE.BACK)
summary(sig.CASE)
summary(sig.SEVERE.ARM)
summary(sig.SEVERE.BACK)
summary(sig.MILD.ARM)
summary(sig.MILD.BACK)
summary(sig.SEVERE.MILD)

## CONFIRM - GET NUMBER OF SIGNIFICANT GENES FOR EACH COMPARISON
sum(sig.ALL$padj < 0.01, na.rm=TRUE)
sum(sig.CASE.ARM$padj < 0.01, na.rm=TRUE)
sum(sig.CASE.BACK$padj < 0.01, na.rm=TRUE)
sum(sig.CASE$padj < 0.01, na.rm=TRUE)
sum(sig.SEVERE.ARM$padj < 0.01, na.rm=TRUE)
sum(sig.SEVERE.BACK$padj < 0.01, na.rm=TRUE)
sum(sig.MILD.ARM$padj < 0.01, na.rm=TRUE)
sum(sig.MILD.BACK$padj < 0.01, na.rm=TRUE)
sum(sig.SEVERE.MILD$padj < 0.01, na.rm=TRUE)

## CREATE A LIST OF SIGNIFICAN GENES (ENSEMBL IDS) FOR EACH COMPARISON
selected.ALL <- row.names(sig.ALL)
selected.CASE.ARM <- row.names(sig.CASE.ARM)
selected.CASE.BACK <- row.names(sig.CASE.BACK)
selected.CASE <- row.names(sig.CASE)
selected.SEVERE.ARM <- row.names(sig.SEVERE.ARM)
selected.SEVERE.BACK <- row.names(sig.SEVERE.BACK)
selected.MILD.ARM <- row.names(sig.MILD.ARM)
selected.MILD.BACK <- row.names(sig.MILD.BACK)
selected.SEVERE.MILD <- row.names(sig.SEVERE.MILD)

## LIST OF UPREGULATED GENES IN EACH GROUP
get.down <- function(sig,lfc=2) {
    down <- row.names(sig[which(sig$log2FoldChange < lfc),])
    down <- get_gene_names(down)
    return(down)
}

down.ALL <- get.down(sig.ALL)
down.CASE.ARM <- get.down(sig.CASE.ARM)
down.CASE.BACK <- get.down(sig.CASE.BACK)
down.CASE <- get.down(sig.CASE,1.1)
down.SEVERE.ARM <- get.down(sig.SEVERE.MILD)
down.SEVERE.BACK <- get.down(sig.SEVERE.BACK)
down.MILD.ARM <- get.down(sig.MILD.ARM)
down.MILD.BACK <- get.down(sig.MILD.BACK)
down.SEVERE.MILD <- get.down(sig.SEVERE.MILD,1.1)

## LIST OF DOWNREGULATED IN EACH GROUP
get.up <- function(sig,lfc=2) {
    up <- row.names(sig[which(sig$log2FoldChange > lfc),])
    up <- get_gene_names(up)
    return(up)
}

up.ALL <- get.up(sig.ALL)
up.CASE.ARM <- get.up(sig.CASE.ARM)
up.CASE.BACK <- get.up(sig.CASE.BACK)
up.CASE <- get.up(sig.CASE,1.1)
up.SEVERE.ARM <- get.up(sig.SEVERE.MILD)
up.SEVERE.BACK <- get.up(sig.SEVERE.BACK)
up.MILD.ARM <- get.up(sig.MILD.ARM)
up.MILD.BACK <- get.up(sig.MILD.BACK)
up.SEVERE.MILD <- get.up(sig.SEVERE.MILD,1.1)

## GET THE NAMES OF THE SIGNIFICANT GENES (ENSEMBL ID, DESCTIPTION, CHROMOSOME, BIOTYPE, HGNC SYMBOL)
genes.ALL <- get_gene_names(selected.ALL)
genes.CASE.ARM <- get_gene_names(selected.CASE.ARM)
genes.CASE.BACK <- get_gene_names(selected.CASE.BACK)
genes.CASE <- get_gene_names(selected.CASE)
genes.SEVERE.ARM <- get_gene_names(selected.SEVERE.ARM)
genes.SEVERE.BACK <- get_gene_names(selected.SEVERE.BACK)
genes.MILD.ARM <- get_gene_names(selected.MILD.ARM)
genes.MILD.BACK <- get_gene_names(selected.MILD.BACK)
genes.SEVERE.MILD <- get_gene_names(selected.SEVERE.MILD)

diff_expr_results <- function(sig.res,genes,dds){
    ## GET DATA FROM SIG
    sig <- as.data.frame(sig.res)[order(sig.res$padj),]
    ## GET DATA FROM DDS (NORMALIZED COUNTS)
    dds <- as.data.frame(counts(dds, normalized=TRUE))
    ## MERGE THE TWO
    expr_table <- merge(sig,dds,by="row.names",sort=FALSE)
    colnames(expr_table)[1] <- "GeneID"
    colnames(expr_table)[3] <- "log2FC"
    ## GET THE DESCRIPTION, CROMOSOME AND HGNC ID FROM THE ANNOTATION TABLE OF THE DATASET
    expr_table$Description <- genes$description[match(expr_table$GeneID, genes$ensembl_id)]
    expr_table$Chr <- genes$chromosome[match(expr_table$GeneID, genes$ensembl_id)]
    expr_table$GeneID <- ifelse(genes$HGNC[match(expr_table$GeneID, genes$ensembl_id)] == "",
                                expr_table$GeneID,
                                genes$HGNC[match(expr_table$GeneID, genes$ensembl_id)])
    ## ORDER THE TABLE
    len <- length(colnames(expr_table))
    expr_table <- expr_table[c(1, (len-1):len, 2:(len-2))]
    ## ## FORMAT THE NUMBERS
    expr_table$pvalue <- formatC(expr_table$pvalue, format = "e", digits = 2)
    expr_table$padj <- formatC(expr_table$padj, format = "e", digits = 2)
    expr_table[,c(4:7,10:len)] <- sapply(expr_table[,c(4:7,10:len)], function(number) format(round(number, 2), nsmall = 2))
    ## STATS ONLY
    stats_table <- expr_table[ ,c(1:9)]
    
    ## WRITE OUT THE TABLE
    name <- deparse(substitute(sig.res))
    out.table <- paste0(tables.out.dir,"stats_table", substr(name,4,nchar(name)),".tex")
    names(stats_table) <- paste0('\\textbf{',colnames(stats_table),"}")

    lft <- ifelse(substr(name,5,nchar(name)) %in% c("CASE","SEVERE.MILD"),1.1,2)
    
    ## x <- "\\begin{tabular}{ L{8.5em} L{28em} R{2em} R{5.1em} R{3.75em} R{3em} R{2.5em} R{3.75em} R{3.75em} }"
    table.head <- "\\begin{longtable}{ L{8.5em} L{28em} R{2em} R{5.1em} R{3.75em} R{3em} R{2.5em} R{3.75em} R{3.75em} }"
    table.caption <- paste0("\\caption[Differentially expressed genes in the ",substr(name,5,nchar(name))," comparison set.]{\\textbf{Significantly differentially expressed genes in the ",substr(name,5,nchar(name))," comparison set, where FDR$ \\leq$ 0.01 (1\\%) \\& {\\textbar}Log$_2$FC{\\textbar} $\\geq$",lft,".}}")
    table.label <- paste0("\\label{tab:sig.",tolower(substr(name,5,nchar(name))),"}\\\\")
    ##
    write.table(table.head, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" & ", eol=" \n")
    write.table(table.caption, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
    write.table(table.label, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
    ##
    write.table("\\toprule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    write.table(paste0("\\rowcolor{gray!25} ", paste(names(stats_table),collapse=" & ")),
                file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
    write.table("\\midrule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    write.table(stats_table, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
    write.table("\\bottomrule\n\\end{longtable}", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    return(expr_table)
}

table.ALL <- (diff_expr_results(sig.ALL,genes.ALL,dds.ALL))
table.CASE.ARM <- (diff_expr_results(sig.CASE.ARM,genes.CASE.ARM,dds.CASE.ARM))
table.CASE.BACK <- (diff_expr_results(sig.CASE.BACK,genes.CASE.BACK,dds.CASE.BACK))
table.CASE <- (diff_expr_results(sig.CASE,genes.CASE,dds.CASE))
table.SEVERE.ARM <- (diff_expr_results(sig.SEVERE.ARM,genes.SEVERE.ARM,dds.SEVERE.ARM))
table.SEVERE.BACK <- (diff_expr_results(sig.SEVERE.BACK,genes.SEVERE.BACK,dds.SEVERE.BACK))
table.MILD.ARM <- (diff_expr_results(sig.MILD.ARM,genes.MILD.ARM,dds.MILD.ARM))
table.MILD.BACK <- (diff_expr_results(sig.MILD.BACK,genes.MILD.BACK,dds.MILD.BACK))
table.SEVERE.MILD <- (diff_expr_results(sig.SEVERE.MILD,genes.SEVERE.MILD,dds.SEVERE.MILD))


## ------------------------------ DESCRIPTIVE VISUALIZATIONS

## ===== PLOT PCA FOR EACH SAMPLE AND PUT IN ONE FILE
pca.ALL <- plot_pca(rld.ALL)
pca.CASE.ARM <- plot_pca(rld.CASE.ARM)
pca.CASE.BACK <- plot_pca(rld.CASE.BACK)
pca.CASE <- plot_pca(rld.CASE)
pca.SEVERE.ARM <- plot_pca(rld.SEVERE.ARM)
pca.SEVERE.BACK <- plot_pca(rld.SEVERE.BACK)
pca.MILD.ARM <- plot_pca(rld.MILD.ARM)
pca.MILD.BACK <- plot_pca(rld.MILD.BACK)
pca.SEVERE.MILD <- plot_pca(rld.SEVERE.MILD)

## CREATE LAYOUT 
lay <- rbind(c(1,2,3), c(4,5,6), c(7,8,9),10)

mylegend <- get_legend(pca.CASE.ARM)

## PLOT ALL IN ONE PAGE
pdf(paste0(figures.out.dir,'plot.pca.pdf'), paper='a4r',width = 0, height = 0, onefile=FALSE)
grid.arrange(pca.ALL + ggtitle(title.ALL) + theme(legend.position="none"),
             pca.CASE.ARM + ggtitle(title.CASE.ARM) + theme(legend.position="none"),
             pca.CASE.BACK + ggtitle(title.CASE.BACK) + theme(legend.position="none"),
             pca.CASE + ggtitle(title.CASE) + theme(legend.position="none"),
             pca.SEVERE.ARM + ggtitle(title.SEVERE.ARM) + theme(legend.position="none"),
             pca.SEVERE.BACK + ggtitle(title.SEVERE.BACK) + theme(legend.position="none"),
             pca.MILD.ARM + ggtitle(title.MILD.ARM) + theme(legend.position="none"),
             pca.MILD.BACK + ggtitle(title.MILD.BACK) + theme(legend.position="none"),
             pca.SEVERE.MILD + ggtitle(title.SEVERE.MILD) + theme(legend.position="none"),
             layout_matrix = lay,
             widths = c(1,1,1),
             heights = c(1,1,1,0.05),
             mylegend
             ## top=textGrob("PCA Plots for the different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()


## PLOT ALL IN ONE PAGE
pdf(paste0(figures.out.dir,'plot.pca.arms.pdf'), paper='a4',width = 0, height = 0, onefile=FALSE)
grid.arrange(pca.ALL + ggtitle(title.ALL) + theme(legend.position="none"),
             pca.CASE.ARM + ggtitle(title.CASE.ARM) + theme(legend.position="none"),
             pca.SEVERE.ARM + ggtitle(title.SEVERE.ARM) + theme(legend.position="none"),
             pca.MILD.ARM + ggtitle(title.MILD.ARM) + theme(legend.position="none"),
             layout_matrix = rbind(c(1,2),c(3,4),5),
             widths = c(1,1),
             heights = c(1,1,0.05,1,1),
             mylegend
             ## top=textGrob("PCA Plots for the different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()


## PLOT ALL IN ONE PAGE
pdf(paste0(figures.out.dir,'plot.pca.backs.pdf'), paper='a4',width = 0, height = 0, onefile=FALSE)
grid.arrange(pca.CASE.BACK + ggtitle(title.CASE.BACK) + theme(legend.position="none"),
             pca.SEVERE.BACK + ggtitle(title.SEVERE.BACK) + theme(legend.position="none"),
             pca.MILD.BACK + ggtitle(title.MILD.BACK) + theme(legend.position="none"),
             pca.CASE + ggtitle(title.CASE) + theme(legend.position="none"),
             pca.SEVERE.MILD + ggtitle(title.SEVERE.MILD) + theme(legend.position="none"),
             layout_matrix = rbind(c(1,2),c(3,4),c(5,5),c(6,6)),
             widths = c(1,1),
             heights = c(1,1,1,0.05,1),
             mylegend
             ## top=textGrob("PCA Plots for the different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()


## ===== CREATE VOLCANO PLOT FROM EACH COMPARISON DATASET AND PUT THEM IN ONE FILE (EASY VISUALISATION)
volcano.ALL <- my_volcano(res.ALL)
volcano.CASE.ARM <- my_volcano(res.CASE.ARM)
volcano.CASE.BACK <- my_volcano(res.CASE.BACK)
volcano.CASE <- my_volcano(res.CASE, 1.1, 0.01)
volcano.SEVERE.ARM <- my_volcano(res.SEVERE.ARM)
volcano.SEVERE.BACK <- my_volcano(res.SEVERE.BACK)
volcano.MILD.ARM <- my_volcano(res.MILD.ARM)
volcano.MILD.BACK <- my_volcano(res.MILD.BACK)
volcano.SEVERE.MILD <- my_volcano(res.SEVERE.MILD, 1.1, 0.01)

## CREATE LAYOUT
lay <- rbind(c(1,2,3), c(4,5,6), c(7,8,9))

## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'plot.volcano.pdf'), paper='a4r', width = 0, height = 0, onefile=FALSE)
grid.arrange(volcano.ALL + ggtitle(title.ALL),
             volcano.CASE.ARM + ggtitle(title.CASE.ARM),
             volcano.CASE.BACK + ggtitle(title.CASE.BACK),
             volcano.CASE + ggtitle(title.CASE),
             volcano.SEVERE.ARM + ggtitle(title.SEVERE.ARM),
             volcano.SEVERE.BACK + ggtitle(title.SEVERE.BACK),
             volcano.MILD.ARM + ggtitle(title.MILD.ARM),
             volcano.MILD.BACK + ggtitle(title.MILD.BACK),
             volcano.SEVERE.MILD + ggtitle(title.SEVERE.MILD),
             layout_matrix = lay,
             widths = c(1,1,1),
             heights = c(1,1,1)
             ## top=textGrob("Volcano Plots for the different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()


## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'plot.volcano.arms.pdf'), paper='a4', width = 0, height = 0, onefile=FALSE)
grid.arrange(volcano.ALL + ggtitle(title.ALL),
             volcano.CASE.ARM + ggtitle(title.CASE.ARM),
             volcano.SEVERE.ARM + ggtitle(title.SEVERE.ARM),
             volcano.MILD.ARM + ggtitle(title.MILD.ARM),
             layout_matrix = rbind(c(1,2),c(3,4),5),
             widths = c(1,1),
             heights = c(1,1,1,1)
             ## top=textGrob("Volcano Plots for the different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()

## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'plot.volcano.backs.pdf'), paper='a4', width = 0, height = 0, onefile=FALSE)
grid.arrange(volcano.CASE.BACK + ggtitle(title.CASE.BACK),
             volcano.SEVERE.BACK + ggtitle(title.SEVERE.BACK),
             volcano.MILD.BACK + ggtitle(title.MILD.BACK),
             volcano.CASE + ggtitle(title.CASE),
             volcano.SEVERE.MILD + ggtitle(title.SEVERE.MILD),
             layout_matrix = rbind(c(1,2),c(3,4),5,6),
             widths = c(1,1),
             heights = c(1,1,1,1)
             ## top=textGrob("Volcano Plots for the different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()

## ===== PLOT HEATMAPS - ADD KNOWN TARGETS AND OTHER DESCRIPTIVE DATA TO THE HEATMAP
## Get the known targets and check how many of my data are targets
## https://www.targetvalidation.org/disease/EFO_0000717/associations?fcts=datatype:genetic_association;known_drug;literature
## 31 October 2018

## GET KNOWN TARGETS
known_targets <- read.csv(file = 'targets_associated_with_systemic_scleroderma.csv', row.names = 1, header = TRUE)

known_in_my_data.ALL <- subset(genes.ALL, HGNC %in% row.names(known_targets))
known_in_my_data.CASE.ARM <- subset(genes.CASE.ARM, HGNC %in% row.names(known_targets))
known_in_my_data.CASE.BACK <- subset(genes.CASE.BACK, HGNC %in% row.names(known_targets))
known_in_my_data.CASE <- subset(genes.CASE, HGNC %in% row.names(known_targets))
known_in_my_data.SEVERE.ARM <- subset(genes.SEVERE.ARM, HGNC %in% row.names(known_targets))
known_in_my_data.SEVERE.BACK <- subset(genes.SEVERE.BACK, HGNC %in% row.names(known_targets))
known_in_my_data.MILD.ARM <- subset(genes.MILD.ARM, HGNC %in% row.names(known_targets))
known_in_my_data.MILD.BACK <- subset(genes.MILD.BACK, HGNC %in% row.names(known_targets))
known_in_my_data.SEVERE.MILD <- subset(genes.SEVERE.MILD, HGNC %in% row.names(known_targets))

## PLOT THE HEATMAPS
pdf(paste0(figures.out.dir,'plot.heatmap.ALL.pdf'), paper='a4',onefile=FALSE)
heat.ALL <- plot_heat(rld.ALL,selected.ALL,known_in_my_data.ALL,genes.ALL,title.ALL)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.CASE.ARM.pdf'), paper='a4', onefile=FALSE)
heat.CASE.ARM <- plot_heat(rld.CASE.ARM,selected.CASE.ARM,known_in_my_data.CASE.ARM,genes.CASE.ARM,title.CASE.ARM)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.CASE.BACK.pdf'), paper='a4', onefile=FALSE)
heat.CASE.BACK <- plot_heat(rld.CASE.BACK,selected.CASE.BACK,known_in_my_data.CASE.BACK,genes.CASE.BACK,title.CASE.BACK)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.CASE.pdf'), paper='a4', onefile=FALSE)
heat.CASE <- plot_heat(rld.CASE,selected.CASE,known_in_my_data.CASE,genes.CASE,title.CASE)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.SEVERE.ARM.pdf'), paper='a4', onefile=FALSE)
heat.SEVERE.ARM <- plot_heat(rld.SEVERE.ARM,selected.SEVERE.ARM,known_in_my_data.SEVERE.ARM,genes.SEVERE.ARM,title.SEVERE.ARM)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.SEVERE.BACK.pdf'), paper='a4', onefile=FALSE)
heat.SEVERE.BACK <- plot_heat(rld.SEVERE.BACK,selected.SEVERE.BACK,known_in_my_data.SEVERE.BACK,genes.SEVERE.BACK,title.SEVERE.BACK)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.MILD.ARM.pdf'), paper='a4', onefile=FALSE)
heat.MILD.ARM <- plot_heat(rld.MILD.ARM,selected.MILD.ARM,known_in_my_data.MILD.ARM,genes.MILD.ARM,title.MILD.ARM)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.MILD.BACK.pdf'), paper='a4', onefile=FALSE)
heat.MILD.BACK <- plot_heat(rld.MILD.BACK,selected.MILD.BACK,known_in_my_data.MILD.BACK,genes.MILD.BACK,title.MILD.BACK)
dev.off()

pdf(paste0(figures.out.dir,'plot.heatmap.SEVERE.MILD.pdf'), paper='a4', onefile=FALSE)
heat.SEVERE.MILD <- plot_heat(rld.SEVERE.MILD,selected.SEVERE.MILD,known_in_my_data.SEVERE.MILD,genes.SEVERE.MILD,title.SEVERE.MILD)
dev.off()

## #####------------------------------  Visualizations of the data for QC ------------------------------#####
## # Number of rejections
pdf(paste0(figures.out.dir,'plot.rejections.pdf'), paper='a4r',width = 0, height = 0, onefile=FALSE)
par(mfrow = c(3, 3))
rej.ALL <- plot_rejections_plot(res.ALL) #+ title(title.1)
rej.CASE.ARM <- plot_rejections_plot(res.CASE.ARM) #+ title(title.2)
rej.CASE.BACK <- plot_rejections_plot(res.CASE.BACK) #+ title(title.3)
rej.SEVERE.ARM <- plot_rejections_plot(res.SEVERE.ARM) #+ title(title.4)
rej.SEVERE.BACK <- plot_rejections_plot(res.SEVERE.BACK) #+ title(title.5)
rej.MILD.ARM <- plot_rejections_plot(res.MILD.ARM) #+ title(title.6)
rej.MILD.BACK <- plot_rejections_plot(res.MILD.BACK) #+ title(title.7)
rej.SEVERE.MILD <- plot_rejections_plot(res.SEVERE.MILD) #+ title(title.8)
rej.CASE <- plot_rejections_plot(res.CASE) #+ title(title.9)
dev.off()
par(mfrow = c(1, 1))

## # Boxplot.
pdf(paste0(figures.out.dir,'plot_boxplot.pdf'), width = 12, height = 7, onefile=FALSE)
par(mfrow = c(3, 3))
b_p.1 <- plot_boxplots(dds.ALL, title.ALL)
b_p.2 <- plot_boxplots(dds.CASE.ARM, title.CASE.ARM)
b_p.3 <- plot_boxplots(dds.CASE.BACK, title.CASE.BACK)
b_p.4 <- plot_boxplots(dds.CASE, title.CASE)
b_p.5 <- plot_boxplots(dds.SEVERE.ARM, title.SEVERE.ARM)
b_p.6 <- plot_boxplots(dds.SEVERE.BACK, title.SEVERE.BACK)
b_p.7 <- plot_boxplots(dds.MILD.ARM, title.MILD.ARM)
b_p.8 <- plot_boxplots(dds.MILD.BACK, title.MILD.BACK)
b_p.9 <- plot_boxplots(dds.SEVERE.MILD, title.SEVERE.MILD)
dev.off()
par(mfrow = c(1, 1))

## ------------------------------ GET ENSEMBL BIOTYPES OF SIGNIFICANT GENES AFTER ANALYSES 

sig.biotypes <- biotypes[which(biotypes$ensembl_id %in% c(selected.ALL, selected.CASE.ARM, selected.CASE.BACK,
                                                          selected.CASE, selected.SEVERE.ARM, selected.SEVERE.BACK,
                                                          selected.MILD.ARM, selected.MILD.BACK, selected.SEVERE.MILD)),]

## GET THE FREQUENCY TABLE
freq.sig <- as.data.frame(table(sig.biotypes$biotype))
freq.sig <- freq.sig[order(freq.sig$Freq, decreasing=TRUE),]
colnames(freq.sig) <- c("Gene Biotype","Number of Genes")

## PLOT
pdf(paste0(figures.out.dir,'plot.biotypes.after.pdf'), paper='a4', onefile=FALSE)
ggplot(data=freq.sig, aes(x=reorder(`Gene Biotype`, -`Number of Genes`), y=`Number of Genes`)) +
    geom_bar(stat="identity", width=0.7) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5)) +
    xlab("Gene Biotypes") +
    ylab("Counts") +
    ggtitle("Categories of Differentially Expressed Genes")
dev.off()


## ------------------------------ UPSET ANALYSIS TO FIND GENES SIGNIFICANT
## CREATE A MATRIX FOR THE SIGNIFICANT 
analysis.genes <- unique(c(selected.ALL,selected.CASE.ARM,selected.CASE.BACK,
                           selected.CASE,selected.SEVERE.ARM,selected.SEVERE.BACK,
                           selected.MILD.ARM,selected.MILD.BACK,selected.SEVERE.MILD))

the_table <- as.data.frame(analysis.genes)

the_table$ALL <- as.numeric(the_table$analysis.genes %in% selected.ALL)
the_table$CASE.ARM <- as.numeric(the_table$analysis.genes %in% selected.CASE.ARM)
the_table$CASE.BACK <- as.numeric(the_table$analysis.genes %in% selected.CASE.BACK)
the_table$CASE <- as.numeric(the_table$analysis.genes %in% selected.CASE)
the_table$SEVERE.ARM <- as.numeric(the_table$analysis.genes %in% selected.SEVERE.ARM)
the_table$SEVERE.BACK <- as.numeric(the_table$analysis.genes %in% selected.SEVERE.BACK)
the_table$MILD.ARM <- as.numeric(the_table$analysis.genes %in% selected.MILD.ARM)
the_table$MILD.BACK <- as.numeric(the_table$analysis.genes %in% selected.MILD.BACK)
the_table$SEVERE.MILD <- as.numeric(the_table$analysis.genes %in% selected.SEVERE.MILD)

rownames(the_table) <- the_table$analysis.genes


the_table <- the_table[,-1]

## ===== WRITE OUT THE UPSET CSV TABLE FOR USE WITH THE INTERACTIVE UPSET
the_names <- get_gene_names(rownames(the_table))
the_table$HGNC <- the_names$HGNC[match(row.names(the_table), the_names$ensembl_id)]
the_table$NAME <- the_names$description[match(row.names(the_table), the_names$ensembl_id)]
write.csv(the_table, file=paste0(tables.out.dir,"scleroderma.csv"))

## ===== SIG GENES
final_genes <- the_table[!(the_table$CASE == 1 | the_table$SEVERE.MILD == 1 |
                           (the_table$CASE == 1 & the_table$SEVERE.MILD == 1 )),]

final_genes.ensembl <- unique(row.names(final_genes))

arm_genes.hgnc <- filter(the_table,(ALL == 1 |
                                    CASE.ARM == 1 |
                                    SEVERE.ARM == 1 |
                                    MILD.ARM == 1 ) &
                                   (CASE != 1 &
                                    CASE.BACK != 1 &
                                    CASE != 1 &
                                    SEVERE.BACK != 1 &
                                    MILD.BACK != 1 &
                                    SEVERE.MILD != 1))$HGNC

most.hgnc <- filter(the_table,(ALL == 1 &
                               CASE.ARM == 1 &
                               CASE.BACK == 1 &
                               SEVERE.ARM == 1 &
                               SEVERE.BACK == 1 &
                               MILD.ARM == 1 &
                               MILD.BACK == 1 ))$HGNC

arm_genes.hgnc <- unique(c(arm_genes.hgnc,most.hgnc))

arm_genes.hgnc <- arm_genes.hgnc[arm_genes.hgnc != ""]

final_genes.hgnc <- unique(the_names[the_names$ensembl_id %in% final_genes.ensembl,]$HGNC)
final_genes.hgnc <- final_genes.hgnc[final_genes.hgnc != ""]
write.csv(final_genes[,colnames(final_genes) %in% c("HGNC","NAME")], file=paste0(tables.out.dir,"final.gene.list.csv"), row.names=TRUE,col.names=TRUE)

the_table$ensembl <- row.names(the_table)

length(row.names(the_table))

## REMOVE NA'S AND ADD 0 INSTEAD
the_table[is.na(the_table)] <- 0

## ------------------------------ THIS SECTION IS FOR UPSETR ANALYSIS
full.set  <- c("ALL","CASE.ARM","CASE.BACK",
               "CASE","SEVERE.ARM","SEVERE.BACK",
               "MILD.ARM","MILD.BACK","SEVERE.MILD")

sets <- names(the_table[1:9])
compare <- c(title.ALL, title.CASE.ARM, title.CASE.BACK, title.CASE, title.SEVERE.ARM, title.SEVERE.BACK, title.MILD.ARM, title.MILD.BACK, title.SEVERE.MILD)
             
metadata <- as.data.frame(cbind(sets, compare))
names(metadata) <- c("sets", "Comparisons")

## FUNCTION FOR PLOTTING VOLCANO PLOTS ALONGSIDE UPSET PLOT
plot_this <- function(upset.data,compare) {
    name <- compare
    compare_name <- paste0("compare.",compare)
    title <- metadata$Comparisons[match(subset(metadata$sets, metadata$sets == compare_name), colnames(the_table))]
    print(title)
    lft <- ifelse(name %in% c("CASE","SEVERE.MILD"),1.1,2)
    print(lft)
    
    dds.set <- eval(as.name(paste0('dds.',name)))

    ## RAW SET
    raw.set <- as.data.frame(results(dds.set))
    raw.set <- raw.set[,colnames(raw.set) %in% c("log2FoldChange","pvalue","padj")]
    raw.set$log.pvalue <- -log10(raw.set[["pvalue"]])
    
    ## PVALUE CUTOFF
    pval.cut <- max(na.omit(raw.set[which(raw.set$padj < 0.01 & abs(raw.set$log2FoldChange) > lft),])$pvalue)
    
    colored <- upset.data[which(upset.data$color != "gray65" & upset.data$ensembl %in% row.names(raw.set)),]
    #print(colored)

    raw.set$color <- colored$color[match(row.names(raw.set), colored$ensembl)]
    #print(raw.set[which(raw.set$color != "NA"),])
    
    color.set <- raw.set[which(raw.set$color != "NA"),]
    
    ## ## THE PLOT ITSELF
    ## plot <- (ggplot(data=raw.set[which(!row.names(raw.set) %in% color.set$ensembl),], aes(x=log2FoldChange, y=log.pvalue)) +
    ## geom_point(size=0.1, alpha=1 , color = "#E0E0E0") +
    plot <- (ggplot(data=color.set, aes(x=log2FoldChange, y=log.pvalue, color=color)) +
             geom_point(size=0.1, alpha=1) +
             scale_color_identity() +
             theme_bw() +
             annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste(unlist(eval(as.name(compare_name))[[1]]), collapse = " "), size=1.6) +
             annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste(unlist(eval(as.name(compare_name))[[2]]), collapse = " "), size=1.6) +
             annotate("rect", xmin = -lft, xmax = lft, ymin = -log10(pval.cut), ymax = Inf, alpha = 0.1) +
             annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(pval.cut), alpha = 0.1) +
             geom_vline(aes(xintercept=-lft), size=0.2, alpha=0.5, linetype="dashed") +
             geom_vline(aes(xintercept=lft), size=0.2, alpha=0.5, linetype="dashed") +
             geom_hline(aes(yintercept=-log10(pval.cut)), size=0.2, alpha=0.5, linetype="dashed") +
             ggtitle(title) +
             xlab("log2(Fold-Change)") +
             ylab("-log10(P-Value)") +
             scale_x_continuous(breaks = seq(-50,50,by=5), limits=c(-15,15)) +
             scale_y_continuous(breaks = seq(0,50,5),limits=c(0,35)) +
             theme(panel.grid.major = element_blank(),
                   legend.position="none",
                   panel.grid.minor = element_blank(),
                   axis.title = element_text(size=5),
                   axis.text = element_text(size=4),
                   plot.title = element_text(size=7, hjust = 0.5),
                   plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm")
                   )
    )
}

## ARMS QUERY

## PLOT THE UPSETR PLOT, VOLCANO PLOT, QUERIES AND INTERSECTIONS
pdf(paste0(figures.out.dir,'plot.upset.intersection.pdf'), paper='a4', width = 0, height = 0, onefile=FALSE)
upset(the_table[,c(1:9,12)],#which(!colnames(the_table) %in% c("HGNC","NAME"))],
      sets = full.set, nintersects = NA,
      mainbar.y.label = "Differentially Expressed Genes",
      sets.x.label = "Number of Differentially Expressed Genes", nsets = 9,
      keep.order = TRUE, query.legend = "none", sets.bar.color = "#56B4E9",
      order.by = 'degree',# decreasing = FALSE,
      matrix.dot.alpha = 0.2, mb.ratio = c(0.65, 0.35),
      text.scale = 0.75,
      queries = list(
          list(query = intersects, params = "ALL",color = "#4363d8", active = T),
          list(query = intersects, params = "CASE.ARM",color = "#e6194B", active = T),
          list(query = intersects, params = "SEVERE.ARM",color = "#f58231", active = T),
          list(query = intersects, params = "MILD.ARM",color = "#911eb4", active = T),
          ## ## ARMS INTERSECTS
          list(query = intersects, params = list("CASE.ARM", "SEVERE.ARM","MILD.ARM"),
               color = "#3cb44b", active = T),
          list(query = intersects, params = list("CASE.ARM", "SEVERE.ARM"),
               color = "#3cb44b", active = T),
          list(query = intersects, params = list("CASE.ARM","MILD.ARM"),
               color = "#3cb44b", active = T),
          list(query = intersects, params = list("CASE.ARM","ALL"),
               color = "#3cb44b", active = T),
          list(query = intersects, params = list("ALL", "CASE.ARM","CASE.BACK", "SEVERE.ARM","SEVERE.BACK","MILD.ARM","MILD.BACK"),
               color = "#800000", active = T)
      ),
      attribute.plots = list(
          gridrows = 200,
          plots = list(
              list(plot = plot_this, x = "ALL", queries=T),
              list(plot = plot_this, x = "CASE.ARM", queries=T),
              list(plot = plot_this, x = "CASE.BACK", queries=T),
              list(plot = plot_this, x = "CASE", queries=T),
              list(plot = plot_this, x = "SEVERE.ARM", queries=T),
              list(plot = plot_this, x = "SEVERE.BACK", queries=T),
              list(plot = plot_this, x = "MILD.ARM", queries=T),
              list(plot = plot_this, x = "MILD.BACK", queries=T),              
              list(plot = plot_this, x = "SEVERE.MILD", queries=T)),
          ncols = 3)
      )
dev.off()


## ENRICHR
dbs <- listEnrichrDbs()$libraryName
dbs <- dbs[dbs != "NCI-Nature_2015"]
my.db.list=c("GO_Molecular_Function_2013",
             "GO_Biological_Process_2013",
             "GO_Cellular_Component_2013",
             "GO_Molecular_Function_2015",
             "GO_Biological_Process_2015",
             "GO_Cellular_Component_2015",
             "GO_Molecular_Function_2017",
             "GO_Biological_Process_2017",
             "GO_Cellular_Component_2017",
             "GO_Molecular_Function_2017b",
             "GO_Biological_Process_2017b",
             "GO_Cellular_Component_2017b",
             "GO_Molecular_Function_2018",
             "GO_Biological_Process_2018",
             "GO_Cellular_Component_2018",
             "WikiPathways_2013",
             "WikiPathways_2015",
             "WikiPathways_2016",
             "KEGG_2013",
             "KEGG_2015",
             "KEGG_2016",
             "BioCarta_2013",
             "BioCarta_2015",
             "BioCarta_2016",
             "Panther_2015",
             "Panther_2016",
             "Reactome_2013",
             "Reactome_2015",
             "Reactome_2016",
             "OMIM_Disease",
             "OMIM_Expanded",
             "Jensen_DISEASES",
             "Jensen_COMPARTMENTS",
             "Jensen_TISSUES")

dbs <- dbs[which(dbs %in% my.db.list)]

## ALL SIGNIFICANT GENES COMBINED!
enriched.ALL.COMBINED <- enrichr(final_genes.hgnc, dbs)

## ARMS ENRICHMENT
enriched.ARMS <- enrichr(arm_genes.hgnc, dbs)

## SIGNIFICANT GENES BY SUBSET
enriched.ALL <- enrichr(unique(genes.ALL$HGNC[genes.ALL$HGNC != ""]), dbs)
enriched.CASE.ARM <- enrichr(unique(genes.CASE.ARM$HGNC[genes.CASE.ARM$HGNC != ""]), dbs)
enriched.CASE.BACK <- enrichr(unique(genes.CASE.BACK$HGNC[genes.CASE.BACK$HGNC != ""]), dbs)
enriched.CASE <- enrichr(unique(genes.CASE$HGNC[genes.CASE$HGNC != ""]), dbs)
enriched.SEVERE.ARM <- enrichr(unique(genes.SEVERE.ARM$HGNC[genes.SEVERE.ARM$HGNC != ""]), dbs)
enriched.SEVERE.BACK <- enrichr(unique(genes.SEVERE.BACK$HGNC[genes.SEVERE.BACK$HGNC != ""]), dbs)
enriched.MILD.ARM <- enrichr(unique(genes.MILD.ARM$HGNC[genes.MILD.ARM$HGNC != ""]), dbs)
enriched.MILD.BACK <- enrichr(unique(genes.MILD.BACK$HGNC[genes.MILD.BACK$HGNC != ""]), dbs)
enriched.SEVERE.MILD <- enrichr(unique(genes.SEVERE.MILD$HGNC[genes.SEVERE.MILD$HGNC != ""]), dbs)

put_enrichment_into_files <- function(set) {
    set.name <- deparse(substitute(set))
    set.name <- substr(set.name,10,nchar(set.name))
    lapply(names(set), function(db) 
        if (nrow(set[[db]]) >= 1 & nrow(set[[db]][which(set[[db]]$Adjusted.P.value <= 0.05),]) >= 1) {
            print(c("There are ",nrow(set[[db]])," in this set"))
            write.csv(as.data.frame(set[[db]][,c("Term","Overlap","P.value","Adjusted.P.value","Combined.Score","Genes")][which(set[[db]]$Adjusted.P.value <= 0.05),]),
                      file=paste0(enrichr.out.dir,"enrichr_annotation_",db,"_",set.name,".csv"),
                      row.names=FALSE,col.names=TRUE)
        } else {}
        )
}


put_enrichment_into_files(enriched.ALL.COMBINED)
put_enrichment_into_files(enriched.ARMS)

put_enrichment_into_files(enriched.ALL)
put_enrichment_into_files(enriched.CASE.ARM)
put_enrichment_into_files(enriched.CASE.BACK)
put_enrichment_into_files(enriched.CASE)
put_enrichment_into_files(enriched.SEVERE.ARM)
put_enrichment_into_files(enriched.SEVERE.BACK)
put_enrichment_into_files(enriched.MILD.ARM)
put_enrichment_into_files(enriched.MILD.BACK)
put_enrichment_into_files(enriched.SEVERE.MILD)

## ATTEMPT AT VISUALIZING GO TERMS

plot.enrichr <- function(set,db.in){
    if (nrow(set[[db.in]]) >= 1 & nrow(set[[db.in]][which(set[[db.in]]$Adjusted.P.value <= 0.05),]) >= 1) {
        the.set <- set[[db.in]][order(set[[db.in]]$Adjusted.P.value, decreasing=FALSE),]
        if (nrow(the.set[which(the.set$Adjusted.P.value < 0.05 & str_detect(the.set$Term,"_Mus musculus") == FALSE),]) >= 1) {
            the.set <- head(the.set[which(the.set$Adjusted.P.value < 0.05 & str_detect(the.set$Term,"_Mus musculus") == FALSE),],20)
            the.set$gene.counts <- count.fields(textConnection(the.set$Genes),sep=";")
            the.set <- the.set[,which(colnames(the.set) %in% c("Term","P.value","Combined.Score","gene.counts"))]
            if (str_detect(db.in,"GO_") == TRUE) {
                ## the.set$Term <- gsub(".*\\(|\\)", "", the.set$Term)
                the.set$Term <- gsub(" \\(GO:.*","",the.set$Term)
            } else if (str_detect(db.in,"WikiPathways_") == TRUE) {
                the.set$Term <- gsub("_Homo.*","",the.set$Term)
            } else if (str_detect(db.in,"Reactome_") == TRUE) {
                the.set$Term <- gsub("_Homo.*","",the.set$Term)
            } else if (str_detect(db.in,"KEGG_")) {
                the.set$Term <- substr(the.set$Term,0,nchar(the.set$Term)-22)
            } else {}
        } else {
            the.set <- data.frame(Term="No Significant Results",P.value=0.1,Combined.Score=1,gene.counts=1)
        }
    } else {
        the.set <- data.frame(Term="NO SIGNIFICANT RESULTS",P.value=0.1,Combined.Score=1,gene.counts=1)
    }
    ## PLOT
    ggplot(data=the.set, aes(x=-log10(P.value), y=gene.counts, color=-log10(P.value), size=Combined.Score)) +
        geom_point(alpha=0.5,shape=1) +
        guides(size = "none", colour = "legend") +
        geom_text_repel(data=head(the.set,5), aes(label=Term),
                        size=1.5, min.segment.length = 0.4,
                        segment.size = 0.1, box.padding = 0.5,
                        ## label.padding = 0.1, 
                        ## nudge_y = 0.1,
                        ## nudge_x = 0.1,
                        direction = "both",
                        hjust = 1
                        ) +
        theme_bw() +
        xlab("-log10(P-Value)") +
        ylab("Gene Counts") +
        scale_x_continuous(breaks = seq(0,15,by=0.5)) +
        scale_y_continuous(breaks = seq(0,15,1)) +
        scale_fill_manual(values="white") +
        scale_color_gradient(low='black',high='red') +
        theme(panel.grid.major = element_blank(),
              legend.position="none",
              ## legend.position = "bottom",
              ## legend.justification = "top",
              ## legend.key = element_rect(fill = "white", size=0.1),
              ## legend.text = element_text(size = 3),
              ## legend.title = element_text(size = 4,face = 'bold'),
              ## legend.key.size = unit(0.25, "cm"),
              ## legend.margin = margin(t = 0, unit='cm'),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=5),
              axis.text = element_text(size=4),
              plot.title = element_text(size=7, hjust = 0.5),
              plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"))
}


## enriched.ARMS
go.biol.ARMS <- plot.enrichr(enriched.ARMS,"GO_Biological_Process_2018")
go.mol.ARMS <- plot.enrichr(enriched.ARMS,"GO_Molecular_Function_2018")
go.cell.ARMS <- plot.enrichr(enriched.ARMS,"GO_Cellular_Component_2018")
WIKI.ARMS <- plot.enrichr(enriched.ARMS,"WikiPathways_2016")
KEGG.ARMS <- plot.enrichr(enriched.ARMS,"KEGG_2016")
REACTOME.ARMS <- plot.enrichr(enriched.ARMS,"Reactome_2016")
pdf(paste0(figures.out.dir,'enrichment.plot.ARMS_Terms.pdf'), paper='a4', width = 0, height = 0, onefile=FALSE)
grid.arrange(go.biol.ARMS + ggtitle("GO Biological Process 2018"),
             go.mol.ARMS + ggtitle("GO Molecular Function 2018"),
             go.cell.ARMS + ggtitle("GO Cellular Component 2018"),
             WIKI.ARMS + ggtitle("WikiPathways 2016"),
             KEGG.ARMS + ggtitle("KEGG 2016"),
             REACTOME.ARMS + ggtitle("Reactome 2016"),
             ## + ggtitle(title.MILD.BACK),
             ## + ggtitle(title.SEVERE.MILD),
             layout_matrix = rbind(c(1,2), c(3,4), c(5,6)),
             widths = c(1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for WIKI Pathways of different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()


##enriched.ALL.COMBINED
## WIKI PATHWAYS
go.biol.ALL.COMBINED <- plot.enrichr(enriched.ALL.COMBINED,"GO_Biological_Process_2018")
go.mol.ALL.COMBINED <- plot.enrichr(enriched.ALL.COMBINED,"GO_Molecular_Function_2018")
go.cell.ALL.COMBINED <- plot.enrichr(enriched.ALL.COMBINED,"GO_Cellular_Component_2018")
WIKI.ALL.COMBINED <- plot.enrichr(enriched.ALL.COMBINED,"WikiPathways_2016")
KEGG.ALL.COMBINED <- plot.enrichr(enriched.ALL.COMBINED,"KEGG_2016")
REACTOME.ALL.COMBINED <- plot.enrichr(enriched.ALL.COMBINED,"Reactome_2016")

## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'enrichment.plot.COMBINED_Terms.pdf'), paper='a4', width = 0, height = 0, onefile=FALSE)
grid.arrange(go.biol.ALL.COMBINED + ggtitle("GO Biological Process 2018"),
             go.mol.ALL.COMBINED + ggtitle("GO Molecular Function 2018"),
             go.cell.ALL.COMBINED + ggtitle("GO Cellular Component 2018"),
             WIKI.ALL.COMBINED + ggtitle("WikiPathways 2016"),
             KEGG.ALL.COMBINED + ggtitle("KEGG 2016"),
             REACTOME.ALL.COMBINED + ggtitle("Reactome 2016"),
              ## + ggtitle(title.MILD.BACK),
              ## + ggtitle(title.SEVERE.MILD),
             layout_matrix = rbind(c(1,2), c(3,4), c(5,6)),
             widths = c(1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for WIKI Pathways of different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()



## WIKI PATHWAYS
WIKI.ALL <- plot.enrichr(enriched.ALL,"WikiPathways_2016")
WIKI.CASE.ARM <- plot.enrichr(enriched.CASE.ARM,"WikiPathways_2016")
WIKI.CASE.BACK <- plot.enrichr(enriched.CASE.BACK,"WikiPathways_2016")
WIKI.CASE <- plot.enrichr(enriched.CASE,"WikiPathways_2016")
WIKI.SEVERE.ARM <- plot.enrichr(enriched.SEVERE.ARM,"WikiPathways_2016")
WIKI.SEVERE.BACK <- plot.enrichr(enriched.SEVERE.BACK,"WikiPathways_2016")
WIKI.MILD.ARM <- plot.enrichr(enriched.MILD.ARM,"WikiPathways_2016")
WIKI.MILD.BACK <- plot.enrichr(enriched.MILD.BACK,"WikiPathways_2016")
WIKI.SEVERE.MILD <- plot.enrichr(enriched.SEVERE.MILD,"WikiPathways_2016")
## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'enrichment.plot.WIKI_Terms.pdf'), paper='a4r', width = 0, height = 0, onefile=FALSE)
grid.arrange(WIKI.ALL + ggtitle(title.ALL),
             WIKI.CASE.ARM + ggtitle(title.CASE.ARM),
             WIKI.CASE.BACK + ggtitle(title.CASE.BACK),
             WIKI.CASE + ggtitle(title.CASE),
             WIKI.SEVERE.ARM + ggtitle(title.SEVERE.ARM),
             WIKI.SEVERE.BACK + ggtitle(title.SEVERE.BACK),
             WIKI.MILD.ARM + ggtitle(title.MILD.ARM),
             WIKI.MILD.BACK + ggtitle(title.MILD.BACK),
             WIKI.SEVERE.MILD + ggtitle(title.SEVERE.MILD),
             layout_matrix = lay,
             widths = c(1,1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for WIKI Pathways of different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()

## REACTOME PATHWAYS
REACTOME.ALL <- plot.enrichr(enriched.ALL,"Reactome_2016")
REACTOME.CASE.ARM <- plot.enrichr(enriched.CASE.ARM,"Reactome_2016")
REACTOME.CASE.BACK <- plot.enrichr(enriched.CASE.BACK,"Reactome_2016")
REACTOME.CASE <- plot.enrichr(enriched.CASE,"Reactome_2016")
REACTOME.SEVERE.ARM <- plot.enrichr(enriched.SEVERE.ARM,"Reactome_2016")
REACTOME.SEVERE.BACK <- plot.enrichr(enriched.SEVERE.BACK,"Reactome_2016")
REACTOME.MILD.ARM <- plot.enrichr(enriched.MILD.ARM,"Reactome_2016")
REACTOME.MILD.BACK <- plot.enrichr(enriched.MILD.BACK,"Reactome_2016")
REACTOME.SEVERE.MILD <- plot.enrichr(enriched.SEVERE.MILD,"Reactome_2016")
## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'enrichment.plot.REACTOME_Terms.pdf'), paper='a4r', width = 0, height = 0, onefile=FALSE)
grid.arrange(REACTOME.ALL + ggtitle(title.ALL),
             REACTOME.CASE.ARM + ggtitle(title.CASE.ARM),
             REACTOME.CASE.BACK + ggtitle(title.CASE.BACK),
             REACTOME.CASE + ggtitle(title.CASE),
             REACTOME.SEVERE.ARM + ggtitle(title.SEVERE.ARM),
             REACTOME.SEVERE.BACK + ggtitle(title.SEVERE.BACK),
             REACTOME.MILD.ARM + ggtitle(title.MILD.ARM),
             REACTOME.MILD.BACK + ggtitle(title.MILD.BACK),
             REACTOME.SEVERE.MILD + ggtitle(title.SEVERE.MILD),
             layout_matrix = lay,
             widths = c(1,1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for REACTOME Pathways of different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()

## KEGG PATHWAYS
KEGG.ALL <- plot.enrichr(enriched.ALL,"KEGG_2016")
KEGG.CASE.ARM <- plot.enrichr(enriched.CASE.ARM,"KEGG_2016")
KEGG.CASE.BACK <- plot.enrichr(enriched.CASE.BACK,"KEGG_2016")
KEGG.CASE <- plot.enrichr(enriched.CASE,"KEGG_2016")
KEGG.SEVERE.ARM <- plot.enrichr(enriched.SEVERE.ARM,"KEGG_2016")
KEGG.SEVERE.BACK <- plot.enrichr(enriched.SEVERE.BACK,"KEGG_2016")
KEGG.MILD.ARM <- plot.enrichr(enriched.MILD.ARM,"KEGG_2016")
KEGG.MILD.BACK <- plot.enrichr(enriched.MILD.BACK,"KEGG_2016")
KEGG.SEVERE.MILD <- plot.enrichr(enriched.SEVERE.MILD,"KEGG_2016")
## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'enrichment.plot.KEGG_Terms.pdf'), paper='a4r', width = 0, height = 0, onefile=FALSE)
grid.arrange(KEGG.ALL + ggtitle(title.ALL),
             KEGG.CASE.ARM + ggtitle(title.CASE.ARM),
             KEGG.CASE.BACK + ggtitle(title.CASE.BACK),
             KEGG.CASE + ggtitle(title.CASE),
             KEGG.SEVERE.ARM + ggtitle(title.SEVERE.ARM),
             KEGG.SEVERE.BACK + ggtitle(title.SEVERE.BACK),
             KEGG.MILD.ARM + ggtitle(title.MILD.ARM),
             KEGG.MILD.BACK + ggtitle(title.MILD.BACK),
             KEGG.SEVERE.MILD + ggtitle(title.SEVERE.MILD),
             layout_matrix = lay,
             widths = c(1,1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for KEGG Pathways of different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()

## CELLULAR COMPARTMENTS
go.cell.ALL <- plot.enrichr(enriched.ALL,"GO_Cellular_Component_2018")
go.cell.CASE.ARM <- plot.enrichr(enriched.CASE.ARM,"GO_Cellular_Component_2018")
go.cell.CASE.BACK <- plot.enrichr(enriched.CASE.BACK,"GO_Cellular_Component_2018")
go.cell.CASE <- plot.enrichr(enriched.CASE,"GO_Cellular_Component_2018")
go.cell.SEVERE.ARM <- plot.enrichr(enriched.SEVERE.ARM,"GO_Cellular_Component_2018")
go.cell.SEVERE.BACK <- plot.enrichr(enriched.SEVERE.BACK,"GO_Cellular_Component_2018")
go.cell.MILD.ARM <- plot.enrichr(enriched.MILD.ARM,"GO_Cellular_Component_2018")
go.cell.MILD.BACK <- plot.enrichr(enriched.MILD.BACK,"GO_Cellular_Component_2018")
go.cell.SEVERE.MILD <- plot.enrichr(enriched.SEVERE.MILD,"GO_Cellular_Component_2018")
## CREATE LAYOUT
lay <- rbind(c(1,2,3), c(4,5,6), c(7,8,9))
## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'enrichment.plot.GO_Cellular_Terms.pdf'), paper='a4r', width = 0, height = 0, onefile=FALSE)
grid.arrange(go.cell.ALL + ggtitle(title.ALL),
             go.cell.CASE.ARM + ggtitle(title.CASE.ARM),
             go.cell.CASE.BACK + ggtitle(title.CASE.BACK),
             go.cell.CASE + ggtitle(title.CASE),
             go.cell.SEVERE.ARM + ggtitle(title.SEVERE.ARM),
             go.cell.SEVERE.BACK + ggtitle(title.SEVERE.BACK),
             go.cell.MILD.ARM + ggtitle(title.MILD.ARM),
             go.cell.MILD.BACK + ggtitle(title.MILD.BACK),
             go.cell.SEVERE.MILD + ggtitle(title.SEVERE.MILD),
             layout_matrix = lay,
             widths = c(1,1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for GO Cellular Compartment Terms",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()
## BIOLOGICAL PROCESSES
go.biol.ALL <- plot.enrichr(enriched.ALL,"GO_Biological_Process_2018")
go.biol.CASE.ARM <- plot.enrichr(enriched.CASE.ARM,"GO_Biological_Process_2018")
go.biol.CASE.BACK <- plot.enrichr(enriched.CASE.BACK,"GO_Biological_Process_2018")
go.biol.CASE <- plot.enrichr(enriched.CASE,"GO_Biological_Process_2018")
go.biol.SEVERE.ARM <- plot.enrichr(enriched.SEVERE.ARM,"GO_Biological_Process_2018")
go.biol.SEVERE.BACK <- plot.enrichr(enriched.SEVERE.BACK,"GO_Biological_Process_2018")
go.biol.MILD.ARM <- plot.enrichr(enriched.MILD.ARM,"GO_Biological_Process_2018")
go.biol.MILD.BACK <- plot.enrichr(enriched.MILD.BACK,"GO_Biological_Process_2018")
go.biol.SEVERE.MILD <- plot.enrichr(enriched.SEVERE.MILD,"GO_Biological_Process_2018")
## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'enrichment.plot.GO_Biological_Terms.pdf'), paper='a4r', width = 0, height = 0, onefile=FALSE)
grid.arrange(go.biol.ALL + ggtitle(title.ALL),
             go.biol.CASE.ARM + ggtitle(title.CASE.ARM),
             go.biol.CASE.BACK + ggtitle(title.CASE.BACK),
             go.biol.CASE + ggtitle(title.CASE),
             go.biol.SEVERE.ARM + ggtitle(title.SEVERE.ARM),
             go.biol.SEVERE.BACK + ggtitle(title.SEVERE.BACK),
             go.biol.MILD.ARM + ggtitle(title.MILD.ARM),
             go.biol.MILD.BACK + ggtitle(title.MILD.BACK),
             go.biol.SEVERE.MILD + ggtitle(title.SEVERE.MILD),
             layout_matrix = lay,
             widths = c(1,1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for GO Biological Terms of different comparisons",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()
## MOLECULAR FUNCTIONS
go.mol.ALL <- plot.enrichr(enriched.ALL,"GO_Molecular_Function_2018")
go.mol.CASE.ARM <- plot.enrichr(enriched.CASE.ARM,"GO_Molecular_Function_2018")
go.mol.CASE.BACK <- plot.enrichr(enriched.CASE.BACK,"GO_Molecular_Function_2018")
go.mol.CASE <- plot.enrichr(enriched.CASE,"GO_Molecular_Function_2018")
go.mol.SEVERE.ARM <- plot.enrichr(enriched.SEVERE.ARM,"GO_Molecular_Function_2018")
go.mol.SEVERE.BACK <- plot.enrichr(enriched.SEVERE.BACK,"GO_Molecular_Function_2018")
go.mol.MILD.ARM <- plot.enrichr(enriched.MILD.ARM,"GO_Molecular_Function_2018")
go.mol.MILD.BACK <- plot.enrichr(enriched.MILD.BACK,"GO_Molecular_Function_2018")
go.mol.SEVERE.MILD <- plot.enrichr(enriched.SEVERE.MILD,"GO_Molecular_Function_2018")
## PLOT ALL IN ONE PAGE  
pdf(paste0(figures.out.dir,'enrichment.plot.GO_Molecular_Terms.pdf'), paper='a4r', width = 0, height = 0, onefile=FALSE)
grid.arrange(go.mol.ALL + ggtitle(title.ALL),
             go.mol.CASE.ARM + ggtitle(title.CASE.ARM),
             go.mol.CASE.BACK + ggtitle(title.CASE.BACK),
             go.mol.CASE + ggtitle(title.CASE),
             go.mol.SEVERE.ARM + ggtitle(title.SEVERE.ARM),
             go.mol.SEVERE.BACK + ggtitle(title.SEVERE.BACK),
             go.mol.MILD.ARM + ggtitle(title.MILD.ARM),
             go.mol.MILD.BACK + ggtitle(title.MILD.BACK),
             go.mol.SEVERE.MILD + ggtitle(title.SEVERE.MILD),
             layout_matrix = lay,
             widths = c(1,1,1),
             heights = c(1,1,1)
             ## top=textGrob("Enrichment for GO Molecular Function Terms",
             ##              gp = gpar(fontsize = 12,fontface = 'bold',vjust = 1))
             )
dev.off()

## ------------------------------ Pathway Analysis ------------------------------#####
## Map ENSEMBL IDs to Entrez for pathway analysis
## Order by P-value

## Detach the PROPRER PACKAGE - INTERFERES WITH PATHWAY ANALYSIS
detach("package:PROPER", unload=TRUE)

data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs = kegg.sets.hs[sigmet.idx.hs]
columns(org.Hs.eg.db)

mapping <- data.frame(row.names=final_genes.hgnc)
mapping$ensembl <- mapIds(org.Hs.eg.db, keys=row.names(mapping), column="SYMBOL", keytype="SYMBOL", multiVals="first")
mapping$entrez <- mapIds(org.Hs.eg.db, keys=row.names(mapping), column="ENTREZID", keytype="SYMBOL", multiVals="first")
mapping$name   <- mapIds(org.Hs.eg.db, keys=row.names(mapping), column="GENENAME", keytype="SYMBOL", multiVals="first")


do_pathway <- function(res) {
    ## GET THE NAME OF THE COMPARISON SET WE ARE WORKING WITH FROM RES NAME
    path.name <- deparse(substitute(res))
    path.name <- substring(path.name,5,nchar(path.name))
    
    ## ORDER THE DIFFERENTIAL EXPRESSION RESULTS BY ADJUSTED PVALUE
    resPath <- res[order(res$padj),]
    
    ## GET THE ENTREZ ID'S FOR MAPPING TO KEGG PATHWAYS
    resPath$symbol = mapIds(org.Hs.eg.db, keys=row.names(resPath), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
    resPath$entrez = mapIds(org.Hs.eg.db, keys=row.names(resPath), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
    resPath$name =   mapIds(org.Hs.eg.db, keys=row.names(resPath), column="GENENAME", keytype="ENSEMBL", multiVals="first")
    
    ## EXTRACT FOLD-CHANGE INFORMATION FROM THE DE SET
    foldchanges = resPath$log2FoldChange
    names(foldchanges) = resPath$entrez
    
    ## RUN PATHWAY ANALYSIS WITH GAGE
    keggres <- gage(foldchanges, gsets=kegg.sets.hs, same.dir=TRUE) #, set.size=c(10,50))
    
    ## GET THE PATHWAYS UP-REGULARED
    select.greater <- keggres$greater[, "q.val"] < 0.05 & !is.na(keggres$greater[,"q.val"])
    up.reg <- rownames(keggres$greater)[select.greater]
    lapply(up.reg, function(x) print(paste0('upreg_',x)))
    
    ## GET THE PATHWAYS DOWN-REGULARED
    select.less <- keggres$less[, "q.val"] < 0.05 & !is.na(keggres$less[,"q.val"])
    down.reg <- rownames(keggres$less)[select.less]
    lapply(down.reg, function(x) print(paste0('downreg_',x)))
    
    ## GET THE TABLE FOR THE VALUES
    if (nrow(keggres$greater[select.greater,]) >=1) {
        ## print(head(keggres$greater[select.greater,]))
        path.table <- as.data.frame(keggres$greater[select.greater,])
        path.table$p.geomean <- formatC(path.table$p.geomean, format = "e", digits = 2)
        path.table$stat.mean <- formatC(path.table$stat.mean, format = "e", digits = 2)
        path.table$p.val <- formatC(path.table$p.val, format = "e", digits = 2)
        path.table$q.val <- formatC(path.table$q.val, format = "e", digits = 2)
        path.table$exp1 <- formatC(path.table$exp1, format = "e", digits = 2)
        print(path.table)
        write.csv(path.table,
                  file=paste0(path.out.dir,"pathway_table_",path.name,"_greater.csv"),
                  row.names=TRUE,col.names=TRUE)

        out.table = paste0(path.out.dir,"pathway_table_",path.name,"_greater.tex")
        table.head <- "\\begin{longtable}{ L{28em} L{5em} R{5em} R{5em} R{5em} R{5em} R{5em} }"
        table.caption <- paste0("\\caption[Up-regulated pathways in the ",path.name," comparison identified by \\gage{}]{\\textbf{Up-regulated pathways in the ",path.name," comparison identified by \\gage{}.}}")
        table.label <- paste0("\\label{tab:gage.up.",tolower(path.name),"}\\\\")
        ##
        write.table(table.head, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" & ", eol=" \n")
        write.table(table.caption, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
        write.table(table.label, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
        ##
        write.table("\\toprule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
        write.table(paste0("\\rowcolor{gray!25} Pathway name & ", paste(names(path.table),collapse=" & ")),
                    file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
        write.table("\\midrule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
        write.table(path.table, file=out.table, row.names=TRUE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
        write.table("\\bottomrule\n\\end{longtable}", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")

    }else {}
    
    if (nrow(keggres$less[select.less,]) >=1) {
        ## print(head(keggres$less[select.less,]))
        path.table <- as.data.frame(keggres$less[select.less,])
        path.table$p.geomean <- formatC(path.table$p.geomean, format = "e", digits = 2)
        path.table$stat.mean <- formatC(path.table$stat.mean, format = "e", digits = 2)
        path.table$p.val <- formatC(path.table$p.val, format = "e", digits = 2)
        path.table$q.val <- formatC(path.table$q.val, format = "e", digits = 2)
        path.table$exp1 <- formatC(path.table$exp1, format = "e", digits = 2)
        print(path.table)
        write.csv(path.table,
                  file=paste0(path.out.dir,"pathway_table_",path.name,"_less.csv"),
                  row.names=TRUE,col.names=TRUE)
        
        out.table = paste0(path.out.dir,"pathway_table_",path.name,"_less.tex")
        table.head <- "\\begin{longtable}{ L{28em} L{5em} R{5em} R{5em} R{5em} R{5em} R{5em} }"
        table.caption <- paste0("\\caption[Down-regulated pathways in the ",path.name," comparison identified by \\gage{}]{\\textbf{Down-regulated pathways in the ",path.name," comparison identified by \\gage{}.}}")
        table.label <- paste0("\\label{tab:gage.down.",tolower(path.name),"}\\\\")
        ##
        write.table(table.head, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" & ", eol=" \n")
        write.table(table.caption, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
        write.table(table.label, file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n")
        ##
        write.table("\\toprule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
        write.table(paste0("\\rowcolor{gray!25} Pathway name & ", paste(names(path.table),collapse=" & ")),
                    file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
        write.table("\\midrule", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
        write.table(path.table, file=out.table, row.names=TRUE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
        write.table("\\bottomrule\n\\end{longtable}", file=out.table, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    }else {}
    
    ## GET THE LIST
    keggresnames <- c(down.reg, up.reg)
    print(keggresnames)
    
    keggresids <- substr(c(up.reg, down.reg), 1, 8)
    keggresids <- substr(keggresnames, 1, 8)
    
    keggresnames.table <- data.frame(row.names=keggresids)
    keggresnames.table$names <- substr(keggresnames,10,nchar(keggresnames))
    
    print(keggresids)
    print(keggresnames.table)
    
    ## PLOT PATHWAYS USING PATHVIEW
    sapply(keggresids, function(pid) pathview(gene.data=foldchanges,
                                              pathway.id=pid,
                                              species="hsa" ,
                                              kegg.dir = kegg.out.dir,
                                              limit = list(gene = 10, cpd = 10),
                                              out.suffix = paste0(gsub(' ','.',keggresnames.table$names[match(pid, row.names(keggresnames.table))]),"_",path.name)
                                              )
           )

    ## MOVE PATHWAYS TO PATHWAY FOLDER
    system(paste0('mv *',path.name,".png ",'"',path.out.dir,'"'))
}

try(do_pathway(res.ALL), TRUE)
try(do_pathway(res.CASE.ARM), TRUE)
try(do_pathway(res.CASE.BACK), TRUE)
try(do_pathway(res.CASE), TRUE)
try(do_pathway(res.SEVERE.ARM), TRUE)
try(do_pathway(res.SEVERE.BACK), TRUE)
try(do_pathway(res.MILD.ARM), TRUE)
try(do_pathway(res.MILD.BACK), TRUE)
try(do_pathway(res.SEVERE.MILD), TRUE)

library(PROPER)

## write(utils::toLatex(sessionInfo()), file=paste0(out.dir,"session.info.tex"))
## write(sessioninfo::session_info(), file="session.info.tex")
## clipr::write_clip(sessioninfo::session_info())

## ============================== THIS SECTION IS FOR FUNCTIONS!! ============================== ##
##
## ===== THE FUNCTION: FOR GETTING SIGNIFICANT GENES
get_sig <- function(res,fdr,lfc){
    sig <- res[!is.na(res$padj) & res$padj < fdr & abs(res$log2FoldChange) >= lfc,]
    return(sig)
}
## ===== THE FUNCTION: ENDS HERE!

## ===== THE FUNCTION: GET THE NAMES OF THE SIGNIFICANT GENES FROM ENSEMBL
get_ensembl_annotations <- function(selected) {
    gene_names <- getBM(attributes=c('ensembl_gene_id', 'description', 'gene_biotype', 'chromosome_name', 'hgnc_symbol'),
                        filters = 'ensembl_gene_id',
                        values = selected,
                        mart = ensembl,
                        uniqueRows = TRUE,
                        verbose=FALSE)
    gene_names <- as.data.frame(gene_names)
    colnames(gene_names) <- c('ensembl_id','description', 'biotype', 'chromosome', 'HGNC')
    gene_names$description <- sapply(gene_names$description, function(x) str_replace(x, "\\[.*\\]", ""))
    return(gene_names)
}

get_gene_names <- function(selected) {
    gene_names <- ensembl.annotations[which(ensembl.annotations$ensembl_id %in% selected),]
    return(gene_names)
}

## ===== THE FUNCTION: ENDS HERE!

## ===== THE FUNCTION: FOR PLOTTING REJECTIONS
plot_rejections_plot <- function(res) {
    plot(metadata(res)$filterNumRej, type="b", ylab="number of rejections", xlab="quantiles of filter")
    lines(metadata(res)$lo.fit, col="red")
    abline(v=metadata(res)$filterTheta)
}
## ===== THE FUNCTION: ENDS HERE!

## ===== THE FUNCTION: FOR PLOTTING TRANSFORMATIONS
plot_transformations <- function(dds, rld, vsd) {
    dds_est <- estimateSizeFactors(dds)
    df <- bind_rows(as_data_frame(log2(counts(dds_est, normalized=TRUE)[, 1:2]+1))
                    %>% mutate(transformation = "log2(x + 1)"), as_data_frame(assay(rld)[, 1:2])
                    %>% mutate(transformation = "rlog"), as_data_frame(assay(vsd)[, 1:2])
                    %>% mutate(transformation = "vst"))
    
    colnames(df)[1:2] <- c("x", "y")
    
    ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 100) +
        coord_fixed() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=7),
              axis.text = element_text(size=5),
              plot.title = element_text(size=8, hjust = 0.5)) +
        facet_grid( . ~ transformation)
}
## ===== THE FUNCTION: ENDS HERE!

## ===== THE FUNCTION: FOR PLOTTING BOXPLOTS
plot_boxplots <- function(dds, title) {
    boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, main=title, cex.axis=1 , cex.names=1)
}
## THE FUNTION: ENDS HERE!

## ===== THE FUNCTION: FOR PLOTTING SAMPLE-DISTANCE MATRICES
plot_distance_matrix <- function(rld) {
    #mycols <- brewer.pal(8, "Dark2")[1:length(unique(condition))]
    sampleDists <- dist(t(assay(rld)))
    sampleDistMatrix <- as.matrix( sampleDists )
    rownames(sampleDistMatrix) <- paste( rownames(sampleDistMatrix), rld$site, sep = " - " )
    colnames(sampleDistMatrix) <- NULL #paste( rownames(sampleDistMatrix), rld$age, sep = " - " )
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    pheatmap(sampleDistMatrix,
             clustering_distance_rows = sampleDists,
             clustering_distance_cols = sampleDists,
             col = colors)
}
## ===== THE FUNCTION: ENDS HERE!

## ===== THE FUNCTION: FOR PRINCIPAL COMPONENT ANALYSIS (GETTING THE PCA DATA - MYPLOTPCA, THEN PLOTTING PCA - PLOT_PCA)
myPlotPCA <- function (object, intgroup = c('condition','site','age','disease.duration'), ntop = 500, returnData = FALSE) {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    
    intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
    group <- if (length(intgroup) > 1) {
                 factor(apply(intgroup.df, 1, paste, collapse = " : "))
             }
             else {
                 colData(object)[[intgroup]]
             }

    d <- data.frame(PC1 = pca$x[, 1],
                    PC2 = pca$x[, 2],
                    PC3 = pca$x[, 3],
                    PC4 = pca$x[, 4],
                    group = group, intgroup.df, name = colnames(object))
    
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
}

plot_pca <- function(rld) {
    pcaData <- myPlotPCA(rld, intgroup = c('condition','site','age','gender','disease.duration'), returnData = TRUE)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    ggplot(pcaData, aes(x=PC1, y=PC2, color=condition, shape=site, size=disease.duration)) +
        geom_point(alpha=0.4) +
        theme_bw() +
        scale_color_manual(values = c(case="red", ctrl="blue")) +
        ##scale_size_continuous(limits = c(1,10), guide="none") +
        scale_size_manual(values=c('0'=1,'12'=1.5,'24'=2,'36'=2.5,'48'=3,'60'=3.5), guide="none") +
        geom_text_repel(aes(label=paste(pcaData$name,'(', pcaData$age, ')', sep='')),
                        size=1.5, segment.size = 0.1) +
        xlab(paste0("PC1: ", percentVar[1], "% variance")) +
        ylab(paste0("PC2: ", percentVar[2], "% variance")) +
        scale_shape_manual(values=c(arm=16,back=1)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              legend.position = "bottom",
              legend.justification = "top",
              legend.key = element_rect(fill = "gray", size=0.1),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 7,face = 'bold'),
              legend.key.size = unit(0.25, "cm"),
              legend.margin = margin(t = 0, unit='cm'),
              axis.title = element_text(size=6),
              axis.text = element_text(size=5),
              plot.title = element_text(size=7, hjust = 0.5))
}
## ===== THE FUNCTION: ENDS HERE!

get_legend<-function(plot){
    tmp <- ggplot_gtable(ggplot_build(plot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
}

## ===== THE FUNCTION: FOR PLOTTING VOLCANO PLOTS USING GGPLOT - BETTER THAN THE OTHER ONE!
my_volcano <- function(res, lfc.cut = 2, sig.cut = 0.01) {
    raw.set <- eval(as.name(paste0('dds',substr(deparse(substitute(res)),4,nchar(deparse(substitute(res)))))))
    filtered.set <- deparse(substitute(res))

    ## RAW SET
    raw.set <- as.data.frame(results(raw.set))
    raw.set <- raw.set[,colnames(raw.set) %in% c("log2FoldChange","pvalue","padj")]
    raw.set$log.pvalue <- -log10(raw.set[["pvalue"]])
    
    ## THE DATASET
    mydata <- as.data.frame(res)[,colnames(res) %in% c("log2FoldChange","pvalue","padj")]
    mydata$log.pvalue <- -log10(mydata[["pvalue"]])
    mydata$sig.thr <- ifelse(mydata$padj < sig.cut, 1, 0)
    mydata$lfc.thr <- ifelse(abs(mydata$log2FoldChange) > lfc.cut, 1, 0)
    mydata$both.thr <- ifelse(mydata$padj < sig.cut & abs(mydata$log2FoldChange) > lfc.cut, 1, 0)

    ## PVALUE CUTOFF
    ## pval.cut <- max(na.omit(mydata[which(mydata$both.thr == 1),])$pvalue)
    pval.cut <- max(na.omit(raw.set[which(raw.set$padj < sig.cut & abs(raw.set$log2FoldChange) > lfc.cut),])$pvalue) # & row.names(raw.set) %in% eval(as.name(set))),]
    ## pval.cut <- sig.cut
    
    ## THE SIGNIFICANT SET
    set <- paste0('selected', substring(filtered.set,4,nchar(filtered.set)))
    
    ## THE COMPARISON NAME
    compare_name <- paste0('compare', substring(filtered.set,4,nchar(filtered.set)))
    
    ## THE MOST HIGHLY EXPRESSED TOP 10 GREEN
    mydata.top <- mydata[which(mydata$both.thr == 1 & row.names(mydata) %in% eval(as.name(set))),]
    mydata.top <- row.names(head(mydata.top[order(mydata.top$padj),],10))

    names <- get_gene_names(mydata.top)$HGNC
    print(names)
    print(mydata.top)

    green  <- row.names(raw.set[which(raw.set$padj < sig.cut & abs(raw.set$log2FoldChange) > lfc.cut & row.names(raw.set) %in% eval(as.name(set))),])
    blue   <- row.names(raw.set[which(raw.set$padj < sig.cut & abs(raw.set$log2FoldChange) > lfc.cut & !row.names(raw.set) %in% eval(as.name(set))),])
    red    <- row.names(raw.set[which(raw.set$padj < sig.cut & abs(raw.set$log2FoldChange) < lfc.cut),])
    gray   <- row.names(raw.set[which(raw.set$padj > sig.cut & abs(raw.set$log2FoldChange) < lfc.cut),])
    orange <- row.names(raw.set[which(raw.set$padj > sig.cut & abs(raw.set$log2FoldChange) > lfc.cut),])

    print(length(green))
    print(length(blue))
    print(length(red))
    print(length(gray))
    print(length(orange))
    
    ggplot(data=raw.set[which(row.names(raw.set) %in% green), ], aes(x=log2FoldChange, y=log.pvalue)) +
        geom_point(size=0.2, color='green') +
        geom_label_repel(data=raw.set[which(row.names(raw.set) %in% mydata.top),],
                         aes(label=names),#row.names(raw.set)[which(row.names(raw.set) %in% mydata.top)]),
                         size=1, min.segment.length = 0.2, segment.size = 0.1, label.padding = 0.1, box.padding = 0.1) +
        geom_point(data=raw.set[which(row.names(raw.set) %in% red),], aes(x=log2FoldChange, y=log.pvalue), colour="red", size=0.2) +
        geom_point(data=raw.set[which(row.names(raw.set) %in% gray),], aes(x=log2FoldChange, y=log.pvalue, colour = "gray"), size=0.2) +
        geom_point(data=raw.set[which(row.names(raw.set) %in% orange),], aes(x=log2FoldChange, y=log.pvalue), colour="orange", size=0.2) +
        geom_point(data=raw.set[which(row.names(raw.set) %in% blue),], aes(x=log2FoldChange, y=log.pvalue), colour="blue", size=0.2) +        
        scale_color_identity() +
        theme_bw() +
        annotate("text", x = -Inf, y = Inf, hjust = 0, vjust = 1, label = paste(unlist(eval(as.name(compare_name))[[1]]), collapse = " "), size=1.6) +
        annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1, label = paste(unlist(eval(as.name(compare_name))[[2]]), collapse = " "), size=1.6) +
        annotate("rect", xmin = -lfc.cut, xmax = lfc.cut, ymin = -log10(pval.cut), ymax = Inf, alpha = 0.1) +
        annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -log10(pval.cut), alpha = 0.1) +
        geom_vline(aes(xintercept=-lfc.cut), size=0.2, alpha=0.5, linetype="dashed") +
        geom_vline(aes(xintercept=lfc.cut), size=0.2, alpha=0.5, linetype="dashed") +
        geom_hline(aes(yintercept=-log10(pval.cut)), size=0.2, alpha=0.5, linetype="dashed") +
        xlab("log2(Fold-Change)") +
        ylab("-log10(P-Value)") +
        scale_x_continuous(breaks = seq(-50,50,by=5)) +
        scale_y_continuous(breaks = seq(0,50,5)) +
        theme(panel.grid.major = element_blank(),
              legend.position="none",
              panel.grid.minor = element_blank(),
              axis.title = element_text(size=5),
              axis.text = element_text(size=4),
              plot.title = element_text(size=7, hjust = 0.5),
              plot.margin = unit(c(0.5, 0.25, 0.25, 0.25), "cm"))
}
## ===== THE FUNCTION: ENDS HERE

## volcano.ALL <- my_volcano(res.ALL)

## ===== THE FUNCTION: FOR CREATING VOLCANO PLOT (VISUAL REPRESENTATION OF SIGNIFICANT GENES IN EACH COMPARISON DATASET)
volcanoplot <- function (res, lfcthresh=0, sigthresh=0.01, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
    with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
    with(subset(res, padj < sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
    with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
    if (labelsig) {
        require(calibrate)
        with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh),
             textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
    }
    legend(legendpos, xjust=1, yjust=1,
           legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"),
           pch=20, col=c("red","orange","green"))
}
## ===== THE FUNCTION: ENDS HERE!

## ===== THE FUNCTION: FOR PLOTTING HEATMAPS (TAKES RLD, SELECTED GENES, GENES KNONT TO BE INVOLVED, GENE NAMES AND TITLE)
plot_heat <- function(rld,selected,known_in_my_data,gene_names,title) {
    ## The data
    mat  <- assay(rld)[ selected, ]
    mat  <- mat - rowMeans(mat)
    
    ## For the colors
    the_col <- data.frame(target = rep('No', nrow(mat)), stringsAsFactors=F)
    rownames(the_col) <- rownames(mat)
    the_col$target[rownames(the_col) %in% known_in_my_data$ensembl_id] <- "Yes"
    ## anno <- as.data.frame(colData(rld)[, 1,drop=FALSE])
    anno <- as.data.frame(colData(rld)[, c('condition','site','gender','lung','myositis','anti.scl'),
                                       drop=FALSE])
    
    rownames(mat) <- paste(ifelse(gene_names$HGNC[match(rownames(mat),gene_names$ensembl_id)] == "",
                                  rownames(mat),
                                  gene_names$HGNC[match(rownames(mat),gene_names$ensembl_id)]),
                           gene_names$chromosome[match(rownames(mat),gene_names$ensembl_id)],
                           gene_names$description[match(rownames(mat),gene_names$ensembl_id)],
                           sep="\t-\t")
    
    rownames(the_col) <- paste(ifelse(gene_names$HGNC[match(rownames(the_col),gene_names$ensembl_id)] == "",
                                      rownames(the_col),
                                      gene_names$HGNC[match(rownames(the_col),gene_names$ensembl_id)]),
                               gene_names$chromosome[match(rownames(the_col),gene_names$ensembl_id)],
                               gene_names$description[match(rownames(the_col),gene_names$ensembl_id)],
                               sep="\t-\t")
   
    ann_colors = list(
        lung = c('1' = "red", '0'="green"),
        myositis = c('1' = "red", '0'="green"),
        anti.scl = c('1' = "red", '0'="green"),
        condition = c(case = "red", ctrl="green"),
        gender = c(F = "pink", M="blue"),
        site = c(arm = "springgreen4", back = "gray80"),
        target = c(Yes = "green4", No = "cornsilk")
    )
    
    pheatmap(mat, annotation_col = anno,  
             fontsize = 6, fontsize_row = 4,
             color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
             border_color = NA, treeheight_col = 30, 
             treeheight_row = 35,fontsize_col = 5, 
             show_rownames = TRUE, cellwidth = 6,
             cellheight = 3.5, cluster_cols = TRUE,
             clustering_distance_rows = "euclidean",
             annotation_row = the_col,
             annotation_names_row = TRUE,
             annotation_legend = FALSE,
             annotation_colors = ann_colors)
    ## ,
    ## main=title)
}

height <- function(genes_len) {
    length(genes_len$ensembl_id) * x + y
}
## ===== THE FUNCTION: ENDS write!

write.latex.table <- function(table,columns,extra=NULL,out.file) {
    write.table(columns, file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" & ", eol=" \n")
    write.table("\\toprule", file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    lapply(extra, function(row) write.table(row, file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n"))
    write.table(paste0("\\rowcolor{gray!25} ", paste(names(table),collapse=" & ")), file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
    write.table("\\midrule", file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    write.table(table, file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
    write.table("\\bottomrule\n\\end{tabular}", file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
}

write.long.table <- function(table,columns,extra=NULL,out.file) {
    write.table(columns, file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" & ", eol=" \n")
    write.table("\\toprule", file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    lapply(extra, function(row) write.table(row, file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE,sep=" & ", eol=" \n"))
    write.table(paste0("\\rowcolor{gray!25} ", paste(names(table),collapse=" & ")), file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
    write.table("\\midrule", file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
    write.table(table, file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \\\\\n")
    write.table("\\bottomrule\n\\end{tabular}", file=out.file, row.names=FALSE, col.names=FALSE, quote=FALSE, append=TRUE, sep=" & ", eol=" \n")
}

## ============================= FUNCTIONS END HERE!!! ============================= ##

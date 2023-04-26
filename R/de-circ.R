#' @import limma
#' @importFrom edgeR DGEList
#' @export
de_circrna = function (circ_mtx, condition, design, contr, 
                       min_num_samples=2, lfc=1, p.value=0.05, 
                       voomWithQualityWeights=FALSE, plotVoomResult=FALSE, ...) {
    #--------------------------------------------------------------------------------
    # DE circRNA 
    #--------------------------------------------------------------------------------
    # Filter by min_num_samples
    circ_mtx <- circ_mtx[rowSums(circ_mtx>0)>=min_num_samples,]

    # Create DGE
    circ_DGE <- edgeR::DGEList(counts = circ_mtx, group = condition)

    # Condition
    if (voomWithQualityWeights) {
        v <- limma::voomWithQualityWeights(circ_DGE, design, plot=plotVoomResult)
    } else {
        v <- limma::voom(circ_DGE, design, plot=plotVoomResult)
    }

    vfit <- limma::lmFit(v, design)
    vfit <- limma::contrasts.fit(vfit, contrasts=contr)
    efit <- limma::eBayes(vfit)
    circ_de <- limma::topTable(efit, p.value=p.value, lfc=lfc, ...)
    circ_de
}

.test_de_circrna = function () {
    #min_num_samples = 2
    #plotVoomResult=FALSE
    #p.value=0.05
    #lfc=1
    sampleFile <- system.file("extdata","samples_U2932_PQR.txt", package="EasyCirc")
    condition <- factor(rep(c("DMSO", "PQR"),3))

    circ <- read_ciri_output(sampleFile)
    names(circ)
    circ_df <- circ$circ_df 
    circ_mtx <- circ$circ_mtx
    head(circ_df)
    head(circ_mtx)
    nrow(circ_mtx)
    nrow(circ_df)
    design <- model.matrix(~0+condition)
    contr <- limma::makeContrasts("PQRvsDMSO"  = conditionPQR - conditionDMSO, levels = colnames(design))
    circ_de <- de_circrna(circ_mtx, condition, design, contr, lfc=0, p.value =0.3 ,voomWithQualityWeights=TRUE)
    nrow(circ_de)
    head(circ_de)
    circ_mtx[rownames(circ_de),]

}

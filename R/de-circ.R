#' @title EasyCircR - DE circRNAs
#'
#' @description Differential expressed analysis of detected circRNAs performed with \code{limma} package.
#' 
#' @author Luca Parmigiani, Antonino Aparo, Simone Avesani
#' 
#' @param circ_mtx the count matrix of circRNAs as results of \code{EasyCircR::read_ciri_output()}.
#'
#' @param condition \code{numeric}, \code{character} or \code{factor} with two levels
#'   that groups the samples into two conditions.
#' @param design design matrix with rows corresponding to samples and columns to coefficients 
#' to be estimated.
#' @param contr contrast matrix computed by \code{limma::makeContrasts()}.
#' @param min_num_samples removing all the circRNAs that are not present in at least the specified number of samples. Default is \code{2}.
#' @param lfc minimum absolute log2-fold-change required. 
#' @param p.value Benjamini-Hochberg adjusted p-values cutoff. Only genes with lower p-values are listed. Default is \code{0.05}.
#' @param voomWithQualityWeights \code{logical}. If \code{TRUE} combine voom observational-level weights with sample-specific quality weights in a designed experiment.
#' @param plotVoomResult \code{logical}. If \code{TRUE} the plot of the mean-variance trend and sample-specific weights is displayed.
#' @param ... other arguments are passed to \code{limma:topTable()}.
#' 
#' @return A \code{dataframe} storing for each circRNA following features: \describe{
#' \item{logFC}{estimate of the log2-fold-change corresponding to the effect or contrast}
#' \item{AveExpr}{average log2-expression of circRNA in all the samples}
#' \item{t}{moderated t-statistic}
#' \item{P.Value}{raw p-value}
#' \item{adj.P.Value}{adjusted Benjamini-Hochberg p-value}
#' \item{B}{log-odds that the circRNA is differential expressed}}
#' 
#' @examples
#' design <- model.matrix(~0+condition)
#' contr <- limma::makeContrasts("PQRvsDMSO"  = conditionPQR - conditionDMSO, levels = colnames(design))
#' 
#' circ_de <- de_circrna(circ_mtx, condition, design, contr, lfc=0, p.value=0.05, 
#'                       voomWithQualityWeights=FALSE, min_num_samples=2)
#'
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

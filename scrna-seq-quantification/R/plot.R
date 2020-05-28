plotInfReps <- function(y, idx, x, cov=NULL,
                        colsDrk=c("dodgerblue","goldenrod4","royalblue4",
                                  "red3","purple4","darkgreen"),
                        colsLgt=c("lightblue1","goldenrod1","royalblue1",
                                  "salmon1","orchid1","limegreen"),
                        xaxis, xlab, ylim,
                        main, mainCol,
                        legend=FALSE,
                        legendPos="topleft",
                        legendTitle=FALSE,
                        legendCex=1,
                        useMean=TRUE,
                        applySF=FALSE,
                        reorder,
                        thin) {
  
  hasInfReps <- any(grepl("infRep", assayNames(y)))
  # logical switch for if variance assay is present
  # if not, just show the count or mean
  showVar <- TRUE
  if (!hasInfReps) {
    showVar <- "variance" %in% assayNames(y)
    if (useMean & !("mean" %in% assayNames(y))) {
      message("using 'counts' assay, as 'mean' is missing, see argument 'useMean'")
      useMean <- FALSE
    }
  }
  stopifnot(x %in% names(colData(y)))
  condition <- colData(y)[[x]]
  # whether x is a factor variable or numeric (e.g. pseudotime)
  xfac <- is(condition, "factor")
  if (!xfac) stopifnot(is(condition, "numeric"))
  stopifnot(length(colsDrk) == length(colsLgt))
  if (xfac) {
    ncond <- nlevels(condition)
    stopifnot(ncond <= length(colsDrk))
    colsDrk <- colsDrk[seq_len(ncond)]
    colsLgt <- colsLgt[seq_len(ncond)]
  } else {
    if (hasInfReps) {
      message("boxplot of inf. reps over numeric 'x' not supported yet")
      hasInfReps <- FALSE
      stopifnot("variance" %in% assayNames(y))
    }
  }
  if (missing(xaxis)) {
    if (xfac) {
      xaxis <- ncol(y) < 30
    } else {
      xaxis <- TRUE
    }
  }
  if (missing(thin)) {
    thin <- if (ncol(y) >= 400) 2 else if (ncol(y) >= 150) 1 else 0
  } else {
    stopifnot(thin >= 0 & thin <= 2)
  }
  if (!is.null(cov)) {
    stopifnot(cov %in% names(colData(y)))
    covariate <- factor(colData(y)[[cov]])
    stopifnot(is(covariate, "factor"))
    ngrp <- nlevels(covariate)
  }
  infRepsScaled <- FALSE
  if (!is.null(metadata(y)$infRepsScaled)) {
    infRepsScaled <- metadata(y)$infRepsScaled
  }
  # single cell?
  sc <- FALSE 
  if (!is.null(metadata(y)$tximetaInfo$type)) {
    if (metadata(y)$tximetaInfo$type == "alevin") {
      sc <- TRUE
    }
  }
  if (missing(xlab)) {
    if (xfac) {
      xlab <- if (sc) "cells" else "samples"
    } else {
      xlab <- x
    }
  }
  if (missing(reorder)) {
    if (xfac) {
      reorder <- sc
    } else {
      reorder <- FALSE
    }
  } else {
    if (!xfac & reorder) stop("reorder not used when 'x' is numeric")
  }
  ylab <- if (infRepsScaled) "scaled counts" else "counts"
  # this is a dummy variable used when making the plot()
  # if we don't put x-axis ticks, then we will move
  # the label up higher using title()
  xlabel <- if (xaxis) xlab else ""
  if (missing(main)) {
    if (missing(mainCol)) {
      if (is.character(idx)) {
        main <- idx
      } else {
        main <- rownames(y)[idx]
      }
    } else {
      stopifnot(mainCol %in% names(mcols(y)))
      main <- mcols(y)[idx,mainCol]
    }
  }
  if (reorder) {
    if (hasInfReps) {
      infReps <- assays(y[idx,])[grep("infRep",assayNames(y))]
      value <- colMeans(unlist(infReps))
    } else {
      which.assay <- if (useMean) "mean" else "counts"
      value <- assays(y)[[which.assay]][idx,]
      if (applySF & !is.null(y$sizeFactor)) {
        value <- value/y$sizeFactor
      }
    }
    if (is.null(cov)) {
      o <- order(condition, value)
    } else {
      o <- order(covariate, condition, value)
    }
    # not reordering:
  } else {
    if (xfac) {
      if (is.null(cov)) {
        o <- order(condition)
      } else {
        o <- order(covariate, condition)
      }
      # numeric 'x', don't reorder
    } else {
      o <- seq_along(condition)
    }
  }
  # col - dark color for boxplot border, points
  # col.hglt - highlight color for inside boxplot, lines when n >= 400
  if (xfac) {
    if (is.null(cov)) {
      samp.nums <- unlist(lapply(table(condition), seq_len))
      col <- rep(colsDrk, table(condition))
      col.hglt <- rep(colsLgt, table(condition))
    } else {
      vec.tab <- as.vector(table(condition, covariate))
      samp.nums <- unlist(lapply(vec.tab, seq_len))
      col <- rep(rep(colsDrk, ngrp), vec.tab)
      col.hglt <- rep(rep(colsLgt, ngrp), vec.tab)
    }
  } else {
    if (is.null(cov)) {
      col <- colsDrk[1]
      col.hglt <- colsLgt[1]
    } else {
      col <- colsDrk[covariate]
      col.hglt <- colsLgt[covariate]
    }
  }
  ### boxplot ###
  if (hasInfReps) {
    infReps <- assays(y[idx,])[grep("infRep",assayNames(y))]
    cts <- unlist(infReps)[,o]
    ymax <- max(cts)
    ymin <- if (is.null(cov)) 0 else -0.02 * ymax
    if (missing(ylim)) {
      ylim <- c(ymin,ymax)
    } else {
      stopifnot(length(ylim) == 2)
    }
    boxplot2(cts, col=col, col.hglt=col.hglt, ylim=ylim,
             xlab=xlabel, ylab=ylab, main=main)
    ### point and line plot by 'x' levels or numeric 'x' ###
  } else {
    which.assay <- if (useMean) "mean" else "counts"
    cts <- assays(y)[[which.assay]][idx,o]
    if (showVar) {
      sds <- sqrt(assays(y)[["variance"]][idx,o])
      Q <- qnorm(.975)
    }
    if (applySF & !is.null(y$sizeFactor)) {
      cts <- cts / y$sizeFactor[o]
      if (showVar) {
        sds <- sds / y$sizeFactor[o]
      }
      ylab <- "scaled counts"
    }
    ymax <- if (showVar) max(cts + Q*sds) else max(cts)
    ymin <- 0
    if (xfac & !is.null(cov)) {
      ymin <- -0.02 * ymax
    }
    if (missing(ylim)) {
      ylim <- c(ymin, ymax)
    } else {
      stopifnot(length(ylim) == 2)
    }
    if (xfac) {
      plot(cts, type="n", main=main,
           xaxt="n", ylim=ylim,
           xlab=xlabel, ylab=ylab)
    } else {
      plot(condition, cts, type="n", main=main,
           xaxt="n", ylim=ylim,
           xlab=xlabel, ylab=ylab)
    }
    seg.lwd <- if (thin == 0) 2 else if (thin == 1) 1 else 3
    seg.col <- if (thin < 2) col else col.hglt
    if (showVar) {
      if (xfac) {
        segments(seq_along(cts), pmax(cts - Q*sds, 0),
                 seq_along(cts), cts + Q*sds,
                 col=seg.col, lwd=seg.lwd)
      } else {
        segments(condition, pmax(cts - Q*sds, 0),
                 condition, cts + Q*sds,
                 col=seg.col, lwd=seg.lwd)
      }
    }
    pts.pch <- if (thin == 0) 22 else 15
    pts.lwd <- if (thin == 0) 1. else 1
    pts.cex <- if (thin == 0) 1 else 0.5
    if (xfac) {
      points(cts,
             col=col, pch=pts.pch, bg=col.hglt,
             cex=pts.cex, lwd=pts.lwd)
    } else {
      points(condition, cts,
             col=col, pch=pts.pch, bg=col.hglt,
             cex=pts.cex, lwd=pts.lwd)
    }
  }
  if (xaxis) {
    if (xfac) {
      axis(1, seq_along(condition), samp.nums)
    } else {
      axis(1)
    }
  }
  if (!xaxis) {
    title(xlab=xlab, mgp=c(1,1,0))
  }
  if (xfac & !is.null(cov)) {
    cuts <- cumsum(table(covariate))
    segments(c(1,cuts[-ngrp]+1),ymin,cuts,ymin,lwd=3,
             col=rep(c("black","grey60"),length=ngrp))
  }
  if (legend & (xfac | !is.null(cov))) {
    group <- if (xfac) condition else covariate
    group.name <- if (xfac) x else cov
    ltitle <- if (legendTitle) group.name else NULL
    legend(legendPos, legend=levels(group), title=ltitle,
           col=colsDrk, pt.bg=colsLgt, pch=22,
           cex=legendCex, bg="white", box.col=NA, inset=.01)
  }
}
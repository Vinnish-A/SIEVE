metaIntegration = R6Class(
  classname = "metaIntegration",
  public = list(
    rules   = NULL,
    exprs   = NULL,
    phen    = NULL,
    omics   = NULL,
    target  = NULL,
    aligned = NULL,
    initialize = function(intes_, target_ = NULL) {

      if (is.null(names(intes_))) stop("Integrations should be presented in a named list.")

      self$exprs = setNames(lapply(intes_, \(ele_) ele_$expr), names(intes_))
      self$phen  = intes_[[1]]$phen
      self$omics = names(self$exprs)

      self$target  = target_
      self$aligned = F

      self$rules = list(
        .getSamples = function (phen_, condition_) {
          conditions_ = c("default", "pairwise", "faraway")
          if (!condition_ %in% conditions_) {
            stop("Condition shoule be limited in default, pairwise or faraway")
          }

          if (condition_ == "default") {
            samples_  = phen_[grepl("01$", sample), list(sample)][[1]]
          } else if (condition_ == "pairwise") {
            patients_ = phen_[grepl("01$|11$", sample), .N, by = patient][N>1][[1]]
            samples_  = phen_[patient %in% patients_ & grepl("01$|11$", sample)][[1]]
          } else if (condition_ == "faraway") {
            patients_ = phen_[grepl("01$|06$", sample), .N, by = patient][N>1][[1]]
            samples_  = phen_[patient %in% patients_ & grepl("01$|06$", sample)][[1]]
          }

          return(samples_)
        }
      )

      private$.checkIt()
      cat("Prepared\n")

    },
    setTarget = function (target_) {

      self$target = target_

    },
    extract = function(feas_, condition_ = "default", fixPhen_ = T) {

      results_ = lapply(feas_, private$.extract, condition_ = condition_)
      names(results_) = feas_

      results_ = results_[!sapply(results_, is.null)]

      if (length(results_) == 0) {
        stop("Non valid features captured!")
      }

      results_ = lapply(results_, \(result__) result__[!private$.allNaRowsOf(result__)])

      return(geneFeature$new(results_, self$phen, names(self$exprs), self$target, fixPhen_))

    },
    align = function(omics_ = self$omics) {

      self$omics = omics_
      self$exprs = self$exprs[self$omics]

      feaToKeep_ = Reduce(\(ele1__, ele2__) intersect(ele1__, ele2__), lapply(self$exprs, colnames))
      smpToKeep_ = Reduce(\(ele1__, ele2__) intersect(ele1__, ele2__), lapply(self$exprs, `[[`, "sample"))

      self$exprs = lapply(self$exprs, \(expr__) expr__[smpToKeep_, ..feaToKeep_])
      self$phen  = self$phen[smpToKeep_]

      self$aligned = T

    },
    apply = function(fun_, omics_ = self$omics, condition_ = "default", target_ = self$target, enableProgressBar_ = F) {

      eachOmic_ = list()
      for (omic_ in omics_) {

        cat("Calculating ", omic_, "...\n", sep = "")
        expr_ = private$.data(omic_, condition_, target_)

        targetFun_ = function (ind__, data__) {
          fun_(data__[[2]], data__[[ind__]])
        }

        result_ = private$.apply(targetFun_, expr_, enableProgressBar_ = enableProgressBar_)

        eachValue_ = vector("list", length(result_[[1]]))
        for (i in seq_along(result_[[1]])) {
          eachValue_[[i]] = sapply(result_, `[`, i)
        }
        eachValue_ = as.data.table(eachValue_)

        if (!is.null(result_[[1]])) names(eachValue_) = names(result_[[1]])
        eachValue_$feature = colnames(expr_)[-c(1, 2)]
        eachOmic_[[omic_]] = eachValue_
      }

      return(eachOmic_)

    },
    transform = function(fun_, omics_ = self$omics, condition_ = "default", target_ = self$target, enableProgressBar_ = F, return_ = T) {

      eachOmic_ = list()
      for (omic_ in omics_) {

        cat("Calculating ", omic_, "...\n", sep = "")
        expr_ = private$.data(omic_, condition_, target_)

        targetFun_ = function (ind__, data__) {
          fun_(data__[[ind__]])
        }

        result_ = private$.apply(targetFun_, expr_, enableProgressBar_ = enableProgressBar_)
        result_ = as.data.table(c(list(expr_[[1]]), result_))
        colnames(result_) = colnames(expr_)[-2]
        eachOmic_[[omic_]] = result_
      }

      if (return_) {
        return(eachOmic_)
      } else {
        self[omics_] = eachOmic_[omics_]
      }

    },
    data = function(omic_, condition_ = "default", target_ = self$target) {

      private$.data(omic_, condition_, target_)

    }
  ),
  private = list(
    .checkIt = function () {

      trans_ = function (df__) {
        colnames(df__)[1] = "sample"
        df__ = as.data.table(df__)
        setkey(df__, "sample")

        df__
      }


      self$exprs = lapply(self$exprs, trans_)
      names(self$exprs) = self$omics

      samplesInExprs_ = unique(unlist(lapply(self$exprs, `[[`, "sample")))

      self$phen = trans_(self$phen)
      self$phen = self$phen[samplesInExprs_]
      setkey(self$phen, "sample")

    },
    .getTargetSamples = function (target_ = self$target) {

      if (!is.null(target_)) {
        toHave_ = c("sample", target_)
        phen_ = self$phen[, ..toHave_]
        phen_ = phen_[!is.na(phen_[[target_]])]
      } else {
        phen_ = self$phen
      }

      return(phen_$sample)

    },
    .getSamples = function (condition_, target_ = self$target) {

      intersect(self$rules$.getSamples(self$phen, condition_), private$.getTargetSamples(target_))

    },
    .extract = function(fea_, condition_ = "default", target_ = self$target) {

      samples_ = private$.getSamples(condition_, target_)

      result_ = vector("list", length(self$exprs))
      namesForResult_ = names(self$exprs)
      for (i in seq_along(self$exprs)) {

        if (fea_ %in% colnames(self$exprs[[i]])) {
          result_[[i]] = self$exprs[[i]][samples_, ..fea_]
        } else {
          result_[[i]] = NA
          cat(fea_, "doesn't exist in", namesForResult_[[i]], "return NA instead\n")
        }

      }
      names(result_) = paste("value", namesForResult_, sep = "_")
      result_ = as.data.table(result_)

      if (all(is.na(result_))) {
        warning(fea_, " doesn't exist in any omic data, return NULL")
        return(NULL)
      }

      result_$sample = samples_
      result_$gene   = fea_

      setkey(result_, "sample")

      return(result_)

    },
    .allNaRowsOf = function(result_) {

      toDrop_ = c("sample", "gene")
      naRows_ = Reduce(`&`, lapply(result_[, !..toDrop_], is.na))

      return(naRows_)

    },
    .data = function(omic_, condition_ = "default", target_ = self$target) {

      expr_    = self$exprs[[omic_]]
      samples_ = intersect(private$.getSamples(condition_, target_), expr_$sample)

      toHave_ = c("sample", target_)
      phen_ = self$phen[samples_, ..toHave_]
      expr_ = merge(phen_, expr_, all.x = T)

      return(expr_)

    },
    .apply = function(warpedFun_, data_, enableProgressBar_ = F) {

      inds_ = 3:ncol(data_)

      if (enableProgressBar_) {
        suppressWarnings(library(progressr))
        handlers(global = TRUE)

        p = progressor(along = inds_)
        finalFun_ = function(ind__, data__) {
          p(sprintf("x=%g", ind__))
          warpedFun_(ind__, data__)
        }
      } else {
        finalFun_ = function(ind__, data__) {
          warpedFun_(ind__, data__)
        }
      }

      result_ = suppressMessages(future.apply::future_lapply(inds_, finalFun_, data__ = data_))

      return(result_)

    }
  )
)

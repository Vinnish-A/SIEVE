geneFeature = R6Class(
  classname = "geneFeature",
  public = list(
    values = NULL,
    phen   = NULL,
    omics  = NULL,
    target = NULL,
    initialize = function(values_, phen_, omics_, target = NULL, fixPhen_ = T) {

      self$values = values_
      self$phen   = phen_
      self$omics  = omics_

      if (fixPhen_) {
        notBadSamples_ = unique(unlist(lapply(self$values, \(value__) value__$sample)))
        self$phen = self$phen[notBadSamples_]
      }

    },
    setTargets = function (target_) {

      self$target = target_

    },
    getValue = function(target_ = self$target, omics_ = self$omics, wayOfRmNa_ = "any") {

      results_ = lapply(self$values, private$.rmNA, way_ = wayOfRmNa_, paste0("value_", omics_))

      if (!is.null(target_)) {
        targets_ = lapply(results_, \(each__) self$phen[each__$sample, ..target_])
        results_ = lapply(1:length(targets_), \(i__) cbind(results_[[i__]], targets_[[i__]]))
      }

      names(results_) = names(self$values)
      return(results_)
    }
  ),
  private = list(
    .rmNA = function(dt_, cols_, way_ = "any") {

      samplesToDrop_ = switch(
        way_,
        all = Reduce(`&`, lapply(dt_[, ..cols_], is.na)),
        any = Reduce(`|`, lapply(dt_[, ..cols_], is.na))
      )

      return(dt_[!samplesToDrop_])

    }
  )
)

loadAll = function (phen_, condition_) {

  samples_ = phen_[, list(sample)][[1]]

  return(samples_)

}

loadTCGA = function (phen_, condition_) {

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

transposition_of = function(dt_) {

  dt_ = dt_ |> dplyr::distinct(sample, .keep_all = T)

  genes = dt_$sample
  samples = colnames(dt_)[-1]

  to_drop = "sample"
  dt_ = transpose(dt_[, !..to_drop])
  dt_ = dt_[, lapply(.SD, \(vec) ifelse(is.na(vec), mean(vec, na.rm = T), vec))]

  colnames(dt_) = genes
  dt_[, sample := samples]

  neworder_ = colnames(dt_)[c(ncol(dt_), 2:ncol(dt_)-1)]
  dt_[, ..neworder_]

}

integration_of = function(expr_, phen_, to_trans_ = T) {

  phen_ = as.data.table(phen_)

  if (to_trans_) {
    expr_ = transposition_of(expr_)
  } else {
    expr_ = as.data.table(expr_)
  }

  con_ = expr_[, c(1, 2)]; colnames(con_)[2] = "test"
  phen_ = merge(con_, phen_, by = "sample", all = T)

  setkey(expr_, "sample"); setkey(phen_, "sample")

  list(expr = expr_, phen = phen_)

}

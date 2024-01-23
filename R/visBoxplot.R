mrnaVsSnp = function(metaInte_, fea_) {

  data_ = metaInte_$extract(fea_)$getValue() |>
    rbindlist() |>
    mutate(value_snp = factor(ifelse(value_snp == 0, "wt", "mt"), levels = c("wt", "mt")))

  library(rstatix)
  library(ggsignif)

  statWilcox_ = data_ %>%
    group_by(gene) %>%
    wilcox_test(value_mRNA ~ value_snp) %>%
    add_xy_position(x = "gene")

  jitterCalculator_ = function(data__, class__, inClass__, col1__ = "dist_cat_n", col2__ = "scat_adj", dodge__ = dodge_) {

    data__[[col1__]] = as.numeric(factor(data__[[class__]]))
    data__[[col2__]] = dodge__*(as.numeric(data__[[inClass__]]) - sum(range(as.numeric(data__[[inClass__]])))/2)/length(levels(data__[[inClass__]]))

    data__

  }

  dodge_ = 1
  data_ |>
    mutate(gene = str_to_title(gene)) |>
    jitterCalculator_("gene", "value_snp") |>
    ggplot(aes(x = gene, y = value_mRNA)) +
    geom_boxplot(aes(color = value_snp), outlier.size = 0, position = position_dodge(width = dodge_), width = 0.5, outlier.stroke = 0) +
    geom_jitter(aes(x = scat_adj + dist_cat_n, color = value_snp), alpha = 0.2, width = 0.1, stroke = 0) +
    scale_color_manual(values = c("firebrick", "navy"), name = "Type") +
    ylim(NA, ceiling(max(statWilcox_$y.position)*1.1)) +
    xlab("") +
    ylab("Expression") +
    theme_bw() +
    theme(legend.position = "top") +
    geom_signif(annotations = statWilcox_$p,
                xmin = statWilcox_$xmin,
                xmax = statWilcox_$xmax,
                y_position = max(statWilcox_$y.position)*1.1,
                color="black")

}

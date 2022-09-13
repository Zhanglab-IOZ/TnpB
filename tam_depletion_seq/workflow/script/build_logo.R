suppressPackageStartupMessages({
    library(tidyverse)
    library(magrittr)
    library(broom)
    library("optparse")
})

args <- OptionParser() |>
    add_option("--in-test") |>
    add_option("--out-logodata") |>
    add_option("--n-boot", type = "integer", default = 1000) |>
    add_option("--seed", type = "integer", default = NULL) |>
    parse_args2()

test_df <- read_tsv(args$options$in_test, show_col_types = FALSE)
control_total <- sum(test_df$control)
test_total <- sum(test_df$test)

set.seed(args$options$seed)

if (!any(test_df$sig.filter)) {
    cat("idx\tA\tC\tG\tT\tentropy\tlower\tupper\n", file = args$options$out_logodata)
    quit()
}

test_df |>
    filter(sig.filter) |>
    # bootstrapped fold change
    mutate(
        boot_logfc = (
            ((rpois(n() * args$options$n_boot, test + 0.5) + 1) / test_total) /
                ((rpois(n() * args$options$n_boot, control + 0.5) + 1) / control_total)
        ) |> log() |> matrix(ncol = args$options$n_boot)
    ) |>
    # weight letters by log fold change
    select(tam, logfc, boot_logfc) |>
    mutate(chars = str_split(tam, "", 5)) |>
    unnest_longer(chars, values_to = "char", indices_to = "idx") |>
    mutate(char = factor(char, levels = c("A", "C", "G", "T"))) |>
    group_by(idx, char, .drop = FALSE) |>
    summarise(
        weight = sum(logfc),
        boot_weight = boot_logfc |> pmin(0) |> colSums() |> t(),
        .groups = "drop_last"
    ) |>
    # calculate frequency and entropy
    mutate(
        weight = weight / sum(weight),
        boot_weight = boot_weight / matrix(colSums(boot_weight), nrow = nrow(boot_weight), ncol = args$options$n_boot, byrow = TRUE)
    ) |>
    summarise(
        weights = setNames(weight, char) |> list(),
        entropy = log(4) + sum(weight * replace(log(weight), weight == 0, 0)),
        boot_entropy = log(4) + colSums(boot_weight * replace(log(boot_weight), boot_weight == 0, 0)) |> t()
    ) |>
    unnest_wider(weights) |>
    # get confident interval from bootstrapped entropy
    mutate(
        entropy_ci = apply(boot_entropy, 1, quantile, probs = c(0.05, 0.95)) |> matrix(nrow = 2) |> `rownames<-`(c("lower", "upper")) |> array_branch(margin = 2)
    ) |>
    select(!boot_entropy) |>
    unnest_wider(entropy_ci) |>
    mutate_all(format, digits = 0, nsmall = 5, trim = TRUE) |>
    write_tsv(file = args$options$out_logodata)
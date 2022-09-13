suppressPackageStartupMessages({
    library(tidyverse)
    library(magrittr)
    library(broom)
    library("optparse")
})

args <- OptionParser() |>
    add_option("--out-test") |>
    add_option("--out-plot") |>
    add_option("--z-score", type = "integer", default = 5) |>
    parse_args2()

test_df <-
    full_join(
        read_tsv(file = args$args[1], col_names = c("tam", "control"), col_types = "ci"),
        read_tsv(file = args$args[2], col_names = c("tam", "test"), col_types = "ci"),
        by = "tam"
    ) %>%
    replace_na(list(control = 0, test = 0)) %T>%
    with({
        control_total <<- sum(control)
        test_total <<- sum(test)
    }) %>%
    mutate(
        enrich_control = (control + 1) / (control_total + 2) * 1024,
        enrich_test = (test + 1) / (test_total + 2) * 1024,
        logfc = log(enrich_test / enrich_control),
        var_sampling = 1 / (test + 1) + 1 / (control + 1)
    ) %T>%
    with({
        var_total_full <<- mean(logfc^2)
    }) %T>% {
        . |>
            filter(logfc > 0) |>
            with({
                var_total <<- mean(logfc[logfc > 0]^2)
                var_amplify <<- mean(logfc[logfc > 0]^2 - var_sampling[logfc > 0])
            })
    } %>%
    mutate(
        sd = sqrt(var_sampling + var_amplify),
        z.score = logfc / sd,
        p.value = pnorm(z.score),
        q.value = p.adjust(p.value, "fdr"),
        sig.filter = z.score < -args$options$z_score
    ) %>%
    arrange(p.value)

write_tsv(test_df, file = args$options$out_test)

sd_df <- tibble(
    y = seq(max(-6, log(1024 / (test_total))), 5, length.out = 101) |> exp(),
    z = list(c(-args$options$z_score)),
) |>
    unnest_longer(z) |>
    rowwise() |>
    mutate(
        solve = optimize(function(x) {
            (z * sqrt(1 / (x * (control_total + 2) / 1024) + 1 / (y * (test_total + 2) / 1024) + var_amplify) - log(y / x))^2
        }, lower = 0, upper = 1e3, tol = 1e-9) |> list(),
    ) |>
    ungroup() |>
    unnest_wider(solve) |>
    filter(objective < 1e-9) |>
    select(x = minimum, y, z)

ggsave(
    args$options$out_plot,
    width = 5,
    height = 4,
    plot = ggplot() +
        scale_x_log10() +
        scale_y_log10() +
        coord_fixed(xlim = range(test_df$enrich_control), ylim = range(test_df$enrich_test)) +
        geom_point(aes(x = enrich_control, y = enrich_test, color = sig.filter), data = test_df, alpha = 0.8) +
        labs(
            x = "Control TAM abundance",
            y = "Experiment TAM abundance",
            color = sprintf("z < %s", -args$options$z_score)
        ) +
        # theme_classic() +
        geom_abline() +
        # geom_abline(slope = 1, intercept = -args$options$z_score / log(10) * sqrt(var_total_full), color = "red", linetype = "dashed") +
        # geom_abline(slope = 1, intercept = -args$options$z_score / log(10) * sqrt(var_total), color = "red", linetype = "solid") +
        geom_path(aes(x = x, y, group = z), data = sd_df, color = "blue", linetype = "solid")
)

library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(igraph)
library(progressr)
library(ggplot2)


meta <- readRDS("data/cdi-metadata-pos_vid.rds") |> as_tibble()
g <- readRDS("data/child_net_graph.rds") |> upgrade_graph()
pure_growth <- readRDS("./pure_growth_v4.rds")

multisample <- function(n, x, fct, simplify2vec = FALSE) {
    k <- split(x, fct)
    y <- mapply(sample, split(x, fct), n, SIMPLIFY = FALSE)
    return(if(simplify2vec) unname(do.call(c, y)) else y)
}


future::plan("multicore")

# POS matched sampling ----
niter <- 1000
tmp <- with_progress(seq_len(niter) |>
    future_map(function(iter, data, g, meta, p) {
        map(data, function(x_seed, g, meta, p) {
            p()
            map(x_seed, function(x_model, g, meta) {
                selected <- logical(vcount(g))
                netstats <- tibble(
                    nproduced = integer(length(x_model)),
                    indegree = numeric(length(x_model)),
                    clustcoef = numeric(length(x_model)),
                    aspl = numeric(length(x_model))
                )
                for (i in seq_along(x_model)) {
                    z <- !selected
                    pos_count <- table(x_model[[i]]$pos)
                    vid_sample <- multisample(pos_count, meta$vid[z], meta$pos[z]) |>
                        unlist()
                    selected[vid_sample] <- TRUE
                    vid_vocab <- which(selected)
                    g_vocab <- induced_subgraph(g, vid_vocab)
                    netstats$nproduced[i] <- sum(selected)
                    netstats$indegree[i] <- median(degree(g_vocab, mode = "in"))
                    netstats$clustcoef[i] <- transitivity(g_vocab, type = "global")
                    netstats$aspl[i] <- mean_distance(g_vocab)
                }
                return(netstats)
            },
            g = g,
            meta = meta)
        },
        g = g,
        meta = meta,
        p = p
        )
    },
    data = pure_growth$pure_growth,
    g = g,
    meta = meta,
    p = progressor(niter * 200),
    .options = furrr_options(seed = TRUE)
    )
)

saveRDS(tmp, "pure-growth-ran-netstats_v4.rds")


# Uniform random sampling ----
niter <- 10000
ra_unif <- with_progress(seq_len(niter) |>
    future_map(function(iter, g, p) {
        p()
        selected <- logical(vcount(g))
        vocab_sizes <- seq(60, 600, by = 20)
        sample_sizes <- diff(c(0, vocab_sizes))
        netstats <- tibble(
            nproduced = integer(length(vocab_sizes)),
            indegree = numeric(length(vocab_sizes)),
            clustcoef = numeric(length(vocab_sizes)),
            aspl = numeric(length(vocab_sizes))
        )
        for (i in seq_along(vocab_sizes)) {
            z <- !selected
            vid_sample <- sample(meta$vid[z], size = sample_sizes[i])
            selected[vid_sample] <- TRUE
            vid_vocab <- which(selected)
            g_vocab <- induced_subgraph(g, vid_vocab)
            netstats$nproduced[i] <- sum(selected)
            netstats$indegree[i] <- median(degree(g_vocab, mode = "in"))
            netstats$clustcoef[i] <- transitivity(g_vocab, type = "global")
            netstats$aspl[i] <- mean_distance(g_vocab)
        }
        return(netstats)
    },
    g = g,
    p = progressor(niter),
    .options = furrr_options(seed = TRUE)
    )
)

saveRDS(ra_unif, "pure-growth-unif-ran-netstats_v4.rds")
#
pure_growth$pure_growth[[1]]$loa[[1]]

#pure_growth$netstats <- pure_growth$pure_growth |>
pure_growth$netstats <- pure_growth$pure_growth |>
    future_map(function(x_seed, g, p) {
        map(x_seed, function(x_model, g, p) {
            p()
            selected <- logical(vcount(g))
            netstats <- tibble(
                nproduced = integer(length(x_model)),
                indegree = numeric(length(x_model)),
                clustcoef = numeric(length(x_model)),
                aspl = numeric(length(x_model))
            )
            for (i in seq_along(x_model)) {
                z <- !selected
                pos_count <- table(x_model[[i]]$pos)
                vid_sample <- x_model[[i]]$ind
                selected[vid_sample] <- TRUE
                g_vocab <- induced_subgraph(g, which(selected))
                netstats$nproduced[i] <- sum(selected)
                netstats$indegree[i] <- median(degree(g_vocab, mode = "in"))
                netstats$clustcoef[i] <- transitivity(g_vocab, type = "global")
                netstats$aspl[i] <- mean_distance(g_vocab)
            }
            return(netstats)
        }, g = g, p = p)
    }, g = g, p = progressor(600)) |> with_progress()


pgnetstats[[1]]$loa
tmp[[1]][[1]]$loa
pure_growth$netstats[[1]]$loa

saveRDS(pure_growth |> select(seed_ind, group, seed_type, netstats), "./pure_growth_netstats_v4.rds")


# Summarize POS MATCHED RANs ----
pure_growth_ran_netstats_summary <- tmp |>
    map(function(x_iter) {
        map(x_iter, function(x_seed) {
            list_rbind(x_seed, names_to = "model") |> mutate(model = as.factor(model))
        }) |> list_rbind(names_to = "seed_ind")
    }) |> list_rbind(names_to = "iter") |>
    group_by(seed_ind, model, nproduced) |>
    summarize(across(c(indegree, clustcoef, aspl), list(mean_ran = mean, sd_ran = sd))) |>
    ungroup()



# Summarize UNIFORM RANs ----
pure_growth_unif_ran_netstats_summary <- ra_unif |>
    list_rbind(names_to = "iter") |>
    group_by(nproduced) |>
    summarize(across(c(indegree, clustcoef, aspl), list(mean_ran = mean, sd_ran = sd))) |>
    ungroup()


pure_growth_netstats <- pure_growth$netstats |>
    map(function(x_seed) {
        list_rbind(x_seed, names_to = "model") |> mutate(model = as.factor(model))
    }) |> list_rbind(names_to = "seed_ind") |>
    left_join(pure_growth_unif_ran_netstats_summary, by = c("nproduced")) |>
    mutate(
        indegree_z = (indegree - indegree_mean_ran) / indegree_sd_ran,
        clustcoef_z = (clustcoef - clustcoef_mean_ran) / clustcoef_sd_ran,
        aspl_z = (aspl - aspl_mean_ran) / aspl_sd_ran
    ) |> left_join(
        pure_growth |> ungroup() |> select(group) |> mutate(seed_ind = 1:200)
    )


saveRDS(pure_growth_netstats, "./pure_growth_netstats_with_ran_z_v4.rds")

# Beyond this point is a mess... but it is where all the figures are being made
tmpplot <- pure_growth_netstats |>
    drop_na() |>
    filter(nproduced <= 600, group == "NA") |>
    group_by(model, nproduced) |>
    summarize(across(ends_with("_z"), list(m=mean,s=sd), .names = "{.col}:{.fn}")) |>
    ungroup() |>
    pivot_longer(
        cols = contains("_z"),
        names_to = c("metric", ".value"),
        names_sep = ":"
    ) |>
    ggplot(aes(x = nproduced, y = m, color = model)) +
        geom_ribbon(aes(ymin=m-s, ymax=m+s, fill = model), alpha = .2) +
        geom_line() +
        facet_wrap(~metric, scale="free_y") +
        theme_bw(base_size = 18)

ggsave("tmpplot_v4.pdf", tmpplot, width = 11, height = 4, unit = "in")


xsim <- pure_growth_netstats |>
    drop_na() |>
    filter(nproduced <= 600, group == "NA") |>
    group_by(model, nproduced) |>
    summarize(across(ends_with("_z"), list(m=mean,s=sd), .names = "{.col}:{.fn}")) |>
    ungroup() |>
    pivot_longer(
        cols = contains("_z"),
        names_to = c("metric", ".value"),
        names_sep = ":"
    )

cdi <- readRDS("data/asd-na_osg-2023_11_20-forms.rds") |>
    inner_join(meta, by = "num_item_id") |>
    filter(group == "NA")


yy <- cdi |>
    select(subjectkey, interview_age, sex, nproduced, form) |>
    distinct() |>
    left_join(
        readRDS("../asd-na_netstats-bs/data/asd-na_netstats-ran.rds") |>
            tibble() |>
            filter(nproduced >= 60, nproduced <= 600, group == "NA") |>
            select(subjectkey, interview_age, sex, nproduced,indegree_z = z_indegree_med, clustcoef_z = z_clust, aspl_z = z_dist)
        )

yy |> filter(nproduced >= 60, nproduced <= 600) |> nrow()


xkids <- readRDS("../asd-na_netstats-bs/data/asd-na_netstats-ran.rds") |>
    tibble() |>
    filter(nproduced >= 60, nproduced <= 600, group == "NA") |>
    rename(indegree_z = z_indegree_med, clustcoef_z = z_clust, aspl_z = z_dist) |>
    pivot_longer(
        cols = ends_with("_z"),
        names_to = "metric",
        values_to = "m"
    ) |>
    select(subjectkey, nproduced, metric, m) |>
    mutate(s = NA, model = "kid")


saveRDS(xkids, "kids-netstats-ran_v4.rds")

tmpplot_kids <- ggplot(xsim, aes(x = nproduced, y = m, color = model)) +
    geom_point(data = xkids, color = "grey") +
    geom_ribbon(aes(ymin=m-s, ymax=m+s, fill = model), alpha = .2) +
    geom_line() +
    xlim(c(0, 600)) +
    facet_wrap(~metric, scale="free_y") +
    theme_bw(base_size = 18)

ggsave("tmpplot_kids_v4.pdf", tmpplot_kids, width = 11, height = 4, unit = "in")

tmpplot_kids_smooth <- ggplot(xsim, aes(x = nproduced, y = m, color = model)) +
    geom_point(data = xkids, color = "grey") +
    geom_smooth(data = xkids, color = "black", fill = "black", alpha = .3) +
    geom_ribbon(aes(ymin=m-s, ymax=m+s, fill = model), alpha = .2) +
    geom_line() +
    xlim(c(0, 600)) +
    facet_wrap(~metric, scale="free_y") +
    theme_bw(base_size = 18)

ggsave("puregrowth_unif-ran_kids_smooth_v4.pdf", tmpplot_kids_smooth, width = 11, height = 4, unit = "in")


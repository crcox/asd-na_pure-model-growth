library(dplyr)
library(tidyr)
library(netgrowr)
library(purrr)
library(furrr)
library(progressr)
library(igraph)
library(ggplot2)

source("./generate-seed-vocabularies.R")

# Load network of child-oriented associations ----
g <- upgrade_graph(readRDS("data/child_net_graph.rds"))
adjmat <- as_adjacency_matrix(g)


# Load WordBank and NDAR CDI expressive vocabularies ----
#   group     n
# 1 NA     1416
# 2 ASD     472
meta <- readRDS("data/cdi-metadata-pos_vid.rds")

cdi <- readRDS("data/asd-na_osg-2023_11_20-forms.rds") |>
    inner_join(meta, by = "num_item_id") |>
    group_by(group, word)

word_fits <- cdi |>
    group_split() |>
    map(\(x) {
        glm(produced ~ nproduced, data = x, family = "binomial")
    }, .progress = TRUE) |>
    tibble(word_fit = _) |>
    bind_cols(group_keys(cdi)) |>
    left_join(meta, by = "word")


seed_vocabs <- generate_seed_vocabularies(word_fits, size = 60) |>
    list_rbind(names_to = "seed_ind") |>
    group_by(group, seed_ind) |>
    nest()


within_group_diffs <- split(seed_vocabs, seed_vocabs$group) |>
    map(\(x) {
        apply(combn(100, 2), 2, function(ij) {
            a <- x$data[[ij[1]]]$num_item_id
            b <- x$data[[ij[2]]]$num_item_id
            sum(a %in% b)
        })
    })

within_group_diffs |> map(summary)


between_group_diffs <- expand.grid(i = 1:100, j = 1:100) |>
    pmap(
        \(i, j, a, b) {
            sum(a[[i]]$num_item_id %in% b[[j]]$num_item_id)
        },
        a = seed_vocabs[seed_vocabs$group == "ASD", "data"]$data,
        b = seed_vocabs[seed_vocabs$group == "NA", "data"]$data
    ) |>
    unlist()

summary(between_group_diffs)


future::plan("multicore")

pure_growth <- map(seed_vocabs$data, function(x, adjmat, p) {
    p()
    map(
        list(
            loa = lure_of_the_associates,
            pac = preferential_acquisition,
            pat = preferential_attachment
        ),
        function(growth_model, known, adjmat) {
            learned <- list(tibble(
                ind = unname(which(known)),
                word = rownames(adjmat)[known]
            ))
            vocab_sizes <- seq(80, 600, by = 20)
            for (i in seq_along(vocab_sizes)) {
                gv <- growth_model(adjmat, known)
                r <- rank(gv, ties.method = "random", na.last = FALSE)
                r <- (max(r) - r)
                k <- r < 20
                learned <- append(learned, list(tibble(ind = unname(which(k)), word = rownames(adjmat)[k])))
                known[k] <- TRUE
            }
            return(learned)
        },
        known = rownames(adjmat) %in% x$word, adjmat = adjmat
    )
}, adjmat = as.matrix(adjmat), p = progressor(200)) |> with_progress()
#, .options = furrr_options(seed = TRUE)

seed_vocabs$pure_growth <- pure_growth
seed_vocabs$seed_type <- ifelse(seed_vocabs$seed_ind > 1, "sampled", "maxprob") |> as.factor()
seed_vocabs <- rename(seed_vocabs, seed_vocab = data)

seed_vocabs$pure_growth <- seed_vocabs$pure_growth |>
    map_depth(.depth=3, ~ left_join(.x, select(meta, ind=vid, pos), by="ind"), .progress=TRUE)

saveRDS(seed_vocabs, "pure_growth_v3.rds")
saveRDS(within_group_diffs, "seed_vocabs_within_group_diffs_v3.rds")
saveRDS(between_group_diffs, "seed_vocabs_between_group_diffs_v3.rds")


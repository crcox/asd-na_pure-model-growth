library(dplyr)
library(tidyr)
library(netgrowr)
library(purrr)
library(readr)
library(igraph)
library(ggplot2)

source("./R/generate-seed-vocabularies.R")

# Load network of child-oriented associations ----
g <- upgrade_graph(readRDS("data/child_net_graph.rds"))
adjmat <- as_adjacency_matrix(g)


# Load WordBank and NDAR CDI expressive vocabularies ----
#   group     n
# 1 NA     1416
# 2 ASD     472
meta <- readRDS("data/cdi-metadata-pos_vid.rds")

wb_to_ndar <- read_csv(
    "data/wb-to-ndar.csv",
    col_types = c(
        word = col_character(),
        ndar_item_id = col_integer(),
        wb_item_id = col_integer()
    )
) |> mutate(across(ends_with("_id"), as.integer))

cdi_bad_asd_id <- readRDS("data/asd-na_osg-2023_11_20-forms.rds") |>
    mutate(num_item_id = as.integer(num_item_id)) |>
    group_by(group)

cdi_split <- group_split(cdi_bad_asd_id)
names(cdi_split) <- group_keys(cdi_bad_asd_id) |> pull(group)


cdi_asd_fixed_id <- cdi_split$ASD |>
    left_join(wb_to_ndar |> select(num_item_id = ndar_item_id, wb_item_id), by = "num_item_id") |>
    mutate(num_item_id = wb_item_id) |>
    select(-wb_item_id)

cdi_split$ASD <- cdi_asd_fixed_id

cdi <- list_rbind(cdi_split) |>
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


pure_growth <- map(seed_vocabs$data, function(x, adjmat) {
    k <- rownames(adjmat) %in% x$word
    map(
        c(loa = lure_of_the_associates,
          pac = preferential_acquisition,
          pat = preferential_attachment),
        function(growth_model, known, adjmat) {
            learned <- list(tibble(
                ind = unname(which(known)),
                word = rownames(adjmat)[known]
            ))
            while (any(!known)) {
                gv <- growth_model(adjmat, known)
                r <- rank(gv, ties.method = "random")
                k <- r <= 20
                learned <- append(learned, list(tibble(ind = unname(which(k)), word = rownames(adjmat)[k])))
                known <- known | k
            }
            return(learned)
        },
        known = k, adjmat = adjmat
    )
}, adjmat = as.matrix(adjmat), .progress = TRUE)

seed_vocabs$pure_growth <- pure_growth
seed_vocabs$seed_type <- ifelse(seed_vocabs$seed_ind > 1, "sampled", "maxprob") |> as.factor()
seed_vocabs <- rename(seed_vocabs, seed_vocab = data)

seed_vocabs$pure_growth <- seed_vocabs$pure_growth |>
    map_depth(.depth=3, ~ left_join(.x, select(meta, ind=vid, pos), by="ind"), .progress=TRUE)

saveRDS(seed_vocabs, "pure_growth_20250722.rds")
saveRDS(within_group_diffs, "seed_vocabs_within_group_diffs_20250722.rds")
saveRDS(between_group_diffs, "seed_vocabs_between_group_diffs_20250722.rds")

seed_vocabs$pure_growth[[1]]$loa |> map_int(nrow)
seed_vocabs$pure_growth[[1]]$pac |> map_int(nrow)
seed_vocabs$pure_growth[[1]]$pat |> map_int(nrow)

diff_counts <- matrix(0, nrow=32, ncol=3, dimnames=list(NULL, c("loa", "pac", "pat")))
for (i in 1:32) {
    diff_counts[i, ] <- c(
        loa=setdiff(pure_model_growth$autistic$loa[[i]]$word, pure_model_growth$nonautistic$loa[[i]]$word) |> length(),
        pac=setdiff(pure_model_growth$autistic$pac[[i]]$word, pure_model_growth$nonautistic$pac[[i]]$word) |> length(),
        pat=setdiff(pure_model_growth$autistic$pat[[i]]$word, pure_model_growth$nonautistic$pat[[i]]$word) |> length()
    )
}


pure_model_growth_cumulative <- map_depth(pure_model_growth, 2, \(x) {
    y <- list(x[[1]])
    for (i in 2:32) {
        y <- append(y, list(bind_rows(y[[i-1]], x[[i]])))
    }
    return(y)
})



netstats <- map_depth(pure_model_growth_cumulative, 3, function(x, g) {
    g <- igraph::induced_subgraph(g, x$ind)
    cc <- igraph::transitivity(g, type = "localaverage")
    aspl <- igraph::mean_distance(g)
    indegree <- median(igraph::degree(g, mode = "in"))
    return(tibble(vocab_size = nrow(x), clust_coef = cc, aspl = aspl, indegree = indegree))
}, g = g, .progress = TRUE) |>
    map_depth(2, bind_rows) |>
    map(\(x) {
        list_rbind(x, names_to = "growth_model")
    }) |>
    list_rbind(names_to = "group") |>
    pivot_longer(
        cols = c(clust_coef, aspl, indegree),
        names_to = "metric",
        values_to = "value"
    ) |>
    mutate(across(c(group, growth_model, metric), as.factor))

p <- ggplot(netstats, aes(x = vocab_size, y = value, color = group)) +
    geom_line() +
    facet_grid(metric ~ growth_model, scales = "free_y")

ggsave("pure_model_growth_netstats_02.png", p)

# Load VSOA and metadata ----
vsoa <- readRDS("data/asd-na_vsoa.rds") |>
    select(num_item_id, vsoa_autistic = "ASD,ASD", vsoa_nonautistic = "NA.NA", vsoa_diff = "ASD-NA.ASD") |>
    as_tibble()

cdi_metadata_preproc <- left_join(readRDS('data/cdi-metadata-preproc.rds'), vsoa)

vsoa_df <- cdi_metadata_preproc |>
    select(num_item_id, cue_CoxHae, vsoa_autistic, vsoa_nonautistic) |>
    pivot_longer(
        cols = c(vsoa_autistic, vsoa_nonautistic),
        names_to = "group",
        values_to = "vsoa",
        names_prefix = "vsoa_"
    ) |>
    mutate(num_item_id = as.integer(num_item_id), group = as.factor(group)) |>
    group_by(group) |>
    mutate(
        vsoa_sample_prob = if_else(vsoa < -100, -100, if_else(vsoa > 200, NA, vsoa)),
        vsoa_sample_prob = (max(vsoa_sample_prob, na.rm = TRUE) - vsoa_sample_prob) / sum(vsoa_sample_prob, na.rm = TRUE)
    ) |>
    ungroup()

# Generate seed vocabularies ----
# Selecting the 60 words with lowest VSOA, without ties, will always return the
# same seed vocabulary for each group. However, it has drawbacks. For
# example, if the 60th and 61st VSOA are 80, both words are equally likely
# to be in a child's early vocabulary but one of them will be included in
# the seed vocabulary and the other will be excluded. To address this, 99
# additional seed vocabularies are selected by sampling words according to
# a probability distribution determined by their VSOA. The procedure is:
#
#  1. Set all negative VSOAs to zero. This is done because model estimates of VSOA become outrageously extreme the most and least likely words to be known, and these extreme values will skew the probability distribution. To
DETERMINISTIC_SEED_VOCAB <- TRUE
seed_vocabs <- vector("list", 100)
for (i in seq_along(seed_vocabs)) {
    if (DETERMINISTIC_SEED_VOCAB) {
        seed_vocabs[[i]] <- vsoa_df |>
            group_by(group) |>
            slice_min(order_by = vsoa, n = 60, with_ties = FALSE)
        DETERMINISTIC_SEED_VOCAB <- FALSE

    } else {
        seed_vocabs[[i]] <- vsoa_df |>
            filter(!is.na(vsoa_sample_prob)) |>
            group_by(group) |>
            slice_sample(n = 60, weight_by = vsoa_sample_prob)
    }
}


seed_index <- map(seed_vocabs, \(x) {
    split(x$cue_CoxHae, x$group)
})

tmp <- vapply(seed_index, \(x) {
    c(
        sum(x$autistic %in% seed_index[[1]]$autistic),
        sum(x$nonautistic %in% seed_index[[1]]$nonautistic)
    )
}, integer(2))

seed_index <- split(seed_vocabs[[1]]$cue_CoxHae, seed_vocabs[[1]]$group)

seed_comp <- list(
    aut_only = seed_index$autistic[!(seed_index$autistic %in% seed_index$nonautistic)],
    non_only = seed_index$nonautistic[!(seed_index$nonautistic %in% seed_index$autistic)],
    shared = seed_index$autistic[seed_index$autistic %in% seed_index$nonautistic],
    diff = setdiff(seed_index$autistic, seed_index$nonautistic)
)
print(seed_comp)

assertthat::are_equal(rownames(assocnet$child_preproc), colnames(assocnet$child_preproc))

seed_known <- map(seed_index, \(x, adjmat) {
    rownames(adjmat) %in% x
}, adjmat = adjmat)


pure_model_growth <- map(seed_known, \(k, net) {
    map(c(loa = lure_of_the_associates, pac = preferential_acquisition, pat = preferential_attachment), \(f, known, net) {
        learned <- list(tibble(ind = unname(which(known)), word = rownames(net)[known]))
        while (any(!known)) {
            gv <- f(net, known)
            r <- rank(gv, ties.method = "random")
            k <- r <= 20
            learned <- append(learned, list(tibble(ind = unname(which(k)), word = rownames(net)[k])))
            known <- known | k
        }
        return(learned)
    }, known = k, net = net)
}, net = as.matrix(adjmat))

lengths(pure_model_growth$nonautistic)

diff_counts <- matrix(0, nrow=32, ncol=3, dimnames=list(NULL, c("loa", "pac", "pat")))
for (i in 1:32) {
    diff_counts[i, ] <- c(
        loa=setdiff(pure_model_growth$autistic$loa[[i]]$word, pure_model_growth$nonautistic$loa[[i]]$word) |> length(),
        pac=setdiff(pure_model_growth$autistic$pac[[i]]$word, pure_model_growth$nonautistic$pac[[i]]$word) |> length(),
        pat=setdiff(pure_model_growth$autistic$pat[[i]]$word, pure_model_growth$nonautistic$pat[[i]]$word) |> length()
    )
}


pure_model_growth_cumulative <- map_depth(pure_model_growth, 2, \(x) {
    y <- list(x[[1]])
    for (i in 2:32) {
        y <- append(y, list(bind_rows(y[[i-1]], x[[i]])))
    }
    return(y)
})



netstats <- map_depth(pure_model_growth_cumulative, 3, function(x, g) {
    g <- igraph::induced_subgraph(g, x$ind)
    cc <- igraph::transitivity(g, type = "localaverage")
    aspl <- igraph::mean_distance(g)
    indegree <- median(igraph::degree(g, mode = "in"))
    return(tibble(vocab_size = nrow(x), clust_coef = cc, aspl = aspl, indegree = indegree))
}, g = g, .progress = TRUE) |>
    map_depth(2, bind_rows) |>
    map(\(x) {
        list_rbind(x, names_to = "growth_model")
    }) |>
    list_rbind(names_to = "group") |>
    pivot_longer(
        cols = c(clust_coef, aspl, indegree),
        names_to = "metric",
        values_to = "value"
    ) |>
    mutate(across(c(group, growth_model, metric), as.factor))

p <- ggplot(netstats, aes(x = vocab_size, y = value, color = group)) +
    geom_line() +
    facet_grid(metric ~ growth_model, scales = "free_y")

ggsave("pure_model_growth_netstats_02.png", p)

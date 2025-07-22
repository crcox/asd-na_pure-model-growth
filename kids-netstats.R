library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(progressr)
library(igraph)
library(ggplot2)

# Load network of child-oriented associations ----
g <- upgrade_graph(readRDS("data/child_net_graph.rds"))

cdi <- readRDS("data/asd-na_osg-2023_11_20-forms.rds") |>
    inner_join(meta, by = "num_item_id") |>
    group_by(group, form, sex, subjectkey, interview_age, nproduced)


kids_netstats <- group_split(cdi) |>
    map(\(x, g) {
        g_vocab <- induced_subgraph(g, x$vid[x$produced])
        tibble(
            indegree = median(degree(g_vocab, mode = "in")),
            clustcoef = transitivity(g_vocab, type = "global"),
            aspl = mean_distance(g_vocab)
        )
    }, g = g, .progress = TRUE) |>
    list_rbind() |>
    bind_cols(group_keys(cdi))

saveRDS(kids_netstats, "kids-netstats-raw-v3.rds")

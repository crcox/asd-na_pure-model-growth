#' Generate seed vocabularies from GLM fits
#'
#' @param glms A numeric vector of VSOA (vocabulary size of acquisition) estimates for each word for each group.
#' @param groups A factor coding each VSOA as belonging to each group.
#' @param size An integer indicating the desired size for seed vocabularies.
#' @param trim Lower and upper bounds for VSOA. Values below the lower bound
#'   will be set to `trim[1]`, and values above the upper bound will be set to
#'   `NA`, ensuring those words will never appear in seed vocabularies.
#'
#' @details
#' Selecting the 60 words with lowest VSOA, without ties, will always return the
#' same seed vocabulary for each group. However, it has drawbacks. For
#' example, if the 60th and 61st VSOA are 80, both words are equally likely
#' to be in a child's early vocabulary but one of them will be included in
#' the seed vocabulary and the other will be excluded. To address this, 99
#' additional seed vocabularies are selected by sampling words according to
#' a probability distribution determined by their VSOA. The procedure is:
#'
#'  1. Set all negative VSOAs to the lower trim bound. This is done because
#'     model estimates of VSOA become extreme for the most and least
#'     likely words to be known, and these extreme values will skew the
#'     probability distribution. To 
generate_seed_vocabularies <- function(word_fits, size = 60, gate = 0) {
    word_fits <- word_fits |>
        mutate(
            p_in_seed = map_dbl(word_fits$word_fit, function(model) {
                predict(model, data.frame(nproduced = size), type = "response")
            })
        ) |>
        drop_na() |>
        filter(p_in_seed > gate) |>
        mutate(p_sample = p_in_seed / sum(p_in_seed))

    DETERMINISTIC_SEED_VOCAB <- TRUE
    seed_vocabs <- vector("list", 100)
    for (i in seq_along(seed_vocabs)) {
        if (DETERMINISTIC_SEED_VOCAB) {
            seed_vocabs[[i]] <- word_fits |>
                group_by(group) |>
                slice_min(order_by = p_in_seed, n = size, with_ties = FALSE)
            DETERMINISTIC_SEED_VOCAB <- FALSE

        } else {
            seed_vocabs[[i]] <- word_fits |>
                group_by(group) |>
                slice_sample(n = size, weight_by = p_in_seed)
        }
    }

    return(seed_vocabs)
}



#' Generate seed vocabularies from VSOA estimates
#'
#' @param words A character vector (or factor).
#' @param vsoa A numeric vector of VSOA (vocabulary size of acquisition) estimates for each word for each group.
#' @param groups A factor coding each VSOA as belonging to each group.
#' @param size An integer indicating the desired size for seed vocabularies.
#' @param trim Lower and upper bounds for VSOA. Values below the lower bound
#'   will be set to `trim[1]`, and values above the upper bound will be set to
#'   `NA`, ensuring those words will never appear in seed vocabularies.
#'
#' @details
#' Selecting the 60 words with lowest VSOA, without ties, will always return the
#' same seed vocabulary for each group. However, it has drawbacks. For
#' example, if the 60th and 61st VSOA are 80, both words are equally likely
#' to be in a child\'s early vocabulary but one of them will be included in
#' the seed vocabulary and the other will be excluded. To address this, 99
#' additional seed vocabularies are selected by sampling words according to
#' a probability distribution determined by their VSOA. The procedure is:
#'
#'  1. Set all negative VSOAs to the lower trim bound. This is done because
#'     model estimates of VSOA become extreme for the most and least
#'     likely words to be known, and these extreme values will skew the
#'     probability distribution. To 
generate_seed_vocabularies_VSOA <- function(words, vsoa, group, size = 60, trim = NULL) {
    if (is.null(trim)) {
        vsoa_trim <- vsoa
    } else {
        vsoa_trim <- if_else(vsoa < trim[1], trim[1], if_else(vsoa > trim[2], NA, vsoa))
    }
    vsoa_df <- tibble(group=group, vsoa=vsoa, vsoa_trim=vsoa_trim) |>
        group_by(group) |>
        mutate(
            vsoa_sample_prob = if_else(vsoa < -100, -100, if_else(vsoa > 200, NA, vsoa)),
            vsoa_sample_prob = (max(vsoa_sample_prob, na.rm = TRUE) - vsoa_sample_prob) / sum(vsoa_sample_prob, na.rm = TRUE)
        ) |>
        ungroup()

    DETERMINISTIC_SEED_VOCAB <- TRUE
    seed_vocabs <- vector("list", 100)
    for (i in seq_along(seed_vocabs)) {
        if (DETERMINISTIC_SEED_VOCAB) {
            seed_vocabs[[i]] <- vsoa_df |>
                group_by(group) |>
                slice_min(order_by = vsoa, n = size, with_ties = FALSE)
            DETERMINISTIC_SEED_VOCAB <- FALSE

        } else {
            seed_vocabs[[i]] <- vsoa_df |>
                filter(!is.na(vsoa_sample_prob)) |>
                group_by(group) |>
                slice_sample(n = size, weight_by = vsoa_sample_prob)
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
}

library(dplyr)

x <- readRDS("data/asd-na_osg-2023_11_20.rds")

y <- x |>
    group_by(group, subjectkey, interview_age, nproduced, sex) |>
    mutate(form = factor(if_else(n()>400, 1, 2), levels = 1:2, labels = c("WS", "WG"))) |>
    ungroup()

y |>
    select(group, subjectkey, form, interview_age, sex) |>
    distinct() |>
    count(group)

saveRDS(y, "data/asd-na_osg-2023_11_20-forms.rds")

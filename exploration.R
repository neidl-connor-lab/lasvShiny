lf_covid <- LF_Covid_Compare %>%
    tibble() %>% filter(Disease %in% c("COVID", "LF")) %>%
    gather(key = "Factor", value = "Value", -Disease, -Subset)

lf_covid %>% ggplot(aes(Subset, Value)) +
    geom_bar(stat = "identity") +
    facet_wrap("Factor", scales = "free_y") +
    scale_y_log10() +
    geom_errorbar(stat = "identity")


adm_all_mat <- cyto_all_patients_dpi_corrected %>%
    select(id, Outcome, ct:IFNy) %>% group_by(id) %>% slice(1) %>%
    column_to_rownames(var = "id") %>%
    na.omit() %>%
    mutate(Outcome = factor(Outcome))

res.pca <- adm_all_mat %>% select(-Outcome) %>%
    as.matrix() %>%
    prcomp(scale = TRUE)

fviz_eig(res.pca)

fviz_pca_var(res.pca,
             col.var = adm_all_mat$Outcome,
             # palette = c("#00AFBB",  "#FC4E07"),
             # col.var = "contrib", # Color by contributions to the PC
             # gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)



adm.umap <- umap(adm_all_mat[,-1])
dims_umap <- adm.umap$layout %>%
    set_colnames(c("UMAP1", "UMAP2")) %>%
    data.frame() %>%
    cbind(Outcome = adm_all_mat$Outcome)

ggplot(dims_umap, aes(UMAP1, UMAP2, color=Outcome)) + geom_point()

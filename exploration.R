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



box.data %>% filter(Factor %in% c("PAI1")) %>%
ggplot(aes(x=Outcome, y = Value)) +
    geom_boxplot() +
    geom_point(aes(color=AgeBin), size=5, position = position_jitterdodge()) +
    scale_x_discrete(labels=c("F", "S", "OFI", "HC")) +
    scale_y_log10() +
    theme(text = element_text(size=24),
          plot.title = element_text(hjust = .5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.spacing = unit(1.25, "lines"),
          axis.title.x = element_blank(),
          axis.text = element_text(size = 20),
          axis.title.y = element_blank()
          )


pre.box.data %>% select(Outcome, ct:RANTES) %>%
    group_by(Outcome) %>% summarise_all(mean, na.rm=T)  %>%
    column_to_rownames(var = "Outcome") %>%
    t() %>% data.frame() %>% rownames_to_column("Factor") %>%
    dplyr::rename(f=Fatal,
                  s=Survivor,
                  ofi = "Other.Febrile.Illness",
                  hc = "Healthy.Control") -> mean.by.outcome


mean.by.outcome  %>% mutate_if(is.numeric, list(log2)) %>%
    mutate(f.s.lfc = f - s,
           f.hc.lfc = f - hc,
           s.hc.lfc = s - hc,
           ofi.hc.lfc = ofi - hc) %>%
    select(Factor, f.s.lfc, f.hc.lfc, s.hc.lfc, ofi.hc.lfc) %>%
    gather(key = "fc", value = "value", -Factor) %>%
    mutate(fc=factor(fc, levels = c("f.s.lfc", "f.hc.lfc", "s.hc.lfc", "ofi.hc.lfc")))-> fc.data

# %>%
#     filter(fc %in% c("f.fc", "s.fc", "ofi.fc"))


ggplot(fc.data, aes(x=Factor, y=value, fill=fc)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap("Factor", scales = "free")



# adm %>% select(ct:RANTES) %>% mutate_if(is.numeric, list(log2)) %>% cor(use = "complete.obs") %>%
#     # mutate_if(is.numeric, list( function (x) ( x - mean(x, na.rm = TRUE) ) / sd(x, na.rm = TRUE) )) %>%
#     # column_to_rownames(Outcome)
#     prcomp() %>% fviz_eig()


library(mclust)

table(adm$Outcome)



# %>%
    # Mclust(outcome)
# summary(mclust.obj)
# plot(mclust.obj, what = "classification")

glm(Outcome ~ PAI1 + IL8 + TM, family = binomial, data = adm) %>% visreg("PAI1", scale="response")

library(ClusterR)
dat <- select(adm, Outcome, ct:aHMCSF) %>%
    drop_na() #%>% write.csv("lasvData.csv", row.names = F)

X = dat %>% select(-Outcome) %>% mutate_all(list(log2))
y = dat$Outcome

dat = center_scale(X, mean_center = T, sd_scale = T)  # centering and scaling the data
gmm = GMM(dat, 8, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 100,
          em_iter = 100, verbose = F)

as.data.frame(dat) %>% gather() %>%
    ggplot(aes(x="value", group="key")) + geom_histogram(stat="count", alpha = .6)
# predict centroids, covariance matrix and weights

pr = predict_GMM(dat, gmm$centroids, gmm$covariance_matrices, gmm$weights)
opt_gmm = Optimal_Clusters_GMM(dat, max_clusters = 180, criterion = "BIC",

                               dist_mode = "maha_dist", seed_mode = "random_subset",

                               km_iter = 100, em_iter = 100, var_floor = 1e-10,

                               plot_data = T)

res = external_validation(as.numeric(y), pr$cluster_labels,

                          method = "adjusted_rand_index", summary_stats = T)

res
plotGMM::plot_GMM(gmm, k=8)

pca_dat = stats::princomp(dat)$scores[, 1:2]
gmm = GMM(pca_dat, 8, dist_mode = "maha_dist", seed_mode = "random_subset", km_iter = 10,
          em_iter = 10, verbose = F)



plot_2d(pca_dat, km$clusters, km$centroids)
gmm$covariance_matrices
pr$cluster_proba
km = KMeans_rcpp(pca_dat, clusters = 2, num_init = 5, max_iters = 100)


dis.lasv <- makemultdata(dat$PAI1, dat$TM,
                         cuts = median(c(dat$PAI1, dat$TM)))


mixmdl <- mixtools::multmixEM(as.matrix(dat[,2:3]), k = 2)
# mixtools::normalmixEM(dat$PAI1, k = 2)
plot_cut_point(mixmdl, plot = T, color = "wesanderson")
plot_mm(mixmdl, k = 2)


A <- c(1.51, 1.92, 1.08, 2.04, 2.14, 1.76, 1.17)
B <- c(1.69, 0.64, .9, 1.41, 1.01, .84, 1.28, 1.59)
C <- c(1.56, 1.22, 1.32, 1.39, 1.33, 1.54, 1.04, 2.25, 1.49)
D <- c(1.3, .75, 1.26, .69, .62, .9, 1.2, .32)
E <- c(.73, .8, .9, 1.24, .82, .72, .57, 1.18, .54, 1.3)

dis.coal <- makemultdata(A, B, C, D, E,
                         cuts = median(c(A, B, C, D, E)))
em.out <- multmixEM(dis.coal)
em.out[1:4]



###--------------------------
###
cases <- cyto_all_patients_dpi_corrected %>%
    select(id, dpi, Outcome, ct, !!cytokines, !!clin.chem) %>%
    filter(Outcome %in% c("Fatal", "Survivor"))

thresh <- topHits %>% select(Threshold = thresh.val, var=Factor) %>%
    rbind(topHitsClin %>%
              select(Threshold = thresh.val, var=Factor) %>%
              mutate(var = str_replace_all(var, pattern = "PAI1", replacement = "PAI1clin"))
    )

##for ease of viewing, cutting off at 12dpi --
##change this to see full range for survivors

long <- cases %>%
    group_by(id) %>%
    gather(key = "var", value = "Value", -Outcome, -id, -dpi)%>%
    left_join(thresh, by="var")%>%
    mutate(Factor = get.pretty.name(var),
           Units = get.pretty.factor(var),
           dpi = floor(dpi)) %>%
    na.omit() %>%
    group_by(id, var) %>%
    filter(dpi < 12, Value >= .1) %>%
    arrange(id)


survive.severe <- long %>% filter(dpi == 0 & Outcome=='Survivor') %>%
    filter(var=="ct" & Value<Threshold | var!="ct" & Value>Threshold)%>%
    select(id, var)%>%
    left_join(long, by = c("id", "var"))%>%
    add_column(Timing = "Predicted Fatal") %>%
    select(id, Outcome, dpi, var, Factor, Value, Threshold, Units, Timing) %>%
    data.frame()

survive.mild <- long %>% filter(dpi == 0 & Outcome=='Survivor') %>%
    filter(var=="ct" & Value>Threshold | var!="ct" & Value<Threshold)%>%
    select(id, var)%>%
    left_join(long, by = c("id", "var"))%>%
    add_column(Timing = "Predicted Survivor") %>%
    select(id, Outcome, dpi, var, Factor, Value, Threshold, Units, Timing) %>%
    data.frame()

fatal.late <- sqldf("SELECT * FROM long WHERE id IN (SELECT DISTINCT id FROM long WHERE dpi > 3 and Outcome='Fatal')") %>%
    add_column(Timing = "Late (>3 days)") %>%
    select(id, Outcome, dpi, var, Factor, Value, Threshold, Units, Timing) %>%
    data.frame()

fatal.early <- long %>% filter(Outcome == "Fatal", !id %in% fatal.late$id) %>%
    add_column(Timing = "Early (<=3 days)") %>%
    select(id, Outcome, dpi, var, Factor, Value, Threshold, Units, Timing) %>%
    data.frame()

cases.plot.survive <- rbind(survive.severe, survive.mild)
cases.plot.fatal <- rbind(fatal.late, fatal.early)
cases.plot.all <- rbind(cases.plot.fatal, cases.plot.survive) %>% ungroup()



time.data <- longitudinalPlots()

time.plot <- function(time.data) {
    ggplot(time.data, aes(x = dpi, y = Value, color = Timing )) +
        color_palette(c("red", "red4", "blue", "darkblue")) +
        theme_pubr() +
        facet_grid(. ~ Timing#, labeller =label_wrap_gen(width = 13, multi_line = T)
        ) +
        geom_point(alpha=0.7) +
        theme(text = element_text(size = 20),
              legend.position="none",
              strip.background = element_blank(),
              strip.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              # axis.ticks.x = element_blank(),
              axis.text.x = element_blank()
        ) +
        xlim(0,11)+
        labs(y=time.data$Units) +
        geom_hline(aes(yintercept = Threshold), linetype="dashed", color="black") +
        geom_smooth(method = "loess", color="black", formula = y~x)
}


l <- cyto_all_patients_dpi_corrected %>% group_by(id) %>%
    gather(key = "var", value = "Value", -Outcome, -id, -dpi)%>%
    # left_join(thresh, by="var")%>%
    mutate(Factor = get.pretty.name(var),
           Units = get.pretty.factor(var),
           dpi = floor(dpi)) %>%
    na.omit() %>%
    group_by(id, var) %>%
    filter(dpi < 12, Value >= .1) %>%
    arrange(id)

fatal.late <- sqldf("SELECT * FROM l WHERE id IN
                    (SELECT DISTINCT id FROM long WHERE dpi > 3 and Outcome='Fatal')") %>%
    # add_column(Timing = "Late (>3 days)") %>%
    select(id, Outcome, dpi, var, Factor, Value, Units) %>%
    mutate(Value = as.numeric(Value)) %>%
    data.frame()


fatal.late %>% filter(var %in% c("PAI1", "Fibrinogen", "DDIMER", "ct", "AntithrombinIII"), !is.na(Value)) %>%
    ggplot(aes(x = dpi, y = Value , color = id)) +
    # color_palette(c("red", "red4", "blue", "darkblue")) +
    theme_pubr() +
    facet_wrap(. ~ Factor, scales = "free"#, labeller =label_wrap_gen(width = 13, multi_line = T)
    ) +
    scale_y_continuous(trans = "log10") +
    geom_point(alpha=0.7) +
    geom_line() +
    # theme(text = element_text(size = 20),
    #       legend.position="none",
    #       strip.background = element_blank(),
    #       strip.text.x = element_blank(),
    #       axis.title.x = element_blank(),
    #       axis.title.y = element_blank(),
    #       # axis.ticks.x = element_blank(),
    #       axis.text.x = element_blank()
    # ) +
    # labs(y=time.data$Units) +
    # geom_hline(aes(yintercept = Threshold), linetype="dashed", color="black") +
# geom_smooth(method = "loess", color="black", formula = y~x)
    xlim(0,11)


#--------------
# Stats
#
#


ttest <- function(x) {apply(x, 2, FUN = t.test) %>%
        sapply(function(x) {
            c(x$estimate[1],
              ci.lower = x$conf.int[1],
              ci.upper = x$conf.int[2])
        }
        ) %>% data.frame()%>%
        mutate_if(is.numeric, exp) %>%
        mutate_if(is.numeric, round, 3)  %>%
        rbind(Missing = missing(x)) %>%
        t() %>% data.frame()%>%
        set_colnames(c("Mean", "CI.low", "CI.high", "Missing"))
}

stat.table <- function(df, cytokines, out1, out2) {

    lc <- log.cyto(df, cytokines) %>% filter(Outcome %in% c(!!out1, !!out2))

    p1 <- lc %>% filter(Outcome %in% out1) %>% select(!!cytokines) %>% ttest()
    p2 <- lc %>% filter(Outcome %in% out2) %>% select(!!cytokines) %>% ttest()

    p.2samp.unequal <- sapply(select(lc, all_of(cytokines)), function(y) wilcox.test(y~lc$Outcome))[3,] %>%
        data.frame() %>%
        mutate_if(is.numeric, round, 4) %>%
        set_rownames(c("p.value")) %>%
        gather(key = "Var", value = "p.value") %>%
        mutate(p.adj = p.adjust(p.value, method = "bonferroni", n = length(cytokines)))

    table1 <- data.frame(factor=cytokines %>% get.pretty.factor(),
                         fatal=sprintf("%.1f (%.1f - %.1f) [%.0f/%.0f]",
                                       p2$Mean,
                                       p2$CI.low,
                                       p2$CI.high,
                                       p2$Missing,
                                       nrow(lc %>% filter(Outcome=="Fatal"))),
                         surv=sprintf("%.1f (%.1f - %.1f) [%.0f/%.0f]",
                                      p1$Mean,
                                      p1$CI.low,
                                      p1$CI.high,
                                      p1$Missing,
                                      nrow(lc %>% filter(Outcome=="Survivor"))),
                         # p.value=p.2samp.unequal$p.value,
                         p.adj=p.2samp.unequal$p.adj
    ) %>%
        set_colnames(c("Factor",
                       paste(out2, " (95% CI) [missing/total]"),
                       paste(out1, " (95% CI) [missing/total]"),
                       # "p.value",
                       "p.adj")) %>%
        mutate(#p.value = sub(pattern="^0$", "<0.0001", p.value),
            p.adj = case_when(p.adj == 0 ~ "<0.0001",
                              # p.adj >= 0.1 ~ "NS",
                              TRUE ~ as.character(p.adj))
            # p.adj = sub(pattern="^0$", "<0.0001", p.adj),
            # p.adj = sub(pattern = "^1$", "NS", p.adj)
        )

    return(table1)
}
# "Other Febrile Illness"
out1 <- "Healthy Control"
out2 <- "Fatal" # wilcox test needs two levels...may need new grouping variable of Lassa vs non-LF
stat.table(pre.box.data, cytokines = c(cytokines), out1, out2)



nas <- pre.box.data %>% filter(Outcome %in% "Other Febrile Illness", !is.na(TNFa))

bun.cre <- adm %>%
    mutate(bun.cre = as.numeric(BUN/CRE)) %>%
    select(Outcome, bun.cre)

# Age---------------------------------------
#
#
adm.over40 <- adm %>% filter(Age >40)



age.sig <- adm.over40 %>%
    stat.table(cytokines = c("ct", cytokines))%>%
    filter(p.adj < 0.05) %>%
    arrange(p.adj)

age.factors <- age.sig$Factor %>% as.character()



glmDat <- glmFits(cases.df = adm.over40, c("ct", cytokines))

glmCo <- glmDat %>%
    glmCoefs(cytokines= c("ct", cytokines))

rocData <- glmDat %>%
    rocFits(cytokines= c("ct", cytokines))

names(rocData) <- c("ct", cytokines)

rocRes <- rocCoords(rocFits = rocData, glmCoefs = glmCo, cytokines= c("ct", cytokines))



# Late Fatals
fatal.late.ids <- sqldf("SELECT DISTINCT id FROM long WHERE dpi > 3 and Outcome='Fatal'") %>% unlist()

fatal.late.adm <- adm_all %>% filter(id %in% fatal.late.ids)


lc <- log.cyto(fatal.late.adm, cytokines)%>% filter(Outcome %in% c("Survivor", "Fatal"))

p.over40 <- lc %>% filter(Age>=40) %>% select(!!cytokines) %>% ttest()
p.under40 <- lc %>% filter(Outcome=="Fatal") %>% select(!!cytokines) %>% ttest()


sapply(select(lc, all_of(cytokines)), function(y) glm(y~lc$Age))


#Pairwise logistic regression ---------------
#
cytokines = colnames(adm %>% select(ct:RANTES))




glmDat <- glmPairFits(adm, cytokines)

glmCo <- glmDat %>%
    glmCoefs(cytokines)

rocData <- glmDat %>%
    rocFits(cytokines)

# names(rocData) <- cytokines


rocRes <- rocCoords(rocFits = rocData, glmCoefs = glmCo, cytokines)

rocResFormat <- rocRes %>% rownames_to_column("Var") %>%
    left_join(units, by = "Var") %>%
    select(Name, AUC, Threshold=thresh.val, Units, p.value,
           specificity, sensitivity, ppv, npv) %>%
    mutate_if(is.numeric, round, 3)




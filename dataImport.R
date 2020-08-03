pacman::p_load(tidyverse, readxl)



##DPI correction
cytokines <- c("TNFa", "sTNFRI", "sTNFRII", "sCD40L", "sFAS",
               "AntithrombinIII", "DDIMER",  "Fibrinogen", "PAI1", "TM", "tPA", "vWF",
               "EGF", "Pecam1", "sESelectin", "PSELECTIN", "siCAM1", "sVCAM1", "SAANG", "VEGF",
               "CFH",  "CRP", "GMCSF", "aHMCSF", "FRACTALINE",  "GRO", "IFNa2", "IFNy",
               "IL1RA", "sIL2Ra", "IL6", "IL8", "IL10", "IP10",
               "MCP1", "MCP2", "MCP3", "MIP1a", "MIP1b", "RANTES")

clin.chem <- c("ALT",  "AST", "BUN", "CRE")
cyto_all_patients_dpi_corrected <- read_xlsx("cyto_all_patients_dpi_corrected.xlsx") %>%
    mutate(Outcome = factor(Outcome, levels = c("Fatal",
                                                "Survivor",
                                                "Other Febrile Illness",
                                                "Healthy Control")))

primaryIDs <- read_rds("primaryIDs.RDS") #saved from 2015_data_20200909

cyto2015 <- cyto_all_patients_dpi_corrected %>%
    filter(id %in% primaryIDs)
# filter(date<as.Date("2016-05-31") | id=="314LV16", !is.na(ct))
#   #patient 314LV16 from this dataset had no date

cyto2017 <- cyto_all_patients_dpi_corrected %>%
    filter(!id %in% primaryIDs)
# filter(date>=as.Date("2016-05-31") | id=="151LV17", !is.na(ct))
#   #patient 151LV17 from this dataset had no date


###-----Select Data Functions----


# Which samples to consider
day0.select <- function(cyto.data, outcomes=c("Survivor", "Fatal") ) {
    cyto.data %>%
        select(id, dpi, date, Age, Outcome, ct, !!cytokines) %>%
        group_by(id) %>%
        filter(Outcome %in% outcomes, dpi==0) %>%
        ungroup() %>%
        as.data.frame()
}

clin.adm.select  <- function(cyto.data, outcomes=c("Survivor", "Fatal") ){
    cyto.data %>%
        select(id, dpi.clin=dpi, Outcome, ct, !!cytokines, !!clin.chem) %>%
        group_by(id) %>%
        filter(!is.na(AST), Outcome %in% outcomes, dpi.clin <= 1) %>%
        slice(which.min(dpi.clin)) %>%
        ungroup() %>%
        as.data.frame()
}

peak.select <- function(cyto.data, outcomes=c("Survivor", "Fatal")) {
    cyto.data %>%
        select(id, date, dpi, Outcome, ct, !!cytokines, !!clin.chem) %>%
        group_by(id)%>%
        filter(Outcome %in% outcomes) %>%
        slice(which.min(ct)) %>%
        ungroup()%>%
        as.data.frame()
}

last.day.select <- function(dataframe, outcomes=c("Survivor", "Fatal")){
    dataframe %>%
        select(id, date, Age, Outcome, ct, !!cytokines, !!clin.chem) %>%
        group_by(id) %>%
        filter(Outcome %in% outcomes, date==dplyr::last(date), !is.na(ct)) %>%
        ungroup() %>%
        as.data.frame()
}

count.outcome <- function(cyto.data) {
    cyto.data %>%
        group_by(Outcome) %>%
        dplyr::count(name = "Samples", Patients = n_distinct(id)) %>%
        arrange(Outcome)
}

missing <- function(df) {
    df %>%
        ungroup() %>%
        summarise_all(list(~sum(is.na(.)))) %>%
        data.frame()
}


log.cyto <- function(cyto.data, cytos )  { cyto.data %>%
        mutate_at(.vars = cytos[!cytos %in% "ct"], .funs = list(log)) %>%
        rationalize()
}


###-----Subset Data----


adm2015_clin <- cyto2015  %>% clin.adm.select(c("Survivor", "Fatal")) %>% arrange(Outcome)

adm2015 <- cyto2015  %>%
    day0.select(c("Survivor", "Fatal")) %>%
    arrange(Outcome) %>%
    left_join(select(adm2015_clin, id, Outcome, !!clin.chem), by = c("id", "Outcome"))

adm2017_clin <- cyto2017  %>% clin.adm.select(c("Survivor", "Fatal")) %>% arrange(Outcome)

adm2017 <- cyto2017 %>%
    day0.select(c("Survivor", "Fatal")) %>%
    arrange(Outcome)%>%
    left_join(select(adm2017_clin, id, Outcome, !!clin.chem), by = c("id", "Outcome"))

adm_all <- rbind(adm2015, adm2017)
adm_clin <- rbind(adm2015_clin, adm2017_clin)


peak_all <- cyto_all_patients_dpi_corrected %>%
    peak.select(c("Survivor", "Fatal")) %>% arrange(Outcome)
peak2015 <- cyto2015  %>% peak.select(c("Survivor", "Fatal")) %>% arrange(Outcome)
peak2017 <- cyto2017 %>% peak.select(c("Survivor", "Fatal")) %>% arrange(Outcome)




###-----Units----


units <- read_xlsx("units.xlsx")

get.pretty.factor <- function(cytokines) {
    lapply(cytokines, function(x) as.character(units$Factor[match(x, units$Var)])) %>%
        unlist()
}
get.pretty.name <- function(cytokines) {
    lapply(cytokines, function(x) as.character(units$Name[match(x, units$Var)])) %>%
        unlist()
}
get.pretty.units <- function(cytokines) {
    lapply(cytokines, function(x) as.character(units$Units[match(x, units$Var)])) %>%
        unlist()
}

ttest <- function(x) {apply(x, 2, FUN = t.test) %>%
        sapply(function(x) {
            c(x$estimate[1],
              ci.lower = x$conf.int[1],
              ci.upper = x$conf.int[2])
        }
        ) %>% data.frame()%>%
        mutate_at(vars(-ct), exp) %>%
        mutate_if(is.numeric, round, 3)  %>%
        rbind(Missing = missing(x)) %>%
        t() %>% data.frame()%>%
        set_colnames(c("Mean", "CI.low", "CI.high", "Missing"))
}



###-----Stats----



# lc <- log.cyto(adm2015, cytokines)%>% filter(Outcome %in% c("Survivor", "Fatal"))
stat.table <- function(df, cytokines) {

    lc <- log.cyto(df, cytokines)%>% filter(Outcome %in% c("Survivor", "Fatal"))

    p.Surv <- lc %>% filter(Outcome=="Survivor") %>% select(!!cytokines) %>% ttest()
    p.Fatal <- lc %>% filter(Outcome=="Fatal") %>% select(!!cytokines) %>% ttest()

    p.2samp.unequal <- sapply(select(lc, all_of(cytokines)), function(y) wilcox.test(y~lc$Outcome))[3,] %>%
        data.frame() %>%
        mutate_if(is.numeric, round, 4) %>%
        set_rownames(c("p.value")) %>%
        gather(key = "Var", value = "p.value") %>%
        mutate(p.adj = p.adjust(p.value, method = "bonferroni", n = length(cytokines)))

    table1 <- data.frame(factor=cytokines,
                         fatal=sprintf("%.1f (%.1f - %.1f) [%.0f/%.0f]",
                                       p.Fatal$Mean,
                                       p.Fatal$CI.low,
                                       p.Fatal$CI.high,
                                       p.Fatal$Missing,
                                       nrow(lc %>% filter(Outcome=="Fatal"))),
                         surv=sprintf("%.1f (%.1f - %.1f) [%.0f/%.0f]",
                                      p.Surv$Mean,
                                      p.Surv$CI.low,
                                      p.Surv$CI.high,
                                      p.Surv$Missing,
                                      nrow(lc %>% filter(Outcome=="Survivor"))),
                         # p.value=p.2samp.unequal$p.value,
                         p.adj=p.2samp.unequal$p.adj
    ) %>%
        set_colnames(c("Factor",
                       "Fatal (95% CI) [missing/total]",
                       "Survivor (95% CI) [missing/total]",
                       # "p.value",
                       "p.adj")) %>%
        mutate(#p.value = sub(pattern="^0$", "<0.0001", p.value),
            p.adj = case_when(p.adj == 0 ~ "<0.0001",
                              p.adj >= 0.1 ~ "NS",
                              TRUE ~ as.character(p.adj))
            # p.adj = sub(pattern="^0$", "<0.0001", p.adj),
            # p.adj = sub(pattern = "^1$", "NS", p.adj)
        )

    return(table1)
}

table1 <- stat.table(adm2015, cytokines = c("ct", cytokines, clin.chem))

supp.table1 <- stat.table(adm2017, cytokines = c("ct", cytokines, clin.chem))

# stat.table(cytoIntake2015, cytokines)


sig.vars.padj <- data.frame(var = table1$Factor[which(table1$p.adj < 0.05 & supp.table1$p.adj < 0.05)])
# adm_all <- adm_all %>% dplyr::rename_with(.cols = ct:CRE, .fn = get.pretty.name)


###-----Plotting----


cytoTop5 <- c('ct', 'PAI1', 'aHMCSF', 'TM', 'IL8', 'sTNFRI')

colMax <- function(cytoData) sapply(cytoData, max, na.rm = TRUE)

colMaxTop5 <- cyto_all_patients_dpi_corrected %>%
    select(!!cytoTop5) %>%
    colMax() %>%
    data.frame() %>%
    rownames_to_column() %>%
    set_colnames(c("Factor", "Max"))


cytoColMax <- cyto_all_patients_dpi_corrected %>%
    select(ct:CRE) %>%
    colMax() %>%
    data.frame() %>%
    rownames_to_column() %>%
    set_colnames(c("Factor", "Max"))


box.data <- cyto_all_patients_dpi_corrected %>%
    filter(Outcome=="Healthy Control" | (Outcome=="Other Febrile Illness" & dpi==0)) %>%
    ungroup() %>%
    select(Outcome, ct, !!cytokines, !!clin.chem) %>%
    rbind(select(adm_all, Outcome, ct, !!cytokines, !!clin.chem)) %>%
    gather(key = "Factor", value = "Value", -Outcome) %>%
    left_join(cytoColMax, by="Factor")%>%
    filter(!is.na(Value)) %>%
    mutate(Name = get.pretty.name(Factor),
           Units = get.pretty.units(Factor),
           Name_Units = get.pretty.factor(Factor),
           Outcome = factor(Outcome, levels = c("Fatal", "Survivor", "Other Febrile Illness", "Healthy Control"))) %>%
    mutate_if(is.numeric, round, 2)

#https://stackoverflow.com/questions/6461209/how-to-round-up-to-the-nearest-10-or-100-or-x
roundUpNice <- function(x, nice=c(1,2,4,5,6,8,10)) {
    if(length(x) != 1) stop("'x' must be of length 1")
    10^floor(log10(x)) * nice[[which(x <= 10^floor(log10(x)) * nice)[[1]]]]
}


my_comparisons <- list(c("Fatal", "Survivor"),
                       c("Fatal", "Other Febrile Illness"),
                       c("Fatal", "Healthy Control"),
                       c("Survivor", "Other Febrile Illness"),
                       c("Survivor", "Healthy Control"))



###-----Modeling----


glmSingle <- function(cases, cytokines) {

    fmla <- sapply(1:length(cytokines), function(x, var = cytokines) {
        fmla <- as.formula(paste0("as.factor(Outcome) ~ ", var[x]))
    })
    fits <- lapply(fmla, function(x, dfcases = cases) {
        glm(x, data=dfcases, family = binomial, na.action = na.omit)
    })

    p.value <- sapply(fits, function(x) summary(x)$coefficients [8])
    # p.age <- sapply(fits, function(x) summary(x)$coefficients[11] )
    roc.fits <- lapply(fits, function(x)
        roc(response = x$y, predictor = x$fitted.values, na.rm = T,
            auc = TRUE, levels = c(0, 1),
            direction = "auto", quiet = T) )

    fit.intercept <- sapply(1:length(cytokines),
                            function(x) (coef(fits[[x]])[1] ) )
    # fit.age <- sapply(1:length(cytokines),
    # function(x) (coef(fits[[x]])[2] ) )
    fit.cytokine <- sapply(1:length(cytokines),
                           function(x) (coef(fits[[x]])[2] ) )

    AUC <- sapply(roc.fits, function(x) x$auc )
    #
    # find optimal accuracy of ROC fit
    roc.coords <- sapply(roc.fits, function(y)
        coords(y, ret=c("threshold", "specificity", "sensitivity",
                        "accuracy", "ppv", "npv", "tpr", "tnr", "fnr", "fpr", "fdr",
                        "tn", "tp", "fn", "fp"), transpose=T,
               x="best", best.method="c"), simplify = T) %>%
        t() %>% data.frame()%>%
        cbind(fit.intercept, fit.cytokine, AUC, p.value) %>%
        mutate(thresh.val = (log(threshold/(1-threshold)) - fit.intercept) / fit.cytokine) %>%
        set_rownames(cytokines)

    # Threshold <- sapply(1:length(cytokines),
    #                     function(x) (log(opt_t[x]/(1-opt_t[x])) - coef(fits[[x]])[1] ) / coef(fits[[x]])[2]
    #                     )

    return(roc.coords)
}
#
sig.cytokines <- sig.vars.padj$var %>% as.character()

#Run glmSingle with correctly formated dataframe
topGLM <- function(cases, sig.cytokines) { cases %>%
        select(Outcome, !!sig.cytokines) %>%
        filter(Outcome %in% c("Survivor", "Fatal")) %>%
        glmSingle(cytokines = sig.cytokines) %>%
        rownames_to_column(var = "Factor") %>%
        mutate(AUC=round(AUC, digits=3),
               p.value = round(p.value, digits = 4),
               p.value = sub(pattern="^0$", "<0.0001", p.value),
               Threshold=round(x = thresh.val, digits = 1),
               Threshold=str_c(Threshold, get.pretty.units(sig.cytokines), sep = " "),
               auc.rank = row_number(desc(AUC)),
               MCC = mcc(TP = tp, FP = fp, TN = tn, FN = fn))  %>%
        arrange(Factor)
}


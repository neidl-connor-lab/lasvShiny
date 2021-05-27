

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
        select(id, dpi, Outcome, !!clin.chem) %>%
        group_by(id) %>%
        filter(!is.na(AST), Outcome %in% outcomes, dpi <= 1) %>%
        slice(which.min(dpi)) %>%
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



###-----Units-----

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


# get.pat.meta <- function(cytokines) {
#     lapply(cytokines, function(x) as.character(units$Factor[match(x, units$Var)])) %>%
#         unlist()
# }

###-----Stats----



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
                              # p.adj >= 0.1 ~ "NS",
                              TRUE ~ as.character(p.adj))
            # p.adj = sub(pattern="^0$", "<0.0001", p.adj),
            # p.adj = sub(pattern = "^1$", "NS", p.adj)
        )

    return(table1)
}


colMax <- function(cytoData) sapply(cytoData, max, na.rm = TRUE)

cyto.corr <- function(cases, rects) {
    correlations <- abs(cor(cases, method = "pearson", use = "complete.obs"))

    corrplot(correlations,
             cl.lim=c(0,max(correlations)),
             method="circle",
             order = "hclust",
             hclust.method = "complete",
             tl.col="black",
             tl.srt=45,
             addrect = rects,
             is.corr=FALSE)
}


###-----Modeling----


glmFits <- function(cases.df, cyto, ...) {

    fmla <- sapply(1:length(cyto), function(x, var = cyto) {
        as.formula(paste0("as.factor(Outcome) ~ ", var[x]))
    })


    glmFits <- lapply(fmla, function(y) {
        glm(y, data=cases.df, family = binomial, na.action = na.omit)
    })

    return(glmFits)
}




glmCoefs <-  function(glmFits, cytokines=cytokines, ...)  {
    # get glm coefficients and significance from fits (glmFits)
    fit.intercept <- sapply(1:length(cytokines),
                            function(x) (coef(glmFits[[x]])[1] ) )
    fit.cytokine <- sapply(1:length(cytokines),
                           function(x) (coef(glmFits[[x]])[2] ) )
    p.value <- sapply(glmFits, function(x) summary(x)$coefficients [8])

    glmCoefs <- cbind(fit.intercept, fit.cytokine, p.value) %>%
        set_rownames(cytokines)

    return(glmCoefs)
    }

# find optimal accuracy of ROC fit
# arguments:
# 1. model fits generated by glmFits
# 2. list of cytokines

rocFits <- function(glmFits, cytokines=cytokines, ...) {

    # calculate the ROC curve
    rocFits <- lapply(glmFits, function(x)
        roc(response = x$y, predictor = x$fitted.values, na.rm = T,
            # auc = TRUE, levels = c(0, 1), direction = "auto",
            quiet = T) ) #%>% names(cytokines)

    return(rocFits)
}

rocCoords <- function(rocFits, glmCoefs, cytokines=cytokines, ...) {
    # AUCs of each ROC curve
    AUC <- sapply(rocFits, function(x) x$auc )

    # coordinates of the best (top-left corner) threshold of ROC curve
    # and measures of prediction accuracy
    rocCoords <- sapply(rocFits, function(y)
        coords(y, ret=c("threshold", "specificity",
                        "sensitivity", "accuracy",
                        "ppv", "npv", "tpr", "tnr",
                        "fnr", "fpr", "fdr",
                        "tn", "tp", "fn", "fp"),
               transpose=T,
               x="best",
               best.method="c",
               best.weights=c(5, .2)), #NEED TO TEST WEIGHTS
         simplify = T) %>%
        data.frame() %>%
        t() %>% data.frame() %>%
        cbind(AUC, glmCoefs) %>%
        mutate(thresh.val = (log(threshold/(1-threshold)) -
                                 fit.intercept) / fit.cytokine) %>%
        set_rownames(cytokines)

    return(rocCoords)
}

# Threshold <- sapply(1:length(cytokines),
#                     function(x) (log(opt_t[x]/(1-opt_t[x])) - coef(fits[[x]])[1] ) / coef(fits[[x]])[2]
#                     )
#


# USAGE---------------------
# glmDat <- glmFits(adm, cytokines)
#
# glmCo <- glmDat %>%
#     glmCoefs(cytokines)
#
# rocData <- glmDat %>%
#     rocFits(cytokines)
#
# names(rocData) <- cytokines
#
# rocRes <- rocCoords(rocFits = rocData, glmCoefs = glmCo, cytokines)
# ---------------------------


#Same as glmFits, but for pairs of cytokines
glmPairFits <- function(cases, cytokines, ...) {
    # lapply(1:(length(cytokines)-1) (x) {
    #     sapply((x+1 : length(cytokines) (y) {

    fmla <- sapply(1 : (length(cytokines)-1), function(x, var = cytokines) {
        sapply( (x+1) : length(cytokines),function(y, var = cytokines) {
            fmla <- as.formula(paste0("as.factor(Outcome) ~ ", var[x], "+", var[y]))
        })
    })%>% unlist()


    glmPairFits <- lapply(fmla, function(x, dfcases = adm_all) {
        glm(x, data=dfcases, family = binomial, na.action = na.omit)
    })

#     return(glmPairFits)
# }
#
#
#
# glmPairCoefs <-  function(glmPairFits, cytokines=cytokines, ...)  {
#     # get glm coefficients and significance from fits (glmFits)
#     # fit.intercept <- sapply(1:length(glmPairFits),
#     #                         function(x) (coef(glmPairFits[[x]])[1] ) )
#     # fit.cytokine1 <- sapply(1:length(glmPairFits),
#     #                        function(x) (coef(glmPairFits[[x]])[2] ) )
#     # fit.cytokine2 <- sapply(1:length(glmPairFits),
#     #                         function(x) (coef(glmPairFits[[x]])[3] ) )


    ## Not working ---- need this for threshold values
    # glmCoefs <- sapply(glmPairFits, function(x) {
        cytokine1 <- sapply(1:length(fmla), function(x) {
                            fmla[[x]][[3]][[2]]}) %>% as.character()  # %>% as.vector()

        cytokine2 <- sapply(1:length(fmla), function(x) {
                            fmla[[x]][[3]][[3]]}) %>% as.character()
    #     fit.intercept <- coef(glmPairFits[[x]])[1]
    #     fit.cytokine1 <- coef(glmPairFits[[x]])[2]
    #     fit.cytokine2 <- coef(glmPairFits[[x]])[3]
    #     })
    #
    #     data.frame() %>%
    #         t() %>%
    #         set_colnames(c("cytokine1", "cytokine2", "fit.intercept", "fit.cytokine1", "fit.cytokine2"))
    #



    # p.value <- sapply(glmPairFits, function(x) summary(x)$coefficients [8])
#
#     return(glmCoefs)
# }
#
#
# rocPairCoords <- function(rocFits, glmCoefs, cytokines=cytokines, ...) {
    # AUCs of each ROC curve
    rocFits <- lapply(glmPairFits, function(x)
        roc(response = x$y, predictor = x$fitted.values, na.rm = T,
            # auc = TRUE, levels = c(0, 1), direction = "auto",
            quiet = T) ) #%>% names(cytokines)

    AUC <- sapply(rocFits, function(x) x$auc )

    # coordinates of the best (top-left corner) threshold of ROC curve
    # and measures of prediction accuracy
    rocCoords <- sapply(rocFits, function(y)
        coords(y, ret=c(#"threshold",
                        "specificity",
                        "sensitivity", "accuracy",
                        "ppv", "npv", "tpr", "tnr",
                        "fnr", "fpr", "fdr",
                        "tn", "tp", "fn", "fp"),
               transpose=T,
               x="best",
               best.method="c"),
        simplify = T) %>%

        t() %>%
        data.frame() %>%
        cbind(AUC, cytokine1, cytokine2)

    return(rocCoords)
}


# glmPairs <- glmPairFits(cases = adm, cytokines = cytokines)



#
# glmPairDat <- glmPairFits(adm, cytokines)
#
# glmPairCo <- glmPairDat %>%
#     glmPairCoefs(cytokines)
#
# rocPairData <- glmPairDat %>%
#     rocFits(cytokines)
#
# # names(rocData) <- cytokines
#
# rocRes <- rocPairCoords(rocFits = rocPairData, glmCoefs = glmPairCo, cytokines)
#



#Run glmSingle with correctly formated dataframe
# topGLM <- function(cases, sig.cytokines) { cases %>%
#         select(Outcome, !!sig.cytokines) %>%
#         filter(Outcome %in% c("Survivor", "Fatal")) %>%
#         glmSingle(cytokines = sig.cytokines) %>%
#         rownames_to_column(var = "Factor") %>%
#         mutate(AUC=round(AUC, digits=3),
#                p.value = round(p.value, digits = 4),
#                p.value = sub(pattern="^0$", "<0.0001", p.value),
#                Threshold=round(x = thresh.val, digits = 1),
#                Threshold=str_c(Threshold, get.pretty.units(sig.cytokines), sep = " "),
#                auc.rank = row_number(desc(AUC)),
#                MCC = mcc(TP = tp, FP = fp, TN = tn, FN = fn))  %>%
#         arrange(Factor)
# }

epi.summary <- function(cyto.data) {
    # cyto.data = cyto2017
    cyto.data <- cyto.data %>%
        filter(Outcome %in% c("Survivor", "Fatal"))
    cyto.data$Outcome <- factor(cyto.data$Outcome, levels = c("Survivor", "Fatal"))

    epi.data <- cyto.data %>%
        select(id, Outcome, Age, Sex, State) %>%
        group_by(id) %>%
        slice(id=1)

    Age <- epi.data %>%
        # filter(!is.na(Age)) %>%
        mutate(Age = case_when(Age < 21 ~ 1,
                               21 <= Age & Age <= 40 ~ 2,
                               Age > 40 ~ 3,
                               is.na(Age) ~ 4)) %>%
        ungroup() %>% #group_by(Age) %>%
        #summarise(list)
        count(Outcome, Age) %>%
        spread(key = Age, value = n) %>%
        column_to_rownames(1) %>%
        t() %>% as.data.frame() %>%
        mutate(Total=Survivor+Fatal, "Percent Survived"=Survivor*100/Total) %>%
        add_column(Group=c("<21", "21-40", ">40", "Unknown"), .before = 1)


    Sex <- epi.data %>%
        # filter(!is.na(Sex)) %>%
        ungroup() %>% #group_by(Age) %>%
        #summarise(list)
        count(Outcome, Sex) %>%
        spread(key = Sex, value = n) %>%
        data.frame(., row.names = 1) %>%
        t() %>% data.frame(.) %>%
        mutate(Total=Survivor+Fatal, "Percent Survived"=Survivor*100/Total) %>%
        add_column(Group=c("Female", "Male", "Unknown"), .before = 1)


    State <- epi.data %>%
        # filter(State %in% c("EBONYI", "EDO", "ONDO")) %>%
        ungroup() %>% #group_by(Age) %>%
        #summarise(list)
        count(Outcome, State) %>%
        spread(key = State, value = n) %>%
        data.frame(., row.names = 1) %>%
        replace(is.na(.), 0) %>%
        mutate(Other = select(., -EDO, -ONDO, -EBONYI) %>% reduce(., `+`)) %>%
        select("Ebonyi"=EBONYI, "Edo"=EDO, "Ondo"=ONDO, "Other/Unknown"=Other) %>%
        t() %>% data.frame(.) %>% set_colnames(c("Survivor", "Fatal"))%>%
        mutate(Total=Survivor+Fatal, "Percent Survived"=Survivor*100/Total) %>%
        add_column(Group=c("Ebonyi", "Edo", "Ondo", "Other/Unknown"), .before = 1)


    results <- rbind(Age, Sex, State) %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        mutate_if(is.numeric, signif, 3)

    return(results)
}


longitudinalPlots <- function(){
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


    return(cases.plot.all)
}



pacman::p_load(tidyverse, readxl, hablar, magrittr)

source("glmFunctions.R")

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

cytoAll <- cyto_all_patients_dpi_corrected %>%
    mutate(Dataset = case_when(
        id %in% primaryIDs ~ "Primary",
        !id %in% primaryIDs ~ "Secondary"
            ) %>% factor(),
        date.onset=lubridate::date(date.onset),
        date=lubridate::date(date)
        ) %>%
    relocate(Dataset, .after=id)


metaPatients <- cytoAll %>%
    select(id, Dataset, Outcome, Age, Sex, State) %>%
    mutate(AgeBin = case_when(Age < 20 ~ "<20",
                           20 <= Age & Age <= 40 ~ "20-40",
                           Age > 40 ~ ">40",
                           is.na(Age) ~ "Unknown") %>%
               factor(levels = c("<20", "20-40", ">40", "Unknown")),
           Sex = case_when(is.na(Sex) ~ "Unknown",
                           TRUE ~ as.character(Sex)) %>% factor(),
           State = case_when(State %in% c("EDO", "ONDO", "EBONYI") ~ str_to_title(State),
                             is.na(State) ~ "Unknown",
                             TRUE ~ "Other") %>% factor())%>%
    distinct() %>% arrange(id, Age) %>% group_by(id) %>% slice(1)



###-----Subset Data----


adm2015_clin <- cyto2015 %>% clin.adm.select(c("Survivor", "Fatal")) %>% arrange(Outcome)

adm2015 <- cyto2015  %>%
    day0.select(c("Survivor", "Fatal")) %>%
    arrange(Outcome) %>%
    left_join(select(adm2015_clin, id, Outcome, !!clin.chem), by = c("id", "Outcome"))

adm2017_clin <- cyto2017 %>% clin.adm.select(c("Survivor", "Fatal")) %>% arrange(Outcome)

adm2017 <- cyto2017 %>%
    day0.select(c("Survivor", "Fatal")) %>%
    arrange(Outcome)%>%
    left_join(select(adm2017_clin, id, Outcome, !!clin.chem), by = c("id", "Outcome"))


adm_all <- rbind(adm2015, adm2017)
# admCrtl <- cyto_all_patients_dpi_corrected %>%
#     filter(Outcome=="Healthy Control" | (Outcome=="Other Febrile Illness" & dpi==0)) %>%
#     ungroup()

# cyto_all_patients_dpi_corrected %>% select(sFAS:IFNy) %>% filter(all(is.na(.)))
#
# ladmLASV <- rbind(adm2015, adm2017)
# admClin <- rbind(adm2015_clin, adm2017_clin)
# admAll <- rbind(admLASV, admCrtl)

peak_all <- cyto_all_patients_dpi_corrected %>%
    peak.select(c("Survivor", "Fatal")) %>% arrange(Outcome)
peak2015 <- cyto2015  %>% peak.select(c("Survivor", "Fatal")) %>% arrange(Outcome)
peak2017 <- cyto2017 %>% peak.select(c("Survivor", "Fatal")) %>% arrange(Outcome)




###-----Units----


units <- read_xlsx("units.xlsx")

table1 <- stat.table(adm2015, cytokines = c("ct", cytokines, clin.chem))

supp.table1 <- stat.table(adm2017, cytokines = c("ct", cytokines, clin.chem))

# stat.table(cytoIntake2015, cytokines)


sig.vars.padj <- data.frame(var = table1$Factor[which(table1$p.adj < 0.05 & supp.table1$p.adj < 0.05)])
# adm_all <- adm_all %>% dplyr::rename_with(.cols = ct:CRE, .fn = get.pretty.name)


###-----Plotting----


cytoTop5 <- c('ct', 'PAI1', 'aHMCSF', 'TM', 'IL8', 'sTNFRI')


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
    select(id, Outcome, ct, !!cytokines, !!clin.chem) %>%
    rbind(select(adm_all, id, Outcome, ct, !!cytokines, !!clin.chem)) %>%
    gather(key = "Factor", value = "Value", -id, -Outcome) %>%
    left_join(cytoColMax, by="Factor")%>%
    filter(!is.na(Value)) %>%
    mutate(Name = get.pretty.name(Factor),
           Units = get.pretty.units(Factor),
           Name_Units = get.pretty.factor(Factor),
           Outcome = factor(Outcome, levels = c("Fatal", "Survivor", "Other Febrile Illness", "Healthy Control"))) %>%
    mutate_if(is.numeric, round, 2)%>%
    left_join(metaPatients %>% select(-Outcome, -Age),
              by = "id")

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


sig.cytokines <- sig.vars.padj$var %>% as.character()



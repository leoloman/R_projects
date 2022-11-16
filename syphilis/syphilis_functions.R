relative_risk_function <- function(ref_pos, ref_neg, case_pos, case_neg){
  
  rr_matrix <-  matrix(c(ref_neg, case_neg, ref_pos,case_pos), nrow=2)
  
  rr_result <- riskratio(rr_matrix)$measure
  
  return(c(rr_result[2,2], rr_result[2,1], rr_result[2,3]))
}


load_clean_uk_pop_data <- function(sheetname){
  uk_age_pop_est <- read_xls('data/syphilis/mid-1991-unformatted-data-file.xls', sheet = sheetname, n_max=1)
  
  # convert transpose dataframe turning the column names into a column - dropping the first 4 rows
  # then convert the age into the age bands which we have in the syphilis dataset
  uk_age_pop_est <- as_tibble(cbind(Age = names(uk_age_pop_est)[4:length(uk_age_pop_est)], population = t(uk_age_pop_est)[4:length(uk_age_pop_est)])) %>% 
    mutate(Age = ifelse(Age <= 14,'<=14',
                 ifelse(Age <=19,'15-19',
                 ifelse(Age<=24,'20-24',
                 ifelse(Age<=29,'25-29',
                 ifelse(Age<=34,'30-34',
                 ifelse(Age<=44,'35-44',
                 ifelse(Age<=54,'45-54', '>55')))))))) %>% 
    group_by(across('Age')) %>% 
    summarise(population = sum(as.numeric(population)))
  return(uk_age_pop_est)
}


relative_risk_age <- function(df, gender = NA, ref_group = NA){
  # load population data by what gender is chosen
  if (is.na(gender)) {
    uk_age_pop_est <- load_clean_uk_pop_data('Mid-1991 Persons')
    plot_title <- 'Relative Risk split by age for all Persons'
  } 
  else if (gender == 'Male') {
    uk_age_pop_est <- load_clean_uk_pop_data('Mid-1991 Males')
    plot_title <- 'Relative Risk split by age for Males'
  } else {
    uk_age_pop_est <- load_clean_uk_pop_data('Mid-1991 Females')
    plot_title <- 'Relative Risk split by age for Females'
  }
  
  if (is.na(gender)){
    grouped_df <- df %>% group_by(across(Age)) %>% summarise(counts= n())
  } else {
    grouped_df <- df %>% filter(Sex == gender) %>% group_by(across(Age)) %>% summarise(counts= n())
  }
  
  uk_rr_age <- as.data.frame(uk_age_pop_est %>% left_join(grouped_df))
  
  uk_rr_age$population <- uk_rr_age$population - uk_rr_age$counts
    
  ref_group <- uk_rr_age[uk_rr_age$counts==max(uk_rr_age$counts),]
  
  # filt_uk_rr_age <- uk_rr_age[uk_rr_age$Age!= ref_group$Age,]
  
  rr_rr <- t(mapply(relative_risk_function, ref_group$counts, ref_group$population, 
                    uk_rr_age$counts, uk_rr_age$population))
  
  colnames(rr_rr) <- c('V1','V2','V3')
  
  uk_rr_age <- cbind(uk_rr_age, rr_rr)
  uk_rr_age[uk_rr_age$counts==max(uk_rr_age$counts),c('V1','V2','V3')] <- NA
  uk_rr_age$Age <- factor(uk_rr_age$Age, 
                          levels = c("<=14" ,"15-19" ,"20-24", "25-29", 
                                     "30-34", "35-44" ,"45-54",">55"))
  
  ggplot(uk_rr_age, aes(Age, V2)) + 
    geom_point(size=2) + 
    geom_errorbar(aes(ymin=V1, ymax=V3), width=.2) +
    geom_hline(yintercept=1) + 
    labs(y='Risk Ratio', x='Age', caption = paste0('Reference Group: ',ref_group$Age,'\nPositive Cases: ',ref_group$counts,'\nPopulation: ',ref_group$population)) + 
    ggtitle(plot_title) + 
    theme_bw()
}


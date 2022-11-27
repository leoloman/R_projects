relative_risk_function <- function(ref_pos, ref_neg, case_pos, case_neg){
  
  rr_matrix <-  matrix(c(ref_neg, case_neg, ref_pos,case_pos), nrow=2)
  
  rr_result <- riskratio(rr_matrix)$measure
  
  return(c(rr_result[2,2], rr_result[2,1], rr_result[2,3]))
}

load_clean_gender_age_pop_data <- function(group_by_vars, gender_filter = NA, eth_filter=NA,
                                           path = 'syphilis/data/gender_age_ethnicity_breakdown.csv'){
  ga_age_pop_est <- read.csv(path, skip  = 1)
  # clean data to match what is present in the syphilis data
  ga_age_pop_est <- ga_age_pop_est %>% mutate(age = str_trim(age), 
                                              age = recode(age, "under 5" = "<=14",
                                  "5 to 9" = "<=14",
                                  "10 to 14"="<=14",
                                  "15 to 19" = "15-19",
                                  "20 to 24" ="20-24",
                                  "25 to 29" = "25-29",
                                  "30 to 34" = "30-34",
                                  "35 to 39" = "35-44",
                                  "40 to 44"="35-44",
                                  "45 to 49"="45-54",
                                  "50 to 54"="45-54",
                                  "55 to 59"=">55" ,
                                  "60 to 64"=">55" ,
                                  "65 to 69"=">55" ,
                                  "70 to 74"=">55" ,
                                  "75 to 79"=">55" ,
                                  "80 to 84"=">55" ,
                                  "85<="=">55" ),
                            race = recode(race, 
                                          " American Indian, eskimo or aleut" = " Other",
                                          " asian or Pacific Islander" = " Other"),
                            race = str_to_title(race),
                            # remove white spaces 
                            race = gsub(' ','',race),
                            Gender = str_to_title(Gender)) %>% 
    group_by(Gender, age, race) %>% summarise(count = sum(count))
  
  if(is.na(gender_filter) != TRUE){
    ga_age_pop_est <- ga_age_pop_est%>% filter(Gender == gender_filter)
  }
  
  if(is.na(eth_filter) != TRUE){
    ga_age_pop_est <- ga_age_pop_est%>% filter(race == eth_filter)
  }
  
  ga_age_pop_est <- ga_age_pop_est %>% group_by(across(group_by_vars)) %>% summarise(count = sum(count))
  

  
  return(ga_age_pop_est)
}

relative_risk_age <- function(df, gender = NA, ref_group = NA, path = 'syphilis/data/gender_age_ethnicity_breakdown.csv'){
  # load population data by what gender is chosen
  if (is.na(gender)) {
    ga_age_pop_est <- load_clean_gender_age_pop_data(c('age'), path = path)
    plot_title <- 'Relative Risk split by age for all Persons'
    names(ga_age_pop_est) <- c('Age','population')
  } 
  else if (is.na(gender) == FALSE) {
    ga_age_pop_est <- load_clean_gender_age_pop_data(c("Gender",'age'),gender, path = path)
    plot_title <- paste0('Relative Risk split by age for ',gender,'s')
    names(ga_age_pop_est) <- c('gender','Age','population')
  } 
  
  if (is.na(gender)){
    grouped_df <- df %>% group_by(across(Age)) %>% summarise(counts= n())
  } else {
    grouped_df <- df %>% filter(Sex == gender) %>% group_by(across(Age)) %>% summarise(counts= n())
  }
  
  ga_rr_age <- as.data.frame(ga_age_pop_est %>% left_join(grouped_df))
  
  ga_rr_age$population <- ga_rr_age$population - ga_rr_age$counts

  ref_group <- ga_rr_age[ga_rr_age$counts==max(ga_rr_age$counts),]
  
  # filt_uk_rr_age <- uk_rr_age[uk_rr_age$Age!= ref_group$Age,]
  
  rr_rr <- t(mapply(relative_risk_function, ref_group$counts, ref_group$population, 
                    ga_rr_age$counts, ga_rr_age$population))
  
  colnames(rr_rr) <- c('V1','V2','V3')
  
  ga_rr_age <- cbind(ga_rr_age, rr_rr)
  ga_rr_age[ga_rr_age$counts==max(ga_rr_age$counts),c('V1','V2','V3')] <- NA
  ga_rr_age$Age <- factor(ga_rr_age$Age, 
                          levels = c("<=14" ,"15-19" ,"20-24", "25-29", 
                                     "30-34", "35-44" ,"45-54",">55"))
  
  ggplot(ga_rr_age, aes(Age, V2)) + 
    geom_point(size=2) + 
    geom_errorbar(aes(ymin=V1, ymax=V3), width=.2) +
    geom_hline(yintercept=1) + 
    labs(y='Risk Ratio', x='Age', caption = paste0('Reference Group: ',ref_group$Age,'\nPositive Cases: ',ref_group$counts,'\nPopulation: ',ref_group$population)) + 
    ggtitle(plot_title) + 
    theme_bw()
}


# The purpose of this script is to prepare data for model stacking

# file paths, custom fucntions and packages
scripts <- 'C:/Users/tommy/Google Drive/Coursework/1machine_learning/scripts/'
data <- "C:/Users/tommy/Google Drive/Coursework/1machine_learning/data/"
raw <- "C:/Users/tommy/Google Drive/Coursework/1machine_learning/data/raw/"
source(paste0(scripts,'helper_functions.R')) # custom functions
LoadPackages(c("dplyr", "ggplot2", "caret", "tibble", 'Boruta', 'stringr'))

# load data and combine it
load(paste0(raw,"MLProjectData.RData"))
train_data <-  MLProjectData %>% mutate(set = 'training')
submit_data <- test.data %>% mutate(set = "submission")
all_data <- bind_rows(train_data, submit_data)

################################################################################
# create features
 
imp_num_vars <- c('num6', 'num5', 'num58', 'num4',
                  'num23', 'num32', "num18", "num58",
                  "num1", "num3") 

# create square terms
num_vars_df <- all_data %>% select(imp_num_vars)

sq_num_vars <- num_vars_df * num_vars_df
new_names <- names(sq_num_vars) %>% str_replace_all("num", 'sq_num')
names(sq_num_vars) <- new_names

all_data <- bind_cols(all_data, sq_num_vars)


# create interaction terms
int_terms <- do.call(cbind, combn(colnames(num_vars_df), 2, FUN= function(x)
  list(setNames(data.frame(num_vars_df[,x[1]]*num_vars_df[,x[2]]),
                paste(x, collapse="_")) )))

all_data <- bind_cols(all_data, int_terms)

# combine certain category features into a new binary feature
all_data <- all_data %>% 
  mutate(imp_cat = ifelse(cat13 == T | 
                            cat16 == T |
                            cat6 == T |
                            cat19 == T |
                            cat11 == T, 1, 0 ))

# combine all logical category features
cat_log <- all_data %>% 
  select(starts_with("cat")) %>% 
  select(-cat1, -cat2) %>%
  mutate(cat_sum = rowSums(.))

all_data$cat_sum <- cat_log$cat_sum

################################################################################
# select features

# read in boruta feature importance df that I saved from explore.Rmd
boruta_df2 <- read.csv(paste0(data,"boruta_df2.csv"))
boruta_df2$col <- as.character(boruta_df2$col)

# create vector of top 40 features from boruta
mod_vars <- boruta_df2 %>%
  slice(1:40) %>% 
  pull(col)

# add some others that i think is important
mod_vars <- c("set", "target", "cat_sum", mod_vars, "num1", "num4" )

select_data <- all_data %>% select(mod_vars)

################################################################################
# Split data

# sumbission data
submit_data <- select_data %>% filter(set == "submission")

# non submission data
select_data2 <- select_data %>% 
  filter(set != "submission") %>% 
  mutate(id = row_number()) %>%
  select(id, everything())

# split to train and test sets
set.seed(888)
train <- select_data2  %>% sample_frac(.8) # random 80% for training
test <- anti_join(select_data2, train, by = "id") # the remaining for test

# save the data sets 
# saveRDS(submit_data, paste0(data,"submit_data.rds"))
# saveRDS(train, paste0(data,"train.rds"))
# saveRDS(test, paste0(data,"test.rds"))

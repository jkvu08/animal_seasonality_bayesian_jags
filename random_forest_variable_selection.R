##############################
#### Random forest models ####
##############################
library(ggplot2)
library(reshape2)
library(doParallel)
library(tidyverse)
library(dplyr)
library(raster)
library(corrplot)
library(corrgram)
library(lubridate)
library(Hmisc)
library(randomForest)
library(ranger)
library(tree)
library(party)
library(vip)
library(caret)
library(utils)

# detect cores for parallel computing
numcores <- makePSOCKcluster(detectCores())

registerDoParallel(numcores) # register cores

# load data
data <- read.csv('data.csv', header = T, stringsAsFactors = F)

# subsample data to values with at least 24 follow length
datasub <- data[which(data$nfollow5 ==24),] 

# convert categorical variables to factors
datasub$Individual <- factor(datasub$Individual)
datasub$reprostate <- factor(datasub$reprostate, levels = c('MEG','G','LAC','NR')) # seet the levels for the reproductive states
datasub$year <- as.factor(datasub$year) # turn into factors, don't need to sete levels cause it'll automatically go from smalleest to largers
datasub$month <- as.factor(datasub$month) # turn into factors
datasub$biweekly <- as.factor(datasub$biweekly) # turn into factors
datasub$Sex <- as.factor(datasub$Sex) # turn into factors

# summarise data to biweekly resolution (temporal resolution of analysis)
bw_sub <- datasub %>% group_by(year,biweekly, reprostate, fruit_prop, flower_prop,
                                    fruit_shannon, flower_shannon) %>% 
                          summarise(Adults = mean(Adults,na.rm= T), groupchange = mean(groupchange, na.rm  =T),
                                    Infants = mean(Infants, na.rm =T), Juveniles = mean(Juveniles, na.rm = T),
                                    Young = mean(Young, na.rm = T), Rain = sum(Rain, na.rm = T), Daylength = mean(Daylength, na.rm = T),
                                    ERA.Temp2m = mean(ERA.Temp2m, na.rm = T), Feed = mean(Feed, na.rm = T), Rest = mean(Rest, na.rm = T),
                                    Travel = mean(Travel,na.rm =T), FrBw = mean(FrBw, na.rm =T), FlBw = mean(FlBw, na.rm = T),
                                    cp_avg = mean(cumpath, na.rm =T), cp_var = var(cumpath,na.rm = T), tc_avg = mean(treechange,na.rm =T),
                                    tc_var = var(treechange,na.rm =T), gc_avg = mean(groupchange,na.rm =T), gc_var = var(groupchange,na.rm =T), 
                                    adults_var = var(Adults, na.rm = T))
  
# visualize the movement distance
cp1 <- ggplot(datasub, aes(x = biweekly, y = cumpath))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('cumulative path 2 hr (m)') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

cp2 <- ggplot(bw_sub, aes(x = biweekly, y = cp_avg))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('mean cumulative path biweekly (m)') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

cp3 <- ggplot(bw_sub, aes(x = biweekly, y = cp_var))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('var cumulative path biweekly (m)') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

grid.arrange(cp1,cp2,cp3,nrow = 3)

# visualize the frequency at which the animals change trees 
tc1 <- ggplot(datasub, aes(x = biweekly, y = treechange))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('tree change 2hr follow') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

tc2 <- ggplot(bw_sub, aes(x = biweekly, y = tc_avg))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('mean tree change biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

tc3 <- ggplot(bw_sub, aes(x = biweekly, y = tc_var))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('mean tree change biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

tc4 <- ggplot(bw_sub, aes(x = biweekly, y = prop_ut))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('prop unique trees biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

grid.arrange(tc1,tc2,tc3,tc4, nrow = 4)

# visualize the frequency at which the group changes
gc1 <- ggplot(datasub, aes(x = biweekly, y = groupchange))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('group composition change 2 hr') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

gc2 <- ggplot(bw_sub, aes(x = biweekly, y = gc_avg))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('mean group composition change biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

gc3 <- ggplot(bw_sub, aes(x = biweekly, y = gc_var))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('var group composition change biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

grid.arrange(gc1,gc2,gc3,nrow = 3)

# visualize the adult counts 
ac1 <- ggplot(datasub, aes(x = biweekly, y = Adults))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('mean adult count 2 hr') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

ac2 <- ggplot(bw_sub, aes(x = biweekly, y = Adults))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('mean adult count biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

ac3 <- ggplot(bw_sub, aes(x = biweekly, y = adults_var))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('var mean adult count biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

grid.arrange(ac1,ac2,ac3,nrow = 3)

# visualize the proportion of time spent feeding
fp <- ggplot(bw_sub, aes(x = biweekly, y = Feed))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('feeding prop biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

tp <- ggplot(bw_sub, aes(x = biweekly, y = Travel))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('travel prop biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

rp <- ggplot(bw_sub, aes(x = biweekly, y = Rest))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('resting prop biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

frp <- ggplot(bw_sub, aes(x = biweekly, y = FrBw))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('fruit diet prop biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

flp <- ggplot(bw_sub, aes(x = biweekly, y = FlBw))+
  geom_boxplot(fill = 'gray75', outlier.colour = "black", outlier.shape = 1) +
  ylab('flower diet prop biweekly') +
  theme_bw(base_size = 10)+
  theme(strip.text = element_blank())

grid.arrange(fp,tp,rp,frp,flp,nrow = 5)

# Random Forest Modeling

#lets look at a single tree for the combined data
# also get rid of the NA values
start.time<-Sys.time() # get time to monitor how long analysis takes 
cptree <- ranger(formula = cumpath~., # enter the formula  which is the response~. the . signifies that the model should take all the other data  as covariates
                 num.trees=1,importance = 'impurity', # indicate number of trees  & set importance to impurity to later quantify which varaibles are most important for predicting the response
                 data=datasub[which(!is.na(datasub$cumpath)),c(1,3,5:32)], # enter in the data, make sure to subset out the single response and the predictors for that response, get rid of all else when inputting data
                 mtry = floor(29/3),  # set the number of variables the algorithm should try for each node (floor rounds the number, default recommended is the number of covariates/3)
                 respect.unordered.factors = 'order',oob.error = T) #  respect the categorical variables as factors, and report the out of bag error
print(Sys.time() - start.time) # print out how long the code took to run

cptree # print out results

# full forest for the cumulative path
start.time<-Sys.time()
cpforest <- ranger(cumpath~.,num.trees=1000,importance = 'impurity',
                   data=datasub[which(!is.na(datasub$cumpath)),c(1,3,5:32)],
                   mtry = floor(29/3), respect.unordered.factors = 'order', oob.error = TRUE)
print(Sys.time() - start.time) 

print(cpforest)

# plot our the variable importance, tells use how much it contributes to reducing the mean squared error of the model 
# greater importance it what we are looking for!
vip(cpforest,29)  # put in the random forest model object and the top number of variables, I just set it to all of them

# the relationship with resting and traveling proportions are kind of obvious, so let's get rid of them and rerun
cpforest2 <- ranger(cumpath~.,num.trees=1000,importance = 'impurity',
                    data=datasub[which(!is.na(datasub$cumpath)),c(1,3,5:24,27,30:32)],
                    mtry = floor(24/3), respect.unordered.factors = 'order', oob.error = TRUE)
cpforest2

# again map the variable importance
vip(cpforest2,24)+
  theme_bw(base_size = 13)

#lets parse down variables a bit:
cpforest3 <- ranger(cumpath~Feed+tod+ERA.Temp2m+Individual+FrBw+Adults,num.trees=5000,
                    importance = 'impurity',data=datasub[which(!is.na(datasub$cumpath)),],
                    mtry = floor(6/3), respect.unordered.factors = 'order', oob.error = TRUE)

cpforest3 # print results 
vip(cpforest3) # get variable importance

# iterate through all combinations
# extract table for predictors for each resposnse and dump the NA values
cp_df <- datasub[which(!is.na(datasub$cumpath)),]
cpa_df <- bw_sub[which(!is.na(bw_sub$cp_avg)),]

# get predictor column names
cp_col <- c('Young','Adults','ERA.Temp2m', 'reprostate','fruit_prop','fruit_shannon','FlBw','FrBw',
            'flower_prop','flower_shannon','Feed','FBw','month','year','tod','Individual') # 15 options

# generate the equations for the random forest
cp_combo <- cp_col
for (i in 2:15){ # range of numbers indicates the number of variables that should be used
  v <- combn(cp_col, i, simplify = F) # get all combinations
  v <- sapply(1:length(v), FUN =function(x) paste0(v[[x]],collapse = '+')) # concatenate with a plus sign
  dup_ct <- as.integer(grepl('Feed', v, fixed = T)) # see if Feed is a variable and convert to integer
  dup_ct <- dup_ct + as.integer(grepl('Rest', v, fixed = T)) # see if rest is a varible cponver to integer and add
  dup_ct <- dup_ct + as.integer(grepl('Travel', v, fixed = T)) # check travel
  dup_ct <- dup_ct + as.integer(grepl('FBw', v, fixed = T)) # check feeding
  dup_ct <- dup_ct + as.integer(grepl('RBw', v, fixed = T)) # check resting
  dup_ct <- dup_ct + as.integer(grepl('TBw', v, fixed = T)) # check traveling
  index <- which(dup_ct > 1) # get indices that have multiple of these variables
  v <- v[-index] # get rid of the indices since these variables are highly correlated
  cp_combo <- c(cp_combo, v) # add the models to the combinations
  print(i)
}
cp_combo <- sapply(1:length(cp_combo), FUN = function(x) paste0('cumpath~', cp_combo[x]))

# run all variable combinations for random forest
start.time<-Sys.time() 
cp_rf <- foreach(i=1:length(cp_combo), .combine=rbind, .packages=c('ranger','caret','randomForest','party')) %dopar% {
  forest <- ranger(eval(parse(text=cp_combo[i])),
                   data=cp_df, num.trees=1000,importance = 'impurity',
                   respect.unordered.fcptors = TRUE, oob.error = TRUE)
  return(c(cp_combo[i], forest$prediction.error, forest$r.squared))
}

# format results
cp_rf <- as.data.frame(cp_rf) # save results as dataframes
colnames(cp_rf) <- c('mode','MSE','r2') # rename columns
cp_rf$MSE <- as.numeric(cp_rf$MSE) # format datatype
cp_rf$r2 <- as.numeric(cp_rf$r2) # format datatype
cp_rf <- arrange(cp_rf, desc(r2))
print('cumulative path done')
print(Sys.time() - start.time)
# Load packages 
library("readxl")
library("stringr")
library("dplyr")
library("stringdist")
library("psych")
library("sets")
library("cluster")
library("stats")

# Load data from already transformed excel file 
#df <- read_excel('/Users/merelklok/Documents/Econometrie/Business analytics & quantitive marketing/Computer Science/Paper/Data/data_excel.xlsx')
df <- data_excel
df[1,5]

# Clean data and create feature list 
features <- list()
titles <- df %>% select(title)
for(i in 1:nrow(df)) {
  titles[i,] <- tolower(titles[i,])
  titles[i,] <- gsub('Inch|inches|"|-"|-inch| inch|inch', 'inch', titles[i,])
  df[i,5] <- gsub('Inch|inches|"|-"|-inch| inch|inch', 'inch', df[i,5])
  titles[i,] <- gsub('Hz|HZ| hz|-hz|hz|hertz|Hertz', 'hertz', titles[i,])
  df[i,30] <- gsub('Hz|HZ| hz|-hz|hz|hertz|Hertz', 'hertz', df[i,30])
  titles[i,] <- gsub('-|(|)|,|diagonal|diag|diag.', "", titles[i,])
  titles[i,] <- gsub('newegg|neweggcom|amazon|amazoncom|bestbuycom|best|buy|bestbuy|com|tv|the|[()]|class|/',"", titles[i,])
  j <- 1
  next_feature <- FALSE 
  while(next_feature == FALSE){
    feature <- word(titles[i,], start = j)
    feature <- gsub('\\[|\\]', "", feature)
    j <- j + 1
    if(is.na(feature)){
      next_feature <- TRUE
    }
    else if (nchar(feature) < 2 || feature %in% features) {}
    else if(grepl("([0-9].*[a-z])", feature)){
      features <- append(features, feature)
    }
    
  }
}
df[,c('title')] <- titles 

#Create useful lists of aspects needed in analysis 
names <- unlist(df[,1])

#Create binary matrix based on title-features 
bi_matrix <- data.frame(matrix(data=0, nrow=length(features), ncol=length(names)+1, byrow=TRUE))
bi_matrix[,1] <- unlist(features)
colnames(bi_matrix) <- c("Features",names)
for(i in 1:nrow(df)){
  for(k in 1:nrow(bi_matrix)){
    #if(bi_matrix[k,1] %in% unlist(strsplit(as.character(df[i,c("title")])," ")) || bi_matrix[k,1] %in% df[i,5] || bi_matrix[k,1] %in% df[i,30]){
      if(grepl(bi_matrix[k,1], df[i,c("title")])){
        bi_matrix[k,i+1] <- 1
      } 
    }
  }
#}

#Apply min-hashing on binary matrix ##DETERMINE THE AMOUNT OF HASHES --> 50% of amount of rows? 
n_hashing <- 600
min_hash <- data.frame(matrix(data=0, nrow=n_hashing, ncol=length(names), byrow=TRUE))
colnames(min_hash) <- c(names)
for(i in 1:n_hashing){
  sequence <- data.frame(matrix(data=unlist(list(sample(1 : nrow(bi_matrix), size = nrow(bi_matrix), replace = F))), nrow=nrow(bi_matrix),ncol=1,byrow=TRUE))
  for(k in 2:ncol(bi_matrix)){
    value_hashing <- min(sequence[which(bi_matrix[,k] == 1),1])
    min_hash[i,k-1] <- value_hashing
  }
}

#Get bootstrap datasets
bootstraps <- list()
train <- list()
test <- list()
for(i in 1:5){
  train[[i]] <- list()
  test[[i]] <- list()
  bootstrap <- sample(1:length(names), replace = TRUE)
  bootstrap_train <- unique(bootstrap)
  train[[i]] <- bootstrap_train
  test[[i]] <- setdiff(c(1:length(names)),bootstrap_train)
}

#Apply LSH on signature matrix and final duplicate detection###############
#met bootstraps
F1_star_save2 <- list()
F1_save2 <- list()
pair_completeness_save2 <- list()
pair_quality_save2 <- list()
pair_comp_save2 <- list()
pair_dupl_save2 <- list()
r_seq <- seq(1,800,2)
for(i in 1:5){
  testbi <- bi_matrix[,test[[i]]]
  testdf <- df[test[[i]],]
  testnames <- names[test[[i]]]
  testsig <- min_hash[,test[[i]]]
  realdupl <- func_realdupl(testdf,testnames)
  realduplam <- length(realdupl)
  F1_star_save2[[i]] <- NaN*r_seq
  F1_save2[[i]] <- NaN*r_seq
  pair_completeness_save2[[i]] <- NaN*r_seq
  pair_quality_save2[[i]] <- NaN*r_seq
  pair_comp_save2[[i]] <- NaN*r_seq
  pair_dupl_save2[[i]] <- NaN*r_seq
  j <- 1
  for(r in r_seq){
  candidates <- func_candidates(testsig, r)
  F1_star_save2[[i]][j] <- func_F1_star(candidates,realduplam,testdf)
  dissim_matrix <- func_jac_sim(testdf,testsig,candidates,testbi)
  pairqual <- func_pairqual(dissim_matrix,0.45)
  pair_completeness <- pairqual$dupl/realduplam
  pair_quality <- pairqual$dupl/pairqual$comp
  pair_completeness_save2[[i]][j] <- pair_completeness
  pair_quality_save2[[i]][j] <- pair_quality
  pair_comp_save2[[i]][j] <- pairqual$comp
  pair_dupl_save2[[i]][j]<- pairqual$dupl
  F1_save2[[i]][j] <- 2*(pair_completeness*pair_quality)/(pair_completeness+pair_quality)
  j <- j + 1
}
}

#################determine threshold on training data######################
trainbi <- bi_matrix[,train[[5]]]
traindf <- df[test[[5]],]
trainnames <- names[train[[5]]]
trainsig <- min_hash[,train[[5]]]
realdupl <- func_realdupl(traindf,trainnames)
realduplam <- length(realdupl)
F1_star_save <- list()
F1_save <- list()
for(r in seq(10,25,1)){
  candidates <- func_candidates(trainsig, r)
  F1_star <- func_F1_star(candidates,realduplam,traindf)
  F1_star_save[r] <- F1_star
  dissim_matrix <- func_jac_sim(traindf,trainsig,candidates,trainbi)
  #for(e in seq(0.01,0.99,0.02)){
  pairqual <- func_pairqual(dissim_matrix,0.45)
  pair_completeness <- pairqual$dupl/realduplam
  pair_quality <- pairqual$dupl/pairqual$comp
  F1_save <- 2*(pair_completeness*pair_quality)/(pair_completeness+pair_quality)
  #}
}

###########################create graphs #####################################
#average results F1*
F1_star_matrix <- matrix(NA,5,400)
for(i in 1:5){
  F1_star_matrix[i,]<- unlist(F1_star_save2[[i]])
}

F1_star_means <- matrix(NA,1,400)
for(i in 1:ncol(F1_star_matrix)){
  F1_star_means[1,i] <- mean(F1_star_matrix[,i])
}

#average results F1
F1_matrix <- matrix(NA,5,400)
for(i in 1:5){
  F1_matrix[i,]<- unlist(F1_save2[[i]])
}

F1_means <- matrix(NA,1,400)
for(i in 1:ncol(F1_matrix)){
  F1_means[1,i] <- mean(F1_matrix[,i])
}

#average results pair quality
PQ_matrix <- matrix(NA,5,400)
for(i in 1:5){
  PQ_matrix[i,]<- unlist(pair_quality_save2[[i]])
}
PQ_means <- matrix(NA,1,400)
for(i in 1:ncol(PQ_matrix)){
  PQ_means[1,i] <- mean(PQ_matrix[,i])
}

#average results pair completeness
PC_matrix<- matrix(NA,5,400)
for(i in 1:5){
  PC_matrix[i,]<- unlist(pair_completeness_save2[[i]])
}

PC_means <- matrix(NA,1,400)
for(i in 1:ncol(PC_matrix)){
  PC_means[1,i] <- mean(PC_matrix[,i])
}

#compute fraction of comparisons
frac_comparisons <- list()
n <- length(names[test[[1]]])
testsig <- min_hash[,test[[1]]]
totalcomparisons <- (n*(n-1))/2
j <- 1
for(i in r_seq){
  candidates <- func_candidates(testsig, i)
  comparisons <- which(candidates == 1)
  num_comparisons <- length(comparisons)
  frac_comparisons[j] <- num_comparisons/totalcomparisons
  j <- j + 1
}

#create graphs
frac_comp_test <- frac_comparisons[-1]
PC_means_test <- PC_means[,-1]
PQ_means_test <- PQ_means[,-1]
F1_means_test <- F1_means[,-1]
F1_star_means_test <- F1_star_means[,-1]
plot(frac_comp_test,PC_means_test, type = "l", lty = 1, xlab = "Fraction of comparisons", ylab = "Pair Completeness")
plot(frac_comp_test,PQ_means_test, type = "l", lty = 1, xlab = "Fraction of comparisons", ylab = "Pair Quality")
plot(frac_comp_test,F1_means_test, type = "l", lty = 1, xlab = "Fraction of comparisons", ylab = "F1")
plot(frac_comp_test,F1_star_means_test, type = "l", lty = 1, xlab = "Fraction of comparisons", ylab = "F1*")

  


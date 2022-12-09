##################### FUNCTIONS ###############################################
func_candidates <- function(signature, r) {
  n <- nrow(signature)
  b <- floor(n/r)
  r_start = seq(1,n, by= r)
  r_end = seq(r, n, by=r)
  r_end[b] <- n
  candidate_pairs = matrix(0,nrow=ncol(signature),ncol=ncol(signature))
  for(i in 1:b) {
    band <- signature[r_start[i]:r_end[i],]
    buckets <- list()
    for(j in 1:ncol(signature)){
      test_band <- paste(toString(band[,j]), collapse="")
      buckets <- append(buckets, test_band)
    }
    for(k in 1:(length(buckets))){
      pairs <- which(buckets == buckets[[k]])
      for(m in 1:length(pairs)){
        if(pairs[m] != k){
          candidate_pairs[k,pairs[m]] <- 1
        }
      }
    }
  }
  return(candidate_pairs)
}

func_realdupl <- function(df,names){
  realdupl <- list()
  duplicates <- names[duplicated(names)]
  duplicates <- duplicates[!duplicated(duplicates)]
  list <- list()
  for(i in 1:length(duplicates)){
    modelnrs <- which(names==duplicates[i])
    list[[i]] <- modelnrs
  }
  return(list)
}

func_F1_star <- function(candidates,realdupltotal,df) {
  TP <- 0
  pairs <- which(candidates == 1, arr.ind=TRUE)
  for(i in 1:nrow(pairs)){
    if(df[pairs[i,1],4] == df[pairs[i,2],4]){
      TP <- TP + 1
    }
  }
  TP <- TP
  compl <- TP/realdupltotal
  qual <- TP/which(candidates==1)
  F1_star <- (2*qual*compl)/(qual + compl)
  return(F1_star)
}

func_jaccard <- function(x,y) {
  inters <- length(intersect(which(x == 1), which(y == 1)))
  uni <- length(which(x==1)) + length(which(y==1)) - inters
  return(inters/uni)
}

func_jac_sim <- function(df,signature,candidates,bi_matrix){
  brands <- unlist(df[,55])
  shops <- unlist(df[,2])
  
  dissim_matrix = matrix(1000,nrow=ncol(signature),ncol=ncol(signature))
  pairs <- which(candidates == 1, arr.ind=TRUE)
  
  for(i in 1:nrow(pairs)){
    if(!is.na(brands[pairs[i,1]]) && !is.na(tolower(brands[pairs[i,2]])) && tolower(brands[pairs[i,1]])!= tolower(brands[pairs[i,2]])) {
      dissim_matrix[pairs[i,1],pairs[i,2]] <- 1000
    } else if(shops[pairs[i,1]] == shops[pairs[i,2]]){
      dissim_matrix[pairs[i,1],pairs[i,2]] <- 1000
    } else {
      dissim <- 1 - func_jaccard(unlist(bi_matrix_test[,(pairs[i,1])]),unlist(bi_matrix_test[,(pairs[i,2])]))
      dissim_matrix[pairs[i,1],pairs[i,2]] <- dissim 
    }
  } 
  
  return(dissim_matrix)
}

func_pairqual <- function(dissim_matrix,epsilon){
  dupl <- 0
  comp <- 0
  real_pairs <- list()
  dupl_candidates <- func_dupl_candidates(dissim_matrix, epsilon)
  if(length(dupl_candidates)>0){
    for(i in 1:length(dupl_candidates)){ 
      list <- list()
      list <- dupl_candidates[[i]]
      dupl_candidates[[i]] <- NaN*seq(length(list))
      while(length(list)>1){
        for(m in 2:length(list)){
          if(df[list[1],4] == df[list[m],4]){
            dupl <- dupl + 1 
            real_pairs[dupl] <- c(list[1],list[m])
          }
        }
        list <- list[-1]
        comp <- comp + 1
      }
    }
  }
  output <- list("dupl" = dupl, "comp" = comp, "realpairs" = real_pairs)
  return(output)
}

func_dupl_candidates <- function(dissim_matrix,epsilon){
  dupl_candidates <- list()
  for(i in 1:nrow(dissim_matrix)){
    indices <- which(dissim_matrix[,i] <= epsilon)
    indices <- indices[which(indices>i)]
    if(length(indices) > 0){
      indices_list <- list()
      indices_list <- c(i,unlist(indices))
      dupl_candidates <- append(dupl_candidates,list(indices_list))
    }
  }
  return(dupl_candidates)
}

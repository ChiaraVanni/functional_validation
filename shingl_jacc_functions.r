# Jaccard similarity
JaccardSimilarity <- function(x, y) {
  non_zero <- which(x | y)
  set_intersect <- sum( x[non_zero] & y[non_zero] )
  set_union <- length(non_zero)
  return(set_intersect / set_union)
}

# Jaccard function in case Shingling is not needed because we have mono-domain annotations
Jaccard <- function(doc_annot){
  doc_annot <- setNames(doc_annot , c("memb","annot","pres"))
  J <- spread(doc_annot,annot,pres, fill=0)
  rownames(J) <- J$memb
  J <- J[-1]
  pr_DB$set_entry( FUN = JaccardSimilarity, names = c("JaccardSimilarity") )
  # jaccard similarity distance matrix
  ds <- dist(J, method = "JaccardSimilarity" )
  # delete the new entry
  pr_DB$delete_entry("JaccardSimilarity")
  return(ds)
}

# Shingling function, specify k-shingles (in case of multi-domain annotations on the same ORFs)
Shingling <- function(document, k) {
  shingles <- character( length = length(document) - k + 1 )
  for( i in 1:( length(document) - k + 1 ) ) {
    shingles[i] <- paste( document[ i:(i + k - 1) ], collapse = " " )
  }
  return( unique(shingles) )
}

# Jaccard function combined with Shingling (multi-domain annotations)
Jaccard_shingl <- function(doc_annot){
  multi.1 <- doc_annot %>% select(-1,-2)
  multi.1 <- split(multi.1, seq(nrow(multi.1)))
  multi.sh <- lapply(multi.1, function(x){
    Shingling(x, k=2) }) # "shingle" our example document, with k = n-grams, we choose a k=2
  # unique shingles sets across all documents
  doc_dict <- unlist(multi.sh) %>% unique()
  # "characteristic" matrix
  M <- lapply(multi.sh, function(set, dict) {
    as.integer(dict %in% set)
  }, dict = doc_dict) %>% data.frame()
  # set the names for both rows and columns
  setnames(M, paste( "doc", 1:length(multi.sh), sep = "_" ) )
  rownames(M) <- doc_dict
  # how similar is two given document, jaccard similarity
  # create a new entry in the registry
  pr_DB$set_entry( FUN = JaccardSimilarity, names = c("JaccardSimilarity") )
  # jaccard similarity distance matrix
  dm <- dist(t(M), method = "JaccardSimilarity" )
  # delete the new entry
  pr_DB$delete_entry("JaccardSimilarity")
  return(dm)
}

# Function to write results
write_res <- function(tbl,type,clstr,size,med,prop){
  res <- data.frame(rep=clstr$rep[1],
                    jacc_median_raw=med,
                    jacc_median_sc=med*prop,
                    type=type,
                    prop_type=prop,
                    prop_partial=(dim(merge(tbl, clstr,by="memb") %>% filter(partial!="00"))[1])/(dim(clstr)[1]),
                    stringsAsFactors =F)
}

# Comprehensive function for the clusters with multi-domain annotations
# Calls Shingling, Jaccard_shingl and write_res
MultiAnnot <- function(annot_data, annot_type, clstr, size, m_type){
  # split multiple annotations (in many columns), fill empty cells with NAs
  m1 <- annot_data %>% dplyr::select(memb,annot_type) %>% setNames(c("memb","annot")) %>% group_by(memb) %>%
    mutate(annot=paste(annot,collapse = "|")) %>% mutate(pres=1) %>%
    distinct %>% ungroup()
  m2 <- cbind(m1,str_split_fixed(m1$annot, "\\|", m_type))
  m2[,2] <- NULL
  empty_as_na <- function(x){
    if("factor" %in% class(x)) x <- as.character(x) 
    ifelse(as.character(x)!="", x, NA)
  }
  m2 <- m2 %>% mutate_all(funs(empty_as_na))
  multi <- m2 %>% filter(is.na(m2[,4])==F)
  m <- dim(multi)[2] - 2
  n.multi <- dim(multi)[1]
  prop.multi <- n.multi/size
  dm <- Jaccard_shingl(multi)
  if(n.multi==1){
    median <- 1
  } else {
    median <- median(as.matrix(dm)[lower.tri(as.matrix(dm))])
  }
  annotated=n.multi
  prop.annot=prop.multi
  res.m <- write_res(multi,paste("Multi_",annot_type,sep = ""),clstr, size, median, prop.annot)
  # In case there are some member with mono-domain annotations (no shingling)
  singl <- m2 %>% filter(is.na(m2[,4])==T) %>% mutate(pres=1) %>% select(memb,`1`,pres) %>%
    setNames(c("memb","annot","pres"))
  if(dim(singl)[1] > 0){
    n.singl <- dim(singl)[1]
    prop.singl <- n.singl/size
    ds <- Jaccard(singl)
    if(n.singl==1){
      median <- 1
    } else {
      median <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
    }
    prop.annot=prop.singl
    res.s <- write_res(singl,paste("Singl_",annot_type,sep = ""),clstr, size, median, prop.annot)
    res <- rbind(res.s,res.m)
  }else{
    res = res.m
  }
  return(res)
}

library(tidyverse)
library(proxy)
library(stringr)
library(data.table)
library(textreuse)
library(vegan)
library(parallel)

#Subset of annotated clusters

cluster_annot <- fread("data/pro_pang_gc_annot_cvl_all.tsv, stringsAsFactors = F, header = F) %>%
  setNames(c("rep","memb","pf","acc","clan","partial")) %>%
  dplyr::select(rep,memb,partial,pf,clan)

# Pfam terminal (C, N) domains of same proteins
pfam_shared_term <- read.table("data/pfam_comm_names.txt",
                               stringsAsFactors = F, header = F) %>% setNames(c("pfam_name","com_name"))
# Cluster lists
cluster.list <- split(cluster_annot, list(cluster_annot$rep), drop=TRUE)

rm(cluster_annot)
gc()

# Main function
shingl.jacc <- function(cluster, pfam){
  # functional annotation data per cluster
  test_multi <- as.data.frame(cluster) %>% setNames(c("rep","memb","partial","pf","clan"))
  # cluster size
  size <- dim(test_multi)[1]
  # remove not annotated members
  clstrnoNA <- test_multi %>% drop_na()
  # number of annotated members (cluster size)
  prop.annot <- (dim(clstrnoNA)[1])/size
  #select members and annotations (pfam domains and clans)
  m1 <- clstrnoNA %>% dplyr::select(memb,pf,clan)
  ma <- max(str_count(m1$pf, "\\|")) + 1 # max number of multiple domains on the same orfs
  mc <- max(str_count(m1$clan, "\\|")) + 1 # max number of multiple clans on the same orfs

  # Homog_pf: all members annotated to the same Pfam domain
  if(length(unique(m1$pf))==1){
      m1 <- m1 %>% dplyr::select(memb,pf) %>% mutate(pres=1)
      ds <- Jaccard(m1)
      if(dim(m1)[1]==1){  #Only one annotated member
        median <- 1
      }else {
        median <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
      }
      res <- write_res(m1,"Homog_pf", clstrnoNA, size, median, prop.annot)
  }
  # Homog_clan: all members annotated to the same Pfam clan
  else if(length(unique(m1$clan))==1 & any(m1$clan!="no_clan")==T){
      m1 <- m1 %>% dplyr::select(memb,clan) %>% mutate(pres=1)
      ds <- Jaccard(m1)
      median <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
      res <- write_res(m1,"Homog_clan", clstrnoNA, size, median, prop.annot)
      }
  # Not homogeneous annotations
  # First: check for Pfam different terminal domains of the same protein
  else{
        pfam_term <- pfam
        # if the cluster contains any member with multiple annotation
        # split the annotation in multiple rows
        if(ma>1){
          multi_annot <- strsplit(m1$pf, split = "\\|")
          m1 <- data.frame(memb = rep(m1$memb, sapply(multi_annot, length)), pf = unlist(multi_annot), clan=rep(m1$clan, sapply(multi_annot, length)), stringsAsFactors = F)
        }
        # annotation matches with the list of Pfam terminal domains of same proteins
        term <- m1 %>% filter(grepl(paste(pfam_term$com_name,collapse="|"),pf))
        # if there is a correspondance, replace terminal-domain names with the common ones
        if(dim(term)[1]>0){
        m1 <- m1 %>% mutate(pf=plyr::mapvalues(as.vector(.$pf), from = pfam_term$pfam_name, to = pfam_term$com_name))
        }
        # Homog_pf_term: using the common name, do we have homogeneous annotations?
        if(length(unique(m1$pf))==1){
            m1 <- m1 %>% dplyr::select(memb,pf) %>% mutate(pres=1) %>% distinct
            ds <- Jaccard(m1)
            median <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
            res <- write_res(m1,"Homog_pf_term", clstrnoNA, size, median, prop.annot)
        }
        # Still not homogeneous annotations..
        else{
          # Only mono-domain annotations
            if(ma==1 & mc==1){
              m1.a <- m1 %>% dplyr::select(memb,pf) %>% mutate(pres=1) %>% distinct
              ds <- Jaccard(m1.a)
              median.a <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
              m1.c <- m1 %>% dplyr::select(memb,clan) %>% mutate(pres=1) %>% distinct
              ds <- Jaccard(m1.c)
              median.c <- median(as.matrix(ds)[lower.tri(as.matrix(ds))])
                if(median.a >= median.c){
                  median <- median.a
                  res <- write_res(m1,"Mono_pf", clstrnoNA, size, median, prop.annot)
                }else{
                  median <- median.c
                  res <- write_res(m1,"Mono_clan", clstrnoNA, size, median, prop.annot)
                }
              }else{ # The cluster contains also multi-domain annotations
                # using clans
                res.c <- MultiAnnot(m1,"clan", clstrnoNA, size, mc)
                #using pfam names
                res.a <- MultiAnnot(m1,"pf", clstrnoNA, size, ma)
                # Choose the best results between pfam name and clan
                if(any(!grepl("no_clan",m1$clan))==T && any(as.numeric(res.a$jacc_median_sc)< as.numeric(res.c$jacc_median_sc))==T){
                  res <- res.c
                }else{
                  res <- res.a
                }
              }
        }
  }
  return(res)
}

source("shingl_jacc_functions.r")
res.list <- mclapply(cluster.list, shingl.jacc, pfam=pfam_shared_term, mc.cores = getOption("mc.cores",28))
results <- plyr::ldply(res.list, data.frame)
res.parsed.1 <- results %>% select(rep,jacc_median_raw,jacc_median_sc,type,prop_type,prop_partial) %>% group_by(rep) %>%
  mutate(count = n_distinct(type)) %>% filter(count==1)
res.parsed.2 <- results %>% select(rep,jacc_median_raw,jacc_median_sc,type,prop_type,prop_partial) %>% group_by(rep) %>%
  mutate(count = n_distinct(type)) %>% filter(count==2) %>%
  filter(prop_type==max(prop_type)) %>%
  filter(prop_partial==min(prop_partial)) %>%
  filter(jacc_median_sc==max(jacc_median_sc)) %>%
  filter(jacc_median_raw==max(jacc_median_raw)) %>% group_by(rep) %>% slice(1)

shingl_jacc.res <- rbind(res.parsed.1,res.parsed.2) %>% select(-count)
write.table(shingl_jacc.res, "results/pro_pan_gc_func_eval.tsv",col.names=T,row.names=F,sep="\t",quote=F)

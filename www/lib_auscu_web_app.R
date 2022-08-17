## name: lib_auscu_web_app.R
##    
## author: Federico M. Stefanini
## email: federico.stefanini@unifi.it
## date: 31 July 2020
## update: 29 August 2020
## update: 4 Sept. 2020 -  rel:  1.9
## update 10 Sept. 2020
## update 11 Sept. 2020 -  rel:  1.11
## update 26 Oct. 2020 -  rel:  1.12
## update 15 August 2022 -  rel:  1.13



#library(BiocManager)
#options(repos = BiocManager::repositories())

library(utils)

options(repos=c(Bioconductor = "https://bioconductor.org/packages/3.11/bioc", 
                BioCsoft = "https://bioconductor.org/packages/3.11/bioc", 
                BioCann = "https://bioconductor.org/packages/3.11/data/annotation", 
                BioCexp = "https://bioconductor.org/packages/3.11/data/experiment", 
                BioCworkflows = "https://bioconductor.org/packages/3.11/workflows", 
                CRAN = "https://cloud.r-project.org"))

library(shiny)
library(shinyjs)
library(stringr)
library(shinyalert)
library(shinybusy)
library(magrittr)
library(Rcpp)
library(igraph)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(graph)
library(Rgraphviz)
library(gRbase)
library(RBGL)
library(gRain)
library(magrittr)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

mapp_aggrega_ILD <- function(
                       dbLS  
){
  res<-list(
    dataNow = dbLS$dataNow,
    myNotes =dbLS$myNotes)
  # component
  tmp_out <- aggregate_ILD(dbLS$tavola1,F)
  res$tavola1 <- tmp_out[[1]]
  res$tavola2 <- tmp_out[[2]]
  # component 5
  tmp1 <- as_tibble(tmp_out[[1]])
  tmp11 <- tmp1 %>% gather( 'CHF', 'Pneumonia', 'ILDs', 'LungScars', key="Disease",value="Prob")
  tmp111 <- tmp11 %>% dplyr::filter(ID == 'Init_prob' | ID =='Final_prob' )
  res$figura1 <-  ggplot(tmp111, aes(fill=ID, y=Prob, x=Disease)) + 
                          geom_bar(position="dodge", stat="identity");
  # component 6    
  tmp2 <- as_tibble(tmp_out[[2]])
  tmp22 <- tmp2[3,] %>% gather( 'CHF', 'Pneumonia', 'ILDs', 'LungScars', key="Disease",value="BFval")
  res$figura2 <-  ggplot(tmp22, aes(  y=BFval, x=Disease)) + 
                          geom_bar(stat="identity");
  #
  res$current_evi <- dbLS$current_evi # 7
  res$spe_sen_LS  <- dbLS$spe_sen_LS  # 8
  # componente 9
  res$effe_size_stat <-  list(res_all_evid=list(tavolaPro=NA,tavolaBF=NA),
                              res_loo=list(),
                              res_switchin=list(),
                              res_missing=list());
  ttmmpp <- aggregate_ILD(dbLS$effe_size_stat$res_all_evid$tavolaPro,T)
  res$effe_size_stat$res_all_evid$tavolaPro <- ttmmpp[[1]]
  res$effe_size_stat$res_all_evid$tavolaBF <- ttmmpp[[2]]
  # cicli
  for(aux in seq_along(dbLS$effe_size_stat$res_loo)){
    
    ta_pro <- dbLS$effe_size_stat$res_loo[[aux]]$tavolaPro
    ttmmpp <- aggregate_ILD(ta_pro,T)
    res$effe_size_stat$res_loo[[aux]] <- list(deletedEvidence=NA)
    res$effe_size_stat$res_loo[[aux]]$deletedEvidence <- dbLS$effe_size_stat$res_loo[[aux]]$deletedEvidence  
    res$effe_size_stat$res_loo[[aux]]$tavolaPro  <- ttmmpp[[1]]
    res$effe_size_stat$res_loo[[aux]]$tavolaBF   <- ttmmpp[[2]]
    }
   
  for(aux in seq_along(dbLS$effe_size_stat$res_switchin)){
    ta_pro <- dbLS$effe_size_stat$res_switchin[[aux]]$tavolaPro
    ttmmpp <- aggregate_ILD(ta_pro,T)
    res$effe_size_stat$res_switchin[[aux]] <-  list(node=NA,state=NA)
    res$effe_size_stat$res_switchin[[aux]]$node <- dbLS$effe_size_stat$res_switchin[[aux]]$node
    res$effe_size_stat$res_switchin[[aux]]$state <- dbLS$effe_size_stat$res_switchin[[aux]]$state
    res$effe_size_stat$res_switchin[[aux]]$tavolaPro  <- ttmmpp[[1]]
    res$effe_size_stat$res_switchin[[aux]]$tavolaBF   <- ttmmpp[[2]]
  }
  #
   
  for(aux in seq_along(dbLS$effe_size_stat$res_missing)){
    ta_pro <- dbLS$effe_size_stat$res_missing[[aux]]$tavolaPro
    ttmmpp <- aggregate_ILD(ta_pro,T)
    res$effe_size_stat$res_missing[[aux]] <-  list()
    res$effe_size_stat$res_missing[[aux]]$node <- dbLS$effe_size_stat$res_missing[[aux]]$node
    res$effe_size_stat$res_missing[[aux]]$state <- dbLS$effe_size_stat$res_missing[[aux]]$state 
    res$effe_size_stat$res_missing[[aux]]$tavolaPro <- ttmmpp[[1]]
    res$effe_size_stat$res_missing[[aux]]$tavolaBF  <- ttmmpp[[2]]
  }
  #  
  res$release_cur <- dbLS$release_cur  # 
  return(res);#  
}


function(){
  debug(mapp_aggrega_ILD)
  
  mapp_aggrega_ILD(params)
  undebug(mapp_aggrega_ILD)
  
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

aggregate_ILD <- function(tavoPro,
                          labels_out= TRUE){
  
  
   if( !("tbl" %in% class(tavoPro))){ 
     tmpdat   <- as_tibble(tavoPro)
     tmpdat <- mutate(tmpdat,ID=dimnames(tavoPro)[[1]],.before=1)
     tavoPro <- tmpdat
   } 
    nuova_tabPro <- cbind((tavoPro[1:3]) %>% mutate(ILDs = tavoPro[[4]]+tavoPro[[5]]+tavoPro[[6]]),
                         tavoPro[7]
                          )
   inipro <- as.numeric(nuova_tabPro[1,-1])
   finpro <-  as.numeric(nuova_tabPro[2,-1])
    
   nuova_tabPro[3,2:5] <-finpro/inipro
   nuova_tabPro[4,2:5] <-inipro/finpro
   
   nuova_tabBF <- nuova_tabPro[1:3,]
   nuova_tabBF[[1]] <- c("iniOdds","finOdds","BFval")
   iniOdds <-inipro/(1-inipro)  
   finOdds <- finpro /(1-finpro )
     
   nuova_tabBF[1,2:5] <-  as.numeric(iniOdds)
   nuova_tabBF[2,2:5] <-  as.numeric(finOdds)
   nuova_tabBF[3,2:5] <-  as.numeric(finOdds/iniOdds)
    
   
   if(labels_out ==TRUE){
     nuova_tabPro <-  as.data.frame(nuova_tabPro[,-1],
                                    row.names = nuova_tabPro[,1]);
     nuova_tabBF <-  as.data.frame(nuova_tabBF[,-1],
                                    row.names = nuova_tabBF[,1]);
     
   }
   
list(nuova_tabPro,nuova_tabBF) 
}

                             
function(){
  aggregate_ILD(params$tavola1, T)
  aggregate_ILD(params$tavola1, F)
  params$tavola1
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

diverging_bchart_sos <- function(effe_size_stat,
                                 type_sco=c("tavolaPro","tavolaBF")[1]){
  
  
  nomi_diseases <- dimnames(effe_size_stat[["res_all_evid"]][["tavolaPro"]])[[2]];
  if(type_sco == "tavolaPro"){
    ve_score <-  effe_size_stat[["res_all_evid"]][["tavolaPro"]]["Final_prob",]
    etichettaY <- "Probability ratio"
    etichetta_main <- " posterior probability (ratio) "
  }else if(type_sco == "tavolaBF"){
    ve_score <-  effe_size_stat[["res_all_evid"]][["tavolaBF"]]["BFval",]
    etichettaY <- "Bayes Factor"
    etichetta_main <-" Bayes Factor (ratio)"
  }else{
    stop("Wrong selection in diverging_bchart_sos.")
  }
  posizBest <-  which(ve_score == max(ve_score))
   
  nameBest <-  nomi_diseases[posizBest]
  
   
  if(type_sco == "tavolaPro"){
    ve_sco_sos <- lapply(effe_size_stat[["res_switchin"]],function(vx){
      list(score= as.numeric(vx[[type_sco]]["Final_prob",][posizBest]),
           disease = nameBest,
           typeScore = type_sco,
           node_sos = vx[["node"]],
           state= vx[["state"]]
      )})
  }else  if(type_sco == "tavolaBF"){
    ve_sco_sos <- lapply(effe_size_stat[["res_switchin"]],function(vx){
      list(score= as.numeric(vx[[type_sco]]["BFval",][posizBest]),
           disease = nameBest,
           typeScore = type_sco,
           node_sos = vx[["node"]],
           state= vx[["state"]]
      )})
  }# end if else
  
  
  myTB <- tibble(
    typeScore = sapply(ve_sco_sos,function(vx){vx$typeScore}),
    node_sos = sapply(ve_sco_sos,function(vx){vx$node_sos}),
    state= sapply(ve_sco_sos,function(vx){vx$state}),
    scoreRatio  = sapply(ve_sco_sos,function(vx){vx$score}) /as.numeric(ve_score[posizBest]), #  
    top_disease = sapply(ve_sco_sos,function(vx){vx$disease}), 
    node_state=NA
  )
  myTB <- mutate(myTB, node_state = paste0(.data$node_sos,"=",.data$state))
  
  myGR <- ggplot(data = myTB,
                 aes(x = reorder(node_state, scoreRatio), y = scoreRatio-1,
                     fill = scoreRatio> 1.0))+
    geom_bar(stat = "identity")+
    coord_flip()+
    labs(x = "Switched  evidence", y = etichettaY,
         title = paste0("Changes of ",etichetta_main," after switching state to one variable"),
         subtitles = paste0("Top disease after collecting evidence: ",nameBest ))+
    scale_y_continuous(labels = function(y) y + 1) +
    theme_minimal()+
    guides(fill = FALSE)
  
  return(myGR)
} 




#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

diverging_bchart_loo <- function(effe_size_stat,
                                  type_sco=c("tavolaPro","tavolaBF")[1]){
  
  
  nomi_diseases <- dimnames(effe_size_stat[["res_all_evid"]][["tavolaPro"]])[[2]]
  if(type_sco == "tavolaPro"){
    ve_score <-  effe_size_stat[["res_all_evid"]][["tavolaPro"]]["Final_prob",]
    etichettaY <- "Probability ratio"
    etichetta_main <- " posterior probability (ratio) "
  }else if(type_sco == "tavolaBF"){
    ve_score <-  effe_size_stat[["res_all_evid"]][["tavolaBF"]]["BFval",]
    etichettaY <- "Bayes Factor"
    etichetta_main <-" Bayes Factor (ratio)"
  }else{
    stop("Wrong selection in diverging_bchart_loo.")
  }
  posizBest <-  which(ve_score == max(ve_score))
  nameBest <-  nomi_diseases[posizBest]
   
  if(type_sco == "tavolaPro"){
    ve_sco_loo <- lapply(effe_size_stat[["res_loo"]],function(vx){
      list(score= as.numeric(vx[[type_sco]]["Final_prob",][posizBest]),
           disease = nameBest,
           typeScore = type_sco,
           deletedEvid = vx[["deletedEvidence"]]
      )})
  }else  if(type_sco == "tavolaBF"){
    ve_sco_loo <- lapply(effe_size_stat[["res_loo"]],function(vx){
      list(score= as.numeric(vx[[type_sco]]["BFval",][posizBest]),
           disease = nameBest,
           typeScore = type_sco,
           deletedEvid = vx[["deletedEvidence"]]
      )})
  }#  
  
  myTB <- tibble(
    typeScore = sapply(ve_sco_loo,function(vx){vx$typeScore}),
    deletedEvid = sapply(ve_sco_loo,function(vx){vx$deletedEvid}),
    scoreRatio  = sapply(ve_sco_loo,function(vx){vx$score}) /as.numeric(ve_score[posizBest]), 
     
    top_disease = sapply(ve_sco_loo,function(vx){vx$disease}) 
     
  )
  
  myGR <- ggplot(data = myTB,
                 aes(x = reorder(deletedEvid, scoreRatio), y = scoreRatio-1,
                     fill = scoreRatio> 1.0))+
    geom_bar(stat = "identity")+
    coord_flip()+
    labs(x = "Deleted evidence", y = etichettaY,
         title = paste0("Changes of ",etichetta_main," after deleting one variable"),
         subtitles = paste0("Top disease after collecting evidence: ",nameBest ))+
    scale_y_continuous(labels = function(y) y + 1) +
    theme_minimal()+
    guides(fill = FALSE)
  
  return(myGR)
} 








#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

diverging_bchart_miss <- function(effe_size_stat,
                                  type_sco=c("tavolaPro","tavolaBF")[1]){
  
   #  
   nomi_diseases <- dimnames(effe_size_stat[["res_all_evid"]][["tavolaPro"]])[[2]]
   if(type_sco == "tavolaPro"){
     ve_score <-  effe_size_stat[["res_all_evid"]][[type_sco]]["Final_prob",]
     etichettaY <- "Probability ratio"
     etichetta_main <- " posterior probability (ratio) "
   }else{
     ve_score <-  effe_size_stat[["res_all_evid"]][[type_sco]]["BFval",]
     etichettaY <- "Bayes Factor"
     etichetta_main <-" Bayes Factor "
   }
   posizBest <-  which(ve_score == max(ve_score))
   # 
   nameBest <-  nomi_diseases[posizBest]
   
   
   if(type_sco == "tavolaPro"){
     ve_sco_miss <- lapply(effe_size_stat[["res_missing"]],function(vx){
         list(score= as.numeric(vx[[type_sco]]["ratioPostPre",][posizBest]),
            variable = names(vx[[type_sco]]["ratioPostPre",][posizBest]),
            state = vx[["state"]],
            node = vx[["node"]]
       )})
     }else{
       ve_sco_miss <- lapply(effe_size_stat[["res_missing"]],function(vx){
         list(score= as.numeric(vx[[type_sco]]["BFval",][posizBest]),
              variable = names(vx[[type_sco]]["BFval",][posizBest]),
              state = vx[["state"]],
              node = vx[["node"]]
         )})
       
    }
       

  myTB <- tibble(
    miss_imputed =sapply(ve_sco_miss,function(vx){vx$node}),
    miss_state = sapply(ve_sco_miss,function(vx){vx$state}),
    score      = sapply(ve_sco_miss,function(vx){vx$score}),
    top_disease=sapply(ve_sco_miss,function(vx){vx$variable}),
    impu_evid=NA
  )
  myTB <- mutate(myTB, impu_evid = paste0(.data$miss_imputed,"=",.data$miss_state))
  
  myGR <- ggplot(data = myTB,
                 aes(x = reorder(impu_evid, score), y = score-1,
                     fill = score> 1.0))+
    geom_bar(stat = "identity")+
    coord_flip()+
    labs(x = "Imputed evidence", y = etichettaY,
         title = paste0("Changes of ",etichetta_main," after imputing one unobserved variable"),
         subtitles = paste0("Top disease after collecting evidence: ",nameBest ))+
    scale_y_continuous(labels = function(y) y + 1) +
    theme_minimal()+
    guides(fill = FALSE)
  
  return(myGR)
} 







#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
reset_evidences <- function(input,session){
  updateSelectInput(session,"Age_above60", selected = "not-informed")
  updateSelectInput(session,"Sex", selected = "not-informed")
  updateSelectInput(session,"Exposures", selected = "not-informed")
  updateSelectInput(session,"ActiEx_Smoke", selected = "not-informed")
  updateSelectInput(session,"GERD", selected = "not-informed")
  updateSelectInput(session,"Histo_autoim_rheuma", selected = "not-informed")
  updateSelectInput(session,"Histo_HF_cardi_dis", selected = "not-informed")
  updateSelectInput(session,"Pneumotox_medicat", selected = "not-informed")
  updateSelectInput(session,"Asthma", selected = "not-informed")
  updateSelectInput(session,"Histo_pulmo_dis", selected = "not-informed")
  updateSelectInput(session,"Familia_histo_FILD", selected = "not-informed")
  updateSelectInput(session,"Orthopnea", selected = "not-informed")
  updateSelectInput(session,"Per_malleo_edema", selected = "not-informed")
  updateSelectInput(session,"Compla_limit_exerc_toler", selected = "not-informed")
  updateSelectInput(session,"Digit_clubb", selected = "not-informed")
  updateSelectInput(session,"Onset_sympto", selected = "not-informed")
  updateSelectInput(session,"Progr_chron_dysp_exer", selected = "not-informed")
  updateSelectInput(session,"Cough", selected = "not-informed")
  updateSelectInput(session,"Fever", selected = "not-informed")
  updateSelectInput(session,"Chest_tight", selected = "not-informed")
  updateSelectInput(session,"Dysphagia", selected = "not-informed")
  updateSelectInput(session,"Lab_test_Whi_Blo_Cel", selected = "not-informed")
  
  updateSelectInput(session,"Raynaud_phenom", selected = "not-informed")
  updateSelectInput(session,"Decreas_DLCO70", selected = "not-informed")
  updateSelectInput(session,"Restr_lung_dis_TLC80", selected = "not-informed")
  updateSelectInput(session,"Oxigen_sat_rest90", selected = "not-informed")
  updateSelectInput(session,"Lips_cianos", selected = "not-informed")
  updateSelectInput(session,"AuscVC_Basal", selected = "not-informed")
  
  updateSelectInput(session,"AuscVC_Lateral", selected = "not-informed")
  updateSelectInput(session,"AuscVC_Higher", selected = "not-informed")
  updateSelectInput(session,"Other_abnor_lun_soun", selected = "not-informed")
  updateSelectInput(session,"Third_heart_sound", selected = "not-informed")
  updateSelectInput(session,"AuscVC_Sym_same_lev", selected = "not-informed")
  updateSelectInput(session,"EE_AuscVC_Basal", selected = "not-informed")
  updateSelectInput(session,"EE_AuscVC_Lateral", selected = "not-informed")
  updateSelectInput(session,"EE_AuscVC_Higher", selected = "not-informed")
  updateSelectInput(session,"EE_AuscVC_Sym_same_lev", selected = "not-informed")
  updateTextInput(session,
                  inputId="propagaStatus",
                  value="NO")
}



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

get_sample_space <- function(){
  
 sam_spa <- list()
  
 sam_spa[["EE_AuscVC_Basal"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["EE_AuscVC_Sym_same_lev"]] <- list(omega = c("Bilateral","Unilateral"), posterior = list())
 sam_spa[["EE_AuscVC_Lateral"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["EE_AuscVC_Higher"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Age_above60"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Sex"]] <- list(omega = c('Male','Female'), posterior = list())
 sam_spa[["Exposures"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["ActiEx_Smoke"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["GERD"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Histo_autoim_rheuma"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Histo_HF_cardi_dis"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Pneumotox_medicat"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Asthma"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Histo_pulmo_dis"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Familia_histo_FILD"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Orthopnea"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Per_malleo_edema"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Compla_limit_exerc_toler"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Digit_clubb"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Onset_sympto"]] <- list(omega = c('FAST','SLOW'), posterior = list())
 sam_spa[["Progr_chron_dysp_exer"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Cough"]] <- list(omega = c('NO','DRY','PRODUCTIVE'), posterior = list())
 sam_spa[["Fever"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Chest_tight"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Dysphagia"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Lab_test_Whi_Blo_Cel"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Raynaud_phenom"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Decreas_DLCO70"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Restr_lung_dis_TLC80"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Oxigen_sat_rest90"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["Lips_cianos"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["AuscVC_Basal"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["AuscVC_Sym_same_lev"]] <- list(omega = c("Bilateral","Unilateral"), posterior = list())
 sam_spa[["AuscVC_Lateral"]] <- list(omega = c('YES','NO'), posterior = list())
 sam_spa[["AuscVC_Higher"]] <- list(omega = c('YES','NO'), posterior = list()) 
 sam_spa[["Other_abnor_lun_soun"]] <- list(omega = c('Wheezing','Squacks','Pleur_frictio','None'), posterior = list())
 sam_spa[["Third_heart_sound"]] <- list(omega = c('YES','NO'), posterior = list())

sam_spa
}





#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

two_tables_summary <- function(
               iniProb,  #  
               finalProb #  
               ){
  
  iniOdds <- iniProb/(1-iniProb)  
  finOdds <-   finalProb / (1-finalProb)
  BFval <- finOdds/iniOdds
  probRatio <- finalProb/iniProb
  invProbRatio <-  1/probRatio
  # tables
  tavolaPro <- rbind(iniProb,
                     finalProb,
                     probRatio, 
                     invProbRatio)
  dimnames(tavolaPro)[[1]] <- c("Init_prob","Final_prob","ratioPostPre","ratioPrePost")
  tavolaBF <- rbind(iniOdds,finOdds,BFval)
  return(list(
    tavolaPro = tavolaPro,
    tavolaBF  = tavolaBF  
  ))
}


tmp <- function(){
  argoIni <- c( d1=0.5 ,d2=0.3, d3=0.2)
  argoFin <-  c(d1=0.1 ,d2=0.2, d3=0.7)
  two_tables_summary(argoIni,argoFin)
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
effe_size_table_IE <-function(
                        p_disease_0,         
                        p_disease_full,      
                        p_disease_full_loo,  
                        p_switch_one_in,     
                        post_miss_one_in     
                        ){                   
  # step 1
  step1 <- two_tables_summary(p_disease_0$Disease,
                              p_disease_full$Disease)
  # step2  
  res_loo <- list()
  nomi_manifest <- names(p_disease_full_loo)
  for(aux in seq_along(p_disease_full_loo)){
    tmpTab <- two_tables_summary( p_disease_full$Disease,
                                  p_disease_full_loo[[aux]]$Disease)
    res_loo[[aux]] <- list(deletedEvidence = nomi_manifest[aux],
                         tavolaPro = tmpTab$tavolaPro,
                         tavolaBF  = tmpTab$tavolaBF)
  }
  # step3  
  res_sw1in <- list()
  punta <- 1
  nomi_nodes <- names(p_switch_one_in)
  for(auxNode in seq_along(p_switch_one_in)){
       for(auxState in seq_along(p_switch_one_in[[auxNode]]$omega)){
         tmpTab <- two_tables_summary(p_disease_full$Disease,
                                      p_switch_one_in[[auxNode]]$posterior[[auxState]]$Disease)
         res_sw1in[[punta]] <- list(node= nomi_nodes[auxNode],
              state = p_switch_one_in[[auxNode]]$omega[auxState],
              tavolaPro = tmpTab$tavolaPro,
              tavolaBF  = tmpTab$tavolaBF
              )
         punta <- punta + 1
       }
    
  }
  #step  
  res_mis1in <- list()
  punta <- 1
  nomi_nodes <- names(post_miss_one_in)
  for(auxNode in seq_along(post_miss_one_in)){
    for(auxState in seq_along(post_miss_one_in[[auxNode]]$omega)){
      tmpTab <- two_tables_summary(p_disease_full$Disease,
                                   post_miss_one_in[[auxNode]]$posterior[[auxState]]$Disease)
      res_mis1in[[punta]] <- list(node= nomi_nodes[auxNode],
                                 state = post_miss_one_in[[auxNode]]$omega[auxState],
                                 tavolaPro = tmpTab$tavolaPro,
                                 tavolaBF  = tmpTab$tavolaBF
      )
      punta <- punta + 1
    }
    
  }
  
  return(
    list(
      #  
      res_all_evid = list(tavolaPro = step1$tavolaPro,
                     tavolaBF = step1$tavolaBF),
      #  
      res_loo = res_loo,
      #  
      res_switchin = res_sw1in,
      #  
      res_missing = res_mis1in
   ))
}



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
effe_size <-function(myFindings,
                     input,session ){
    #  
    show_modal_spinner( spin = "fulfilling-bouncing-circle", 
                        color = "#0e4aa2",
                        text = "Propagation in progress: please wait ...") 
     current_evi <- myFindings[!is.na(myFindings)]
     names_context <- c("EE_AuscVC_Sym_same_lev","EE_AuscVC_Basal")
    names_restri <-  setdiff(names(current_evi), c("EE_AuscVC_Sym_same_lev","EE_AuscVC_Basal"))
     bn_3_3c <-  compile_my_bn (input,session)
     compNet <- retractEvidence(bn_3_3c)
    
    if(is.na(myFindings$EE_AuscVC_Sym_same_lev) || is.na(myFindings$EE_AuscVC_Basal )){
        
       return(NULL);
    }
    #  
    compNet <- setEvidence(compNet, evidence= myFindings[names_context])
    post_disease_0 <- querygrain( compNet, nodes=c("Disease"))  
     
    compNet <- setEvidence(compNet, evidence= myFindings[names_restri])
    post_disease_full <- querygrain( compNet, nodes=c("Disease"))  
    
    post_disease_full_loo <- list() 
    for(aux in seq_along(names_restri)){
       nome_cor <- names_restri[aux]
       compNet <- retractEvidence(compNet, nodes =  nome_cor, propagate = TRUE)
        
       post_disease_full_loo[[nome_cor]] <- querygrain( compNet, nodes=c("Disease"))
       compNet <- setEvidence(compNet, evidence= current_evi[nome_cor], propagate = TRUE)
    }
    
    post_switch_one_in <- get_sample_space()
     
    post_switch_one_in <- post_switch_one_in[names_restri]
    for(aux in names_restri){ 
      
      post_switch_one_in[[aux]]$omega <- setdiff(post_switch_one_in[[aux]]$omega, 
                                                myFindings[[aux]])
      #  
      for(auxalte in post_switch_one_in[[aux]]$omega){
        compNet <- retractEvidence(compNet, nodes = aux, propagate = TRUE)
        fakeEvidence <-  list();fakeEvidence[[aux]] <- auxalte;
        compNet <- setEvidence(compNet, evidence= fakeEvidence, propagate = TRUE)
        post_switch_one_in[[aux]]$posterior[[auxalte]] <- querygrain( compNet, nodes=c("Disease"))
      }
      compNet <- retractEvidence(compNet, nodes = aux, propagate = TRUE)
     
      compNet <- setEvidence(compNet, evidence= current_evi[aux], propagate = TRUE)
    }# 
    #
     
    names_missing <- setdiff(names(myFindings),names(current_evi))
    #  
    names_missing <- setdiff(names_missing, c("AuscVC_Basal","AuscVC_Lateral",
                               "AuscVC_Higher","AuscVC_Sym_same_lev"));
    post_miss_one_in <- get_sample_space()
    post_miss_one_in <- post_miss_one_in[names_missing]
    for(aux in names_missing){#
      for(auxalte in post_miss_one_in[[aux]]$omega){
         
        fakeEvidence <-  list();fakeEvidence[[aux]] <- auxalte;
        compNet <- setEvidence(compNet, evidence= fakeEvidence, propagate = TRUE)
        post_miss_one_in[[aux]]$posterior[[auxalte]] <- querygrain( compNet, nodes=c("Disease"))
        compNet <- retractEvidence(compNet, nodes = aux, propagate = TRUE)
      }
       
    }# 
    
  
    efSiTableIE <- effe_size_table_IE( post_disease_0,
                                     post_disease_full,
                                     post_disease_full_loo,
                                     post_switch_one_in,
                                     post_miss_one_in)
     
    remove_modal_spinner()
    return(efSiTableIE )
}#  





 
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
propaga_ev <-function(
                    myFindings,
                    session ){
  #  
  show_modal_spinner( spin = "fulfilling-bouncing-circle", 
                      color = "#0e4aa2",
                      text = "Propagation in progress: please wait ..."
  ) 
   net_3_3.c <-   get_userData("compiled_net",session = session)
  #
 
  findiLS <- myFindings[!is.na(myFindings)]
  res.Disease.evid <-  NA
  compiNet.r <- retractEvidence(net_3_3.c)
  res.Disease.0 <- querygrain( compiNet.r,
                               nodes=c("Disease"),
                               type="marginal")  
  
  if(length(findiLS) >0){
    net.fin.p <- setEvidence(compiNet.r, evidence= findiLS)
    res.Disease.evid <- querygrain( net.fin.p, nodes=c("Disease"))  
  }else{
    res.Disease.evid <- querygrain( compiNet.r, nodes=c("Disease"))
  }   
  
  iniOdds <- res.Disease.0$Disease/(1-res.Disease.0$Disease)  
  finOdds <-   res.Disease.evid$Disease / (1-res.Disease.evid$Disease)
  BFval <- finOdds/iniOdds
  probRatio <- res.Disease.evid$Disease/res.Disease.0$Disease
  invProbRatio <-  1/probRatio
  # tables
  tavolaPro <- rbind(res.Disease.0$Disease,res.Disease.evid$Disease,
                     probRatio, invProbRatio)
  dimnames(tavolaPro)[[1]] <- c("Init_prob","Final_prob",
                                "ratioPostPre","ratioPrePost")
  tavolaBF <- rbind(iniOdds,finOdds,BFval)
  remove_modal_spinner()
  return(
    list(initial= res.Disease.0,
         conditioned = res.Disease.evid,
         ratioPostPre = probRatio,
         ratioPrePost = invProbRatio,
         iniOdds = iniOdds,
         finalOdds = finOdds,
         valBF = BFval,
         tavolaPro= tavolaPro,
         tavolaBF =tavolaBF )
  )  
}# end of function


 

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
set_userData <- function( name,
                          value,
                          session){
  assign(name,value,envir=session$userData )
}






#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

get_userData <- function( name,
                          session){
   
  if(exists(name,envir=session$userData)) {
    return(get(name,envir=session$userData))
  }else{
    return(NA)
  }
}




#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#  

read_evidences <-  function(input, session){
  res <- list(
    Age_above60= c(NA),
    Sex= c(NA),
    Exposures= c(NA),
    ActiEx_Smoke= c(NA),
    GERD= c(NA),
    Histo_autoim_rheuma= c(NA),
    Histo_HF_cardi_dis= c(NA),
    Pneumotox_medicat= c(NA),
    Asthma= c(NA),
    Histo_pulmo_dis= c(NA),
    Familia_histo_FILD= c(NA),
    #
    Orthopnea= c(NA),
    Per_malleo_edema= c(NA),
    Compla_limit_exerc_toler= c(NA),
    Digit_clubb= c(NA),
    Onset_sympto= c(NA),
    Progr_chron_dysp_exer= c(NA),
    Cough= c(NA),
    Fever= c(NA),
    Chest_tight= c(NA),
    Dysphagia= c(NA),
    Lab_test_Whi_Blo_Cel= c(NA),
    Raynaud_phenom= c(NA),
    Decreas_DLCO70= c(NA),
    Restr_lung_dis_TLC80= c(NA),
    Oxigen_sat_rest90= c(NA),
    Lips_cianos= c(NA),
    #
    AuscVC_Basal     = c(NA),
    AuscVC_Lateral   = c(NA),
    AuscVC_Higher    = c(NA),
    Other_abnor_lun_soun = c(NA),
    Third_heart_sound = c(NA),
    AuscVC_Sym_same_lev = c(NA),
    #
    EE_AuscVC_Basal     = c(NA),
    EE_AuscVC_Lateral   = c(NA),
    EE_AuscVC_Higher    = c(NA),
    EE_AuscVC_Sym_same_lev = c(NA)
    );
    
  if(input$Age_above60 != "not-informed"){
         res$Age_above60 <-  input$Age_above60
    } 
  if(input$Sex  != "not-informed"){
     res$Sex  <-  input$Sex
  } 
  if(input$Exposures  != "not-informed"){
    res$Exposures  <-  input$Exposures
  } 
  if(input$ActiEx_Smoke  != "not-informed"){
    res$ActiEx_Smoke  <-  input$ActiEx_Smoke
  } 
  if(input$GERD  != "not-informed"){
    res$GERD  <-  input$GERD
  } 
  if(input$Histo_autoim_rheuma  != "not-informed"){
    res$Histo_autoim_rheuma  <-  input$Histo_autoim_rheuma
  } 
  if(input$Histo_HF_cardi_dis  != "not-informed"){
    res$Histo_HF_cardi_dis  <-  input$Histo_HF_cardi_dis
  } 
  if(input$Pneumotox_medicat  != "not-informed"){
    res$Pneumotox_medicat  <-  input$Pneumotox_medicat
  } 
  if(input$Asthma  != "not-informed"){
    res$Asthma  <-  input$Asthma
  } 
  if(input$Histo_pulmo_dis  != "not-informed"){
    res$Histo_pulmo_dis  <-  input$Histo_pulmo_dis
  } 
  if(input$Familia_histo_FILD  != "not-informed"){
    res$Familia_histo_FILD  <-  input$Familia_histo_FILD
  } 
  if(input$Orthopnea  != "not-informed"){
    res$Orthopnea  <-  input$Orthopnea
  } 
  if(input$Per_malleo_edema  != "not-informed"){
    res$Per_malleo_edema  <-  input$Per_malleo_edema
  } 
  if(input$Compla_limit_exerc_toler  != "not-informed"){
    res$Compla_limit_exerc_toler  <-  input$Compla_limit_exerc_toler
  } 
  if(input$Digit_clubb  != "not-informed"){
    res$Digit_clubb  <-  input$Digit_clubb
  } 
  if(input$Onset_sympto  != "not-informed"){
    res$Onset_sympto  <-  input$Onset_sympto
  } 
  if(input$Progr_chron_dysp_exer  != "not-informed"){
    res$Progr_chron_dysp_exer  <-  input$Progr_chron_dysp_exer
  } 
  if(input$Cough  != "not-informed"){
    res$Cough  <-  input$Cough
  } 
  if(input$Fever  != "not-informed"){
    res$Fever  <-  input$Fever
  } 
  if(input$Chest_tight  != "not-informed"){
    res$Chest_tight  <-  input$Chest_tight
  } 
  if(input$Dysphagia  != "not-informed"){
    res$Dysphagia  <-  input$Dysphagia
  } 
  if(input$Lab_test_Whi_Blo_Cel  != "not-informed"){
    res$Lab_test_Whi_Blo_Cel  <-  input$Lab_test_Whi_Blo_Cel
  } 
  if(input$Raynaud_phenom  != "not-informed"){
    res$Raynaud_phenom  <-  input$Raynaud_phenom
  } 
  if(input$Decreas_DLCO70  != "not-informed"){
    res$Decreas_DLCO70  <-  input$Decreas_DLCO70
  } 
  if(input$Restr_lung_dis_TLC80  != "not-informed"){
    res$Restr_lung_dis_TLC80  <-  input$Restr_lung_dis_TLC80
  } 
  if(input$Oxigen_sat_rest90  != "not-informed"){
    res$Oxigen_sat_rest90  <-  input$Oxigen_sat_rest90
  } 
  if(input$Lips_cianos  != "not-informed"){
    res$Lips_cianos  <-  input$Lips_cianos
  } 
  if(input$AuscVC_Basal  != "not-informed"){
    res$AuscVC_Basal  <-  input$AuscVC_Basal
  } 
  if(input$AuscVC_Lateral  != "not-informed"){
    res$AuscVC_Lateral  <-  input$AuscVC_Lateral
  } 
  if(input$AuscVC_Higher  != "not-informed"){
    res$AuscVC_Higher  <-  input$AuscVC_Higher
  } 
  if(input$Other_abnor_lun_soun  != "not-informed"){
    res$Other_abnor_lun_soun  <-  input$Other_abnor_lun_soun
  } 
  if(input$Third_heart_sound  != "not-informed"){
    res$Third_heart_sound  <-  input$Third_heart_sound
  } 
  if(input$AuscVC_Sym_same_lev  != "not-informed"){
    res$AuscVC_Sym_same_lev  <-  input$AuscVC_Sym_same_lev
  } 
  if(input$EE_AuscVC_Basal  != "not-informed"){
    res$EE_AuscVC_Basal  <-  input$EE_AuscVC_Basal
  } 
  if(input$EE_AuscVC_Lateral  != "not-informed"){
    res$EE_AuscVC_Lateral  <-  input$EE_AuscVC_Lateral
  } 
  if(input$EE_AuscVC_Higher  != "not-informed"){
    res$EE_AuscVC_Higher  <-  input$EE_AuscVC_Higher
  } 
  if(input$EE_AuscVC_Sym_same_lev  != "not-informed"){
    
    if(input$EE_AuscVC_Sym_same_lev  == "Yes-Bilateral"){
      res$EE_AuscVC_Sym_same_lev  <-  "Bilateral"
    }else if(input$EE_AuscVC_Sym_same_lev  == "No-Bilateral"){
      res$EE_AuscVC_Sym_same_lev  <-  "Unilateral"
    }else{
      stop("Error in translating evidence on VC detectors");
    }
  } 

 return(res)
}






#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#  

compile_my_bn <-  function(input,session){
  show_modal_spinner( spin = "fulfilling-bouncing-circle", 
                      color = "#0e4aa2",
                      text = "Compilation in progress: please wait ..."
                     ) 
  cpt_all_LS <- get_userData("cpt_all_LS",session)
  
  
  spe_sen_LS <-list(
       EE_AuscVC_Basal = c(speci= input$speciBasal,
                           sensi=input$sensiBasal), 
       EE_AuscVC_Lateral= c(speci=input$speciLateral,
                            sensi=input$sensiLateral),       
       EE_AuscVC_Higher= c(speci=input$speciHigher,
                           sensi=input$sensiHigher),       
       EE_AuscVC_Sym_same_lev= c(speci=input$speciSym,
                                 sensi=input$sensiSym
       ))
   
  set_userData("spe_sen_LS",spe_sen_LS,session)
   
  
  EE_AuscVC_Basal <- cptable( c("EE_AuscVC_Basal","AuscVC_Basal"),
        values = c(# 
          spe_sen_LS[["EE_AuscVC_Basal"]]["sensi"], 1-spe_sen_LS[["EE_AuscVC_Basal"]]["sensi"],#          
          1-spe_sen_LS[["EE_AuscVC_Basal"]]["speci"],  spe_sen_LS[["EE_AuscVC_Basal"]]["speci"] #   
        ), levels= c('YES','NO'));
  #  
  EE_AuscVC_Lateral <- cptable(c("EE_AuscVC_Lateral","AuscVC_Lateral"),
       values = c(# YES, NO
         spe_sen_LS[["EE_AuscVC_Lateral"]]["sensi"], 1-spe_sen_LS[["EE_AuscVC_Lateral"]]["sensi"],#          
         1-spe_sen_LS[["EE_AuscVC_Lateral"]]["speci"], spe_sen_LS[["EE_AuscVC_Lateral"]]["speci"] #   
                           ), levels= c('YES','NO'));
  # 
  EE_AuscVC_Higher <- cptable( c( "EE_AuscVC_Higher","AuscVC_Higher"),
         values = c(#  YES, NO
           spe_sen_LS[["EE_AuscVC_Higher"]]["sensi"], 1- spe_sen_LS[["EE_AuscVC_Higher"]]["sensi"],#           
           1-spe_sen_LS[["EE_AuscVC_Higher"]]["speci"],    spe_sen_LS[["EE_AuscVC_Higher"]]["speci"]  #   
         ), levels= c("YES","NO"));
  #  
  EE_AuscVC_Sym_same_lev <- cptable( c( "EE_AuscVC_Sym_same_lev","AuscVC_Sym_same_lev"),
         values = c(# YES, NO
           spe_sen_LS[["EE_AuscVC_Sym_same_lev"]]["sensi"], 1-spe_sen_LS[["EE_AuscVC_Sym_same_lev"]]["sensi"],#           
           1-spe_sen_LS[["EE_AuscVC_Sym_same_lev"]]["speci"], spe_sen_LS[["EE_AuscVC_Sym_same_lev"]]["speci"]  # 
         ), levels= c("Bilateral","Unilateral"));
  #
  cpt_all_LS[['EE_AuscVC_Basal']] <- EE_AuscVC_Basal
  cpt_all_LS[['EE_AuscVC_Lateral']] <- EE_AuscVC_Lateral
  cpt_all_LS[['EE_AuscVC_Higher']] <- EE_AuscVC_Higher
  cpt_all_LS[['EE_AuscVC_Sym_same_lev']] <- EE_AuscVC_Sym_same_lev
  #  
  plist_3_3 <- compileCPT(cpt_all_LS)
  #  
  net_3_3.0 <- grain(plist_3_3)
  net_3_3.c <-  compile(net_3_3.0)
  remove_modal_spinner()
  net_3_3.c 
}



#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


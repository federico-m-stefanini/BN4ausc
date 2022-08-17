# title: "BN for auscultation"
# author: "Federico M. Stefanini"
# date: "4/8/2020 - rel. 1.3"
# date: "31/8/2020 - rel. 1.4"
# date: "1-2/9/2020 - rel. 1.5"
# date: "1-4/9/2020 - rel. 1.6"
# date: "5/9/2020 - rel. 1.7"   
# date: "10/9/2020 - rel. 1.8"   
# date: "10/9/2020 - rel. 1.9"   
# date: "10/9/2020 - rel. 1.10"   
# date: "11/9/2020 - rel. 1.11"   
# date: "26/10/2020 - rel. 1.12"   
# date: "15/08/2022 - rel. 1.13"   


CURRENT_RELEASE <-  "1.13"
CURRENT_DATE <-  "15/08/2022"#as.character(Sys.Date())


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




source("www/lib_auscu_web_app.R")



shinyApp(
    # UI ---------------
    ui = fluidPage(
        titlePanel("A Bayesian Network for classifying patients after auscultation of velcro crackles"),
        HTML("<h4><strong>Team:</strong>",
             "Auth1, Auth2, Auth3 </h4>"),
         h5(paste0(CURRENT_DATE," ", "- rel. ",CURRENT_RELEASE)),
        
        br(),
        br(),
        shinyjs::useShinyjs(),  
        shinyalert::useShinyalert(),
        # 
        sidebarLayout(
            
        ## Sidebar -----
        shinyjs::hidden(
            div(id = "Sidebar", sidebarPanel(
                    
                     
                    conditionalPanel(
                        condition = "input.tabselected > 1",
                        
                        textInput(inputId ="bnStatus", 
                                  label ="Bayesian Network compiled?", 
                                  value = "NO", 
                                  width = '80%', 
                                  placeholder = NULL),
                        #br(),
                        actionButton("goCompile", 
                                     "Compile BN",
                                     icon= icon(name="hand-right",
                                                lib="glyphicon"),
                                     style="font-weight:bold;",
                                     width="80%"),
                      
                        br(),br(),br(),#br(),br(),
                       textInput(inputId ="propagaStatus", 
                                  label ="Evidence in the Bayesian Network propagated?", 
                                  value = "NO", 
                                  width = '80%', 
                                  placeholder = NULL),
                        #br(),
                        actionButton("goPropagate", 
                                     "Propagate evidence.",
                                     icon= icon(name="hand-right",
                                                lib="glyphicon"),
                                     style="font-weight:bold;",
                                     width="80%"),
                       br(),br(),br(),
                       HTML("<strong>All data provided as input can be saved by clicking on:</strong>"),
                       actionButton("makeSnapshot", 
                                    "Make snapshot.",
                                    icon= icon(name="hand-down",
                                               lib="glyphicon"),
                                    style="font-weight:bold;",
                                    width="80%"),
                       br(),
                       downloadButton("downloadSnapshot", "Download last snapshot",
                                      style="font-weight:bold;",
                                      width="80%"),
                       br(),br(),br(),
                       fileInput("uploadSnapshot",
                                 "A previously saved snapshot  can be loaded by uploading a file .RData",
                                 width="80%",
                                 multiple = FALSE,
                                 accept = ".RData"),
                       HTML("<strong>(wait for the upload to be complete)</strong>"),
                       actionButton("installSnapshot", "Install snapshot",
                                    icon= icon(name="hand-up",
                                               lib="glyphicon"),
                                    style="font-weight:bold;",
                                    width="80%")
                    )# end  
            
                ))), #  
            
            # Main-Panel ------
            mainPanel(
                
                tabsetPanel(
                    id = "tabselected", type = "pills",
                    # Login -------
                    tabPanel("Authentication",
                             value = 1,
                             br(),
                             textInput("username", "Username"),
                             passwordInput("password", label = "Passwort"),
                             actionButton("login", "Login"),
                             textOutput("pwd")
                             
                    )  # closes 
                )  # closes       
            )  # closes                       
        ) # closes 
    ), # closes 
    
    
    
    
    
    
    
  # SERVER --------------------------------------------
  server = function(input, output, session){
        
     load("www/bncomponents_1_1.RData",envir = session$userData)
      allNodes <- c("Disease", "Age_above60", "Sex", "Exposures", "ActiEx_Smoke", 
                    "GERD", "Histo_autoim_rheuma", "Histo_HF_cardi_dis", "Pneumotox_medicat", 
                    "Asthma", "Histo_pulmo_dis", "Familia_histo_FILD", "Orthopnea", 
                    "Per_malleo_edema", "Compla_limit_exerc_toler", "Digit_clubb", 
                    "Onset_sympto", "Progr_chron_dysp_exer", "Cough", "Fever", "Chest_tight", 
                    "Dysphagia", "Lab_test_Whi_Blo_Cel", "Raynaud_phenom", "Decreas_DLCO70", 
                    "Restr_lung_dis_TLC80", "Oxigen_sat_rest90", "Lips_cianos", 
                    "Other_abnor_lun_soun", "Third_heart_sound",
                    "AuscVC_Basal", "AuscVC_Sym_same_lev",
                    "AuscVC_Lateral", "AuscVC_Higher", 
                    #
                    "EE_AuscVC_Basal","EE_AuscVC_Sym_same_lev",
                    "EE_AuscVC_Lateral","EE_AuscVC_Higher")
      set_userData("allNodes",allNodes,session)
      nameNodesActi <- c("Disease",   "EE_AuscVC_Basal", "EE_AuscVC_Sym_same_lev",
                         "AuscVC_Basal","AuscVC_Sym_same_lev",
                         "AuscVC_Lateral", "AuscVC_Higher")
      set_userData("nameNodesActi",nameNodesActi,session)
      nameNodesNoActi <- c("Disease", "Age_above60", "Sex", "Exposures", "ActiEx_Smoke", 
                           "GERD", "Histo_autoim_rheuma", "Histo_HF_cardi_dis", "Pneumotox_medicat", 
                           "Asthma", "Histo_pulmo_dis", "Familia_histo_FILD", "Orthopnea", 
                           "Per_malleo_edema", "Compla_limit_exerc_toler", "Digit_clubb", 
                           "Onset_sympto", "Progr_chron_dysp_exer", "Cough", "Fever", "Chest_tight", 
                           "Dysphagia", "Lab_test_Whi_Blo_Cel", "Raynaud_phenom", "Decreas_DLCO70", 
                           "Restr_lung_dis_TLC80", "Oxigen_sat_rest90", "Lips_cianos", 
                            "Other_abnor_lun_soun", 
                           "Third_heart_sound",
                           #
                           "EE_AuscVC_Lateral", "EE_AuscVC_Higher")
      set_userData("nameNodesNoActi",nameNodesNoActi,session)
      
      ## passwords -----------------------------------
      user_vec <- c(  "golem" = "golem" )
        
 
      
      ## saveSnapshot ----------------------------------------
      observeEvent(input$makeSnapshot, {
           evidenza_corrente <- read_evidences(input,session)
           spe_sen_LS <-list(
             EE_AuscVC_Basal = c(speci= input$speciBasal,
                                 sensi=input$sensiBasal), 
             EE_AuscVC_Lateral= c(speci=input$speciLateral,
                                  sensi=input$sensiLateral),       
             EE_AuscVC_Higher= c(speci=input$speciHigher,
                                 sensi=input$sensiHigher),       
             EE_AuscVC_Sym_same_lev= c(speci=input$speciSym,
                                       sensi=input$sensiSym))
           myNotes <-  input$testoNote
           last_snap_file <-  tempfile(pattern= paste0("wapp_snapshot_",
                                       gsub(":","_",gsub(" ","_",Sys.time())),"_"),
                                       fileext=".RData")
           save(evidenza_corrente,spe_sen_LS,myNotes,
                file=last_snap_file)
           set_userData("last_snap_file", last_snap_file, session)
           #
      })
  
      ## downloadSnapshot ----------------------------------
      output$downloadSnapshot <- downloadHandler(
        filename = function() {
          pathFiDown <- paste0("snapshot",gsub(":","_",gsub(" ","_",Sys.time())),".RData")
          return(pathFiDown)
        },
        content = function(file) {
          last_snap_file <- get_userData("last_snap_file",session)
          return(file.copy(from= last_snap_file, to= file))
        }
      )
      
      
          
       
      observeEvent(input$installSnapshot, {
           tmpLoadEnv <- new.env()
            load(input$uploadSnapshot$datapath, envir = tmpLoadEnv)
          
           reset_evidences(input,session)
           evidenza_corrente <- get("evidenza_corrente",envir=tmpLoadEnv)
           evidenza_corrente <- evidenza_corrente[!is.na(evidenza_corrente)]
           nomi_corrente <- names(evidenza_corrente)
           for(aux in seq_along(evidenza_corrente)){
             if(nomi_corrente[aux] == "EE_AuscVC_Sym_same_lev"){# gestione eccezione 
                  if(evidenza_corrente[[aux]] =="Bilateral"){
                    messaggio <-  "Yes-Bilateral"
                  }else if(evidenza_corrente[[aux]] =="Unilateral"){
                    messaggio <-  "No-Bilateral"
                  }
              evidenza_corrente[[aux]] <- messaggio
              }
              updateSelectInput(session,inputId = nomi_corrente[aux],
                                selected = evidenza_corrente[[aux]])
           }
           #
           spe_sen_LS<- get("spe_sen_LS",envir=tmpLoadEnv)
           spe_sen_LS <- lapply(spe_sen_LS,function(vx)as.numeric(vx))
           updateSliderInput(session,inputId = 'speciBasal', 
                             value=  spe_sen_LS$EE_AuscVC_Basal[1]);
           updateSliderInput(session,inputId = 'speciLateral',
                             value=  spe_sen_LS$EE_AuscVC_Lateral[1]);
           updateSliderInput(session,inputId = 'speciHigher',
                             value= spe_sen_LS$EE_AuscVC_Higher[1] );
           updateSliderInput(session,inputId = 'speciSym',
                             value=  spe_sen_LS$EE_AuscVC_Sym_same_lev[1]);
           updateSliderInput(session,inputId = 'sensiBasal',
                             value= spe_sen_LS$EE_AuscVC_Basal[2]);
           updateSliderInput(session,inputId = 'sensiLateral',
                             value= spe_sen_LS$EE_AuscVC_Lateral[2]);
           updateSliderInput(session,inputId = 'sensiHigher',
                             value= spe_sen_LS$EE_AuscVC_Higher[2]);
           updateSliderInput(session,inputId = 'sensiSym',
                             value=spe_sen_LS$EE_AuscVC_Sym_same_lev[2] );
           updateTextInput(session,
                           inputId="bnStatus",
                           value="NO")
           #
           myNotes<- get("myNotes",envir=tmpLoadEnv)
           updateTextAreaInput(session,inputId="testoNote",value=myNotes)
           
      })
      
      
        
      
      
      
      
      ## saveB button ----------------------------------------
        observeEvent(input$saveB, {
          # check for evidence introduced
          evid_informed <-  read_evidences(input,session)
          evid_informed <- evid_informed[!is.na(evid_informed)]
          if(input$propagaStatus =="NO"){
            shinyalert::shinyalert(title = "Warning!",
                                   text= paste0("Warning: Before preparing the extended report ",
                                                " you must compile the network and propagate evidences! "),
                                   type = "info",
                                   showConfirmButton=TRUE,
                                   timer=0) 
            return(FALSE);
          }
          if(length(evid_informed)<3){
            shinyalert::shinyalert(title = "Warning!",
                                   text= paste0("Warning: the minimum of evidence to prepare an extended report ",
                                                " was not provided. ",
                                                " Besides VCBasal and VCSymmetry nodes at least one more node must be ",
                                                " informed before propagation!"),
                                   type = "info",
                                   showConfirmButton=TRUE,
                                   timer=0) 
            return(FALSE);
          }
          #
          show_modal_spinner( spin = "fulfilling-bouncing-circle", 
                              color = "#0e4aa2",
                              text = "Preparation in progress: please wait ...") 
          
          tavola_1 <- get_userData("tavola_1", session)
          tavola_2 <- get_userData("tavola_2", session)
          figura_1 <- get_userData("figura_1", session)
          figura_2 <- get_userData("figura_2", session) 
          current_evi <- get_userData("current_evi", session) 
          spe_sen_LS <- get_userData("spe_sen_LS",session)
          effe_size_stat <- get_userData("effe_size_stat",session)
         repo_dir <- tempdir(check=TRUE)
          tmp_folder <-  paste0("tmp_",substr(as.character(runif(1)),3,300))
          repo_dir_tmp <- paste0(repo_dir,"/",tmp_folder)
          flag_created <- dir.create(repo_dir_tmp)
          name_BNrep <- paste0("BNrep_",substr(as.character(runif(1)),3,300))
          set_userData("repo_dir",repo_dir,session)
          set_userData("repo_dir_tmp",repo_dir_tmp,session)
          set_userData("name_BNrep",name_BNrep,session)
          set_userData("tmp_folder",tmp_folder,session)
          # 
          argo_report_LS_ini <- list(
            dataNow =  paste0(date()),
            myNotes =  input$testoNote,
            tavola1  = tavola_1,
            tavola2  = tavola_2,
            figura1 =  figura_1,
            figura2 =  figura_2,
            current_evi = current_evi,
            spe_sen_LS = spe_sen_LS,
            effe_size_stat= effe_size_stat,
            release_cur = CURRENT_RELEASE
          ) 
           if(TRUE){# 26 ottobre 2020
              argo_report_LS <-  mapp_aggrega_ILD(argo_report_LS_ini)
          }else{
              argo_report_LS <-  argo_report_LS_ini
          }
          rmarkdown::render(input = "www/BN_auscult_report.Rmd", 
                  output_format = "html_document", 
                  params = argo_report_LS, #  
                  intermediates_dir = repo_dir_tmp,
                   output_dir = repo_dir, 
                    output_file = name_BNrep)
          # output  
          rmarkdown::render(input = "www/BN_auscult_report.Rmd", 
                              output_format = "word_document", 
                            params = argo_report_LS, 
                            intermediates_dir = repo_dir_tmp,
                            output_dir = repo_dir, 
                             
                            output_file =  name_BNrep)

      remove_modal_spinner()
      })
      
      
      ## Download docx ---------------------
      output$downloadDocx <- downloadHandler(
        filename = function() {
            pathFiDown <- paste0("BNreport",gsub(":","_",gsub(" ","_",Sys.time())),".docx")
          return(pathFiDown)
        },
        content = function(file) {
          repo_dir   <- get_userData("repo_dir",session)
          name_BNrep <- get_userData("name_BNrep",session)
          file_da_input <- file.path(repo_dir,paste0(name_BNrep,".docx"))
            return(file.copy( 
                 from= file_da_input, to= file))
        }
      )

      
      ## download HTML -------------------------------
      output$downloadHTML <- downloadHandler(
        filename = function() {
          pathFiDown <- paste0("BNreport",gsub(":","_",gsub(" ","_",Sys.time())),".html")
          return(pathFiDown)
        },
        content = function(file) {
          repo_dir   <- get_userData("repo_dir",session)
          name_BNrep <- get_userData("name_BNrep",session)
          file_da_input <- file.path(repo_dir,paste0(name_BNrep,".html"))
           return(file.copy( from= file_da_input, to= file))
        }
      )
      
      
      
        
        ## go compile -----------------------
        observeEvent(input$goCompile,{
             updateTextInput(session,
                            inputId="bnStatus",
                            value="YES")
            
            compiled_net <- compile_my_bn(input,session)
            set_userData("compiled_net",compiled_net,session = session)
        })
        
        
        
      
      
      
      
      
      
      
      
      
      
        ## goPropagate EVIDENCE  ----------------------
        observeEvent(input$goPropagate, {
              
         
            if(input$bnStatus == "YES"){
                 current_evi <- read_evidences(input, session)
                 if(is.na(current_evi$EE_AuscVC_Sym_same_lev) || is.na(current_evi$EE_AuscVC_Basal)){
                  shinyalert::shinyalert(title = "Warning!",
                           text= paste0("Warning: the minimum of evidence was not provided.",
                                        " Inform VCBasal and VCSymmetry nodes before propagation!"),
                           type = "info",
                           showConfirmButton=TRUE,
                           timer=0)  
                  
                }else{
                     # 
                    propa_stat <- propaga_ev(current_evi,session)
                    updateTextInput(session,
                                inputId="propagaStatus",
                                value="YES")
                     effe_size_stat <- effe_size(current_evi,input, session)
                    set_userData("effe_size_stat", effe_size_stat,session)
                    
                    tmp1 <- as_tibble(propa_stat$tavolaPro)
                    tmp1 <- mutate(tmp1,ID= row.names(propa_stat$tavolaPro),.before=1)
                    tmp11 <- tmp1 %>% gather( 'CHF', 'Pneumonia', 'ILD_IPF', 'ILD_CTD', 'ILD_other', 'LungScars', key="Disease",value="Prob")
                    tmp111 <- tmp11 %>% dplyr::filter(ID == 'Init_prob' | ID =='Final_prob' )
                    mygg1 <- ggplot(tmp111, aes(fill=ID, y=Prob, x=Disease)) + 
                                  geom_bar(position="dodge", stat="identity");
                    tmp2 <- as_tibble(propa_stat$tavolaBF)
                    tmp2 <- mutate(tmp2,ID= row.names(propa_stat$tavolaBF),.before=1)
                    tmp22 <- tmp2[3,] %>% gather( 'CHF', 'Pneumonia', 'ILD_IPF', 'ILD_CTD', 'ILD_other', 'LungScars', key="Disease",value="BFval")
                    mygg2 <- ggplot(tmp22, aes(  y=BFval, x=Disease)) + 
                               geom_bar(stat="identity");
                    
                    set_userData("tavola_1", tmp1,session)
                    set_userData("tavola_2", tmp2,session)
                    set_userData("figura_1", mygg1,session)
                    set_userData("figura_2", mygg2,session) 
                    set_userData("current_evi", current_evi,session) 
                     spe_sen_LS <- list(
                      EE_AuscVC_Basal = c(speci= input$speciBasal,
                                          sensi=input$sensiBasal), 
                      EE_AuscVC_Lateral= c(speci=input$speciLateral,
                                           sensi=input$sensiLateral),       
                      EE_AuscVC_Higher= c(speci=input$speciHigher,
                                          sensi=input$sensiHigher),       
                      EE_AuscVC_Sym_same_lev= c(speci=input$speciSym,
                                                sensi=input$sensiSym))
                    argo_report_LS_ini <- list(
                      dataNow =  paste0(date()),
                      myNotes =  "<none>",
                      tavola1  = tmp1,
                      tavola2  = tmp2,
                      figura1 =  mygg1,
                      figura2 =  mygg2,
                      current_evi = current_evi,
                      spe_sen_LS = spe_sen_LS,
                      effe_size_stat= effe_size_stat,
                      release_cur = CURRENT_RELEASE
                    ) 
                    argo_report_LS <-  mapp_aggrega_ILD(argo_report_LS_ini)
                    # OUTPUT  
                    output$block1 <- renderTable(argo_report_LS$tavola1,digits=4)
                    output$block2 <- renderTable(argo_report_LS$tavola2,digits=4)
                    output$mygg1 <-  renderPlot(argo_report_LS$figura1 )
                    output$mygg2 <-  renderPlot(argo_report_LS$figura2 )
                    
                    
              }# end else
            }else{
                shinyalert::shinyalert(title = "Warning!",
                           text="Warning: click on compile button before propagating evidence.",
                           type = "info",
                           showConfirmButton=TRUE,
                           timer=0)   
                }    
        })
        
        
        ## speci&sensi input ----------------------------------
        observeEvent(input$speciAll, {
          
            valore_spec <- input$speciAll[1]
            updateSliderInput(session,
                              inputId = 'speciBasal',
                              value= valore_spec);
            updateSliderInput(session,
                              inputId = 'speciLateral',
                              value= valore_spec);
            updateSliderInput(session,
                              inputId = 'speciHigher',
                              value= valore_spec);
            updateSliderInput(session,
                              inputId = 'speciSym',
                              value= valore_spec);
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
            updateTextInput(session,
                            inputId="propagaStatus",
                            value="NO")
        })
        observeEvent(input$sensiAll, {
            
            valore_sensi <- input$sensiAll[1]
            updateSliderInput(session,
                              inputId = 'sensiBasal',
                              value= valore_sensi);
            updateSliderInput(session,
                              inputId = 'sensiLateral',
                              value= valore_sensi);
            updateSliderInput(session,
                              inputId = 'sensiHigher',
                              value= valore_sensi);
            updateSliderInput(session,
                              inputId = 'sensiSym',
                              value= valore_sensi);
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
            updateTextInput(session,
                            inputId="propagaStatus",
                            value="NO")
        })        
        observeEvent(input$speciBasal, {    
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
          
        })
        observeEvent(input$speciLateral, {    
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        observeEvent(input$speciHigher, {    
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        observeEvent(input$speciSym, {    
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        #
        observeEvent(input$sensiBasal, {    
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        observeEvent(input$sensiLateral, {    
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        observeEvent(input$sensiHigher, {    
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        observeEvent(input$sensiSym, {    
            updateTextInput(session,
                            inputId="bnStatus",
                            value="NO")
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        
      
        
        ## input EE VC BASAL& EE Sym Lev--------------------------------------------
        observeEvent(input$EE_AuscVC_Basal, {
           
          if(input$EE_AuscVC_Basal  != 'not-informed'  & input$EE_AuscVC_Sym_same_lev != 'not-informed'){
            for(aux in nameNodesNoActi){
              shinyjs::enable(id=aux)
            }  
          }else{
            for(aux in nameNodesNoActi){
              shinyjs::disable(id=aux)
            } 
          }
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        
        observeEvent(input$EE_AuscVC_Sym_same_lev, {
            if(input$EE_AuscVC_Basal  != 'not-informed'  & input$EE_AuscVC_Sym_same_lev != 'not-informed'){
            for(aux in nameNodesNoActi){
              shinyjs::enable(id=aux)
            }  
          }else{
            for(aux in nameNodesNoActi){
              shinyjs::disable(id=aux)
            } 
          }
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")
        })
        
        ## input al other evidences -------------------------------------------
        observeEvent(input$EE_AuscVC_Lateral,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$EE_AuscVC_Higher,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        # 
        observeEvent(input$AuscVC_Basal,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$AuscVC_Lateral,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$AuscVC_Higher,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$AuscVC_Sym_same_lev,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        #
        observeEvent(input$Age_above60,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Sex,{updateTextInput(session, inputId="propagaStatus", value="NO")})  
        observeEvent(input$Exposures,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$ActiEx_Smoke,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$GERD,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Histo_autoim_rheuma,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Histo_HF_cardi_dis,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Pneumotox_medicat,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Asthma,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Histo_pulmo_dis,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Familia_histo_FILD,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Orthopnea,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Per_malleo_edema,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Compla_limit_exerc_toler,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Digit_clubb,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Onset_sympto,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Progr_chron_dysp_exer,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Cough,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Fever,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Chest_tight,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Dysphagia,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Lab_test_Whi_Blo_Cel,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Raynaud_phenom,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Decreas_DLCO70,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Restr_lung_dis_TLC80,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Oxigen_sat_rest90,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Lips_cianos,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Other_abnor_lun_soun,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        observeEvent(input$Third_heart_sound,{updateTextInput(session, inputId="propagaStatus", value="NO")})
        
        
        ## Clear form---------------------
        
        observeEvent(input$clearForm,{ 
           reset_evidences(input,session)
          updateTextInput(session,
                          inputId="propagaStatus",
                          value="NO")  
       })
          
        observeEvent(input$login, {
            
           
                
            shinyjs::show(id = "Sidebar")
 
                  
 
            appendTab(inputId = "tabselected",
                              
                     tabPanel("VC detection",
                          value = 2,
                          br(),
                          p("Select values of sensitivity and specificity ",
                            " while detecting velcro crackles."),
                          p("Changes performed on 'all detections' affect all nodes below"),
                          br(),
                          sliderInput(inputId="speciAll",
                                      label="Specificity for all detections",
                                      min = 0,
                                      max = 1,
                                      value= 0.73,
                                      step=0.01,
                                     
                                      width='75%',
                                      sep=""
                                      
                          ),
                          sliderInput(inputId="sensiAll",
                                      label="Sensitivity for all detections",
                                      min = 0,
                                      max = 1,
                                      value= 0.78,
                                      step=0.01,
                                     
                                      width='75%',
                                      sep=""
                                      
                          ),
                          br(),
                          hr(),
                          p("Node by node alternative specification."),
                          br(),
                          sliderInput(inputId="speciBasal",
                                      label="Specificity for EE_AuscVC_Basal ",
                                      min = 0,
                                      max = 1,
                                      value= 0.73,
                                      step=0.01,
                                     
                                      width='75%',
                                      sep=""
                                      
                          ),
                          sliderInput(inputId="sensiBasal",
                                      label="Sensitivity for EE_AuscVC_Basal ",
                                      min = 0,
                                      max = 1,
                                      value= 0.78,
                                      step=0.01,
                                     
                                      width='75%',
                                      sep=""
                                      
                          ),
                          br(),  br(),      
                          sliderInput(inputId="speciLateral",
                                      label="Specificity for EE_AuscVC_Lateral ",
                                      min = 0,
                                      max = 1,
                                      value= 0.73,
                                      step=0.01,
                                     
                                      width='75%',
                                      sep=""
                                      
                          ),
                          sliderInput(inputId="sensiLateral",
                                      label="Sensitivity for EE_AuscVC_Lateral ",
                                      min = 0,
                                      max = 1,
                                      value= 0.78,
                                      step=0.01,
                                      
                                      width='75%',
                                      sep=""
                                      
                          ),
                          br(),  br(),      
                          sliderInput(inputId="speciHigher",
                                      label="Specificity for EE_AuscVC_Higher ",
                                      min = 0,
                                      max = 1,
                                      value= 0.73,
                                      step=0.01,
                                     
                                      width='75%',
                                      sep=""
                                      
                          ),
                          sliderInput(inputId="sensiHigher",
                                      label="Sensitivity for EE_AuscVC_Higher ",
                                      min = 0,
                                      max = 1,
                                      value= 0.78,
                                      step=0.01,
                                      
                                      width='75%',
                                      sep=""
                                      
                          ),
                          br(),br(),      
                          sliderInput(inputId="speciSym",
                                      label="Specificity for EE_AuscVC_Sym_same_lev ",
                                      min = 0,
                                      max = 1,
                                      value= 0.73,
                                      step=0.01,
                                      
                                      width='75%',
                                      sep=""
                                      
                          ),
                          sliderInput(inputId="sensiSym",
                                      label="Sensitivity for EE_AuscVC_Sym_same_lev ",
                                      min = 0,
                                      max = 1,
                                      value= 0.78,
                                      step=0.01,
                                      
                                      width='75%',
                                      sep=""
                                      
                          ),
                          br(),      
                          hr()
                      ) # closes 
                    )
            
 
    ## Evidence UI -------------------------
            
            appendTab(inputId = "tabselected",
                 tabPanel("Evidence",
                     value = 3,
                    br(),
                    p('Select evidences to update class attribution. ',
                      'Please note that the minimum information is made by:'),					
					          HTML('<ol>',
                          '<li>Auscultation VC Symmertry on same level: Bilateral-YES or Bilateral-NO.</li>',
                          '<li>Auscultation VC Basal: YES or NO. </li>',
                          '</ol>',
          					     '<p> In the specific case in which Auscultation VC Basal = NO, ',
          					     ' the network should be informed either by setting :',
          					     '<ul>',
                         '<li>[EITHER] Auscultation VC Lateral = YES/NO.</li>',
          					     '<li>[OR] Auscultaiton VC Higher = YES/NO.	</li>',
                                   '</ul></p>',
					               '<p>Please start entering evidence with Basal and Symmetry nodes to activate ',
					               ' all other nodes.</p>'
					               ),
					          HTML('<br>To reset the form, please click the "Clear form" button:'),
					          actionButton(inputId="clearForm",
					                       label="Clear form (evidence)",
					                       icon= icon(name="hand-down",
					                                  lib="glyphicon"),
					                       style="font-weight:bold;",
					                       width="200px"),
					          br(),br(),br(),
					selectInput(inputId="EE_AuscVC_Basal",
					            label="EE_AuscVC_Basal",
					            choices = c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE), 
					selectInput(inputId="EE_AuscVC_Sym_same_lev",
					            label="EE_AuscVC_Sym_same_lev",
					            choices =  c('not-informed',"Yes-Bilateral","No-Bilateral"),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),
					br(),
					shinyjs::disabled(
  					# detector 
  					selectInput(inputId="EE_AuscVC_Lateral",
  					            label="EE_AuscVC_Lateral",
  					            choices =   c('not-informed','YES','NO'),
  					            selected = "not-informed",
  					            width="420%",
  					            multiple = FALSE,
  					            selectize = TRUE),  
  					selectInput(inputId="EE_AuscVC_Higher",
  					            label="EE_AuscVC_Higher",
  					            choices =    c('not-informed','YES','NO'),
  					            selected = "not-informed",
  					            width="420%",
  					            multiple = FALSE,
  					            selectize = TRUE)
          ),# end 
					br(),
					br(),
					hr(),
					shinyjs::disabled(
					  selectInput(inputId="Age_above60",
					            label="Age_above60",
					            choices = c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),
					selectInput(inputId="Sex",
					            label="Sex",
					            choices = c('not-informed','Male','Female'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Exposures",
					            label="Exposures",
					            choices = c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="ActiEx_Smoke",
					            label="ActiEx_Smoke",
					            choices = c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="GERD",
					            label="GERD",
					            choices =   c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Histo_autoim_rheuma",
					            label="Histo_autoim_rheuma",
					            choices =   c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Histo_HF_cardi_dis",
					            label="Histo_HF_cardi_dis",
					            choices = c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Pneumotox_medicat",
					            label="Pneumotox_medicat",
					            choices = c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Asthma",
					            label="Asthma",
					            choices = c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Histo_pulmo_dis",
					            label="Histo_pulmo_dis",
					            choices =  c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Familia_histo_FILD",
					            label="Familia_histo_FILD",
					            choices =  c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Orthopnea",
					            label="Orthopnea",
					            choices = c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Per_malleo_edema",
					            label="Per_malleo_edema",
					            choices =  c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Compla_limit_exerc_toler",
					            label="Compla_limit_exerc_toler",
					            choices =  c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Digit_clubb",
					            label="Digit_clubb",
					            choices =  c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Onset_sympto",
					            label="Onset_sympto",
					            choices =  c('not-informed','FAST','SLOW') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Progr_chron_dysp_exer",
					            label="Progr_chron_dysp_exer",
					            choices = c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Cough",
					            label="Cough",
					            choices =  c('not-informed','NO','DRY','PRODUCTIVE'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Fever",
					            label="Fever",
					            choices =  c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Chest_tight",
					            label="Chest_tight",
					            choices =  c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Dysphagia",
					            label="Dysphagia",
					            choices = c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Lab_test_Whi_Blo_Cel",
					            label="Lab_test_Whi_Blo_Cel",
					            choices =   c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Raynaud_phenom",
					            label="Raynaud_phenom",
					            choices =   c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Decreas_DLCO70",
					            label="Decreas_DLCO70",
					            choices =   c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Restr_lung_dis_TLC80",
					            label="Restr_lung_dis_TLC80",
					            choices =   c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Oxigen_sat_rest90",
					            label="Oxigen_sat_rest90",
					            choices =    c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Lips_cianos",
					            label="Lips_cianos",
					            choices =  c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE)   
					),
					#
					# 
					shinyjs::hidden(
					  selectInput(inputId="AuscVC_Basal",
					            label="AuscVC_Basal",
					            choices =    c('not-informed','YES','NO'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),
					  selectInput(inputId="AuscVC_Sym_same_lev",
					              label="AuscVC_Sym_same_lev",
					              choices =  c('not-informed',"Bilateral","Unilateral"),
					              selected = "not-informed",
					              width="420%",
					              multiple = FALSE,
					              selectize = TRUE),
					selectInput(inputId="AuscVC_Lateral",
					            label="AuscVC_Lateral",
					            choices =  c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="AuscVC_Higher",
					            label="AuscVC_Higher",
					            choices =   c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE)
					),# 
					#
					shinyjs::disabled(
					selectInput(inputId="Other_abnor_lun_soun",
					            label="Other_abnor_lun_soun",
					            choices =    c('not-informed','Wheezing','Squacks',
                                                'Pleur_frictio','None'),
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE),  
					selectInput(inputId="Third_heart_sound",
					            label="Third_heart_sound",
					            choices =   c('not-informed','YES','NO') ,
					            selected = "not-informed",
					            width="420%",
					            multiple = FALSE,
					            selectize = TRUE)  
					#
					#
					 
					),  # 
					br(),
					br(),
					br(),
					br(),
					br(),
					br() 
					)# 
        )
            
            
            
            
            ## Output tab UI -------------------------------------
            appendTab(inputId = "tabselected",
                 tabPanel("Output",
                     value = 4,
                     br(),
                     p("Here initial probability values are defined with respect to no  evidence ",
                     " (null evidence)."),
                     br(),
                     p("Prior and posterior probability values."),
                     tableOutput("block1"),
                     br(),
                     p("Bayes Factor values."),
                      tableOutput("block2"),
                     br(),
                     br(),
                     plotOutput(outputId="mygg1", width='500px', height ='400px'),
                     br(),
                     br(),
                     plotOutput(outputId="mygg2", width='500px', height ='400px')
                 ) #       
            )# end 
            
            
            ##   UI --------------------------------------
            appendTab(inputId = "tabselected",
                      tabPanel("Save",
                               value = 6,
                               br(),
                                p("Before saving results of the analysis, be sure that you compiled and propagated the desired evidences ",
                                 "so that in the 'Output' tab results are ready to be saved. "),
                               p("If you want to add some notes, please fill the space below."),
                               textAreaInput("testoNote",
                                             label="Your notes:",
                                             placeholder = "Fill with your notes before saving...",
                                             width="500pt",
                                             height = "200pt",
                                             resize="both"),
                               br(),
                               p("Click on 'Prepare summary' and then on 'Download' buttons below."),
                               br(),
                              
                               actionButton(inputId="saveB",
                                            label="Prepare summary",
                                             icon= icon(name="envelope",
                                                        lib="glyphicon"),
                                             style="font-weight:bold;",
                                             width="200px"),
                               br(),br(),
                              
                               downloadButton("downloadHTML", "Download last HTML",
                                              style="font-weight:bold;",
                                              width="300px"),
                               downloadButton("downloadDocx", "Download last docx",
                                              style="font-weight:bold;",
                                              width="300px"),
                               br(),br(),br()
                      ) # closes tabPanel 
            )
            
            
            ## Help  UI ---------------------------------------
            appendTab(inputId = "tabselected",
                      tabPanel("Help",
                               value = 5,
                               br(),
                               h4("Instructions."),
                               p("(1) Select the 'VC detection' tab by clicking on it. Use sliders 'Specificity for all detections' ",
                                 "and 'Sensitivity for all detections' to define those parameters for all detectors at once. ",
                                 "Revise those parameters node by node as needed using sliders from ",
                                 "'Specificity for EE_AuscVC_Basal' down to the end of this first tab."),
                               p("(2) Compile the BN by clicking on the 'Compile BN' button on the left panel."),
                               p("(3) Click on the 'Evidence' tab and enter known values at each node."),
                               p("(4) Click on 'Propagate evidence' button on the left panel to calculate results given observed evidence."),
                               p('(5) Click on "Output" tab and read  results with reference "no evidence".'),
                               p('(6) Make notes and extend the analysis by clicking on "Prepare Summary". '),
                               p("(7) Save the extended report in docx and html formats clicking on buttons.")
                      ) # closes          
            )# end 
            
            
                    
            removeTab(inputId = "tabselected",target = "1")
            shinyjs::disable(id="bnStatus")
            shinyjs::disable(id="propagaStatus")
           
            
        }) # 
        
        
    } # 
) # Closes 

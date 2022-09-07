# BN4ausc

This open source  Shiny App  has been developed for classifying patients
after auscultation of velcro crackles using a Bayesian Network.

The context information came from several experts, and some of them
also took care of developing     an automated web interview for 
the validation of the Bayesian Network.


I developed all the code uploaded here, but I also collaborated with
the research team on  the statistical side of this project, 
for example by  preparing a training seminar to introduce physicians
to the Bayesian elicitation of probability values.
The preliminary elicitation of structural features of a Bayesian graphical model
was performed with the expert prof.Dr.Wuertemberger.
I recently developed a Bayesian  approach for the statistical analysis of
validation data produced by the planned validation study: 

*  Stefanini F.M., 2021, Evaluating heterogeneity of agreement with strong
prior information, 
pp 1630-1635,  In: SIS2021, ed. Perna C., Salvati N., F. Schirripa Spagnolo.
PEARSON publisher,
(WWW.PEARSON.COM), ISBN 9788891927361        
   

I want to thank the following members of the research team that provided
  their expertise in many steps  of this work:   
       

* Prof.Dr.Gebhard Wuertemberger:  Evaluation of  contextual medical
 information   covering epidemiological and clinical  data; definition of
 evidences and categories, including conversion  to discrete scales
 for the Bayesian Network;  elicitation of   large   Directed Acyclic Graphs  and
selection of  candidate reductions suited for  classification purposes;   development
 of 10 patient-cases  added to the validation dataset after refinement.  

* Dr.Norbert Metzdorf: Reviewed the epidemiological data and the conversion of patient-cases to categories and associated state of evidences.    

* Dr.Armin Furtwaengler: Reviewed the patient-cases.   

* Dr.Giacomo Sgalla: Designed 10 patient-cases and reviewed the classes and states.  

* Dr.Shaun Bender: Co-developed the validation step and  the statistical analysis
 of  validation data.    

* Prof.Dr.Giancarlo Crocetti:  Implemented a companion Website for validation,
  which is not covered in the development of the BN, but enables its ongoing 
  validation.    

* Roberto Slepetys: Project lead in Boehringer Ingelheim International GmbH, 
proponent of the application of Bayesian Network as explanatory tool for 
velcro crackle findings in lung auscultation,
  performed the literature search covering epidemiological data, clinical
 evidences and categories, conversion of patient-cases to categories and 
associated state of evidences and
  co-designed the validation framework and web based tool.    
  
Thanks are also due to  Boehringer Ingelheim International GmbH, 
Germany, for funding this part
of  a larger  project dealing with velcro crackles.      



Please note that the included License applies only to the code   I have developed, 
not  to other R packages needed to run it. 
Last, this is  R code developed for research purposes,
and   is not yet ready (validated) for any use in the medical field.


 

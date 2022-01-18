cols <- colnames(Drugresponse_DepMapID_HistSubtype_SampleInfo)
drugs <- cols[2:4687]
epidrugs <- read.csv("epidrugs.csv", col.names = FALSE)
epidrugs <- epidrugs$FALSE.

compare <- function(dataset, label, level1, level2) {
  label <-substitute(label)
  
  #subset datamatrix to only include the two levels of interest
  data <- subset(dataset, eval(label) == level1 | eval(label) ==  level2) 
  
  #initialize empty dataframe to store test results
  results <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(results) <- c('Drug', 'p-value', "testType")
  
  #runs the below code for every drug in the data matrix
  for (drug in drugs) {
    
    #see if data is normally distributed
    normal_test <- shapiro.test(data[[drug]]) #p>0.05, t-test; p<0.05, wilcoxon test
    
    #if p-value from above test is less than 0.05, run a t-test between the two labels of interest, and append results to dataframe
    if (normal_test$p.value < 0.05) {
      test = t.test(as.formula(paste(drug, "~ HistSubtype")), data = data, var.equal = TRUE)
      testType = "t-test"
      results[nrow(results) + 1,] = c(drug, test$p.value, testType)
      
    }
    
    #otherwise, do a wilcoxon test, and append results to dataframe
    else {
      test = wilcox.test(as.formula(paste(drug, "~ HistSubtype")), data = data)
      testType = "wilcoxon"
      results[nrow(results) + 1,] = c(drug, test$p.value, testType)
      
    }
    
  }
  
  #create folder to store results
  dir.create(file.path("results", paste(sep = "", label, "-", level1, "vs", level2 )), showWarnings = FALSE)
  
  #load in drug annotation file
  info <- read.csv("primary-screen-replicate-treatment-info_DrugIDs.csv")
  
  
  #create dataframe displaying results only for epigenetic drugs of interest
  z = results
  epiresults <- as.data.frame(z[(z$Drug %in% epidrugs),])
  
  #format drug names to match drug info file
  results$Drug <- gsub('\\.', '-', results$Drug)
  results$Drug <- gsub('--.*', '', results$Drug)
  
  
  epiresults$Drug <- gsub('\\.', '-', epiresults$Drug)
  epiresults$Drug <- gsub('--.*', '', epiresults$Drug)
  
  #Merge drug test tables and drug annotation table by broad ID
  info$Drug <- info$broad_id
  results<-merge(results, info[, c("Drug", "name", "dose", "perturbation_type", "screen_id", "purity", "moa", "target", "disease.area", "indication", "smiles", "phase")], by="Drug")
  results<-results[!duplicated(results), ]
  epiresults<-merge(epiresults, info[, c("Drug", "name", "dose", "perturbation_type", "screen_id", "purity", "moa", "target", "disease.area", "indication", "smiles", "phase")], by="Drug")
  epiresults<-epiresults[!duplicated(epiresults), ]
  
  
  #create object (list) containing a dataframe displaying test results for all drugs, and a dataframe of test results only for epigenetic drugs
  res = list(results, epiresults)
  
  #save results to csv files within respective folder
  write.csv(res[[1]], file.path("results", paste(sep = "", label, "-", level1, "vs", level2 ), "ALL.csv"))
  write.csv(res[[2]], file.path("results", paste(sep = "", label, "-", level1, "vs", level2 ), "EPI.csv"))
  
  
  #return the object for downstream analysis
  return(res)
  
}

#runs the above function. In this test: 
#Drugresponse_DepMapID_HistSubtype_SampleInfo is the dataset to be analyzed
#HistSubtype is the variable label to be tested within
#"LinkerMut" and "CoreMut" are the two groups within that variable that will be tested for differences. (IMPT: make sure these two parameters are passed in as strings)
linkervscore <- compare(Drugresponse_DepMapID_HistSubtype_SampleInfo, HistSubtype, "LinkerMut", "CoreMut")



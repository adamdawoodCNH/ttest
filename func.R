cols <- colnames(Drugresponse_DepMapID_HistSubtype_SampleInfo)
drugs <- cols[2:4687]
epidrugs <- read.csv("epidrugs.csv", col.names = FALSE)
epidrugs <- epidrugs$FALSE.

compare <- function(dataset, label, level1, level2) {
  label <-substitute(label)
  
  #subset datamatrix to only include the two levels of interest
  data <- subset(dataset, eval(label) == level1 | eval(label) ==  level2) 
  
  #initialize empty dataframe to store test results
  results <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(results) <- c('Drug', 'p-value')
  
  #runs the below code for every drug in the data matrix
  for (drug in drugs) {
    
    #see if data is normally distributed
    normal_test <- shapiro.test(data[[drug]]) #p>0.05, t-test; p<0.05, wilcoxon test
    
    #if p-value from above test is less than 0.05, run a t-test between the two labels of interest, and append results to dataframe
    if (normal_test$p.value < 0.05) {
      test = t.test(as.formula(paste(drug, "~ HistSubtype")), data = data, var.equal = TRUE)
      results[nrow(results) + 1,] = c(drug, test$p.value)
      
    }
    
    #otherwise, do a wilcoxon test, and append results to dataframe
    else {
      test = wilcox.test(as.formula(paste(drug, "~ HistSubtype")), data = data)
      results[nrow(results) + 1,] = c(drug, test$p.value)
      
    }
    
  }
  
  #create dataframe displaying results only for epigenetic drugs of interest
  z = results
  epiresults <- as.data.frame(z[(z$Drug %in% epidrugs),])
  
  #create object (list) containing a dataframe displaying test results for all drugs, and a dataframe of test results only for epigenetic drugs
  res = list(results, epiresults)
  
  #return the object
  return(res)
  
}


linkervscore <- compare(Drugresponse_DepMapID_HistSubtype_SampleInfo, HistSubtype, "LinkerMut", "CoreMut")


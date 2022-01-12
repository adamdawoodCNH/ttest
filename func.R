cols <- colnames(Drugresponse_DepMapID_HistSubtype_SampleInfo)
drugs <- cols[2:4687]
epidrugs <- read.csv("epidrugs.csv", col.names = FALSE)
epidrugs <- epidrugs$FALSE.

compare <- function(dataset, label, level1, level2) {
  label <-substitute(label)
  
  data <- subset(dataset, eval(label) == level1 | eval(label) ==  level2)
  results <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(results) <- c('Drug', 'p-value')
  
  for (drug in drugs) {
    
    normal_test <- shapiro.test(data[[drug]]) #p>0.05, t-test; p<0.05, wilcoxon test
    
    if (normal_test$p.value < 0.05) {
      test = t.test(as.formula(paste(drug, "~ HistSubtype")), data = data, var.equal = TRUE)
      results[nrow(results) + 1,] = c(drug, test$p.value)
      
    }
    else {
      test = wilcox.test(as.formula(paste(drug, "~ HistSubtype")), data = data)
      results[nrow(results) + 1,] = c(drug, test$p.value)
      
    }
  }
  
  epiresults <- as.data.frame(z[ , (colnames(z) %in% epidrugs)])
  
  res = list(results, epiresults)
  
  return(res)
  
}

linkervscore <- compare(HistSubtype, "LinkerMut", "CoreMut")


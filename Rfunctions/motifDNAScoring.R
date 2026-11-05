####DNA binding energy calculation according to ProBound/PWM motif model####

#How to use:
#install RJSONIO with install.packages('RJSONIO')
#safe this file and keep its path
#At the top of whichever script the function is used, add the line:
#source('PATH_TO_THIS_SCRIPT/motifDNAScoring.R)
#remember to delete the test code before using

#Mandatory input:
#DNAList: a List of DNA sequences with length larger than the motif length, and only containing 'ATGC'
#PWM / modelFile: a PWM with rows representing ddG or a JSON file output from probound usually named fit.models.consensus.json
#k: length of the motif/kmer

#Optional input:
#lflank/rflank: added flank to the DNAList, works the best with the least favorable bases at the ends of the motif,
#               used for frame shifting optimization
#mode: the index of model to be fit if a multi-model probound model is input

#Work flow
#????????delete all codes below this line and above the other indicator line after testing????????#
#ProBoundModel <- '~/fit.models.consensus.json'#change path to wherever the JSON file is located
#testSeqs <- c('CAGTTAAGCTAAATCATTTGTGCTGTCAATGGGCGGATTTGTGGGTCCTCGCAAGGGTTA', #arbitrary DNA List for testing
#              'TTGCCCAAAAGTGACTGCTTTTGTGGTTATGAAAATTAAATTTGTTTAATTTCCGTTTCATTGGTTTTT',
#              'TCAGTCAGCCATGCGGCTTCCGAGCGAAAAGTCCAAGAGGCAGT',
#              'ATTCGCATTGAGCTTAGACCCACATGGAC',
#              'ATCACAATTCGACTAACGAGAGCCCCACTGGGGCCAAACGTCATAAATCCACTAATAAACAGTAACCACAGTAACCCGGCCAGTATG',
#              'CGGAAAGAGCAGCCGAGGGAGAAGAGAACAGCAGAATAGCAGAACAGCCCGCCAGAGGTAATTTGAGATTTTC')
#pwm <- pb2pwm(ProBoundModel, mode = 1)
#seqsTable <- PWM2Score(testSeqs, pwm) #output has 3 cols. seq: the input sequence:
                                      #kmer: the highest scoring kmer frame; scores: the energy score calculated
#To compare scores with different pwms, divide the output scores with the maximum score that can be output by a model
#To get the maximum score of a model, run pbMaxScore for a ProBound model, pwmMaxScore for a pwm model
#pbMax <- pbMaxScore(ProBoundModel, mode = 1)
#pwmMax <- pwmMaxScore(pwm)
#the two values should be the same in this case

#To normalize the energy scores
#seqsTable$NormalizedScores <- seqsTable$scores/pbMax
#????????delete all codes above this line and below the other indicator line after testing????????#

####make PWM####
pb2pwm <- function(modelFile, mode = 1){
  JSON <- readLines(modelFile)
  JSON_matrix <- JSON2Matrix(JSON_text = JSON, expo = exp(1), mode = mode)
  JSON_matrix <- log(JSON_matrix)
  pwm <- apply(JSON_matrix, 2, function(column) column - mean(column))
  return(pwm)
}

####make seqTable####
makeSeqTable <- function(DNAList, rflank = 'GA',
                         lflank = 'CC',
                         k = 12){
  AllSeqs <- DNAList
  exSeqs <- paste(rflank, AllSeqs, lflank, sep = '')
  seqsTable <- matrix(ncol = 2, nrow = sum(nchar(exSeqs)-k+1), data = 'none')
  for(i in 1:length(exSeqs)){
    seq <- exSeqs[i]
    seqList <- c()
    for(j in 1:(nchar(seq)-k+1)){
      seqList <- c(seqList, substr(seq, j, j+k-1))
    }
    add <- min(which(seqsTable[,1] == 'none'))
    seqsTable[add:(add-1+length(seqList)),1] <- rep(seq, length(seqList))
    seqsTable[add:(add-1+length(seqList)),2] <- seqList
  }
  #print('seqsTable Ready')
  seqsTable <- as.data.frame(seqsTable)
  colnames(seqsTable) <- c('seq','kmer')
  return(seqsTable)
}

####PWM2Score####
PWM2Score <- function(DNAList,
                      pwm,
                      rflank = 'GA',
                      lflank = 'CC',
                      k = 12,
                      includeSeq = FALSE){
  seqsTable <- makeSeqTable(DNAList, rflank, lflank, k)
  scores <- pwmScoring(pwm,seqList = seqsTable$kmer)
  seqsTable$scores <- scores
  if(!includeSeq){
    seqsTable <- seqsTable[,2:3]
  }
  #for(i in 1:length(seqs)){
  #  seqsTable[seqsTable$seq == seqs[i],'scoresMod'] <- seqsTable[seqsTable$seq == seqs[i],'scores'] - max(seqsTable[seqsTable$seq == seqs[i],'scores'])
  #}
  #seqsTable <- seqsTable[seqsTable$scoresMod == 0,]
  #seqsTable <- seqsTable[,-4]
  return(seqsTable)
}

####pwmScoring####
pwmScoring <- function(pwm, seqList){
  scores <- c(rep(0,length(seqList)))
  for(i in 1:length(seqList)){
    scores[i] <- scoreSeq(seqList[i], pwm)
  }
  return(scores)
}
####pbScoring####
pbScoring <- function(modelFile, seqList, mode = 1){
  JSON <- readLines(modelFile)
  JSON_matrix <- JSON2Matrix(JSON_text = JSON, expo = exp(1), mode = mode)
  JSON_matrix <- log(JSON_matrix)
  pwm <- apply(JSON_matrix, 2, function(column) column - mean(column))
  scores <- c(rep(0,length(seqList)))
  for(i in 1:length(seqList)){
    scores[i] <- scoreSeq(seqList[i], pwm)
  }
  return(scores)
}

####pbMaxScore####
pbMaxScore <- function(modelFile, mode = 1){
  JSON <- readLines(modelFile)
  JSON_matrix <- JSON2Matrix(JSON_text = JSON, expo = exp(1), mode = mode)
  JSON_matrix <- log(JSON_matrix)
  pwm <- apply(JSON_matrix, 2, function(column) column - mean(column))
  maxScore <- 0
  for(i in 1:ncol(pwm)){
    maxScore <- maxScore + max(pwm[,i])
  }
  return(maxScore)
}
####pwmMaxScore####
pwmMaxScore <- function(pwm){
  maxScore <- 0
  for(i in 1:ncol(pwm)){
    maxScore <- maxScore + max(pwm[,i])
  }
  return(maxScore)
}


####function: scoreSeq####
scoreSeq <- function(seq, pwm){
  if(nchar(seq) != ncol(pwm)){
    return(print('not same length'))
  }
  out <- 0
  for(i in 1:nchar(seq)){
    base <- substr(seq, i,i)
    if(base == 'N'){
      out <- out + min(pwm[,i])
    }else{
      out <- out + pwm[base,i]
    }
  }
  return(out)
}



####function: JSON2Matrix (read JSON file of motif energy matrix and output as matrix)####
JSON2Matrix <- function(JSON_text, expo = exp(1), mode = 1){
  library(RJSONIO)
  scoring <- RJSONIO::fromJSON(JSON_text)
  scoring_matrix <- scoring$coefficients$bindingModes[[mode + 1]]$mononucleotide
  DN_seq <- strsplit(scoring$modelSettings$letterOrder, "")[[1]]
  N_pos <- scoring$modelSettings$bindingModes[[mode+1]]$size
  if(expo < 0){
    out_matrix <- matrix(data = scoring_matrix, nrow = 4, ncol = N_pos, byrow = FALSE)
  }else{
    out_matrix <- matrix(data = expo^scoring_matrix, nrow = 4, ncol = N_pos, byrow = FALSE)
  }
  colnames(out_matrix) <- paste('P',c(1:N_pos),sep = '')
  rownames(out_matrix) <- DN_seq
  return(out_matrix)
}

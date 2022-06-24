#1. RNA folding

seq="AAUACUCCGUUGCAGCAU"
RNAseq = strsplit(seq, "")[[1]] #Create vector of the sequence

path=c()
index_pairs=c()

#Find the WC pair and its index 
get_pair <- function(first, second, rna_df){
  if (rna_df[first, second] != 0){
    pair = paste(rnames[first], "-", rnames[second])
    index_pair = paste(first, second)
    return(c(index_pair, pair))
  }
}

#Traceback function
traceback <- function(trace_df, rna_df, row = 1, col = N, index_pairs = c(), path=c()){
  trace = as.integer(strsplit(trace_df[row,col], ",")[[1]])
  
  for (tr in trace){
    if (identical(tr, numeric(0))){
      return(NULL) 
    }
    if (tr == 100){
      print(index_pairs)
      print(path)
      print("Alternate path###########")
      
    } else if (tr == 1){ #First case: left box
      index_pair = paste(row, col)
      if (!(index_pair %in% index_pairs)){
        index_pairs = c(index_pairs, index_pair)}
      traceback(trace_df, rna_df, row, col-1, index_pairs, path)
      next
      
    } else if (tr == 2){ #Second case: down box
      index_pair = paste(row, col)
      if (!(index_pair %in% index_pairs)){
        index_pairs = c(index_pairs, index_pair)}
      traceback(trace_df, rna_df, row+1, col, index_pairs, path)
      next 
      
    } else if ((tr == 3) ) { #Third case: down-left box
      a = get_pair(row, col, rna_df)
      index_pair = a[1]
      if (!(index_pair %in% index_pairs)){
        index_pairs = c(index_pairs, index_pair)}
      pair = a[2]
      if (!(pair %in% path)){
        path = c(path, pair)}
      traceback(trace_df, rna_df, row+1, col-1, index_pairs, path)
      next
      
    }
    
    
    else if (tr == 4 ) { #??? Forth case
      
      k = as.double(strsplit(trace_df[row,col], " ")[[1]][2])
      get_pair(k, col, rna_df)
      traceback(trace_df, rna_df, k, col)
      get_pair(row, k-1, rna_df)
      traceback(trace_df, rna_df, row, k-1, index_pairs, path)
    }
  }
}



N = length(RNAseq) 


restrain = 5 #our restrain
wc_bonds = c('GC', 'CG', 'AU', 'UA') #WC bonds

#Create the column/row names of the matrix
rnames=c()
for (row in 1:N){
  rnames = c(rnames, paste(RNAseq[row], row, sep=""))
}
cnames = rnames

#Initiate the crude energy matrix and the traceback matrix with a value >> 0 (100 in our case) where is needed
rna_matrix = matrix(100, N, N) 
traceback_matrix = matrix("100", N, N)
##Convert crude energy matrix to dataframe and add col and row names
rna_dataframe = data.frame(rna_matrix)
colnames(rna_dataframe) <- cnames
rownames(rna_dataframe) <- cnames

traceback_df = rna_dataframe

for (row in ((N)-restrain):1){
  for (col in (row+restrain):(N)){
    rna_dataframe[row, col]=0
  }
}

traceback_df = traceback_matrix
rna_dataframe



#Fill the crude energy table and the traceback table
for (row in ((N)-restrain):1){
  #col <- columns
  for (col in (row+restrain):(N)){
    k_vector = 0
    k_indeces = 0

    #To check k when it is valid
    if ((row) < (col)){
      k_vector = c() #Initialize k_vector of length 4
      k_indeces = c() #
      for (k in (row+2):(col-1)){
        k_vector = c(k_vector, rna_dataframe[k, col] + rna_dataframe[row, (k-1)])
        k_indeces = c(k_indeces, k)
      }
    }
  
  #Find the energy of the bond
  en_bond = 4
  pair = paste(RNAseq[col], RNAseq[row], sep="")
  if (pair %in% wc_bonds){
    en_bond = -4
  }
  if (pair %in% c("GU", "UG")){
    en_bond = 0
  }
  
  #Find the 
  min_k = min(k_vector) #Which is the min k | if k was not valid then =0
  k_index = which.min(k_vector) #Find index of min k
  k_num = k_indeces[k_index] 
  vector = c() #Initialize vector with the 4 B(row,col)
    #Find the min B(row,col)
  vector = c(rna_dataframe[row, (col-1)], rna_dataframe[(row+1), col], rna_dataframe[(row+1), (col-1)] + en_bond, min_k)
  rna_dataframe[row,col] = min(vector)
  #For the traceback table...
  index = which(vector == min(vector))


  if ("4" %in% index){ #If min is with the 4th case
    new_4 = paste("4", k_num)
    replace(index, index=="4", new_4)
    traceback_df[row,col] =  index #Store 4 (the method number) of min k and ???????????
  }
  
  else {
    traceback_df[row,col] = toString(index)
  }
}
}



traceback(traceback_df, rna_dataframe)

rna_dataframe

traceback_df



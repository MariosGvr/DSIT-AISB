#2. c-RMSD and d-RMSD


# Loading package
library(wordspace)
library(hash)

#c-RMSD
datasets_path = "~/MASTER/2nd Semester/Alg. in Structural Biology"  #Change this if necessary
setwd(datasets_path) #Change working directory

h <- hash() #Create a hash table instance in which we are going to store every conformation
conf_counter = 1
atom_counter = 0
line_counter = 1
count=0
#Read the file and create the hash table
processFile = function(filename) {
  total_atoms = Inf
  total_conformations = Inf
  con = file(filename, "r")
  while ( TRUE ) {
    if (atom_counter==total_atoms){ #Check if the conformation is completed, add it to the hash table and initiate a new one
      count = count +1
      .set(h, keys=conf_counter, values=conformation)
      conformation <- matrix(nrow=total_atoms, ncol=3)
      conf_counter = conf_counter +1
      atom_counter = 0
      line_counter = 3
    }
    line = readLines(con, n = 1)
    
    #Find how many conformations exist in this file
    if (line_counter==1){
      total_conformations = as.double(line)
    }
    
    #Find how many atoms per conformation
    if (line_counter==2){
      total_atoms = as.double(line)
      conformation <- matrix(nrow=total_atoms, ncol=3)
    }
    
    #Add coordinates to the conformation table
    if (line_counter>=3){
      line = unlist(lapply(strsplit(line, split = "\t"), as.numeric))
      if ( conf_counter>total_conformations ) {
        break
      }
      conformation[line_counter-2,] <- line
            atom_counter = atom_counter + 1
    }
    line_counter = line_counter +1
    
  }
close(con)
return(h)
}
h = processFile("80_conformations.txt")


#c-RMSD

##Calculate non-optimal c-RMSD  between all 3160 combinations of conformation
result_cRMSD = c()
for (conf1 in 1:(length(ls(h))-1)){
  conf1 = as.character(conf1)
  X <- h[[conf1]]

  for (conf2 in (as.integer(conf1)+1):length(ls(h))){
    conf2 = as.character(conf2)
    sum_norm=0
    norm_row = 0
    if (conf1 == conf2){
      next
    }

    Y <- h[[conf2]]
    substr <- X-Y #Subtract the two conformations

    for (row in substr){
      norm_row = 0
      for (i in row){
        norm_row = norm_row + i**2
      }
      sum_norm = sum_norm + norm_row
    }

    cRMSD = sqrt(sum_norm/369) #Calculate the c-RMSD of the two conformations and...
    result_cRMSD = c(result_cRMSD, cRMSD) #...append it to the final vector in which we store all c-RMSD values
  }
}
mean_cRMSD <- mean(result_cRMSD)
median_cRMSD <-  median(result_cRMSD)
range = (max(result_cRMSD) - min(result_cRMSD))/10
i = min(result_cRMSD)
breaks_vector = c(i, i+1*range, i+2*range, i+3*range, i+4*range,
                      i+5*range, i+6*range, i+7*range, i+8*range,
                      i+9*range, i+10*range)
hist(result_cRMSD,
     main="Histogram of all non-optimal c-RMSD",
     xlab="c-RMSD distances",
     col="darkmagenta",
     breaks = breaks_vector)
print(paste("Mean non-optimal c-RMSD:", mean_cRMSD))
print(paste("Median non-optima c-RMSD:", median_cRMSD))

##Calculate optimal c-RMSD  between all 3160 combinations of conformation
opt_result_cRMSD = c()
for (conf1 in 1:(length(ls(h))-1)){
  conf1 = as.character(conf1)
  X <- h[[conf1]]
  Xc <- colMeans(X) #Find the Xc coordinates
  X_xc <- sweep(X, 2, Xc) #Subtract from the conformation the Xc
  for (conf2 in (as.integer(conf1)+1):length(ls(h))){
    conf2 = as.character(conf2)
    sum_norm=0

    Y <- h[[conf2]]
    Yc <- colMeans(Y) #Find the Yc coordinates
    Y_yc <- sweep(Y, 2, Yc) #Subtract from the conformation the Yc
    
    svd_output = svd(t(X_xc)%*%Y_yc) #Calculate the svd
    S <- matrix(diag(svd_output$d), ncol=3)
    U <- svd_output$u
    V <- svd_output$v
    
    #Calculate the Q matrix
    Q <- U%*%t(V) 
    
    if (det(Q)<0){
      U <- sweep(U, 2, c(1,1,-1), `*`)
      Q <- U%*%t(V)
    }
    
    #Calculate the c-RMSD
    for (i in 1:dim(X)[1]){
      row_x <- X[i,]
      Q_x <- row_x%*%Q
      row_y <- Y[i,]
      substr <- Q_x-row_y
      
      norm_row = 0
      for (i in substr){
        norm_row = norm_row + i**2
      }
      sum_norm = sum_norm + norm_row
      
    }
    
    cRMSD = sqrt(sum_norm/dim(X)[1]) #Calculate the c-RMSD of the two conformations and...
    opt_result_cRMSD = c(opt_result_cRMSD, cRMSD) #...append it to the final vector in which we store all c-RMSD values
  }
}
mean_opt_cRMSD <- mean(opt_result_cRMSD)
median_opt_cRMSD <-  median(opt_result_cRMSD)

range = (max(opt_result_cRMSD) - min(opt_result_cRMSD))/10
i = min(opt_result_cRMSD)
breaks_vector_opt = c(i, i+1*range, i+2*range, i+3*range, i+4*range,
                      i+5*range, i+6*range, i+7*range, i+8*range,
                      i+9*range, i+10*range)
hist(opt_result_cRMSD,
     main="Histogram of all optimal c-RMSD distances",
     xlab="c-RMSD distances",
     col="darkmagenta",
     breaks = breaks_vector_opt)
print(paste("Mean c-RMSD:", mean_opt_cRMSD))
print(paste("Median c-RMSD:", median_opt_cRMSD))





#d-RMSD

##Calculate the distances and store them in the d hash table
temp_vector_dist = c()
d <- hash()
for (conf1 in 1:(length(ls(h)))){
  conf1 = as.character(conf1)
  X <- h[[conf1]] #Take each conformation 369x3 table

  for (row1 in 1:dim(X)[1]){ #dim(X)[1] = 369
    if (row1 == dim(X)[1]){break}
    for (row2 in (row1+1):dim(X)[1]){
      temp_vector_dist <- c(temp_vector_dist, dist(rbind(X[row1,], X[row2,])))
    }
  }
  .set(d, keys=conf1, values=temp_vector_dist)
  temp_vector_dist <- c()
  print(conf1)
}
  
  
result_dRMSD <- 0
for (conf1 in 1:(length(ls(d))-1)){
  if (conf1 == length(ls(d))){break}
  for (conf2 in (conf1+1):length(ls(d))){
       conf1 = as.character(conf1)
       X <- d[[conf1]]
       conf2 = as.character(conf2)
       Y <- d[[conf2]]
       dRMSD <- sqrt(sum((X-Y)**2)/length(X))
       result_dRMSD = c(result_dRMSD, dRMSD)
       #print(paste("Conf. ", conf1, " VS Conf. ", conf2, ": ", dRMSD))
    
  }
  print(conf1)
}


mean_dRMSD <- mean(result_dRMSD)
median_dRMSD <-  median(result_dRMSD)
print(paste("Mean d-RMSD:", mean_dRMSD))
print(paste("Median d-RMSD:", median_dRMSD))

range = (max(result_dRMSD) - min(result_dRMSD))/10
i = min(result_dRMSD)
breaks_vector_dRMSD = c(i, i+1*range, i+2*range, i+3*range, i+4*range,
                      i+5*range, i+6*range, i+7*range, i+8*range,
                      i+9*range, i+10*range)

hist(result_dRMSD, breaks = breaks_vector_dRMSD,
     main="Histogram of all d-RMSD distances",
     xlab="d-RMSD distances",
     col="darkmagenta"
)
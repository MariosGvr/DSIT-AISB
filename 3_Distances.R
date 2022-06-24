#3. Distances

# Loading package
library(hash)


datasets_path = "~/MASTER/2nd Semester/Alg. in Structural Biology/GavrielatosMarios_Assignment1/" #Change this if necessary
filename = "covid.txt" #Name of the file in which we have stored the 50 coordinates of the Ca atoms
setwd(datasets_path) #Change working directory

h <- hash() #Initiate hash table in which we are going to store the conforamtions


cov_matrix = read.csv(filename, sep=",", header=F) #Read the csv file...
.set(h, keys=1, values=cov_matrix) #...and add it to the hash table


M_emb <- matrix(nrow=50, ncol=50)
M <- matrix(nrow=50, ncol=50)

for (row1 in 1:dim(M)[1]){
  for (row2 in (row1):dim(M)[1]){
    a = dist(rbind(cov_matrix[row1,], cov_matrix[row2,])) #Find the distance
    
    #Distance matrix with embedding
    M_emb[row1, row2] = (a**2)/2
    M_emb[row2, row1] = (a**2)/2
    
    #Distance matrix without embedding
    M[row1, row2] = a 
    M[row2, row1] = a
    }
}

#Create the Cayley-Menger (border) matrix
cm_matrix_emb <- rbind(rep(1, times=dim(M_emb)[1]), M_emb)
cm_matrix_emb <- cbind(rep(1, times=dim(M_emb)[1]+1), cm_matrix_emb)
cm_matrix_emb[1,1] = 0
rank_B <- qr(cm_matrix_emb)$rank
print(paste("Rank(B_emb):", rank_B))


#cm_matrix <- rbind(rep(1, times=dim(M)[1]), M)
#cm_matrix <- cbind(rep(1, times=dim(cm_matrix)[1]), cm_matrix)
#cm_matrix[1,1] = 0
#rank_B <- qr(cm_matrix)$rank
#print(paste("Rank(B):", rank))

#Create the perpetuated B matrices
##5%
cm_matrix_5 = cm_matrix_emb
for (row in 2:50){
  for (col in (row+1):51){
    #Find the error
    value = cm_matrix_emb[row, col]*0.05
    #Choose randomly if we are going to add or subtract error
    values = c(value, -value) 
    error = sample(values, 1)
    #Create the error symmetrically
    cm_matrix_5[row, col] = cm_matrix_5[row, col] + error
    cm_matrix_5[col, row] = cm_matrix_5[col, row] + error
  }
}
rank_5 = qr(cm_matrix_5)$rank
print(paste("The rank of the 5% pertubated B matrix:", rank_5))

##10%
cm_matrix_10 = cm_matrix_emb
for (row in 2:50){
  for (col in (row+1):51){
    value = cm_matrix_emb[row, col]*0.1
    values = c(value, -value)
    error = sample(values, 1)
    cm_matrix_10[row, col] = cm_matrix_10[row, col] + error
    cm_matrix_10[col, row] = cm_matrix_10[col, row] + error
  }
}
rank_10 = qr(cm_matrix_10)$rank
print(paste("The rank of the 10% pertubated B matrix:", rank_10))


#Find the gram matrices

##5%
gram_matrix_5 = matrix(0, 50, 50)
for (i in 1:50){
  for (j in 1:50){
    i_j = (cm_matrix_5[i+1,2] - cm_matrix_5[i+1,j+1] + cm_matrix_5[j+1,2])
    gram_matrix_5[i,j] = i_j
    gram_matrix_5[j,i] = i_j
  }
}
#cm_matrix_5 <- cm_matrix_5[-1,-1]
#cm_matrix_5 <- sweep(cm_matrix_5, 2, cm_matrix_5[2,])

##10%
gram_matrix_10 = matrix(0, 50, 50)
for (i in 1:50){
  for (j in 1:50){
    i_j = (cm_matrix_10[i+1,2] - cm_matrix_10[i+1,j+1] + cm_matrix_10[j+1,2])
    gram_matrix_10[i,j] = i_j
    gram_matrix_10[j,i] = i_j
  }
}




#Perform SVD
##5%
udv_5 = svd(gram_matrix_5, nu = 3, nv = 3)
S_5 <- udv_5$d[1:3] ##Find the 3 largest singular values and force rank(G)=3
S_5 <- matrix(diag(S_5), ncol=3)
U_5 <- udv_5$u
V_5 <- udv_5$v

p_matrix_05 <- t(sqrt(S_5)%*%t(V_5)) #Get 3D coordinates
.set(h, keys=2, values=p_matrix_05) #Add 3D coordinates to the hash table


##10%
udv_10 = svd(gram_matrix_10, nu = 3, nv = 3)
S_10 <- udv_10$d[1:3] #Find the 3 largest singular values and force rank(G)=3
S_10 <- matrix(diag(S_10), ncol=3)
U_10 <- udv_10$u
V_10 <- udv_10$v

p_matrix_10 <- t(sqrt(S_10)%*%t(V_10)) #Get 3D coordinates
.set(h, keys=3, values=p_matrix_10) #Add 3D coordinates to the hash table

#Find the c-RMSD values
result_cRMSD = c()
for (conf in 2:3){
  X <- h[["1"]]
  
  sum_norm=0
  norm_row = 0

  Y <- h[[as.character(conf)]]
  substr <- X-Y
    
  for (row in substr){
    norm_row = 0
    for (i in row){
      norm_row = norm_row + i**2
    }
    sum_norm = sum_norm + norm_row
  }
    
  cRMSD = sqrt(sum_norm/50)
  result_cRMSD = c(result_cRMSD, cRMSD)
}
print(paste("c-RMSD of 5% perturbation vs original structure:", result_cRMSD[1]))
print(paste("c-RMSD of 10% perturbation vs original structure:", result_cRMSD[2]))



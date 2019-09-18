# this is the code I used to split the liver dataset into test and training. 

print(load("GTEx.liver.geneSymbols.matrix.asinh.RData"))

# for 175 samples, 80% is 140 samples
num.test <- sample(1:175, 35)
mtx.test <- mtx[, num.test]
mtx.train <- mtx[, -num.test]

dim(mtx.test)
dim(mtx.train)

save(mtx.test, file = "GTEx.liver.test.RData")
save(mtx.train, file = "GTEx.liver.train.RData")

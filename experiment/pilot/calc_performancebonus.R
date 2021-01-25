#Perfomance bonus payments pilot 2




set.seed(12345)
randinfs <- round(runif(5)*25) 
selectedinfs <- randinfs+seq(0, 96, 24)

matrix(datlong[which(datlong$inforder %in% selectedinfs),]$acc, nrow=5) # each column is a participants acc on chosen infs.
#1 pound within 5%, .5 within 10, .25 within 15
# pp1: .5, 0, .5, 0, 1
# pp2: 0,0,0,0,0
# pp3: .5, .5, 0, .5 .5 
# pp4: .5, .5, .5, 
# pp5: 5 * .5
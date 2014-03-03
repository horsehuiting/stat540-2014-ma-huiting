p = 0.2 # prob to get head
n<-c(100,1000,10000,100000)
n1 = n[1]
e<-c()
y<-rbinom(n1,1,0.2)
e[1]<-sum(y)/length(y)
p-e[1]

n2 = n[2]
y<-rbinom(n2,1,0.2)
e[2]<-sum(y)/length(y)
p-e[2]

n3 = n[3]
y<-rbinom(n3,1,0.2)
e[3]<-sum(y)/length(y)
p-e[3]

n4 = n[4]
y<-rbinom(n4,1,0.2)
e[4]<-sum(y)/length(y)
p-e[4]

plot(e~n)
abline(0.2,0)
## yes, when the sample size is getting bigger, the fraction of tosses which are heads
## is getting closer to real p

nn<-100 
avg<-c()
MAD<-c()
variance<-c()
iqr<-c()
Y <- matrix(NA,length(n),nn)
for (j in 1:length(n)){
  X<-matrix(NA,nn,n[j])
  for (i in 1:nn){
    X[i,1:n[j]] <- rbinom(n[j],1,p) 
    Y[j,i] <- sum(X[i,1:n[j]])/n[j]
    
  }
}
rownames(Y, do.NULL = FALSE)
rownames(Y) <- c("100","1000","10000","100000")
boxplot(Y,use.cols=FALSE) #p=0.2


##p=0.4

p=0.4
nn<-100 
avg<-c()
MAD<-c()
variance<-c()
iqr<-c()
Y <- matrix(NA,length(n),nn)
for (j in 1:length(n)){
  X<-matrix(NA,nn,n[j])
  for (i in 1:nn){
    X[i,1:n[j]] <- rbinom(n[j],1,p) 
    Y[j,i] <- sum(X[i,1:n[j]])/n[j]
    
  }
}

rownames(Y, do.NULL = FALSE)
rownames(Y) <- c("100","1000","10000","100000")
boxplot(Y,use.cols=FALSE)




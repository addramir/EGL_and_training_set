nlop=1000
n_sample=10000
out=array(NA,c(n_sample,2))
for (i in 1:nlop){
  g=rnorm(n_sample)
  liab=g*runif(1,-1,1)*rnorm(n_sample)
  liab=liab-mean(liab)
  liab=liab/sd(liab)
  y=rep(0,n_sample)
  y[liab>=runif(1,0.001,3)]=1
  
  prev=sum(y==1)/n_sample
  l1=summary(glm(y~g,family = "binomial"))
  
  z=l1$coefficients["g",3]
  varg=var(g)
  se=1/sqrt(n_sample*varg*prev*(1-prev))
  beta=z*se
  
  out[i,1]=beta
  out[i,2]=l1$coefficients["g",1]
}

plot(out[,1],out[,2])
abline(0,1)

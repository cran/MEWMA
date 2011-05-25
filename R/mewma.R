mewma <-
function(x,n,lambda){
h4<-matrix(c(
8.6336,
9.6476,
10.083,
10.3114,
10.4405,
10.5152,
10.5581,
10.5816,
10.5932,
10.814,
11.8961,
12.3505,
12.5845,
12.7143,
12.788,
12.8297,
12.8524,
12.8635,
12.7231,
13.8641,
14.3359,
14.576,
14.7077,
14.7818,
14.8234,
14.846,
14.857,
14.5363,
15.7293,
16.217,
16.4629,
16.5965,
16.6711,
16.7127,
16.7352,
16.7463,
16.2634,
17.5038,
18.0063,
18.2578,
18.3935,
18.4687,
18.5105,
18.5331,
18.5442,
17.9269,
19.2113,
19.7276,
19.9845,
20.1223,
20.1982,
20.2403,
20.2631,
20.2743,
19.541,
20.8665,
21.396,
21.6581,
21.798,
21.8747,
21.9171,
21.9401,
21.9515,
21.1152,
22.4796,
23.0217,
23.2887,
23.4307,
23.5082,
23.551,
23.5742,
23.5858,
22.6565,
24.0579,
24.6119,
24.8838,
25.0278,
25.1062,
25.1493,
25.1728,
25.1846),nrow=9)



require(MASS)


if(class(x)!="matrix")(cat("x must be a matrix           "))

a<-dim(x)
x1<-matrix(0,a[1]/n,a[2])
z<-matrix(0,a[1]/n,a[2])
y<-matrix(0,a[1]/n,a[2])
t2<-matrix(0,a[1]/n,1)
t3<-matrix(0,a[1]/n,1)

for(j in 1:a[2]){
  for(i in 1:(a[1]/n)){
   y[i,j]<-mean(x[(i*n-(n-1)):(i*n),j])
 }
}

for(i in 1:(a[1]/n)){
  for(j in 1:a[2]){
   x1[i,j]<-y[i,j]-mean(y[,j])
   ifelse(i==1,z[i,j]<-lambda*x1[i,j],z[i,j]<-lambda*x1[i,j]+(1-lambda)*z[i-1,j])     
  }
}
cova<-covariance(x,n)

for(i in 1:(a[1]/n)){
weights<-cova*(lambda*(1-((1-lambda)^(2*i)))/(2-lambda))
inv<-ginv(weights)
za<-matrix(z[i,])

t2[i,]<-t(za)%*%inv%*%za

  } 
 
rownames(h4)<-c(seq(0.1,0.9, by =0.1))
colnames(h4)<-c(1:9)

m1<-rownames(h4)
m2<-colnames(h4)
l<-lambda*10
ucl<-h4[m1[l],m2[a[2]-1]]
 
ifelse(max(t2)>ucl,s3<-max(t2),s3<-ucl )

plot(t2,ylim=c(0,ceiling(s3)),type="o") 
s1<-c(ucl,ucl)
s2<-c(0,0)
abline(s1,s2)

cat("Multivariate Exponentially Weighted Moving Average (MEWMA) Control Chart    ")
cat("Upper Control Limits(UCL)" )
print(ucl)

k<-1
for(i in 1:(a[1]/n)){
 if(t2[i]>ucl) {
 t3[k,1]<-i
   k<- k+1 }
 } 
 
 if(k>1){
 cat("The following(s) point(s) fall outside of the control limits" )
  for(i in 1:(a[1]/n)){
    if(t3[i]!=0)(print(t3[i]))
     }
 }

outList = list ("Multivariate Exponentially Weighted Moving Average (MEWMA) Control Chart",t2=t2,covariance=cova)
    invisible(outList)

}


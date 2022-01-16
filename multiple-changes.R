### Simulations. Generate data with THREE transient changes. 
### Study FAR, FRR, the distribution of estimates.
### Thresholds are already calculated, we just download our results.

n = 1000; alpha = 0.05; 
Nruns = 20000; NrunsH = 200000; 
mu0 = 0; sigma0 = 1; sigma1 = 1;
a = floor(c(0.15,0.45,0.75)*n); b = a+0.1*n;		
MU = c(seq(0.1,1.5,0.1),2);  Nmu = length(MU); H = rep(0,Nmu); 

# MU = means, H = thresholds for them

### Thresholds
# Input = read.csv("C:\\Users\\baron\\Documents\\Research\\Malov\\simulation.results.mean.exact.csv")
# MU = Input$Mu1[1:11]; H = Input$Threshold[1:11];

### Function estimating the threshold h as the quantile from the null dist of Wmax

threshold = function(n,mu,alpha,NrunsH){
  # n = sample size; mu = disturbed distribution mean. 
  # We are actually calculating the null distribution of the test statistic Wn
  Wmax = rep(0,NrunsH); 
  for (run in 1:NrunsH){
   	 X = rnorm(n,0,1);      # Generate a sequence from the base distr.
   	 z = log( dnorm(X,mu,1)/dnorm(X,0,1) );	    # log-likelihood ratios
   	 S = cumsum(z);   							    # Partial sums of z
   	 W = S-cummin(S);				  		        # Cusum process
     Wmax[run] = max(W); 
  } # end for-loop by run
  h = quantile(Wmax,1-alpha);
  return(h);
}  # End of function "threshold"


# threshold = function(n,mu,alpha,NrunsH)
# Example:
# threshold(100, 0.3, 0.05, 20000)



### Define a function for the unknown number of transient change-point estimation
### for known Normal distributions

MLE.multiple = function(X,m0,s0,m1,s1,h){  # Estimation of a single transient change-point
  # X = data sequence, m and s are means and sd of pre- and post-change distributions
  n = length(X);  
  z = log( dnorm(X,m1,s1)/dnorm(X,m0,s0) );	# log-likelihood ratios

  t = 0;    # t = last estimated change-point
  Mode = 0; # Under Mode = 0, seek a change from F to G
            # Under Mode = 1, seek a change from G to F
  a.hat = rep(0,0); b.hat = rep(0,0);

  while (t < n){ 
    S = cumsum(z[(t+1):n]);                     # LLR random walk, after t
    if (Mode == 0){
       W = S - cummin(S);               			# Cusum process
       t.exceed = which(W >= h); 
       if (length(t.exceed)>0){ 
          tau = min(t.exceed);                  # Detection time
          KerW = which(W[1:tau]==0);            # Kernel, zeros of W
          if (length(KerW)>0){
              a.new = max(KerW);                # MLE
              } else { a.new=0; } # end of "if KerW non-empty"
          a.hat = c(a.hat, t+a.new);
          Mode = 1;
          t = t+a.new;
          } else { t = n;                       # what if tau is empty
       } # end of "if tau is not empty"
       } else {                                 # what if Mode==1
       S = -S;
       W = S - cummin(S);               			# Cusum process
       t.exceed = which(W >= h); 
       if (length(t.exceed)>0){ 
          tau = min(t.exceed);                  # Detection time
          KerW = which(W[1:tau]==0);            # Kernel, zeros of W
          if (length(KerW)>0){
              b.new = max(KerW);                # MLE
              } else { b.new=0; } # end of "if KerW non-empty"
          b.hat = c(b.hat, t+b.new);
          Mode = 0;
          t = t+b.new;
          } else { t = n; b.hat = c(b.hat, n);  # what if tau is empty
       } # end of "if tau is not empty"
    } # end of "if Mode==0"
  } # end of while(t<n) loop
  Na = length(a.hat); Nb = length(b.hat);
  return(data.frame(a.hat,b.hat));               
}   # end of function MLE.single


### MC study ###

Nmu = length(MU);                  
FAR = rep(0,Nmu); FRR = FAR; 
Detected = matrix(rep(0,Nmu*7),Nmu,7); 
a1.hat.mean = FAR; b1.hat.mean = FAR; a1.hat.sd = FAR; b1.hat.sd = FAR;
a2.hat.mean = FAR; b2.hat.mean = FAR; a2.hat.sd = FAR; b2.hat.sd = FAR;
a3.hat.mean = FAR; b3.hat.mean = FAR; a3.hat.sd = FAR; b3.hat.sd = FAR;

InstReg = c(rep(0,a[1]), rep(1,b[1]-a[1]),
rep(0,a[2]-b[1]), rep(1,b[2]-a[2]), rep(0,a[3]-b[2]), 
rep(1,b[3]-a[3]), rep(0,n-b[3]));      # Actual instability regions


for (k in 1:Nmu){ mu1 = MU[k]; 
   OUTPUT = data.frame(Nmu,k);
   write.table(OUTPUT,file="C:\\Users\\baron\\Documents\\Research\\Malov\\output.txt")

   h = threshold(n,mu1,alpha,NrunsH); H[k] = h;
   
   a1=rep(0,Nruns); a2=a1; a3=a1; b1=a1; b2=a1; b3=a1;
   Ndetected=a1;
   false.alarm = a1; false.readjustment = a1;

   for (run in 1:Nruns){  
     X = c( rnorm(a[1],mu0,sigma0), rnorm(b[1]-a[1],mu1,sigma1), 
            rnorm(a[2]-b[1],mu0,sigma0), rnorm(b[2]-a[2],mu1,sigma1), 
            rnorm(a[3]-b[2],mu0,sigma0), rnorm(b[3]-a[3],mu1,sigma1), 
            rnorm(n-b[3],mu0,sigma0) );

     MLE.result = MLE.multiple(X,mu0,sigma0,mu1,sigma1,h);

     a.hat = MLE.result$a.hat;  b.hat = MLE.result$b.hat;  
     Na = length(a.hat);
     Ndetected[run] = Na;
     if (Na==3){
        a1[run]=a.hat[1]; b1[run]=b.hat[1];
        a2[run]=a.hat[2]; b2[run]=b.hat[2];
        a3[run]=a.hat[3]; b3[run]=b.hat[3];
               }

     # False alarm and false readjustment
     if (Na > 1){
        for (chp in 1:(Na-1)){
            GoodDetect = sum( InstReg[(a.hat[chp]+1) : b.hat[chp]] ); # false alarms
              if (GoodDetect == 0){ false.alarm[run]=1; }
            GoodNoDetect = sum( 1-InstReg[(b.hat[chp]+1) : a.hat[chp+1]] ); # false readjustment    
              if (GoodNoDetect == 0){ false.readjustment[run]=1; }
           } # end of the chp-loop
        # The last re-adjustment period is up to n
            GoodDetect = sum( InstReg[(a.hat[Na]+1) : b.hat[Na]] ); # false alarms
              if (GoodDetect == 0){ false.alarm[run]=1; }
            if (b.hat[Na] < n){
            GoodNoDetect = sum( 1-InstReg[(b.hat[Na]+1) : n] ); # false readjustment    
              if (GoodNoDetect == 0){ false.readjustment[run]=1; }
            } # if the last b.hat < n
     } # End of "if Na > 0"

     if (Na == 1){
            GoodDetect = sum( InstReg[(a.hat[1]+1) : b.hat[1]] ); # false alarms
              if (GoodDetect == 0){ false.alarm[run]=1; }
            if (b.hat[1] < n){
            GoodNoDetect = sum( 1-InstReg[(b.hat[1]+1) : n] ); # false readjustment    
              if (GoodNoDetect == 0){ false.readjustment[run]=1; }
            } # if the last b.hat < n
     } # End of "if Na == 0"

   }  # End of MC runs for the given parameters
   
   FAR[k] = mean(false.alarm); 
   FRR[k] = mean(false.readjustment);

   for (Ndet in 1:6){ Detected[k,Ndet] = sum(Ndetected==(Ndet-1)); }
      Detected[k,7] = sum(Ndetected>=6);

   Ndet3 = which(Ndetected==3);
   if (length(Ndet3)>0){ 
     a1.hat.mean[k] = mean(a1[Ndet3]); 
     a1.hat.sd[k]   = sd(  a1[Ndet3]); 
     a2.hat.mean[k] = mean(a2[Ndet3]); 
     a2.hat.sd[k]   = sd(  a2[Ndet3]); 
     a3.hat.mean[k] = mean(a3[Ndet3]); 
     a3.hat.sd[k]   = sd(  a3[Ndet3]); }
   Ndet3 = which(Ndetected==3 & b3 < n);
   if (length(Ndet3)>0){ 
     b1.hat.mean[k] = mean(b1[Ndet3]); 
     b1.hat.sd[k]   = sd(  b1[Ndet3]); 
     b2.hat.mean[k] = mean(b2[Ndet3]); 
     b2.hat.sd[k]   = sd(  b2[Ndet3]); 
     b3.hat.mean[k] = mean(b3[Ndet3]); 
     b3.hat.sd[k]   = sd(  b3[Ndet3]); }  
} # End of the k-loop, over mu1

ProbDetected0 = Detected[,1]/Nruns;
ProbDetected1 = Detected[,2]/Nruns;
ProbDetected2 = Detected[,3]/Nruns;
ProbDetected3 = Detected[,4]/Nruns;
ProbDetected4 = Detected[,5]/Nruns;
ProbDetected5 = Detected[,6]/Nruns;
ProbDetected6Plus = Detected[,7]/Nruns;

Output = data.frame(MU, H, FAR, FRR, 
ProbDetected0, ProbDetected1, ProbDetected2, ProbDetected3, 
ProbDetected4, ProbDetected5, ProbDetected6Plus, 
a1.hat.mean, b1.hat.mean, a2.hat.mean, b2.hat.mean, a3.hat.mean, b3.hat.mean,
a1.hat.sd, b1.hat.sd, a2.hat.sd, b2.hat.sd, a3.hat.sd, b3.hat.sd
);

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\simulation.multiple.csv")

print(Output);


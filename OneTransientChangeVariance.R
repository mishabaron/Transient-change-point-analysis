### Simulations. Generate data with one transient change in VARIANCE.
### Estimate threshold h WITHOUT DOOB.
### Normal(0,1) changes to Normal(0,sigma), known parameters

n = 1000; alpha = 0.05; 
Nruns = 50000; NrunsH = 200000;  # MC runs for the threshold h and for all simulations 
mu0 = 0; mu1 = 0; sigma1 = 1;
SIGMA = c(0.5,0.75,0.9,0.95,1.05,1.1,1.25,1.5,2);
a = floor(0.5*n); b = floor(0.7*n);		# Transient change period

### Function estimating the threshold h as the quantile from the null dist of Wmax

threshold = function(n,sigma,alpha,NrunsH){
  # n = sample size; sigma = disturbed distribution std. 
  Wmax = rep(0,NrunsH); Wn = rep(0,NrunsH);
  for (run in 1:NrunsH){
   	 X = rnorm(n,0,1);      # Generate a sequence from the base distr.
   	 z = log( dnorm(X,0,sigma)/dnorm(X,0,1) );	# log-likelihood ratios
   	 S = cumsum(z);   							    # Partial sums of z
   	 W = S-cummin(S);				  		        # Cusum process
     Wmax[run] = max(W); 
     Wn[run] = W[n];
  } # end for-loop by run
  # hist(Wmax,30); 
  h = quantile(Wmax,1-alpha);      
  return(h);
  }  # End of function "threshold"


### Define a function for a single transient change-point estimation
### for known Normal distributions

MLE.single = function(X,sigma,h){  # Estimation of a single transient change-point
  # X = data sequence, mu is the post-change mean. N(0,1) -> N(mu,1)
  n = length(X); m1 = 0; m0 = 0; s1 = sigma; s0 = 1;
  z = log( dnorm(X,m1,s1)/dnorm(X,m0,s0) );	# log-likelihood ratios
  S = cumsum(z);                            # LLR random walk
  W = S - cummin(S);               			# Cusum process
  b.hat = which.max(W); a.hat = which.min(S[1:b.hat]);  # MLE
  detection = 1*( W[b.hat] >= h ); # If h is exceeded, ch.pt. is detected
  return(c(a.hat,b.hat,detection));
}   # end of function MLE.single


### MC study with different means, N(0,1) -> N(0,sigma)

Nscena = length(SIGMA);                    # Number of scenarious
H = rep(0,Nscena); prob.false.alarm = H; prob.detection = H; 
a.hat.mean = H; b.hat.mean = H; a.hat.sd = H; b.hat.sd = H;
Threshold.Doob = H;

for (k in 1:Nscena){ sigma1 = SIGMA[k];
   OUTPUT = data.frame(Nscena,k);
   write.table(OUTPUT,file="C:\\Users\\baron\\Documents\\Research\\Malov\\output.txt")

   h = threshold(n,sigma1,alpha,NrunsH);  # Estimate the threshold h = - log( alpha / Eexp(W_n) )
   false.alarm.run = rep(0,Nruns);
   detection.run = rep(0,Nruns);
   a.hat = rep(0,Nruns);
   b.hat = rep(0,Nruns);
   for (run in 1:Nruns){
     # The case of no change
     X = rnorm(n,mu0,sigma0); 
     MLE.result = MLE.single(X,sigma1,h);     # Result = c(a.hat,b.hat,detection)  
     false.alarm.run[run] = MLE.result[3];
     
     # The case of a transient change
     X = c( rnorm(a,mu0,sigma0), rnorm(b-a,mu1,sigma1), rnorm(n-b,mu0,sigma0) );
     MLE.result = MLE.single(X,sigma1,h);
     a.hat[run] = MLE.result[1];
     b.hat[run] = MLE.result[2];
     detection.run[run] = MLE.result[3];
                    }  # End of MC runs for the given parameters

   H[k] = h; 
   prob.false.alarm[k] = mean(false.alarm.run); # False alarm rate | no change
   prob.detection[k] = mean(detection.run);     # Probability of detection | change
   a.hat.mean[k] = mean(a.hat);
   a.hat.sd[k] = sd(a.hat);
   b.hat.mean[k] = mean(b.hat);
   b.hat.sd[k] = sd(b.hat);
} # End of the k-loop, over scenarios

Threshold = H;

Mu0=rep(0,Nscena); Mu1=Mu0; Sigma0=rep(sigma0,Nscena); Sigma1=SIGMA; 
Params = data.frame(Mu0,Mu1,Sigma0,Sigma1);

Output = data.frame(Params, Threshold, 
         prob.false.alarm, prob.detection,
         a.hat.mean, a.hat.sd, b.hat.mean, b.hat.sd);

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\simulation.variance.csv")

print(Output);

### Simulations. Generate data with one transient change.
### Estimate threshold h. Detect the change-point. 
### Study FAR, FRR, the distribution of estimates.

### One example: Normal(0,1) changes to Normal(0.3,1), known parameters

n = 100; alpha = 0.05; 
Nruns = 1000; Nruns.h = 50000;  # MC runs for the threshold h and for all simulations 
mu0 = 0; sigma0 = 1; sigma1 = 1;
a = floor(0.5*n); b = floor(0.7*n);		# Transient change period
plot = 0;     # Print and plot the last results, if wanted, let plot=1


### Function estimating the threshold h = - log( alpha / Eexp(W_n) )

threshold = function(n,m0,s0,m1,s1,alpha,Nruns.h){
  # n = sample size; m0,s0 = base parameters; m1,s1 = disturbed parameters 
  Wn = rep(0,Nruns.h);        # Initialize CUSUM
  for (run in 1:Nruns.h){
     X = rnorm(n,m0,s0);      # Generate a sequence from the base distr.
     z = log( dnorm(X,m1,s1)/dnorm(X,m0,s0) );	# log-likelihood ratios
     S = cumsum(z);   							    # Partial sums of z
     Wn[run] = S[n]-min(S);  		      # Cusum process
  } # end for-loop by run
  h = log( mean(exp(Wn)) / alpha );  # Threshold
  return(h);
  }  # End of function "threshold"

# threshold = function(n,m0,s0,m1,s1,alpha,Nruns.h)
# Example:
# threshold(1000, 0, 1, 0.3, 1, 0.05, 20000)


### Define a function for a single transient change-point estimation
### for known Normal distributions

MLE.single = function(X,m0,s0,m1,s1,h){  # Estimation of a single transient change-point
  # X = data sequence, m and s are means and sd of pre- and post-change distributions
  n = length(X);
  z = log( dnorm(X,m1,s1)/dnorm(X,m0,s0) );	# log-likelihood ratios
  S = cumsum(z);                            # LLR random walk
  W = S - cummin(S);               			# Cusum process
  b.hat = which.max(W); a.hat = which.min(S[1:b.hat]);  # MLE
  detection = 1*( W[b.hat] >= h ); # If h is exceeded, ch.pt. is detected
  return(c(a.hat,b.hat,detection));
}   # end of function MLE.single


### MC study with different means, N(0,1) -> N(mu1,1)

MU = c( seq(0.2,1,0.2),seq(1.5,8,0.5)); 
Nscena = length(MU);                    # Number of scenarious
H = rep(0,Nscena); prob.false.alarm = H; prob.detection = H; 
a.hat.mean = H; b.hat.mean = H; a.hat.sd = H; b.hat.sd = H;

for (k in 1:Nscena){ mu1 = MU[k];
   OUTPUT = data.frame(Nscena,k);
   write.table(OUTPUT,file="C:\\Users\\baron\\Documents\\Research\\Malov\\output.txt")

   h = threshold(n,mu0,sigma0,mu1,sigma1,alpha,Nruns.h);  # Estimate the threshold h = - log( alpha / Eexp(W_n) )
   false.alarm.run = rep(0,Nruns);
   detection.run = rep(0,Nruns);
   a.hat = rep(0,Nruns);
   b.hat = rep(0,Nruns);
   for (run in 1:Nruns){
     # The case of no change
     X = rnorm(n,mu0,sigma0); 
     MLE.result = MLE.single(X,mu0,sigma0,mu1,sigma1,h); # Result = c(a.hat,b.hat,detection)  
     false.alarm.run[run] = MLE.result[3];
     
     # The case of a transient change
     X = c( rnorm(a,mu0,sigma0), rnorm(b-a,mu1,sigma1), rnorm(n-b,mu0,sigma0) );
     MLE.result = MLE.single(X,mu0,sigma0,mu1,sigma1,h);
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
Threshold.lower.bound = MU*qnorm((1-alpha)^(1/n)) - MU^2/2;
Mu0=rep(0,Nscena); Mu1=MU; Sigma0=rep(sigma0,Nscena); Sigma1=rep(sigma1,Nscena); 
Params = data.frame(Mu0,Mu1,Sigma0,Sigma1);
Violates = 1*(Threshold < Threshold.lower.bound);

Output = data.frame(Params, Threshold, Threshold.lower.bound, Violates,
         prob.false.alarm, prob.detection,
         a.hat.mean, a.hat.sd, b.hat.mean, b.hat.sd);

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\simulation.results.csv")

print(Output);
sum(Violates);

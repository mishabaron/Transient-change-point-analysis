### Simulations. Generate data with one transient change.
### Estimate threshold h WITHOUT DOOB.
### Study FAR, FRR, the distribution of estimates.
### Normal(0,1) changes to Normal(mu1,1), known parameters

n = 1000; alpha = 0.05; 
Nruns = 50000; NrunsH = 200000;  # MC runs for the threshold h and for all simulations 
mu0 = 0; sigma0 = 1; sigma1 = 1;
MU = c( seq(0.05,0.4,0.05),seq(0.6,1,0.2)); 
a = floor(0.5*n); b = floor(0.7*n);		# Transient change period

### Function estimating the threshold h as the quantile from the null dist of Wmax

threshold = function(n,mu,alpha,NrunsH){
  # n = sample size; mu = disturbed distribution mean. 
  # We are actually calculating the null distribution of the test statistic Wn
  Wmax = rep(0,NrunsH); Wn = rep(0,NrunsH);
  for (run in 1:NrunsH){
   	 X = rnorm(n,0,1);      # Generate a sequence from the base distr.
   	 z = log( dnorm(X,mu,1)/dnorm(X,0,1) );	# log-likelihood ratios
   	 S = cumsum(z);   							    # Partial sums of z
   	 W = S-cummin(S);				  		        # Cusum process
     Wmax[run] = max(W); 
     Wn[run] = W[n];
  } # end for-loop by run
  # hist(Wmax,30); 
  h = quantile(Wmax,1-alpha);
  h.Doob = log(mean(exp(Wn))/alpha);  
  # h = max(quantile(Wmax,1-alpha),-log(alpha));           
  return(c(h,h.Doob));
  }  # End of function "threshold"


# threshold = function(n,mu,alpha,NrunsH)
# Example:
# threshold(100, 0.3, 0.05, 20000)


### Define a function for a single transient change-point estimation
### for known Normal distributions

MLE.single = function(X,mu,h){  # Estimation of a single transient change-point
  # X = data sequence, mu is the post-change mean. N(0,1) -> N(mu,1)
  n = length(X); m1 = mu; m0 = 0; s1 = 1; s0 = 1;
  z = log( dnorm(X,m1,s1)/dnorm(X,m0,s0) );	# log-likelihood ratios
  S = cumsum(z);                            # LLR random walk
  W = S - cummin(S);               			# Cusum process
  b.hat = which.max(W); a.hat = which.min(S[1:b.hat]);  # MLE
  detection = 1*( W[b.hat] >= h ); # If h is exceeded, ch.pt. is detected
  return(c(a.hat,b.hat,detection));
}   # end of function MLE.single


### MC study with different means, N(0,1) -> N(mu1,1)

Nscena = length(MU);                    # Number of scenarious
H = rep(0,Nscena); prob.false.alarm = H; prob.detection = H; 
a.hat.mean = H; b.hat.mean = H; a.hat.sd = H; b.hat.sd = H;
Threshold.Doob = H;

for (k in 1:Nscena){ mu1 = MU[k];
   OUTPUT = data.frame(Nscena,k);
   write.table(OUTPUT,file="C:\\Users\\baron\\Documents\\Research\\Malov\\output.txt")

   hresult = threshold(n,mu1,alpha,NrunsH);  # Estimate the threshold h = - log( alpha / Eexp(W_n) )
   h = hresult[1];
   false.alarm.run = rep(0,Nruns);
   detection.run = rep(0,Nruns);
   a.hat = rep(0,Nruns);
   b.hat = rep(0,Nruns);
   for (run in 1:Nruns){
     # The case of no change
     X = rnorm(n,mu0,sigma0); 
     MLE.result = MLE.single(X,mu1,h);     # Result = c(a.hat,b.hat,detection)  
     false.alarm.run[run] = MLE.result[3];
     
     # The case of a transient change
     X = c( rnorm(a,mu0,sigma0), rnorm(b-a,mu1,sigma1), rnorm(n-b,mu0,sigma0) );
     MLE.result = MLE.single(X,mu1,h);
     a.hat[run] = MLE.result[1];
     b.hat[run] = MLE.result[2];
     detection.run[run] = MLE.result[3];
                    }  # End of MC runs for the given parameters

   H[k] = h; Threshold.Doob[k] = hresult[2];
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

Output = data.frame(Params, Threshold, Threshold.lower.bound, Threshold.Doob,
         Violates,
         prob.false.alarm, prob.detection,
         a.hat.mean, a.hat.sd, b.hat.mean, b.hat.sd);

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\simulation.results.csv")

plot(MU,Threshold,col="blue",type="l",ylim=c(0,max(Threshold)));
lines(MU,Threshold.lower.bound,col="red");
#lines(MU,Threshold.Doob,col="black");
text(8,8,"Exact threshold",col="blue")
text(8,7,"Lower bound",col="red")
#text(8,6,"Doob threshold",col="black")



print(Output);
sum(Violates);

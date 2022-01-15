### Simulations. Generate data with one transient change.
### Estimate threshold h WITHOUT DOOB.
### Study the POWER for different mu1 and (b-a)
### Normal(0,1) changes to Normal(mu1,1), known parameters

n = 1000; alpha = 0.05; 
Nruns = 10000; NrunsH = 100000;  # MC runs for the threshold h and for all simulations 
mu0 = 0; sigma0 = 1; sigma1 = 1;
MU = seq(0.1,1,0.1); 
Nmu = length(MU);
change.interval = c(seq(0.05,0.5,0.05))*n;
a = floor(0.3*n); b = a + change.interval;		# Transient change period
Nb = length(b);


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
  h = quantile(Wmax,1-alpha);         
  return(h);
  }  # End of function "threshold"


### Define a function for a single transient change-point estimation

MLE.single = function(X,mu,h){  # Estimation of a single transient change-point
  # X = data sequence, mu is the post-change mean. N(0,1) -> N(mu,1)
  n = length(X); m1 = mu; m0 = 0; s1 = 1; s0 = 1;
  z = log( dnorm(X,m1,s1)/dnorm(X,m0,s0) );	# log-likelihood ratios
  S = cumsum(z);                            # LLR random walk
  W = S - cummin(S);               			# Cusum process
  detection = 1*( max(W) >= h ); 			# If h is exceeded, ch.pt. is detected
  return(detection);
}   # end of function MLE.single


### MC study with different means, N(0,1) -> N(mu1,1)

Prob.Detection = matrix(rep(0,Nmu*Nb),Nmu,Nb); 

for (k in 1:Nb){
  for (scena in 1:Nmu){ mu1 = MU[scena];
     OUTPUT = data.frame(Nmu,Nb,scena,k);
     write.table(OUTPUT,file="C:\\Users\\baron\\Documents\\Research\\Malov\\output.txt")

     h = threshold(n,mu1,alpha,NrunsH);  # Estimate the threshold h = - log( alpha / Eexp(W_n) )
     detection.run = rep(0,Nruns);
     for (run in 1:Nruns){
       X = c( rnorm(a,mu0,sigma0), rnorm(b[k]-a,mu1,sigma1), rnorm(n-b[k],mu0,sigma0) );
     MLE.result = MLE.single(X,mu1,h);
     detection.run[run] = MLE.result;
                    }  # End of MC runs for the given parameters
  Prob.Detection[scena,k] = mean(detection.run);     # Probability of detection | change
  } # End of the k-loop
} # End of scena-loop, over scenarios

Output = data.frame(MU, Prob.Detection);

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\power-analysis-mean.csv")

print(Output);
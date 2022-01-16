### Simulations. Generate data with one transient change.
### Estimate threshold h WITHOUT DOOB.
### Study the POWER for different sigma1 and (b-a)
### Normal(0,1) changes to Normal(0,sigma1), known parameters

n = 1000; alpha = 0.05; 
Nruns = 20000; NrunsH = 100000;    # MC runs for the threshold h and for all simulations 
mu0 = 0; mu1 = 0; sigma0 = 1; 
SIGMA = c(0.5,0.75,0.9,0.95,1.05,1.1,1.25,1.5,2);
Ns = length(SIGMA);
change.interval = c(seq(0.05,0.5,0.05))*n;
a = floor(0.3*n); b = a + change.interval;		# Transient change period
Nb = length(b);


### Function estimating the threshold h as the quantile from the null dist of Wmax

threshold = function(n,sigma,alpha,NrunsH){
  # n = sample size; mu = disturbed distribution mean. 
  # We are actually calculating the null distribution of the test statistic Wn
  Wmax = rep(0,NrunsH); 
  for (run in 1:NrunsH){
   	 X = rnorm(n,0,1);      # Generate a sequence from the base distr.
   	 z = log( dnorm(X,0,sigma)/dnorm(X,0,1) );	    # log-likelihood ratios
   	 S = cumsum(z);   							    # Partial sums of z
   	 W = S-cummin(S);				  		        # Cusum process
     Wmax[run] = max(W); 
  } # end for-loop by run
  h = quantile(Wmax,1-alpha);         
  return(h);
}  # End of function "threshold"


### Define a function for a single transient change-point estimation

MLE.single = function(X,sigma,h){  # Estimation of a single transient change-point
  # X = data sequence, mu is the post-change mean. N(0,1) -> N(0,s)
  n = length(X); m1 = 0; m0 = 0; s1 = sigma; s0 = 1;
  z = log( dnorm(X,m1,s1)/dnorm(X,m0,s0) );	# log-likelihood ratios
  S = cumsum(z);                               # LLR random walk
  W = S - cummin(S);               			    # Cusum process
  detection = 1*( max(W) >= h ); 			    # If h is exceeded, ch.pt. is detected
  return(detection);
}   # end of function MLE.single


### MC study with different means, N(0,1) -> N(0,s)

Prob.Detection = matrix(rep(0,Ns*Nb),Ns,Nb); 

for (k in 1:Nb){
  for (scena in 1:Ns){ sigma1 = SIGMA[scena];
     OUTPUT = data.frame(Ns,Nb,scena,k);
     write.table(OUTPUT,file="C:\\Users\\baron\\Documents\\Research\\Malov\\output.txt")

     h = threshold(n,sigma1,alpha,NrunsH);  
     detection.run = rep(0,Nruns);
     for (run in 1:Nruns){
       X = c( rnorm(a,0,1), rnorm(b[k]-a,0,sigma1), rnorm(n-b[k],0,1) );
     MLE.result = MLE.single(X,sigma1,h); 
     detection.run[run] = MLE.result; 
                    }  # End of MC runs for the given parameters
  Prob.Detection[scena,k] = mean(detection.run);     # Probability of detection | change
  } # End of the k-loop
} # End of scena-loop, over scenarios

Output = data.frame(SIGMA, Prob.Detection);

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\power-analysis-variance.csv")

print(Output);
### Simulations. Generate data with one transient change.
### Estimate threshold h WITHOUT DOOB.
### Study FAR, FRR, the distribution of estimates.
### Normal(0,1) changes to Normal(mu1,1), known parameters

n = 1000; alpha = 0.05; 
NrunsH = 100000;  # MC runs for the threshold h and for all simulations 
mu0 = 0; sigma0 = 1; sigma1 = 1;
MU = c( seq(0.02,0.1,0.02),seq(0.15,1,0.05),seq(1,2,0.2),seq(2.5,7.5,0.5) ); 
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
  h = quantile(Wmax,1-alpha);         
  return(h);
  }  # End of function "threshold"


Nscena = length(MU); H = rep(0,Nscena);

for (k in 1:Nscena){ mu = MU[k];
   OUTPUT = data.frame(Nscena,k);
   write.table(OUTPUT,file="C:\\Users\\baron\\Documents\\Research\\Malov\\output.txt")

   H[k] = threshold(n,mu,alpha,NrunsH);  
} # End of the k-loop, over scenarios

Threshold = H;
MU1 = seq(0,7.5,0.05);
Threshold.lower.bound = MU1*qnorm((1-alpha)^(1/n)) - MU1^2/2;

plot(MU,Threshold,col="blue",type="l",lwd=2,ylim=c(0,max(Threshold)),
  xlab="Magnitude of change");
lines(MU1,Threshold.lower.bound,lwd=2,col="red");
text(6.9,7.8,"Exact threshold",col="blue")
text(6.9,7,"Lower bound",col="red")


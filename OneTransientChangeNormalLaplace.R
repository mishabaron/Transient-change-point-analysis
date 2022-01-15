### Simulations. Change from N(0,1) to Laplace(0,1/sqrt(2))

n = 1000; alpha = 0.05; 
Nruns = 20000; NrunsH = 100000;  # MC runs for the threshold h and for all simulations 
a = floor(0.5*n); b = floor(0.7*n);		# Transient change period

### Function estimating the threshold h as the quantile from the null dist of Wmax

threshold = function(n,NrunsH){
  Wmax = rep(0,NrunsH); Wn = rep(0,NrunsH);
  for (run in 1:NrunsH){
   	 X = rnorm(n,0,1);   
   	 z = log( sqrt(2)/2*exp(-sqrt(2)*abs(X)) / dnorm(X,0,1) );	
   	 S = cumsum(z);   							    # Partial sums of z
   	 W = S-cummin(S);				  		        # Cusum process
     Wmax[run] = max(W); 
  } # end for-loop by run
  h = quantile(Wmax,1-alpha);      
  return(h);
  }  # End of function "threshold"


### Define a function for a single transient change-point estimation

MLE.single = function(X,h){  
  n = length(X); 
  z = log( sqrt(2)/2*exp(-sqrt(2)*abs(X)) / dnorm(X,0,1) );	
  S = cumsum(z);                            # LLR random walk
  W = S - cummin(S);               			# Cusum process
  b.hat = which.max(W); a.hat = which.min(S[1:b.hat]);  # MLE
  detection = 1*( W[b.hat] >= h ); # If h is exceeded, ch.pt. is detected
  return(c(a.hat,b.hat,detection));
}   # end of function MLE.single


### MC study with N(0,1) -> Laplace(0,1/sqrt(2))

   h = threshold(n,NrunsH);  
   detection.run = rep(0,Nruns);
   a.hat = rep(0,Nruns);
   b.hat = rep(0,Nruns);
   for (run in 1:Nruns){
     # Generate the middle sample from the Double Exponential distribution

     Bern = rbinom(b-a,1,0.5)*2-1;      # Sign
     Y = rexp(b-a,rate=sqrt(2));

     X = c( rnorm(a,0,1), Bern*Y, rnorm(n-b,0,1) );

     MLE.result = MLE.single(X,h);
     a.hat[run] = MLE.result[1];
     b.hat[run] = MLE.result[2];
     detection.run[run] = MLE.result[3];
                    }  # End of MC runs for the given parameters

   H = h; 
   prob.detection = mean(detection.run);     # Probability of detection | change
   a.hat.mean = mean(a.hat);
   a.hat.sd = sd(a.hat);
   b.hat.mean = mean(b.hat);
   b.hat.sd = sd(b.hat);


Threshold = H;

Output = data.frame(Threshold, prob.detection,
         a.hat.mean, a.hat.sd, b.hat.mean, b.hat.sd);

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\simulation.laplace.csv")

print(Output);

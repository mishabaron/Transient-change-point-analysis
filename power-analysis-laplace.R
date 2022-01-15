### Simulations. Generate data with one transient change.
### Estimate threshold h WITHOUT DOOB.
### Study the POWER for different sigma1 and (b-a)
### Normal(0,1) changes to Laplace(0,1/sqrt(2))

n = 1000; alpha = 0.05; 
Nruns = 10000; NrunsH = 100000;  # MC runs for the threshold h and for all simulations 
change.interval = c(seq(0.05,0.5,0.05))*n;
a = floor(0.3*n); b = a + change.interval;		# Transient change period
Nb = length(b);


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
  detection = 1*( max(W) >= h ); # If h is exceeded, ch.pt. is detected
  return(detection);
}   # end of function MLE.single


### MC study with N(0,1) -> Laplace(0,1/sqrt(2))
Ns=1;

Prob.Detection = matrix(rep(0,Ns*Nb),Ns,Nb); 

h = threshold(n,NrunsH);  

for (k in 1:Nb){
  for (scena in 1:Ns){ 
     detection.run = rep(0,Nruns);
     for (run in 1:Nruns){
 
       Bern = rbinom(b[k]-a,1,0.5)*2-1;      # Sign
       Y = rexp(b[k]-a,rate=sqrt(2));

       X = c( rnorm(a,0,1), Bern*Y, rnorm(n-b,0,1) );

       MLE.result = MLE.single(X,h); 
       detection.run[run] = MLE.result; 
                    }  # End of MC runs for the given parameters
  Prob.Detection[scena,k] = mean(detection.run);     # Probability of detection | change
  } # End of the k-loop
} # End of scena-loop, over scenarios

Output = data.frame(Prob.Detection);

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\power-analysis-laplace.csv")

print(Output);
### Simulations. ROC as a function of threshold h

n = 100; alpha = 0.05; beta = 0.05;
Nruns = 10000; 
mu0 = 0; sigma0 = 1; sigma1 = 1;   

### Define a function for a single transient change-point estimation for known Normal distributions

MLE.single = function(X,m1,h){  # Estimation of a single transient change-point
  n = length(X);
  z = log( dnorm(X,mu1,sigma1)/dnorm(X,mu0,sigma0) );	# log-likelihood ratios
  S = cumsum(z);  			   							    # Partial sums of z
  W = S - cummin(S);			
  b.hat = which.max(W); a.hat = which.min(S[1:b.hat]);  # MLE
  detection = 1*( W[b.hat] >= h ); # If h is exceeded, ch.pt. is detected
  return(c(a.hat,b.hat,detection));
}   # end of function MLE.single


### MC study, N(0,1) -> N(mu1,1)

mu1 = c(0.2, 0.2, 0.5, 0.5, 1, 1); 
change.interval = c(0.1,0.2,0.1,0.2,0.1,0.2)*n;
a = 0.5*n; b = a+change.interval;
Nscena = length(mu1);
H = seq(0,10,0.4); NH = length(H);

Prob.False.Alarm = matrix(rep(0,NH*Nscena),NH,Nscena);
Prob.Detection = matrix(rep(0,NH*Nscena),NH,Nscena);

for (scena in 1:Nscena){
  for (k in 1:NH){ h = H[k];
     OUTPUT = data.frame(Nscena,scena,NH,k);
     write.table(OUTPUT,file="C:\\Users\\baron\\Documents\\Research\\Malov\\output.txt")
 
     false.alarm.run = rep(0,Nruns); detection.run = rep(0,Nruns);
     for (run in 1:Nruns){
       # The case of no change
       X = rnorm(n,mu0,sigma0); 
       MLE.result = MLE.single(X,mu1[scena],h); # Result = c(a.hat,b.hat,detection)  
       false.alarm.run[run] = MLE.result[3];
     
       # The case of a transient change
       X = c( rnorm(a,mu0,sigma0), rnorm(change.interval[scena],mu1[scena],sigma1), rnorm(n-b[scena],mu0,sigma0) );
       MLE.result = MLE.single(X,mu1[scena],h);
       detection.run[run] = MLE.result[3];
                    }  # End of MC runs for the given parameters
     Prob.False.Alarm[k,scena] = mean(false.alarm.run); # False alarm rate | no change
     Prob.Detection[k,scena] = mean(detection.run);     # Probability of detection | change
  } # End of the k-loop, over thresholds
} # End of the scena-loop, over (mu1, b-a)-scenarios

Output = data.frame(H, Prob.False.Alarm, Prob.Detection);
print(Output)

write.csv(Output, file = "C:\\Users\\baron\\Documents\\Research\\Malov\\ROC-results.csv")

color = c("red","blue","green","purple","orange","black");

plot( Prob.False.Alarm[,1], Prob.Detection[,1], lwd=2, col=color[1],
xlim=c(0,1), ylim=c(0,1), type="l",
xlab="Probability of a false alarm", ylab="Probability of detection",
main = "ROC curves for detecting a change in the mean" );
for (k in 2:Nscena){
lines( Prob.False.Alarm[,k], Prob.Detection[,k], lwd=2, col=color[k])
                    }
  

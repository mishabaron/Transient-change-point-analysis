Transient-change-point-analysis
R codes of Michael Baron and Sergey Malov for the detection and estimation of single and multiple transient changes

List of codes and their description:

OneTransientChangeExactThreshold.R
 Simulations. Generate data with one transient change.
 Estimate threshold h WITHOUT DOOB.
 Study FAR, FRR, the distribution of estimates.
 Normal(0,1) changes to Normal(mu1,1), known parameters

OneTransientChangeMean.R
 Simulations. Generate data with one transient change.
 Estimate threshold h. Detect the change-point. 
 Study FAR, FRR, the distribution of estimates.

OneTransientChangeNormalLaplace.R
 Simulations. Change from N(0,1) to Laplace(0,1/sqrt(2))

OneTransientChangeVariance.R
 Simulations. Generate data with one transient change in VARIANCE.
 Estimate threshold h WITHOUT DOOB.
 Normal(0,1) changes to Normal(0,sigma), known parameters

ROC.R
 Simulations. ROC as a function of threshold h

multiple-changes.R
 Simulations. Generate data with THREE transient changes. 
 Study FAR, FRR, the distribution of estimates.
 Thresholds are already calculated, we just download our results.

power-analysis-laplace.R
 Simulations. Generate data with one transient change.
 Estimate threshold h WITHOUT DOOB.
 Study the POWER for different sigma1 and (b-a)
 Normal(0,1) changes to Laplace(0,1/sqrt(2))

power-analysis-mean.R
 Simulations. Generate data with one transient change.
 Estimate threshold h WITHOUT DOOB.
 Study the POWER for different mu1 and (b-a)
 Normal(0,1) changes to Normal(mu1,1), known parameters

power-analysis-variance.R
 Simulations. Generate data with one transient change.
 Estimate threshold h WITHOUT DOOB.
 Study the POWER for different sigma1 and (b-a)
 Normal(0,1) changes to Normal(0,sigma1), known parameters

threshold-low-bound.R
 Simulations. Generate data with one transient change.
 Estimate threshold h WITHOUT DOOB.
 Study FAR, FRR, the distribution of estimates.
 Normal(0,1) changes to Normal(mu1,1), known parameters

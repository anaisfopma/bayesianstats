Hello!

Besides this readme file, you will find the following in the repository:
Report_Bayesian: This is the report.
Rcode_Bayesian: This is the R code to replicate the results presented in the report.

Below are instructions on which lines of code to execute to obtain the results presented in the report. 
Note 1: Please note that the code needs to be executed in this order to obtain the right results. 
Note 2: The packages "NHANES", "dplyr", "ggplot2", "gridEXtra", "rjags", and "bain" need to be installed.
Note 3: For some code (especially for the trace plots) it might be better to execute the code line by line instead of all at once. At least for me, executing it all at once seemed to be too much to handle for R.

TO OBTAIN:        trace plots for model 1
EXECUTE LINES:    1:193

TO OBTAIN:        parameter estimates, 95% CCIs, and MC errors for model 1
EXECUTE LINES:    195:216

TO OBTAIN:        trace plots for model 2
EXECUTE LINES:    218:305

TO OBTAIN:        parameter estimates, 95% CCIs, and MC errors for model 2
EXECUTE LINES:    307:328

TO OBTAIN:        posterior predictive p-value
EXECUTE LINES:    102:154 THEN 199:201 THEN 330:362

TO OBTAIN:        DIC for model 1
EXECUTE LINES:    364:397

TO OBTAIN:        DIC for model 2
EXECUTE LINES:    399:426

TO OBTAIN:        Bayes Factor
EXECUTE LINES:    428:443

TO OBTAIN:        frequentist parameter estimates and 95% CIs for model 1
EXECUTE LINES:    445:452

TO OBTAIN:        frequentist parameter estimates and 95% CIs for model 2
EXECUTE LINES:    454:457
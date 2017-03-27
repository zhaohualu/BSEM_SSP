These R scripts are provided under the GNU General Public License Version 3 (GPL V3)
A copy of the license can be obtained from 
https://www.gnu.org/licenses/gpl-3.0.en.html

Please cite the following papers if you used to R scripts in your research or publications.

Lu, Z.-H., Chow, S.-M., and Loken, E. (2016). Bayesian factor analysis as a
variable-selection problem: Alternative priors and consequences. Multivariate Behavioral
Research, 51(4):519–539.

Lu, Z.-H., Chow, S.-M., and Loken, E. (2016+). A Comparison of Bayesian and Frequentist Model Selection Methods for Factor Analysis Models. Underrevision.




main1.R is an example code. Procedure and explanation of the code are described below.

1. Simulation.R is a script for generating data based on a factor analysis model. 
Change the script for different setting. If empirical data are used, comment out the source("Simulation.R") and change the command for data input, i.e., "DataMat<-read.table("SimulatedData.txt",header=F)".

2. Change the dimension setting in main1.R 

3. Prior.R is a script for setting the hyperparameters in the prior distributions. Change the hyperparameters as needed. 

4. ind.R is a script for setting if the elements in the loading matrix and intercept to be fixed, estimated like common Bayesian CFA, or with SSP. Change the setting as needed. 

5. init1.R sets the initial values of the parameters. Change the settings and the values as needed.

6. The "BFA" function conduct the Bayesian factor analysis, the MCMC samples are returned as a list.

7. Calculate the posterior mean, SE and posterior model inclusion probabilities (if SSP is used) with the "MCMCSummary" function

8. If the model is a confirmatory model (no 2 in ind.R), use the "BayesianMCC" function to calculate the Bayesian model comparison criteria, including the marginal model probability for Bayes factor, BIC, DIC, and LOO-PSIS.

9. If the model uses SSP for some loading, 
9.1 In CandidateModelsIdentification.R, Describe the loading structure of the candidate models
9.2 use MMP_SSP function to calculate the marginal model probability for Bayes factor

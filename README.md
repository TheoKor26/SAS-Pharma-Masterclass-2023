# SAS Pharma Masterclass 2023
Project assignment on a RWE (Real-World Evidence) scenario about Asthma-exacerbation rates and biologics recommendation in random patients.
Following analysis steps are performed:
Simulation of routine data and generation of random patients;
Covariates preparation prior propensity score matching (PSM);
Propensity Score Modeling + Matching (model approach: logistic regression, matching criteria: 1:1 Nearest Neighbor);
Carry forward of matched observations;
Computing exacerbation rate after matching via Poisson regression for count data on matched dataset (model approach: generalized linear model, link function: Log);
Estimating the time to first exacerbation using a stratified Cox model (stratification by matched pairs) and tests the proportional hazard assumption;
Producing the Kaplan-Meier survival with treatment-effect stratification (biologics vs. standard);
Predicting patient characteristics for a guideline-compliant biologic recommendation through logistic regression (optimal decision cut-off via Youden's index);
Decision tree to translate the logistic model into clinically readable if-then rules and export scoring code to be embedded into practice management systems.

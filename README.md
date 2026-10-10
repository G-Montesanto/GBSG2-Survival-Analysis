Survival Analysis in R: GBSG2 Dataset

Overview

This project presents a survival analysis of the German Breast Cancer Study Group 2 (GBSG2) dataset using R.

The main objective is to investigate the relationship between clinical and biological characteristics and patient survival, with particular attention to the number of involved lymph nodes.

The analysis combines exploratory data analysis, Kaplan–Meier estimation, Cox proportional hazards regression, functional form assessment, model diagnostics, and sensitivity analysis.

Project Objectives

The main objectives of this project are:

* Explore the dataset and describe the main variables.
* Investigate the relationship between prognostic factors and survival.
* Estimate survival probabilities using Kaplan–Meier curves.
* Fit Cox proportional hazards models and assess their assumptions.
* Compare alternative functional forms for continuous covariates.
* Evaluate the potential benefits of logarithmic transformations and spline functions.
* Investigate influential observations and assess the stability of the results.
* Compare candidate models using statistical criteria and concordance indices.

Dataset

The analysis uses the GBSG2 dataset available through the TH.data R package.

The dataset contains information on 686 patients with breast cancer, including survival time, censoring status, treatment, tumor characteristics, hormonal receptor measurements, age, and the number of involved lymph nodes.

The main variables considered include:

* time: observed follow-up time.
* cens: event indicator used to identify whether the event was observed.
* horTh: hormonal therapy indicator.
* age: patient’s age.
* menostat: menopausal status.
* tumsize: tumor size.
* tgrade: tumor grade.
* pnodes: number of involved lymph nodes.
* progrec: progesterone receptor level.
* estrec: estrogen receptor level.

The observed follow-up time is converted from days to months for interpretability.

The dataset contains 686 observations and no missing values in the variables examined.

Methodology

1. Data Preparation and Exploratory Data Analysis

The first stage involves loading the dataset, checking its structure, examining missing values, and preparing the survival variables.

Descriptive statistics and graphical summaries are used to investigate the distributions of the main covariates.

Particular attention is given to the distribution of the number of involved lymph nodes, which is positively skewed.

The relationship between progesterone and estrogen receptor levels is also examined using Spearman’s rank correlation.

2. Kaplan–Meier Survival Analysis

Kaplan–Meier estimators are used to estimate the survival function.

Survival curves are examined for the overall sample and across selected clinical groups, including hormonal therapy, tumor grade, and categories of involved lymph nodes.

Log-rank tests are used to investigate differences between survival curves.

These analyses provide an initial, unadjusted assessment of the relationship between clinical characteristics and survival.

3. Cox Proportional Hazards Models

Cox proportional hazards regression is used to investigate associations between patient characteristics and the hazard of the event.

The analysis includes univariable models and multivariable models that account for multiple covariates simultaneously.

Hazard ratios and their confidence intervals are used to quantify the estimated associations.

The proportional hazards assumption is assessed using Schoenfeld residuals and the cox.zph() function.

Where appropriate, stratified Cox models are considered to accommodate covariates for which a common baseline hazard may be unsuitable.

4. Functional Form Assessment

The functional relationship between continuous covariates and the log hazard is an important aspect of Cox model specification.

Particular attention is given to pnodes, the number of involved lymph nodes.

Alternative specifications are investigated, including:

* A linear functional form.
* A logarithmic transformation using log1p().
* Natural spline functions.

Martingale residuals are used as an exploratory tool to investigate potential departures from linearity.

The candidate models are compared using the Akaike Information Criterion (AIC), while likelihood ratio tests are used when the models being compared are nested.

This approach helps evaluate whether a more flexible functional form provides a meaningful improvement over a simpler specification.

5. Hazard Ratio Curves

Hazard ratio curves are constructed to examine how the estimated hazard changes with the number of involved lymph nodes.

The hazard ratios are expressed relative to a reference value of three involved lymph nodes.

Pointwise confidence intervals are calculated to represent uncertainty in the estimated hazard ratios at individual covariate values.

These curves help illustrate potential nonlinear associations that may not be adequately represented by a single constant hazard ratio per additional lymph node.

6. Model Diagnostics and Influential Observations

Model diagnostics are used to investigate the stability and adequacy of the fitted Cox models.

The analysis includes:

* Assessment of the proportional hazards assumption.
* Examination of influential observations using DFBETA diagnostics.
* A sensitivity analysis in which the most influential observations are temporarily excluded.
* Comparison of model estimates before and after the sensitivity analysis.

The purpose is to assess whether a small number of observations substantially affect the estimated associations.

The primary analysis retains the original observations. The sensitivity analysis is exploratory and does not, by itself, justify excluding observations from the main model.

Main Findings

Survival Patterns

Kaplan–Meier analyses indicate differences in survival across selected clinical groups.

In particular, survival differs across categories defined by the number of involved lymph nodes, suggesting that this variable is an important prognostic factor in the dataset.

These unadjusted comparisons describe survival patterns but do not account for potential confounding by other covariates.

Number of Involved Lymph Nodes

The number of involved lymph nodes shows a strong association with the hazard of the event.

In the univariable Cox model with a linear term for pnodes, the estimated hazard ratio is approximately 1.06 per additional involved lymph node.

This corresponds to an estimated increase of approximately 6% in the hazard per additional lymph node, under the model’s linear functional form and proportional hazards assumptions.

The estimate is an association, not a causal effect, and its interpretation depends on the adequacy of the model specification.

Functional Form Comparison

The analysis compares linear, logarithmic, and spline-based specifications for the number of involved lymph nodes.

The likelihood ratio test comparing the linear model with the specified spline model provides evidence that the spline offers a better fit than the linear specification.

However, the logarithmic model has a slightly lower AIC than the particular spline model initially considered.

This illustrates the importance of distinguishing between comparisons of nested models using likelihood ratio tests and comparisons of candidate models using information criteria.

The findings do not establish that one functional form is universally preferable. Additional sensitivity analyses using alternative spline specifications would help assess the robustness of this conclusion.

Proportional Hazards Assumption

The proportional hazards assumption is assessed using Schoenfeld residuals.

For the stratified spline model examined, the global test does not provide statistically significant evidence against the proportional hazards assumption.

A nonsignificant test does not prove that the assumption is satisfied, but it indicates that the diagnostic test did not detect a clear violation.

Sensitivity Analysis

The sensitivity analysis examines the effect of temporarily excluding the five most influential observations identified by the diagnostic procedure.

The estimated concordance index changes only modestly, and the hazard ratio curves remain broadly similar across much of the central range of the number of involved lymph nodes.

Greater differences appear in the upper tail, where relatively few observations are available.

This suggests that the central part of the estimated relationship is reasonably stable under this particular sensitivity analysis, while estimates at extreme covariate values should be interpreted more cautiously.

Model Comparison and Interpretation

The logarithmic and spline-based models have similar concordance indices.

The difference between their in-sample C-indices is approximately 0.002, with the spline model showing a slightly higher value.

Such a small difference should not be interpreted as evidence that the spline model has superior predictive performance.

Because the concordance indices are calculated on the same observations used to estimate the models, they may be optimistic. Bootstrap validation or cross-validation would be useful for a more reliable assessment of predictive discrimination.

Overall, the results support an association between the number of involved lymph nodes and the hazard of the event. The analysis also suggests that the functional form used to represent this relationship deserves careful consideration.

Model selection should account for fit, complexity, diagnostic results, interpretability, and the stability of the estimates rather than relying on a single metric.

Methodological Limitations

Several methodological limitations should be considered when interpreting the results.

Spline Specification

The initial spline model uses manually selected knots at 3, 7, and 15 involved lymph nodes.

The sensitivity of the results to alternative knot locations and spline degrees of freedom has not yet been fully assessed.

Additional analyses could compare natural splines with different degrees of freedom and alternative knot placements, including placements based on the observed covariate distribution.

In-Sample Concordance Index

The reported concordance indices are calculated using the same data used to fit the models.

Consequently, they may overestimate performance on new observations.

Small differences between the C-indices should not be interpreted as conclusive evidence of superior predictive performance.

Bootstrap validation or cross-validation would provide a more informative assessment of model discrimination and optimism.

Functional Form of Other Covariates

The functional form assessment focuses primarily on the number of involved lymph nodes.

Further checks within the multivariable model are needed to evaluate whether alternative functional forms for age and tumor size would improve the model specification.

Confidence Intervals for Hazard Ratio Curves

The confidence intervals displayed for hazard ratio curves are pointwise rather than simultaneous.

They quantify uncertainty at individual values of the covariate and should not be interpreted as providing a specified coverage probability for the entire curve.

Influential Observations

The sensitivity analysis provides information about the stability of the results when influential observations are temporarily excluded.

However, influential observations are not necessarily erroneous. Their exclusion should not be justified solely by their influence on model estimates.

The primary analysis therefore retains the original observations.

Observational Nature of the Data

The dataset is observational, and the estimated associations do not establish causal relationships.

Potential confounding, model misspecification, and the limitations of the available variables should be considered when interpreting the findings.

Tools and Technologies

The project uses R and the following packages:

* TH.data for the GBSG2 dataset.
* survival for Kaplan–Meier estimation and Cox regression.
* survminer for survival analysis visualization.
* dplyr for data manipulation.
* ggplot2 for data visualization.
* splines for flexible modeling of continuous covariates.
* car for additional hypothesis testing.

How to Run

To reproduce the analysis:

1. Install R.
2. Install the required packages if they are not already available.
3. Download or clone this repository.
4. Open analysis.R in RStudio or another R-compatible environment.
5. Run the script from beginning to end.

The dataset is loaded from the TH.data package, so a separate local copy of the dataset is not required.

The script contains comments describing the main analytical steps and model comparisons.

Project Structure

The repository is organized as follows:

GBSG2-Survival-Analysis/
├── README.md
└── analysis.R

Additional folders, such as figures/, may be added to store exported plots and other project outputs.

Conclusion

This project demonstrates the application of survival analysis techniques to a public clinical dataset using R.

By combining Kaplan–Meier estimation, Cox proportional hazards regression, alternative functional forms, model diagnostics, and sensitivity analysis, the study investigates the relationship between clinical characteristics and survival.

The analysis highlights the prognostic relevance of the number of involved lymph nodes and illustrates the importance of evaluating functional form assumptions and model stability.

Further work should focus on assessing sensitivity to spline specifications, evaluating the functional forms of additional covariates in the multivariable model, and validating predictive performance using bootstrap or cross-validation.

The project provides a foundation for further development in statistical modeling, survival analysis, and reproducible data analysis.
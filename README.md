Survival Analysis in R: GBSG2 Dataset

A survival analysis project using R to investigate prognostic factors in breast cancer patients, with a particular focus on the functional relationship between the number of positive lymph nodes and the hazard of the event.

The project combines exploratory data analysis, Kaplan–Meier estimation, Cox proportional hazards regression, nonlinear modeling with natural splines, model diagnostics, and sensitivity analysis.

Note: The R script is commented in Italian.

Overview

Survival analysis is used to study the time until an event occurs while accounting for observations whose event time is not fully observed.

This project applies survival analysis techniques to the German Breast Cancer Study Group 2 (GBSG2) dataset. The main objective is to investigate the association between patient characteristics and survival outcomes and to assess how the specification of continuous predictors affects Cox model estimates.

Particular attention is given to pnodes, the number of positive lymph nodes, whose relationship with the hazard may not be adequately represented by a simple linear term.

The project emphasizes statistical interpretation, model diagnostics, and the distinction between goodness of fit and predictive performance.

Project Objectives

The main objectives are to:

* Explore the dataset and describe its main clinical and biological characteristics.
* Estimate survival functions using Kaplan–Meier methods.
* Compare survival patterns across selected patient groups.
* Fit univariable and multivariable Cox proportional hazards models.
* Assess the proportional hazards assumption.
* Investigate alternative functional forms for the number of positive lymph nodes.
* Compare linear, logarithmic, and natural spline specifications.
* Visualize hazard ratios relative to a clinically interpretable reference value.
* Examine influential observations using DFBETA diagnostics.
* Assess the sensitivity of model estimates to the temporary exclusion of influential observations.
* Compare candidate models using the Akaike Information Criterion (AIC), likelihood ratio tests, and concordance indices.

Dataset

The analysis uses the GBSG2 dataset provided by the R package TH.data.

The dataset contains information on 686 breast cancer patients, including follow-up time, event status, treatment, tumor characteristics, hormonal receptor measurements, age, and the number of positive lymph nodes.

The main variables include:

* time: observed follow-up time.
* cens: event indicator.
* horTh: hormonal therapy status.
* age: patient age.
* menostat: menopausal status.
* tumsize: tumor size.
* tgrade: tumor grade.
* pnodes: number of positive lymph nodes.
* progrec: progesterone receptor level.
* estrec: estrogen receptor level.

The follow-up time is converted from days to months to make the results easier to interpret.

The dataset contains 686 observations, with 299 observed events. No missing values were identified in the variables examined for this analysis.

Methodology

1. Data Preparation and Exploratory Data Analysis

The analysis begins by loading the dataset, examining its structure, checking for missing values, and preparing the variables required for survival modeling.

Descriptive statistics and graphical summaries are used to explore the distributions of the main covariates.

Particular attention is given to the distribution of pnodes, which is positively skewed. This feature motivates the investigation of transformations and flexible functional forms.

The association between progesterone and estrogen receptor levels is also explored using Spearman’s rank correlation.

Exploratory analyses provide context for the subsequent models but do not establish causal relationships.

2. Kaplan–Meier Survival Analysis

Kaplan–Meier estimators are used to estimate the survival function for the overall sample and for selected clinical groups.

Survival curves are examined according to:

* Hormonal therapy status.
* Tumor grade.
* Categories of positive lymph nodes.

Log-rank tests are used to investigate differences between selected survival curves.

These comparisons provide an initial, unadjusted assessment of survival patterns. They do not account for potential confounding by other patient characteristics.

3. Cox Proportional Hazards Regression

Cox proportional hazards models are fitted to investigate the association between patient characteristics and the hazard of the event.

Both univariable and multivariable specifications are considered. The multivariable models account for multiple covariates simultaneously and use transformations for selected variables where specified in the analysis.

Hazard ratios and confidence intervals are used to quantify the estimated associations.

The proportional hazards assumption is assessed using Schoenfeld residuals and the cox.zph() function.

Stratified Cox models are also considered to allow different baseline hazard functions across selected categorical variables.

A nonsignificant proportional hazards test is not proof that the assumption holds exactly; it indicates only that the diagnostic test did not detect statistically significant evidence against it.

4. Functional Form Assessment

A central part of the project is the investigation of how pnodes should enter the Cox model.

Three functional specifications are compared:

* Linear: the log hazard changes linearly with the number of positive lymph nodes.
* Logarithmic: the predictor is transformed using log1p(pnodes).
* Natural spline: a flexible function allows the association to depart from a simple linear relationship.

Martingale residuals are examined as an exploratory diagnostic for possible departures from linearity.

Candidate models are compared using the AIC. A likelihood ratio test is used to compare the nested linear and spline specifications.

These criteria answer different questions. The likelihood ratio test evaluates whether the additional flexibility of the specified spline improves model fit relative to the linear specification, whereas the AIC balances fit against model complexity.

Consequently, a significant likelihood ratio test does not necessarily imply that the spline has the lowest AIC among all candidate models.

5. Hazard Ratio Curves

Hazard ratio curves are used to visualize the estimated association between pnodes and the hazard.

The curves are expressed relative to a reference value of three positive lymph nodes. This makes it possible to interpret the estimated hazard at other values relative to the same reference.

The logarithmic and spline specifications are compared to illustrate how the choice of functional form affects the estimated relationship.

Pointwise 95% confidence intervals are used to represent uncertainty at individual values of pnodes. They are not simultaneous confidence bands for the entire curve.

Estimates near the extremes of the predictor distribution should be interpreted cautiously because relatively few observations may be available in those regions.

6. Model Diagnostics and Influential Observations

Model diagnostics are used to investigate the adequacy and stability of the fitted models.

The analysis includes:

* Assessment of the proportional hazards assumption using Schoenfeld residuals.
* Examination of influential observations using DFBETA residuals.
* Identification of observations with the largest absolute DFBETA values.
* A sensitivity analysis that temporarily excludes the five most influential observations.
* Comparison of model summaries and concordance indices before and after this exclusion.

The sensitivity analysis is exploratory. An observation can be influential without being erroneous, and influence alone is not sufficient justification for removing a patient from the primary analysis.

The original observations are retained in the primary model.

Main Findings

Survival Patterns

The Kaplan–Meier analyses describe differences in estimated survival across selected clinical groups, including categories defined by the number of positive lymph nodes.

These comparisons are unadjusted and should be interpreted as descriptive evidence of differences in survival patterns rather than as estimates of independent or causal effects.

Association Between Positive Lymph Nodes and Hazard

In the univariable Cox model with a linear term for pnodes, the estimated hazard ratio is approximately 1.06 per additional positive lymph node.

Under the linear specification and proportional hazards assumption, this corresponds to an estimated increase of approximately 6% in the hazard per additional positive lymph node.

This interpretation depends on the assumption that the log hazard changes linearly with pnodes. If the true relationship is nonlinear, a single hazard ratio per additional lymph node may not adequately describe the association.

The estimate represents an association, not a causal effect.

Comparison of Functional Forms

The comparison of candidate specifications produces the following approximate results:

* The logarithmic model has an AIC of 2597.77.
* The specified spline model has an AIC of 2600.18.
* The likelihood ratio test comparing the linear and spline models gives a p-value of approximately 0.0001.
* The in-sample concordance index is approximately 0.682 for the logarithmic model and 0.684 for the spline model.

The likelihood ratio test provides evidence that the specified spline improves fit relative to the linear specification.

However, the logarithmic model has a slightly lower AIC than the spline model considered here. The spline model has a marginally higher in-sample concordance index, but the difference is approximately 0.002.

These results are not contradictory because the criteria measure different aspects of model performance.

The findings do not establish that one functional form is universally preferable. The choice should also consider model complexity, diagnostic results, interpretability, and stability.

Proportional Hazards Assumption

The proportional hazards assumption is assessed using Schoenfeld residuals and cox.zph().

For the stratified spline model examined in the analysis, the global test does not provide statistically significant evidence against the proportional hazards assumption.

This result should not be interpreted as proof that proportional hazards holds exactly. It indicates that the test did not detect a statistically significant violation.

Sensitivity to Influential Observations

The DFBETA diagnostics identify observations that have relatively large influence on the estimated model coefficients.

The five most influential observations are temporarily excluded in a sensitivity analysis, and the model is refitted using the remaining observations.

The analysis script reports an approximate change in the concordance index from 0.684 in the original model to 0.687 in the sensitivity model.

This relatively small change suggests that the concordance index is not substantially altered by this particular exclusion. Nevertheless, a similar concordance index does not guarantee that all regression coefficients or hazard ratio curves remain unchanged.

Differences may be more pronounced in the upper tail of pnodes, where observations are less numerous and flexible functional forms can be less precisely estimated.

The sensitivity analysis does not justify excluding these patients from the primary analysis. It illustrates why influential observations and sparse regions of the predictor distribution deserve careful examination.

Interpretation of Predictive Performance

The logarithmic and spline models have similar in-sample concordance indices.

The difference of approximately 0.002 is too small to support a meaningful claim that the spline model has superior predictive performance.

Moreover, these concordance indices are calculated on the same observations used to estimate the models. They may therefore be optimistic estimates of performance on new patients.

Bootstrap validation or cross-validation would be useful to assess predictive discrimination more reliably.

Methodological Limitations

Several limitations should be considered when interpreting the findings.

Spline Specification

The spline model uses manually specified knots at 3, 7, and 15 positive lymph nodes, with boundary knots at 1 and 51.

The sensitivity of the findings to alternative knot locations and degrees of freedom has not been fully assessed. Further work could compare alternative spline specifications and knot placements based on the observed distribution of pnodes.

In-Sample Concordance

The reported concordance indices are calculated on the data used to fit the models.

They should not be interpreted as validated estimates of predictive performance on unseen observations. Bootstrap validation or cross-validation would provide a more reliable assessment.

Functional Forms of Other Covariates

The functional form investigation focuses primarily on pnodes.

Although selected transformations are used for other continuous covariates, additional systematic checks within the multivariable model could assess whether alternative specifications for age and tumor size improve model adequacy.

Confidence Intervals for Hazard Ratio Curves

The confidence intervals displayed for hazard ratio curves are pointwise rather than simultaneous.

They quantify uncertainty at individual covariate values and do not provide a specified coverage probability for the entire estimated curve.

Influential Observations

The sensitivity analysis examines the effect of temporarily excluding the five observations with the largest values of the chosen influence measure.

This is one exploratory diagnostic, not a universal procedure for identifying invalid observations. Influential cases may be legitimate patients whose characteristics are uncommon.

The original data are therefore retained in the primary analysis.

Observational Data

The dataset is observational. The estimated associations do not establish causal relationships.

Potential confounding, model misspecification, and limitations of the available variables should be considered when interpreting the findings.

Figures

The repository’s figures/ directory contains the graphical outputs associated with the analysis.

These figures support the exploration of survival patterns, comparisons between candidate models, hazard ratio estimation, and model diagnostics.

See the figures/ directory for the available plots.

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
2. Clone or download this repository.
3. Open analysis.R in RStudio or another R-compatible environment.
4. Install any required packages that are not already available.
5. Run the script from beginning to end.

The dataset is loaded from the TH.data package, so a separate local copy is not required.

The script contains comments in Italian describing the main analytical steps, model comparisons, and interpretation of results.

Project Structure

The repository is organized as follows:

GBSG2-Survival-Analysis/
├── README.md
├── analysis.R
└── figures/
    └── Graphical outputs from the analysis

Conclusion

This project demonstrates the application of survival analysis techniques to a clinical dataset using R.

By combining Kaplan–Meier estimation, Cox proportional hazards regression, alternative functional forms, model diagnostics, and sensitivity analysis, the project investigates the association between patient characteristics and the hazard of the event.

The number of positive lymph nodes emerges as an important prognostic variable in the analysis. The comparison of linear, logarithmic, and spline specifications also illustrates how model specification can influence the interpretation of this association.

The results highlight the importance of distinguishing model fit from predictive performance and of evaluating assumptions, influential observations, and uncertainty before drawing conclusions.

Further work should assess sensitivity to alternative spline specifications, investigate the functional forms of additional covariates, and validate predictive performance using bootstrap methods or cross-validation.

Overall, the project provides a foundation for further work in survival analysis, statistical modeling, and reproducible data analysis.

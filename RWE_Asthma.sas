/* ==============================================================================
   10/11/2023-Case Study RWE-Asthma: Building Cohorts + Propensity Score Matching
   ============================================================================== */

/* 1) Generate sample dataset (simulated routine data)
   Variables:
- PatientID
- age, sex
- asthma_sev_prev (exacerbations in the pre-period, 0–5)
- eosinophils (blood eosinophils, 10^9/L)
- comorbid_rhino (rhinosinusitis, 0/1)
- comorbid_adip (obesity, 0/1)
- ics_dose_prev (mean ICS dose in the pre-period, numeric)
- ode_cure_prev (OCS courses in the pre-period, count)
- biologikum_prev (biologics in the pre-period, 0/1)
- treatment (0=ICS/LABA, 1=biologic)
- index_date (index date as SAS date)
- exacerb_12m (exacerbations in the 12-month follow-up)
- hosp_12m    (1=asthma-related hospitalization, 0=no)
*/

data asthma_raw;
  length PatientID $10;
  format index_date date9.;
  do i = 1 to 3000;
  do month = 3, 6, 9, 12;
    quarter = int((month - 1) / 3) + 1;  /* 1,2,3,4 */

/* Quarterly exacerbations (simulated) */
    /* Assumption: Exacerb_12m is over 12 months, we divide roughly by quarter */
    rate_per_quarter = exacerb_12m / 4;
/* Random variation around this value */
    exacerb_quarter = max(0, int(rate_per_quarter + rand("normal")*0.5));

/* Quarterly hospitalization (0/1, simplified) */
    hosp_quarter = (rand("uniform") < hosp_12m / 4);
    PatientID = cat("P", put(i, 3.));
    age       = 12 + int(60*rand("uniform"));   /* 12–71 Jahre */
    sex       = ifc(rand("uniform") < 0.5, "M", "F");

    asthma_sev_prev = min(5, int(3*rand("uniform") + rand("normal")*1));
    if asthma_sev_prev < 0 then asthma_sev_prev = 0;

    eosinophils     = 0.1 + 0.8*rand("uniform") + 0.3*rand("normal");
    if eosinophils < 0.05 then eosinophils = 0.05;

    comorbid_rhino  = (rand("uniform") < 0.4);
    comorbid_adip   = (rand("uniform") < 0.3);

    ics_dose_prev   = 200 + 600*rand("uniform");
    ode_cure_prev   = min(4, int(2*rand("uniform") + rand("normal")*0.7));
    if ode_cure_prev < 0 then ode_cure_prev = 0;

    biologikum_prev = (rand("uniform") < 0.15);

/* Treatment: more severe patients more likely to receive a biologic */
    logit_p = -1.5
            + 0.7*(asthma_sev_prev >= 2)
            + 0.6*(eosinophils > 0.3)
            + 0.4*ode_cure_prev
            + 0.5*biologikum_prev
            + 0.3*comorbid_rhino;
    p_treat = exp(logit_p) / (1 + exp(logit_p));
    treatment = (rand("uniform") < p_treat);

/* Index date: randomly between 2018-01-01 and 2023-12-31 */
    index_date = '01JAN2018'd + int(365*6*rand("uniform"));

/* Outcomes (simulated, with a slight effect of treatment) */
    base_exac = 0.5 + 0.4*asthma_sev_prev + 0.3*(eosinophils > 0.3)
                + 0.2*ode_cure_prev - 0.25*treatment;
    exacerb_12m = max(0, int(base_exac + rand("normal")*0.8));

    base_hosp = -1.5 + 0.6*asthma_sev_prev + 0.4*(eosinophils > 0.3)
                + 0.3*ode_cure_prev - 0.4*treatment;
    logit_h = base_hosp;
    phosp   = exp(logit_h) / (1 + exp(logit_h));
    hosp_12m = (rand("uniform") < phosp);

    output;
    end;
  end;
run;

data asthma_raw2;
  set asthma_raw;
  format study_end_date event_date censor_date date9.;
  study_end_date = intnx('month', index_date, 12, 'same');

  if exacerb_quarter = 1 then do;
    event = 1;
/* Approximated event time, e.g., end of month 3 */
    event_date = intnx('month', index_date, 3, 'same');
    time_event = event_date - index_date;
    censor_date = event_date;
  end;
  else do;
    event = 0;
    censor_date = study_end_date;
    time_event = censor_date - index_date;
  end;
  drop i;
run;


/* 2) Prepare data: Model covariates for PSM
   We use only variables from the pre-period (prior to index_date).
*/

data asthma_pre;
  set asthma_raw2;
  /* Binary covariates for PSM */
  sev2      = (asthma_sev_prev >= 2);
  eos_high  = (eosinophils > 0.3);
  rhino     = comorbid_rhino;
  adip      = comorbid_adip;

  /* Treatment as binary: 1 = Biologikum, 0 = ICS/LABA */
  treat_bin = treatment;
run;

/* 3) PROC PSMATCH: Propensity Score Modeling + Matching
- PSMODEL: logistic regression for P(Treatment=1 | covariates)
- MATCH: 1:1 Nearest Neighbor, caliper = 0.2 SD of the logit
- ASSESS: Balance check (standardized differences)
- OUTPUT: matched dataset with _PS_, _LPS_, _MatchID_
*/

ods graphics on;

proc psmatch data=asthma_pre region=cs;
  /* classification variables */
  class sex treat_bin sev2 eos_high rhino adip;

  /* Propensity Score Modell: logistic regression */
  psmodel treat_bin(Treated='1') = sex age sev2 eos_high rhino adip ics_dose_prev ode_cure_prev;

  /* Matching: 1:1 Nearest Neighbor, Kaliper=0.2 */
  match method=GREEDY(k=1) stat=lps caliper=0.2;

  /* Balance-Check */
  assess lps var=(sex age sev2 eos_high rhino adip ics_dose_prev ode_cure_prev)
         / plots=(boxplot barchart stddiff);

  /* Output: matched dataset */
  output out(obs=match)=asthma_matched
         lps=_LPS_
         matchid=_MatchID_;
run;

ods graphics off;

/* 4) describing matched dataset */

proc sort data=asthma_matched out=asthma_matched_cls;
  by _MatchID_;
run;

proc print data=asthma_matched_cls(obs=20);
  var PatientID sex age month sev2 eos_high rhino adip ics_dose_prev ode_cure_prev
      treat_bin _PS_ _LPS_ _MatchID_;
run;

/* 5) Check balance: standardized differences before/after matching
   We calculate these explicitly for several key variables.
*/

/* Unmatched: all observations */
proc means data=asthma_pre noprint;
           output out=means_pre
           var=age month ics_dose_prev ode_cure_prev asthma_sev_prev eosinophils;
run;

proc means data=asthma_pre noprint;
           output out=means_pre_t
           var=age month ics_dose_prev ode_cure_prev asthma_sev_prev eosinophils;
  where treat_bin = 1;
run;

proc means data=asthma_pre noprint;
           output out=means_pre_c
           var=age month ics_dose_prev ode_cure_prev asthma_sev_prev eosinophils;
  where treat_bin = 0;
run;

/* Matched:  matched observations only */
proc means data=asthma_matched_cls noprint;
           output out=means_match
           var=age month ics_dose_prev ode_cure_prev asthma_sev_prev eosinophils;
run;

proc means data=asthma_matched_cls noprint;
           output out=means_match_t
           var=age month ics_dose_prev ode_cure_prev asthma_sev_prev eosinophils;
  where treat_bin = 1;
run;

proc means data=asthma_matched_cls noprint;
           output out=means_match_c
           var=age month ics_dose_prev ode_cure_prev asthma_sev_prev eosinophils;
  where treat_bin = 0;
run;

/* 6) Example analysis: Exacerbation rate after matching
   Poisson regression on matched dataset, without further adjustment
   (since balance is good).
   The result is incidence rate ratios, which directly indicate by what percentage the biologic
   reduces the exacerbation rate.
*/

/* Mixed Model: Exacerbation rate (count data) */
/* ============================================================
   Generalized Linear Mixed Model (GLMM)
   Outcome:    Exacerbation count (count → Poisson/NegBin)
   Link:       Log
   ============================================================ */

proc glimmix  data=asthma_matched_cls method=rspl;
  class   PatientID treat_bin sex sev2;
  model exacerb_12m = treat_bin|month age sev2  comorbid_adip  comorbid_rhino/ dist=poisson link=log
  																		 SOLUTION ODDSRATIO;
  RANDOM INTERCEPT  / SUBJECT=PatientID type=vc;
    /* Inzidenzraten-Ratios (Kontraste) ausgeben */
  ESTIMATE "Biologikum vs. Standard"
           treat_bin 1 -1  / EXP CL;
  LSMEANS treat_bin / ILINK PDIFF ADJUST=BON;  /* Back-transformiert auf Skala */
  TITLE "GLMM Poisson: Exacerbation rate over time";
run;


/* 7) – Cox estimates the time to first exacerbation using a stratified
   Cox model (stratification by matched pairs) and tests the PH assumption.
   In addition, KM curves are included for visualization.
*/
/* Cox regression: Time to first exacerbation */
/* ============================================================
   Cox model on matched dataset
   Stratification by matched pairs (paired Cox)
   ============================================================ */

ods graphics on;
PROC PHREG DATA = asthma_matched_cls COVS(AGGREGATE);

  CLASS treat_bin(REF='1') sex sev2(REF='1') comorbid_adip  comorbid_rhino
        / PARAM=REF;

/* Stratification by Match-ID = conditional analysis */
  MODEL time_event * event(0) = treat_bin age sex sev2 comorbid_adip comorbid_rhino  eosinophils ics_dose_prev ode_cure_prev
        / TIES=EFRON RL;

  STRATA  _MatchID_;                  /* looking at matched-Pairs */
  hazardratio treat_bin / diff=ref;
  baseline out=cox_surv
           survival=surv
           cumhaz=cumhaz
           / rowid=treat_bin;
  TITLE "Stratified Cox: Time to exacerbation (PSM pairs) ";
RUN;
ods graphics off;
/* Then you get an output table with estimated survival and cumulative hazard,
   which you can plot later. */

/* Kaplan-Meier on matched dataset */
ods graphics on;
PROC LIFETEST DATA  = asthma_matched_cls
              PLOTS = survival(atrisk cl);
  TIME   time_event * event(0);
  STRATA treat_bin / TEST=LOGRANK;
  TITLE  "KM curves after PSM: Biologic vs. Standard ";
RUN;
ods graphics off;

/* 8) Logistic regression identifies which patient characteristics predict a guideline-compliant
   biologic recommendation, and calculates the optimal decision cut-off via Youden's index.
*/
/* PART B – Guideline Recommendations: Biologic Prescription Decision */

/* B1 – Logistic regression: Who should receive a biologic? */
/* ============================================================
   Objective: Identification of clinical thresholds
   Outcome: Biologic suitability (1=yes, 0=no)
   Basis: Clinical guidelines + biomarkers
   ============================================================ */

/* Derive recommendation variable from clinical criteria */
DATA asthma_raw3;
  SET asthma_raw2;

  /* Guideline-based criteria (e.g., GINA Step 4/5): */
  bio_suitab = 0;
  IF asthma_sev_prev >= 3 AND
     (eosinophils >= 0.30 OR base_exac >=2)
                          THEN bio_suitab = 1;

  /* Biomarker categories for subgroup analysis*/
  IF      eosinophils <  0.15 THEN eos_group = "Low";
  ELSE IF eosinophils <  0.40 THEN eos_group = "Average";
  ELSE                          eos_group = "High";

  IF      base_exac < 0.5  THEN exac_group = "Low";
  ELSE IF base_exac <  2  THEN exac_group = "Average";
  ELSE                         exac_group = "High";
RUN;

/* ── Logistic model: Predictors of the recommendation ────── */
PROC LOGISTIC DATA = asthma_raw3 DESCENDING PLOTS(ONLY)=(ROC EFFECT);

  CLASS sex   (REF='F')
        asthma_sev_prev (REF='1')
        eos_group (REF='Low')
        exac_group (REF='Low')
        / PARAM=REF;

  MODEL bio_suitab =
        age sex base_exac
        eosinophils
        asthma_sev_prev comorbid_adip  comorbid_rhino
        / FIRTH           /* Penalized MLE – stable under separation */
          LACKFIT                   /* Hosmer-Lemeshow */
          CLPARM=PL       /* Profil-Likelihood KI*/
          CTABLE PPROB=0.5          /* Classification Table            */
          RISKLIMITS                /* Odds Ratios                       */
          OUTROC=roc_data;          /* Export ROC-curve           */

  /* Optimal Cut-off via Youden-Index */
  OUTPUT OUT  = logit_scored
         PRED = p_suitab
         XBETA= logit_score;

  TITLE "/* Logistic regression: Predictors of biologic recommendation */";
RUN;

/* ── Optimal threshold ─────────────────────────────── */
DATA roc_optimal;
  SET roc_data;
  youden = _SENSIT_ + (1 - _1MSPEC_) - 1;
RUN;

PROC SORT DATA=roc_optimal; BY DESCENDING youden; RUN;

TITLE "Optimal Cut-off (Youden-Index)";
PROC PRINT DATA=roc_optimal (OBS=1);
  VAR _PROB_ _SENSIT_ _1MSPEC_ youden;
RUN;
title;

/* 9) – Decision tree translates the logistic model into clinically readable if-then rules and
   exports scoring code that can be embedded directly into practice management systems.
*/

/* B3 – Decision tree: Clinical prescription rule */
/* ============================================================
   PROC HPSPLIT: Decision tree for prescription logic
   Provides interpretable if-then rules for clinicians
   ============================================================ */

PROC HPSPLIT DATA = asthma_raw3 seed=12345;

  CLASS  sex;
  MODEL bio_suitab= sex age base_exac
        eosinophils asthma_sev_prev comorbid_adip  comorbid_rhino ;

  GROW VARIANCE;
  PRUNE COSTCOMPLEXITY;

  CODE FILE="/home/u50238372/SAS Pharma Master Class/scoring_rule.sas"; /* export rules */

  TITLE " Decision tree: Biologic prescription rule ";
RUN;

# CLIF ICU Respiratory Failure Trajectories, Air Pollution Exposure, and Outcomes

## CLIF VERSION
2.1

---

# Objective

This project evaluates whether **long-term ambient air pollution exposure (PM₂.₅ and NO₂)** is associated with **early physiologic trajectories of acute respiratory failure (ARF)** and clinical outcomes among critically ill patients admitted to the ICU.

Using high-resolution physiologic data from the **Critical Care Informatics Framework (CLIF)**, the study reconstructs **hourly respiratory and physiologic trajectories during the first 72 hours of ICU care**. These trajectories are then grouped into phenotypic clusters representing distinct patterns of respiratory failure evolution.

The analysis evaluates whether chronic environmental exposure is associated with:

- Distinct early respiratory failure trajectory phenotypes
- Increased severity of respiratory failure
- Greater risk of invasive mechanical ventilation
- Higher ICU mortality
- Longer ICU length of stay

Because CLIF data are distributed across multiple institutions, this project is designed as a **federated analysis pipeline** in which each site executes the analysis locally and exports **PHI-safe cluster summaries and model results**.

---

# Study Design Overview

The analytic workflow consists of four major stages:

### 1. ICU Cohort Construction

Eligible ICU hospitalizations are identified using CLIF administrative tables.

### 2. Hourly Physiologic Trajectory Reconstruction

During the first **0–72 hours after ICU admission**, the following physiologic domains are reconstructed hourly:

- Respiratory support state
- Acute respiratory failure subtype
- Vasopressor exposure

### 3. Trajectory Phenotyping

Sequence clustering is performed using **TraMineR optimal matching distances** to identify distinct respiratory trajectory phenotypes.

### 4. Environmental Exposure Modeling

Multivariable models evaluate whether long-term exposure to:

- PM₂.₅
- NO₂

is associated with trajectory phenotype membership and clinical outcomes.

---

# Required CLIF Tables and Fields

Please refer to the official resources for detailed table specifications:

- [CLIF Data Dictionary](https://clif-icu.com/data-dictionary)  
- [CLIF Tools](https://clif-icu.com/tools)  
- [CLIF ETL Guide](https://clif-icu.com/etl-guide)

The following CLIF tables are required.

---

# 1. `patient`

Required fields:

- `patient_id`
- `birth_date`
- `sex_category`
- `race_category`
- `ethnicity_category`

Used to determine demographic characteristics and age.

---

# 2. `hospitalization`

Required fields:

- `patient_id`
- `hospitalization_id`
- `admission_dttm`
- `discharge_dttm`
- `age_at_admission`
- `county_code`
- `zipcode_five_digit`
- `zipcode_nine_digit`
- `census_tract`

Used for:

- cohort construction
- geographic exposure assignment
- outcome ascertainment

---

# 3. `adt`

Required fields:

- `hospitalization_id`
- `in_dttm`
- `out_dttm`
- `location_category`

Used to:

- identify ICU stays
- determine ICU admission time
- define the **index time (t0)** for trajectory reconstruction.

---

# 4. `respiratory_support`

Required fields:

- `hospitalization_id`
- `recorded_dttm`
- `device_category`
- `mode_category`
- `fio2_set`

Optional fields:

- `peep_set`
- `resp_rate_set`
- `tidal_volume_set`
- `plateau_pressure_obs`
- `peak_inspiratory_pressure_obs`

Used to:

- determine respiratory support state
- calculate FiO₂ exposure
- define ventilatory trajectories

Respiratory support states are collapsed into the following tiers:

- ROOM AIR
- LOW_O2
- NIV
- IMV
- OTHER

---

# 5. `vitals`

Required fields:

- `hospitalization_id`
- `recorded_dttm`
- `vital_category`
- `vital_value`

Used to extract:

- SpO₂

These measurements are paired with FiO₂ to detect **hypoxemic respiratory failure**.

---

# 6. `labs`

Required fields:

- `hospitalization_id`
- `lab_result_dttm`
- `lab_category`
- `lab_value_numeric`

Used to extract:

- PaO₂
- PaCO₂
- arterial pH

These measurements define physiologic respiratory failure subtypes.

---

# 7. `medication_admin_continuous`

Required fields:

- `hospitalization_id`
- `admin_dttm`
- `med_group`
- `mar_action_group`

Used to detect **vasoactive medication administration**, which indicates circulatory failure.

---

# 8. `patient_assessments`

Required fields:

- `hospitalization_id`
- `recorded_dttm`
- `assessment_category`
- `numerical_value`

Used to extract:

- Glasgow Coma Scale (GCS)

These data contribute to **SOFA score calculation**.

---

# 9. `hospital_diagnosis`

Required fields:

- `hospitalization_id`
- `diagnosis_code`
- `diagnosis_code_format`
- `poa_present`

Used to derive:

- Charlson comorbidity score
- metastatic cancer indicators
- neurologic comorbidity flags

---

# Cohort Identification

Eligible hospitalizations must meet:

1. Age ≥ 18 years  
2. ICU admission recorded in the `adt` table  
3. Available respiratory support data  
4. Available geographic identifier for exposure assignment  

Only the **first ICU admission per hospitalization** is used.

---

# Index Time (t0)

- t0 = first ICU admission time


All physiologic trajectories are anchored relative to this time.

---

# Trajectory Construction

Hourly trajectories are constructed for **0–72 hours following ICU admission**.

Each hour is assigned a composite physiologic state composed of:

### Respiratory Support Tier

Derived from `respiratory_support.device_category`.

### Acute Respiratory Failure Subtype

Defined using physiologic criteria.

Hypoxemic ARF:

- SpO₂ < 90% on room air  
- PaO₂ ≤ 60 mmHg on room air  
- PaO₂ / FiO₂ ≤ 300  

Hypercapnic ARF:

- PaCO₂ ≥ 45 mmHg **AND** pH < 7.35

Combined states:

- `ARF_HYPOX`
- `ARF_HYPER`
- `ARF_MIXED`
- `NO_ARF`

### Vasopressor Exposure

Binary indicator of vasoactive medication administration.

---

# Sequence Clustering

Hourly states are combined into patient-level sequences.

Sequence clustering is performed using:

- **Optimal Matching distance**
- **Transition-rate substitution matrix**
- **Ward hierarchical clustering**

The optimal number of clusters is determined using:

- silhouette width
- cluster interpretability

Clusters represent distinct **respiratory failure trajectory phenotypes**.

---

# Environmental Exposure Variables

Long-term ambient exposure estimates are assigned using patient geography.

Primary exposures:

- **PM₂.₅** — long-term fine particulate matter exposure
- **NO₂** — long-term nitrogen dioxide exposure

Exposure values are standardized (z-scores) prior to modeling.

---

# Outcomes

Primary outcome:

- In-hospital mortality or hospice discharge

Secondary outcomes:

- ICU length of stay
- Need for invasive mechanical ventilation
- Timing of mechanical ventilation
- Respiratory failure severity

---

# Statistical Modeling

Each site performs the following models locally.

### Trajectory Phenotype Model
 
Multinomial logistic regression: 
All physiologic trajectories are anchored relative to this time.

---

# Trajectory Construction

Hourly trajectories are constructed for **0–72 hours following ICU admission**.

Each hour is assigned a composite physiologic state composed of:

### Respiratory Support Tier

Derived from `respiratory_support.device_category`.

### Acute Respiratory Failure Subtype

Defined using physiologic criteria.

Hypoxemic ARF:

- SpO₂ < 90% on room air  
- PaO₂ ≤ 60 mmHg on room air  
- PaO₂ / FiO₂ ≤ 300  

Hypercapnic ARF:

- PaCO₂ ≥ 45 mmHg **AND** pH < 7.35

Combined states:

- `ARF_HYPOX`
- `ARF_HYPER`
- `ARF_MIXED`
- `NO_ARF`

### Vasopressor Exposure

Binary indicator of vasoactive medication administration.

---

# Sequence Clustering

Hourly states are combined into patient-level sequences.

Sequence clustering is performed using:

- **Optimal Matching distance**
- **Transition-rate substitution matrix**
- **Ward hierarchical clustering**

The optimal number of clusters is determined using:

- silhouette width
- cluster interpretability

Clusters represent distinct **respiratory failure trajectory phenotypes**.

---

# Environmental Exposure Variables

Long-term ambient exposure estimates are assigned using patient geography.

Primary exposures:

- **PM₂.₅** — long-term fine particulate matter exposure
- **NO₂** — long-term nitrogen dioxide exposure

Exposure values are standardized (z-scores) prior to modeling.

---

# Outcomes

Primary outcome:

- In-hospital mortality or hospice discharge

Secondary outcomes:

- ICU length of stay
- Need for invasive mechanical ventilation
- Timing of mechanical ventilation
- Respiratory failure severity

---

# Statistical Modeling

Each site performs the following models locally.

### Trajectory Phenotype Model

Multinomial logistic regression:

All physiologic trajectories are anchored relative to this time.

---

# Trajectory Construction

Hourly trajectories are constructed for **0–72 hours following ICU admission**.

Each hour is assigned a composite physiologic state composed of:

### Respiratory Support Tier

Derived from `respiratory_support.device_category`.

### Acute Respiratory Failure Subtype

Defined using physiologic criteria.

Hypoxemic ARF:

- SpO₂ < 90% on room air  
- PaO₂ ≤ 60 mmHg on room air  
- PaO₂ / FiO₂ ≤ 300  

Hypercapnic ARF:

- PaCO₂ ≥ 45 mmHg **AND** pH < 7.35

Combined states:

- `ARF_HYPOX`
- `ARF_HYPER`
- `ARF_MIXED`
- `NO_ARF`

### Vasopressor Exposure

Binary indicator of vasoactive medication administration.

---

# Sequence Clustering

Hourly states are combined into patient-level sequences.

Sequence clustering is performed using:

- **Optimal Matching distance**
- **Transition-rate substitution matrix**
- **Ward hierarchical clustering**

The optimal number of clusters is determined using:

- silhouette width
- cluster interpretability

Clusters represent distinct **respiratory failure trajectory phenotypes**.

---

# Environmental Exposure Variables

Long-term ambient exposure estimates are assigned using patient geography.

Primary exposures:

- **PM₂.₅** — long-term fine particulate matter exposure
- **NO₂** — long-term nitrogen dioxide exposure

Exposure values are standardized (z-scores) prior to modeling.

---

# Outcomes

Primary outcome:

- In-hospital mortality or hospice discharge

Secondary outcomes:

- ICU length of stay
- Need for invasive mechanical ventilation
- Timing of mechanical ventilation
- Respiratory failure severity

---

# Statistical Modeling

Each site performs the following models locally.

### Trajectory Phenotype Model

Multinomial logistic regression: trajectory_cluster ~ PM2.5 + NO2 + age + sex + race

---

### Mortality Model

Logistic regression: death_or_hospice ~ PM2.5 + NO2 + trajectory_cluster + covariates

---

### Interaction Models

death_or_hospice ~ NO2 * trajectory_cluster


These models test whether pollution exposure modifies mortality risk differently across trajectory phenotypes.

---

### Landmark Survival Model

A **72-hour landmark survival analysis** evaluates mortality risk beyond the initial trajectory window.

---

# Federated Analysis Framework

Because patient-level data cannot be shared across sites, the project uses a **federated analysis approach**.

Each site executes the full pipeline locally and exports only:

- cluster summaries
- trajectory signatures
- regression coefficients
- model diagnostics

No patient-level data are shared.

Central analysis then performs:

- cross-site cluster harmonization
- pooled meta-analysis of exposure associations

---

# Project Outputs

Each site produces PHI-safe outputs including:

### Trajectory Outputs

- trajectory cluster assignments
- representative sequence plots
- cluster signature plots

### Cluster Phenotype Summaries

- cluster severity
- SOFA scores
- comorbidity burden
- vasopressor exposure

### Exposure Association Results

- trajectory association models
- mortality models
- interaction models

### Figures

- trajectory signatures
- predicted mortality by exposure
- cluster severity profiles

Final outputs are saved in: output/final

---

# Running the Project

## 1 Update Configuration

Edit: config/config.json

Follow instructions in: config/README.md

---

## 2 Run the Cohort Identification App

Run: 01_luncx_cohort.R

This interactive tool generates the configuration file and verifies table availability.

---

## 3 Execute the Analysis Pipeline

Run the full analysis pipeline: 02_federated_clusters.R

This script performs:

- trajectory construction
- sequence clustering
- exposure modeling
- figure generation
- federated export creation







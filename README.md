# *CLIF Post–Lung Cancer Resection Respiratory Recovery Project (PLCR3)*

## CLIF VERSION 

2.1

---

## Objective

To evaluate whether long-term ambient air pollution exposure (PM₂.₅ and NO₂) is associated with postoperative respiratory recovery severity among patients undergoing lung cancer resection who require ICU care.

This project integrates high-resolution CLIF ICU physiologic data with preoperative environmental exposure estimates to quantify:

- Respiratory support at ICU admission
- Early ventilatory trajectory and escalation within 24 hours
- Duration of invasive mechanical ventilation
- Development of postoperative pulmonary complications (PPCs)
- ICU length of stay
- In-hospital mortality

We hypothesize that greater chronic pollution burden is associated with more severe postoperative respiratory recovery.

---

## Required CLIF Tables and Fields

Please refer to the [CLIF data dictionary](https://clif-icu.com/data-dictionary), [CLIF Tools](https://clif-icu.com/tools), [ETL Guide](https://clif-icu.com/etl-guide), and [specific table contacts](https://github.com/clif-consortium/CLIF?tab=readme-ov-file#relational-clif) for more information on constructing the required tables and fields.

The following tables are required:

### 1. **patient**
- `patient_id`
- `birth_date`
- `sex_category`
- `race_category`
- `ethnicity_category`

Used to determine age and demographic covariates.

---

### 2. **hospitalization**
- `patient_id`
- `hospitalization_id`
- `admission_dttm`
- `discharge_dttm`
- `age_at_admission`
- `zipcode_five_digit`
- `zipcode_nine_digit`
- `census_tract`
- `county_code`

Used for:
- Study period restriction
- Age calculation
- Geographic linkage for air pollution exposure assignment

---

### 3. **adt**
- `hospitalization_id`
- `in_dttm`
- `out_dttm`
- `location_category`

Used to:
- Identify ICU admissions
- Determine ICU admission occurring after lung resection
- Define ICU length of stay
- Anchor index time (t0)

---

### 4. **hospital_diagnosis**
- `hospitalization_id`
- `icd_code`
- `poa_present`
- `diagnosis_dttm`

Used to:
- Identify lung cancer diagnosis
- Require present-on-admission flag (`poa_present == 1`)
- Ensure cancer predates hospitalization

---

### 5. **patient_procedures**
- `hospitalization_id`
- `procedure_code`
- `procedure_code_type`
- `procedure_billed_dttm`
- `procedure_dttm`

Used to:
- Identify lung cancer resection procedures (lobectomy, segmentectomy, pneumonectomy, wedge resection)
- Establish temporal sequence (procedure must occur before ICU admission)

---

### 6. **respiratory_support**
- `hospitalization_id`
- `recorded_dttm`
- `device_category`
- `mode_category`
- `fio2_set`
- `peep_set`
- `resp_rate_set`
- `tidal_volume_set`
- `plateau_pressure_obs`
- `peak_inspiratory_pressure_obs`

Used to:
- Characterize respiratory support at ICU admission
- Quantify escalation within 24 hours
- Model ventilatory trajectory
- Calculate duration of invasive mechanical ventilation

---

### 7. (Recommended) **labs**
- `hospitalization_id`
- `lab_result_dttm`
- `lab_category`
- `lab_value`

Recommended categories include:
- PaO2
- Lactate
- Creatinine
- WBC

Used to define postoperative pulmonary complications and physiologic severity.

---

## Cohort Identification

### Inclusion Criteria

A hospitalization is included if all of the following are met:

1. Age ≥ 18 years at admission  
2. Lung cancer diagnosis present on admission  
   - ICD-9-CM or ICD-10-CM lung cancer codes  
   - `poa_present == 1`  
3. Lung resection procedure during hospitalization  
4. ICU admission occurring after the procedure  
   - `procedure_time < icu_in_time`  
   - Closest ICU admission after resection selected  
   - Optional restriction: ICU admission within 48 hours of surgery  
5. Valid geographic identifier available for exposure linkage  

### Index Time (t0)

`t0 = ICU admission time following lung resection`

All recovery and complication metrics are anchored to this time.

---

## Outcome Framework

The primary construct is **Postoperative Respiratory Recovery Severity**, defined using:

- Respiratory support category at ICU admission
- Escalation to invasive mechanical ventilation within 24 hours
- Duration of invasive ventilation
- Development of acute respiratory failure
- Reintubation
- Postoperative pneumonia
- ICU length of stay
- In-hospital mortality

These outcomes represent a continuum of postoperative respiratory vulnerability.

---

## Exposure Variables

Preoperative long-term ambient exposure estimates:

- Annual average PM₂.₅
- Annual average NO₂

Exposure is assigned using available geographic identifiers (census tract, ZIP, or county).

---

## Expected Results

The project will produce:

1. A curated analytic cohort of post–lung cancer resection ICU admissions  
2. Summary distributions of respiratory support at ICU entry  
3. Ventilatory escalation and trajectory analyses  
4. Incidence of postoperative pulmonary complications  
5. Multivariable models estimating associations between air pollution exposure and:
   - Ventilation duration  
   - Escalation risk  
   - PPC incidence  
   - ICU length of stay  
   - Mortality  

Final project results should be saved in the [`output/final`](output/README.md) directory.

---

## Detailed Instructions for Running the Project

### 1. Update `config/config.json`

Follow instructions in the [config/README.md](config/README.md) file for detailed configuration steps.

**Note: if using the `01_run_cohort_id_app.R` file, this step is not necessary as the app will create the config file for the user.**

---

### 2. Set Up the Project Environment

#### R

Run:


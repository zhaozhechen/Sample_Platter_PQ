# Sampler Platter transfer-entropy results

## Methods

Each analysis tests lags 0 through 72 in the dataset's native observation steps. `Lag_time` and every plot x-axis convert those native steps to hours using the modal timestamp interval. Thus the physical horizon varies with native resolution.

The generalized workflow preserves the project histogram estimator, quantile-folded extremes, optional structural-zero bins, and 300 shuffled null realizations. Lag windows are retained only when timestamps are contiguous at the detected native resolution.

Temperature analyses use five-day same-clock anomalies, consistent with the prior project workflow. Hubbard Brook growing-season runs retain April through October. Casper is omitted as requested. No stream-temperature-to-discharge direction is analyzed.

Cumulative TE is defined here exactly as requested: the sum of TE values that exceed their shuffled TE threshold, divided by 72 configured lag steps. The lag-zero position can contribute to the numerator, but the fixed denominator remains 72.

## Results summary

| Site | Variable pair | Native resolution | N | Start | End | Peak TE | Peak lag (steps and hours) | Peak significant | Cumulative TE | Notes |
|---|---|---:|---:|---|---|---:|---:|---:|---:|---|
| Angelo DRY | Air temperature anomaly -> Stream temperature anomaly | 60 min | 40592 | 2015-12-03 | 2020-07-20 | 0.060601 | 1 (1 h) | yes | 0.0037244 | Same hourly file and anomaly choice as prior v3 analysis |
| Angelo DRY | Discharge -> Stream temperature anomaly | 60 min | 40603 | 2015-12-03 | 2020-07-20 | 0.010592 | 1 (1 h) | no | 0 | Same hourly file and sink anomaly choice as prior v3 analysis No TE value exceeded the shuffled threshold. |
| Angelo ELDER | Air temperature anomaly -> Stream temperature anomaly | 30 min | 172456 | 2013-05-10 | 2023-04-24 | 0.03759 | 2 (1 h) | yes | 0.005885 | Same Elder_subdaily file and anomaly choice as prior v3 analysis |
| Angelo ELDER | Discharge -> Stream temperature anomaly | 30 min | 172798 | 2013-05-10 | 2023-04-24 | 0.0010205 | 46 (23 h) | no | 0 | Same Elder_subdaily file and sink anomaly choice as prior v3 analysis No TE value exceeded the shuffled threshold. |
| Angelo ELDER | Precipitation -> Discharge | 30 min | 172800 | 2013-05-08 | 2023-04-26 | 0.0052341 | 3 (1.5 h) | yes | 0.0013967 | Same Elder_subdaily file as prior v3 analysis |
| Atlanta No Business | Precipitation -> Discharge | 15 min | 298268 | 2009-10-01 | 2025-09-02 | 0.0055075 | 2 (0.5 h) | yes | 0.0012122 | Three-site split and P-Q file retained from the 2026-07-10 update |
| Atlanta Richland Creek | Precipitation -> Discharge | 15 min | 323921 | 2015-10-01 | 2025-09-02 | 0.0062358 | 4 (1 h) | yes | 0.0015133 | Three-site split and P-Q file retained from the 2026-07-10 update |
| Atlanta Crooked Creek | Precipitation -> Discharge | 15 min | 319431 | 2015-10-01 | 2025-09-02 | 0.012469 | 4 (1 h) | yes | 0.0019074 | Three-site split and P-Q file retained from the 2026-07-10 update |
| Atlanta No Business | Discharge -> Stream temperature anomaly | 15 min | 517247 | 2009-10-03 | 2025-09-29 | 0.00044719 | 11 (2.75 h) | no | 0 | Q and temperature merged by UTC timestamp; no reverse temperature-to-Q analysis No TE value exceeded the shuffled threshold. |
| Atlanta Richland Creek | Discharge -> Stream temperature anomaly | 15 min | 543940 | 2009-10-03 | 2025-09-29 | 0.0006679 | 0 (0 h) | no | 0 | Q and temperature merged by UTC timestamp; no reverse temperature-to-Q analysis No TE value exceeded the shuffled threshold. |
| Atlanta Crooked Creek | Discharge -> Stream temperature anomaly | 15 min | 495499 | 2009-10-03 | 2025-09-29 | 0.001842 | 0 (0 h) | yes | 0.000046769 | Q and temperature merged by UTC timestamp; no reverse temperature-to-Q analysis |
| Atlanta No Business | Precipitation -> Specific conductivity | 15 min | 500869 | 2009-10-01 | 2025-10-01 | 0.0094219 | 2 (0.5 h) | yes | 0.0015295 | P and conductivity from the combined Atlanta Q-P-conductivity file |
| Atlanta Richland Creek | Precipitation -> Specific conductivity | 15 min | 514417 | 2009-10-01 | 2025-10-01 | 0.010084 | 1 (0.25 h) | yes | 0.0015911 | P and conductivity from the combined Atlanta Q-P-conductivity file |
| Atlanta Crooked Creek | Precipitation -> Specific conductivity | 15 min | 520925 | 2009-10-01 | 2025-10-01 | 0.012406 | 5 (1.25 h) | yes | 0.0029776 | P and conductivity from the combined Atlanta Q-P-conductivity file |
| Atlanta No Business | Specific conductivity -> Stream temperature anomaly | 15 min | 538378 | 2009-10-03 | 2025-09-29 | 0.0015084 | 22 (5.5 h) | no | 0 | Conductivity and temperature merged by UTC timestamp No TE value exceeded the shuffled threshold. |
| Atlanta Richland Creek | Specific conductivity -> Stream temperature anomaly | 15 min | 540483 | 2009-10-03 | 2025-09-29 | 0.00098603 | 0 (0 h) | no | 0 | Conductivity and temperature merged by UTC timestamp No TE value exceeded the shuffled threshold. |
| Atlanta Crooked Creek | Specific conductivity -> Stream temperature anomaly | 15 min | 537370 | 2009-10-03 | 2025-09-29 | 0.0013706 | 0 (0 h) | no | 0 | Conductivity and temperature merged by UTC timestamp No TE value exceeded the shuffled threshold. |
| Baltimore | Precipitation -> Discharge | 5 min | 852425 | 2009-05-11 | 2017-12-22 | 0.0034606 | 18 (1.5 h) | yes | 0.0011057 | Same BES_USGS_QP_0917 file as prior v3 analysis |
| Hubbard Brook W3 (all series) | Precipitation -> Discharge | 15 min | 432998 | 2011-08-26 | 2023-12-31 | 0.004545 | 0 (0 h) | yes | 0.00096091 | Full-series version requested; same file as prior v3 analysis |
| Hubbard Brook W3 (growing season) | Precipitation -> Discharge | 15 min | 252902 | 2011-08-26 | 2023-10-31 | 0.0050174 | 0 (0 h) | yes | 0.0010372 | April through October; same file as prior v3 analysis |
| Hubbard Brook W8 (all series) | Precipitation -> Discharge | 15 min | 1919853 | 1968-05-01 | 2023-12-31 | 0.0035481 | 2 (0.5 h) | yes | 0.00091778 | Full-series version requested; same file as prior v3 metadata |
| Hubbard Brook W8 (growing season) | Precipitation -> Discharge | 15 min | 1123147 | 1968-05-01 | 2023-10-31 | 0.0040209 | 2 (0.5 h) | yes | 0.00095008 | April through October; same file as prior v3 metadata |
| Konza | Precipitation -> Discharge | 15 min | 516106 | 2010-03-19 | 2024-12-31 | 0.0016572 | 0 (0 h) | yes | 0.0003018 | Updated interpolated-Q file specified in the 2026-07-10 analysis log |
| Starkey | Precipitation -> Discharge | 15 min | 425068 | 2011-07-11 | 2024-05-01 | 0.00023382 | 29 (7.25 h) | no | 0 | Same merged P-Q file used in the prior Starkey v3 analysis No TE value exceeded the shuffled threshold. |
| H.J. Andrews GSWS02 | Precipitation -> Discharge | 30 min | 397989 | 1957-10-01 | 2018-10-01 | 0.0081952 | 1 (0.5 h) | yes | 0.0016199 | GSWS02 was the H.J. Andrews reference site in the earlier common P-Q selection; new Keith processed file |

## Analysis details

### Angelo_DRY_Tair_Tstream

- Site and pair: Angelo DRY; Air temperature anomaly -> Stream temperature anomaly
- Record: 40592 complete observations from 2015-12-03 through 2020-07-20 23:00:00.
- Native resolution: 60 minutes; tested lags: 0-72 steps (0-72 hours).
- Peak TE: 0.060601 at 1 steps (1 hours); above shuffled threshold: yes.
- Significant lag positions: 6; cumulative TE: 0.0037244.
- Notes: Same hourly file and anomaly choice as prior v3 analysis
- [Full lag table](../04_Results/Tables/TE_df_Angelo_DRY_Tair_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Angelo_DRY_Tair_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Angelo_DRY_Tair_Tstream.pdf)

### Angelo_DRY_Q_Tstream

- Site and pair: Angelo DRY; Discharge -> Stream temperature anomaly
- Record: 40603 complete observations from 2015-12-03 through 2020-07-20 23:00:00.
- Native resolution: 60 minutes; tested lags: 0-72 steps (0-72 hours).
- Peak TE: 0.010592 at 1 steps (1 hours); above shuffled threshold: no.
- Significant lag positions: 0; cumulative TE: 0.
- Notes: Same hourly file and sink anomaly choice as prior v3 analysis No TE value exceeded the shuffled threshold.
- [Full lag table](../04_Results/Tables/TE_df_Angelo_DRY_Q_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Angelo_DRY_Q_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Angelo_DRY_Q_Tstream.pdf)

### Angelo_ELDER_Tair_Tstream

- Site and pair: Angelo ELDER; Air temperature anomaly -> Stream temperature anomaly
- Record: 172456 complete observations from 2013-05-10 through 2023-04-24 21:00:00.
- Native resolution: 30 minutes; tested lags: 0-72 steps (0-36 hours).
- Peak TE: 0.03759 at 2 steps (1 hours); above shuffled threshold: yes.
- Significant lag positions: 20; cumulative TE: 0.005885.
- Notes: Same Elder_subdaily file and anomaly choice as prior v3 analysis
- [Full lag table](../04_Results/Tables/TE_df_Angelo_ELDER_Tair_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Angelo_ELDER_Tair_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Angelo_ELDER_Tair_Tstream.pdf)

### Angelo_ELDER_Q_Tstream

- Site and pair: Angelo ELDER; Discharge -> Stream temperature anomaly
- Record: 172798 complete observations from 2013-05-10 through 2023-04-24 21:00:00.
- Native resolution: 30 minutes; tested lags: 0-72 steps (0-36 hours).
- Peak TE: 0.0010205 at 46 steps (23 hours); above shuffled threshold: no.
- Significant lag positions: 0; cumulative TE: 0.
- Notes: Same Elder_subdaily file and sink anomaly choice as prior v3 analysis No TE value exceeded the shuffled threshold.
- [Full lag table](../04_Results/Tables/TE_df_Angelo_ELDER_Q_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Angelo_ELDER_Q_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Angelo_ELDER_Q_Tstream.pdf)

### Angelo_ELDER_P_Q

- Site and pair: Angelo ELDER; Precipitation -> Discharge
- Record: 172800 complete observations from 2013-05-08 through 2023-04-26 21:00:00.
- Native resolution: 30 minutes; tested lags: 0-72 steps (0-36 hours).
- Peak TE: 0.0052341 at 3 steps (1.5 hours); above shuffled threshold: yes.
- Significant lag positions: 46; cumulative TE: 0.0013967.
- Notes: Same Elder_subdaily file as prior v3 analysis
- [Full lag table](../04_Results/Tables/TE_df_Angelo_ELDER_P_Q.csv) | [PNG figure](../04_Results/Figures/TE_lag_Angelo_ELDER_P_Q.png) | [PDF figure](../04_Results/Figures/TE_lag_Angelo_ELDER_P_Q.pdf)

### Atlanta_NB_P_Q

- Site and pair: Atlanta No Business; Precipitation -> Discharge
- Record: 298268 complete observations from 2009-10-01 04:00:00 through 2025-09-02 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0055075 at 2 steps (0.5 hours); above shuffled threshold: yes.
- Significant lag positions: 73; cumulative TE: 0.0012122.
- Notes: Three-site split and P-Q file retained from the 2026-07-10 update
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_NB_P_Q.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_NB_P_Q.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_NB_P_Q.pdf)

### Atlanta_RC_P_Q

- Site and pair: Atlanta Richland Creek; Precipitation -> Discharge
- Record: 323921 complete observations from 2015-10-01 04:00:00 through 2025-09-02 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0062358 at 4 steps (1 hours); above shuffled threshold: yes.
- Significant lag positions: 73; cumulative TE: 0.0015133.
- Notes: Three-site split and P-Q file retained from the 2026-07-10 update
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_RC_P_Q.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_RC_P_Q.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_RC_P_Q.pdf)

### Atlanta_CC_P_Q

- Site and pair: Atlanta Crooked Creek; Precipitation -> Discharge
- Record: 319431 complete observations from 2015-10-01 04:00:00 through 2025-09-02 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.012469 at 4 steps (1 hours); above shuffled threshold: yes.
- Significant lag positions: 73; cumulative TE: 0.0019074.
- Notes: Three-site split and P-Q file retained from the 2026-07-10 update
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_CC_P_Q.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_CC_P_Q.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_CC_P_Q.pdf)

### Atlanta_NB_Q_Tstream

- Site and pair: Atlanta No Business; Discharge -> Stream temperature anomaly
- Record: 517247 complete observations from 2009-10-03 04:00:00 through 2025-09-29 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.00044719 at 11 steps (2.75 hours); above shuffled threshold: no.
- Significant lag positions: 0; cumulative TE: 0.
- Notes: Q and temperature merged by UTC timestamp; no reverse temperature-to-Q analysis No TE value exceeded the shuffled threshold.
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_NB_Q_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_NB_Q_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_NB_Q_Tstream.pdf)

### Atlanta_RC_Q_Tstream

- Site and pair: Atlanta Richland Creek; Discharge -> Stream temperature anomaly
- Record: 543940 complete observations from 2009-10-03 04:00:00 through 2025-09-29 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0006679 at 0 steps (0 hours); above shuffled threshold: no.
- Significant lag positions: 0; cumulative TE: 0.
- Notes: Q and temperature merged by UTC timestamp; no reverse temperature-to-Q analysis No TE value exceeded the shuffled threshold.
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_RC_Q_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_RC_Q_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_RC_Q_Tstream.pdf)

### Atlanta_CC_Q_Tstream

- Site and pair: Atlanta Crooked Creek; Discharge -> Stream temperature anomaly
- Record: 495499 complete observations from 2009-10-03 04:00:00 through 2025-09-29 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.001842 at 0 steps (0 hours); above shuffled threshold: yes.
- Significant lag positions: 2; cumulative TE: 0.000046769.
- Notes: Q and temperature merged by UTC timestamp; no reverse temperature-to-Q analysis
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_CC_Q_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_CC_Q_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_CC_Q_Tstream.pdf)

### Atlanta_NB_P_Cond

- Site and pair: Atlanta No Business; Precipitation -> Specific conductivity
- Record: 500869 complete observations from 2009-10-01 04:00:00 through 2025-10-01 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0094219 at 2 steps (0.5 hours); above shuffled threshold: yes.
- Significant lag positions: 26; cumulative TE: 0.0015295.
- Notes: P and conductivity from the combined Atlanta Q-P-conductivity file
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_NB_P_Cond.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_NB_P_Cond.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_NB_P_Cond.pdf)

### Atlanta_RC_P_Cond

- Site and pair: Atlanta Richland Creek; Precipitation -> Specific conductivity
- Record: 514417 complete observations from 2009-10-01 04:00:00 through 2025-10-01 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.010084 at 1 steps (0.25 hours); above shuffled threshold: yes.
- Significant lag positions: 38; cumulative TE: 0.0015911.
- Notes: P and conductivity from the combined Atlanta Q-P-conductivity file
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_RC_P_Cond.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_RC_P_Cond.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_RC_P_Cond.pdf)

### Atlanta_CC_P_Cond

- Site and pair: Atlanta Crooked Creek; Precipitation -> Specific conductivity
- Record: 520925 complete observations from 2009-10-01 05:15:00 through 2025-10-01 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.012406 at 5 steps (1.25 hours); above shuffled threshold: yes.
- Significant lag positions: 36; cumulative TE: 0.0029776.
- Notes: P and conductivity from the combined Atlanta Q-P-conductivity file
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_CC_P_Cond.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_CC_P_Cond.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_CC_P_Cond.pdf)

### Atlanta_NB_Cond_Tstream

- Site and pair: Atlanta No Business; Specific conductivity -> Stream temperature anomaly
- Record: 538378 complete observations from 2009-10-03 04:00:00 through 2025-09-29 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0015084 at 22 steps (5.5 hours); above shuffled threshold: no.
- Significant lag positions: 0; cumulative TE: 0.
- Notes: Conductivity and temperature merged by UTC timestamp No TE value exceeded the shuffled threshold.
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_NB_Cond_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_NB_Cond_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_NB_Cond_Tstream.pdf)

### Atlanta_RC_Cond_Tstream

- Site and pair: Atlanta Richland Creek; Specific conductivity -> Stream temperature anomaly
- Record: 540483 complete observations from 2009-10-03 04:00:00 through 2025-09-29 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.00098603 at 0 steps (0 hours); above shuffled threshold: no.
- Significant lag positions: 0; cumulative TE: 0.
- Notes: Conductivity and temperature merged by UTC timestamp No TE value exceeded the shuffled threshold.
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_RC_Cond_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_RC_Cond_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_RC_Cond_Tstream.pdf)

### Atlanta_CC_Cond_Tstream

- Site and pair: Atlanta Crooked Creek; Specific conductivity -> Stream temperature anomaly
- Record: 537370 complete observations from 2009-10-03 04:00:00 through 2025-09-29 03:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0013706 at 0 steps (0 hours); above shuffled threshold: no.
- Significant lag positions: 0; cumulative TE: 0.
- Notes: Conductivity and temperature merged by UTC timestamp No TE value exceeded the shuffled threshold.
- [Full lag table](../04_Results/Tables/TE_df_Atlanta_CC_Cond_Tstream.csv) | [PNG figure](../04_Results/Figures/TE_lag_Atlanta_CC_Cond_Tstream.png) | [PDF figure](../04_Results/Figures/TE_lag_Atlanta_CC_Cond_Tstream.pdf)

### Baltimore_P_Q

- Site and pair: Baltimore; Precipitation -> Discharge
- Record: 852425 complete observations from 2009-05-11 23:00:00 through 2017-12-22 23:55:00.
- Native resolution: 5 minutes; tested lags: 0-72 steps (0-6 hours).
- Peak TE: 0.0034606 at 18 steps (1.5 hours); above shuffled threshold: yes.
- Significant lag positions: 73; cumulative TE: 0.0011057.
- Notes: Same BES_USGS_QP_0917 file as prior v3 analysis
- [Full lag table](../04_Results/Tables/TE_df_Baltimore_P_Q.csv) | [PNG figure](../04_Results/Figures/TE_lag_Baltimore_P_Q.png) | [PDF figure](../04_Results/Figures/TE_lag_Baltimore_P_Q.pdf)

### HB_w3_P_Q_full

- Site and pair: Hubbard Brook W3 (all series); Precipitation -> Discharge
- Record: 432998 complete observations from 2011-08-26 14:30:00 through 2023-12-31 23:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.004545 at 0 steps (0 hours); above shuffled threshold: yes.
- Significant lag positions: 72; cumulative TE: 0.00096091.
- Notes: Full-series version requested; same file as prior v3 analysis
- [Full lag table](../04_Results/Tables/TE_df_HB_w3_P_Q_full.csv) | [PNG figure](../04_Results/Figures/TE_lag_HB_w3_P_Q_full.png) | [PDF figure](../04_Results/Figures/TE_lag_HB_w3_P_Q_full.pdf)

### HB_w3_P_Q_GS

- Site and pair: Hubbard Brook W3 (growing season); Precipitation -> Discharge
- Record: 252902 complete observations from 2011-08-26 14:30:00 through 2023-10-31 23:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0050174 at 0 steps (0 hours); above shuffled threshold: yes.
- Significant lag positions: 63; cumulative TE: 0.0010372.
- Notes: April through October; same file as prior v3 analysis
- [Full lag table](../04_Results/Tables/TE_df_HB_w3_P_Q_GS.csv) | [PNG figure](../04_Results/Figures/TE_lag_HB_w3_P_Q_GS.png) | [PDF figure](../04_Results/Figures/TE_lag_HB_w3_P_Q_GS.pdf)

### HB_w8_P_Q_full

- Site and pair: Hubbard Brook W8 (all series); Precipitation -> Discharge
- Record: 1919853 complete observations from 1968-05-01 05:45:00 through 2023-12-31 23:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0035481 at 2 steps (0.5 hours); above shuffled threshold: yes.
- Significant lag positions: 73; cumulative TE: 0.00091778.
- Notes: Full-series version requested; same file as prior v3 metadata
- [Full lag table](../04_Results/Tables/TE_df_HB_w8_P_Q_full.csv) | [PNG figure](../04_Results/Figures/TE_lag_HB_w8_P_Q_full.png) | [PDF figure](../04_Results/Figures/TE_lag_HB_w8_P_Q_full.pdf)

### HB_w8_P_Q_GS

- Site and pair: Hubbard Brook W8 (growing season); Precipitation -> Discharge
- Record: 1123147 complete observations from 1968-05-01 05:45:00 through 2023-10-31 23:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0040209 at 2 steps (0.5 hours); above shuffled threshold: yes.
- Significant lag positions: 73; cumulative TE: 0.00095008.
- Notes: April through October; same file as prior v3 metadata
- [Full lag table](../04_Results/Tables/TE_df_HB_w8_P_Q_GS.csv) | [PNG figure](../04_Results/Figures/TE_lag_HB_w8_P_Q_GS.png) | [PDF figure](../04_Results/Figures/TE_lag_HB_w8_P_Q_GS.pdf)

### Konza_P_Q

- Site and pair: Konza; Precipitation -> Discharge
- Record: 516106 complete observations from 2010-03-19 17:45:00 through 2024-12-31 00:45:00.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.0016572 at 0 steps (0 hours); above shuffled threshold: yes.
- Significant lag positions: 62; cumulative TE: 0.0003018.
- Notes: Updated interpolated-Q file specified in the 2026-07-10 analysis log
- [Full lag table](../04_Results/Tables/TE_df_Konza_P_Q.csv) | [PNG figure](../04_Results/Figures/TE_lag_Konza_P_Q.png) | [PDF figure](../04_Results/Figures/TE_lag_Konza_P_Q.pdf)

### Starkey_P_Q

- Site and pair: Starkey; Precipitation -> Discharge
- Record: 425068 complete observations from 2011-07-11 23:30:00 through 2024-05-01.
- Native resolution: 15 minutes; tested lags: 0-72 steps (0-18 hours).
- Peak TE: 0.00023382 at 29 steps (7.25 hours); above shuffled threshold: no.
- Significant lag positions: 0; cumulative TE: 0.
- Notes: Same merged P-Q file used in the prior Starkey v3 analysis No TE value exceeded the shuffled threshold.
- [Full lag table](../04_Results/Tables/TE_df_Starkey_P_Q.csv) | [PNG figure](../04_Results/Figures/TE_lag_Starkey_P_Q.png) | [PDF figure](../04_Results/Figures/TE_lag_Starkey_P_Q.pdf)

### HJA_GSWS02_P_Q

- Site and pair: H.J. Andrews GSWS02; Precipitation -> Discharge
- Record: 397989 complete observations from 1957-10-01 through 2018-10-01.
- Native resolution: 30 minutes; tested lags: 0-72 steps (0-36 hours).
- Peak TE: 0.0081952 at 1 steps (0.5 hours); above shuffled threshold: yes.
- Significant lag positions: 65; cumulative TE: 0.0016199.
- Notes: GSWS02 was the H.J. Andrews reference site in the earlier common P-Q selection; new Keith processed file
- [Full lag table](../04_Results/Tables/TE_df_HJA_GSWS02_P_Q.csv) | [PNG figure](../04_Results/Figures/TE_lag_HJA_GSWS02_P_Q.png) | [PDF figure](../04_Results/Figures/TE_lag_HJA_GSWS02_P_Q.pdf)

## Interpretation notes

A peak at lag zero describes same-timestamp predictive information under this estimator and should not be interpreted as a delayed causal response. A peak marked as not significant did not exceed the shuffled threshold. Cumulative TE is a comparative summary across the common 72-step lag range; because physical horizons differ by native resolution, cross-resolution comparisons should retain that distinction.

## Reproducibility and verification

Configuration and dataset decisions are in `00_Data/Site_info.csv`. The complete input audit is in `04_Results/Data_Audit/Input_data_summary.csv`, and the verification manifest is in `04_Results/verification_manifest.csv`.

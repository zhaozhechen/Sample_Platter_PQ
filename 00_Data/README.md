# Data location

Raw Sampler Platter datasets remain in `../../Data` and are treated as read-only.
`Site_info.csv` records the exact file or timestamp-matched file pair, columns,
filters, transformations, and TE settings used for each analysis. The configured
maximum lag is 72 native observation steps; output lag labels are converted to
hours from the detected modal resolution.

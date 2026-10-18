# Reproducibility and evidence review

## Checked during portfolio preparation

- Read all top-level MATLAB sources and the complete 56-page report text.
- Inspected embedded Simulink functions and model constant/workspace references.
- Confirmed that all 13 paths loaded by the main script exist.
- Preserved included source/data bytes; SHA-256 hashes are recorded in the review manifest outside this repository.
- Visually reviewed selected rendered presentation pages.
- Transcribed the report's p. 26 circle-fit table to `results/reported-circle-fit.csv`; these are reported values, not rerun output.

## Runtime limitation

A batch run was attempted with the installed MATLAB R2025b. MATLAB terminated at startup with `System Error: File system inconsistency`, before executing the analysis. No regenerated numerical results or successful Simulink simulation are claimed.

## Issues to resolve before calling the analysis independently reproduced

1. Run `Experimental_project_latest` in MATLAB with the required toolboxes and check fitted frequencies, damping, coherence, and optimization convergence.
2. Reconcile clamp cancellation: the report describes subtraction, while the script adds 0.160 kg to the apparent mass. Confirm sensor polarity and the physical force convention.
3. Check the 10-400 Hz excitation range against the 10-500 Hz processing range.
4. Retain the full p. 26 NFQ ranges; the conclusion's blanket `NFQ > 0.999` statement does not match that table. The current circle-fit implementation uses an uncentred radial score labelled `R^2`; it should not be described as conventional regression R-squared.
5. The normalized FRF-vector correlation called MAC in the script is not a spatial mode-shape MAC.
6. Validate inverse-polynomial branches and bounds over measured data. The mapping is non-monotonic, and selection from setpoint frequency alone does not establish the actual thermal state.
7. Validate the model's thermal observer against temperature or hardware response. Observer-derived feedback is not independent physical feedback validation. The supplied model uses embedded functions and workspace quantities including inverse polynomial coefficients and `fn1_fit` / `fn2_fit` endpoints.
8. Report confidence intervals, repeated-test variation, calibration records, and the construction of merged low-current records if available. No supplied workflow generates those merged records.
9. Obtain the analytical-model implementation and the reference paper's complete citation before asserting that part is reproducible. Neither is supplied here.

The source is deliberately unchanged; these notes document issues rather than silently changing the experiment or controller.

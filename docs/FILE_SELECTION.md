# File selection review

## Included without modification

- `Experimental_project_latest.m`: primary analysis and inverse-polynomial fitting workflow.
- `autocross.m`: directly called by the primary analysis.
- `plotNyquistCircle.m`: directly called by the primary analysis.
- `circleFit.m`: called by `plotNyquistCircle`; this indirect dependency must be retained.
- `controller.slx`: supplied hybrid-control model, using embedded MATLAB functions.
- The 13 exact `Selected_data` files loaded by the primary script, with original relative paths.
- The original report, retaining its author list, and selected rendered pages.



## Scope of dependency review

Static review of all top-level `.m` files and the `.slx` XML package. The retained workflow is `Experimental_project_latest -> autocross` and `Experimental_project_latest -> plotNyquistCircle -> circleFit`. Controller coefficients are created in the primary script, and the model uses embedded functions. No file was excluded solely because the main script did not call it: indirect calls and embedded model references were also checked. Runtime dependency validation remains pending because MATLAB failed at startup.

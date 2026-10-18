# File selection review

## Included without modification

- `Experimental_project_latest.m`: primary analysis and inverse-polynomial fitting workflow.
- `autocross.m`: directly called by the primary analysis.
- `plotNyquistCircle.m`: directly called by the primary analysis.
- `circleFit.m`: called by `plotNyquistCircle`; this indirect dependency must be retained.
- `controller.slx`: supplied hybrid-control model, using embedded MATLAB functions.
- The 13 exact `Selected_data` files loaded by the primary script, with original relative paths.
- The original report, retaining its author list, and selected rendered pages.

## Excluded from this proposed repository

- `ATMD_Supervisory_Controller.m`: Not called by the main script or referenced by name in the supplied model XML.
- `buildLSCEmodel.m`: Not called by the main script or referenced by name in the supplied model XML.
- `extractPoles.m`: Not called by the main script or referenced by name in the supplied model XML.
- `LSCE_B.m`: Not called by the main script or referenced by name in the supplied model XML. Also calls `preprocessFRF`, which is not provided.
- `LSCE_from_FRF.m`: Not called by the main script or referenced by name in the supplied model XML. Contains an explicitly simplified pole approximation; not evidence of validated LSCE identification.
- `plotStab.m`: Not called by the main script or referenced by name in the supplied model XML.
- `plotStabilization.m`: Not called by the main script or referenced by name in the supplied model XML.
- `pole2modal.m`: Not called by the main script or referenced by name in the supplied model XML.
- `SMA_Controller_Nearest.m`: Not called by the main script or referenced by name in the supplied model XML. Contains example breakpoint values with comments asking for replacement; not treated as measured heating/cooling data.
- `svdReduction.m`: Not called by the main script or referenced by name in the supplied model XML.

Also excluded: the `.asv` backup, `.slxc` cache, `slprj/` generated products, MATLAB Drive tags, unused data records, duplicate/raw data trees, and two loose JPEGs superseded here by presentation figures. Nothing was removed from or changed in the original archive.

## Scope of dependency review

Static review of all top-level `.m` files and the `.slx` XML package. The retained workflow is `Experimental_project_latest -> autocross` and `Experimental_project_latest -> plotNyquistCircle -> circleFit`. Controller coefficients are created in the primary script, and the model uses embedded functions. No file was excluded solely because the main script did not call it: indirect calls and embedded model references were also checked. Runtime dependency validation remains pending because MATLAB failed at startup.

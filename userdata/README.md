# dev_varRBErobOpt Userdata Workspace

This folder is the versioned matRad userdata workspace used by the
`dev_varRBErobOpt` robust optimization workflows in this integration branch. It
keeps workflow macros, lightweight user resources, and selected patient data
inside the checkout that validates the workflows.

## Layout

- `macros`: executable workflow entrypoints organized by anatomical site
  (`prostate`, `breast`, `h&n`) plus shared local helpers.
- `patients`: patient geometry and DICOM/MAT input data location. The prostate
  case `3482` used by the workflow macros is tracked with Git LFS.
- `hluts`: Hounsfield-unit lookup tables.
- `machines`: user-configured machine data.
- `output`: generated workflow runs, cached dose influence matrices, figures,
  results, and performance records. This folder is intentionally not tracked and
  is created by workflow runs when needed.
- `scripts`: standard matRad user-script folder, kept available for standalone
  scripts.

## Workflows

Current workflow macros are:

- `macros/openWorkflowMacroBuilder.m`: GUI-only macro builder.
- `macros/prostate/runProstateInterval2Workflow.m`: prostate photon workflow
  using the `interval2_001` template and INTERVAL2 robust optimization.
- `macros/prostate/runProstateMultipleWorkflow.m`: prostate photon workflow
  using the `comparison_001` template with reference plus multiple robust plans.
- `macros/h&n/runHeadNeckInterval2Workflow.m`: H&N proton workflow using the
  `interval2_001` template and the `2F` beam set.

The macros use the `planWorkflow` submodule from the active matRad checkout. After
`matRad_rc` initializes matRad, the workflows construct `planWorkflow.Workflow`
and execute the staged workflow:

```matlab
workflow.prepare(...)
workflow.precompute(...)
workflow.pullDose(...)
workflow.optimize(...)
workflow.sample(...)
workflow.analyze(...)
workflow.save()
```

Run a macro from MATLAB after matRad has been initialized:

```matlab
cd(fullfile(userDataRoot,'macros'))
openWorkflowMacroBuilder

cd(fullfile(userDataRoot,'macros','prostate'))
runProstateInterval2Workflow
runProstateMultipleWorkflow

cd(fullfile(userDataRoot,'macros','h&n'))
runHeadNeckInterval2Workflow
```

## Patient Data

The prostate macros currently use `patients/prostate/3482`; the H&N macro
uses `patients/h&n/2`. DICOM cases are stored with Git LFS. After cloning or
updating this repository, run:

```bash
git lfs pull
```

## planWorkflow Path

In this integration checkout, `planWorkflow` is resolved from the matRad submodule:

```matlab
<matRadRoot>/submodules/planWorkflow
```

If the package is not available, initialize the submodule from the matRad root:

```bash
git submodule update --init submodules/planWorkflow
```

## Tests

Fast unit tests for `+planWorkflow` are included with the external toolbox. Run them
from MATLAB with:

```matlab
setupPlanWorkflow('/path/to/matRad');
planWorkflow.runTests('/path/to/matRad');
```

The tests cover strict configuration handling, analysis dose scaling,
robustness strategies, and workflow artifact persistence with synthetic data.

## Outputs

Workflow outputs are written below `output`, which is created by the workflow
when needed. Each run folder stores separated MAT files:

- `workflow_state.mat`: lightweight stage/state manifest.
- `workflow_data.mat`: heavy intermediate workflow data.
- `workflow_results.mat`: final optimization, sampling, and analysis results.
- `workflow_performance.mat`: timing and memory records per workflow stage.

Dose influence matrix cache files are stored below `output/cache`.

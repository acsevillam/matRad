# Workflow Macros

This folder contains executable workflow entrypoints for the external userdata
workspace, organized first by anatomy and then by radiation mode.

## Layout

- `breast/photons/`: photon breast workflow macros.
- `prostate/photons/`: photon prostate workflow macros.
- `prostate/protons/`: proton prostate workflow macros.
- `shared/prostate/`: shared prostate workflow implementations used by
  modality-specific wrappers.
- `helpers/`: shared local helper functions used by macros.

General macros:

- `openWorkflowMacroBuilder.m`: opens the planWorkflow GUI only, so templates
  and workflow macros can be exported without executing workflow stages.

Current breast photon macros:

- `breast/photons/runPhotonBreastInterval2Workflow.m`: breast workflow using the
  `interval2_001` template and INTERVAL2 robust optimization.

Current prostate photon macros:

- `prostate/photons/runPhotonProstate1MctRobustPtvWorkflow.m`: prostate workflow using the
  `PTV_001` template and robust PTV optimization for the `1_mct` case.
- `prostate/photons/runPhotonProstateInterval2Workflow.m`: prostate workflow using the
  `interval2_001` template and INTERVAL2 robust optimization.
- `prostate/photons/runPhotonProstateMultipleWorkflow.m`: prostate workflow using the
  `comparison_001` template with reference plus multiple robust plans.

Current prostate proton macros:

- `prostate/protons/runProtonProstateMultipleWorkflow.m`: prostate workflow using the
  `comparison_001` template with `protons`, `Generic`/`constRBE`, `RBExD`, and
  the `2F` beam set.

Shared prostate implementations:

- `shared/prostate/runProstateMultipleWorkflowCore.m`: shared prostate
  multiple-plan workflow. It receives `workflowConfig.prepare` from the caller,
  while the photon and proton wrappers provide their modality-specific prepare
  settings.

These macros use the active matRad checkout's `submodules/planWorkflow` package.
Each macro declares all effective settings in `workflowConfig` up front before
any `workflow.*` stage executes, using nested structs per stage such as
`workflowConfig.prepare`, `workflowConfig.precompute`,
`workflowConfig.pullDose`, and `workflowConfig.sampling`. The stage methods
are called without config structs, and `workflow.gui()` is called before
`workflow.prepare()` so the plan editor can modify the in-memory configuration
when MATLAB is running with UI support. Stage-level config arguments remain
supported by the toolbox API, but these macros do not use them for now.
`workflowConfig.prepare` owns the workflow anatomy (`description`), radiation
mode, anatomical template selection, and beam set. The plan target is defined by
the selected template objectives, not by the macro config.
Scenario settings are explicit per surface: the reference plan lives in
`workflowConfig.precompute.reference`, robust plans live in
`workflowConfig.precompute.robustPlans`, and sampling uncertainty settings live
under `workflowConfig.sampling`.

All macros accept optional overrides for `caseID`, `rootPath`, and
`cacheRootPath`. They can be passed as name-value pairs, a scalar struct, or
positional values in that order:

```matlab
cd(fullfile(userDataRoot,'macros'))
openWorkflowMacroBuilder()

cd(fullfile(userDataRoot,'macros','breast','photons'))
runPhotonBreastInterval2Workflow('caseID','4136','rootPath','/tmp/userdata')

cd(fullfile(userDataRoot,'macros','prostate','photons'))
runPhotonProstateInterval2Workflow('caseID','4136','rootPath','/tmp/userdata')
runPhotonProstateInterval2Workflow(struct('caseID','4136'))
runPhotonProstateInterval2Workflow('4136','/tmp/userdata','/tmp/cache')

cd(fullfile(userDataRoot,'macros','prostate','protons'))
runProtonProstateMultipleWorkflow('caseID','4136','rootPath','/tmp/userdata')
```

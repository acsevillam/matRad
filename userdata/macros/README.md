# Workflow Macros

This folder contains executable workflow entrypoints for the external userdata
workspace, organized by anatomical site.

## Layout

- `prostate/`: prostate workflow macros.
- `breast/`: breast workflow macros.
- `h&n/`: head-and-neck workflow macros.
- `helpers/`: shared local helper functions used by macros.

General macros:

- `openWorkflowMacroBuilder.m`: opens the planWorkflow GUI only, so templates
  and workflow macros can be exported without executing workflow stages.

Current prostate macros:

- `prostate/runProstateInterval2Workflow.m`: prostate workflow using the
  `interval2_001` template and INTERVAL2 robust optimization.
- `prostate/runProstateMultipleWorkflow.m`: prostate workflow using the
  `comparison_001` template with reference plus multiple robust plans.

Current head-and-neck macros:

- `h&n/runHeadNeckInterval2Workflow.m`: H&N workflow using the
  `interval2_001` template with `protons`, `Generic`/`constRBE`, and the
  `2F` beam set.

These macros use the active matRad checkout's `submodules/planWorkflow` package.
Each macro declares all effective settings in `workflowConfig` up front before
any `workflow.*` stage executes, using nested structs per stage such as
`workflowConfig.prepare`, `workflowConfig.precompute`,
`workflowConfig.dosePulling`, and `workflowConfig.sampling`. The stage methods
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

cd(fullfile(userDataRoot,'macros','prostate'))
runProstateInterval2Workflow('caseID','4136','rootPath','/tmp/userdata')
runProstateInterval2Workflow(struct('caseID','4136'))
runProstateInterval2Workflow('4136','/tmp/userdata','/tmp/cache')

cd(fullfile(userDataRoot,'macros','h&n'))
runHeadNeckInterval2Workflow('caseID','2')
```

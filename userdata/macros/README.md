# Workflow Macros

This folder contains executable workflow entrypoints for the external userdata
workspace, organized first by anatomy and then by radiation mode.

## Layout

- `breast/photons/`: photon breast workflow macros.
- `prostate/carbon/`: carbon prostate workflow macros.
- `prostate/helium/`: helium prostate workflow macros.
- `prostate/photons/`: photon prostate workflow macros.
- `prostate/protons/`: proton prostate workflow macros.
- `shared/breast/`: shared breast workflow implementations used by
  breast photon wrappers.
- `shared/prostate/`: shared prostate workflow implementations used by
  modality-specific wrappers.
- `helpers/`: shared local helper functions used by macros.

General macros:

- `openWorkflowMacroBuilder.m`: opens the planWorkflow GUI only, so templates
  and workflow macros can be exported without executing workflow stages.

Current breast photon macros:

- `breast/photons/runPhotonBreastMctPTVWorkflow.m`: breast workflow using the
  `PTV_001` template and PTV optimization for the `4136_mct` case.
- `breast/photons/runPhotonBreastInterval2Workflow.m`: breast workflow using the
  `interval2_001` template and INTERVAL2 robust optimization.
- `breast/photons/runPhotonBreastMultipleWorkflow.m`: breast workflow using the
  `comparison_001` template with reference plus PTV and INTERVAL2 robust plans.
- `breast/photons/runPhotonBreastMctMultipleWorkflow.m`: breast workflow using the
  `comparison_001` template with reference plus PTV and INTERVAL2 robust plans
  for the `4136_mct` case.

Current prostate photon macros:

- `prostate/photons/runPhotonProstateMctPTVWorkflow.m`: prostate workflow using the
  `PTV_001` template and PTV optimization for the `1_mct` case.
- `prostate/photons/runPhotonProstateInterval2Workflow.m`: prostate workflow using the
  `interval2_001` template and INTERVAL2 robust optimization.
- `prostate/photons/runPhotonProstateMultipleWorkflow.m`: prostate workflow using the
  `comparison_001` template with reference plus multiple robust plans.

Current prostate proton macros:

- `prostate/protons/runProtonProstateMctPTVWorkflow.m`: prostate workflow using the
  `PTV_001` template with `protons`, `Generic`/`constRBE`, `RBExD`, and the
  `2F` beam set for the `1_mct` case.
- `prostate/protons/runProtonProstateMultipleWorkflow.m`: prostate workflow using the
  `comparison_001` template with `protons`, `Generic`/`constRBE`, `RBExD`, and
  the `2F` beam set.

Current prostate carbon macros:

- `prostate/carbon/runCarbonProstateMctPTVWorkflow.m`: prostate workflow using the
  `PTV_001` template with `carbon`, `Generic`/`LEM`, `RBExD`, and the `2F`
  beam set for the `1_mct` case.

Current prostate helium macros:

- `prostate/helium/runHeliumProstateMctPTVWorkflow.m`: prostate workflow using the
  `PTV_001` template with `helium`, `Generic`/`HEL`, `RBExD`, and the `2F`
  beam set for the `1_mct` case.

Shared implementations:

- `shared/breast/runBreastMultipleWorkflowCore.m`: shared breast multiple-plan
  workflow. It receives `workflowConfig.prepare` from the caller, while the
  photon wrappers provide their case-specific prepare settings.
- `shared/prostate/runProstatePTVWorkflowCore.m`: shared prostate PTV workflow.
  It receives `workflowConfig.prepare` from the caller, while wrappers provide
  their case-specific prepare settings.
- `shared/prostate/prostateMctPrepareConfig.m`: common prostate `1_mct` case
  settings, including patient IDs, anatomy, template, and DICOM metadata.
- `shared/prostate/ensureProstatePrepareConfig.m`: shared prostate prepare
  config validator used by modality helpers.
- `shared/prostate/photonProstatePrepareConfig.m`: photon prostate modality
  settings.
- `shared/prostate/protonProstatePrepareConfig.m`: proton prostate modality
  settings.
- `shared/prostate/carbonProstatePrepareConfig.m`: carbon prostate modality
  settings.
- `shared/prostate/heliumProstatePrepareConfig.m`: helium prostate modality
  settings.
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
runPhotonBreastMctPTVWorkflow('rootPath','/tmp/userdata')
runPhotonBreastInterval2Workflow('caseID','4136','rootPath','/tmp/userdata')
runPhotonBreastMultipleWorkflow('caseID','4136','rootPath','/tmp/userdata')
runPhotonBreastMctMultipleWorkflow('rootPath','/tmp/userdata')

cd(fullfile(userDataRoot,'macros','prostate','photons'))
runPhotonProstateMctPTVWorkflow('rootPath','/tmp/userdata')
runPhotonProstateInterval2Workflow('caseID','4136','rootPath','/tmp/userdata')
runPhotonProstateInterval2Workflow(struct('caseID','4136'))
runPhotonProstateInterval2Workflow('4136','/tmp/userdata','/tmp/cache')

cd(fullfile(userDataRoot,'macros','prostate','protons'))
runProtonProstateMctPTVWorkflow('rootPath','/tmp/userdata')
runProtonProstateMultipleWorkflow('caseID','4136','rootPath','/tmp/userdata')

cd(fullfile(userDataRoot,'macros','prostate','carbon'))
runCarbonProstateMctPTVWorkflow('rootPath','/tmp/userdata')

cd(fullfile(userDataRoot,'macros','prostate','helium'))
runHeliumProstateMctPTVWorkflow('rootPath','/tmp/userdata')
```

# Workflow Macros

This folder contains local workflow entrypoints for planWorkflow.

Macros now use a single architecture:

- `breast/`, `prostate/`, `head_and_neck/`: workflow wrappers grouped only by
  anatomy and modality.
- `shared/specs/`: local `MacroSpec` catalog for anatomy, modality, case,
  template, robust plans, scenarios, and stage defaults.
- `helpers/`: thin helpers that load the spec catalog and the
  `planWorkflow` submodule.

Generic macro behavior lives in `planWorkflow.macros`:

- `planWorkflow.macros.MacroSpec.normalize` validates the public spec
  contract.
- `planWorkflow.macros.MacroRunner.build` builds `workflowConfig` without
  instantiating `planWorkflow.Workflow`; it is a support API for tests and the
  GUI builder, not a user macro mode.
- `planWorkflow.macros.MacroRunner.run` always executes the complete workflow
  stages.

## Spec IDs

Specs are resolved by `shared/specs/macroSpecCatalog.m` using:

```text
<site>.<modality>.<caseID>.<plan>
```

The execution profile is selected separately with `profile = 'prod'` or
`profile = 'testing'`. `prod` is the base behavior; `testing` applies only
explicit profile overrides.

Examples:

```matlab
macroSpecCatalog('breast.photons.4136_mct.COWC')
macroSpecCatalog('prostate.protons.1_mct.COWC_noAngles','profile','testing')
macroSpecCatalog('head_and_neck.photons.2.INTERVAL2','profile','prod')
```

Supported sites are:

- `breast`
- `prostate`
- `head_and_neck`, mapped internally to planWorkflow description `h&n`

Plan IDs stay canonical across anatomies when the template supports them:

```text
PTV, COWC, STOCH, cCOWC, PROB2, INTERVAL2, INTERVAL3
```

## Defaults

`prod` profile uses:

- `executionMode = 'run'`
- `openGui = true`
- `lockPlanSet = true`
- `allowCustomRobustPlans = false`
- `precompute.doseResolution = [3 3 3]`

`testing` profile inherits the same behavior and overrides only:

- `precompute.doseResolution = [5 5 5]`

Each wrapper declares its GUI default explicitly:

```matlab
% Execution defaults.
% profile: prod uses production defaults; testing applies explicit overrides.
% openGui: true opens the planWorkflow GUI before running stages.
profile = 'prod';
openGui = true;

% MacroSpec selectors.
% site: anatomy/catalog key (breast, prostate, head_and_neck).
% particleType: radiation mode used by prepare.radiationMode.
% caseID: patient/case identifier resolved under userdata.
% robustness: canonical robust plan key or multiple-plan selector.
% samplingProfile: sampling profile (default, noAngles).
site = 'breast';
particleType = 'photons';
caseID = '4136_mct';
robustness = 'cCOWC';
samplingProfile = 'default';

% Optimization scenario toggles.
% Controls robust scenarios used for precompute and optimization plans.
optimizationScenario = struct( ...
    'ctActive', false, ...
    'setupActive', true, ...
    'rangeActive', false, ...
    'gantryActive', false, ...
    'couchActive', false);

% Sampling scenario toggles.
% Controls active uncertainty dimensions used by the sampling stage.
samplingScenario = struct( ...
    'ctActive', false, ...
    'setupActive', true, ...
    'rangeActive', false, ...
    'gantryActive', false, ...
    'couchActive', false);

% Resolve the MacroSpec and run the complete workflow.
specId = macroSpecId( ...
    site,particleType,caseID,robustness,samplingProfile);
```

Change `profile` to `testing` in a specific macro if its default should be a
local testing run. Change `openGui` to `false` when the default run should be
non-interactive, for example on a cluster. Change `particleType`, `caseID`, or
`robustness` to point the macro to a different catalog spec, such as
`photons`, `4136_mct`, and `cCOWC`. Use `samplingProfile = 'noAngles'` for
the no-angle prostate profiles; this keeps `robustness` as the robust plan
selector instead of encoding sampling behavior in that field. Call-time
`profile` and `openGui` arguments still override the wrapper defaults.

`optimizationScenario` controls the robust scenario used to build the
optimization/precompute robust plans. `samplingScenario` controls the scenario
dimensions used by the sampling stage. Both can also be overridden at call time.

The same options are accepted by all wrappers:

```text
profile, caseID, rootPath, cacheRootPath, planKeys, planTemplate,
randomSeed, openGui, optimizationScenario, samplingScenario, lockPlanSet,
overrides, allowCustomRobustPlans
```

Stage overrides may be passed as `prepare`, `precompute`, `pullDose`,
`optimize`, `sampling`, or `analysis` structs.

## Serial Jobs

Serial jobs run macro parameter sets one after another. They are useful for
local testing or cluster submission when the same macro needs to be repeated
with changed parameters, such as robustness, patient case, profile, GUI mode,
or scenario toggles.

Jobs support two equivalent shapes:

- `runs`: explicit fully resolved run entries.
- `base` or `bases` plus `parameterSets`: each parameter set is merged over
  each base, then executed in series.

Example with one base and several robustness parameter sets:

```matlab
job = struct();
job.id = 'breast.photons.4136_mct.robustnessSweep';
job.profile = 'testing';
job.openGui = false;
job.site = 'breast';
job.particleType = 'photons';
job.caseID = '4136_mct';
job.samplingProfile = 'default';
job.optimizationScenario = optimizationScenario;
job.samplingScenario = samplingScenario;
job.parameterSets = repmat(struct('label','','robustness',''),1,3);
job.parameterSets(1).label = 'PTV';
job.parameterSets(1).robustness = 'PTV';
job.parameterSets(2).label = 'COWC';
job.parameterSets(2).robustness = 'COWC';
job.parameterSets(3).label = 'cCOWC';
job.parameterSets(3).robustness = 'cCOWC';

jobResult = runWorkflowMacroJob(job,'rootPath',userDataRoot);
```

Example with several bases and shared parameter sets:

```matlab
job.bases = repmat(struct('label','','site','','particleType','', ...
    'caseID','','samplingProfile','default'),1,2);
job.bases(1).label = 'breast';
job.bases(1).site = 'breast';
job.bases(1).particleType = 'photons';
job.bases(1).caseID = '4136_mct';
job.bases(2).label = 'hn';
job.bases(2).site = 'head_and_neck';
job.bases(2).particleType = 'photons';
job.bases(2).caseID = '2';
job.parameterSets = repmat(struct('robustness',''),1,2);
job.parameterSets(1).robustness = 'PTV';
job.parameterSets(2).robustness = 'INTERVAL2';
```

Example with the same plan and different patients:

```matlab
job = struct();
job.id = 'breast.photons.multiple.caseSweep';
job.profile = 'testing';
job.openGui = false;
job.site = 'breast';
job.particleType = 'photons';
job.robustness = 'multiple';
job.samplingProfile = 'default';
job.parameterSets = repmat(struct('caseID',''),1,2);
job.parameterSets(1).caseID = '4136_mct';
job.parameterSets(2).caseID = '4136';
```

Use `resolveWorkflowMacroJob(job)` to inspect the resolved run list without
executing stages. Use `runWorkflowMacroJob(job)` to execute all resolved runs
serially. A ready-to-edit example is available under `jobs/`.

## Examples

```matlab
cd(fullfile(userDataRoot,'macros','breast','photons'))
runBreastPhotonMctCOWCWorkflow('rootPath',userDataRoot)
runBreastPhotonMctCOWCWorkflow('profile','testing', ...
    'rootPath',userDataRoot)
runBreastPhotonMctCOWCWorkflow('openGui',false, ...
    'rootPath',userDataRoot)
runBreastPhotonMctCOWCWorkflow( ...
    'optimizationScenario',struct('ctActive',false), ...
    'samplingScenario',struct('gantryActive',false,'couchActive',false), ...
    'rootPath',userDataRoot)

cd(fullfile(userDataRoot,'macros','prostate','protons'))
result = runProstateProtonMctCOWCWorkflow( ...
    'rootPath',userDataRoot, ...
    'sampling',struct('sampling_shiftSD',[3 6 3]));

cd(fullfile(userDataRoot,'macros','head_and_neck','photons'))
runHeadNeckPhoton2Interval2Workflow('rootPath',userDataRoot)
```

The macro builder can open any spec in the GUI without running stages:

```matlab
cd(fullfile(userDataRoot,'macros'))
openWorkflowMacroBuilder('prostate.photons.3482.INTERVAL2', ...
    'profile','testing')
openWorkflowMacroBuilder('specId','breast.photons.4136_mct.COWC')
```

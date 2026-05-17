# robOpt planWorkflow Userdata Workspace

This folder contains the local userdata used by the `robOpt_planWorkflow`
integration checkout.

## Layout

- `macros`: workflow entrypoints, jobs, helpers, and local macro specs for
  `planWorkflow`.
- `hluts`: local Hounsfield-unit lookup tables.
- `machines`: local machine data.
- `patients`: local clinical cases required by the integration workflows. This
  folder is ignored.
- `output`: generated workflow outputs and caches. This folder is ignored.
- `scripts`: standalone user scripts.

## Usage

Initialize matRad from the integration checkout, then run macros from the
`userdata/macros` tree:

```matlab
cd('/path/to/robOpt_planWorkflow')
matRad_rc
cd(fullfile(userDataRoot,'macros'))
openWorkflowMacroBuilder
```

The `planWorkflow` package is loaded from:

```text
submodules/planWorkflow
```

## Patient Data

The integration clinical cases are expected locally under `userdata/patients`.
They are intentionally ignored and are not tracked as a submodule or as Git
files. Populate this folder from the local clinical cases checkout before
running workflows.

Generated outputs and local patient data stay local to the integration
environment and should not be committed.

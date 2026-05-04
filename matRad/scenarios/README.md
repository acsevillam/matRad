Scenario module layout
======================

This folder exposes the public scenario creation entry points:

- `matRad_createScenarioModel.m`: preferred factory for scenario models.
- `matRad_multScen.m`: deprecated compatibility wrapper.

Implementation files are grouped by responsibility:

- `base/`: abstract base classes and shared scenario interfaces.
- `models/`: concrete scenario model implementations.
- `applicators/`: applicators that apply uncertainty realizations to CT,
  setup, range, gantry, and couch dimensions.
- `helpers/`: shared helper functions used by scenario model construction.

Naming conventions
------------------

- `scenarioId`: stable public id of a scenario realization.
- `scenarioRowIx`: row position in the active scenario realization table.
- `fullScenIx` or `scenarioDijIx`: linear index into DIJ/scenario-mask cell
  arrays.
- `<dimension>ScenId`: stable id when a dimension references an external
  discrete object. Example: `ctScenId` indexes CT/cube/DVF data.
- `<dimension>ScenIx`: local position or subscript in a dimension grid/list,
  only when the code explicitly resolves that position to a dimension value
  or id. Examples: `ctScenIx`, `setupScenIx`, `rangeScenIx`.
- `<dimension>Value(s)`: realization values consumed by applicators. Existing
  names such as `setupShift` and `rangeShift` are compatibility helpers for
  current applicators; generic scenario code should use `scenarioValues` and
  realization component names.

Dimension activity
------------------

- `scenarioDimensionActive` is the public list of active uncertainty
  dimensions/applicators:
  `ct`, `setup`, `range`, `gantry`, and `couch`.
- Internal component names such as `setup.x` and `range.absolute` are
  realization components, not scenario dimensions.
- Angular component names such as `gantry.beam1` and `couch.beam1` are
  per-beam realization components generated from the public `gantry` and
  `couch` dimensions. Angular values are stored in degrees.
- Active continuous components must have strictly positive uncertainty scales.
- Inactive continuous components remain present in `scenarioComponents` and
  `scenarioValues` for compatibility, but are generated at nominal value and
  ignored by probability evaluation.
- CT scenarios are controlled by `ctScenProb`: include only CT scenario ids
  that should be active. Included rows must have positive probability.

Model support
-------------

- `nomScen` generates nominal CT/setup/range/gantry/couch realizations.
- `rndScen` supports random setup, range, gantry, and couch uncertainty. Angular
  dimensions require `numOfBeams` to be set before scenarios are built; each
  sample stores one gantry and/or couch offset per beam.
- Grid-based models (`wcScen`, `impScen`, and `truncatedImpScen`) currently
  support setup and range dimensions only. They intentionally reject active
  `gantry` or `couch` dimensions until a gridded angular scenario convention is
  implemented.
- Dose engines that do not support angular uncertainty must continue to reject
  active `gantry` or `couch` dimensions explicitly.

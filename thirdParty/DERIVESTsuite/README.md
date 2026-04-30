DERIVESTsuite
=============

This folder contains the minimal subset of John D'Errico's DERIVESTsuite
required by matRad's c-COWC Cheap-Minimax robust optimization model:

- `gradest.m`
- `derivest.m`

Source
------

John D'Errico, Adaptive Robust Numerical Differentiation, MATLAB Central
File Exchange, file ID 13490.

https://www.mathworks.com/matlabcentral/fileexchange/13490-adaptive-robust-numerical-differentiation

Usage in matRad
---------------

`gradest` is used to estimate the scenario-level gradient of the
Cheap-Minimax rank aggregation. The radiotherapy model implemented in matRad
is based on:

Sevilla AC, Cabal G, Wahl N, Puerta ME, Rivera JC. A robust optimization model
for intensity-modulated radiotherapy: Cheap-Minimax. Medical Physics.
2025;52:3360-3376. doi:10.1002/mp.17709.

Maintenance scope
-----------------

This is pure MATLAB code with no MEX build step, no network access, and no
external runtime beyond the two files vendored here. matRad only calls it when
an optimization objective explicitly uses the `c-COWC` robustness setting.

License
-------

DERIVESTsuite is distributed under the BSD license used by MATLAB Central File
Exchange submissions. See `LICENSE.md` in this folder.

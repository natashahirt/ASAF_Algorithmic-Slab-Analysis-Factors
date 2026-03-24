# AsapToolkit (Vendored)

This directory contains a vendored copy of the `AsapToolkit`. 

## Purpose
`AsapToolkit` provides geometric generators, force analysis utilities, and IO functions for interacting with structural models. It is included locally to ensure compatibility and allow for specific structural synthesis modifications.

## Structure
- `AsapSections/`: Section property calculations.
- `ForceAnalysis/`: Internal force and envelope calculations.
- `Generation/`: Parametric structural generators (Trusses, Frames, Domes).
- `Geometry/`: Displacement and model geometry utilities.
- `SteelSections/`: Database and types for standard steel shapes.

## Maintenance
Changes to this toolkit should be made carefully, as it serves as the geometric backbone for `StructuralSynthesizer`.

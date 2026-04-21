# vsclib

Internal ViennaCS lookup database files used by the Python workflow.

## Layout

- `implant/*.json`: recipe defaults and table bindings for implant/damage models.
- `anneal/*.json`: dopant/defect/cluster Arrhenius defaults for anneal model setup.

## Design intent

- Runtime should consume this normalized database format, not raw external source files.
- Database entries keep source provenance fields so values can be traced to original calibration files.
- Config keys can override any value at runtime.

## Current files

- `implant/si_taurus_defaults.json`
  - Maps species/material/substrate to internal implant and damage table copies in `implant/tables/`.
- `anneal/si_score_defaults.json`
  - Provides Si/Ge anneal defaults derived from `data/lib/score/Params/*`.

## Example integration keys (`examples/implantMask2D/implantMask2D.py`)

- Implant DB:
  - `implantUseInternalDb=1`
  - `implantDbFile=../../vsclib/implant/si_taurus_defaults.json`
- Anneal DB:
  - `annealUseInternalDb=1`
  - `annealDbFile=../../vsclib/anneal/si_score_defaults.json`

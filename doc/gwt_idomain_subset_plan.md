# Transport on a subset of the flow model grid

## Purpose

Allow the IDOMAIN array of a GWT model to define a transport domain that is
smaller than the active domain of the GWF model that supplies its flows. Before
this change, a GWF-GWT exchange required the two models to have identical
IDOMAIN arrays, so a transport simulation always had to carry every active flow
cell.

The GWT IDOMAIN must be a *subset* of the GWF IDOMAIN: a cell that is active in
the transport model must be active in the flow model, but the transport model
may exclude cells that are active in the flow model.

## Implementation

### Grid maps

Both models read the same user grid, so the two reduced node numberings differ
only by the excluded cells. The FMI package (`FlowModelInterfaceType`) stores
the maps between them and is put into a *mapped* mode by `igwfmapped`:

| array        | length          | contents                                        |
| ------------ | --------------- | ----------------------------------------------- |
| `gwfnodemap` | transport nodes | flow model node for each transport node         |
| `gwfnodeinv` | flow nodes      | transport node for each flow node, 0 if excluded |
| `gwfjamap`   | transport nja   | flow connection for each transport connection    |
| `gwfdropia`  | transport nodes + 1 | index into `gwfdropja` for each transport node |
| `gwfdropja`  | dropped conns   | flow connections that are excluded from transport |

`map_gwf_grid` builds the maps in `exg_df`, before the transport model arrays
are allocated, using `get_nodeuser` and `get_nodenumber` so that DIS, DISV, and
DISU are all handled without a `select type`. A transport connection with no
counterpart in the flow model (for example, a vertical pass-through cell that
exists only in the transport model) is an error.

### Array transfer

In mapped mode the FMI arrays `gwfhead`, `gwfsat`, `gwfspdis`, `gwfflowja`,
`gwfstrgss`, and `gwfstrgsy` are owned by FMI and sized on the transport grid,
exactly as they are when flows are read from a budget file. The exchange stores
pointers to the corresponding flow model arrays through `set_gwf_sources` and
`set_gwf_storage`, and `map_gwf_values` refills the transport arrays in
`fmi_ad`. The transfer is done at the start of the transport model advance,
after the flow solution for the time step is complete, so downstream packages
(ADV, DSP, MST, SSM, ...) index the arrays with transport node numbers and need
no changes.

Boundary package flows reach the transport model through `PackageBudgetType`.
In mapped mode `nodelist` stays owned by the package budget object, `nbound`,
`flow`, and `auxvar` still point into the flow model package, and
`map_nodelist` refills `nodelist` from the flow model nodelist each time step.
Entries in excluded cells are set to zero, which SSM already skips.

### Flow across the transport boundary

`gwfflowja(idiag)` carries the flow residual for a cell. Water that moves across
a connection the transport model does not have is subtracted from the diagonal
in `map_gwf_values`, so the residual seen by the transport model is the residual
of the reduced domain rather than the residual of the flow model.

That residual is handled by the FMI flow imbalance correction: outflow leaves at
the concentration of the cell and inflow enters at zero concentration. For a
purely advective problem with upstream weighting this reproduces the full-domain
solution exactly in the retained cells. FLOW\_IMBALANCE\_CORRECTION is therefore
effectively required, and the exchange writes a warning when a reduced transport
domain is used without it. A specified concentration on the perimeter of the
transport domain (CNC, or an SSM boundary) is the way to bring mass in.

## Unsupported combinations

Each of these is reported as an error rather than silently producing a wrong
answer:

- a transport cell that is inactive in the flow model;
- a transport connection that does not exist in the flow model;
- an advanced transport package (LKT, SFT, MWT, UZT) connected to a cell outside
  the transport domain, checked in `exg_df` from the flow package budget object
  because the advanced package equations are assembled from flow model cell
  numbers before the exchange is allocated and read;
- a reduced transport domain in a model that is also coupled through a GWT-GWT
  exchange;
- a reduced transport domain in a parallel simulation.

## Not done yet

- GWE. The FMI machinery is shared, so enabling GWE is the same change in
  `exg-gwfgwe.f90` plus the GWE discretization definition files.
- Advanced packages that straddle the transport boundary.
- GWT-GWT exchanges and parallel runs, which need the interface model and the
  virtual model maps to carry the same node and connection maps.

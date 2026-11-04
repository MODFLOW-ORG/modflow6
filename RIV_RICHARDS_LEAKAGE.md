# Uncapped river leakage for Richards flow

## Physical interpretation

A Richards cell can have negative water pressure without losing its ability
to receive water. Water moves according to total hydraulic head:

$$
h = z + \psi,
$$

where $z$ is elevation and $\psi$ is pressure head, which can be negative in
unsaturated soil.

For a wetted riverbed, the river provides head $H_r = STAGE$. Representing the
bed as a hydraulic resistance gives:

$$
Q = C(H_r - h),
$$

where $Q$ is positive into groundwater and $C$ is the effective
riverbed/interface conductance. Negative soil pressure increases the head
difference and draws water from the river. This represents capillary suction,
not an invalid head or an inactive cell.

## Why conventional RIV caps leakage

The conventional RIV relationship is:

$$
Q = C[H_r - \max(h, RBOT)].
$$

It prevents the downstream head seen by the bed from falling below `RBOT`.
This corresponds to assuming atmospheric pressure at the bed underside once
groundwater head falls below that elevation.

Richards flow does not require that assumption: water below the bed may remain
under suction rather than at atmospheric pressure. Removing the cap allows
that suction to influence exchange.

## Assumptions and limitations

Being a Richards cell does not, by itself, make constant-conductance uncapped
leakage physically valid. The approximation additionally assumes:

- A wetted river supplies water continuously.
- The bed maintains liquid continuity and approximately constant conductance.
- The cell head adequately represents the downstream interface head, or
  `COND` includes the unresolved resistance between them.

If the bed desaturates, its conductivity decreases. If an atmospheric-pressure
separation develops, suction is no longer transmitted through the interface
as assumed. Uncapped constant-conductance leakage can then overestimate
infiltration.

Richards conductivity inside the aquifer does not automatically correct a
missing or misrepresented interface resistance.

The defensible description is therefore:

> A connected, riverbed-resistance boundary that permits suction-driven
> infiltration into Richards soil.

It is not a universally correct river boundary for Richards flow.

## Minimal, formulation-independent design

Expose the choice as a boundary leakage law, rather than making RIV detect
whether VSF is active:

| Leakage law | Relationship | Interpretation |
|---|---|---|
| Capped (existing default) | $Q = C[H_r - \max(h, RBOT)]$ | Conventional disconnected-river approximation |
| Uncapped (opt-in) | $Q = C(H_r - h)$ | Connected, resistive riverbed with suction-driven infiltration |

Both laws use total hydraulic head. Neither requires RIV to inspect VSF soil
models, relative permeability, or cell formulation flags. The uncapped law
can also be meaningful for saturated flow.

For a first implementation, preserve the existing default and introduce an
explicit package-level leakage-law option. In uncapped mode, retain `RBOT`
for input compatibility and applicable riverbed geometry checks, but document
that it does not cap leakage.

Do not automatically multiply riverbed conductance by the receiving cell's
relative permeability: that would conflate riverbed resistance with aquifer
resistance. Explicit treatment of bed desaturation or aquifer half-cell
resistance would require a more detailed interface model.

This is a conceptual proposal only; no implementation is included.

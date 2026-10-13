# Semi-Lagrangian transport proofs

This project contains Lean proofs for the semi-Lagrangian transport method.

Lean and Mathlib use version 4.34.0.

## Build

1. Enter `formal/SlTransportProofs`.
2. Run `lake update` after a dependency change.
3. Run `lake exe cache get` after a Mathlib change.
4. Run `lake build` to check all accepted proofs.

The accepted proof files must not contain `sorry`, `admit`, or custom axioms.

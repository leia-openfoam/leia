/-
Copyright (c) 2026 Tomislav Marić. All rights reserved.
Released under GNU General Public License version 3 or later.
Authors: Tomislav Marić
-/
import Mathlib

set_option linter.style.header false

/-!
# Configuration smoke test

This module checks the Lean and Mathlib project configuration.
-/

namespace SlTransportProofs

theorem smoke_test : (1 : ℕ) + 1 = 2 := by
  norm_num

end SlTransportProofs

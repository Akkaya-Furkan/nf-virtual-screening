Smina custom scoring files.

- `vinardo_like.txt`: the four Vinardo terms and weights printed by `smina --scoring vinardo`. It is NOT
  identical to the built-in `vinardo` (acetazolamide crystal-frame pose: -5.63 vs -4.88 kcal/mol; atom typing
  or intramolecular handling differs), so it is only compared with its own zinc variants.
- `vinardo_like_zn1.txt`, `_zn2.txt`: add a Gaussian reward (weight -1 or -2) for a donor nitrogen
  (`NitrogenXSDonor`, i.e. a sulfonamide NH2) at 2.0 Å (width 0.4) from the zinc. This is a pharmacophore-informed
  restraint-like term, not a physical metal-coordination model; weights are arbitrary and were not fitted.

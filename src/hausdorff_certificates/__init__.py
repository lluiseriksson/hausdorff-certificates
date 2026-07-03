"""hausdorff-certificates: finite PSD certificates H_N, L_N for Hausdorff
moment sequences, with an exact rational backend and an interval backend.

Implementation satellite of the physmath-knowledge-tree node
``bridge.hausdorff_moment_certificates``.  Scope contract: certifies finite
matrices only; PSD outcomes are consistency checks, certified negative
eigenvalues are finite falsifiers.  Nothing here proves RH, coercivity in
infinite volume, hRpoly, or anything Clay-adjacent.
"""

__version__ = "0.1.0"

from .certify import (  # noqa: F401
    Certificate,
    VERDICT_INCONCLUSIVE,
    VERDICT_NOT_PSD,
    VERDICT_PSD,
    certify_exact,
    certify_interval,
    coercivity_bundle_exact,
    hausdorff_pair_exact,
    hausdorff_pair_interval,
)
from .hankel import (  # noqa: F401
    complete_monotonicity_table,
    cm_violations,
    hankel_H,
    hankel_L,
    hankel_L_theta,
    shifted_hankel_S,
)
from .intervals import Interval  # noqa: F401
from .moments import (  # noqa: F401
    lambda_min_lower_bound_from_theta,
    moments_from_atoms,
    moments_lebesgue,
    resolvent_trace_moments,
    theta_for_coercivity,
)
from .rational import exact_psd, quadratic_form  # noqa: F401

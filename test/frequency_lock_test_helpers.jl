# Shared test helper. The including module must have `Tracking`, `Hz` and `s`
# in scope.

# A frequency-locked indicator.
_locked() = Tracking.FrequencyLockIndicator(0.0Hz * 0.0s, 0.0s, true)

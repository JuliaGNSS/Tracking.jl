struct NoSatPostProcess <: AbstractSatPostProcess end
struct NoSystemPostProcess <: AbstractSystemPostProcess end
struct NoTrackingPostProcess <: AbstractTrackingPostProcess end

post_process(
    track_state::TrackState{
        <:MultipleSystemSatsState,
        <:AbstractTrackingDopplerEstimator,
        <:AbstractDownconvertAndCorrelator,
        <:NoTrackingPostProcess,
    },
) = track_state
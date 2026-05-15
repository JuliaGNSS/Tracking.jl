module TrackingAcquisitionExt

using Tracking: Tracking, TrackedSat
using Acquisition: AcquisitionResults

function Tracking.TrackedSat(acq::AcquisitionResults; args...)
    TrackedSat(acq.system, acq.prn, acq.code_phase, acq.carrier_doppler; args...)
end

end

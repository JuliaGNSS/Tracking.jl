module TrackingAcquisitionExt

using Tracking: Tracking, SatState
using Acquisition: AcquisitionResults

function Tracking.SatState(acq::AcquisitionResults; args...)
    SatState(acq.system, acq.prn, acq.code_phase, acq.carrier_doppler; args...)
end

end

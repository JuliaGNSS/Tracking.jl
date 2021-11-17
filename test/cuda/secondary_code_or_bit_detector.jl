@testset "CUDA: Secondary code or bit detector" begin
    detector = SecondaryCodeOrBitDetector()
    gpsl1 = GPSL1(use_gpu = Val(true))

    @test Tracking.get_buffer(detector) == 0
    @test Tracking.length(detector) == 0
    @test Tracking.found(detector) == false

    next_detector = Tracking.find(gpsl1, detector, -1.0 + 0.0im)
    @test Tracking.get_buffer(next_detector) == 0
    @test Tracking.length(next_detector) == 1
    @test Tracking.found(next_detector) == false

    next2_detector = Tracking.find(gpsl1, next_detector, 1.0 + 0.0im)
    @test Tracking.get_buffer(next2_detector) == 1
    @test Tracking.length(next2_detector) == 2
    @test Tracking.found(next2_detector) == false

    next3_detector = Tracking.find(gpsl1, next2_detector, 1.0 + 0.0im)
    @test Tracking.get_buffer(next3_detector) == 3
    @test Tracking.length(next3_detector) == 3
    @test Tracking.found(next3_detector) == false

    detector = SecondaryCodeOrBitDetector(1, 1, true)

    @test Tracking.get_buffer(detector) == 1
    @test Tracking.length(detector) == 1
    @test Tracking.found(detector) == true

    next_detector = Tracking.find(gpsl1, detector, 1.0 + 0.0im)
    @test Tracking.get_buffer(next_detector) == 1
    @test Tracking.length(next_detector) == 1
    @test Tracking.found(next_detector) == true

    detector = SecondaryCodeOrBitDetector(0x7ffff80000, 39, false)
    next_detector = Tracking.find(gpsl1, detector, -1.0 + 0.0im)
    @test Tracking.get_buffer(next_detector) == 0x000000fffff00000
    @test Tracking.length(next_detector) == 40
    @test Tracking.found(next_detector) == true
end

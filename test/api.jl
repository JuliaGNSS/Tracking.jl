module APITest

using Test: @test, @testset
import Tracking

# The API is what the manual's `@docs` blocks document: their bare names are the exports,
# their `Tracking.`-qualified ones the unexported names a caller writes qualified, and
# their `TrackingLoops.`- and `Tracking.PositionVelocityTime.`-qualified ones the methods
# this package adds to those packages' functions. The TrackingLoops readers Tracking
# re-exports are documented in TrackingLoops' manual.
@testset "the manual documents the exports bare and the rest qualified" begin
    bare, qualified = Symbol[], Symbol[]
    for page in readdir(joinpath(@__DIR__, "..", "docs", "src"); join = true)
        endswith(page, ".md") || continue
        in_docs = false
        for line in eachline(page)
            startswith(line, "```") && (in_docs = startswith(line, "```@docs"); continue)
            (in_docs && !isempty(strip(line))) || continue
            name = first(split(strip(line), '('))
            startswith(name, r"TrackingLoops\.|Tracking\.PositionVelocityTime\.") &&
                continue
            if startswith(name, "Tracking.")
                push!(qualified, Symbol(chopprefix(name, "Tracking.")))
            else
                push!(bare, Symbol(name))
            end
        end
    end
    # Without `all`, `names` lists the exports; the re-exports are TrackingLoops'.
    exported = filter(names(Tracking)) do name
        name !== :Tracking && Base.binding_module(Tracking, name) === Tracking
    end
    @test sort(unique(bare)) == sort(exported)
    @test all(name -> isdefined(Tracking, name) && name ∉ names(Tracking), qualified)
    @test all(name -> Base.Docs.hasdoc(Tracking, name), [exported; qualified])
end

end

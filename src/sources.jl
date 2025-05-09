module sources

using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using CoolProp.CoolProp
using ..utils: FluidPort, SISO

export Source, Sink, Massflow

## End-points
@mtkmodel Source begin
    @parameters begin
        h
        p
    end

    @components begin
        outlet = FluidPort()
    end

    @equations begin
        outlet.h ~ h
        outlet.p ~ p
    end
end

@mtkmodel Sink begin
    @components begin
        inlet = FluidPort()
    end
end

## In the middle
@mtkmodel Massflow begin
    @extend inlet, outlet = siso = SISO()

    @parameters begin
        ṁ
    end

    @equations begin
        outlet.ṁ ~ ṁ
    end
end
end
 
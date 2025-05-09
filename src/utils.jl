module utils

using ModelingToolkit, OrdinaryDiffEq
using ModelingToolkit: t_nounits as t, D_nounits as D
using CoolProp.CoolProp
using CoolProp_jll

export calc_fluid_state, calc_fluid_property, dp_dT_calc, FluidPort, FluidProperties

## Helper functions
function calc_fluid_state(in1, in2, input_pair::String, fluid_model)
    inputs = get_input_pair_index(input_pair)
    AbstractState_update(fluid_model, inputs, in1, in2)
    return fluid_model
end

function calc_fluid_property(property::String, fluid_state)
    property_index = get_param_index(property)
    return AbstractState_keyed_output(fluid_state, property_index)
end

"""
function dp_dT_calc(fluid_state, fluid_model)
    p1 = calc_fluid_property("P", fluid_state)
    T1 = calc_fluid_property("T", fluid_state)
    ρ = calc_fluid_property("D", fluid_state)

    T2 = T1 + 1e-2
    fluid_state_2 = calc_fluid_state(ρ, T2, "DmassT_INPUTS", fluid_model)
    p2 = calc_fluid_property("P", fluid_state_2)

    dp_dT = (p2 - p1)/(T2 - T1)
    return dp_dT
end
"""

errcode = Ref{Clong}(0)
const message_buffer = Array{UInt8}(undef, 262144)

function AbstractState_first_partial_deriv_mine(handle::Clong, of::Clong, wrt::Clong, constant::Clong)
    buffer_length = length(message_buffer)
    output = ccall( (:AbstractState_first_partial_deriv, libcoolprop), Cdouble, (Clong, Clong, Clong, Clong, Ref{Clong}, Ptr{UInt8}, Clong), handle, of, wrt, constant, errcode, message_buffer::Array{UInt8, 1}, buffer_length)
    if output == -Inf
        error("CoolProp: no correct state has been set with AbstractState_update")
    end
    return output
end

function dp_dT_calc(fluid_state)
    return AbstractState_first_partial_deriv_mine(
        fluid_state, get_param_index("P"), get_param_index("T"), get_param_index("D"))
end

@register_symbolic calc_fluid_state(in1::Real, in2::Real, input_pair::String, fluid_model::Int32)
@register_symbolic calc_fluid_property(property::String, fluid_state::Int32)
@register_symbolic dp_dT_calc(fluid_state::Int32)

## Connectors
@connector FluidPort begin
    @parameters begin
        fluid_model::Int32
    end

    @variables begin
        h(t), [guess=0.0]
        p(t), [guess=0.0]
        ṁ(t), [guess=0.0, connect = Flow]
    end
end

@connector FluidProperties begin
    @parameters begin
        fluid_model::Int32
    end

    @variables begin
        ṁ(t), [guess=0.0, connect = Flow]
    end

    @equations begin
        ṁ ~ 0
    end
end

# Base class
@mtkmodel SISO begin
    @components begin
        inlet = FluidPort()
        outlet = FluidPort()
    end

    @equations begin
        inlet.ṁ + outlet.ṁ ~ 0
        domain_connect(inlet, outlet)
    end
end

end

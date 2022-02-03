using FFTW
using Scapin
using Scapin.Elasticity
using Scapin.Bri17

const T = Float64
const d = 2

C = Hooke{d,T}(1.0, 0.3)
N = (16, 8)
h = (1.0, 1.0)
Γ = DiscreteGreenOperatorBri17{d,T}(C, N, h)

τ = rand(div(d * (d+1), 2), N...)

ε_exp = apply(Γ, τ)

plan = plan_rfft(τ, 2:(d+1))
ε_act = apply(Γ, τ, plan)

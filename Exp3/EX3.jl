using DSP
using Plots
using FixedPointNumbers
using SampledSignals
using Statistics
using PyCall
using SciPy
using Polynomials
plotlyjs()

fa = 40000;
f0 = 100π;
t = 0:1/fa:3;
s = 0.5.*cos.(2π*f0.*t) + 0.3.*cos.(4π*300 .*t) + 0.15.*cos.(6π*f0.*t);
plot(t, s)

Psinal = ((0.5^2)/2)+((0.3^2)/2)+((0.15/2)/2);
B0 = 5;
σ20 = (2.0^(-2*B0))/3;
SNR = pow2db(Psinal/σ20);
println(SNR, " dB")

A = 40;
δr = 10^(-A/20);
δp = 0.05;
ωr = 0.2π;
ωp = 2π*3*f0/fa;
sig = pyimport("scipy.signal");
N, Wn = sig.ellipord(ωp/π, ωr/π, -20*log10(1-δp), -20*log10(δr));
zpkellip = digitalfilter(Lowpass(Wn), Elliptic(N, -20*log10(1-δp), -20*log10(δr)));
a = coefa(zpkellip);
b = coefb(zpkellip);
ω = 0:π;
H = freqz(zpkellip, ω);
plot(ω*fa/(2π), amp2db.(abs.(H)))
grid()
xlabel("f (Hz)");






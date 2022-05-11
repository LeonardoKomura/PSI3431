using DSP
using PlotlyJS
plotlyjs()

Ω = 2*π*500;
τ = 0.5;
f_a = 40000;
t = 0:1/f_a:2;
x0 = (sin.(Ω.*t).^3).*exp.(-t./τ);
display(plot(t, x0))

#revisar amplitude do ruído
wgn = randn((2*f_a)+1);
wgn_10dB = wgn ./ (10^(1/4));
x = x0 .+ wgn_10dB;
#display(plot(t, x))

n = 0:1:100;
h = 0.1.*sinc.(0.1.*(n.-50));
#display(plot(n, h))
hz = PolynomialRatio(h, [1]);
ω = range(0, π, length=500);
H = freqz(hz, ω);
#display(plot(ω/π, abs.(H)))

y = filt(hz, x)
y0 = filt(hz, x0)
n1 = 0:1:80000;
#fazer esse plot no msm gráfico
#display(plot(n1, y0))
#display(plot(n1, y))

#e) X[n] e y[n] tem a amplitude decrescente -> autocorrelação n é constante // ruído tem amplitude constante (WSS)
#f) DSP (Sx): calcular rx[l] = E{X[n+l]X[n]} (autocorrelação) -> Somatória{rx[l]exp(-jwl)}
#   Pm = rx[0] = (1/2pi)*integral{-pi->pi}{Sx}

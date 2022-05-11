using DSP
using Plots
plotlyjs()

N = 101;
n = 0:N-1;
L = (N-1)/2;
ωc = π/4;
h = (ωc/π) * sinc.((ωc/π)*(n.-L));
plot(n,h, line = :stem, marker = (:circle, 3), xlabel = "n", label = "h[n]")

hdf = PolynomialRatio(h, [1])
ω = range(0, π, 2500)
Hd = freqresp(hdf, ω)
plot(ω / π, abs.(Hd), label = "|Hd|", xlabel = "ω/π")

h2 = ((π-ωc)/π) * ((-1).^n) .* sinc.(((π-ωc)/π)*(n.-L));
plot(n,h2, line = :stem, marker = (:circle, 3), xlabel = "n", label = "h2[n]")

hdf2 = PolynomialRatio(h2, [1])
Hd2 = freqresp(hdf2, ω)
plot(ω / π, abs.(Hd2), label = "|Hd2|", xlabel = "ω/π")

function filtro(h, x)
    K = length(x);
    y_loop = zeros(K);
    for i in 1:K
        for j in 0:i-1
            y_loop[i] = (h[j+1]*x[i-j]) + y_loop[i];
        end
    end
    return y_loop
end

x = cos.(π.*n/20) .+ cos.(π.*n/3);
plot(n,x, xlabel = "n", label = "x[n]")

y_pb = filt(hdf, x);
y_pb_f = filtro(h, x);
plot(n,y_pb, label = "y_filt")
plot!(n, y_pb_f, label = "y_filt_2.")

y_pa = filt(hdf2, x);
y_pa_f = filtro(h2, x)
plot(n, y_pa, label = "y_filt")
plot!(n, y_pa_f, label = "y_filtro_2.")

#KAISER
δp = 0.005; 
δr = 0.001;
A = -20*log10(min(δp, δr)); # A = 60
β = 0.1102*(A − 8.7);

ωp = π/4;
ωr = 5π/16;
ωc = (ωp_pb + ωr_pb)/2;
Δω = abs(ωp - ωr);

N = ((A-8)/(2.285*Δω)) + 1;
N = Int64(ceil(N));
L=(N-1)/2;

n=0:N-1;
hk_pb=(ωc/π)*sinc.((ωc/π)*(n.-L)).*kaiser(N,β/π)
plot(n,hk_pb, line = :stem, marker = (:circle, 3), xlabel = "n", label = "hk_pb[n]")

hdfk_pb = PolynomialRatio(hk_pb, [1])
ω = range(0, π, 2500)
Hd = freqresp(hdfk_pb, ω)
plot(ω / π, abs.(Hd), label = "|Hd|", xlabel = "ω/π")

hk_pa=((π-ωc)/π) * ((-1).^n) .* sinc.(((π-ωc)/π)*(n.-L)).*kaiser(N,β/π)
plot(n,hk_pa, line = :stem, marker = (:circle, 3), xlabel = "n", label = "hk_pa[n]")

hdfk_pa = PolynomialRatio(hk_pa, [1])
Hd = freqresp(hdfk_pa, ω)
plot(ω / π, abs.(Hd), label = "|Hd|", xlabel = "ω/π")

N = 101;
n = 0:N-1;
yk_pb = filt(hdfk_pb, x);
plot(n, yk_pb)

yk_pa = filt(hdfk_pa, x);
plot(n, yk_pa)

# Mín-Máx
N = ceil(Int,(-10*log10(δp*δr)-13)/(2.324*Δω))+3;
hmm_pb = remez(N, [(0, ωp/2π) => (1, 1), (ωr/2π, 0.5) => (0, δp/δr)]);
L = (N-1)/2;
n = 0:N-1;
plot(n,hmm_pb, line = :stem, marker = (:circle, 3), xlabel = "n", label = "hmm_pb[n]")

hdfmm_pb = PolynomialRatio(hmm_pb, [1])
Hd = freqresp(hdfmm_pb, ω)
plot(ω / π, abs.(Hd), label = "|Hd|", xlabel = "ω/π")

N = ceil(Int,(-10*log10(δp*δr)-13)/(2.324*Δω))+3;
hmm_pa = remez(N, [(0, (ωp)/2π) => (0, δp/δr), ((ωr)/2π, 0.5) => (1, 1)]);
plot(n,hmm_pa, line = :stem, marker = (:circle, 3), xlabel = "n", label = "hmm_pa[n]")

hdfmm_pa = PolynomialRatio(hmm_pa, [1])
Hd = freqresp(hdfmm_pa, ω)
plot(ω / π, abs.(Hd), label = "|Hd|", xlabel = "ω/π")

N = 101;
n = 0:N-1;
ymm_pb = filt(hdfmm_pb, x);
plot(n, ymm_pb)

ymm_pa = filt(hdfmm_pa, x);
plot(n, ymm_pa)



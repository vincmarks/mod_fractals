
include("setup.jl")

N = 512
KMIN = 4
KMAX = N ÷ 4          # Nyquist-Bereich meiden
OUT = "out"
MAKE_AUDIO = true          

isdir(OUT) || mkdir(OUT)

#######
# default farben
C_BLUE, C_ORANGE, C_RED = "#2a78d6", "#eb6834", "#e34948"
INK2, MUTED, SURFACE    = "#52514e", "#898781", "#fcfcfb"

default(fontfamily = "sans-serif", background_color = SURFACE,
        foreground_color_axis = "#c3c2b7", foreground_color_text = INK2,
        gridcolor = "#e1e0d9", gridalpha = 0.9, framestyle = :axes,
        legendfontsize = 7, titlefontsize = 9, guidefontsize = 8,
        tickfontsize = 7, linewidth = 2, legend_foreground_color = :transparent)

zscore(x) = (x .- mean(x)) ./ std(x)

###############
# Wie fBm bei verschiedenen β aussieht

betas_demo = [1.0, 2.0, 3.0]
fields = Dict(b => synth_fbm(N, b; rng = MersenneTwister(7)) for b in betas_demo)

panels = Plots.Plot[]
for b in betas_demo
    f = zscore(fields[b])
    push!(panels, heatmap(clamp.(f, -2.6, 2.6); c = :grays, aspect_ratio = :equal,
                          axis = false, ticks = false, colorbar = false,
                          title = "β = $(Int(b))", xlims = (1, N), ylims = (1, N)))
end

# 1D ausschnitt
for b in betas_demo
    f = zscore(fields[b])
    push!(panels, plot(f[N ÷ 2, :]; color = C_BLUE, linewidth = 1.1,
                       label = "", ylims = (-4, 4), xlims = (0, N),
                       xlabel = "Position (Pixel)",
                       ylabel = b == betas_demo[1] ? "Amplitude (σ)" : ""))
end

fig1 = plot(panels...; layout = grid(2, 3, heights = [0.62, 0.38]),
            size = (1400, 900),
            plot_title = "FBS für unterschiedliches beta ",
            plot_titlefontsize = 12)
savefig(fig1, joinpath(OUT, "01_fbm_felder.png"))

##########
# Leistungsspektrum eines einzelnen Feldes + Fit, dazu die Summe/Mittelwert-Falle
field = fields[2.0]

# rohes Periodogramm als Punktwolke -- nur eine Stichprobe der Punkte,
# sonst wird der Plot unbrauchbar (N² Punkte bei N=512!)
let
    f = field .- mean(field)
    P = abs2.(fft(f)) ./ length(f)
    kall = kgrid(N)
    sel = (kall .>= 1) .& (kall .<= N ÷ 2)
    kv, Pv = kall[sel], P[sel]
    idx = randperm(MersenneTwister(3), length(kv))[1:12_000]

    k, S = radial_spectrum(field)
    β, icept, r2 = fit_beta(k, S; kmin = KMIN, kmax = KMAX)

    p1 = scatter(kv[idx], Pv[idx]; markersize = 0.8, markeralpha = 0.16,
                 markerstrokewidth = 0, color = MUTED,
                 label = "rohes Periodogramm |F(k)|²",
                 xscale = :log10, yscale = :log10)
    plot!(p1, k, S; color = C_BLUE, label = "radial gemittelt S(k)")
    kf = [KMIN, KMAX]
    plot!(p1, kf, (10 ^ icept) .* kf .^ (-β); color = C_ORANGE, linestyle = :dash,
          label = @sprintf("Fit: β = %.3f  (R² = %.4f)", β, r2))
    plot!(p1; xlabel = "Wellenzahl k (Moden pro Bildbreite)", ylabel = "Leistung S(k)",
          legend = :bottomleft)



    fig2 = plot(p1;
                plot_title = "Vom Periodogramm zum Exponenten", plot_titlefontsize = 12)
    savefig(fig2, joinpath(OUT, "02_spektrum_fit.png"))


end

######
# Allgemeiner Test, ob das passiert was passieren soll

# Kalibrierdiagramm (kommt β_ein wieder raus?) + Leckage-Test
betas_in = 0.5:0.25:4.0
n_real   = 30 # anzahl der mittelungen pro wert (Monte Carlo artiger ansatz) 
mu = Float64[]
sd = Float64[]
for b in betas_in
    vals = [measure(synth_fbm(N, b; rng = MersenneTwister(1000 + round(Int, b * 100) + r));
                    kmin = KMIN, kmax = KMAX)[1] for r in 1:n_real]
    push!(mu, mean(vals))
    push!(sd, std(vals))
end
mae = mean(abs.(mu .- collect(betas_in)))

# Leckage-Test: echte Daten sind nicht periodisch, unsere synth_fbm-Felder
# schon (sie kommen ja direkt aus einer FFT). Um das nachzustellen, schneiden
# wir ein Stück aus einem größeren Feld heraus -- am Rand passt dann nichts
# mehr nahtlos zusammen, genau wie bei einem echten Satellitenbild-Ausschnitt.
betas_leak = 1.0:0.5:5.0
with_w    = Float64[]
without_w = Float64[]
for b in betas_leak
    big  = synth_fbm(1536, b; rng = MersenneTwister(500 + round(Int, b * 10)))
    crop = big[301:300+N, 701:700+N]
    push!(with_w,    measure(crop; kmin = KMIN, kmax = KMAX, window = true)[1])
    push!(without_w, measure(crop; kmin = KMIN, kmax = KMAX, window = false)[1])
end

p3 = plot([0.3, 4.2], [0.3, 4.2]; color = MUTED, linestyle = :dot, linewidth = 1.2,
          label = "ideal");
scatter!(p3, collect(betas_in), mu; yerror = sd, color = C_BLUE, markersize = 4,
         markerstrokecolor = C_BLUE, label = "gemessen (n = $n_real)");
plot!(p3; xlims = (0.3, 4.2), ylims = (0.3, 4.2), xlabel = "β hineingesteckt",
      ylabel = "β zurückgemessen", legend = :topleft,
)

p4 = plot(collect(betas_leak), collect(betas_leak); color = MUTED, linestyle = :dot,
          linewidth = 1.2, label = "ideal");
plot!(p4, collect(betas_leak), with_w;    color = C_BLUE, marker = :circle,
      markersize = 4, label = "mit Hann-Fenster");
plot!(p4, collect(betas_leak), without_w; color = C_RED, marker = :square,
      markersize = 4, label = "ohne Fenster");
plot!(p4; xlabel = "β hineingesteckt", ylabel = "β zurückgemessen", legend = :topleft,
      title = "Nicht-periodischer Ausschnitt: spektrale Leckage")

fig3 = plot(p3, p4; layout = (1, 2), size = (1400, 620),
            plot_title = "Verifikation der Pipeline", plot_titlefontsize = 12)
savefig(fig3, joinpath(OUT, "03_kalibrierung.png"))

######### 
# 3D-Ansicht
surfs = Plots.Plot[]
for b in (2.0, 3.0)
    Z = zscore(fields[b])[1:4:end, 1:4:end]
    push!(surfs, surface(Z; c = :terrain, colorbar = false, camera = (-58, 48),
                         axis = false, ticks = false, title = "β = $(Int(b))",
                         linewidth = 0))
end
fig4 = plot(surfs...; layout = (1, 2), size = (1400, 600),
            plot_title = "", plot_titlefontsize = 12)
savefig(fig4, joinpath(OUT, "04_landschaft_3d.png"))

########
# Ton Bsp.
if MAKE_AUDIO
    sr = 44100 # samplerate
    for (b, name) in ((0.0, "weiss"), (1.0, "rosa"), (2.0, "braun"))
        x = synth_fbm1d(sr * 3, b; rng = MersenneTwister(9)) # 3 Sekunden
        x = x ./ maximum(abs.(x)) .* 0.75 # normierung
        nf = sr ÷ 20  # fade in/out über 50 ms
        x[1:nf]         .*= range(0, 1; length = nf)
        x[end-nf+1:end] .*= range(1, 0; length = nf)
        wavwrite(x, joinpath(OUT, @sprintf("klang_beta%d_%s.wav", Int(b), name));
                 Fs = sr)
    end
end

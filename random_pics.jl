include("setup.jl")
using JpegTurbo
using ImageCore: Gray, gray

TAR   = joinpath(@__DIR__, "random_pics.tar")
OUT   = "out"
N     = 512                # Kantenlänge des analysierten Ausschnitts
KMIN  = 10
KMAX  = N ÷ 4              # Nyquist-nahen, verrauschten Bereich meiden
NPICS = 3_500                  # zum Anschauen der Verteilung einfach hochsetzen
SEED  = 42

VERBOSE = NPICS ≤ 50       # bei großen Stichproben nur noch Fortschritt melden

isdir(OUT) || mkdir(OUT)

########
# Tar-Archiv lesen, ohne es auszupacken
#
# Das Archiv ist 2.3 GB groß und enthält 36 500 Bilder -- alles auszupacken
# (oder pro Bild einmal komplett durchzuscannen) wäre Verschwendung, wenn man
# am Ende nur ein paar zufällige Bilder braucht.
#
# tar ist zum Glück trivial aufgebaut: vor jedem Datei-Inhalt steht ein
# 512-Byte-Header, der Inhalt selbst ist auf ein Vielfaches von 512 Bytes
# aufgefüllt. Man kann also einmal von Header zu Header springen (liest nur
# ~18 MB statt 2.3 GB) und sich merken, wo welche Datei liegt. Danach ist
# jedes einzelne Bild ein seek + read.

struct TarEntry
    name::String
    offset::Int    # Byte-Position, an der der Datei-Inhalt beginnt
    size::Int
end

"Feld fester Breite aus einem tar-Header, ohne die auffüllenden Nullbytes."
field(hdr, range) = rstrip(String(@view hdr[range]), ['\0', ' '])

"""
    tar_index(path; ext) -> Vector{TarEntry}

Inhaltsverzeichnis eines (unkomprimierten) tar-Archivs: Name, Offset und
Größe jeder regulären Datei mit der passenden Endung.
"""
function tar_index(path::AbstractString; ext::AbstractString = ".jpg")
    entries = TarEntry[]
    hdr = Vector{UInt8}(undef, 512)

    open(path) do io
        while !eof(io)
            pos = position(io)
            readbytes!(io, hdr, 512) == 512 || break
            all(iszero, hdr) && break          # zwei Nullblöcke = Archivende

            name = field(hdr, 1:100)
            pre  = field(hdr, 346:500)         # ustar-Präfix bei langen Pfaden
            isempty(pre) || (name = pre * "/" * name)

            szf  = field(hdr, 125:136)         # Größe steht oktal im Header
            sz   = isempty(szf) ? 0 : parse(Int, szf; base = 8)
            typ  = Char(hdr[157])              # '0'/'\0' = reguläre Datei

            if (typ == '0' || typ == '\0') && endswith(name, ext)
                push!(entries, TarEntry(name, pos + 512, sz))
            end

            skip(io, cld(sz, 512) * 512)       # Datei-Inhalt überspringen
        end
    end
    return entries
end

"Rohe Bytes eines Archiv-Eintrags."
function read_entry(io::IO, e::TarEntry)
    seek(io, e.offset)
    return read(io, e.size)
end

########
# Bild -> analysierbares Feld

"""
    to_gray(bytes) -> Matrix{Float64}

JPEG-Bytes als Graustufenbild in [0, 1]. `Gray` gewichtet die Kanäle nach
Rec.601 (0.299 R + 0.587 G + 0.114 B), also nach wahrgenommener Helligkeit
statt als naiver Kanalmittelwert.
"""
to_gray(bytes::Vector{UInt8}) = Float64.(gray.(jpeg_decode(Gray, bytes)))

"""
    center_crop(A, n) -> Matrix oder nothing

Mittiger quadratischer Ausschnitt der Kantenlänge `n`, `nothing` falls das
Bild dafür zu klein ist.

Bewusst ein Ausschnitt und kein Skalieren: Herunterskalieren ist ein
Tiefpassfilter, das genau den hochfrequenten Teil des Spektrums verbiegt,
den wir eigentlich messen wollen. Die Places365-Bilder haben ohnehin eine
kurze Kante von 512 px, der Ausschnitt kostet uns hier also nichts.
"""
function center_crop(A::AbstractMatrix, n::Integer)
    h, w = size(A)
    (h ≥ n && w ≥ n) || return nothing
    i0 = (h - n) ÷ 2 + 1
    j0 = (w - n) ÷ 2 + 1
    return A[i0:i0+n-1, j0:j0+n-1]
end

########
# Auswertung

entries = tar_index(TAR)
@printf("Archiv indiziert: %d JPEGs\n", length(entries))

order = randperm(MersenneTwister(SEED), length(entries))

names   = String[]
betas   = Float64[]
r2s     = Float64[]
spectra = Tuple{Vector{Int}, Vector{Float64}, Float64, Float64}[]  # k, S, β, icept
crops   = Matrix{Float64}[]
skipped = 0                # Bilder, die für einen 512er-Ausschnitt zu klein sind

open(TAR) do io
    for idx in order
        length(betas) < NPICS || break
        e = entries[idx]

        img  = to_gray(read_entry(io, e))
        crop = center_crop(img, N)
        if crop === nothing
            global skipped += 1
            continue
        end

        β, icept, r2, k, S = measure(crop; kmin = KMIN, kmax = KMAX, window = true)

        push!(names, basename(e.name))
        push!(betas, β)
        push!(r2s, r2)
        if length(crops) < 6                # nur die ersten paar Bilder plotten
            push!(crops, crop)
            push!(spectra, (k, S, β, icept))
        end
        if VERBOSE
            @printf("%-32s  %4d x %4d px  ->  Ausschnitt %d²  ->  β = %.3f  (R² = %.4f)\n",
                    basename(e.name), size(img, 1), size(img, 2), N, β, r2)
        elseif length(betas) % 250 == 0
            @printf("  %d / %d Bilder ausgewertet\n", length(betas), NPICS)
        end
    end
end

########
# Bild + Spektrum der analysierten Fotos

for (i, (crop, (k, S, β, icept))) in enumerate(zip(crops, spectra))
    p_img = heatmap(crop; c = :grays, aspect_ratio = :equal, yflip = true,
                    axis = false, ticks = false, colorbar = false,
                    xlims = (1, N), ylims = (1, N),
                    title = "Ausschnitt $(N)×$(N) px\n$(names[i])",
                    titlefontsize = 8)

    p_spec = plot(k, S; seriestype = :scatter, xscale = :log10, yscale = :log10,
                  markersize = 2, markerstrokewidth = 0, label = "S(k)")
    kf = [Float64(KMIN), Float64(KMAX)]
    plot!(p_spec, kf, (10 ^ icept) .* kf .^ (-β); linestyle = :dash,
          label = @sprintf("Fit: β = %.3f  (R² = %.4f)", β, r2s[i]))
    plot!(p_spec; xlabel = "Wellenzahl k", ylabel = "Leistung S(k)",
          title = "Powerspektrum", legend = :bottomleft)

    fig = plot(p_img, p_spec; layout = (1, 2), size = (1200, 520),
               bottom_margin = 6Plots.mm, top_margin = 4Plots.mm)
    savefig(fig, joinpath(OUT, @sprintf("pics_%02d_spektrum.png", i)))
end

########
# Verteilung über viele Bilder 

if NPICS > 1
    @printf("\nβ über %d Bilder (alle %d×%d px): Mittel %.3f, Median %.3f, σ = %.3f, R² im Mittel %.4f\n",
            length(betas), N, N, mean(betas), median(betas), std(betas), mean(r2s))
    @printf("übersprungen (kurze Kante < %d px): %d\n", N, skipped)

    p_hist = histogram(betas; bins = 50, label = "", normalize = :pdf,
                       size = (900, 560), left_margin = 6Plots.mm,
                       bottom_margin = 5Plots.mm, xlabel = L"\beta", ylabel = "Dichte",
                       title = "")
    label_str = @sprintf("\\langle \\beta \\rangle = %.2f", mean(betas))
    vline!(p_hist, [mean(betas)]; linestyle = :dash, linewidth = 2,
       label = latexstring(label_str))
vline!(p_hist, [2.0]; linestyle = :dot, linewidth = 2,
       label = L"\beta" * " = 2 (Toralba et al.)")
    savefig(p_hist, joinpath(OUT, "pics_verteilung.png"))

    open(joinpath(OUT, "pics_beta.csv"), "w") do io
        println(io, "datei,beta,r2")
        for (n, b, r) in zip(names, betas, r2s)
            @printf(io, "%s,%.6f,%.6f\n", n, b, r)
        end
    end
end

# Reproducibility Repo: Statistische Analyse natürlicher Daten

Erste Experimente zu Thema 1 (fraktale Datenanalyse): 2D-fBm-Felder mit
bekanntem Spektralexponenten β synthetisieren, β über das radial gemittelte
Leistungsspektrum zurückmessen und die Pipeline mit einem Kalibrierdiagramm
gegenprüfen, bevor es an echte Beobachtungsdaten geht.

- [`fbm_pipeline.jl`](fbm_pipeline.jl) -- Synthese, Spektrum, Fit
- [`run_analysis.jl`](run_analysis.jl) -- erzeugt alle Abbildungen in `out/`

```
julia run_analysis.jl
```


# Functions for ISSR processing

# Constants
Md = 0.0289644 # Molar Masse Dry Air
Mv = 0.01801528 # Molar Mass Water Vapour

# build rhi from q and T

function rhi_calc(q, T, p)
    epsi = Mv/Md
    #Saturation pressure over ice as in Murphy & Koop:
    psi = exp.(9.550426 .- 5723.265 ./ T .+ 3.53068 .* log.(T) .- 0.00728332 .* T)
    Si = (p .* q ) ./ (epsi .* psi ) # Saturation ratio over ice
    return Si
end


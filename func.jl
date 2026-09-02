# Functions for ISSR processing

# Constants
Md = 0.0289644 # Molar Masse Dry Air
Mv = 0.01801528 # Molar Mass Water Vapour

# build rhi from q and T

function rhi_calc(q, T)
    epsi = Mv/Md
    #Saturation pressure over ice as in Murphy & Koop:
    psi = e^(9.550426 - 5723.265/T + 3.53068*ln(T) -0.00728332*T)
    Si = (p*s)/(epsi*psi) # Saturation ratio over ice
    return Si
end
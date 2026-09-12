import math

class DoublyReinforcedBeam:
    def __init__(self, b: float, d: float, d_prime: float, fc: float, fy: float, As: float, As_prime: float):
        """
        Parameters:
        b        : Beam width (in)
        d        : Depth from top compression fiber to centroid of tension steel (in)
        d_prime  : Depth from top compression fiber to centroid of compression steel (in)
        fc       : Concrete compressive strength (psi)
        fy       : Yield strength of steel (psi)
        As       : Total area of tension reinforcement (in^2)
        As_prime : Area of compression reinforcement (in^2)
        """
        self.b = b
        self.d = d
        self.d_prime = d_prime
        self.fc = fc
        self.fy = fy
        self.As = As
        self.As_prime = As_prime
        
        self.Es = 29000000.0  # Modulus of elasticity of steel (psi)
        self.ey = fy / self.Es # Yield strain of steel
        self.beta1 = self._get_beta1()

    def _get_beta1(self) -> float:
        """Calculate ACI 318 beta_1 factor based on f'c."""
        if self.fc <= 4000:
            return 0.85
        elif self.fc >= 8000:
            return 0.65
        else:
            return 0.85 - (0.05 * (self.fc - 4000) / 1000)

    def analyze_section(self) -> dict:
        """Solves equilibrium using quadratic formula to determine neutral axis (c)."""
        # Assumptions: Tension steel yields (fs = fy)
        # Equilibrium: T = Cc + Cs
        # As*fy = 0.85*fc*beta1*c*b + As'*Es*[(c - d')/c * 0.003 - 0.85*fc/Es]
        # Formulating quadratic equation: A*c^2 + B*c + C = 0
        
        A = 0.85 * self.fc * self.beta1 * self.b
        B = self.As_prime * (0.003 * self.Es - 0.85 * self.fc) - (self.As * self.fy)
        C = -self.As_prime * 0.003 * self.Es * self.d_prime

        # Quadratic formula root
        c = (-B + math.sqrt(B**2 - 4 * A * C)) / (2 * A)
        a = self.beta1 * c

        # Calculate strains
        e_s_prime = 0.003 * (c - self.d_prime) / c  # Compression steel strain
        e_t = 0.003 * (self.d - c) / c               # Tension steel strain

        # Calculate stresses
        fs_prime = min(e_s_prime * self.Es, self.fy) # Compression steel stress capped at fy
        fs = min(e_t * self.Es, self.fy)             # Tension steel stress

        # Re-evaluate forces for precise moment calculation
        Cc = 0.85 * self.fc * a * self.b
        Cs = self.As_prime * (fs_prime - 0.85 * self.fc)
        T = self.As * fs

        # Nominal Moment Capacity (Mn) taking moments about tension steel
        Mn_lb_in = Cc * (self.d - a / 2.0) + Cs * (self.d - self.d_prime)
        Mn_kip_ft = (Mn_lb_in / 1000.0) / 12.0

        # Strength reduction factor (phi) under ACI 318-19
        if e_t >= (self.ey + 0.003):
            phi = 0.90 # Tension-controlled
        elif e_t <= self.ey:
            phi = 0.65 # Compression-controlled
        else:
            phi = 0.65 + (e_t - self.ey) * (0.25 / 0.003) # Transition zone

        phi_Mn_kip_ft = phi * Mn_kip_ft

        return {
            "c_in": round(c, 2),
            "a_in": round(a, 2),
            "strain_tension_et": round(e_t, 5),
            "strain_comp_es_prime": round(e_s_prime, 5),
            "comp_steel_yielded": e_s_prime >= self.ey,
            "phi": round(phi, 3),
            "Mn_kip_ft": round(Mn_kip_ft, 2),
            "phi_Mn_kip_ft": round(phi_Mn_kip_ft, 2)
        }

# Example Usage:
if __name__ == "__main__":
    beam = DoublyReinforcedBeam(
        b=12.0,          # Width = 12 inches
        d=21.5,          # Depth to tension steel = 21.5 inches
        d_prime=2.5,     # Depth to compression steel = 2.5 inches
        fc=4000.0,       # Concrete strength = 4000 psi
        fy=60000.0,      # Steel yield strength = 60000 psi
        As=5.08,         # Tension steel area (4 x #10 bars) = 5.08 in^2
        As_prime=1.57    # Compression steel area (2 x #8 bars) = 1.57 in^2
    )
    
    results = beam.analyze_section()
    
    print("--- DESIGN RESULTS ---")
    for key, value in results.items():
        print(f"{key}: {value}")

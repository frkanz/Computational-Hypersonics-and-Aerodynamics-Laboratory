### A Pluto.jl notebook ###
# v0.12.17

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ f1295720-40c7-11eb-12be-8bd0509deca0
begin
	using Markdown
	using InteractiveUtils
	using LinearAlgebra
end

# ╔═╡ 0a3822a0-40c8-11eb-342c-a77110d06445
begin
    using PlutoUI
    using Plots
    using Printf
    using Interact
end

# ╔═╡ 142d35f0-4192-11eb-244d-8db5ddc08a84
md"_Compressible Similarity (Interactive), version 1_"

# ╔═╡ 024ab760-40c8-11eb-3343-2f5d5cd7535e
md"""

## **Compressible Similarity Solution** 

#### **Description:**

This notebook computes the compressible version of Falkner-Skan Similarity Solution coupled with energy equation

		selfsimilar(M∞, T∞, ηmax, N, itermax, ϵProfile, ϵBC)

If the arguments are missing, it will use the default values.
    
		selfsimilar(M∞=1, T∞=300, ηmax=10, N=50, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)

#### **Compressible Similarity Equation**
Boundary-layer velocity and temperature profiles on the flat plate can be projected onto single profile wich is self-similar profile. It can be represented using the ordinary differential equations (ODEs) below:

$$(cf'')'+ff'' =0$$
$$(a_1g'+a_2f'f'')'+fg'=0$$
    
where 

$$f'=\frac{u}{u_e}$$
$$c=\frac{\rho \mu}{\rho_e \mu_e}$$
$$g=\frac{H}{H_e}$$ 
$$a_1=\frac{c}{\sigma}$$
$$a_2=\frac{(\gamma-1)M^2}{1+(\frac{\gamma-1}{2})M^2}\left(1-\frac{1}{\sigma}\right)c$$


and H is the enthalpy, γ is the ratio of specific heats, M is the edge Mach number, and σ is the Prandtl number. σ and M can be defined as

$$M=\frac{u_e}{\sqrt{\gamma \mathfrak{R}T_e}}$$
$$\sigma=\frac{\mu c_p}{k}$$
    
In this code, σ is assumed as 0.72. The viscosity μ is a function of T and it is calculated as

\begin{equation}
    μ = c₁\frac{T^{3/2}}{T+c₂}    
\end{equation}

c₂ is 110.4 Kelvin. c₁ is disappearing on the nondimensionalizing process. The boundary conditions for the system of ODEs are
    
$$y=0;  f=f'=0$$
$$y\rightarrow \infty;  f',g \rightarrow 0$$

The resultant equations along with the boundary conditions are solved with the Runge-Kutta Fehlberg scheme with Newton's iteration method for missing boundary condition.

Details of RK Fehlberg:
Numerical Recipes, Cambridge

Details of Similarity solution formulation:
Boundary-Layer Theory, 7ᵗʰ edition, Schlichting

Feel free to ask questions!

*Furkan Oz*

foz@okstate.edu
"""

# ╔═╡ 295debb0-40c8-11eb-1d32-a37eb331b81c
function Y1(η,y₁,y₂,y₃,y₄,y₅)
    return y₂
end

# ╔═╡ 32e77ac0-40c8-11eb-3f9a-313f3cb78e79
function Y2(η,y₁,y₂,y₃,y₄,y₅)
    return y₃
end

# ╔═╡ 373ec6f0-40c8-11eb-06c9-df8d801f6678
function Y3(η,y₁,y₂,y₃,y₄,y₅,cμ,T∞)
    return -y₃*((y₅/(2*(y₄)))-(y₅/(y₄+cμ/T∞)))-y₁*y₃*((y₄+cμ/T∞)/(sqrt(y₄)*(1+cμ/T∞)));
end

# ╔═╡ 3c949760-40c8-11eb-1ccd-75955be21825
function Y4(η,y₁,y₂,y₃,y₄,y₅)
    return y₅
end

# ╔═╡ 41c04a90-40c8-11eb-23a8-13a2fef3fb74
function Y5(η,y₁,y₂,y₃,y₄,y₅,cμ,T∞,Pr,γ,M∞)
    return -y₅^2*((0.5/y₄)-(1/(y₄+cμ/T∞)))-Pr*y₁*y₅/sqrt(y₄)*(y₄+cμ/T∞)/(1+cμ/T∞)-(γ-1)*Pr*M∞^2*y₃^2;
end

# ╔═╡ 132d8800-40c8-11eb-04f8-350ebed9c8ee
"""
Computes the compressible Similarity Solution
This notebook computes the compressible Similarity Solution

		selfsimilar(M∞, T∞, ηmax, N, itermax, ϵProfile, ϵBC)

If the arguments are missing, it will use the default values.
    
		selfsimilar(M∞=1, T∞=300, ηmax=10, N=50, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)

Furkan Oz,
foz@okstate.edu, 
    
"""
function selfsimilar(M∞=1, T∞=300, ηmax=10, N=50, itermax=40, ϵProfile=1e-6, ϵBC=1e-6)
    Δη = ηmax/N
    Δη²= Δη^2
    γ  = 1.4   # Ratio of specific heat
    cμ = 110.4 # Sutherland law coefficient for [Kelvin]
    Pr = 0.72  # Prandtl Number
    
    Δ = 1e-7     # Small number for shooting method (Decrease when ϵ is decreased)

    # Initializing the solution vectors
    y₁ = zeros(N+1)   # f
    y₂ = zeros(N+1)   # f'
    y₃ = zeros(N+1)   # f''
    y₄ = zeros(N+1)   # ρ(η)
    y₅ = zeros(N+1)   # ρ(η)'
    η = [(i-1)*Δη for i=1:N+1]

    # Boundary Conditions for Isothermal Case
    y₁[1] = 0
    y₂[1] = 0
    y₅[1] = 0

    α₀ = 0.1 	  # Change that if code does not converge
    β₀ = 3.0      # Change that if code does not converge
    y₃[1] = α₀    # Initial Guess
    y₄[1] = β₀    # Initial Guess
    
    # We want
    # y2[N+1] = 1
    # y4[N+1] = 1
    
	# RK Fehlberg Coefficients see reference for details
    b₂ = 0.2; b₃ = 0.3; 
    b₄ = 0.6; b₅ = 1.0; b₆ = 7.0/8.0;
    a₂₁= 0.2; 
    a₃₁= 3.0/40.0; a₃₂= 9.0/40.0;
    a₄₁= 0.3; a₄₂=-0.9; a₄₃= 6.0/5.0;
    a₅₁=-11.0/54.0;a₅₂= 2.5;  a₅₃=-70.0/27.0; a₅₄= 35.0/27.0;
    a₆₁= 1631.0/55296.0; a₆₂= 175.0/512.0; a₆₃= 575.0/13824.0; 
	a₆₄= 44275.0/110592.0; a₆₅= 253.0/4096.0;
    c₁ = 37.0/378.0; c₂ = 0; c₃ = 250.0/621.0; c₄ = 125.0/594.0;
    c₅ = 0; c₆ = 512.0/1771.0; 
    
    iter = 0
    errorProfile = 1.
    errorBC = 1.
    normₒ = 0
    normₙ = 0
    
    while ϵProfile<=errorProfile && iter<itermax
        y₃[1] = α₀;    # Initial Guess
        y₄[1] = β₀;    # Initial Guess
        
		# First solution for Newton's iteration
        for i=1:N
            
            w₁ = Y5(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i],cμ,T∞,Pr,γ,M∞);
            w₂ = Y5(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₂₁*w₁),cμ,T∞,Pr,γ,M∞);
            w₃ = Y5(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₃₁*w₁+a₃₂*w₂),cμ,T∞,Pr,γ,M∞);
            w₄ = Y5(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃),cμ,T∞,Pr,γ,M∞);
            w₅ = Y5(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄),cμ,T∞,Pr,γ,M∞);
            w₆ = Y5(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅),cμ,T∞,Pr,γ,M∞);

            y₅[i+1] = y₅[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆); 
            
            w₁ = Y4(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y4(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₂₁*w₁), y₅[i]);
            w₃ = Y4(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₅[i]);
            w₄ = Y4(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₅[i]);
            w₅ = Y4(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₅[i]);
            w₆ = Y4(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₅[i]);

            y₄[i+1] = y₄[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

            w₁ = Y3(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i],cμ,T∞);
            w₂ = Y3(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i]+Δη*(a₂₁*w₁), y₄[i], y₅[i],cμ,T∞);
            w₃ = Y3(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₄[i], y₅[i],cμ,T∞);
            w₄ = Y3(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₄[i], y₅[i],cμ,T∞);
            w₅ = Y3(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₄[i], y₅[i],cμ,T∞);
            w₆ = Y3(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₄[i], y₅[i],cμ,T∞);

            y₃[i+1] = y₃[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆); 

            w₁ = Y2(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y2(η[i]+Δη*b₂, y₁[i], y₂[i]+Δη*(a₂₁*w₁), y₃[i], y₄[i], y₅[i]);
            w₃ = Y2(η[i]+Δη*b₃, y₁[i], y₂[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₃[i], y₄[i], y₅[i]);
            w₄ = Y2(η[i]+Δη*b₄, y₁[i], y₂[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₃[i], y₄[i], y₅[i]);
            w₅ = Y2(η[i]+Δη*b₅, y₁[i], y₂[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₃[i], y₄[i], y₅[i]);
            w₆ = Y2(η[i]+Δη*b₆, y₁[i], y₂[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₃[i], y₄[i], y₅[i]);

            y₂[i+1] = y₂[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

            w₁ = Y1(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y1(η[i]+Δη*b₂, y₁[i]+Δη*(a₂₁*w₁), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₃ = Y1(η[i]+Δη*b₃, y₁[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₄ = Y1(η[i]+Δη*b₄, y₁[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₅ = Y1(η[i]+Δη*b₅, y₁[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₆ = Y1(η[i]+Δη*b₆, y₁[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₂[i], y₃[i], y₄[i], y₅[i]);

            y₁[i+1] = y₁[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

        end
        
        # Storing the freestream values for Newton's iteration method
        y₂ₒ = y₂[N+1];
        y₄ₒ = y₄[N+1];
        
        # Small number addition for Newton's iteration method
        y₃[1] = α₀+Δ;  # Initial Guess + Small number
        y₄[1] = β₀;    # Initial Guess
        
		# Second solution for Newton's iteration
        for i=1:N
            
            w₁ = Y5(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i],cμ,T∞,Pr,γ,M∞);
            w₂ = Y5(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₂₁*w₁),cμ,T∞,Pr,γ,M∞);
            w₃ = Y5(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₃₁*w₁+a₃₂*w₂),cμ,T∞,Pr,γ,M∞);
            w₄ = Y5(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃),cμ,T∞,Pr,γ,M∞);
            w₅ = Y5(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄),cμ,T∞,Pr,γ,M∞);
            w₆ = Y5(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅),cμ,T∞,Pr,γ,M∞);

            y₅[i+1] = y₅[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆); 

            w₁ = Y4(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y4(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₂₁*w₁), y₅[i]);
            w₃ = Y4(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₅[i]);
            w₄ = Y4(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₅[i]);
            w₅ = Y4(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₅[i]);
            w₆ = Y4(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₅[i]);

            y₄[i+1] = y₄[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

            w₁ = Y3(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i],cμ,T∞);
            w₂ = Y3(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i]+Δη*(a₂₁*w₁), y₄[i], y₅[i],cμ,T∞);
            w₃ = Y3(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₄[i], y₅[i],cμ,T∞);
            w₄ = Y3(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₄[i], y₅[i],cμ,T∞);
            w₅ = Y3(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₄[i], y₅[i],cμ,T∞);
            w₆ = Y3(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₄[i], y₅[i],cμ,T∞);

            y₃[i+1] = y₃[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆); 

            w₁ = Y2(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y2(η[i]+Δη*b₂, y₁[i], y₂[i]+Δη*(a₂₁*w₁), y₃[i], y₄[i], y₅[i]);
            w₃ = Y2(η[i]+Δη*b₃, y₁[i], y₂[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₃[i], y₄[i], y₅[i]);
            w₄ = Y2(η[i]+Δη*b₄, y₁[i], y₂[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₃[i], y₄[i], y₅[i]);
            w₅ = Y2(η[i]+Δη*b₅, y₁[i], y₂[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₃[i], y₄[i], y₅[i]);
            w₆ = Y2(η[i]+Δη*b₆, y₁[i], y₂[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₃[i], y₄[i], y₅[i]);

            y₂[i+1] = y₂[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

            w₁ = Y1(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y1(η[i]+Δη*b₂, y₁[i]+Δη*(a₂₁*w₁), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₃ = Y1(η[i]+Δη*b₃, y₁[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₄ = Y1(η[i]+Δη*b₄, y₁[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₅ = Y1(η[i]+Δη*b₅, y₁[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₆ = Y1(η[i]+Δη*b₆, y₁[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₂[i], y₃[i], y₄[i], y₅[i]);

            y₁[i+1] = y₁[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

        end
        
        # Storing the freestream values for Newton's iteration method
        y₂ₙ₁ = y₂[N+1];
        y₄ₙ₁ = y₄[N+1];
        
        # Small number addition for Newton's iteration method
        y₃[1] = α₀;    # Initial Guess
        y₄[1] = β₀+Δ;  # Initial Guess + Small number
        
		# Third solution for Newton's iteration
        for i=1:N
            
            w₁ = Y5(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i],cμ,T∞,Pr,γ,M∞);
            w₂ = Y5(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₂₁*w₁),cμ,T∞,Pr,γ,M∞);
            w₃ = Y5(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₃₁*w₁+a₃₂*w₂),cμ,T∞,Pr,γ,M∞);
            w₄ = Y5(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃),cμ,T∞,Pr,γ,M∞);
            w₅ = Y5(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄),cμ,T∞,Pr,γ,M∞);
            w₆ = Y5(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅),cμ,T∞,Pr,γ,M∞);

            y₅[i+1] = y₅[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆); 

            w₁ = Y4(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y4(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₂₁*w₁), y₅[i]);
            w₃ = Y4(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₅[i]);
            w₄ = Y4(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₅[i]);
            w₅ = Y4(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₅[i]);
            w₆ = Y4(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i], y₄[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₅[i]);

            y₄[i+1] = y₄[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

            w₁ = Y3(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i],cμ,T∞);
            w₂ = Y3(η[i]+Δη*b₂, y₁[i], y₂[i], y₃[i]+Δη*(a₂₁*w₁), y₄[i], y₅[i],cμ,T∞);
            w₃ = Y3(η[i]+Δη*b₃, y₁[i], y₂[i], y₃[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₄[i], y₅[i],cμ,T∞);
            w₄ = Y3(η[i]+Δη*b₄, y₁[i], y₂[i], y₃[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₄[i], y₅[i],cμ,T∞);
            w₅ = Y3(η[i]+Δη*b₅, y₁[i], y₂[i], y₃[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₄[i], y₅[i],cμ,T∞);
            w₆ = Y3(η[i]+Δη*b₆, y₁[i], y₂[i], y₃[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₄[i], y₅[i],cμ,T∞);

            y₃[i+1] = y₃[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆); 

            w₁ = Y2(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y2(η[i]+Δη*b₂, y₁[i], y₂[i]+Δη*(a₂₁*w₁), y₃[i], y₄[i], y₅[i]);
            w₃ = Y2(η[i]+Δη*b₃, y₁[i], y₂[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₃[i], y₄[i], y₅[i]);
            w₄ = Y2(η[i]+Δη*b₄, y₁[i], y₂[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₃[i], y₄[i], y₅[i]);
            w₅ = Y2(η[i]+Δη*b₅, y₁[i], y₂[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₃[i], y₄[i], y₅[i]);
            w₆ = Y2(η[i]+Δη*b₆, y₁[i], y₂[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₃[i], y₄[i], y₅[i]);

            y₂[i+1] = y₂[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

            w₁ = Y1(η[i], y₁[i], y₂[i], y₃[i], y₄[i], y₅[i]);
            w₂ = Y1(η[i]+Δη*b₂, y₁[i]+Δη*(a₂₁*w₁), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₃ = Y1(η[i]+Δη*b₃, y₁[i]+Δη*(a₃₁*w₁+a₃₂*w₂), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₄ = Y1(η[i]+Δη*b₄, y₁[i]+Δη*(a₄₁*w₁+a₄₂*w₂+a₄₃*w₃), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₅ = Y1(η[i]+Δη*b₅, y₁[i]+Δη*(a₅₁*w₁+a₅₂*w₂+a₅₃*w₃+a₅₄*w₄), y₂[i], y₃[i], y₄[i], y₅[i]);
            w₆ = Y1(η[i]+Δη*b₆, y₁[i]+Δη*(a₆₁*w₁+a₆₂*w₂+a₆₃*w₃+a₆₄*w₄+a₆₅*w₅), y₂[i], y₃[i], y₄[i], y₅[i]);

            y₁[i+1] = y₁[i] + Δη*(c₁*w₁+c₂*w₂+c₃*w₃+c₄*w₄+c₅*w₅+c₆*w₆);

        end
        
        # Storing the freestream values for Newton's iteration method
        y₂ₙ₂ = y₂[N+1];
        y₄ₙ₂ = y₄[N+1];

        # Calculation of the next initial guess with Newton's iteration method
        p₁₁ = (y₂ₙ₁-y₂ₒ)/Δ;
        p₂₁ = (y₄ₙ₁-y₄ₒ)/Δ;
        p₁₂ = (y₂ₙ₂-y₂ₒ)/Δ;
        p₂₂ = (y₄ₙ₂-y₄ₒ)/Δ;
        r₁ = 1-y₂ₒ;
        r₂ = 1-y₄ₒ;
        Δα = (p₂₂*r₁-p₁₂*r₂)/(p₁₁*p₂₂-p₁₂*p₂₁);
        Δβ = (p₁₁*r₂-p₂₁*r₁)/(p₁₁*p₂₂-p₁₂*p₂₁);
        α₀ = α₀ + Δα;
        β₀ = β₀ + Δβ;
        
        # Profile change between iteration
        normₙ = norm(y₂)
        errorProfile = maximum(abs.(normₙ-normₒ))
		
		# Convergence of boundary condition
        errorBC = abs(y₂[N+1]-1.)
        iter += 1
        normₒ = normₙ
        @printf("%4.4d %16.6e %16.6e \n", iter, errorProfile, errorBC)
        
    end
    
    if errorProfile<=ϵProfile 
        println("")
        println("Solution converged!")
        println("The maximum change between consecutive profiles is less than the error criteria ϵProfile=$ϵProfile.")
    end

    if errorBC<=ϵBC
        println("")
        println("Solution for the boundary condition converged!")
        println("The difference between h(N) and h(N+1) is less than the error criteria ϵBC=$ϵBC.")
    end
    
	# Copying values for logical names
    U = y₂
    T = y₄
    
    # Integration for η --> y transformation
    y = zeros(N+1);
    for i=2:N+1
       y[i] = y[i-1] + y₄[i]*(η[i]-η[i-1]);  
    end
    y = y*sqrt(2);
    
    return η,y,U,T,N
end

# ╔═╡ 730c0e5e-4193-11eb-25d0-8394adde63bb
md"""
-------------------------------------------------------------------------
"""

# ╔═╡ 5a236880-4193-11eb-3d79-d5668e165e21
md"""
Please use the sliders to observe the change with the Mach number and temperature.
"""

# ╔═╡ 5fb1b250-40c8-11eb-10e3-1bff8856c8ae
md"""
M∞:
$(@bind Mtest Slider(0.1:0.1:10; default=1, show_value=true))

"""

# ╔═╡ b3831040-40c8-11eb-20ad-f1c8c2f59e96
md"""
T∞:
$(@bind Ttest Slider(80:20:500; default=300, show_value=true))

"""

# ╔═╡ a838d380-40ca-11eb-22c3-13a7d164bf17
ηᵢ,yᵢ,Uᵢ,Tᵢ,Nᵢ = selfsimilar(Mtest,Ttest);

# ╔═╡ f0cf33be-40c8-11eb-341d-13c3c8b2b99b
begin
	fp = plot(Uᵢ,yᵢ,
		
		label = "U",
		legend = :topleft,
		ylims = (0,40),
		xlabel = "U",
		ylabel = "y/√(νx/U)",
		linewidth = 2,
		linecolor = :black,
		markershape = :circle,
		markercolor = :red,
		)
	sp = plot(Tᵢ,yᵢ,
		
		label = "T",
		legend = :topright,
		xlims = (0,18),
		ylims = (0,40),
		xlabel = "T",
		ylabel = "y/√(νx/U)",
		linewidth = 2,
		linecolor = :black,
		markershape = :circle,
		markercolor = :red,
		)
	plot(fp,sp,title = "Similarity (M=$Mtest - T=$Ttest)",layout=2)
end

# ╔═╡ Cell order:
# ╟─142d35f0-4192-11eb-244d-8db5ddc08a84
# ╟─024ab760-40c8-11eb-3343-2f5d5cd7535e
# ╟─f1295720-40c7-11eb-12be-8bd0509deca0
# ╟─0a3822a0-40c8-11eb-342c-a77110d06445
# ╟─132d8800-40c8-11eb-04f8-350ebed9c8ee
# ╟─295debb0-40c8-11eb-1d32-a37eb331b81c
# ╟─32e77ac0-40c8-11eb-3f9a-313f3cb78e79
# ╟─373ec6f0-40c8-11eb-06c9-df8d801f6678
# ╟─3c949760-40c8-11eb-1ccd-75955be21825
# ╟─41c04a90-40c8-11eb-23a8-13a2fef3fb74
# ╟─730c0e5e-4193-11eb-25d0-8394adde63bb
# ╟─5a236880-4193-11eb-3d79-d5668e165e21
# ╟─5fb1b250-40c8-11eb-10e3-1bff8856c8ae
# ╟─b3831040-40c8-11eb-20ad-f1c8c2f59e96
# ╟─a838d380-40ca-11eb-22c3-13a7d164bf17
# ╟─f0cf33be-40c8-11eb-341d-13c3c8b2b99b

### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 0b868b30-6420-11eb-18f9-1592f7542fe6
begin
    using PlutoUI
    using Plots
    using Printf
	using SpecialFunctions
end

# ╔═╡ 4c4fff40-64c8-11eb-2c45-57f4793ee997
md"_Parabolic Partial Differential Equation, version 1_"

# ╔═╡ 56d1f810-64c8-11eb-3e0d-81c0bbafdc62
md"""
## **Parabolic Partial Differential Equation** 
#### **Description:**
This notebook computes the Parabolic Partial Differential Equation with 4 different discretization methods (2 explicit - 2 implicit)

		FTCS(Δt,Δy,ν)
		DuFortFrankel(Δt,Δy,ν)
		Laasonen(Δt,Δy,ν)
		CrankNicolson(Δt,Δy,ν)

If the arguments are missing, it will use the default values.
    
		FTCS(Δt=0.002,Δy=0.001,ν=0.000217)
		DuFortFrankel(Δt=0.002,Δy=0.001,ν=0.000217)
		Laasonen(Δt=0.002,Δy=0.001,ν=0.000217)
		CrankNicolson(Δt=0.002,Δy=0.001,ν=0.000217)

#### **Parabolic Partial Differential Equation**
In this problem, velocity profile between 2 infinite plates will be investigated. The bottom plates is moving with U₀ and the upper plate is stationary. The equation for the system can be represented by:

$$\frac{∂u}{∂t}=α\frac{∂²u}{∂x²}$$

Boundary conditions for the system:

$$u=U₀=1 \hspace{12pt} y=0$$
$$u=0    \hspace{36pt} y=h$$

where h is the distance between two plates. In this problem, it is taken as 0.04m. The partial system of equation is solved using following discritizations:

#### **Forward Time Central Space (FTCS) Method**

It is a first-order in time and second-order in space. The discritization for this partial equation:

$$u_{i}^{n+1}=u_{i}^{n}+d(u_{i+1}^{n}-2u_{i}^{n}+u_{i-1}^{n})$$

where

$$d=\frac{α(Δt)}{Δx^{2}}$$

For the stability condition:

$$\frac{α(Δt)}{Δx^{2}}=d\leqslant0.5$$

#### **DuFort-Frankel Method**

It is a second-order in time and second-order in space. The discritization for this partial equation:

$$(1+2d)u_{i}^{n+1}=(1-2d)u_{i}^{n-1}+2d(u_{i+1}^{n}+u_{i-1}^{n})$$

where

$$d=\frac{α(Δt)}{Δx^{2}}$$

For the stability condition:

$$\frac{α(Δt)}{Δx^{2}}\leqslant0.5$$

This method requires 2 sets of initial data so in this notebook, FTCS is used to obtain data set at t=Δt.

#### **Laasonen Method**

It is a first-order in time and second-order in space. The discritization for this partial equation:

$$du_{i-1}^{n+1}-(1+2d)u_{i}^{n+1}+du_{i+1}^{n+1}=-u_{i}^{n}$$

where

$$d=\frac{α(Δt)}{Δx^{2}}$$

This method is unconditionally stable.

#### **Crank Nicholson Method**

It is a second-order in time and second-order in space. The discritization for this partial equation:

$$-\frac{1}{2}du_{i-1}^{n+1}-(1+d)u_{i}^{n+1}-\frac{1}{2}du_{i+1}^{n+1}=\frac{1}{2}du_{i-1}^{n}+(1-d)u_{i}^{n}+\frac{1}{2}du_{i+1}^{n}$$

where

$$d=\frac{α(Δt)}{Δx^{2}}$$

This method is unconditionally stable.

Details of schemes:
Computational Fluid Dynamics, Volume I, Hoffmann - Chiang


Feel free to ask questions!


*Furkan Oz*
foz@okstate.edu
"""

# ╔═╡ 29683fb0-6419-11eb-05d0-2d7392bc8afa
"""
Computes the Forward Time Central Space for the Parabolic Partial Differential Equation
This notebook computes the Parabolic Partial Differential Equation
		FTCS(Δt,Δy,ν)
If the arguments are missing, it will use the default values.
    
		FTCS(Δt=0.002,Δy=0.001,ν=0.000217)
Furkan Oz,
foz@okstate.edu, 
    
"""
function FTCS(Δt=0.002,Δy=0.001,ν=0.000217)
	#ν=0.000217 # viscosity [m²/s]
	h = 0.04 # Distance between 2 plate [m]
	T = 1.25 # Total time [s]
	#Δt = 0.002 # Time increment
	Nt = T/Δt+1 # Number of time step
	Nt = round(Int32,Nt)
	
	#Δy = 0.001 # Spacing between node points
	y = 0:Δy:h # y-coordinates
	Nx = h/Δy+1 # Number of Node points
	Nx = round(Int32,Nx)
	
	U₀ = 1 # Plate speed
	d = ν*Δt/(Δy)^2
	#Initializing
	u = zeros(Nx,Nt) 
	
	# Boundary Conditions
	u[1,:] = [U₀ for i=1:Nt]
	u[Nx,:] = [0 for i=1:Nt]
	
	# Initial Conditions
	u[:,1] = [0 for i=1:Nx]
	u[1,1] = U₀
	for ii=1:Nt-1
		for i=2:Nx-1
			u[i,ii+1] = u[i,ii]+d*(u[i+1,ii]-2*u[i,ii]+u[i-1,ii])
		end
	end
	
	return u,y
	
end

# ╔═╡ 36be1c60-6442-11eb-204a-b37677c1acfd
"""
Computes the DuFort Frankel method for the Parabolic Partial Differential Equation
This notebook computes the Parabolic Partial Differential Equation
		DuFortFrankel(Δt,Δy,ν)
If the arguments are missing, it will use the default values.
    
		DuFortFrankel(Δt=0.002,Δy=0.001,ν=0.000217)
Furkan Oz,
foz@okstate.edu, 
    
"""
function DuFortFrankel(Δt=0.002,Δy=0.001,ν=0.000217)
	#ν=0.000217 # viscosity [m²/s]
	h = 0.04 # Distance between 2 plate [m]
	T = 1.25 # Total time [s]
	#Δt = 0.002 # Time increment
	Nt = T/Δt+1 # Number of time step
	Nt = round(Int32,Nt)
	
	#Δy = 0.001 # Spacing between node points
	y = 0:Δy:h # y-coordinates
	Nx = h/Δy+1 # Number of Node points
	Nx = round(Int32,Nx)
	
	U₀ = 1 # Plate speed
	d = ν*Δt/(Δy)^2
	#Initializing
	u = zeros(Nx,Nt) 
	
	# Boundary Conditions
	u[1,:] = [U₀ for i=1:Nt]
	u[Nx,:] = [0 for i=1:Nt]
	
	# Initial Conditions
	u[:,1] = [0 for i=1:Nx]
	u[1,1] = U₀
	for ii=1:1
		for i=2:Nx-1
			u[i,ii+1] = u[i,ii]+d*(u[i+1,ii]-2*u[i,ii]+u[i-1,ii])
		end
	end
	for ii=2:Nt-1
		for i=2:Nx-1
			u[i,ii+1] = ((1-2*d)*u[i,ii-1]+2*d*(u[i+1,ii]+u[i-1,ii]))/(1+2*d)
		end
	end
	return u,y
	
end

# ╔═╡ 899d6a70-6448-11eb-2870-692c6e93c421
"""
Computes the Laasonen method for the Parabolic Partial Differential Equation
This notebook computes the Parabolic Partial Differential Equation
		Laasonen(Δt,Δy,ν)
If the arguments are missing, it will use the default values.
    
		Laasonen(Δt=0.002,Δy=0.001,ν=0.000217)
Furkan Oz,
foz@okstate.edu, 
    
"""
function Laasonen(Δt=0.002,Δy=0.001,ν=0.000217)
	#ν=0.000217 # viscosity [m²/s]
	h = 0.04 # Distance between 2 plate [m]
	T = 1.25 # Total time [s]
	#Δt = 0.002 # Time increment
	Nt = T/Δt+1 # Number of time step
	Nt = round(Int32,Nt)
	
	#Δy = 0.001 # Spacing between node points
	y = 0:Δy:h # y-coordinates
	Nx = h/Δy+1 # Number of Node points
	Nx = round(Int32,Nx)
	
	U₀ = 1 # Plate speed
	d = ν*Δt/(Δy)^2
	
	#Initializing
	A = zeros(Nx,Nx)
	b = zeros(Nx)
	u = zeros(Nx,Nt)
	
	# Boundary Conditions
	A[1,1] = 1
	b[1] = U₀
	A[Nx,Nx] = 1
	b[Nx] = 0
	
	# Initial Conditions
	u[:,1] = [0 for i=1:Nx]
	u[1,1] = U₀
	
	for ii=2:Nt-1
		# Boundary Conditions
		A[1,1] = 1
		b[1] = U₀
		A[Nx,Nx] = 1
		b[Nx] = 0
		for i=2:Nx-1
			A[i,i-1] = d
			A[i,i+0] =-1-2*d
			A[i,i+1] = d
			b[i] =-u[i,ii]
		end
		ulin = A\b
		u[1,ii+1] = b[1]
		u[Nx,ii+1] = b[Nx]
		u[2:Nx-1,ii+1] = [ulin[i] for i=2:Nx-1]
	end
	return u,y
end

# ╔═╡ 9c410110-644c-11eb-3bd9-93935242ac98
"""
Computes the Crank Nicolson method for the Parabolic Partial Differential Equation
This notebook computes the Parabolic Partial Differential Equation
		CrankNicolson(Δt,Δy,ν)
If the arguments are missing, it will use the default values.
    
		CrankNicolson(Δt=0.002,Δy=0.001,ν=0.000217)
Furkan Oz,
foz@okstate.edu, 
    
"""
function CrankNicolson(Δt=0.002,Δy=0.001,ν=0.000217)
	#ν=0.000217 # viscosity [m²/s]
	h = 0.04 # Distance between 2 plate [m]
	T = 1.25 # Total time [s]
	#Δt = 0.002 # Time increment
	Nt = T/Δt+1 # Number of time step
	Nt = round(Int32,Nt)
	
	#Δy = 0.001 # Spacing between node points
	y = 0:Δy:h # y-coordinates
	Nx = h/Δy+1 # Number of Node points
	Nx = round(Int32,Nx)
	
	U₀ = 1 # Plate speed
	d = ν*Δt/(Δy)^2
	
	#Initializing
	A = zeros(Nx,Nx)
	b = zeros(Nx)
	u = zeros(Nx,Nt)
	
	# Boundary Conditions
	A[1,1] = 1
	b[1] = U₀
	A[Nx,Nx] = 1
	b[Nx] = 0
	
	# Initial Conditions
	u[:,1] = [0 for i=1:Nx]
	u[1,1] = U₀
	
	for ii=3:Nt-1
		# Boundary Conditions
		A[1,1] = 1
		b[1] = U₀
		A[Nx,Nx] = 1
		b[Nx] = 0
		for i=2:Nx-1
			A[i,i-1] =-0.5*d
			A[i,i+0] = 1+d
			A[i,i+1] =-0.5*d
			b[i] = 0.5*d*u[i+1,ii]+(1-d)*u[i,ii]+0.5*d*u[i-1,ii]
		end
		ulin = A\b
		u[1,ii+1] = b[1]
		u[Nx,ii+1] = b[Nx]
		u[2:Nx-1,ii+1] = [ulin[i] for i=2:Nx-1]
	end
	return u,y
end

# ╔═╡ 8939c23e-642a-11eb-2de6-ade5ae2f4b4d
function Analytical(t,Δy=0.001,ν=0.000217)
	#ν=0.000217 # viscosity [m²/s]
	h = 0.04 # Distance between 2 plate [m]
	y = 0:Δy:h
	Nx = h/Δy+1 # Number of Node points
	Nx = round(Int32,Nx)
	U₀ = 1 # Plate speed
	η = [y[i]/(2*sqrt(ν*t)) for i=1:Nx]
	η₁= h/(2*sqrt(ν*t))
	u = [U₀*(erfc(η[i])-erfc(2*η₁-η[i])+erfc(2*η₁+η[i])-erfc(4*η₁-η[i])+erfc(4*η₁+η[i])-erfc(8*η₁-η[i])+erfc(8*η₁+η[i])-erfc(16*η₁-η[i])) for i=1:Nx]
	return u,y
end

# ╔═╡ 70d8c3ee-641f-11eb-1d57-a943a361c397
begin
	uFTCS,yFTCS = FTCS()
	plot([uFTCS[:,1],uFTCS[:,90],uFTCS[:,180],uFTCS[:,270],uFTCS[:,360],uFTCS[:,450],uFTCS[:,540]],yFTCS,
			title = "FTCS - \\Deltat = 0.002 - \\Deltay = 0.001",
			label = ["t=0s" "t=0.18s" "t=0.36s" "t=0.54s" "t=0.72s" "t=0.90s" "t=1.08s"],
			legend = :topright,
			ylims = (0,0.04),
			xlabel = "U",
			ylabel = "y",
			linewidth = 2,
			linecolor = :black,
			markershape = :auto,
			markercolor = :red,
			)
end

# ╔═╡ 95697270-6427-11eb-02be-bf3b75fbbd33
begin
	u2,y2 = FTCS(0.005,0.001)
	
	plot([u2[:,1],u2[:,50],u2[:,100],u2[:,150],u2[:,200]],y2,
			title = "FTCS - \\Deltat = 0.005 - \\Deltay = 0.001",
			label = ["t=0s" "t=0.25s" "t=0.50s" "t=0.75s" "t=1.00s"],
			legend = :topright,
			ylims = (0,0.04),
			xlabel = "U",
			ylabel = "y",
			linewidth = 2,
			linecolor = :black,
			markershape = :auto,
			markercolor = :red,
			)
end

# ╔═╡ b1e5b270-643f-11eb-1a37-e3227c598177
begin
	ua1,ya1 = Analytical(0.18,0.001)
	uDF,yDF = DuFortFrankel()
	uLaasonen,yLaasonen = Laasonen()
	uCrankNicolson,yCrankNicolson = CrankNicolson()
	plot([ua1-uFTCS[:,90],ua1-uDF[:,90],ua1-uLaasonen[:,90],ua1-uCrankNicolson[:,90]],[ya1],
				title = "Error - t = 0.18 - \\Deltay = 0.001",
				label = ["FTCS" "D/F" "Laasonen" "CrankNicolson"],
				legend = :topright,
				xlims = (-0.005,0.01),
				ylims = (0,0.04),
				xlabel = "U",
				ylabel = "y",
				linewidth = 2,
				linecolor = :black,
				markershape = :auto,
				markercolor = :red,
				)
end

# ╔═╡ 8c11e680-6440-11eb-076a-0b1e169ef37a
begin
	ua2,ya2 = Analytical(1.08,0.001)
	plot([ua2-uFTCS[:,540],ua2-uDF[:,540],ua2-uLaasonen[:,540],ua2-uCrankNicolson[:,540]],[ya2],
				title = "Error - t = 1.08 -\\Deltay = 0.001",
				label = ["FTCS" "D/F" "Laasonen" "CrankNicolson"],
				legend = :topright,
				xlims = (-0.0005,0.0015),
				ylims = (0,0.04),
				xlabel = "U",
				ylabel = "y",
				linewidth = 2,
				linecolor = :black,
				markershape = :auto,
				markercolor = :red,
				)
end

# ╔═╡ 33a79720-644e-11eb-0eca-e794b9b87ccf
begin
	ua3,ya3 = Analytical(1.00,0.001)
	uLaasonen005,yLaasonen005 = Laasonen(0.005)
	uLaasonen01,yLaasonen01 = Laasonen(0.01)
	uLaasonen1,yLaasonen1 = Laasonen(0.1)
	uLaasonen2,yLaasonen2 = Laasonen(0.2)
	plot([ua3-uLaasonen005[:,200],ua3-uLaasonen01[:,100],ua3-uLaasonen1[:,10],ua3-uLaasonen2[:,5]],[ya2],
					title = "Error - Laasonen - t = 1.00 -\\Deltay = 0.001",
					label = ["\\Deltat=0.005" "\\Deltat=0.01" "\\Deltat=0.1" "\\Deltat=0.2"],
					legend = :topright,
					xlims = (0,0.2),
					ylims = (0,0.04),
					xlabel = "U",
					ylabel = "y",
					linewidth = 2,
					linecolor = :black,
					markershape = :auto,
					markercolor = :red,
					)
end

# ╔═╡ Cell order:
# ╟─0b868b30-6420-11eb-18f9-1592f7542fe6
# ╟─4c4fff40-64c8-11eb-2c45-57f4793ee997
# ╟─56d1f810-64c8-11eb-3e0d-81c0bbafdc62
# ╟─29683fb0-6419-11eb-05d0-2d7392bc8afa
# ╟─36be1c60-6442-11eb-204a-b37677c1acfd
# ╟─899d6a70-6448-11eb-2870-692c6e93c421
# ╟─9c410110-644c-11eb-3bd9-93935242ac98
# ╟─8939c23e-642a-11eb-2de6-ade5ae2f4b4d
# ╟─70d8c3ee-641f-11eb-1d57-a943a361c397
# ╟─95697270-6427-11eb-02be-bf3b75fbbd33
# ╟─b1e5b270-643f-11eb-1a37-e3227c598177
# ╟─8c11e680-6440-11eb-076a-0b1e169ef37a
# ╟─33a79720-644e-11eb-0eca-e794b9b87ccf

### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# ╔═╡ 57aa3170-6e51-11eb-2670-9b77065bc4da
begin
    using PlutoUI
    using Plots
    using Printf
	using SpecialFunctions
end

# ╔═╡ ad3d9470-6e82-11eb-3072-43decfc60ef9
md"_MAC Method for Lid-driven Cavity Problem, version 1_"

# ╔═╡ 484fcd20-6f23-11eb-0dc8-5d49530983a1
md"""
## **Lid-driven Cavity Problem with MAC Method** 
#### **Description:**
This notebook uses MAC method to solve incompressible 2D Navier-Stokes equations for the solution of the lid-driven cavity problem. 

		MAC(Re,h)

If the arguments are missing, it will use the default values.
    
		MAC(Re=0.002,h=0.002)

#### **Lid-driven Cavity Problem**
In this problem, there is a moving lid over a cavity. The lid is moving with a constant speed. The governing equations for the system can be written in the conservative form: 

$$\frac{∂u}{∂x}+\frac{∂u}{∂y}=0$$
$$\frac{∂u}{∂t}+\frac{∂}{∂x}(u^2+\frac{p}{ρ})+\frac{∂(uv)}{∂y}=ν(\frac{∂^2u}{∂x^2}+\frac{∂^2u}{∂y^2})$$
$$\frac{∂v}{∂t}+\frac{∂(uv)}{∂x}(v^2+\frac{p}{ρ})+\frac{∂(uv)}{∂y}=ν(\frac{∂^2v}{∂x^2}+\frac{∂^2v}{∂y^2})$$

These equations can be nondimensionalized with the following values:

$$Re=\frac{ρ_\infty u_\infty L}{μ_\infty}$$
$$t^*=\frac{t u_\infty}{L}$$
$$x^*=\frac{x}{L}$$
$$y^*=\frac{y}{L}$$
$$u^*=\frac{u}{u_\infty}$$
$$v^*=\frac{v}{u_\infty}$$
$$p^*=\frac{p}{ρ_\infty u^2_\infty}$$

After the nondimensionalization * is ignored for simplicity:

$$\frac{∂u}{∂x}+\frac{∂u}{∂y}=0$$
$$\frac{∂u}{∂t}+\frac{∂}{∂x}(u^2+p)+\frac{∂(uv)}{∂y}=\frac{1}{Re}(\frac{∂^2u}{∂x^2}+\frac{∂^2u}{∂y^2})$$
$$\frac{∂v}{∂t}+\frac{∂(uv)}{∂x}(v^2+p)+\frac{∂(uv)}{∂y}=\frac{1}{Re}(\frac{∂^2v}{∂x^2}+\frac{∂^2v}{∂y^2})$$

Since continuity term does not include a time dependent term, there is not communication mechanism between continuity and momentum equations. To overcome this problem, possion equation will be solved. In this notebook, calculations are done in three steps:

1) Calculation of u*.
2) Calculation of p.
3) Calculation of corrected u.

This approach is called pressure correction. The discritization for the system is:

$$I   = \frac{(u_{i+1,j}+u_{i,j})^2-(u_{i,j}+u_{i-1,j})^2}{4Δx}$$
$$II  = \frac{(u_{i,j+1}+u_{i,j})(v_{i+1,j}+v_{i,j})}{4Δy}-\frac{(u_{i,j}+u_{i,j-1})(v_{i+1,j-1}+v_{i,j-1})}{4Δy}$$
$$III = \frac{1}{Re}(\frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{Δx^2}+\frac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{Δy^2})$$
$$\frac{u_{i,j}^*-u_{i,j}^n}{Δt}=-I-II+III$$

For the pressure correction:

$$I=\frac{p_{i+1,j}-2p_{i,j}+p_{i-1,j}}{Δx^2}+\frac{p_{i,j+1}-2p_{i,j}+p_{i,j-1}}{Δy^2}$$
$$II=\frac{\frac{u_{i,j}^*-u_{i-1,j}^*}{Δx}+\frac{v_{i,j}^*-v_{i,j-1}^*}{Δy}}{Δt}$$
$$p_{i,j}=p_{i,j}+ω\frac{Δx^2+Δy^2}{2}(I-II)$$

In this notebook, ω is taken as 1 and this correction method is applied 10 times. Once the pressure is calculated, correct velocity can be calculated for the next time step as:

$$\frac{u^{n+1}_{i,j}-u^*_{i,j}}{Δt}=\frac{p_{i+1,j}-p_{i,j}}{Δx}$$

Same procedure can be applied for the v velocity. Boundary conditions for the system is:

$$u(0,:)=0$$
$$u(N,:)=0$$
$$u(:,0)=0$$
$$u(:,N)=1$$
$$v(0,:)=0$$
$$v(N,:)=0$$
$$v(:,0)=0$$
$$v(:,N)=0$$

For the pressure Nuemann boundary condition is applied as:

$$\frac{∂p}{∂η}=0$$

where η is the wall normal direction.

Details of schemes:
Computational Fluid Dynamics, Volume I, Hoffmann - Chiang


Feel free to ask questions!


*Furkan Oz*
foz@okstate.edu
"""

# ╔═╡ 6bc808d0-6e51-11eb-0794-b9585413d41b

"""
This notebook uses MAC method to solve incompressible 2D Navier-Stokes equations for the solution of the lid-driven cavity problem. 

		MAC(Re,h)

If the arguments are missing, it will use the default values.
    
		MAC(Re=0.002,h=0.002)
		
Furkan Oz,
foz@okstate.edu, 
    
"""
function MAC(Re=100, h=0.02)
	# Initialiazing for domain setup
	xdom = [0 1.0] # Domain limits for x
	ydom = [0 1.0] # Domain limits for y
	Δx = h
	Δy = h
	Δx²= Δx^2
	Δy²= Δy^2 
	nex = floor(Int32,(xdom[2]-xdom[1])/h) #Number of elements in x-direction
	ney = floor(Int32,(xdom[2]-xdom[1])/h) #Number of elements in y-direction
	nnx = nex + 1 #Number of nodes in x-direction
	nny = ney + 1 #Number of nodes in y-direction
	
	# Initialiazing for solver setup
	
	CFL = 1 # CFL condition (less than 1 preferred)
	Δtc = CFL*min(Δx,Δy) # Time-step based on Convective Term
	Δtv = 0.25*Re*min(Δx²,Δy²) # Time-step based on Viscous Term
	σ = 0.5
	a =-2.0/(Δx²) - 2.0/(Δy²)
	Δt =  σ*min(Δtc,Δtv) # Time-step based on lowest condition
	ϵ = 1e-5 # Convergence limit
	
	ttmax = 30000 # Maximum iteration limit
	iimax = 10 # Pressure correction iteration limit
	
	u  = zeros(nnx,nny+1)
	v  = zeros(nnx+1,nny)
	uˢ = zeros(nnx,nny+1)
	vˢ = zeros(nnx+1,nny)
	u₀ = zeros(nnx,nny+1)
	v₀ = zeros(nnx+1,nny)
	
	p  = zeros(nex+2,ney+2)
	p₀ = zeros(nex+2,ney+2)
	
	rᵤ = 0
	rᵥ = 0
	rₚ = 0
	ttend = 0
	
	rᵤₜ = zeros(ttmax)
	rᵥₜ = zeros(ttmax)
	rₚₜ = zeros(ttmax)
	
	#Solver starts
	
	for tt=1:ttmax
		for j=1:nny+1
			for i=1:nnx
				u₀[i,j] = u[i,j]
			end
		end
		for j=1:nny
			for i=1:nnx+1
				v₀[i,j] = v[i,j]
			end
		end
		for j=1:ney+2
			for i=1:nex+2
				p₀[i,j] = p[i,j]
			end
		end

		# Boundary conditions for velocity
		for i=1:nnx
			u[i,1] =-u[i,2]
			u[i,nny+1] = 2.0-u[i,nny]
		end
		for i=1:nny
			v[1,i] =-v[2,i]
			v[nnx+1,i] =-v[nnx,i]
		end
		
		# uˢ calculation
		for j=2:nny
			for i=2:nnx-1
  
  				c₃ = (1.0/Re)*((u[i+1,j]-2.0*u[i,j]+u[i-1,j])/(Δx²)  				                      + (u[i,j+1]-2.0*u[i,j]+u[i,j-1])/(Δy²) )
                   
  				c₁ = (0.25/Δx)*((u[i+1,j]+u[i,j])^2 - (u[i,j]+u[i-1,j])^2)
  
  				c₂ = (0.25/Δy)*((u[i,j+1]+u[i,j])*(v[i+1,j]+v[i,j])                                      - (u[i,j]+u[i,j-1])*(v[i+1,j-1]+v[i,j-1]))
  
  
  				uˢ[i,j] = u[i,j] + Δt*(-c₁-c₂+c₃)
  
			end
		end
		
		# vˢ calculation
		for j=2:nny-1
			for i=2:nnx
  
  				c₃ = (1.0/Re)*((v[i+1,j]-2.0*v[i,j]+v[i-1,j])/(Δx²)  				                      + (v[i,j+1]-2.0*v[i,j]+v[i,j-1])/(Δy²) )
                   
  				c₂ = (0.25/Δy)*((v[i,j+1]+v[i,j])^2 - (v[i,j]+v[i,j-1])^2)
  
  				c₁ = (0.25/Δx)*((u[i,j+1]+u[i,j])*(v[i+1,j]+v[i,j])                                      - (u[i-1,j+1]+u[i-1,j])*(v[i,j]+v[i-1,j]))
  
  
  				vˢ[i,j] = v[i,j] + Δt*(-c₁-c₂+c₃)
  
			end
		end
		
		# Boundary conditions for velocity correction
		for i=1:nnx
			uˢ[i,1] =-uˢ[i,2]
			uˢ[i,nny+1] = 2.0-uˢ[i,nny]
		end
		for i=1:nny
			vˢ[1,i] =-vˢ[2,i]
			vˢ[nnx+1,i] =-vˢ[nnx,i]
		end
		
		# Boundary conditions for pressure
		for i=1:nex+2
			p[i,1] = p[i,2]
			p[i,ney+2] = p[i,ney+1]
		end
		for i=1:ney+2
			p[1,i] = p[2,i]
			p[nex+2,i] = p[nex+1,i]
		end
		
		# Poisson Solver
		ω = 1.0
		for ii = 1:iimax
			for j=2:ney+1
				for i=2:nex+1
  
  					f = ((uˢ[i,j]-uˢ[i-1,j])/Δx + (vˢ[i,j]-vˢ[i,j-1])/Δy)/Δt
                   
  					r = (f - (p[i+1,j]-2.0*p[i,j]+p[i-1,j])/(Δx²)                                           -      (p[i,j+1]-2.0*p[i,j]+p[i,j-1])/(Δy²) )
  
  
  					p[i,j] = p[i,j] + ω*r/a
  
				end
			end
		end
		
		# Correct velocity calculation
		for j=2:nny
			for i=2:nnx-1
  				u[i,j] = uˢ[i,j] - Δt*(p[i+1,j]-p[i,j])/Δx
			end
		end

		for j=2:nny-1
			for i=2:nnx
  				v[i,j] = vˢ[i,j] - Δt*(p[i,j+1]-p[i,j])/Δy
			end
		end
		
		#check residual
		rᵤ = 0.0
		rᵥ = 0.0
		rₚ = 0.0
		for j=2:nny
			for i=2:nnx
  				rᵤ = rᵤ + (u[i,j] - u₀[i,j])^2
  				rᵥ = rᵥ + (v[i,j] - v₀[i,j])^2
  				rₚ = rₚ + (p[i,j] - p₀[i,j])^2
			end
		end
		rᵤ = sqrt(rᵤ/(nex*ney))
		rᵥ = sqrt(rᵥ/(nex*ney))
		rₚ = sqrt(rₚ/(nex*ney))
		
		rᵤₜ[tt] = rᵤ
		rᵥₜ[tt] = rᵥ
		rₚₜ[tt] = rₚ
		ttend = tt
		@show tt rᵤ rᵥ rₚ ttend rᵤₜ[tt]	
		if (rᵤ<=ϵ && rᵥ<=ϵ && rₚ<=ϵ)
			break
		end
		
	end
	ttend = 1:ttend
	
	return u,v,p,nnx,ttend,rᵤₜ,rᵥₜ,rₚₜ
end

# ╔═╡ 9be09b70-6f2a-11eb-0f64-018b8c98e73f
begin
	u100c,v100c,p100c = MAC(100)
	contourf(u100c',levels=50,c= :redsblues, clims=(-1, 1),linewidth = 0,title = "Re=100 - u Contour")
end

# ╔═╡ 6c9f9cc0-6e54-11eb-380f-9f8f2d11343c
begin
	Re100data = [1	1	0.9766	0.84123	0.9688	0.78871	0.9609	0.73722	0.9531	0.68717	0.8516	0.23151	0.7344	0.00332	0.6172	-0.13641	0.5000 -0.20581	0.4531 -0.21090	0.2813	-0.15662	0.1719	-0.1015	0.1016	-0.06434	0.0703	-0.04775	0.0625	-0.04192	0.0547	-0.03717	0	0]
	y100data = [Re100data[i] for i=1:2:34]
	u100data = [Re100data[i] for i=2:2:34]

	u100_0_02,v100_0_02,p100_0_02,nnx100_0_02 = MAC(100,0.02)
	y100_0_02 = -0.01:0.02:1+0.01
	nnx100_0_02m = floor(Int32,(nnx100_0_02-1)/2)
		
	u100_0_01,v100_0_01,p100_0_01,nnx100_0_01 = MAC(100,0.01)
	y100_0_01 = -0.005:0.01:1+0.005
	nnx100_0_01m = floor(Int32,(nnx100_0_01-1)/2)
	
	u100_0_005,v100_0_005,p100_0_005,nnx100_0_005 = MAC(100,0.005)
	y100_0_005 = -0.0025:0.005:1+0.0025
	nnx100_0_005m = floor(Int32,(nnx100_0_005-1)/2)
	
	plot([u100_0_02[nnx100_0_02m,:], u100_0_01[nnx100_0_01m,:], u100_0_005[nnx100_0_005m,:], u100data],[y100_0_02,y100_0_01,y100_0_005,y100data],
				title = "Re=100 - Center u Velocity",
				label = ["MAC h=0.02" "MAC h=0.01" "MAC h=0.005" "Ghia"],
				legend = :bottomright,
				#ylims = (0,0.04),
				xlabel = "U",
				ylabel = "y",
				linewidth = 2,
				linecolor = [:black :black :black :transparent],
				markershape = :auto,
				markercolor = [:red :red :red :green],
				)
end

# ╔═╡ 07e46510-6f1f-11eb-21a0-f50e1e1e67b8
begin
	
	 x,x,x,x,tt100_0_02,ru100_0_02,rv100_0_02,rp100_0_02 = MAC(100,0.02)
	
	plot([tt100_0_02],[log.(ru100_0_02[1:tt100_0_02[end]]),log.(rv100_0_02[1:tt100_0_02[end]]),log.(rp100_0_02[1:tt100_0_02[end]])],
				title = "Re=100 - h=0.02 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				xlims = (0,1500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
end

# ╔═╡ a71f5160-6f26-11eb-1ffb-15467ef1e901
begin
	
	q,q,q,q,tt100_0_01,ru100_0_01,rv100_0_01,rp100_0_01 = MAC(100,0.01)
	
	plot([tt100_0_01],[log.(ru100_0_01[1:tt100_0_01[end]]),log.(rv100_0_01[1:tt100_0_01[end]]),log.(rp100_0_01[1:tt100_0_01[end]])],
				title = "Re=100 - h=0.01 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				#xlims = (0,500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
end

# ╔═╡ 13382bd0-6f2a-11eb-392a-6fc93389bda6
begin
	
	w,w,w,w,tt100_0_005,ru100_0_005,rv100_0_005,rp100_0_005 = MAC(100,0.005)
	
	plot([tt100_0_005],[log.(ru100_0_005[1:tt100_0_005[end]]),log.(rv100_0_005[1:tt100_0_005[end]]),log.(rp100_0_005[1:tt100_0_005[end]])],
				title = "Re=100 - h=0.005 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				#xlims = (0,500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
end

# ╔═╡ 4c463b00-6f2b-11eb-2840-7384b467a27b
begin
	u400c,v400c,p400c = MAC(400)
	contourf(u400c',levels=50,c= :redsblues, clims=(-1, 1),linewidth = 0,title = "Re=400 - u Contour")
end

# ╔═╡ bfad9c20-6e7f-11eb-2fc9-49f1df970d96
begin
	Re400data = [1	1	0.9766	0.75837	0.9688	0.68439	0.9609	0.61756	0.9531	0.55892	0.8516	0.29093	0.7344	0.16256	0.6172	0.02135	0.5000 -0.11477	0.4531 -0.17119	0.2813	-0.32726	0.1719	-0.24299	0.1016	-0.14612	0.0703	-0.10338	0.0625	-0.09266	0.0547	-0.08186	0	0
]
	y400data = [Re400data[i] for i=1:2:34]
	u400data = [Re400data[i] for i=2:2:34]

	u400_0_02,v400_0_02,p400_0_02,nnx400_0_02 = MAC(400,0.02)
	y400_0_02 = -0.01:0.02:1+0.01
	nnx400_0_02m = floor(Int32,(nnx400_0_02-1)/2)
		
	u400_0_01,v400_0_01,p400_0_01,nnx400_0_01 = MAC(400,0.01)
	y400_0_01 = -0.005:0.01:1+0.005
	nnx400_0_01m = floor(Int32,(nnx400_0_01-1)/2)
	
	u400_0_005,v400_0_005,p400_0_005,nnx400_0_005 = MAC(400,0.005)
	y400_0_005 = -0.0025:0.005:1+0.0025
	nnx400_0_005m = floor(Int32,(nnx400_0_005-1)/2)
	
	plot([u400_0_02[nnx400_0_02m,:], u400_0_01[nnx400_0_01m,:], u400_0_005[nnx400_0_005m,:], u400data],[y400_0_02,y400_0_01,y400_0_005,y400data],
				title = "Re=400 - Center Velocity",
				label = ["MAC h=0.02" "MAC h=0.01" "MAC h=0.005" "Ghia"],
				legend = :bottomright,
				#ylims = (0,0.04),
				xlabel = "U",
				ylabel = "y",
				linewidth = 2,
				linecolor = [:black :black :black :transparent],
				markershape = :auto,
				markercolor = [:red :red :red :green],
				)
end

# ╔═╡ 26043960-6f2b-11eb-2fdd-33f7de7fc9e3
begin
	
	 t,t,t,t,tt400_0_02,ru400_0_02,rv400_0_02,rp400_0_02 = MAC(400,0.02)
	
	plot([tt400_0_02],[log.(ru400_0_02[1:tt400_0_02[end]]),log.(rv400_0_02[1:tt400_0_02[end]]),log.(rp400_0_02[1:tt400_0_02[end]])],
				title = "Re=400 - h=0.02 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				xlims = (0,1500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
end

# ╔═╡ 03cf5b90-6f2b-11eb-0dd8-7f2061f9d9d7
begin
	
	r,r,r,r,tt400_0_01,ru400_0_01,rv400_0_01,rp400_0_01 = MAC(400,0.01)
	
	plot([tt400_0_01],[log.(ru400_0_01[1:tt400_0_01[end]]),log.(rv400_0_01[1:tt400_0_01[end]]),log.(rp400_0_01[1:tt400_0_01[end]])],
				title = "Re=400 - h=0.01 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				#xlims = (0,500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
end

# ╔═╡ e680ea90-6f2a-11eb-1ca2-89a396add38a
begin
	
	e,e,e,e,tt400_0_005,ru400_0_005,rv400_0_005,rp400_0_005 = MAC(400,0.005)
	
	plot([tt400_0_005],[log.(ru400_0_005[1:tt400_0_005[end]]),log.(rv400_0_005[1:tt400_0_005[end]]),log.(rp400_0_005[1:tt400_0_005[end]])],
				title = "Re=400 - h=0.005 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				#xlims = (0,500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
end

# ╔═╡ 57ed04c2-6f2b-11eb-20c8-1725c926dd2a
begin
	u1000c,v1000c,p1000c = MAC(1000)
	contourf(u1000c',levels=50,c= :redsblues, clims=(-1, 1),linewidth = 0,title = "Re=1000 - u Contour")
end

# ╔═╡ c8ece650-6e80-11eb-2882-7fed169d61ed
begin
	Re1000data = [1	1	0.9766	0.65928	0.9688	0.57492	0.9609	0.51117	0.9531	0.46604	0.8516	0.33304	0.7344	0.18719	0.6172	0.05702	0.5	-0.0608	0.4531	-0.10648	0.2813	-0.27805	0.1719	-0.38289	0.1016	-0.2973	0.0703	-0.2222	0.0625	-0.20196	0.0547	-0.18109	0	0]
	y1000data = [Re1000data[i] for i=1:2:34]
	u1000data = [Re1000data[i] for i=2:2:34]
	
	u1000_0_02,v1000_0_02,p1000_0_02,nnx1000_0_02 = MAC(1000,0.02)
	y1000_0_02 = -0.01:0.02:1+0.01
	nnx1000_0_02m = floor(Int32,(nnx1000_0_02-1)/2)
		
	u1000_0_01,v1000_0_01,p1000_0_01,nnx1000_0_01 = MAC(1000,0.01)
	y1000_0_01 = -0.005:0.01:1+0.005
	nnx1000_0_01m = floor(Int32,(nnx1000_0_01-1)/2)
	
	u1000_0_005,v1000_0_005,p1000_0_005,nnx1000_0_005 = MAC(1000,0.005)
	y1000_0_005 = -0.0025:0.005:1+0.0025
	nnx1000_0_005m = floor(Int32,(nnx1000_0_005-1)/2)
	
	plot([u1000_0_02[nnx1000_0_02m,:], u1000_0_01[nnx1000_0_01m,:], u1000_0_005[nnx1000_0_005m,:], u1000data],[y1000_0_02,y1000_0_01,y1000_0_005,y1000data],
				title = "Re=1000 - Center Velocity",
				label = ["MAC h=0.02" "MAC h=0.01" "MAC h=0.005" "Ghia"],
				legend = :bottomright,
				#ylims = (0,0.04),
				xlabel = "U",
				ylabel = "y",
				linewidth = 2,
				linecolor = [:black :black :black :transparent],
				markershape = :auto,
				markercolor = [:red :red :red :green],
				)
end

# ╔═╡ 94b982e0-6f29-11eb-2f68-839a418ad1f8
begin
	
	 c,c,c,c,tt1000_0_02,ru1000_0_02,rv1000_0_02,rp1000_0_02 = MAC(1000,0.02)
	
	plot([tt1000_0_02],[log.(ru1000_0_02[1:tt1000_0_02[end]]),log.(rv1000_0_02[1:tt1000_0_02[end]]),log.(rp1000_0_02[1:tt1000_0_02[end]])],
				title = "Re=1000 - h=0.02 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				#xlims = (0,1500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
	
end

# ╔═╡ b5637060-6f28-11eb-0850-5dff9a108c4b
begin
	
	 b,b,b,b,tt1000_0_01,ru1000_0_01,rv1000_0_01,rp1000_0_01 = MAC(1000,0.01)
	
	plot([tt1000_0_01],[log.(ru1000_0_01[1:tt1000_0_01[end]]),log.(rv1000_0_01[1:tt1000_0_01[end]]),log.(rp1000_0_01[1:tt1000_0_01[end]])],
				title = "Re=1000 - h=0.01 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				#xlims = (0,1500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
	
end

# ╔═╡ 2e0fecf0-6f29-11eb-1d52-71a1c891ed5a
begin
	
	 a,a,a,a,tt1000_0_005,ru1000_0_005,rv1000_0_005,rp1000_0_005 = MAC(1000,0.005)
	
	plot([tt1000_0_005],[log.(ru1000_0_005[1:tt1000_0_005[end]]),log.(rv1000_0_005[1:tt1000_0_005[end]]),log.(rp1000_0_005[1:tt1000_0_005[end]])],
				title = "Re=1000 - h=0.005 - Residual",
				label = ["u residual" "v residual" "p residual"],
				legend = :topright,
				#xlims = (0,1500),
				ylims = (-15,-2.5),
				xlabel = "Iteration",
				ylabel = "log(Residual)",
				linewidth = 2,
				linecolor = [:auto],
				#markershape = :auto,
				#markercolor = [:red],
				)
	
end

# ╔═╡ Cell order:
# ╟─57aa3170-6e51-11eb-2670-9b77065bc4da
# ╟─ad3d9470-6e82-11eb-3072-43decfc60ef9
# ╟─484fcd20-6f23-11eb-0dc8-5d49530983a1
# ╟─6bc808d0-6e51-11eb-0794-b9585413d41b
# ╟─9be09b70-6f2a-11eb-0f64-018b8c98e73f
# ╟─6c9f9cc0-6e54-11eb-380f-9f8f2d11343c
# ╟─07e46510-6f1f-11eb-21a0-f50e1e1e67b8
# ╟─a71f5160-6f26-11eb-1ffb-15467ef1e901
# ╟─13382bd0-6f2a-11eb-392a-6fc93389bda6
# ╟─4c463b00-6f2b-11eb-2840-7384b467a27b
# ╟─bfad9c20-6e7f-11eb-2fc9-49f1df970d96
# ╟─26043960-6f2b-11eb-2fdd-33f7de7fc9e3
# ╟─03cf5b90-6f2b-11eb-0dd8-7f2061f9d9d7
# ╟─e680ea90-6f2a-11eb-1ca2-89a396add38a
# ╟─57ed04c2-6f2b-11eb-20c8-1725c926dd2a
# ╟─c8ece650-6e80-11eb-2882-7fed169d61ed
# ╟─94b982e0-6f29-11eb-2f68-839a418ad1f8
# ╟─b5637060-6f28-11eb-0850-5dff9a108c4b
# ╟─2e0fecf0-6f29-11eb-1d52-71a1c891ed5a

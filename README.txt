In this assignment I have done all the normal part. And the ambient temperature in my system is 0(since only the difference between ambient temperature and the grid temperature matters.) And you can change the coefficient of bouyancy force and vorticity in related function.
I have attached two video, one with low vorticity and another one with high vorticity.


Questions for you to answer (10 points)

1.Describe the difference between Lagrangian and Eulerian viewpoints for simulation.  Why is the approach used in this assignment called Semi-Lagragian?
Lagrangian approach treats the continuum as a particle system. Each point in the fluid or solid is labeled as a separate
particle, with a position and a velocity .But in Eulerian approach,instead of tracking each particle, we instead look at fixed points in space and see how the fluid quantities (such as density, velocity, temperature, etc.) measured at those points change in time. The approach used in this assignment called Semi-Lagragian is because we use both methods to simulate the fluids. 

2.Smoke and water can both be simulated with fluids.  Briefly explain how they are similar and how they are different. 
They are similar because they both has Incompressibility feature. The difference is that smoke and water has different pressure in the system and density. Since water has large thickness than smoke, so the Semi-Lagragian system in water should handle larger pressure than in smoke. So smoke simlation should consider buoyancy force and Vorticity.While water simulation should consider Signed Distance.

3.List one advantage and one disadvantage to simulating our fluid on a grid.  Describe two other techniques for simulating fluids and the advantages and disadvantages of each.
Advantage:It¡¯s easier analytically to work with the spatial derivatives like the pressure gradient and viscosity in the
Eulerian viewpoint(in a grid)
Disadvantage:Adding Eulerian approach(the grid approach) in our simulation method is more complicated than pure particle system approach.
Other techniques:
1>.Smoothed-particle hydrodynamics (SPH):advantage: It is a mesh-free Lagrangian method (where the coordinates move with the fluid), and the resolution of the method can easily be adjusted with respect to variables such as the density.
					disadvantage:SPH has zero dissipation.The disadvantage of this in the scheme is that,where dissipative terms are required physically, they must be explicitly added. In particular this is the case for shock-capturing, since shocks lead to a physical increase in entropy. Shocks and other kinds of discontinuities are not adequately captured by the Hamiltonian formulation of SPH since in employing the Euler-Lagrange equations we have assumed that the quantities in the Lagrangian (i.e. thermal energy and velocity) are differentiable, implying that discontinuities in those variables need special treatment.

2>.Procedural Water: advantage:An advantage of procedural animation is its controllability, an important feature in
games. 
		    disadvantage:The main disadvantage of simulating water procedurally is the difficulty of getting the interaction of the water with the boundary or with immersed bodies right

4.How do level set surfaces differ from mesh surfaces? List one advantage and disadvantage of using a mesh to represent a surface.  List one advantage and disadvantage to using a Level Set Surface. 
A¡°level set¡±is simply an implicit surface which is define by function like {x | Q(x) = 0};while mesh surface is defined by faces and vertices. 
mesh surface:advantage:a more easy and intuitive to define a surface and easier to render.
	    disadvantage: when doing intersection tests based on mesh surface, we should loop around all the surfaces and related vertices which is very expensive. 



level set: advantage: The chief advantage  is that there is no connectivity to worry about, making operations such as topology changes (e.g. two water drops merging together, or a sheet breaking apart into droplets) trivial and it  can very easily model smooth surfaces.
          disadvantage: harder to render a level set surface.


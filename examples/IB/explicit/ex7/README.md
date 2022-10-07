This is a tutorial example used to track the IB coordinates of a single post which undergoes prescribed motion. In this example, the IB coordinates of the post are prescribed to lie and move along the level set of a cone who's radius and length is set by the post_deflection_radius and post_length either in the input file, or using the default values as specified in the example.cpp file. After running the executable associated with this example, a "IBCoordinates.txt" file is generate whose columns correspond to the relevant Lagrangian index as specified in the .vertex file and the Eulerian x,y,z coordinates. The rows identify a given IB point at a particular timestep. For example, if the post generated contains 30 IB points and one would like know the location of the Mth IB point at the Nth timestep, one would consult the 30*(N)+M row of the IBCoordinates.txt file. Note that the initial configuration corresponds to the zeroth timestep.

The MATLAB script "generate_posts.m" is used to generate the IB points of the posts located at the center of the square domain. 
The remaining problem parameters are set in an input file (currently `input3d.NFINEST=128` is the only one in the repository).
Notice that a number of problem parameters are currently "hard coded" in this script.
If the size of the domain (`L`, `ASPECT_RATIO_X`, `ASPECT_RATIO_Y`, `ASPECT_RATIO_Z`) or the domain grid resolution (`NFINEST`) are changed in either file, they should be changed in both files.
E.g. to increase the grid resolution, update `NFINEST` in both files and regenerate the posts via the "generate_posts.m" file. 
Notice also that `NFINEST` is a derived quantity --- it is controlled by the number of grid cells on the coarsest level of the AMR grid (`N`), the number of grid levels (`MAX_LEVELS`) and the refinement ratio between levels (`REF_RATIO`).

The post properties are set only in the MATLAB script, and can be modified independently of the rest of the problem settings (except that if the post stiffness is increased, it may be necessary to decrease the value of the maximum time step size, `DT_MAX`, in the `input3d` file).

The fluid properties are set only in the input file.
Currently they correspond to water.



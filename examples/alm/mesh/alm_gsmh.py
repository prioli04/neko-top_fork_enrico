import gmsh
import sys

import numpy as np

# Follows from geometric progression sum formula: Sn = a1 * (q^n - 1) / (q - 1)
# With Sn = L, a1 = dx_refine * q, q = growth_factor
# Substitute in the sum formula and solve for n
def n_progression(L, dx_refine, growth_factor):
    # L -> line length
    # dx_refine -> mesh size in the refinement region
    # growth_factor -> geometric progression ratio for growing the mesh size
    a1 = dx_refine * growth_factor
    n = np.log(1.0 + (growth_factor - 1.0) * L / a1) / np.log(growth_factor)
    return int(np.ceil(n))

# Compute the vector of normalized height (see the extrude function documentation for a definition of normalized height)
def extrude_heights_progression(L, h_refine, dx_refine, growth_factor):
    # L -> extrude length
    # h_refine -> normalized height at the end of the refinement region
    # dx_refine -> mesh size in the refinement region
    # growth_factor -> geometric progression ratio for growing the mesh size

    # Ensure L is positive
    L = np.abs(L)

    # Initialize normalized height vector and loop variables
    h_vec = [] 
    dx, h = dx_refine, h_refine

    # Loop until h gets to 1 (end of the extrude)
    while h < 1.0:
        dx *= growth_factor
        h += dx / L
        h_vec.append(h)

    # Scale h_vec such that the last element is 1
    h_vec = np.array(h_vec) / h_vec[-1] 
    h_vec = h_vec[1:] if h_vec[0] <= 0.5 else h_vec
    return h_vec.tolist()

# Parameters
geo_file_name = "alm.geo"
b = 2.0 # Wing span
domain_spans = [-2.5, 5.0, -2.5, 2.5, -2.5, 2.5] # Domain coordinates in number of spans [xmin, xmax, ymin, ymax, zmin, zmax]
al_refine_spans = [-0.25, 0.25, -1.25, 1.25, -0.25, 0.25] # Refinement around the actuator in number of spans [xmin, xmax, ymin, ymax, zmin, zmax]
al_refine_size_spans = 1.0 / 4.0 # Cell size inside the refinement region, in number of spans
growth_factor = 1.1 # Geometric progression ratio for expanding cell sizes

# Dimensionalize lengths
domain_coords = np.array(domain_spans) * b
al_refine_coords = np.array(al_refine_spans) * b
al_refine_size = al_refine_size_spans * b

# Initialize gmsh model
gmsh.initialize(sys.argv)
gmsh.model.add("alm")
## Create model in the 2D plane at y_min of the refinement box ##
refine_point_ids, refine_line_ids = [], []
left_point_ids, left_line_ids = [], []
right_point_ids, right_line_ids = [], []
top_point_ids, top_line_ids = [], []
bottom_point_ids, bottom_line_ids = [], []
corner_point_ids, corner_line_ids = [], []

top_left_line_ids, top_right_line_ids, bottom_left_line_ids, bottom_right_line_ids = [], [], [], []

# Refinement box points
refine_point_ids.append(gmsh.model.geo.add_point(al_refine_coords[0], 0.0, al_refine_coords[4], 1.0))
refine_point_ids.append(gmsh.model.geo.add_point(al_refine_coords[0], 0.0, al_refine_coords[5], 1.0))
refine_point_ids.append(gmsh.model.geo.add_point(al_refine_coords[1], 0.0, al_refine_coords[4], 1.0))
refine_point_ids.append(gmsh.model.geo.add_point(al_refine_coords[1], 0.0, al_refine_coords[5], 1.0))

# Refinement box lines
refine_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[0], refine_point_ids[2]))
refine_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[2], refine_point_ids[3]))
refine_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[3], refine_point_ids[1]))
refine_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[1], refine_point_ids[0]))

# Refinement box face
refine_curve_loop_id = gmsh.model.geo.add_curve_loop(refine_line_ids[0:4])
refine_face_id = gmsh.model.geo.add_plane_surface([refine_curve_loop_id])

# Left points
left_point_ids.append(gmsh.model.geo.add_point(domain_coords[0], 0.0, al_refine_coords[4], 1.0))
left_point_ids.append(gmsh.model.geo.add_point(domain_coords[0], 0.0, al_refine_coords[5], 1.0))

# Right points 
right_point_ids.append(gmsh.model.geo.add_point(domain_coords[1], 0.0, al_refine_coords[4], 1.0))
right_point_ids.append(gmsh.model.geo.add_point(domain_coords[1], 0.0, al_refine_coords[5], 1.0))

# Top points
top_point_ids.append(gmsh.model.geo.add_point(al_refine_coords[0], 0.0, domain_coords[5], 1.0))
top_point_ids.append(gmsh.model.geo.add_point(al_refine_coords[1], 0.0, domain_coords[5], 1.0))

# Bottom points
bottom_point_ids.append(gmsh.model.geo.add_point(al_refine_coords[0], 0.0, domain_coords[4], 1.0))
bottom_point_ids.append(gmsh.model.geo.add_point(al_refine_coords[1], 0.0, domain_coords[4], 1.0))

# Corner points
corner_point_ids.append(gmsh.model.geo.add_point(domain_coords[0], 0.0, domain_coords[4], 1.0))
corner_point_ids.append(gmsh.model.geo.add_point(domain_coords[0], 0.0, domain_coords[5], 1.0))
corner_point_ids.append(gmsh.model.geo.add_point(domain_coords[1], 0.0, domain_coords[4], 1.0))
corner_point_ids.append(gmsh.model.geo.add_point(domain_coords[1], 0.0, domain_coords[5], 1.0))

# Left box lines
left_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[0], left_point_ids[0]))
left_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[1], left_point_ids[1]))
left_line_ids.append(gmsh.model.geo.add_line(left_point_ids[1], left_point_ids[0]))

# Right box lines
right_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[2], right_point_ids[0]))
right_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[3], right_point_ids[1]))
right_line_ids.append(gmsh.model.geo.add_line(right_point_ids[0], right_point_ids[1]))

# Top box lines
top_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[1], top_point_ids[0]))
top_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[3], top_point_ids[1]))
top_line_ids.append(gmsh.model.geo.add_line(top_point_ids[0], top_point_ids[1]))

# Bottom box lines
bottom_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[0], bottom_point_ids[0]))
bottom_line_ids.append(gmsh.model.geo.add_line(refine_point_ids[2], bottom_point_ids[1]))
bottom_line_ids.append(gmsh.model.geo.add_line(bottom_point_ids[0], bottom_point_ids[1]))

# Top left box lines
top_left_line_ids.append(gmsh.model.geo.add_line(left_point_ids[1], corner_point_ids[1]))
top_left_line_ids.append(gmsh.model.geo.add_line(top_point_ids[0], corner_point_ids[1]))

# Top right box lines
top_right_line_ids.append(gmsh.model.geo.add_line(right_point_ids[1], corner_point_ids[3]))
top_right_line_ids.append(gmsh.model.geo.add_line(top_point_ids[1], corner_point_ids[3]))

# Bottom left box lines
bottom_left_line_ids.append(gmsh.model.geo.add_line(left_point_ids[0], corner_point_ids[0]))
bottom_left_line_ids.append(gmsh.model.geo.add_line(bottom_point_ids[0], corner_point_ids[0]))

# Bottom right box lines
bottom_right_line_ids.append(gmsh.model.geo.add_line(right_point_ids[0], corner_point_ids[2]))
bottom_right_line_ids.append(gmsh.model.geo.add_line(bottom_point_ids[1], corner_point_ids[2]))

# Top left box face
top_left_curve_loop_id = gmsh.model.geo.add_curve_loop([top_line_ids[0], top_left_line_ids[1], -top_left_line_ids[0], -left_line_ids[1]])
top_left_face_id = gmsh.model.geo.add_plane_surface([top_left_curve_loop_id])

# Left box face
left_curve_loop_id = gmsh.model.geo.add_curve_loop([-refine_line_ids[3], -left_line_ids[0], left_line_ids[2], left_line_ids[1]])
left_face_id = gmsh.model.geo.add_plane_surface([left_curve_loop_id])

# Bottom left box face
bottom_left_curve_loop_id = gmsh.model.geo.add_curve_loop([-bottom_line_ids[0], left_line_ids[0], bottom_left_line_ids[0], -bottom_left_line_ids[1]])
bottom_left_face_id = gmsh.model.geo.add_plane_surface([bottom_left_curve_loop_id])

# Bottom box face
bottom_curve_loop_id = gmsh.model.geo.add_curve_loop([-refine_line_ids[0], bottom_line_ids[0], bottom_line_ids[2], -bottom_line_ids[1]])
bottom_face_id = gmsh.model.geo.add_plane_surface([bottom_curve_loop_id])

# Bottom right box face
bottom_right_curve_loop_id = gmsh.model.geo.add_curve_loop([bottom_line_ids[1], bottom_right_line_ids[1], -bottom_right_line_ids[0], -right_line_ids[0]])
bottom_right_face_id = gmsh.model.geo.add_plane_surface([bottom_right_curve_loop_id])

# Right box face
right_curve_loop_id = gmsh.model.geo.add_curve_loop([-refine_line_ids[1], right_line_ids[0], right_line_ids[2], -right_line_ids[1]])
right_face_id = gmsh.model.geo.add_plane_surface([right_curve_loop_id])

# Top right box face
top_right_curve_loop_id = gmsh.model.geo.add_curve_loop([-top_line_ids[1], right_line_ids[1], top_right_line_ids[0], -top_right_line_ids[1]])
top_right_face_id = gmsh.model.geo.add_plane_surface([top_right_curve_loop_id])

# Top box face
top_curve_loop_id = gmsh.model.geo.add_curve_loop([-refine_line_ids[2], top_line_ids[1], -top_line_ids[2], -top_line_ids[0]])
top_face_id = gmsh.model.geo.add_plane_surface([top_curve_loop_id])
gmsh.model.geo.synchronize()

# Cell sizes
# Ceil ensures sizes can't be larger than requested. The +1 is due to Gmsh requiring number of points, not cells
dx_refine, dz_refine = al_refine_coords[1] - al_refine_coords[0], al_refine_coords[5] - al_refine_coords[4]
nx_refine = int(np.ceil(dx_refine / al_refine_size) + 1) + 1
nz_refine = int(np.ceil(dz_refine / al_refine_size) + 1)

# Refinement size in the x direction
gmsh.model.geo.mesh.set_transfinite_curve(refine_line_ids[0], nx_refine)
gmsh.model.geo.mesh.set_transfinite_curve(refine_line_ids[2], nx_refine)
gmsh.model.geo.mesh.set_transfinite_curve(top_line_ids[2], nx_refine)
gmsh.model.geo.mesh.set_transfinite_curve(bottom_line_ids[2], nx_refine)

# Refinement size in the z direction
gmsh.model.geo.mesh.set_transfinite_curve(refine_line_ids[1], nz_refine)
gmsh.model.geo.mesh.set_transfinite_curve(refine_line_ids[3], nz_refine)
gmsh.model.geo.mesh.set_transfinite_curve(left_line_ids[2], nz_refine)
gmsh.model.geo.mesh.set_transfinite_curve(right_line_ids[2], nz_refine)

# Transition size in the x direction to the left
nx_left = n_progression(np.abs(domain_coords[0] - al_refine_coords[0]), al_refine_size, growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(left_line_ids[0], nx_left, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(left_line_ids[1], nx_left, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(top_left_line_ids[1], nx_left, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(bottom_left_line_ids[1], nx_left, coef=growth_factor)

# Transition size in the x direction to the right
nx_right = n_progression(np.abs(domain_coords[1] - al_refine_coords[1]), al_refine_size, growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(right_line_ids[0], nx_right, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(right_line_ids[1], nx_right, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(top_right_line_ids[1], nx_right, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(bottom_right_line_ids[1], nx_right, coef=growth_factor)

# Transition size in the z direction to the bottom
nz_bottom = n_progression(np.abs(domain_coords[4] - al_refine_coords[4]), al_refine_size, growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(bottom_line_ids[0], nz_bottom, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(bottom_line_ids[1], nz_bottom, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(bottom_left_line_ids[0], nz_bottom, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(bottom_right_line_ids[0], nz_bottom, coef=growth_factor)

[gmsh.model.geo.mesh.set_transfinite_surface(id + 1) for id in range(gmsh.model.geo.get_max_tag(2))]
[gmsh.model.geo.mesh.set_recombine(2, id + 1) for id in range(gmsh.model.geo.get_max_tag(2))]

# Transition size in the z direction to the top
nz_top = n_progression(np.abs(domain_coords[5] - al_refine_coords[5]), al_refine_size, growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(top_line_ids[0], nz_top, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(top_line_ids[1], nz_top, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(top_left_line_ids[0], nz_top, coef=growth_factor)
gmsh.model.geo.mesh.set_transfinite_curve(top_right_line_ids[0], nz_top, coef=growth_factor)

[gmsh.model.geo.mesh.set_transfinite_surface(id + 1) for id in range(gmsh.model.geo.get_max_tag(2))]
[gmsh.model.geo.mesh.set_recombine(2, id + 1) for id in range(gmsh.model.geo.get_max_tag(2))]
gmsh.model.geo.synchronize()

# Generate 2D mesh (must happen before the extrude)
gmsh.model.mesh.generate(2)

# Extrude 2D plane to y+
extrude_tags = [(2, i+1) for i in range(gmsh.model.geo.get_max_tag(2))] # Include all surfaces
dy_plus = domain_coords[3]

height_refine_plus = al_refine_coords[3] / dy_plus
height_farfield_plus = extrude_heights_progression(dy_plus, height_refine_plus, al_refine_size, growth_factor)

ny_refine_plus = int(np.ceil(abs(al_refine_coords[3]) / al_refine_size)) + 1
ny_farfield_plus = [1] * len(height_farfield_plus) # Grow farfield layer by layer

gmsh.model.geo.extrude(extrude_tags, 0.0, dy_plus, 0.0, numElements=[ny_refine_plus] + ny_farfield_plus, heights=[height_refine_plus] + height_farfield_plus, recombine=True)

# Extrude 2D plane to y-
dy_minus = domain_coords[2]

height_refine_minus = al_refine_coords[2] / dy_minus
height_farfield_minus = extrude_heights_progression(dy_minus, height_refine_minus, al_refine_size, growth_factor)

ny_refine_minus = int(np.ceil(abs(al_refine_coords[2]) / al_refine_size)) + 1
ny_farfield_minus = [1] * len(height_farfield_minus) # Grow farfield layer by layer

gmsh.model.geo.extrude(extrude_tags, 0.0, dy_minus, 0.0, numElements=[ny_refine_minus] + ny_farfield_minus, heights=[height_refine_minus] + height_farfield_minus, recombine=True)

gmsh.model.geo.synchronize()

# Set physical names
inlet_ids = [63, 85, 107, 261, 283, 305, # Front Face
             310, 332, 354, 376, 398, 244, 420, 266, 288, # Left Face
             200, 178, 156, 134, 46, 222, 112, 68, 90, # Right Face
             59, 217, 199, 397, 415, 257, # Top Face
             111, 129, 147, 345, 327, 309 # Bottom Face
             ]

outlet_ids = [195, 173, 151, 349, 371, 393] # Back Face
volume_ids = list(range(1, gmsh.model.geo.getMaxTag(3) + 1)) # All volumes

gmsh.model.geo.add_physical_group(2, inlet_ids, tag=1, name="inlet") # Add inlet
gmsh.model.geo.add_physical_group(2, outlet_ids, tag=2, name="outlet") # Add outlet
gmsh.model.geo.add_physical_group(3, volume_ids, tag=3, name="fluid") # Add fluid
gmsh.model.geo.synchronize()

gen_mesh = True

if gen_mesh:
    # Generate 3D mesh
    gmsh.option.set_number("Mesh.MshFileVersion", 2.2)
    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.set_order(2)
    gmsh.write("alm.msh")

else:
    gmsh.write("alm.geo_unrolled")

gmsh.finalize()
 

# importing gmsh
import gmsh
import os.path

# this is the class to generate mesh file for square lamina
class SquareLamina:

	def __init__(self, h, file_name, save_path):
		self.h = h                  # mesh size
		self.file_name = file_name  # name of the mesh file
		self.save_path = save_path  # path to save the mesh file

	def msh(self):
		# Gmsh must be initialized:
		gmsh.initialize()

		# Next we add a new model name
		gmsh.model.add("{}".format(self.file_name))

		# mesh size
		lc = self.h

		# points
		gmsh.model.geo.addPoint(0.0, 0.0, 0, lc/5.0, 1)
		gmsh.model.geo.addPoint(0.0, 1.0, 0, lc, 2)
		gmsh.model.geo.addPoint(1.0, 1.0, 0, lc, 3)
		gmsh.model.geo.addPoint(1.0, 0.0, 0, lc/5.0, 4)

		# lines
		gmsh.model.geo.addLine(1, 2, 1)
		gmsh.model.geo.addLine(2, 3, 2)
		gmsh.model.geo.addLine(3, 4, 3)
		gmsh.model.geo.addLine(4, 1, 4)

		# curve loop
		gmsh.model.geo.addCurveLoop([4, 1, 2, 3], 1)

		# surface
		gmsh.model.geo.addPlaneSurface([1], 1)

		# synchronize
		gmsh.model.geo.synchronize()

		# Boundary conditions
		#gmsh.model.addPhysicalGroup(1, [3], name="Fixed Edge")
		#gmsh.model.addPhysicalGroup(1, [5], name="Displaced Edge")

		# generate a 2D mesh
		gmsh.model.mesh.generate(2)

		# save it to disk
		os.chdir(self.save_path)                     # change directory
		gmsh.option.setNumber("Mesh.SaveAll", 1)
		gmsh.write("{}.msh".format(self.file_name))  # write mesh file in the changed directory

		# just to check the geometry created (write .geo file)
		gmsh.write("{}.geo_unrolled".format(self.file_name))

		# end of the Gmsh Python API:
		gmsh.finalize()

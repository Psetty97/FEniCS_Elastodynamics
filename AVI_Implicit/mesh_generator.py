# import the test cases here
#from shear_test import ShearTest
#from peel_test import PeelTest
from square_lamina import SquareLamina

# calling the shear test class
#              |mesh size
#                       |name of the mesh file
#                                                       | path to save the mesh file
#p1 = ShearTest(0.002, "shear_test_0.002_eps_0.015", '/home/fenics/shared/meshes')
#p1.msh()

#p1 = PeelTest(0.0075, "peel_test_0.0075_eps_adaptive", '/home/fenics/shared/meshes')
#p1.msh()

#p1=SquareLamina(0.1, "square_lamina_0.1", '/home/FAUAD/ag33ineh/AVI_Implicit/meshes')
p1=SquareLamina(1.5, "square_lamina_smaller", '/home/FAUAD/ag33ineh/AVI_Implicit/meshes')
p1.msh()



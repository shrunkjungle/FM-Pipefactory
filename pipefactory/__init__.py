import importlib.util

from .UtilityFunctions import rot_vec, get_orthogonal_inplane, find_section, rotate_point_about_point, rotate_triad, rotate_vector, cylinder_geodesic_distance

from .FEM.Mesh import Element, Node
from .FEM.shapefunctions import LinearHex, LinearQuad
from .FEM.intpoints import IntPoint, GaussQuadrature, GaussQuadrature2D
from .FEM.LinearElastic import LinearElasticity
from .FEM.Loads import CollarLoad
from .FEM.FEM import FEM

from .Defects.Defects import Hole, Dimple, Radial_Slit, Weld, RadialCrack
from .Defects.Refinement import AxialRefinement, Ramp

from .Pipe import Pipe
from .BuildTools import PipeParam, PartitionROM

from .UnitBrick import UnitBrick, TwoBricks

def is_extra_installed(package_name):
    return importlib.util.find_spec(package_name) is not None


PARALLEL_ENABLED = is_extra_installed("mpi4py")

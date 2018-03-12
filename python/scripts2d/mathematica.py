from scripts2d.utils import dist_from_code
import sys

d = dist_from_code(sys.argv[1])
print(d.mathematica())

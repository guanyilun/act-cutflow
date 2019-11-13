"""This script aims to implement the cuts in the moby2 framework
without removing the pre-selection. """

from todloop import TODLoop

from routines.tod import LoadTOD, CheckTODLength
from routines.cuts import CutMCE

# initialize loop
loop = TODLoop()
# add list of tod
loop.add_tod_list("./inputs/mr3_pa2_s16_test.txt")
# load tod
loop.add_routine(LoadTOD())
loop.add_routine(CheckTODLength(fmin=10, min_periods=10))  # placeholder values
# cut mce error
loop.add_routine(CutMCe())

# run loop
loop.run(0,20)

"""This script aims to implement the cuts in the moby2 framework
without removing the pre-selection. """

from todloop import TODLoop

from routines.tod import LoadTOD, CheckTODLength, TransformTOD
from routines.cuts import CutMCE, CutSources, CutPlanets, \
    RemoveSyncPickup, CutPartial, FindPathologies

# initialize loop
loop = TODLoop()
# add list of tod
loop.add_tod_list("./inputs/mr3_pa2_s16_test.txt")
loop.add_routine(LoadTOD())
loop.add_routine(CheckTODLength(fmin=10, min_periods=10))  # placeholder values
loop.add_routine(CutMCE())
loop.add_routine(CutSources())
loop.add_routine(CutPlanets())
loop.add_routine(RemoveSyncPickup())
loop.add_routine(CutPartial())
loop.add_routine(TransformTOD())
loop.add_routine(FindPathologies())

# run loop
loop.run(0,20)

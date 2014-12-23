import run_CUBIT_lib as rcl 
import cubit
reload(rcl)

cubit.init([''])

configfile = '/home/psakicki/THESE/MODEL_GWADA/PYTHON_FUNCTIONS/configfiles/run_CUBIT.cfg'
#cubit.cmd("graphics window create")

## GENERATE CONFIG FILE FOR EACH COMBINATION
#rcl.multi_config_generator4cubit(configfile)
## RUN MULTI 
#rcl.geom_cubit_multi(configfile)

# ============== OR ==============

## RUN SOLO
rcl.geom_cubit_solo('/home/psakicki/THESE/MODEL_GWADA/PYTHON_FUNCTIONS/configfiles/run_cubit_solo_exemple.cfg')

## RUN SOLO
#rcl.geom_cubit_solo('/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/output_CUBIT_mk10/meshdir_12/mesh_12.cfg')

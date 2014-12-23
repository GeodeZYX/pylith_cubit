import run_pylith_lib as rpl 
reload(rpl)
import numpy as np

configfile = '/home/psakicki/THESE/MODEL_GWADA/PYTHON_FUNCTIONS/configfiles/run_pylith.cfg'

# Generate the configfiles / spatialdbs files for each combination 
rpl.main_experience_generator(configfile)
# Run in pylith the combinations
rpl.main_run_pylith_multi(configfile)

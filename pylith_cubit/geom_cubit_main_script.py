import geom_cubit_lib as gcl 
import os
import glob
reload(gcl)

#cfpath='/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/output_CFs'
#cfwildcard='*.cf'
#
#os.chdir(cfpath)
#for cf in sorted(glob.glob(cfwildcard)):
#    print cf
#    gcl.geom_cubit_main(cfpath + '/' + cf)

# QUI MARCHE TOUJOURS
gcl.geom_cubit_main('/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/configfile_3a.cf')
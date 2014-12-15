import ConfigParser
import numpy as np
import itertools
import sys

cf = ConfigParser.ConfigParser()
if sys.argv[0] != '':
#    cf.read(sys.argv[0])
    cf.read('/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/bigconfig.multi.cf')
else:
    cf.read('/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/bigconfig.multi.cf')

dic = {}
cleans = cf.sections()
cleans.remove('Input')
cleans.remove('Output')
cleans.remove('backstop_bool')

for s in cleans:
    if not cf.getboolean(s,'variable'):
        val = [cf.getfloat(s,'fixed')]
    else:
        truemin = np.min(cf.getfloat(s,'min'),cf.getfloat(s,'max'))
        truemax = np.max(cf.getfloat(s,'min'),cf.getfloat(s,'max'))
        val = list(np.arange(truemin,truemax+0.1, cf.getfloat(s,'delta')))
    dic[s] = val
 
pstr = 'Parameters'
ostr = 'Output'
istr = 'Input'
                    
for i,V in enumerate(list(itertools.product(*dic.values()))):
    outcf = ConfigParser.SafeConfigParser() 
    outcf.add_section('Input')
    outcf.add_section('Output')
    outcf.add_section('Parameters')   
    # setting Inputs in the new configfile
    outcf.set(istr,'bathyfile',cf.get(istr,'bathyfile'))
    outcf.set(istr,'slabfile',cf.get(istr,'slabfile'))
    # setting Outputs in the new configfile    
    outcf.set(ostr,'output_path',cf.get(ostr,'output_path_meshs'))
    outcf.set(ostr,'output_prefix',cf.get(ostr,'output_prefix'))
    # setting parameters in the new configfile 
    outcf.set(pstr,'backstop_bool','value')
    for j,v in enumerate(list(V)):
        outcf.set(pstr,dic.keys()[j],str(v))
        outcf.set(ostr,'output_suffix',str(i))

        outfile = open(cf.get(ostr,'output_path_configfiles')  + '/' + cf.get(ostr,'output_prefix') + str(i) + '.cf','w')
        outcf.write(outfile)
        outfile.close()
        
                        
                        
                    

    

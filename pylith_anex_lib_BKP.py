import StringIO
from collections import OrderedDict
import ConfigParser
import itertools
import matplotlib.pyplot as plt
import numpy as np
import os
import netCDF4

# ================== README ==================
#A sdb dico have 3 keys : 
#        'name' : STRING, ex  :'x' 'y' 'z' 'density' ...
#        'unit' : STRING, ex : 'kg/m**3'  'm/s'  'm/s' ... 
#        'values' : LIST of the values 
#        
#        ONE value type per dico (but severals values in the values list)
#
# a val dico have n keys ('x','y','density'...) + a entity_name key :
#        param1 : [val1,val2 ...]
#        param2 : [val1,val2 ...]
#        ....
#        paramn : [val1,val2 ...]
#         +++       
#        entity_name : the name of the entity (thank you captain obvious)
#        ONE block/interface per dico
#        SEVERALS value type per dico
#
# a combi dico contains all information about ONE combinaition,
# the keys are like in the configfile '<name of the block/interface>.<name of the parameter>
#
# the discrimination of specific cases is in the configdic2valdic

def open_exodus(exodusin):
    if type(exodusin) is str:
        exodusout = netCDF4.Dataset(filenameExodus, 'a')
    else:
        exodusout = exodusin
    return exodusout
    
def blocksname_list(exodus):
    exodus = open_exodus(exodus)
    nbblocks = len(exodus.dimensions['num_el_blk'])
    blocksnamelist = []
    for i in range(nbblocks):
        blocksnamelist.append(''.join(exodus.variables['eb_names'][i,:]))
    return blocksnamelist

def nodesetname_list(exodus):
    exodus = open_exodus(exodus)
    nbnodesets = len(exodus.dimensions['num_node_sets'])
    nodesetnamelist = []
    for i in range(nbnodesets):
        nodesetnamelist.append(''.join(exodus.variables['ns_names'][i,:]))
    return nodesetnamelist

def get_XY_of_nodeset(exodus,nsdesc):
    # from a EXODUS file imported with the NETCDF4 module
    # and a nodeset descriptor (name or ID)
    # retrun X & Y of the nodes, CUBIT's ID of the nodes and 
    # 'local' indices of the nodes in a pythonnic convention (start at 0)
    exodus = open_exodus(exodus)
    X = exodus.variables['coordx'][:]
    Y = exodus.variables['coordy'][:]
    mapp1 = exodus.variables['node_num_map'][:]
    if type(nsdesc) is int:
        nsid = nsdesc
    elif type(nsdesc) is str:
        nsid = nodesetname_list(exodus).index(nsdesc) + 1      
    nodesidlist = exodus.variables['node_ns'+str(nsid)][:] - 1
    return X[nodesidlist],Y[nodesidlist],mapp1[nodesidlist],nodesidlist

def sort_XY_nodeset(x,y):
    xy = np.vstack((x,y))
    xy2 = np.sort(xy,axis=1)
    xy3 = np.fliplr(xy2)
    return xy3[0,:] , xy3[1,:]
    
def dist(xa,ya,xb,yb):
    return np.sqrt((xa-xb)**2 + (ya-yb)**2)
    
def cumul_dist(X,Y,i1,i2):
    finaldist = 0
    for i in range(i1,i2):
        finaldist = finaldist + dist(X[i],Y[i],X[i+1],Y[i+1])
    return finaldist
               
def unit_giver(namein):
    if namein == 'x' or namein == 'y' or namein == 'z':
            return 'km'
    elif namein == 'density':
        return 'kg/m**3'
    elif namein == 'vs':
        return 'm/s'
    elif namein == 'vp':
        return 'm/2'
    elif namein == 'viscosity':
        return 'Pa*s'
    elif namein == 'left-lateral-slip':
        return 'cm/year'
    elif namein == 'fault-opening':
        return 'cm/year'
    else:
        return None

def make_combi_4_lockzone(Xin,Yin,ymax,lmin,inunit='km',outunit='km'):
    # ymax : the ultimate depth where is it supposed to be a locking
    # lmin : minimal length of a locking zone in the same combi
    #        AND b/w 2 successive top point in 2 different combis
    if inunit == 'km':
        Xin = Xin * 1000
        Yin = Yin * 1000
    if outunit == 'km':
        k = 1000
    else:
        k = 1

    X,Y = sort_XY_nodeset(Xin,Yin)
    imax = np.min(np.flatnonzero(Y < ymax))
    rangei = np.arange(imax)
    xybound_list = []
    i_lasttoppt = 0
    Ncombi_tot = 0
    Ncombi_valid = 0
    for c in list(itertools.combinations(rangei,2)):
        Ncombi_tot += 1
        # excluding locking zone with a too small 
        if cumul_dist(X,Y,c[0],c[1]) < lmin:
            continue
        # excluding locking zone with a top point too close of the last one
        if i_lasttoppt != c[0] and cumul_dist(X,Y,i_lasttoppt,c[0]) < lmin:
            continue
        xybound_list.append(((X[c[0]]/k,Y[c[0]]/k),(X[c[1]]/k,Y[c[1]]/k)))
        i_lasttoppt = c[0]
        Ncombi_valid += 1
        
    print 'end of locking combinations, valid / total :' ,Ncombi_valid, '/', Ncombi_tot
    return xybound_list
     
def check_exodus_configdic(exodus,configfile):
    exodus = open_exodus(exodus)
    blkns_in_exodus = blocksname_list(exodus) + nodesetname_list(exodus)
    
    cf = ConfigParser.ConfigParser()
    cf.read(configfile)
    sections = [s for s in cf.sections() if s.find('.') >= 0 ]
    blkns_in_cf = list(set([ s.split('.')[0] for s in sections ]))
    
    missing_e = []
    for e in blkns_in_exodus:
        if not e in blkns_in_cf:
            missing_e.append(e)
    
    if missing_e != []:
        print " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
        print "WARN : missing parameters for blocs /nodesets in "
        print  configfile
        print " "
        for e in missing_e:
            print e 
        print " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
        
            
    return None
    

def product_on_dic(dicin):
    varNames = sorted(dicin)
    return [dict(zip(varNames, prod)) for prod in itertools.product(*(dicin[varName] for varName in varNames))]


def configfile2combidic_list(exodus,path,faultlist=['fault_top']):
    exodus = open_exodus(exodus)
    cf = ConfigParser.ConfigParser()
    cf.read(path)
    configdic = OrderedDict()
    # Removing section with no variables intervals => with no '.'
    sections = cf.sections()
    for s in sections:
        if s.find('.') == -1 or not cf.has_option(s,'variable'):
            continue
        if not cf.getboolean(s,'variable'):
            val = [cf.getfloat(s,'fixed')]
        elif cf.getfloat(s,'delta') == 0:
            val = [cf.getfloat(s,'fixed')]
            print "WARN : !!! variable parameter but delta = 0 !!!"
        else:
            truemin = np.min([cf.getfloat(s,'min'),cf.getfloat(s,'max')])
            truemax = np.max([cf.getfloat(s,'min'),cf.getfloat(s,'max')])
            val = list(np.arange(truemin,truemax +0.1, cf.getfloat(s,'delta')))
        configdic[s] = val
        
    # MAKING THE COMBINATIONS OF THE LOCKING ZONE
    for faultname in  faultlist:
        print faultname
        ymax = cf.getfloat(faultname,'maxi_locked_depth')
        lmin = cf.getfloat(faultname,'mini_locked_length')
        X,Y,_,_ = get_XY_of_nodeset(exodus,faultname)
        configdic[faultname + '.combiXY'] = make_combi_4_lockzone(X/1000,Y/1000,ymax,lmin,inunit='km',outunit='km')
        # SPECIFIC
        for badfield in [ k for k in configdic.keys() if 'maxi_locked_depth' in k ]:
            print badfield , 'aaa'
            configdic.pop(badfield,None)
        # CARTESIAN PRODUCT OF A CONFIG DIC => MANY COMBIDIC
        combidiclis = product_on_dic(configdic)
        print "total combinations : ", len(combidiclis)
        
    return combidiclis,configdic
    
def combidic2valdic_list(combidic):
    entity_name_lis = []
    # getting the name of all entity (block/nodesets) in the combidico
    # (for the moment the keys are <entity>.<parameter>)
    for k,v in combidic.iteritems():
        entity_name_lis.append(k.split('.',1)[0])
    list(set(entity_name_lis))
    # iterating over all the entities
    valdic_lis = []
    for e in entity_name_lis:
        valdic = OrderedDict()
        valdic['entity_name'] = e
        for k,v in combidic.iteritems():
            if e in k: # here is a parameter of the entity
                p = k.split('.',1)[1]
                if p == 'combiXY': #SPECIFIC CASES
                    delta = 0.1
                    valdic['x'] = [0,v[0][0],v[0][0],v[1][0],v[1][0],0]
                    valdic['y'] = [99, v[0][1] ,v[0][1]-delta ,v[1][1],v[1][1]-delta,-999]
                elif p == 'fault-opening':
                    valdic[p] = [v] * 6
                elif p == 'left-lateral-slip':
                    valdic[p] = [v,v,0,0,v,v]     
                else: #GENERIC CASE
                    valdic['x'] = [0]
                    valdic['y'] = [0]
                    if type(v) is int or float:
                        valdic[p] = [v]
                    else:
                        valdic[p] = v

        valdic_lis.append(valdic)  
    return valdic_lis
 
def valdic2sdbdic_list(valdicin):
    outlistofdic = []
    for k,v in valdicin.iteritems():
        if k == 'entity_name':
            continue
        dic = OrderedDict()
        dic['name'] = k
        dic['unit'] = unit_giver(dic['name'])
        if type(v) is list:
            dic['values'] = v
        else:
            dic['values'] = [v]
                    
        outlistofdic.append(dic)
    
    entity_name = valdicin['entity_name']
    return outlistofdic,entity_name

   
def write_spatialdb2(listofdic,outpath,outname):
    ''' This function take a LIST of sdb dico as argument '''
    strm =  StringIO.StringIO()
    numlocs = 1 # ???????? ADEFINIR !!!!!!!!
    datadim = 0 # ???????? ADEFINIR !!!!!!!!
    if sum([d['name'] == 'z' for d in listofdic]) != 0:
        spacedim = 3
    else:
        spacedim = 2 
    # finding the X value for the unit :
        unitcoef = 1
        for d in listofdic:
            if d['name'] == 'x':
                if d['unit'] == 'km':
                    unitcoef = 1000
    # write usual header    
    strm.write('#SPATIAL.ascii 1\n')
    strm.write('SimpleDB {\n')
    strm.write('  num-values = {}\n'.format(len(listofdic) - spacedim))
    namestr = ''
    unitstr = ''
    for d in listofdic:
        # we must exclude the x,y,z name in the header
        if d['name'] in ('x','y','z'):
            continue
        namestr = namestr + ' ' + d['name']
        namestr = namestr.replace('_','-')
        unitstr = unitstr + ' ' + d['unit']  
    strm.write('  num-names = {}\n'.format(namestr))
    strm.write('  num-units = {}\n'.format(unitstr))
    
    strm.write('  num-locs = {}\n'.format(numlocs))
    strm.write('  data-dim = {}\n'.format(datadim))
    
    strm.write('  space-dim = {}\n'.format(spacedim))
 
    strm.write('  cs-data = cartesian {\n')
    strm.write('    to-meters = {}\n'.format(unitcoef))
    strm.write('    space-dim = {}\n'.format(spacedim))
    strm.write('  } \n')
    strm.write('} \n')
    
    stdformat = '{:8.1f} '
    
    for line in zip(*[ d['values'] for d in listofdic ]):
        tempstr = stdformat * len(listofdic) + '\n'
        strm.write(tempstr.format(*line))
            
    fil = open(outpath + '/' + outname + '.spatialdb','w')
    
    fil.write(strm.getvalue())
    fil.close()
            
    return strm

    
        
def main_config2multiexp(exodus,configfile,outpath,outdirprefix='exp'):
    exodus = open_exodus(exodus)
    check_exodus_configdic(exodus,configfile)       
    combidiclis,confdic = configfile2combidic_list(exodus,configfile) 
    i = 0
    for combidic in combidiclis:
        i = i+1 
        valdiclis = combidic2valdic_list(combidic)
        for valdic in valdiclis:
            sdbdiclis,name = valdic2sdbdic_list(valdic)
            outdirectory = outpath + outdirprefix + '_' + str(i)
            if not os.path.exists(outdirectory):
                os.makedirs(outdirectory)
            write_spatialdb2(sdbdiclis, outdirectory ,name)
      
# ZONE DE TEST
prepath=''
#prepath='/run/user/1000/gvfs/sftp:host=calipso.univ-lr.fr,user=psakicki'

filenameExodus = prepath + "/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/output_MESH/mesh__solo1.exo"
cfpath = prepath + '/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/pylith_bigrun.cf'
outpath = '/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/test_sdb/'

main_config2multiexp(filenameExodus,cfpath,outpath)  
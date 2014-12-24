import StringIO
from progressbar import ProgressBar
from collections import OrderedDict
import ConfigParser
import itertools
import numpy as np
import os
import netCDF4
import sys
import matplotlib.pyplot as plt
import time
import copy
import multiprocessing
import glob
import subprocess
import re
import csv

# ================== README ==================
#A sdb dico have 3 keys : 
#        'name' : STRING, ex  :'x' 'y' 'z' 'density' ...
#        'unit' : STRING, ex : 'kg/m**3'  'm/s'  'm/s' ... 
#        'values' : LIST of the values 
#        
#        ONE value type per dico (but severals values in the values list)
#
# a val dico have n keys ('x','y','density'...) + a entity_name key + entity_type key:
#        param1 : [val1,val2 ...]
#        param2 : [val1,val2 ...]
#        ....
#        paramn : [val1,val2 ...]
#         +++       
#        entity_name : the name of the entity (thank you captain obvious)
#        (all the annex keys must be removed in valdic2sdbdic_list)
#
#        ONE block/interface per dico
#        SEVERALS value type per dico
#
# a combi dico contains all information about ONE combinaition,
# the keys are like in the configfile '<name of the block/interface>.<name of the parameter>
#
# the discrimination of specific cases is in the configdic2valdic

#  ________   ______  _____  _    _  _____   ______                _   _                  
# |  ____\ \ / / __ \|  __ \| |  | |/ ____| |  ____|              | | (_)                 
# | |__   \ V / |  | | |  | | |  | | (___   | |__ _   _ _ __   ___| |_ _  ___  _ __  ___  
# |  __|   > <| |  | | |  | | |  | |\___ \  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __| 
# | |____ / . \ |__| | |__| | |__| |____) | | |  | |_| | | | | (__| |_| | (_) | | | \__ \ 
# |______/_/ \_\____/|_____/ \____/|_____/  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/ 
#                                                                                         
#                                                                                         

def open_exodus(exodusin):
#    print 'open a exodus'
    if type(exodusin) is str:
        exodusout = netCDF4.Dataset(exodusin, 'a')
    else:
        exodusout = exodusin
    return exodusout
    
def make_exodusdic(exodusin):
    if type(exodusin) is dict:
        return exodusin # assuming the dico is a exodusdic
    exodus = open_exodus(exodusin)
    exodus_dic = {}
    exodus_dic['title'] = str(exodus.title)
    if type(exodusin) is str:
        exodus_dic['filename'] = os.path.basename(exodusin)
    else:
        exodus_dic['filename'] = os.path.basename(str(exodus.title).split('(')[1].split(')')[0])
    exodus_dic['blocks'] = blocksname_list(exodus)
    exodus_dic['nodesets'] = nodesetname_list(exodus)
    exodus_dic['X'] = exodus.variables['coordx'][:]
    exodus_dic['Y'] = exodus.variables['coordy'][:]
    exodus_dic['map'] = exodus.variables['node_num_map'][:]
    
    for i in range(len(exodus_dic['nodesets'])):
        exodus_dic['node_ns'+str(i+1)] = exodus.variables['node_ns'+str(i+1)][:]
    
    return exodus_dic

def open_configfile(configfilein):
    if type(configfilein) is str:    
#        print 'open a cf'
        cfout = ConfigParser.ConfigParser()
        cfout.optionxform = str 
        cfout.read(configfilein)
    else:
        cfout  = configfilein
    return cfout
    
def blocksname_list(exodus):
    if type(exodus) is dict:
        blocksnamelist = exodus['blocks']
    elif type(exodus) is str or netCDF4.Dataset :      
        exodus = open_exodus(exodus)
        nbblocks = len(exodus.dimensions['num_el_blk'])
        blocksnamelist = []
        for i in range(nbblocks):
            blocksnamelist.append(''.join([x for x in exodus.variables['eb_names'][i,:] if type(x) is np.string_]))
    return blocksnamelist

def nodesetname_list(exodus):
    if type(exodus) is dict:
        nodesetnamelist = exodus['nodesets']   
    elif type(exodus) is str or netCDF4.Dataset :      
        exodus = open_exodus(exodus)
        nbnodesets = len(exodus.dimensions['num_node_sets'])
        nodesetnamelist = []
        for i in range(nbnodesets):
            nodesetnamelist.append(''.join([x for x in exodus.variables['ns_names'][i,:] if type(x) is np.string_]))
    return nodesetnamelist
    
    
def check_exodus_configfile(exodus,configfile):
#    exodus = open_exodus(exodus)
    exodus = make_exodusdic(exodus)
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
        print ''
        print " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
        print "WARN : missing parameters for blocs /nodesets in "
        print  configfile
        print " "
        for e in missing_e:
            print e 
        print " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! "
        print " "
        
    for e in blkns_in_cf:
        if not e in blkns_in_exodus:
            raise Exception(e + 'is not present in the exodus file !!!')
            
    return None

def get_XY_of_nodeset_from_NETCDFexodus(exodus,nsdesc):
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
    
def get_XY_of_nodeset_from_exodusdic(exodusdic,nsdesc):
    X = exodusdic['X']
    Y = exodusdic['Y']
    mapp1 = exodusdic['map']
    
    if type(nsdesc) is int:
        nsid = nsdesc
    elif type(nsdesc) is str:
        nsid = exodusdic['nodesets'].index(nsdesc) + 1 
        
    nodesidlist = exodusdic['node_ns'+str(nsid)] - 1
    return X[nodesidlist],Y[nodesidlist],mapp1[nodesidlist],nodesidlist


def get_XY_of_nodeset(exodusin,nsdesc):
    if type(exodusin) is dict:
        x,y,a,b = get_XY_of_nodeset_from_exodusdic(exodusin,nsdesc)
    elif type(exodusin) is str or netCDF4.Dataset :
        x,y,a,b = get_XY_of_nodeset_from_NETCDFexodus(exodusin,nsdesc)
    else:
        raise Exception('exodus is not a NETCDF or a dict' )
    return x,y,a,b


#   _____                 _ _   ______                _   _                 
#  / ____|               | | | |  ____|              | | (_)                
# | (___  _ __ ___   __ _| | | | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
#  \___ \| '_ ` _ \ / _` | | | |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
#  ____) | | | | | | (_| | | | | |  | |_| | | | | (__| |_| | (_) | | | \__ \
# |_____/|_| |_| |_|\__,_|_|_| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                                                                           
      
def sort_XY_nodeset(x,y):    
    xy = np.vstack((x,y)).T
    xy2 = (xy[xy[:,1].argsort()])
    return np.flipud(xy2[:,0]), np.flipud(xy2[:,1].T)


def check_if_combiXY(combidic,entity):
    for k,v in combidic.items():
        if entity in k:
            if 'combiXY' in k:
                return True , v
    return False , None

def product_on_dic(dicin):
    varNames = sorted(dicin)
    return [dict(zip(varNames, prod)) for prod in itertools.product(*(dicin[varName] for varName in varNames))]

   
def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)

def listofstr2str(listin):
    outstr = '['
    for s in listin:
        outstr  = outstr + s + ','
    outstr = outstr[:-1] + ']'
    return outstr

def flipper(X,Y):
    # flip X,Y array in order to get the deepest point at the end of the array
    if Y[-1] > Y[1]:
        X = np.flipud(X)
        Y = np.flipud(Y)
    return X,Y   
    
def dist(xa,ya,xb,yb):
    return np.sqrt((xa-xb)**2 + (ya-yb)**2)
    
def cumul_dist(X,Y,i1,i2):
    finaldist = 0
    for i in range(i1,i2):
        finaldist = finaldist + dist(X[i],Y[i],X[i+1],Y[i+1])
    return finaldist

#
#   _____                _     _   ______                _   _                 
#  / ____|              | |   (_) |  ____|              | | (_)                
# | |     ___  _ __ ___ | |__  _  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
# | |    / _ \| '_ ` _ \| '_ \| | |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |___| (_) | | | | | | |_) | | | |  | |_| | | | | (__| |_| | (_) | | | \__ \
#  \_____\___/|_| |_| |_|_.__/|_| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
#                                                                              
#                                                                              


def unit_giver(namein,representative_value = 1 , limit_m_km = 100000 ):
    if namein == 'x' or namein == 'y' or namein == 'z':
        if representative_value < limit_m_km:
            return 'km'
        else:
            return 'm'
    elif namein == 'density':
        return 'kg/m**3'
    elif namein == 'vs':
        return 'm/s'
    elif namein == 'vp':
        return 'm/s'
    elif namein == 'viscosity':
        return 'Pa*s'
    elif namein == 'left-lateral-slip':
        return 'cm/year'
    elif namein == 'fault-opening':
        return 'cm/year'
    else:
        return 'NoneUnit'
        
def make_combi_4_lockzone(Xin,Yin,ymax,lmin,stepmin,inunit='km',outunit='km'):
    # ymax : the ultimate depth where is it supposed to be a locking
    # lmin : minimal length of a locking zone in the same combi
    #        AND b/w 2 successive top point in 2 different combis
    if inunit == 'km':
        Xin = Xin * 1000
        Yin = Yin * 1000
        ymax = ymax * 1000
        lmin = lmin * 1000
        stepmin = stepmin * 1000
    if outunit == 'km':
        k = 1000
    else:
        k = 1
    X,Y = Xin,Yin
    # search if the first elt of X,Y is the top or the bottom of slab
    X,Y = flipper(X,Y)
    imax = np.min(np.flatnonzero(Y < ymax))
    print 'making locking combinations'
    print "ymax = " , ymax , " imax = " , imax
    rangei = np.arange(imax)
    xybound_list = []
    ijbound_list = []
    i_lasttoppt = 0
    Ncombi_tot = 0
    Ncombi_valid = 0
    for c in list(itertools.combinations(rangei,2)):
        Ncombi_tot += 1
        # excluding locking zone with a too small l
        if cumul_dist(X,Y,c[0],c[1]) < lmin:
            continue
        # excluding locking zone with a top point too close of the last one
        if i_lasttoppt != c[0] and cumul_dist(X,Y,i_lasttoppt,c[0]) < stepmin:
            continue
        xybound_list.append(((X[c[0]]/k,Y[c[0]]/k),(X[c[1]]/k,Y[c[1]]/k)))
        ijbound_list.append((c[0],c[1]))

        i_lasttoppt = c[0]
        Ncombi_valid += 1
        
    print 'end of locking combinations, valid / total :' ,Ncombi_valid, '/', Ncombi_tot
    return ijbound_list


def configfile2combidic_list(path,exodusin):
    cf = ConfigParser.ConfigParser()
    cf.read(path)
    sections = cf.sections()   
    exodus_param = cf.get('Parameters','exodus_name')
    
    exodus = open_exodus(exodusin)    
    check_exodus_configfile(exodus,path)
    exodusdic = make_exodusdic(exodusin)
    
    entity_name_lis = []
    for s in sections:
        if '.' in s:
            entity_name_lis.append(s.split('.',1)[0])
    
    entity_name_lis = list(set(entity_name_lis))
    
    paramdic = {}
    for e in entity_name_lis:
        if cf.has_section(e + '.parameters'):
            for param,val in cf.items(e + '.parameters'):
                try:
                    val = eval(val) 
                except NameError:
                    val = val
                if type(val) is list:
                    truemin = min(val[0:2])
                    truemax = max(val[0:2])
                    delta = val[2]
                    paramdic[e + '.' +  param] = list(np.arange(truemin,truemax,delta))           
                    if cf.has_option(e + '.general','spacing_mode'):  # try for the special case of a exponetial progression
                        if cf.get(e + '.general','spacing_mode') == 'exp':
                            paramdic[e + '.' +  param] = [ eval('1e+' + str(x)) for x in list(np.arange(truemin,truemax,delta))]           
                elif type(val) is str:
                    if val.strip()[0] == '-': # using a minus at the beginning of the string as marker for sign inversion
                        val2 = val.strip()[1:].replace(' ','')
                        paramdic[e + '.' +  param] = [val2 + '.' + param +'-'] # putting the minus at the end
                    else:
                        paramdic[e + '.' +  param] = [val.strip() + '.' + param]
                    #temporary value when the current val have to have another one
                    # change a few lines after
                elif type(val) is float or int : # this test MUST be the last one
                    paramdic[e + '.' +  param] = [val]   
                    
    # MAKING THE COMBINATIONS OF THE LOCKING ZONE
    faultlist = [s.split('.')[0] for s in cf.sections() if cf.has_option(s,'locked') and cf.getboolean(s,'locked') ]
    for faultname in  faultlist:
        print faultname
        ymax = cf.getfloat(faultname +'.general','maxi_locked_depth')
        lmin = cf.getfloat(faultname +'.general','mini_locked_length')
        stepmin = cf.getfloat(faultname +'.general','mini_locked_step')
        out_dist_units = cf.get('Parameters','out_dist_units')
        inp_dist_units = cf.get('Parameters','inp_dist_units')
        
        if inp_dist_units == 'km':
            kk = 1000
        else:
            kk = 1
            
        X,Y,_,_ = get_XY_of_nodeset(exodus,faultname)
        paramdic[faultname + '.combiXY'] = make_combi_4_lockzone(X/kk,Y/kk, \
        ymax,lmin,stepmin,inunit=inp_dist_units,outunit=out_dist_units)
        
    # SPECIFIC
#    for badfield in [ k for k in configdic.keys() if 'maxi_locked_depth' in k ]:
#        configdic.pop(badfield,None)
        
    # ===== CARTESIAN PRODUCT OF A CONFIG DIC => MANY COMBIDIC =====
    combidiclis = product_on_dic(paramdic)
    # ==============================================================
    print "total combinations for " + exodusdic['filename'] + ": ", len(combidiclis)
    if len(combidiclis) == 0:
        print "if total combinations == 0, check if imax == 0 too" 
        print "in this case ymax badly defined : problem of units (km <=> m) ?"
        
    # REPLACING TEMP STR VAL BY THE FINAL ONE
    for combidic in combidiclis:
        for k,v in combidic.items():
            if type(v) is str:
                try:
                    if v[-1] == '-': # discrimination in case of a minus inversion
                        combidic[k] = - combidic[v[:-1]]
                    else:
                        combidic[k] = combidic[v]     
                except:
                    print "replacement of the " + k + " key per " + v + " failed in configfile2combidic_list"               
    return combidiclis,paramdic             
    
def combidic2valdic2_list(combidic,configfile,exodusin):
    start = time.time()
    cf = open_configfile(configfile)
    exodus = make_exodusdic(exodusin)
    
    outunit = cf.get('Parameters','out_dist_units')

    if outunit == 'm':
        kk = 1
    elif outunit == 'km':
        kk = 1000
    else:
        kk = 1000
    entity_name_lis = []
    # getting the name of all entity (block/nodesets) in the combidico
    # (for the moment the keys are <entity>.<parameter>)
    # NB : Parameters without dot are naturally excluded
    for k,v in combidic.iteritems():
        entity_name_lis.append(k.split('.',1)[0])
    list(set(entity_name_lis))
    # iterating over all the entities
    valdic_lis = []
    for e in entity_name_lis:
        tempdico = {}
        for p,v in combidic.items() :
            if e in p:
                tempdico[p.split('.')[1]] = v
        valdic = OrderedDict()
        valdic['entity_name'] = e
        valdic['x'] = []
        valdic['y'] = []
        
        # TYPE BLOCK => only ONE line
        if cf.get(e + '.general','type') in ('bloc','block'):
            valdic['x'] = [0]
            valdic['y'] = [0]
            for p,v in tempdico.items():
                valdic[p] = [v]
        # TYPE INTERFACE
        elif cf.get(e + '.general','type') in ('interface'):
            # LOCKED CASE => speed = nominal speed and 0
            if cf.has_option(e + '.general','locked') and cf.getboolean(e + '.general','locked'):
                if not 'combiXY' in tempdico.keys():
                    raise Exception('ERR : ' + ' is a locked interface but have no locking boundary (combiXY)')
                X,Y,_,_ = get_XY_of_nodeset(exodus,e)
                X,Y = flipper(X,Y)
                X = X/kk
                Y = Y/kk
                valdic['x'] = list(X)
                valdic['y'] = list(Y)
                for p,v in tempdico.items():
                    if p == 'combiXY':
                        continue
                    elif p == 'left-lateral-slip':
                        valdic[p] = []
                        i = combidic[e + '.combiXY'][0]
                        j = combidic[e + '.combiXY'][1]

                        for n,(x,y) in enumerate(zip(X,Y)):
                            if i <= n <= j:
                                valdic[p].append(0)
                            else:
                                valdic[p].append(v)
                    else:
                        valdic[p] = [v] * len(valdic['x'])    
            # NOT LOCKED CASE => speed = nominal speed
            else:
                # case where the left lateral slip == 0
                if tempdico['left-lateral-slip'] == 0:
                    valdic['x'] = [0,0]
                    valdic['y'] = [9999,-999999]
                    for p,v in tempdico.items():
                        valdic[p] = [v] * len(valdic['x'])  
                else:   # case where the left lateral slip != 0
                    X,Y,_,_ = get_XY_of_nodeset(exodus,e)
                    X,Y = flipper(X,Y)
                    X = X/kk
                    Y = Y/kk
                    valdic['x'] = list(X)
                    valdic['y'] = list(Y)
                    for p,v in tempdico.items():
                        valdic[p] = [v] * len(valdic['x'])
                    

        valdic_lis.append(valdic) 
#    print "exec. time combidic2valdic :" , time.time() - start
    return valdic_lis
    
    
def valdic2sdbdic_list(valdicin):
    outlistofdic = []
    for k,v in valdicin.iteritems():
        if k == 'entity_name':
            continue
        if k == 'entity_type':
            continue
        dic = OrderedDict()
        if type(v) is list:
            dic['values'] = v
        else:
            dic['values'] = [v]      
        dic['name'] = k
        dic['unit'] = unit_giver(dic['name'],np.max(np.abs(dic['values'])))

        outlistofdic.append(dic)
    
    entity_name = valdicin['entity_name']
    return outlistofdic,entity_name

def write_spatialdb(valdicin,outpath):
    
    listofdic , _ = valdic2sdbdic_list(valdicin)
    outname =  valdicin['entity_name']
#    entitytype = valdicin['entity_name']
    
    strm =  StringIO.StringIO()
    # ----- numloc -----
    protonumloc = list(set([len(d['values']) for d in listofdic ]))    
    if len(protonumloc) != 1:
        print "WARN : write_spatialdb : differents numlocs"
        print "e.g. x : [1,2,3,4], vp : [10,20]"
        for d in listofdic:
            print d['name'],d['values'] 
    numlocs = protonumloc[0]
    # ----- datadim -----
    if len( valdicin['x'] ) == 1 and len( valdicin['y'] ) == 1:
        datadim = 0
    elif len(list(set(valdicin['x']))) != 1 and len(list(set(valdicin['x']))) != 1:
        datadim = 2
    else:
        datadim = 1
    # ----- spacedim ------
    if sum([d['name'] == 'z' for d in listofdic]) != 0:
        spacedim = 3
    else:
        spacedim = 2 
    # ----- unitcoef ------
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
    strm.write('  value-names = {}\n'.format(namestr))
    strm.write('  value-units = {}\n'.format(unitstr))
    
    strm.write('  num-locs = {}\n'.format(numlocs))
    strm.write('  data-dim = {}\n'.format(datadim))
    
    strm.write('  space-dim = {}\n'.format(spacedim))
 
    strm.write('  cs-data = cartesian {\n')
    strm.write('    to-meters = {}\n'.format(unitcoef))
    strm.write('    space-dim = {}\n'.format(spacedim))
    strm.write('  } \n')
    strm.write('} \n')
    
    stdformat = '{:12.3f} '
    
    for line in zip(*[ d['values'] for d in listofdic ]):
        tempstr = stdformat * len(listofdic) + '\n'
        strm.write(tempstr.format(*line))
            
    fil = open(outpath + '/' + outname + '.spatialdb','w')
    
    fil.write(strm.getvalue())
    fil.close()
            
    return strm



#  __  __       _         ______                _   _                 
# |  \/  |     (_)       |  ____|              | | (_)                
# | \  / | __ _ _ _ __   | |__ _   _ _ __   ___| |_ _  ___  _ __  ___ 
# | |\/| |/ _` | | '_ \  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# | |  | | (_| | | | | | | |  | |_| | | | | (__| |_| | (_) | | | \__ \
# |_|  |_|\__,_|_|_| |_| |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
                                                                                                                                     

      
def main_experience_generator(configfile):  
    cf_objt = open_configfile(configfile)    
    
    outpath           = cf_objt.get('Parameters','output_path') 
    local_cfg_generic = cf_objt.get('Parameters','local_generic')
    main_cfg_generic  = cf_objt.get('Parameters','main_generic')
    output_cfg_local  = cf_objt.get('Parameters','output_local')
    output_cfg_main   = cf_objt.get('Parameters','output_main')
    outdirprefix      = cf_objt.get('Parameters','output_dir_prefix')
 
    # Loading config files
    if not os.path.isfile(local_cfg_generic):
            raise Exception('generic config file ' + local_cfg_generic + " don't exist")
    if not os.path.isfile(main_cfg_generic):         
            raise Exception('generic config file ' + main_cfg_generic + " don't exist")

    local_cfg_generic_objt = open_configfile(local_cfg_generic)
    main_cfg_generic_objt = open_configfile(main_cfg_generic) 
    
    exodus_path = cf_objt.get('Parameters','exodus_name')
    
    if os.path.isfile(exodus_path):
        exodusdic = make_exodusdic(exodus_path)
        exodusdiclis = [exodusdic]
    elif os.path.isdir(exodus_path):
        exodusdiclis = []
        for protoexodus in glob.glob( exodus_path + '/' + '*.exo'):
            exodusdic = make_exodusdic(protoexodus)
            exodusdiclis.append(exodusdic)
    else:
        raise Exception('exodus file/folder ' + main_cfg_generic + " don't exist/is empty !")

    i = 0            
    for exodusdic in exodusdiclis:
        print ''
        print ' ==============================================================='
        print "processing combinations for " + exodusdic['filename']
        combidiclis,confdic = configfile2combidic_list(configfile,exodusdic) 
        # L is the number of combis, for folders creation
        L = str(len(str(len(combidiclis))))
        frmtstr = '{:0>' + L + 'd}'
        pbar = ProgressBar(maxval=len(combidiclis)+1)
        pbar.start()
        ipbar = 0
        for combidic in combidiclis:
            start = time.time()
            i = i+1 
            ipbar = ipbar + 1
    #        print i
            valdiclis = combidic2valdic2_list(combidic,cf_objt,exodusdic)
            for valdic in valdiclis:
                outdirectory = outpath + '/' + outdirprefix + '_' + frmtstr.format(i) + '/'
                if not os.path.exists(outdirectory):
                    os.makedirs(outdirectory)
                write_spatialdb(valdic, outdirectory)
            # writing the pylith config files for each combi
            write_configfile(cf_objt,copy.copy(local_cfg_generic_objt),copy.copy(main_cfg_generic_objt),outdirectory,output_cfg_local,output_cfg_main,exodusdic)
            # writing a summary file of the combi
            sumfile = open(outdirectory + '/' + outdirprefix + '_' + str(i) +'.sum', "w")
            sumfile.write('exodus_title , ' + exodusdic['title'] + '\n')            
            sumfile.write('exodus_filename , ' + exodusdic['filename'] + '\n')            
            w = csv.writer(sumfile)
            for key, val in OrderedDict(sorted(combidic.items())).items():
                w.writerow([key, val])
            sumfile.close()
            pbar.update(ipbar)
        pbar.finish()
    #        print 'execution combi :' , time.time() - start
    
    print ''
    print "FINISHING COMBINATION & CONFIGFILES WRITING PROCESS"
    print "total combinations for " + outdirprefix + ' : ' + str(i)
    
    return None

def main_run_pylith_multi(configfile):
    if not os.path.isfile(configfile) :
        raise Exception("ERR: CONFIG FILE DON'T EXIST !!!")
    cf_objt = ConfigParser.ConfigParser()
    cf_objt.read(configfile)
    maindir       = cf_objt.get('Parameters','output_path')  
    outdirprefix  = cf_objt.get('Parameters','output_dir_prefix')
    output_local  = cf_objt.get('Parameters','output_local')
    output_h5_name  = cf_objt.get('Parameters','output_h5_name')
    erase_existing  = cf_objt.getboolean('Parameters','erase_existing')
    
    # defining a valid path
    pylith_path= cf_objt.get('Parameters','pylith_path')  
    dic_env = os.environ
    dic_env['PATH'] = pylith_path+'/bin:'+ dic_env['PATH'] 
    dic_env['PYTHONPATH'] = pylith_path+'/lib/python2.7/site-packages:'+ pylith_path +'/lib64/python2.7/site-packages:'#+ dic_env['PYTHONPATH'] 
    dic_env['LD_LIBRARY_PATH'] = pylith_path+'lib:'+pylith_path+'lib64:' #+os.environ['LD_LIBRARY_PATH']
    
    dirlist = natural_sort([ d for d in os.listdir(maindir) if d.startswith(outdirprefix) ])
    cfwildcard = output_local
    
    errbig_file_path = maindir + "/err_summary.log"
    if os.path.isfile(errbig_file_path):
		os.remove(errbig_file_path)
    err_iter_list = []
    process_args_list = []
    pool = multiprocessing.Pool(processes=8)

    for d in dirlist:
        curdir = maindir + '/' + d
        os.chdir(curdir)
        if not os.path.exists('output'):
            os.makedirs('output')
        for cf in sorted(glob.glob(cfwildcard)):
            args = [cf,curdir,output_h5_name,pylith_path,dic_env,erase_existing,errbig_file_path]
            process_args_list.append(args)


        #  ===== MULTI PROCESS PART

    start_pool = time.time()            
#    results = [ pool.apply(process_unit,args = (a,)) for a in process_args_list ]
    results = pool.map(process_unit,process_args_list)
    results = sorted(results)
    
    errbig_file = open( errbig_file_path , "a")
    
#for f,std,err in results:
#if err != '':
#errbig_file.write('error log not empty for : '+ f +'\n')

    errbig_file.write('total exec. time : ' + str(time.time() - start_pool))
            

        #  ===== CLASSIC MONO PART


#            print '----------------------------------'
#            print "running : " + cf + ' in ' + curdir
#            if os.path.isfile('./output/' + output_h5_name + '.h5') and not erase_existing:
#                print './output/' + output_h5_name + '.h5 exists , skiping ...'
#                continue
#            start = time.time()
#            p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE ,env=dic_env)
#            command = pylith_path +"/bin/pylith " + curdir + '/' + cf
#            stdout,stderr = p.communicate( command )
#            std_file = open("std.log", "w")
#            err_file = open("err.log", "w")
#            std_file.write(stdout)
#            err_file.write(stderr)
#            std_file.close()
#            err_file.close()                
#            print "exec. time : " + str(time.time() - start) + ' s.'
#            if stderr != '':
#                print "err.log is not empty, must be checked !" 
#            err_iter_list.append(curdir)
#    
#    for e in err_iter_list:
#        errbig_file.write(e+'\n')
#    errbig_file.close()       
                
    return None
    
def process_unit(ain):
    cf,curdir,output_h5_name,pylith_path,dic_env,erase_existing,errbig_file_path = ain
    print '---------------- START ------------------'
    os.chdir(curdir)
    print "running : " + cf + ' in ' + curdir
    if os.path.isfile( curdir + '/output/' + output_h5_name + '.h5') and not erase_existing:
        print './output/' + output_h5_name + '.h5 exists , skiping ...'
        return curdir ,'',''
    print ''
    start = time.time()
    p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE , stdout=subprocess.PIPE , stderr=subprocess.PIPE ,env=dic_env)
    command = pylith_path +"/bin/pylith " + curdir + '/' + cf
    stdout,stderr = p.communicate( command )
    std_file = open("std.log", "w")
    err_file = open("err.log", "w")
    std_file.write(stdout)
    err_file.write(stderr)
    std_file.close()
    err_file.close()                
    print '----------------  END  ------------------'
    print "exec. time : " + str(time.time() - start) + ' s. for ' + os.path.basename(curdir)
    if stderr != '':
        print "err.log is not empty, must be checked !" 
        print stderr
        errbig_file = open(errbig_file_path,'a')
        errbig_file.write('error log not empty for : '+ curdir +'\n')
        errbig_file.close()        
    print ''
    return curdir , stdout , stderr
    

#                _ _                          __ _        __ _ _           
#               (_) |                        / _(_)      / _(_) |          
# __      ___ __ _| |_ ___    ___ ___  _ __ | |_ _  __ _| |_ _| | ___  ___ 
# \ \ /\ / / '__| | __/ _ \  / __/ _ \| '_ \|  _| |/ _` |  _| | |/ _ \/ __|
#  \ V  V /| |  | | ||  __/ | (_| (_) | | | | | | | (_| | | | | |  __/\__ \
#   \_/\_/ |_|  |_|\__\___|  \___\___/|_| |_|_| |_|\__, |_| |_|_|\___||___/
#                                                   __/ |                  
#                                                  |___/        


def remove_sections_beg(substring,config_obj):
    list_remove=[string for string in config_obj.sections() if string.startswith(substring)]
    [config_obj.remove_section(section) for section in list_remove]


def nodeset_type_from_configfile(configfilein,typein):
    cf = open_configfile(configfilein)
    outlist = []
    for s in cf.sections():
        if ( cf.has_option(s,'type') and cf.get(s,'type')==typein):
            outlist.append(s.split('.')[0])
    return outlist
    
def write_configfile(input_cfgfile,local_generic,main_generic,outpath,output_local,output_main,exodusdic):
    # This function must have many args (and not only the configfile)
    # because of speed of execution (openning many files with only the config file)
    # here the files are opened before
    
    ######## Checking and opening blocks
    
    ### Check if output path exist and open file for writing
    
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    if type(local_generic) is str:
        if not os.path.isfile(local_generic):
            raise Exception('generic config file ' + local_generic + " don't exist")

    if type(main_generic) is str:
        if not os.path.isfile(main_generic):         
            raise Exception('generic config file ' + main_generic + " don't exist")
        
    out_local=outpath + '/' + output_local # for the local.cfg file
    fic_out_local=open(out_local,'w')
    out_main=outpath  + '/' + output_main # for the main.cfg file
    fic_out_main=open(out_main,'w')
    
    ### Check and open generic files
    
    if type(local_generic) is str and type(main_generic) is str:
        if ( not os.path.isfile(local_generic) or not os.path.isfile(main_generic)):
            print 'WARNING: generic configuration files were not found, change path below and try again! \n' + local_generic
            raise Exception  
    
    ### Read generic files and Exodus
    
    config_local = open_configfile(local_generic)
    config_main = open_configfile(main_generic)
    cfginput = open_configfile(input_cfgfile)
    
    #exodus_name=cfginput.get('Parameters','exodus_name')
    exodus_name= cfginput.get('Parameters','exodus_name') + '/' + exodusdic['filename']
     
    ### Change settings
    
    list_of_blocks = blocksname_list(exodusdic)

    total_time=float(cfginput.get('Parameters','total_time'))
    dt=float(cfginput.get('Parameters','dt'))
    config_local.set('pylithapp.timedependent.formulation.time_step','total_time',"%.1f*year" %total_time)
    config_local.set('pylithapp.timedependent.formulation.time_step','dt',"%.1f*year" %dt)
    
    ### Get lists for bc and interfaces nodesets
        
    nodeset_bc=nodeset_type_from_configfile(cfginput,'bc')
    nodeset_interface=nodeset_type_from_configfile(cfginput,'interface')
    config_local.set('pylithapp.timedependent','bc',listofstr2str(nodeset_bc))
    config_local.set('pylithapp.timedependent','interfaces',listofstr2str(nodeset_interface))
    
    ### Remove unwanted sections
    
    remove_sections_beg('pylithapp.timedependent.bc',config_local)
    remove_sections_beg('pylithapp.timedependent.interfaces',config_local)
    remove_sections_beg('pylithapp.problem',config_local)
    remove_sections_beg('pylithapp.timedependent.materials',config_local)
        
    ### Start defining
        
    for nodeset in nodeset_bc:
        section_name='pylithapp.timedependent.bc.'+ nodeset
        config_local.add_section(section_name)
        config_local.set(section_name,'bc_dof',cfginput.get(nodeset + '.general','bc_dof'))
        config_local.set(section_name,'label',nodeset)
        config_local.set(section_name,'db_initial.label','Dirichlet BC')
    
    ### Set the Type 
    
    k=0
    config_local.add_section('pylithapp.timedependent.interfaces')
    for nodeset in nodeset_interface:
        config_local.set('pylithapp.timedependent.interfaces',nodeset,'pylith.faults.FaultCohesiveKin')
        section_name='pylithapp.timedependent.interfaces.'+ nodeset
        config_local.add_section(section_name)
        config_local.set(section_name,'label',nodeset)
        config_local.set(section_name,'id',100+k)
        k=k+1
        config_local.set(section_name,'quadrature.cell','pylith.feassemble.FIATSimplex')
        config_local.set(section_name,'quadrature.cell.dimension',1)
        
        config_local.add_section(section_name+'.eq_srcs.rupture')
        config_local.set(section_name+'.eq_srcs.rupture','slip_function','pylith.faults.ConstRateSlipFn')
        config_local.add_section(section_name+'.eq_srcs.rupture.slip_function') 
        config_local.set(section_name+'.eq_srcs.rupture.slip_function','slip_time','spatialdata.spatialdb.UniformDB')
        config_local.set(section_name+'.eq_srcs.rupture.slip_function','slip_time.label','Slip time')
        config_local.set(section_name+'.eq_srcs.rupture.slip_function','slip_time.values','[slip-time]')
        config_local.set(section_name+'.eq_srcs.rupture.slip_function','slip_time.data','[0.0*year]')
         
        if cfginput.get(nodeset + '.general','spatialdb')=='file':
            spatialdb_name=nodeset+'.spatialdb'
            config_local.set(section_name+'.eq_srcs.rupture.slip_function','slip_rate.iohandler.filename',spatialdb_name)
            config_local.set(section_name+'.eq_srcs.rupture.slip_function','slip_rate.query_type','linear')
            config_local.set(section_name+'.eq_srcs.rupture.slip_function','slip_rate.label','Final Slip')
        elif cfginput.get(nodeset + '.general','spatialdb')=='uniform':
            config_local.set(section_name+'.eq_srcs.rupture','slip_rate.label','Final Slip')
            config_local.set(section_name+'.eq_srcs.rupture','slip_rate.values','[left-lateral-slip,fault-opening]')
            left_lateral_slip=cfginput.getfloat(nodeset + '.parameters','left-lateral-slip')
            fault_opening=cfginput.getfloat(nodeset + '.parameters','fault-opening')
            config_local.set(section_name+'.eq_srcs.rupture','slip_rate.data','[%.1f*cm/year, %.1f*cm/year]'%(left_lateral_slip,fault_opening))
    
    ####### OUTPUTS
    
    #### Write output files for domains and subdomains
    
    section_name='pylithapp.problem.formulation.output.domain'
    config_local.add_section(section_name)
    config_local.set(section_name,'writer.filename','output/'+cfginput.get('Parameters','output_h5_name')+'.h5')
    
    section_name='pylithapp.problem.formulation.output.subdomain'
    config_local.add_section(section_name)
    config_local.set(section_name,'writer.filename','output/'+cfginput.get('Parameters','surface_nodeset_name')+'.h5')
             
    #### Write output files for faults
    
    for nodeset in nodeset_interface:
        section_name='pylithapp.problem.interfaces.'+nodeset+'.output'
        config_local.add_section(section_name)
        config_local.set(section_name,'writer','pylith.meshio.DataWriterHDF5')    
        config_local.set(section_name,'writer.filename','output/'+nodeset+'.h5')
        
    for block in list_of_blocks:
        section_name='pylithapp.timedependent.materials.'+block+'.output'
        config_local.add_section(section_name) 
        config_local.set(section_name,'writer.filename','output/'+block+'.h5')
             
    config_local.write(fic_out_local)
    fic_out_local.close()
    
    ################################################
    ######## Do the main generic file 
    ################################################
    
    ### Define the exodus file
    
    config_main.set('pylithapp.mesh_generator.reader','filename',exodus_name)
    config_main.set('pylithapp.mesh_generator.reader','coordsys.space_dim',2)
    config_main.set('pylithapp.timedependent','materials',listofstr2str(list_of_blocks))
    remove_sections_beg('pylithapp.timedependent.materials',config_main)
    config_main.add_section('pylithapp.timedependent.materials')
    
    for block in list_of_blocks:
        
        if cfginput.has_option(block + '.general','type_material'):
            if cfginput.get(block + '.general','type_material')=='elastic':
                config_main.set('pylithapp.timedependent.materials',block,'pylith.materials.ElasticPlaneStrain')
            elif cfginput.get(block + '.general','type_material')=='maxwell':
                config_main.set('pylithapp.timedependent.materials',block,'pylith.materials.MaxwellPlaneStrain')
            else:
                print "WARN : no type_material defined for " + block + ", elastic by default"
                config_main.set('pylithapp.timedependent.materials',block,'pylith.materials.ElasticPlaneStrain')
        else:
            print "WARN : no type_material defined for " + block + ", elastic by default"
            config_main.set('pylithapp.timedependent.materials',block,'pylith.materials.ElasticPlaneStrain')
    
        section_name='pylithapp.timedependent.materials.'+block
        id_block = list_of_blocks.index(block) + 1
        config_main.add_section(section_name)
        config_main.set(section_name,'label',block)
        config_main.set(section_name,'id',id_block)
        config_main.set(section_name,'db_properties.label',block+' properties')
        config_main.set(section_name,'db_properties.iohandler.filename',block+'.spatialdb')
        config_main.set(section_name,'quadrature.cell','pylith.feassemble.FIATSimplex')
        config_main.set(section_name,'quadrature.cell.dimension','2')
        section_output='pylithapp.timedependent.materials.'+block+'.output'
        config_main.add_section(section_output)
        config_main.set(section_output,'cell_filter','pylith.meshio.CellFilterAvg')
        config_main.set(section_output,'output_freq','time_step')
        config_main.set(section_output,'time_step','9.99999*year')
        config_main.set(section_output,'writer','pylith.meshio.DataWriterHDF5')
       
    surf_nodeset=cfginput.get('Parameters','surface_nodeset_name')
    config_main.set('pylithapp.problem.formulation.output.subdomain','label',surf_nodeset)
    
    #### Write PETSC
     # Nothing to write
    
    ### Write to output
    
    config_main.write(fic_out_main)
    fic_out_main.close()
    
    return None
    
#  _______ _    _ ______   ______ _   _ _____  
# |__   __| |  | |  ____| |  ____| \ | |  __ \ 
#    | |  | |__| | |__    | |__  |  \| | |  | |
#    | |  |  __  |  __|   |  __| | . ` | |  | |
#    | |  | |  | | |____  | |____| |\  | |__| |
#    |_|  |_|  |_|______| |______|_| \_|_____/    
#    
# ZONE DE TEST
cfpath = '/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/pylith_bigrun.cfg'

#import h5py
#f = h5py.File('/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/experience_mk2/exp_1/output/ground_surface.h5')
#
#f.keys()
#f['geometry'].items()
#arr = np.array(f['vertex_fields']['displacement'])

#clf()
##plt.plot(arr[1,:,0])
##plt.plot(arr[2,:,0],'r')
#plt.plot(arr[3,:,0] - arr[2,:,0],'k')
#
#
#            
#    

#import subprocess 
#p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE ,env=dic_env)
#command = pylith_path +"/bin/pylith" + ' /home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/experience_mk3/exp_993/interseismic.cfg'
#p.communicate( command )
#        dic['PATH'] = '/home/psakicki/THESE/CODES/pylith-2.0.3-linux-x86_64/bin:/opt/goa-6.3/bin:/opt/goa-6.3/bin/Linux-x86_64:/opt/goa-6.3/doris/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/psakicki/THESE/CODES/cTraceo/WORK/bin/:/home/psakicki/ginspc/bin:/home/psakicki/scripts:/home/psakicki/scripts:/home/psakicki/ginspc/bin::/home/psakicki/scripts:/home/psakicki/dynamo/scripts:/home/psakicki/GRGS_utilities/bin:/home/psakicki/INFO:/home/psakicki/GINO:.:/opt/gamit_globk/gamit/bin:/opt/gamit_globk/kf/bin:/opt/gamit_globk/com:/home/psakicki/THESE/CODES/SCRIPTS/'
#        p = Popen('',executable='/bin/bash', stdin=PIPE,env=dic)
#        p.communicate('/home/psakicki/THESE/CODES/pylith-2.0.3-linux-x86_64/bin/pylith ' + cf , shell = True, executable='/bin/bash')            
#        
#        except:
#            print 'hello world !'
#return None


#from subprocess import Popen, PIPE
#
#dic['PATH'] = '/home/psakicki/THESE/CODES/pylith-2.0.3-linux-x86_64/bin:/opt/goa-6.3/bin:/opt/goa-6.3/bin/Linux-x86_64:/opt/goa-6.3/doris/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/psakicki/THESE/CODES/cTraceo/WORK/bin/:/home/psakicki/ginspc/bin:/home/psakicki/scripts:/home/psakicki/scripts:/home/psakicki/ginspc/bin::/home/psakicki/scripts:/home/psakicki/dynamo/scripts:/home/psakicki/GRGS_utilities/bin:/home/psakicki/INFO:/home/psakicki/GINO:.:/opt/gamit_globk/gamit/bin:/opt/gamit_globk/kf/bin:/opt/gamit_globk/com:/home/psakicki/THESE/CODES/SCRIPTS/'
#
#p = Popen('',executable='/bin/bash', stdin=PIPE)#,env=dic)
#p = subprocess.Popen(['/bin/bash'], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
#p.communicate('source /home/psakicki/.bashrc') # pass commands to the opened shell
#p.communicate("echo $PATH") # pass commands to the opened shell
#
#p.communicate('cd $HOME') # pass commands to the opened shell
#p.communicate("pwd") # pass commands to the opened shell
#
#
#
#os.environ['PATH']=pylith_path+'/bin:'+os.environ['PATH']
#os.environ['PYTHONPATH']=pylith_path+'/lib/python2.7/site-packages:'+ pylith_path +'/lib64/python2.7/site-packages:'+ os.environ['PYTHONPATH']
#os.environ['LD_LIBRARY_PATH']=pylith_path+'lib:'+pylith_path+'lib64:'#+os.environ['LD_LIBRARY_PATH']
#

#dic_env = os.environ
#pylith_path='/home/psakicki/THESE/CODES/pylith-2.0.3-linux-x86_64/'
#dic_env['PATH'] = pylith_path+'/bin:'+ dic_env['PYTHONPATH'] 
#dic_env['PYTHONPATH'] = pylith_path+'/lib/python2.7/site-packages:'+ pylith_path +'/lib64/python2.7/site-packages:'+ dic_env['PYTHONPATH'] 
#dic_env['LD_LIBRARY_PATH'] = pylith_path+'lib:'+pylith_path+'lib64:' #+os.environ['LD_LIBRARY_PATH']
#
#
#import subprocess 
#p = subprocess.Popen('',executable='/bin/bash', stdin=subprocess.PIPE ,env=dic_env)
#command = pylith_path +"/bin/pylith" + ' /home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/experience_mk3/exp_993/interseismic.cfg'
#p.communicate( command )

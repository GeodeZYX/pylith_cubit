import StringIO
from collections import OrderedDict
import ConfigParser
import itertools
import numpy as np
import os
import netCDF4
import sys
import matplotlib.pyplot as plt

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
#        entity_name : the name of the entity
#        entity_type : the type of the entity (thank you captain obvious)
#        (all the annex keys must be removed in valdic2sdbdic_list)
#
#        ONE block/interface per dico
#        SEVERALS value type per dico
#
# a combi dico contains all information about ONE combinaition,
# the keys are like in the configfile '<name of the block/interface>.<name of the parameter>
#
# the discrimination of specific cases is in the configdic2valdic

def open_exodus(exodusin):
    if type(exodusin) is str:
        exodusout = netCDF4.Dataset(exodusin, 'a')
    else:
        exodusout = exodusin
    return exodusout
    
def open_configfile(configfilein):
    if type(configfilein) is str:    
        cfout = ConfigParser.ConfigParser()
        cfout.read(configfilein)
    else:
        cfout  = configfilein
    return cfout
    
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
    xy = np.vstack((x,y)).T
    xy2 = (xy[xy[:,1].argsort()])
    return np.flipud(xy2[:,0]), np.flipud(xy2[:,1].T)

    
def dist(xa,ya,xb,yb):
    return np.sqrt((xa-xb)**2 + (ya-yb)**2)
    
def cumul_dist(X,Y,i1,i2):
    finaldist = 0
    for i in range(i1,i2):
        finaldist = finaldist + dist(X[i],Y[i],X[i+1],Y[i+1])
    return finaldist
               
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
        
def flipper(X,Y):
    # flip X,Y array in order to get the deepest point at the end of the array
    if Y[-1] > Y[1]:
        X = np.flipud(X)
        Y = np.flipud(Y)
    return X,Y


def make_combi_4_lockzone_bkp(Xin,Yin,ymax,lmin,stepmin,inunit='km',outunit='km'):
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
    X,Y = sort_XY_nodeset(Xin,Yin)
    imax = np.min(np.flatnonzero(Y < ymax))
    print "ymax = " , ymax , " imax = " , imax
    rangei = np.arange(imax)
    xybound_list = []
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
        i_lasttoppt = c[0]
        Ncombi_valid += 1
        
    print 'end of locking combinations, valid / total :' ,Ncombi_valid, '/', Ncombi_tot
    return xybound_list

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


#def configfile2combidic_list(exodus,path,faultlist=['fault_top']):
#    exodus = open_exodus(exodus)
#    cf = ConfigParser.ConfigParser()
#    cf.read(path)
#    configdic = OrderedDict()
#    # Removing section with no variables intervals => with no '.'
#    sections = cf.sections()
#    for s in sections:
#        if s.find('.') == -1 or not cf.has_option(s,'variable'):
#            continue
#        if not cf.getboolean(s,'variable'):
#            val = [cf.getfloat(s,'fixed')]
#        elif cf.getfloat(s,'delta') == 0:
#            val = [cf.getfloat(s,'fixed')]
#            print "WARN : !!! variable parameter but delta = 0 !!!"
#        else:
#            truemin = np.min([cf.getfloat(s,'min'),cf.getfloat(s,'max')])
#            truemax = np.max([cf.getfloat(s,'min'),cf.getfloat(s,'max')])
#            val = list(np.arange(truemin,truemax +0.1, cf.getfloat(s,'delta')))
#        configdic[s] = val
#        
#    # MAKING THE COMBINATIONS OF THE LOCKING ZONE
#    for faultname in  faultlist:
#        print faultname
#        ymax = cf.getfloat(faultname,'maxi_locked_depth')
#        lmin = cf.getfloat(faultname,'mini_locked_length')
#        X,Y,_,_ = get_XY_of_nodeset(exodus,faultname)
#        configdic[faultname + '.combiXY'] = make_combi_4_lockzone(X/1000,Y/1000,ymax,lmin,inunit='km',outunit='km')
#        # SPECIFIC
#        for badfield in [ k for k in configdic.keys() if 'maxi_locked_depth' in k ]:
#            print badfield , 'aaa'
#            configdic.pop(badfield,None)
#        # CARTESIAN PRODUCT OF A CONFIG DIC => MANY COMBIDIC
#        combidiclis = product_on_dic(configdic)
#        print "total combinations : ", len(combidiclis)
#        
#    return combidiclis,configdic

def configfile2combidic_list(path):
    cf = ConfigParser.ConfigParser()
    cf.read(path)
    sections = cf.sections()    
    exodus = open_exodus(cf.get('Parameters','exodus_name'))
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
        X,Y,_,_ = get_XY_of_nodeset(exodus,faultname)
        paramdic[faultname + '.combiXY'] = make_combi_4_lockzone(X/1000,Y/1000,ymax,lmin,stepmin,inunit='km',outunit='m')
    # SPECIFIC
#    for badfield in [ k for k in configdic.keys() if 'maxi_locked_depth' in k ]:
#        print badfield , 'aaa'
#        configdic.pop(badfield,None)
    # CARTESIAN PRODUCT OF A CONFIG DIC => MANY COMBIDIC
    combidiclis = product_on_dic(paramdic)
    print "total combinations : ", len(combidiclis)
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
        valdic['x'] = []
        valdic['y'] = []
        for k,v in combidic.iteritems():
            if e in k: # here is a parameter of the entity
                p = k.split('.',1)[1]
                if p == 'combiXY': #SPECIFIC CASES
                    delta = 0.1
                    valdic['x'] = [0,v[0][0],v[0][0]-delta,v[1][0],v[1][0]-delta,0]
                    valdic['y'] = [99, v[0][1] ,v[0][1]-delta ,v[1][1],v[1][1]-delta,-999]
                elif p == 'fault-opening':
                    if check_if_combiXY(combidic,e)[0]:
                        valdic[p] = [v] * 6
                    else:
                        valdic['x'] = [0,0]
                        valdic['y'] = [99,-999]
                        valdic[p] = [v,v]   
                elif p == 'left-lateral-slip':
                    if check_if_combiXY(combidic,e)[0]:
                        valdic[p] = [v,v,0,0,v,v]     
                    else:
                        valdic['x'] = [0,0]
                        valdic['y'] = [99,-999]
                        valdic[p] = [v,v]
                else: #GENERIC CASE
                    valdic['x'] = [0]
                    valdic['y'] = [0]
                    if type(v) is int or float:
                        valdic[p] = [v]
                    else:
                        valdic[p] = v

        valdic_lis.append(valdic)  
    return valdic_lis
    
    
def combidic2valdic2_list(combidic,configfile,outunit='m'):
    cf = open_configfile(configfile)
    exodus = open_exodus(cf.get('Parameters','exodus_name'))
    if outunit == 'm':
        kk = 1
    elif outunit == 'km':
        kk = 1000
    else:
        kk = 1000
    entity_name_lis = []
    # getting the name of all entity (block/nodesets) in the combidico
    # (for the moment the keys are <entity>.<parameter>)
    for k,v in combidic.iteritems():
        entity_name_lis.append(k.split('.',1)[0])
    list(set(entity_name_lis))
    # iterating over all the entities
    valdic_lis = []
    for e in entity_name_lis:
        param_lis = []
        val_lis = []
        for p,v in combidic.items() :
            if e in p:
                param_lis.append(p.split('.')[1])
                val_lis.append(v)
        valdic = OrderedDict()
        valdic['entity_name'] = e
        valdic['x'] = []
        valdic['y'] = []
        
        # TYPE BLOCK => only ONE line
        if cf.get(e + '.general','type') in ('bloc','block'):
            valdic['x'] = [0]
            valdic['y'] = [0]
            for p,v in zip(param_lis,val_lis):
                valdic[p] = [v]
        # TYPE INTERFACE
        elif cf.get(e + '.general','type') in ('interface'):
            # locked case => speed = nominal speed and 0
            if cf.has_option(e + '.general','locked') and cf.getboolean(e + '.general','locked'):
                if not 'combiXY' in param_lis:
                    raise Exception('ERR : ' + ' is a locked interface but have no locking boundary (combiXY)')
                X,Y,_,_ = get_XY_of_nodeset(exodus,e)
#                X,Y = sort_XY_nodeset(X,Y)
                X,Y = flipper(X,Y)
                X = X/kk
                Y = Y/kk
                valdic['x'] = list(X)
                valdic['y'] = list(Y)
                for p,v in zip(param_lis,val_lis):
                    if p == 'combiXY':
                        continue
                    elif p == 'left-lateral-slip':
                        valdic[p] = []
#                        x1 = combidic[e + '.combiXY'][0][0]
#                        y1 = combidic[e + '.combiXY'][0][1]
#                        x2 = combidic[e + '.combiXY'][1][0]
#                        y2 = combidic[e + '.combiXY'][1][1]
                        i = combidic[e + '.combiXY'][0]
                        j = combidic[e + '.combiXY'][1]

                        for n,(x,y) in enumerate(zip(X,Y)):
                            if i <= n <= j:
                                valdic[p].append(0)
                            else:
                                valdic[p].append(v)
                    else:
                        valdic[p] = [v] * len(X)    
            # not locked case => speed = nominal speed
            else:
                X,Y,_,_ = get_XY_of_nodeset(exodus,e)
#                X,Y = sort_XY_nodeset(X,Y)
                X,Y = flipper(X,Y)
                X = X/kk
                Y = Y/kk
                valdic['x'] = list(X)
                valdic['y'] = list(Y)
                for p,v in zip(param_lis,val_lis):
                    valdic[p] = [v] * len(X)

        valdic_lis.append(valdic)  
    return valdic_lis
 
def check_if_combiXY(combidic,entity):
    for k,v in combidic.items():
        if entity in k:
            if 'combiXY' in k:
                return True , v
    return False , None
        

 
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
   
#def write_spatialdb2(listofdic,outpath,outname):
#    ''' This function take a LIST of sdb dico as argument '''
#    strm =  StringIO.StringIO()
#    protonumloc = list(set([len(d['values']) for d in listofdic ]))    
#    if len(protonumloc) != 1:
#        print "WARN : write_spatialdb : differents numlocs "
#    numlocs = protonumloc[0]
#    datadim = 0 # ???????? ADEFINIR !!!!!!!!
#    if sum([d['name'] == 'z' for d in listofdic]) != 0:
#        spacedim = 3
#    else:
#        spacedim = 2 
#    # finding the X value for the unit :
#        unitcoef = 1
#        for d in listofdic:
#            if d['name'] == 'x':
#                if d['unit'] == 'km':
#                    unitcoef = 1000
#    # write usual header    
#    strm.write('#SPATIAL.ascii 1\n')
#    strm.write('SimpleDB {\n')
#    strm.write('  num-values = {}\n'.format(len(listofdic) - spacedim))
#    namestr = ''
#    unitstr = ''
#    for d in listofdic:
#        # we must exclude the x,y,z name in the header
#        if d['name'] in ('x','y','z'):
#            continue
#        namestr = namestr + ' ' + d['name']
#        namestr = namestr.replace('_','-')
#        unitstr = unitstr + ' ' + d['unit']  
#    strm.write('  value-names = {}\n'.format(namestr))
#    strm.write('  value-units = {}\n'.format(unitstr))
#    
#    strm.write('  num-locs = {}\n'.format(numlocs))
#    strm.write('  data-dim = {}\n'.format(datadim))
#    
#    strm.write('  space-dim = {}\n'.format(spacedim))
# 
#    strm.write('  cs-data = cartesian {\n')
#    strm.write('    to-meters = {}\n'.format(unitcoef))
#    strm.write('    space-dim = {}\n'.format(spacedim))
#    strm.write('  } \n')
#    strm.write('} \n')
#    
#    stdformat = '{:9.3f} '
#    
#    for line in zip(*[ d['values'] for d in listofdic ]):
#        tempstr = stdformat * len(listofdic) + '\n'
#        strm.write(tempstr.format(*line))
#            
#    fil = open(outpath + '/' + outname + '.spatialdb','w')
#    
#    fil.write(strm.getvalue())
#    fil.close()
#            
#    return strm

        
def main_config2multiexp(exodus,configfile,outpath,local_cfg_generic,main_cfg_generic,output_cfg_local,output_cfg_main,outdirprefix='exp'):
    exodus = open_exodus(exodus)
    check_exodus_configdic(exodus,configfile)       
    combidiclis,confdic = configfile2combidic_list(configfile) 
    i = 0
    for combidic in combidiclis:
        i = i+1 
        print i
        valdiclis = combidic2valdic2_list(combidic,configfile)
        for valdic in valdiclis:
            outdirectory = outpath + '/' + outdirprefix + '_' + str(i)
            if not os.path.exists(outdirectory):
                os.makedirs(outdirectory)
            write_spatialdb(valdic, outdirectory)
            write_configfile(configfile,local_cfg_generic,main_cfg_generic,outdirectory,output_cfg_local,output_cfg_main)

   
# ================== WRITE CONFIG FILES ==================


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
    
def listofstr2str(listin):
    outstr = '['
    for s in listin:
        outstr  = outstr + s + ','
    outstr = outstr[:-1] + ']'
    return outstr
   
#    return [s for s in cf.sections() if ( cf.has_option(s,'type') and cf.get(s,'type')==typein) ]            

def write_configfile(input_cfgfile,local_generic,main_generic,outpath,output_local,output_main):
    
#    input_cfgfile='/Users/baillard/_Moi/Programmation/Post_Doc_LR/Pylith/input_param.cfg'
#    local_generic='/Users/baillard/_Moi/Programmation/Post_Doc_LR/Pylith/interseismic.cfg'
#    main_generic='/Users/baillard/_Moi/Programmation/Post_Doc_LR/Pylith/pylithapp.cfg'
#    outpath='/Users/baillard/_Moi/Programmation/Post_Doc_LR/Pylith/test/'
#    output_local='test.cfg'
#    output_main='pylithtest.cfg'
        
    config_local = ConfigParser.ConfigParser()
    config_main = ConfigParser.ConfigParser()
    cfginput = ConfigParser.ConfigParser()
    config_local.optionxform = str 
    config_main.optionxform = str 
    cfginput.optionxform = str 
    
    ######## Checking and opening blocks
    
    ### Check if output path exist and open file for writing
    
    if not os.path.isdir(outpath):
        os.mkdir(outpath)
    if not os.path.isfile(local_generic) or not os.path.isfile(main_generic):
        raise Exception('generic config file ' + local_generic + " don't exist")
        
    out_local=outpath + '/' + output_local # for the local.cfg file
    fic_out_local=open(out_local,'w')
    out_main=outpath  + '/' + output_main # for the main.cfg file
    fic_out_main=open(out_main,'w')
    
    ### Check and open generic files
    
    if ( not os.path.isfile(local_generic) or not os.path.isfile(main_generic)):
        print 'WARNING: generic configuration files were not found, change path below and try again! \n' + local_generic
        raise Exception  
    
    ### Read generic files and Exodus
    
    config_local.read(local_generic)
    config_main.read(main_generic)
    cfginput.read(input_cfgfile)
    
    exodus_name=cfginput.get('Parameters','exodus_name')
     
    ### Change settings
    
    list_of_blocks = blocksname_list(exodus_name)
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
    config_local.set(section_name,'writer.filename','output/'+cfginput.get('Parameters','h5_output_name')+'.h5')
    
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
    
# ZONE DE TEST
prepath=''
#prepath='/run/user/1000/gvfs/sftp:host=calipso.univ-lr.fr,user=psakicki'

filenameExodus = prepath + "/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/output_MESH/mesh__solo1.exo"
cfpath = prepath + '/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/pylith_bigrun.cfg'
outpath = '/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/test_sdb/'

input_cfgfile='/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/pylith_bigrun.cfg'
outpath = '/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/experience_mk2'
output_main='pylithapp.cfg'
output_local='interseismic.cfg'
local_generic='/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/generic_pylith_config/interseismic.cfg'
main_generic='/home/psakicki/THESE/MODEL_GWADA/WORKING_DIR/generic_pylith_config/pylithapp.cfg'
 
cf = open_configfile(input_cfgfile)

a,b = configfile2combidic_list(input_cfgfile)
combidic2valdic2_list(a[10],input_cfgfile)

faultlist = [s.split('.')[0] for s in cf.sections() if cf.has_option(s,'locked') and cf.getboolean(s,'locked') ]

X,Y,_,_ = get_XY_of_nodeset(filenameExodus,'bottom_litho_DP')
X2,Y2 = sort_XY_nodeset(X,Y)

plt.plot(X,Y)


Xft,Yft ,_,_= get_XY_of_nodeset(filenameExodus,'ground_surface')

plt.plot(Xft,Yft,'k')

nodesetname_list(filenameExodus)

main_config2multiexp(filenameExodus,cfpath,outpath,local_generic,main_generic,output_local,output_main) 

#X,Y,_,_ = get_XY_of_nodeset(filenameExodus,'fault_top')

#write_configfile(input_cfgfile,local_generic,main_generic,outpath,output_local,output_main)
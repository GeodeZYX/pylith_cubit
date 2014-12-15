import geom_cubit_lib as gcl 
import cubit
import numpy as np
import ConfigParser
import os
import sys
import subprocess
from StringIO import StringIO 


# VERSION 141209a
#
#===================================================================
#             _                  _           _             
#            (_)                | |         (_)            
#  _ __  _ __ _  ___  _ __    __| | ___  ___ _  __ _ _ __  
# | '_ \| '__| |/ _ \| '__|  / _` |/ _ \/ __| |/ _` | '_ \ 
# | |_) | |  | | (_) | |    | (_| |  __/\__ \ | (_| | | | |
# | .__/|_|  |_|\___/|_|     \__,_|\___||___/_|\__, |_| |_|
# | |                                           __/ |      
# |_|                                          |___/       
#
#===================================================================

def common_pt(X1,Y1,X2,Y2): 
    # for a boolean (good or not), check the length of the return "common"
    XY1_pts = zip(X1,Y1)
    XY2_pts = zip(X2,Y2)
    common=list(set(XY1_pts).intersection(XY2_pts))  
    return common

def check_if_common_pt(X1,Y1,X2,Y2):
    common = common_pt(X1,Y1,X2,Y2)
    if len(common) > 0:
        print "OK : ", len(common), " common point(s)"
        if len(common) > 1:
            print "but WARN : more that 1 point !!!"
    else:
        print "WARN : no common point in the bathy and the slab file !!!"
    return None
    
def shift_thickness(X,Y,thick):
    alpha = np.arctan2(np.diff(X),np.diff(Y))
    alpha = np.insert(alpha,0,1)
    x_new=[]
    y_new=[]  
    for a,x,y in zip(alpha,X,Y):
        if a <=0:
            dx=np.cos(a)*thick
            x_new.append(x-dx)
            dy=np.sin(a)*thick
            y_new.append(y+dy)
        else:
            dx=np.cos(a)*thick
            x_new.append(x+dx)
            dy=np.sin(a)*thick
            y_new.append(y-dy)  
    return x_new,y_new,alpha
    
def slab_direction(Xslab,Yslab):
    ixl = Xslab.argmin()
    ixr = Xslab.argmax()
    yl = Yslab[ixl]
    yr = Yslab[ixr]
    if yl > yr:
        return 1
    else:
        return -1

def is_sorted(listin):
    return (np.sort(listin) == np.array(listin)).all()
    
def sort_B_according_A(A,B):
    Bout = [b for (a,b) in sorted(zip(A,B))]
    Aout = sorted(A)
    return Aout,Bout

def make_bounding_points(Xbathy,Ybathy,Xslab,Yslab,OPepaiss,DPepaiss,bottom,leftext,rightext):
    slabdir = slab_direction(Xslab,Yslab)
    if slabdir == -1:
        leftepaiss = OPepaiss
        rightepaiss = DPepaiss
    else:
        leftepaiss = DPepaiss
        rightepaiss = OPepaiss   
    if is_sorted(Xbathy):
        print "WARN : Xbathy not sorted !!!"
        Xbathy,Ybathy= sort_B_according_A(Xbathy,Ybathy)
    xl = np.min(Xbathy)
    yl = Ybathy[Xbathy.argmin()]
    xlnew = xl - np.abs(leftext)
    xr = np.max(Xbathy)
    yr = Ybathy[Xbathy.argmax()]
    xrnew = xr + np.abs(rightext)
    left_bathy = [xlnew,yl]  
    left_lit_ast = [xlnew,-np.abs(leftepaiss)]
    left_ast_bot = [xlnew,-np.abs(bottom)]
    right_bathy = [xrnew,yr]    
    right_lit_ast = [xrnew,-np.abs(rightepaiss)]
    right_ast_bot = [xrnew,-np.abs(bottom)]
    if slabdir == -1:
        leftplate ='OP'
        rightplate='DP'
    else:
        leftplate ='DP'
        rightplate='OP'
    outdic = {} 
    outdic['vBathy_'+leftplate] = left_bathy
    outdic['vBottom_Litho_'+leftplate] = left_lit_ast
    outdic['vBottom_Astheno_'+leftplate] = left_ast_bot
    outdic['vBathy_'+rightplate] = right_bathy
    outdic['vBottom_Litho_'+rightplate] = right_lit_ast
    outdic['vBottom_Astheno_'+rightplate] = right_ast_bot             
    common = common_pt(Xbathy,Ybathy,Xslab,Yslab)
    outdic['vscratch'] = np.array([common[0][0],outdic['vBottom_Litho_OP'][1]])    
    return outdic
    
def make_backstop(Xbathy,Ybathy,Xslab,Yslab,distance,angle,dupdown=100,upsafety=0.05):
    direction = slab_direction(Xslab,Yslab)
    common = common_pt(Xbathy,Ybathy,Xslab,Yslab)
    xup = common[0][0] + direction*distance
    yup = np.interp(xup,Xbathy,Ybathy) + upsafety 
    ydown = - dupdown * np.sin(( np.pi / 180 ) * angle) + yup
    xdown = -direction * dupdown *np.cos( (np.pi / 180 ) * angle) + xup
    Xbackstop = [xup,xdown]
    Ybackstop = [yup,ydown]
    return Xbackstop,Ybackstop

def bathyslab_2_OPDPtop(Xbathy,Ybathy,Xslab,Yslab):
    direction = slab_direction(Xslab,Yslab)
    common = common_pt(Xbathy,Ybathy,Xslab,Yslab)
    Xop , Yop = [],[]
    Xdp , Ydp = [],[]
    for xb,yb in zip(Xbathy,Ybathy):
        if direction == -1:
            if xb <= common[0][0]:
                Xop.append(xb)
                Yop.append(yb)
            else:
                Xdp.append(xb)
                Ydp.append(yb)
        else:
            if xb >= common[0][0]:
                Xop.append(xb)
                Yop.append(yb)
            else:
                Xdp.append(xb)
                Ydp.append(yb)
    Xdp = list(Xdp) + list(Xslab)
    Ydp = list(Ydp) + list(Yslab)  
    Xop2,Yop2 = sort_B_according_A(Xop,Yop)
    Xdp2,Ydp2 = sort_B_according_A(Xdp,Ydp)
    return Xop2,Yop2,Xdp2,Ydp2
    
def add_point_in_list(Xlis,Ylis,point):
    Xlis.append(point[0])
    Ylis.append(point[1])
    Xlis,Ylis = sort_B_according_A(Xlis,Ylis)
    return Xlis,Ylis

#===================================================================
#          _____ _    _ ____ _____ _______ 
#         / ____| |  | |  _ \_   _|__   __|
#        | |    | |  | | |_) || |    | |   
#        | |    | |  | |  _ < | |    | |   
#        | |____| |__| | |_) || |_   | |   
#         \_____|\____/|____/_____|  |_|   
#
#===================================================================

def read_file_2_listofvertices(path,unit='*km'):
    filein = open(path)
    Vlis_out = []
    for line in filein:
        fields = line.split()
#        v = cubit.create_vertex(float(fields[0]),float(fields[1]),0)    
        komand = "create vertex x { "+str(float(fields[0]))+ unit + "} y { "+str(float(fields[1])) + unit + "}"
        cubit.silent_cmd(komand)
        v = cubit.vertex(cubit.get_last_id('vertex'))
        Vlis_out.append(v)  
    return Vlis_out 
    
def read_file_2_curve(path,name=''):
    Vlist = read_file_2_listofvertices(path)
    finalcmd = "create curve spline from vertex " + str(Vlist[0].id()) + " to " + str(Vlist[-1].id())
    cubit.silent_cmd(finalcmd)
    C = cubit.curve(cubit.get_last_id('curve'))
    if name == '':
        C.entity_name(os.path.basename(path))
    else:
        C.entity_name(name)
    return C
    
def list_2_listofvertices(X,Y,unit='*km'):
    Vlis_out = []
    for x,y in zip(X,Y):
        komand = "create vertex x { "+str(x)+ unit + "} y { "+str(y) + unit + "}"
        cubit.silent_cmd(komand)
        v = cubit.vertex(cubit.get_last_id('vertex'))
        Vlis_out.append(v)  
    return Vlis_out 
    
def list_2_curve(X,Y,name='',unit='*km'):
    Vlist = list_2_listofvertices(X,Y,unit)
    finalcmd = "create curve spline from vertex " + str(Vlist[0].id()) + " to " + str(Vlist[-1].id())
    cubit.silent_cmd(finalcmd)
    C = cubit.curve(cubit.get_last_id('curve'))
    if name != '':
        C.entity_name(name)
    return C
    
def dico_2_listofvertices(dico,unit='*km'):
    Vlis_out = []
    for name , xy in dico.iteritems():
        komand = "create vertex x { "+str(xy[0])+ unit + "} y { "+str(xy[1]) + unit + "}"
        cubit.silent_cmd(komand)
        v = cubit.vertex(cubit.get_last_id('vertex'))
        v.entity_name(name)
        Vlis_out.append(v)  
    return Vlis_out 
    
def vertices_2_curve(Va,Vb,name=''):
    C = cubit.create_curve(Va,Vb)
    C.entity_name(name)
    return C

def find_list(Clist,desc,desctype='ID'):
    for c in Clist:
        if desctype == 'ID' and c.id() == desc:
            return c
        if desctype == 'name' and c.entity_name() == desc:
            return c
    print "if no return, check the Descriptor Type"
    return None 

def find_extrema(E,extr = 'u'):
    Vlis = []
    extlist = E.bounding_box()    
    if extr == 'u':
        yok = extlist[4]
    elif extr == 'd':
        yok = extlist[1]
    elif extr == 'r':
        xok = extlist[3]
    elif extr == 'l':
        xok = extlist[0]   
    for v in E.vertices():
        if extr == 'u' or extr =='d':
            if v.coordinates()[1] == yok:
                Vlis.append(v)
        if extr == 'l' or extr =='r':
            if v.coordinates()[0] == xok:
                Vlis.append(v)              
    if len(Vlis) == 1:
        return Vlis[0]
    else:
        return Vlis  

def split(Ca,Cb,Ca1name='_1',Ca2name='_2'):
    # Ca : crossing curve (the cutted one) 
    # Cb : crossed curve (the unchanged one)
    # Ca & Cb are objects
    # If Ca[1/2]name starts with a underscore, we keep the original name as a prefix
    #SPLIT
    Caname = Ca.entity_name()
    command = 'split curve ' + str(Ca.id()) + ' crossing curve ' + str(Cb.id())
    print command 
    cubit.cmd(command)
    # GET IDs
    lastID2 = cubit.get_last_id('curve')
    lastID1 = lastID2 - 1
    Ca1 = cubit.curve(lastID1)
    # RENAMING
    if Ca1name[0] == '_':
#        Ca1.entity_name(Caname + Ca1name)
        command = 'curve ' + str(lastID1) + ' rename "' + Caname + Ca1name +'"'        
    else:
#        Ca1.entity_name(Ca1name)
        command = 'curve ' + str(lastID1) + ' rename "' + Ca1name +'"'        
    cubit.cmd(command)    
    Ca2 = cubit.curve(lastID2)
    if Ca2name[0] == '_':
#        Ca2.entity_name(Caname + Ca2name)
        command = 'curve ' + str(lastID2) + ' rename "' + Caname + Ca2name +'"'        
    else:
#        Ca2.entity_name(Ca2name)
        command = 'curve ' + str(lastID2) + ' rename "' + Ca2name +'"'     
    cubit.cmd(command)    
    return Ca1,Ca2

    
def double_split(Ca,Cb,Ca1name='_1',Ca2name='_2',Cb1name='_1',Cb2name='_2'):
    Ca1,Ca2 = split(Ca,Cb,Ca1name,Ca2name)
    Cb1,Cb2 = split(Cb,Ca1,Cb1name,Cb2name)   
    return Ca1,Ca2,Cb1,Cb2

def double_split_desc(Ca,Cb,Ca1name='_1',Ca2name='_2',Cb1name='_1',Cb2name='_2'):
    if type(Ca) is int:
        Ca1,Ca2,Cb1,Cb2 = double_split(cubit.curve(Ca),cubit.curve(Cb),Ca1name,Ca2name,Cb1name,Cb2name)
    elif type(Ca) is str:
        Ca1,Ca2,Cb1,Cb2 = double_split(cubit.curve(cubit.get_id_from_name(Ca)),cubit.curve(cubit.get_id_from_name(Cb)),Ca1name,Ca2name,Cb1name,Cb2name)
    return Ca1,Ca2,Cb1,Cb2  
   
def create_surface_desc(Dlist,name=''):
    # Dlist : liste des DESCRIPTEURS des courbes qui composent la future surface ([int] ou [str])
    OKlist = []
    if type(Dlist[0]) is int:
        for d in Dlist:
            OKlist.append(cubit.curve(d))
    elif type(Dlist[0]) is str:
        for d in Dlist:
            OKlist.append(cubit.curve(cubit.get_id_from_name(d)))
    S = cubit.create_surface(OKlist)
    if name != '':
        S.entity_name(name)
    return S
    
def curver_desc(Cdesc,dx,unit='*km'):
    if type(Cdesc) is int:
        komand = "curve " + str(Cdesc) + " size {" + str(dx) + unit + "}"
    elif type(Cdesc) is str:
        komand = "curve " + str(cubit.get_id_from_name(Cdesc)) + " size {" + str(dx) + unit +"}"    
    print komand
    cubit.silent_cmd(komand)
    return None
    
def curver_bias(C,dx,biasfactor,direct='l',unit='*km'):    
    Vtup = cubit.get_relatives("curve", C.id(), "vertex")
    if cubit.vertex(Vtup[1]).coordinates()[0] > cubit.vertex(Vtup[0]).coordinates()[0]:
        vr = cubit.vertex(Vtup[1])
        vl = cubit.vertex(Vtup[0])
    else:
        vr = cubit.vertex(Vtup[0])
        vl = cubit.vertex(Vtup[1])
    if direct == 'l':
        vstart = vr
    elif direct == 'r':
        vstart = vl
    komand = "curve " + str(C.id()) + " scheme bias fine size {" + str(dx) + unit + "} factor {" + str(biasfactor) + "} start vertex " + str(vstart.id())
    print komand
    cubit.silent_cmd(komand)
    return None
    
def curver_start_end(C,dxstrt,dxend,direct='l',unit='*km'):
    Vtup = cubit.get_relatives("curve", C.id(), "vertex")
    if cubit.vertex(Vtup[1]).coordinates()[0] > cubit.vertex(Vtup[0]).coordinates()[0]:
        vr = cubit.vertex(Vtup[1])
        vl = cubit.vertex(Vtup[0])
    else:
        vr = cubit.vertex(Vtup[0])
        vl = cubit.vertex(Vtup[1])
    if direct == 'l':
        vstart = vr
    elif direct == 'r':
        vstart = vl       
    komand = "curve " + str(C.id()) + " scheme bias fine size {" + str(dxstrt) + unit + "} coarse size  {" + str(dxend) + unit + "} start vertex " + str(vstart.id())
    print komand
    cubit.silent_cmd(komand)
    return None
   
def curver_bias_desc(Cdesc,dx,biasfactor,direction='l'):
    if type(Cdesc) is int:
        curver_bias(cubit.curve(Cdesc),dx,biasfactor,direction)
    elif type(Cdesc) is str:
        curver_bias(cubit.curve(cubit.get_id_from_name(Cdesc)),dx,biasfactor,direction)
    return None
    
def curver_start_end_desc(Cdesc,dx1,dx2,direction='l'):
    if type(Cdesc) is int:
        curver_start_end(cubit.curve(Cdesc),dx1,dx2,direction)
    elif type(Cdesc) is str:
        curver_start_end(cubit.curve(cubit.get_id_from_name(Cdesc)),dx1,dx2,direction)
    return None    

def create_group_nodeset_desc(Clist,grpname,nodesetID=0):
    while nodesetID == 0 or (nodesetID in cubit.get_nodeset_id_list()):
        nodesetID = nodesetID +1
    for c in Clist:
        cubit.silent_cmd('group "' + grpname + '" add node in ' + c) 
    cubit.silent_cmd('nodeset ' + str(nodesetID) + ' group ' + grpname)
    cubit.silent_cmd('nodeset ' + str(nodesetID) + ' name "' + grpname + '"' )


# Fonctions qu'il faut redefinir via des cmds APREPRO
# parce que les methodes Cubit marchent pas ... (?!?)
    
def destroy_curve(Cin):
    # C is a object
    command = 'delete curve ' + str(Cin.id())  
    cubit.silent_cmd(command)
    return None 
    
def destroy_curve_desc(Cdesc):
    if type(Cdesc) is int:
        destroy_curve(cubit.curve(Cdesc))
    elif type(Cdesc) is str:
        destroy_curve(cubit.curve(cubit.get_id_from_name(Cdesc)))
    return None   
    
def rename_curve(C,newname):
    # C is a object
    command = 'curve ' + str(C.id()) + ' rename "' + newname + '"'
    cubit.silent_cmd(command)
    return None  
    
def rename_curve_desc(Cdesc,newname):
    if type(Cdesc) is int:
        rename_curve(cubit.curve(Cdesc),newname)
    elif type(Cdesc) is str:
        rename_curve(cubit.curve(cubit.get_id_from_name(Cdesc)),newname)
    return None
    
# =====================================================================   
#       __  __          _____ _   _ 
#      |  \/  |   /\   |_   _| \ | |
#      | \  / |  /  \    | | |  \| |
#      | |\/| | / /\ \   | | | . ` |
#      | |  | |/ ____ \ _| |_| |\  |
#      |_|  |_/_/    \_\_____|_| \_|
#    
# =====================================================================   
                                                                    
    
def geom_cubit_main(configfile):
    if not os.path.isfile(configfile) :
        raise Exception("ERR: CONFIG FILE DON'T EXIST !!!")
    #  READ CONFIGFILE
    cf = ConfigParser.ConfigParser()
    cf.read(configfile)
    bathyfile = cf.get('Input','bathyfile')
    slabfile = cf.get('Input','slabfile')
    left_extension=cf.getfloat('Parameters','left_extension') # (always positive)
    right_extension=cf.getfloat('Parameters','right_extension')
    down_limit=cf.getfloat('Parameters','down_limit') # (always positive)
    EET_OP=cf.getfloat('Parameters','eet_op') #  Thickness of the OP Litho
    EET_DP=cf.getfloat('Parameters','eet_dp') #  Thickness of the DP Litho
#    backstop_bool=cf.getboolean('Parameters','backstop_bool') # True or False
    x_backstop=cf.getfloat('Parameters','backstop_distance')
    dip_backstop=cf.getfloat('Parameters','backstop_angle')
    output_path=cf.get('Output','output_path')
    output_file_prefix=cf.get('Output','output_prefix')
    output_file_suffix=cf.get('Output','output_suffix')
    output_pathfile = output_path + '/'  + output_file_prefix + '_'  + output_file_suffix + '.exo'
    dx1=cf.getfloat('Parameters','dx_fine')
    dx2=cf.getfloat('Parameters','dx_coarse')
    # =======================================================================
    #  PRIOR DESIGN
    # =======================================================================
    # Loading the BATHY and SLAB files
    XYbathy = np.loadtxt(bathyfile)
    XYslab = np.loadtxt(slabfile)
    # Creating the bottom of the slab
    XslabBot,YslabBot,_ = gcl.shift_thickness(XYslab[:,0],XYslab[:,1],EET_DP)
    # Creating the bounding points in a dictionnary
    outpts = gcl.make_bounding_points(XYbathy[:,0],XYbathy[:,1],XYslab[:,0],XYslab[:,1],EET_OP,EET_DP,down_limit,left_extension,right_extension)
    # Creating the backstop
    Xbackstop,Ybackstop = gcl.make_backstop(XYbathy[:,0],XYbathy[:,1],XYslab[:,0],XYslab[:,1],x_backstop,dip_backstop,dupdown=150)
    # Adding the 'Bottom_Litho_DP' point to the Litho_DP Bottom
    gcl.add_point_in_list(XslabBot,YslabBot,outpts['vBottom_Litho_DP'])
    # Changing the description concept : Bathy & Slab ==> OP & DP
    Xop,Yop,Xdp,Ydp = gcl.bathyslab_2_OPDPtop(XYbathy[:,0],XYbathy[:,1],XYslab[:,0],XYslab[:,1])
    # Adding the limit points of OP & DP
    Xop,Yop = gcl.add_point_in_list(Xop,Yop,outpts['vBathy_OP'])
    Xdp,Ydp = gcl.add_point_in_list(Xdp,Ydp,outpts['vBathy_DP'])
    # PLOT FOR SECURITY
    #plt.clf()
    #plt.axis('equal')
    #plt.plot(XYbathy[:,0],XYbathy[:,1],'+b')
    #plt.plot(XYslab[:,0],XYslab[:,1],'r+-')
    #plt.plot(XslabBot,YslabBot,'k+-')
    #plt.plot(Xbackstop,Ybackstop,'*-y')
    #for p in outpts.viewvalues():
    #    plt.plot(p[0],p[1],'*k')
    #plt.plot(Xop,Yop,'b-')
    #plt.plot(Xdp,Ydp,'r-')
    # =======================================================================
    #  CUBIT MESHING
    # =======================================================================
    # Initalisation
    cubit.cmd('reset')  
    cubit.cmd("#{Units('si')}")  
    # Chargement des courbes
    gcl.list_2_curve(Xop,Yop,'Top_Litho_OP')
    gcl.list_2_curve(Xdp,Ydp,'Top_Litho_limit_DP')
    gcl.list_2_curve(XslabBot,YslabBot,'Bottom_Litho_DP')
    gcl.list_2_curve(Xbackstop,Ybackstop,'Backstop')
    # Chargement des Points Isoles
    gcl.dico_2_listofvertices(outpts)
    ## fabrication de "courbes-segments" a partir de vertex
    v_Bottom_Litho_OP = cubit.vertex(cubit.get_id_from_name('vBottom_Litho_OP'))
    v_Bathy_DP = gcl.find_extrema(cubit.curve(cubit.get_id_from_name("Top_Litho_limit_DP")),'r')
    v_Bathy_OP = gcl.find_extrema(cubit.curve(cubit.get_id_from_name('Top_Litho_OP')),'l')
    v_Bottom_Litho_DP = gcl.find_extrema(cubit.curve(cubit.get_id_from_name("Bottom_Litho_DP")),'r')
    v_scratch = cubit.vertex(cubit.get_id_from_name('vscratch'))
    v_Bottom_Astheno_OP = cubit.vertex(cubit.get_id_from_name('vBottom_Astheno_OP'))
    v_Bottom_Astheno_DP = cubit.vertex(cubit.get_id_from_name('vBottom_Astheno_DP'))
    # 
    gcl.vertices_2_curve(v_Bottom_Litho_OP,v_Bottom_Astheno_OP,'Edge_Astheno_OP')
    gcl.vertices_2_curve(v_Bottom_Litho_OP,v_scratch,'Bottom_Litho_OP_PROTO')
    gcl.vertices_2_curve(v_Bottom_Litho_OP,v_Bathy_OP,'Edge_Litho_OP')
    gcl.vertices_2_curve(v_Bottom_Astheno_OP,v_Bottom_Astheno_DP,'Bottom_PROTO')
    gcl.vertices_2_curve(v_Bathy_DP,v_Bottom_Litho_DP,'Edge_Litho_DP')
    gcl.vertices_2_curve(v_Bottom_Litho_DP,v_Bottom_Astheno_DP,'Edge_Astheno_DP') 
    ## Split des courbes 
    gcl.double_split_desc('Bottom_Litho_OP_PROTO','Top_Litho_limit_DP')
    gcl.double_split_desc("Bottom_PROTO","Top_Litho_limit_DP_1")
    gcl.double_split_desc("Bottom_PROTO_2","Bottom_Litho_DP")
    gcl.double_split_desc("Backstop","Top_Litho_limit_DP_2")
    gcl.double_split_desc("Top_Litho_OP","Backstop_1")
    gcl.double_split_desc("Top_Litho_limit_DP_2_2","Top_Litho_OP_2")
    #    ## suppression des petits morceaux
    gcl.destroy_curve_desc("Top_Litho_limit_DP_1_1")
    gcl.destroy_curve_desc("Bottom_Litho_DP_1")
    gcl.destroy_curve_desc("Backstop_2")
    gcl.destroy_curve_desc("Bottom_Litho_OP_PROTO_2")
    gcl.destroy_curve_desc("Backstop_1_1")
#    ## renommage des courbes
    gcl.rename_curve_desc(7,'Edge_Litho_OP')
    gcl.rename_curve_desc(5,'Edge_Astheno_OP') 
    gcl.rename_curve_desc(9,'Edge_Litho_DP')
    gcl.rename_curve_desc(10,'Edge_Astheno_DP')    
    gcl.rename_curve_desc("Top_Litho_limit_DP_2_1",'Contact_OP_DP')
    gcl.rename_curve_desc("Top_Litho_limit_DP_2_2_1" ,'Contact_Prism_DP')
    gcl.rename_curve_desc("Backstop_1_2",'Contact_Prism_OP')
    gcl.rename_curve_desc("Top_Litho_OP_2_1",'Top_Litho_DP')
    gcl.rename_curve_desc("Top_Litho_OP_2_2",'Top_Prism')
    gcl.rename_curve_desc("Top_Litho_OP_1",'Top_Litho_OP')
    gcl.rename_curve_desc("Bottom_PROTO_2_2",'Bottom_Astheno_DP')
    gcl.rename_curve_desc("Bottom_Litho_DP_2",'Bottom_Litho_DP')
    gcl.rename_curve_desc("Bottom_PROTO_1",'Bottom_Astheno_OP')
    gcl.rename_curve_desc("Bottom_Litho_OP_PROTO_1" ,'Bottom_Litho_OP')
    gcl.rename_curve_desc("Bottom_PROTO_2_1" ,'Front_Litho_DP')
    gcl.rename_curve_desc("Top_Litho_limit_DP_1_2",'Astheno_Litho_Contact')
#    # creation des surfaces
    gcl.create_surface_desc(["Edge_Astheno_DP","Bottom_Astheno_DP","Bottom_Litho_DP"])
    gcl.create_surface_desc(["Front_Litho_DP","Bottom_Litho_DP","Edge_Litho_DP","Top_Litho_DP","Contact_Prism_DP","Contact_OP_DP","Astheno_Litho_Contact"])
    gcl.create_surface_desc(["Bottom_Astheno_OP", "Astheno_Litho_Contact", "Bottom_Litho_OP", "Edge_Astheno_OP"])
    gcl.create_surface_desc(["Bottom_Litho_OP","Edge_Litho_OP","Top_Litho_OP","Contact_Prism_OP","Contact_OP_DP"])
    gcl.create_surface_desc(["Top_Prism","Contact_Prism_DP","Contact_Prism_OP"])   
    # renomage des surfaces
    cubit.surface(1).entity_name('Astheno_DP')
    cubit.surface(2).entity_name('Litho_DP')
    cubit.surface(3).entity_name('Astheno_OP')
    cubit.surface(4).entity_name('Litho_OP')
    cubit.surface(5).entity_name('Prisme')
    # Fusion des surfaces
    cubit.cmd('delete vertex all')
    cubit.cmd('imprint all')
    cubit.cmd('merge all')
    cubit.cmd('stitch volume all')
    cubit.cmd('surface all scheme trimesh')
    cubit.cmd('curve all scheme default')
    cubit.cmd('surface all sizing function none')
    gcl.rename_curve_desc(45,'Astheno_Litho_Contact')
    # nouveau decoupage des courbes
    gcl.curver_desc("Contact_Prism_DP",dx1)
    gcl.curver_desc("Contact_Prism_OP",dx1)
    gcl.curver_desc("Top_Prism",dx1)
    gcl.curver_desc("Bottom_Litho_DP",dx2)
    gcl.curver_desc("Bottom_Astheno_OP",dx2)
    gcl.curver_desc("Front_Litho_DP",dx2)
    gcl.curver_desc("Bottom_Astheno_DP",dx2)
    gcl.curver_desc("Edge_Astheno_OP",dx2)
    gcl.curver_desc("Edge_Litho_OP",dx2)
    gcl.curver_desc("Edge_Astheno_DP",dx2)
    gcl.curver_desc("Edge_Litho_DP",dx2)
    gcl.curver_desc("Contact_OP_DP",dx1)
    gcl.curver_start_end_desc("Bottom_Litho_OP",dx1,dx2,'l')
    gcl.curver_start_end_desc("Top_Litho_OP",dx1,dx2,'l')
    gcl.curver_start_end_desc("Top_Litho_DP",dx1,dx2,'r')
    gcl.curver_start_end_desc("Astheno_Litho_Contact",dx1,dx2,'l')  
    # fabrication du mesh
    cubit.cmd('mesh surface all')
    cubit.cmd('surface all smooth scheme condition number beta 1.7 cpu 10')
    cubit.cmd('smooth surface all')
#    cubit.cmd('surface 1 size auto factor 5')
    # Fabrication de groupe et de nodeset
    for i,s in enumerate(cubit.get_entities("surface")):
        S = cubit.surface(s)
        cubit.cmd('block ' + str(i+1) + ' surface ' + S.entity_name())
        cubit.cmd('block ' + str(i+1) + ' name "' + S.entity_name() + ' "' )
    gcl.create_group_nodeset_desc(['Astheno_Litho_Contact','Contact_OP_DP','Contact_Prism_DP'],"fault_top",20)     
    gcl.create_group_nodeset_desc(['Top_Litho_OP','Top_Prism','Top_Litho_DP'],"ground_surface",20)     
    gcl.create_group_nodeset_desc(['Bottom_Litho_OP'],"bottom_litho_OP",20)     
    gcl.create_group_nodeset_desc(['Bottom_Litho_DP'],"bottom_litho_DP",20)     
    gcl.create_group_nodeset_desc(['Bottom_Astheno_DP','Bottom_Astheno_OP'],"bottom_astheno",20) 
    gcl.create_group_nodeset_desc(['Edge_Litho_OP'],"edge_litho_OP",20) 
    gcl.create_group_nodeset_desc(['Edge_Litho_DP'],"edge_litho_DP",20) 
    gcl.create_group_nodeset_desc(['Edge_Astheno_OP'],"edge_astheno_OP",20) 
    gcl.create_group_nodeset_desc(['Edge_Astheno_DP'],"edge_astheno_DP",20) 
    gcl.create_group_nodeset_desc(['Front_Litho_DP'],"front_litho_DP",20) 
    gcl.create_group_nodeset_desc(['Contact_Prism_OP'],"contact_prism_OP",20)
    # ecriture fichier final
    cubit.cmd('export mesh "' + output_pathfile + '" dimension 2 overwrite')
    return None

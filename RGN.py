#!/bin/python
import igraph
import subprocess
import os
import pdb_atom_renumber
import pdb_residue_renumber
import pdb_download
import sys
import pdb_splitnmr
import networkx as nx
import matplotlib.pyplot as plt
from operator import itemgetter
import re
import time
import datetime
import k_2_communities
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import *
import imp
def format_pdb(pdb_file):
    os.system("mkdir -p clo_cent bet_cent deg_cent")
    if provide_chain == 'n':
        chain = pdb_file[5]
    pdb=pdb_file[0:4]              # pdb chain reference (lowercase)
    if dist_matrix_restart == 'n' and evt_matrix_restart == 'n':
        os.system("rm -r %s EVT_%s"%(pdb,pdb))
    any_residue_numb_additions = 'False'
    if RESTART != 'y' and provided != 'y':
        os.system("wget http://www.rcsb.org/pdb/files/%s.pdb.gz" %pdb) #download pdb
        os.system("gzip -d %s.pdb.gz" %pdb) #unpack pdb
    elif RESTART == 'y':
        os.system("mv %s/%s.pdb . ;rm -r %s"%(pdb,pdb,pdb))
        os.system("mkdir -p %s" %pdb)
    elif provided == 'y':
        print 'Using pdb in current directory'
    os.system("mkdir -p %s" %pdb)
    f=open('%s.pdb' %pdb)      #open pdb to check wether X-RAY or NMR structure
    c=f.readlines()
    b = []
    for line in c:
        if line[0:4] !='ATOM':
                b.append(line)
        elif line[0:4] == 'ATOM' and provide_add_chain == 'y' or line[0:4] == 'HETA' and provide_add_chain == 'y':
            line4=line[0:21]+'A'+line[22:]
            line1=line4[:16]
            line2=line4[17:]
            line3=line1+space+line2
            b.append(line3)
        elif line[0:4]=='ATOM' and line[16] == 'A':
            line1=line[:16]
            line2=line[17:]
            line3=line1+space+line2
            b.append(line3)
        elif line[0:4]=='ATOM' and line[16]==space:
            b.append(line)
    ssbond = []
    hettype = []
    link = []
    for line in b:
        a=line.split()
        if a[0] == 'SSBOND' and a[3] in chain:
            disulfide_res1=a[3],int(a[4]),a[2]
            disulfide_res2 = a[6],int(a[7]),a[5]
            disulfide_res =disulfide_res1,disulfide_res2
            ssbond.append(disulfide_res)
        if 'EXPDTA' in a[0] :
            check=a[1]
        if 'SUBMITTED' in line and 'CONFORMERS' in line :
                last = int(len(a)) - 1
                if a[last] == 'NULL':
                        num_confs = int(1)
                else:
                        num_confs = int(a[last])
        if a[0] == 'HET' and line[12] in chain:
            hettype.append(line[7:10])
            for comp in het_remove:
                if comp in hettype and comp in het_remove:
                    hettype.remove(comp)
        if a[0] == 'LINK':
            if  a[2] != water and a[6] != water:
                link_res1=line[21],int(line[22:27]),str(line[17:20])
                link_res2 = line[51],int(line[52:56]),str(line[47:50])
                link_res =link_res1,link_res2
                if link_res1[0] and link_res2[0] in chain:
                    link.append(link_res)
    if check == 'SOLUTION' and num_confs > 1 and num_confs != 'NULL':       # If NMR, then separate all structures
        for line in b:
                a=line.split()
                if 'REPRESENTATIVE' in line and 'ENSEMBLE' in line:
                        last = int(len(a)) - 1
                        if a[last] == 'NULL':
                            best_conf = 1
                        else:
                            best_conf = int(a[last])-1

                            nmr_pdb=[]
        g = open('%s.pdb' %pdb,'r')
        nmr_pdb = g.readlines()
        g.close()
        models=pdb_splitnmr.splitNMR(nmr_pdb)
        short_pdb = os.path.split('%s.pdb' %pdb)[-1][:-4]
        for index, model in enumerate(models):
            g = open("%s_%i.pdb" % (short_pdb,index),"w")
            g.writelines(model)
            g.close()
        atom_file=[]
        f.close()
        f=open('%s_%s.pdb' %(pdb,best_conf)) # only analyse the first model
        for line in f :
                a=line.split()
                if a[0]=='ATOM' and line[21] in chain:
                        atom_file.append(line.strip('\n'))
        f.close()
    else :
        atom_file=[]
        for line in b:
                a=line.split()
                if a[0]=='ATOM' and line[21] in chain or len(hettype) > 0 and line[0:4]=='HETA' and line[17:20].strip(' ') in hettype and line[21] in chain:
                        atom_file.append(line.strip('\n'))
    g=open('%s_atom.pdb' %pdb,'w')
    for item in atom_file:
            print>>g, item
    g.close()
    if h_added != 'y':
        os.system("reduce -BUILD -QUIET  %s_atom.pdb > %sH.pdb" %(pdb,pdb)) # add hydrogens
    elif h_added == 'y':
        os.system("mv  %s_atom.pdb  %sH.pdb" %(pdb,pdb))
    f.close()
    f=open('%sH.pdb'%pdb)
    pdbH_file=[]
    for line in f:
            a=line.split()
            b=line.strip('\n')
            if a[0]=='ATOM' and line[21] in chain: # prepare pdb chain
                    pdbH_file.append(b)
            elif len(hettype) > 0 and line[0:4]=='HETA' and line[17:20].strip(' ') in hettype and line[21] in chain:
                line1='ATOM     '
                line2=line[9:]
                line3=line1+line2
                b=line3.strip('\n')
                pdbH_file.append(b)
    pdbH_file1 = []
    for line in pdbH_file:
        if line[16] == 'A' or space:
            pdbH_file1.append(line)
    f.close()
    f=open('%s_H.pdb' %pdb,'w')
    for item in pdbH_file1:
        print >> f, item

    f.close()
    f=open('%s_H.pdb' %pdb,'r')
    pdb_file=[]
    for line in f:
        a = line.strip('\n')
        pdb_file.append(a)
    f.close()
    f=open('%s_renumber.pdb' %pdb,'w')
    res_renum=pdb_atom_renumber.pdbAtomRenumber(pdb_file)
    for line in res_renum:
            print >> f, line
    f.close()
    if dilution_plot == 'y':
            os.system("FIRST %s_renumber.pdb -non  -dil 1 -E 0.0  -hbout -phout -L /u/tcmsf1/asf40/bin/FIRST-\
6.2.1-bin-64-gcc3.4.3-O3" %pdb) # identify noncovalent interaction network
            os.system("mv %s_renumber.ps %s"%(pdb,pdb))
    else:
        os.system("FIRST %s_renumber.pdb -non  -E 0.0  -covout -hbout -phout -L /u/tcmsf1/asf40/bin/FIRST-\
6.2.1-bin-64-gcc3.4.3-O3" %pdb) # identify noncovalent interaction network
###############################################################
##################### REFORMAT  ###############################
###############################################################
    reformat_list = []
    for line in res_renum:
        if line[22]== '1' or '2' or '3' or '4' or '5' or '6' or '7' or '8' or '9':
                line1=line[:22]
                line2=line[22:]
                line3=line1+space+line2
                reformat_list.append(line3)
        else:
                reformat_list.append(line)
    f=open('%s_renumber_reform.pdb' %pdb,'w')
    for item in reformat_list:
            f.write("%s \n" % item)
    f.close()
    return pdb,ssbond, link,chain

def hbond_cutoff(deg_cutoff,bet_cutoff,clo_cutoff):
    for cutoff in (deg_cutoff,bet_cutoff,clo_cutoff):
            f = open('hbonds.out')
            cutoff_file=[]
            for line in f :
                    a=line.split()
                    b=line.strip('\n')
                    if float(a[2]) <= cutoff :     # if the cutoff is less than or equal cutoff then include
                            cutoff_file.append(b)
            cutoff_print=open('hb_cutoff.%s' %cutoff, 'w')
            for line in cutoff_file:
                    cutoff_print.write('%s\n'%line)
            cutoff_print.close()
            f.close()
    if matrix_set == 'y':
        cutoff = mat_cutoff
        f = open('hbonds.out')
        cutoff_file=[]
        for line in f :
                a=line.split()
                b=line.strip('\n')
                if float(a[2]) <= cutoff :     # if the cutoff is less than or equal cutoff then include
                        cutoff_file.append(b)
        cutoff_print=open('hb_cutoff.%s' %cutoff, 'w')
        for line in cutoff_file:
                cutoff_print.write('%s\n'%line)
        cutoff_print.close()
        f.close()
    if quality == 'y':
            clo_cutoff=float('-3.0')
            for cutoff in (float(-0.1),float(-0.25),float(-0.5),float(-0.75),float(-1.0),float(-1.5),float(-2.0),float(-2.5),clo_cutoff, float(-3.5),float(-4.0),float(-4.5),float(-5.0),float(-6.0),float(-7.0)):
                    f = open('hbonds.out')
                    cutoff_file=[]
                    for line in f :
                            a=line.split()
                            b=line.strip('\n')
                            if float(a[2]) <= cutoff :     # if the cutoff is less than or equal cutoff then include
                                    cutoff_file.append(b)
                    cutoff_print=open('hb_cutoff.%s' %cutoff, 'w')
                    for line in cutoff_file:
                            cutoff_print.write('%s\n'%line)
                    cutoff_print.close()
                    f.close()

def hbond_atom2res():
    print 'Representing Hydrogen Bond Atom interactions As Residue Interactions'
    for cutoff in (deg_cutoff,bet_cutoff,clo_cutoff):
        atomlinks=[]
        cutoff_file=open('hb_cutoff.%s' %cutoff)

        for line in cutoff_file:
                    l = line.strip('\n').split()
                    atomlinks.append((l[0],l[1],-float(l[2]))) #first atom, second atom, energy
                    del l
        f=open('%s_renumber_reform.pdb' %pdb)
        atomres = {}      #atomres is a dictionary with the key as the
        for line in f:    #atom number and value as residue number type and chain
            l = line.strip().split()
            check_int = any(c.isalpha() for c in l[5])
            if l[0] == 'ATOM' and check_int == False :
                    atomres[l[1]] = l[4],int(l[5]),l[3]  #create dictionary according to residue
            elif check_int == True :
                    print >> error_file, 'HB E: Invalid INT in',l[5],'\n'
                    any_residue_numb_additions = 'True'
        f.close()
        g=[]
        hb_reslinks = []
        for i in atomlinks:
                    if i[0] in atomres and i[1] in atomres:
                        if atomres[i[0]][2] and atomres[i[1]][2] != water: #LINKS DO NOT INCLUDE WATER WATER INTERACTIONS
                            a=atomres[i[0]],atomres[i[1]],i[2],hb_color
                            hb_reslinks.append(a)
                    else:
                            print >> error_file, 'HB E : One of atomes', i[0],' or ', i[1], 'had error  with check_int variable'
        output_file=open('hb_reslinks.%s'%cutoff,'w')
        for reslink in hb_reslinks:
            a=str(reslink)
            output_file.write('%s\n'%a)
        output_file.close()
    if matrix_set == 'y':
        cutoff = mat_cutoff
        atomlinks=[]
        cutoff_file=open('hb_cutoff.%s' %cutoff)

        for line in cutoff_file:
                    l = line.strip('\n').split()
                    atomlinks.append((l[0],l[1],-float(l[2]))) #first atom, second atom, energy
                    del l
#        if man_atom_int == 'y':
#            new_atom_links = atomlinks + man_atom_add
#            atomlinks = new_atom_links
        f=open('%s_renumber_reform.pdb' %pdb)
        atomres = {}      #atomres is a dictionary with the key as the
        for line in f:    #atom number and value as residue number type and chain
            l = line.strip().split()
            check_int = any(c.isalpha() for c in l[5])
            if l[0] == 'ATOM' and check_int == False :
                    atomres[l[1]] = l[4],int(l[5]),l[3]  #create dictionary according to residue
            elif check_int == True :
                    print >> error_file, 'HB E: Invalid INT in',l[5],'\n'
                    any_residue_numb_additions = 'True'
        f.close()
        g=[]
        hb_reslinks = []
        for i in atomlinks:
                    if i[0] in atomres and i[1] in atomres:
                        if atomres[i[0]][2] and atomres[i[1]][2] != water: #LINKS DO NOT INCLUDE WATER WATER INTERACTIONS
                            a=atomres[i[0]],atomres[i[1]],i[2],hb_color
                            hb_reslinks.append(a)
                    else:
                            print >> error_file, 'HB E : One of atomes', i[0],' or ', i[1], 'had error  with check_int variable'
        output_file=open('hb_reslinks.%s'%cutoff,'w')
        for reslink in hb_reslinks:
            a=str(reslink)
            output_file.write('%s\n'%a)
        output_file.close()

        if quality == 'y':
            for cutoff in (clo_cutoff, float(-0.1),float(-0.25),float(-0.5),float(-0.75),float(-1.0),float(-1.5),float(-2.0),float(-2.5),clo_cutoff,float(-3.5),float(-4.0),float(-4.5),float(-5.0),float(-6.0),float(-7.0))\
 :
                atomlinks=[]
                cutoff_file=open('hb_cutoff.%s' %cutoff)
                for line in cutoff_file:
                            l = line.strip('\n').split()
                            atomlinks.append((l[0],l[1],-float(l[2]))) #first atom, second atom, energy
                            del l
                            f=open('%s_renumber_reform.pdb' %pdb)
                            atomres = {}      #atomres is a dictionary with the key as the
                            for line in f:    #atom number and value as residue number type and chain
                                l = line.strip().split()
                                check_int = any(c.isalpha() for c in l[5])
                                if l[0] == 'ATOM' and check_int == False :
                                        atomres[l[1]] = l[4],int(l[5]),l[3]  #create dictionary according to residue
                                elif check_int == True :
                                        print >> error_file, 'HB E: Invalid INT in',l[5],'\n'
                                        any_residue_numb_additions = 'True'
                            f.close()
                            g=[]
                            for i in atomlinks:
                                hb_reslinks = []
                                for i in atomlinks:
                                        if i[0] in atomres and i[1] in atomres:
                                                a=atomres[i[0]],atomres[i[1]],i[2],hb_color
                                                hb_reslinks.append(a)
                                        else:
                                                print >> error_file, 'HB E : One of atomes', i[0],' or ', i[1]\
, 'had error  with check_int variable'
                output_file=open('hb_reslinks.%s'%cutoff,'w')
                for reslink in hb_reslinks:
                        a=str(reslink)
                        output_file.write('%s\n'%a)
                output_file.close()
    return hb_reslinks

def hphobe_atom2res():
###############################################################
##################### HPHOBE choose  ##########################
###############################################################
    f = open('hphobes.out')
    phobe=[]
    for line in f :
        a=line.split()
        b= a[0],a[1]
        phobe.append(b)
###############################################################
##################### HPHOBE ATOM TO RES ######################
###############################################################
    print 'Representing Hydrophobic Interaction Atoms As Residue Interactions'
    phobeatoms = []
    for l in phobe:
        phobeatoms.append((l[0],l[1]))
    del l
    f=open('%s_renumber_reform.pdb' %pdb)
    atomres = {}
    for line in f:
        l = line.strip().split()
        check_int = any(c.isalpha() for c in l[5])
        if l[0] == 'ATOM' and check_int == False:
            atomres[l[1]] = l[4],int(l[5]),l[3]
        elif check_int == True :
            print >> error_file, 'HPHOBE E : Invalid INT in',l[5],'\n'
    f.close()
    g=[]
    for i in phobeatoms:
        phobe_reslinks = []
        for i in phobeatoms:
                if i[0] in atomres and i[1] in atomres:
                        a=atomres[i[0]],atomres[i[1]],phobe_strength, phobe_color
                        phobe_reslinks.append(a)
                else :
                        print >> error_file, 'HPHOBE E: One of atomes', i[0],' or ', i[1], 'had e\
rror  with check_int variable'
    return phobe_reslinks

def cov_atom2res():
###############################################################
##################### COVALENT LINKS ##########################
###############################################################
    f = open('cov.out')
    cov=[]
    for line in f :
        a=line.split()
        b= a[0],a[1]
        cov.append(b)
###############################################################
##################### HCOV ATOM TO RES ######################
###############################################################
    print 'Representing Covalent Interaction Atoms As Residue Interactions'
    covatoms = []
    for l in cov:
        covatoms.append((l[0],l[1]))
    del l
    f=open('%s_renumber_reform.pdb' %pdb)
    atomres = {}
    for line in f:
        l = line.strip().split()
        check_int = any(c.isalpha() for c in l[5])
        if l[0] == 'ATOM' and check_int == False:
            atomres[l[1]] = l[4],int(l[5]),l[3]
        elif check_int == True :
            print >> error_file, 'HCOV E : Invalid INT in',l[5],'\n'
    f.close()
    g=[]
    for i in covatoms:
        cov_links = []
        for i in covatoms:
                if i[0] in atomres and i[1] in atomres:
                        a=atomres[i[0]],atomres[i[1]],cov_strength, cov_color
                        cov_links.append(a)
                else :
                        print >> error_file, 'COV E: One of atomes', i[0],' or ', i[1], 'had error  with check_int variable'
    if len(ssbond) > 0:
        ssbond_orig = {}
        f=open('%s.pdb'%pdb)
        for line in f:
            a=line.split()
            if line[0:4] == 'ATOM' and line[17:20] == 'CYS' and line[12:16].strip() == 'CA':
                coord = float(line[31:39].strip()),float(line[39:47].strip()),float(line[47:55].strip())
                res = line[21],int(line[22:27]),line[17:20]
                ssbond_orig[res] = coord
        f.close()
        ssbond_renum = {}
        f=open('%s_renumber_reform.pdb'%pdb)
        for line in f:
            a=line.split()
            if a[0] == 'ATOM' and a[3] == 'CYS' and a[2].strip() == 'CA':
                coord = float(line[31:39].strip()),float(line[39:47].strip()),float(line[47:55].strip())
                res = line[21],int(line[22:27]),line[17:20]
                ssbond_renum[coord] = res
        f.close()
        for link1 in ssbond:
            if link1[0] != link1[1] and link1[0][0] in chain and link1[1][0] in chain:
                c = ssbond_renum[ssbond_orig[link1[0]]],ssbond_renum[ssbond_orig[link1[1]]],cov_strength,cov_color
                cov_links.append(c)
    if len(link) > 0:
        for link1 in link:
            c = link1[0],link1[1],cov_strength,cov_color
            cov_links.append(c)    
    cov_links_set = set(cov_links)
    del cov_links
    cov_links = list(cov_links_set)
    return cov_links

def make_network():
    if matrix_set == 'y':
        mat=open('hb_reslinks.%s'%mat_cutoff)
        hb_reslinks_mat=[]
        for line in mat:
            hb_reslinks_mat.append(eval(line.strip('\n')))
        mat.close()
    bet=open('hb_reslinks.%s'%bet_cutoff)
    hb_reslinks_bet=[]
    for line in bet:
            hb_reslinks_bet.append(eval(line.strip('\n')))
    bet.close()
    deg=open('hb_reslinks.%s'%deg_cutoff)
    hb_reslinks_deg=[]
    for line in deg:
            hb_reslinks_deg.append(eval(line.strip('\n')))
    deg.close()
    clo=open('hb_reslinks.%s'%clo_cutoff)
    hb_reslinks_clo=[]
    for line in clo:
            hb_reslinks_clo.append(eval(line.strip('\n')))
    clo.close()
    network_list_deg=cov_links+phobe_reslinks+hb_reslinks_deg
    network_list_bet=cov_links+phobe_reslinks+hb_reslinks_bet
    if no_phobe == 'y':
        if matrix_set == 'y':
            network_list_mat = cov_links+hb_reslinks_mat
        network_list_clo = cov_links+hb_reslinks_clo
    elif no_phobe != 'y':
        if matrix_set == 'y':
            network_list_mat = cov_links+phobe_reslinks+hb_reslinks_mat
        network_list_clo = cov_links+phobe_reslinks+hb_reslinks_clo
    print "Running Network Analysis at cutoff of -3.0 kcal/mol (unweighted)"
    print "NB: Clo = -3.0, Bet = -2.5 (WEIGHTED), Deg= -2.0"
    g = nx.Graph()
    for l in network_list_clo:
        g.add_edge(l[0],l[1],{'w': 1,'c':l[3]})
    b = nx.Graph()
    for l in network_list_bet:
        b.add_edge(l[0],l[1],{'w': l[2],'c':l[3]})
    d = nx.Graph()
    for l in network_list_deg:
        d.add_edge(l[0],l[1],{'w': l[2],'c':l[3]})
    m = nx.Graph()
    if matrix_set == 'y':
        for l in network_list_mat:
            m.add_edge(l[0],l[1],{'w': l[2],'c':l[3]})
        mod_dict = []
        res_dict = {}
        res_dict_back={}
        count=0
        for line in sorted(g.nodes()):
            res_dict[line]=count
            res_dict_back[count]=line
            count=count+1
        for line in network_list_mat:
            a=res_dict[line[0]]
            x=res_dict[line[1]]
            add=a,x,line[2]
            mod_dict.append(add)
        num = len(res_dict)
    if matrix_set == 'y':
        return b,d,m,network_list_mat,g
    else:
        return b,d,m,g
def modular(mod):
        color_list =['red','blue','green','cyan','pink','orange','grey','yellow','white','black','purple','magenta',
                     'red','blue','green','cyan','pink','orange','grey','yellow','white','black','purple','magenta']
        print 'Identifying modules in RGN'
    ################## COMMUNITED EDGE BETWEENESSS ###############
        dendrogram_mod = igraph.Graph.community_edge_betweenness(mod, clusters=None, directed=False, weights=True)
        clusters = dendrogram_mod.as_clustering()
        membership = clusters.membership
        import csv
        index_list = []
        count=0
        for i in range(num):
            index_list.append(mod.vs()[count].index)
            count=count+1
        mod_list = []
        for index, membership in zip(index_list, clusters.membership):
            a=index,membership
            mod_list.append(a)
        res_mod_list = []
        for ind,mem in mod_list:
            a=res_dict_back[ind],mem
            res_mod_list.append(a)
        mod_list_sorted = sorted(res_mod_list, key=lambda res_mod_list : res_mod_list[1])
        partition = {}
        for line in mod_list_sorted:
            partition[line[0]] = line[1]
        for com in set(partition.values()) :
            list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]
            p=open('module_pymol_%s.pml'%com,'w')
            print >> p, 'select mod_',com,','
            for line in list_nodes:
                print >> p,'chain', line[0],' and resi',line[1]
            p.close()
            os.system("sed -i 's/_ /_/g' module_pymol_%s.pml "%com)
        a=open("syn.x",'w')
        print >> a, "for x in $(ls module_pymol_*)\ndo\nsed -i ':a;N;$!ba;s/\\n/ + /g' $x\ndone \ncat module* > pymol_module.pml ; rm module*;sed -i 's/, + / , /g' pymol_module.pml;"
        t=open('tmp','w')
        print >>t, 'from pymol import cmd\nfrom pymol.cgo import *\nload %s.pdb\nbg_color white\nhide lines, %s\nshow cartoon, %s'%(pdb,pdb,pdb)
        t.close()
        a.close()
        layout = mod.layout("kk")
        igraph.plot(mod, "modular.png", layout=layout, vertex_color=[color_list[x] for x in clusters.membership], vertex_size=5)
        os.system("chmod +x syn.x;./syn.x")
        os.system("mv pymol_module.pml tmp1; cat tmp tmp1 > pymol_module.pml")
        os.system("rm syn.x tmp*;mkdir %s/Module;cp %s.pdb %s/Module ;mv pymol_module.pml %s/Module;mv modular.png %s/Module"%(pdb,pdb,pdb,pdb,pdb))

        if network_x_commun == 'y':
            partition = community.best_partition(m)
            mod=community.modularity(partition,m)
            for com in set(partition.values()) :
                list_nodes = [nodes for nodes in partition.keys() if partition[nodes] == com]
                p=open('module_pymol_%s.pml'%com,'w')
                print >> p, 'select mod_',com,','
                for line in list_nodes:
                    print >> p,'chain', line[0],' and resi',line[1]
                p.close()
                os.system("sed -i 's/_ /_/g' module_pymol_%s.pml "%com)
            a=open("syn.x",'w')
            print >> a, "for x in $(ls module_pymol_*)\ndo\nsed -i ':a;N;$!ba;s/\\n/ + /g' $x\ndone \ncat module* > pymol_module.pml ; rm module*;sed -i 's/, + / , /g' pymol_module.pml;"
            t=open('tmp','w')
            print >>t, 'from pymol import cmd\nfrom pymol.cgo import *\nload %s.pdb\nbg_color white\nhide lines, %s\nshow cartoon, %s'%(pdb,pdb,pdb)
            t.close()
            a.close()
            os.system("chmod +x syn.x;./syn.x")
            os.system("mv pymol_module.pml tmp1; cat tmp tmp1 > pymol_module.pml")
            os.system("rm syn.x tmp*;mkdir %s/Module;cp %s.pdb %s/Module ;mv pymol_module.pml %s/Module"%(pdb,pdb,pdb,pdb))

def interaction_matrix(b,d,m,network_list_mat):
        print 'Initiating EVT analysis Matrices'
        if dist_matrix_restart == 'y':
            os.system("mv EVT_%s/matrix_dist_bet_full.txt ."%pdb )
        if evt_matrix_restart == 'y':
            os.system("mv EVT_%s/full_EVT_matrix.txt ."%pdb)
        os.system("rm -r EVT_%s"%pdb)
        os.system("mkdir -p EVT_%s" %pdb)
        nodes_list = sorted(b.nodes())
        n=len(nodes_list)
        matrix = [[0 for x in range(n)] for x in range(n)]
        dict1 = {}
        dict2 = {}
        count = 0
        nodes_dict = {}
        nodes_dict1 = {}
        nodes_dict_node = {}
        print 'Building Interaction Matrix'
        for i in range(n):
            label = nodes_list[i][0]+str(nodes_list[i][1])
            nodes_dict1[label]=i
            nodes_dict[i]=label
            row = i+1
            nodes_dict_node[label] = row
        for line in network_list_mat:
            dict2[count] = line[:]
            count=count+1
        len_int = len(network_list_mat)
        for i in range(len_int):
            label1 = dict2[i][0][0]+str(dict2[i][0][1])
            label2 = dict2[i][1][0]+str(dict2[i][1][1])
            dict1[i]=label1,label2,dict2[i][2]
        for i in dict1.keys():
            col1=nodes_dict1[dict1[i][0]]
            col2=nodes_dict1[dict1[i][1]]
            matrix[col1][col2] = dict1[i][2]
            matrix[col2][col1] = dict1[i][2]
        np.savetxt('EVT_%s/matrix_weighted_bet.txt' %pdb, matrix, delimiter=' ')
        del i
        if dist_matrix_restart == 'y':
            print 'Loading EVT and Distance matrices'
            os.system("mv matrix_dist_bet_full.txt EVT_%s"%pdb)
            dist_matrix_full = np.loadtxt("EVT_%s/matrix_dist_bet_full.txt"%pdb)
        elif dist_matrix_restart != 'y':
            print 'Generating Distance Matrix'
            dist_matrix_full = [[0 for x in range(n)] for x in range(n)]
            for node1 in nodes_list:
                for node2 in nodes_list:
                    dist=nx.shortest_path_length(m,source=node1,target=node2, weight='w')
                    col1=nodes_dict1[node1[0]+str(node1[1])]
                    col2=nodes_dict1[node2[0]+str(node2[1])]
                    dist_matrix_full[col1][col2] = dist
            np.savetxt('EVT_%s/matrix_dist_bet_full.txt' %pdb, dist_matrix_full, delimiter=' ')
        num_relation = []
        for i in range(n):
            relate = "residue number is:",nodes_list[i],"matrix column is:",nodes_dict1[nodes_list[i][0]+str(nodes_list[i][1])], nodes_list[i][0]+str(nodes_list[i][1])
            num_relation.append(relate)
        f=open('EVT_%s/matrix_residue_key.txt' %pdb,'w')
        for item in num_relation:
            print>>f, item
        f.close()
        return n, nodes_list, nodes_dict, nodes_dict1, dist_matrix_full

def EVT(n,nodes_list,nodes_dict,nodes_dict1,dist_matrix_full):
    ###############################################################
    ##################### EVT Analysis with Mathematica ###########
    ###############################################################
        if evt_matrix_restart == 'y':
            print 'Running EVT analysis Matrix'
            os.system("mv full_EVT_matrix.txt EVT_%s"%pdb)
            f=open('EVT_%s/EVT.nb'%pdb,'w')
            print >>f,'Import["%s/EVT_%s/full_EVT_matrix.txt", mymatrix ,"Table"];'%(work_dir,pdb)
            print >>f, 'Do[EVTprof = mymatrix[[n]];Export["%s/EVT_%s/EVT_raw_data_"<>ToString[n]<>".txt", EVTprof,"Text"],{n,%s}]'%(work_dir,pdb,n)
            f.close()
        elif evt_matrix_restart == 'n':
            f=open('EVT_%s/EVT.nb'%pdb,'w')
            print >>f, '#!/usr/local/shared/bin/MathematicaScript'
            print >>f, 'affinitymatrix =   Import["%s/EVT_%s/matrix_weighted_bet.txt","Table"];'%(work_dir,pdb)
            print >>f, 'd = Table[Sum[affinitymatrix[[i, j]], {j, 1, Length[affinitymatrix]}], {i, 1, Length[affinitymatrix]}];'
            print >>f, 'T = Table[If[d[[i]] == 0, 0.0, affinitymatrix[[i, j]]/d[[i]]], {i, 1,Length[affinitymatrix]}, {j, 1, Length[affinitymatrix]}];'
            print >>f, 'RedT = Table[Drop[T, {k}, {k}], {k, 1, Length[T]}];'
            print >>f, 'F[k_, Nmax_] := IdentityMatrix[Length[T] - 1] + Sum[MatrixPower[RedT[[k]], l], {l, 1, Nmax}];'
            print >>f, 'row = ConstantArray[0, %s];' %(n-1)
            print >>f,'column = ConstantArray[1, %s];' %n
            print >>f,'F1[k_, Nmax_] := Insert[F[k, Nmax], row, k];'
            print >>f,'F2[k_, Nmax_] := Transpose[Insert[Transpose[F1[k, Nmax]], column, k]];'
            print >>f, 'MEVT[Nmax_] := ParallelSum[F2[k, Nmax], {k, 1, Length[T]}]/Length[T];'
            print >>f, 'max = %d;' %Nmax
            print >>f,'mymatrix=MEVT[max];'
            print >>f,'Export["%s/EVT_%s/full_EVT_matrix.txt", mymatrix ,"Table"];'%(work_dir,pdb)
            print >>f, 'Do[EVTprof = mymatrix[[n]];Export["%s/EVT_%s/EVT_raw_data_"<>ToString[n]<>".txt", EVTprof,"Text"],{n,%s}]'%(work_dir,pdb,n)
            f.close()
        print 'Running Estimated Visiting Time Analysis '
        os.system("math -script EVT_%s/EVT.nb"%pdb)
        print 'EVT Analysis Complete'
        os.system("sed -i 's/*^/e/g' EVT_%s/EVT_raw_data*"%pdb)
        for which in range(n):
            filename = 'EVT_%s/EVT_raw_data_%s.txt'%(pdb,which+1)
            label = nodes_list[which][0]+str(nodes_list[which][1])
            os.system('mkdir EVT_%s/%s'%(pdb,label))
            rename = 'EVT_%s/%s/EVT_raw_data_%s.txt'%(pdb,label,label)
            os.rename(filename, rename)
    #######################################################################
    ##################### EVT ANALYSIS ####################################
    #######################################################################
        def getTopDictionaryKeys(dictionary,number):
            topList = []
            a = dict(dictionary)
            for i in range(number):
                m = max(a, key=a.get)
                topList.append([i,m,a[m]])
                del a[m]
            return topList

        full_evt_matrix = np.loadtxt("EVT_%s/full_EVT_matrix.txt"%pdb)
        ####### IMPROVED FIT #######
        full_evt_matrix_fit = [[0 for x in range(n)] for x in range(n)]
        from scipy.optimize import curve_fit
        def line(x, a, b, c):
            return  a * np.exp(-b * x) + c
        for node in range(n):
            x1=[]
            y1=[]
            for node1 in range(n):
                x1 = dist_matrix_full[node][node1]
                y1 = full_evt_matrix[node][node1]
                value = x1*y1
                full_evt_matrix_fit[node][node1] = value
        #######################
        full_evt_matrix_scaled = full_evt_matrix_fit
        for i in range(n):
            label = nodes_list[i][0]+str(nodes_list[i][1])
            np.savetxt('EVT_%s/%s/EVT_scaled_data_%s.txt' %(pdb,label,label),full_evt_matrix_scaled[i],fmt='%.14f')
            f=open('EVT_%s/%s/EVT_scaled_data_%s.txt'%(pdb,label,label))
            count=1
            tmp_dict={}
            for line in f:
                a=line.strip('\n')
                tmp_dict[nodes_dict[count-1]] = float(a)
                count=count+1
            z=open('EVT_%s/%s/EVT_scaled_ordered.txt'%(pdb,label),'w')
            for line in getTopDictionaryKeys(tmp_dict,n):
                z.write('%s\n'%line)
            z.close()
            f.close()

            z=open('EVT_%s/%s/EVT_raw_ordered.txt'%(pdb,label),'w')
            f=open('EVT_%s/%s/EVT_raw_data_%s.txt'%(pdb,label,label))
            count=1
            for line in f:
                a=line.strip('\n')
                tmp_dict[nodes_dict[count-1]] = float(a)
                count=count+1
            for line in getTopDictionaryKeys(tmp_dict,n):
                z.write('%s\n'%line)
            z.close()
            f.close()
            f=open('EVT_%s/%s/EVT_scaled_data_%s.txt'%(pdb,label,label))
            count=0
            dict_scaled = {}
            scaled_data = []
            for line in f:
                dict_scaled[count]=float(line.strip('{}\n'))
                data_scaled = count+1,float(line.strip('{}\n'))
                scaled_data.append(data_scaled)
                count = count+1
            f.close()
            z=open('EVT_%s/%s/EVT_scaled.data'%(pdb,label),'w')
            for item in scaled_data:
                value = ' '.join(map(str, item))
                z.write("%s\n"%value)
            z.close()
            pdb_scaled = []
            plot = open('EVT_%s/%s/plot.gp' %(pdb,label),'w')
            plot_input = 'set title "EVT Distance Scaled Data";set mxtics 10; set grid ytics xtics mxtics;unset key; set xlabel "Sequence Number";set ylabel "Distance Scaledalised EVT";set terminal postscript eps enhanced color font "Helvetica,25";set output "EVT_%s/%s/%s_EVT_scaled.eps" ; plot "EVT_%s/%s/EVT_scaled.data" using 1:2 with line lt -1' %(pdb,label,label,pdb,label)
            plot.write(plot_input)
            plot.close()
            os.system("gnuplot EVT_%s/%s/plot.gp"%(pdb,label))
            os.system("rm EVT_%s/%s/plot.gp"%(pdb,label))

            f=open('EVT_%s/%s/EVT_raw_data_%s.txt'%(pdb,label,label))
            count=0
            dict_raw = {}
            raw_data = []
            for line in f:
                   dict_raw[count]=float(line.strip('{}\n'))
                   data_raw = count+1,float(line.strip('{}\n'))
                   raw_data.append(data_raw)
                   count = count+1
            f.close()

            pdb_raw = []
            z=open('EVT_%s/%s/EVT.data'%(pdb,label),'w')
            for item in raw_data:
                value = ' '.join(map(str, item))
                z.write("%s\n"%value)
            z.close()
            pdb_scaled = []
            plot = open('EVT_%s/%s/plot.gp' %(pdb,label),'w')
            plot_input = 'set title "EVT Raw Data";set mxtics 10; set grid xtics ytics mxtics; unset key; set xlabel "Sequence Number";set ylabel "Raw  EVT";set terminal \
postscript eps enhanced color font "Helvetica,25";set output "EVT_%s/%s/%s_EVT_raw.eps" ; plot "EVT_%s/%s/EVT.data" using 1:2 with line lt -1 ' %(pdb,label,label,pdb,label)
            plot.write(plot_input)
            plot.close()
            os.system("gnuplot EVT_%s/%s/plot.gp"%(pdb,label))
            os.system("rm EVT_%s/%s/plot.gp"%(pdb,label))
            from sympy import N
            f=open("%s_renumber.pdb"%pdb)
            for line in f:
                a=line[0:78]
                res=int(a[22:27])
                letter = a[21]
                if letter+str(res) not in nodes_dict1.keys():
                    print '%s%s NOT IN GRAPH'%(letter,str(res))
                else:
                    res_add = dict_scaled[nodes_dict1[letter+str(res)]]
                    res_add1 = format(res_add,'.3f')
                    c=a[:60]+res_add1+' '+a[76:78]
                    pdb_scaled.append(c)

                    res_raw=letter+str(int(a[22:27]))
                    res_add_raw = dict_raw[nodes_dict1[res_raw]]
                    res_add_raw1 = format(res_add_raw,'.3f')
                    raw_line=a[:60]+res_add_raw1+' '+a[76:78]
                    pdb_raw.append(raw_line)
            f.close()
            max_evt_raw=max(dict_raw.values())
            min_evt_raw=min(dict_raw.values())
            middle_evt_raw = max_evt_raw - ((abs(max_evt_raw) + abs(min_evt_raw))/2)
            max_evt=max(dict_scaled.values())
            min_evt=min(dict_scaled.values())
            middle_evt = max_evt - ((abs(max_evt) + abs(min_evt))/2)
            a=open('EVT_%s/%s/%s_scaled.pdb'%(pdb,label,pdb),'w')
            for line in pdb_scaled:
                a.write("%s\n" %line)
            a.close()
            a=open('EVT_%s/%s/%s_raw.pdb'%(pdb,label,pdb),'w')
            for line in pdb_raw:
                a.write("%s\n" %line)
            a.close()
            f=open('EVT_%s/%s/scaled_putty.pml'%(pdb,label),'w')
            print >>f, 'from pymol import cmd'
            print >>f, 'from pymol.cgo import *'
            print >>f, 'load %s_scaled.pdb'%pdb
            print >>f, 'bg_color white'
            print >>f, 'spectrum b, blue_white_red, minimum=%s, maximum=%s'%(min_evt,max_evt)
            print >>f, 'as cartoon'
            print >>f, 'cartoon putty'
            print >>f, 'ramp_new colorbar, none, [%s,%s, %s], [blue,white,red]'%(min_evt,middle_evt,max_evt)
            print >>f, 'sele initiation_site, chain %s and resi %s'%(label[0],int(label[1:]))
            print >>f, 'color yellow, initiation_site; center initiation_site'
            print >>f, ' show sticks, initiation_site'
            f.close()
            f=open('EVT_%s/%s/raw_putty.pml'%(pdb,label),'w')
            print >>f, 'from pymol import cmd'
            print >>f, 'from pymol.cgo import *'
            print >>f, 'load %s_raw.pdb'%pdb
            print >>f, 'bg_color white'
            print >>f, 'spectrum b, blue_white_red, minimum=%s, maximum=%s'%(min_evt_raw,max_evt_raw)
            print >>f, 'as cartoon'
            print >>f, 'cartoon putty'
            print >>f, 'ramp_new colorbar, none, [%s,%s, %s], [blue,white,red]'%(min_evt_raw,middle_evt_raw,max_evt_raw)
            print >>f, 'sele initiation_site, chain %s and resi %s'%(label[0],int(label[1:]))
            print >>f, 'color yellow, initiation_site; center initiation_site'
            print >>f, ' show sticks, initiation_site'
            f.close()
            f=open('EVT_%s/%s/to_signal_raw.txt'%(pdb,label),'w')
            to_signal = {}
            values = []
            dists = []
            for node_to in range(n):
                value = full_evt_matrix[node_to][nodes_dict1[label]]
                dist = dist_matrix_full[node_to][nodes_dict1[label]]
                dists.append(dist)
                values.append(value)
            def line(x, a, b, c):
                return  a * np.exp(-b * x) + c
            x = np.asarray(dists)
            y = np.asarray(values)
            popt, pcov = curve_fit(line, x, y, maxfev = 2000)
            y_new = []
            for num in range(len(y)):
                y_fit = y[num] - (line(x[num], *popt))
                y_new.append(y_fit)
            y_fit = np.asarray(y_new)

            for node in range(n):
                to_signal[nodes_dict[node]]=y_fit[node]
                f.write('%s\n'%y_fit[node])
            f.close()

            to_signal_pdb = []
            f=open("%s_renumber.pdb"%pdb)
            for line in f:
                a=line[0:78]
                res=int(a[22:27])
                letter = a[21]
                if letter+str(res) not in nodes_dict1.keys():
                    print '%s%s NOT IN GRAPH'%(letter,str(res))
                else:
                    res_add = to_signal[letter+str(res)]
                    res_add1 = format(res_add,'.3f')
                    c=a[:60]+res_add1+' '+a[76:78]
                    to_signal_pdb.append(c)
            f.close()
            max_evt=max(to_signal.values())
            min_evt=min(to_signal.values())
            middle_evt = max_evt - ((abs(max_evt) + abs(min_evt))/2)
            a=open('EVT_%s/%s/%s_signal_to_fit.pdb'%(pdb,label,label),'w')
            for line in to_signal_pdb:
                a.write("%s\n" %line)
            a.close()

            f=open('EVT_%s/%s/%s_signal_to_fit.pml'%(pdb,label,label),'w')
            print >>f, 'from pymol import cmd'
            print >>f, 'from pymol.cgo import *'
            print >>f, 'load %s_signal_to_fit.pdb'%label
            print >>f, 'bg_color white'
            print >>f, 'spectrum b, blue_white_red, minimum=%s, maximum=%s'%(min_evt,max_evt)
            print >>f, 'as cartoon'
            print >>f, 'cartoon putty'
            print >>f, 'ramp_new colorbar, none, [%s,%s, %s], [blue,white,red]'%(min_evt,middle_evt,max_evt)
            print >>f, 'sele initiation_site, chain %s and resi %s'%(label[0],int(label[1:]))
            print >>f, 'color yellow, initiation_site; center initiation_site'
            print >>f, ' show sticks, initiation_site'
            f.close()

            z=open('EVT_%s/%s/EVT_to_signal_ordered.txt'%(pdb,label),'w')
            for line in getTopDictionaryKeys(to_signal,n):
                z.write('%s\n'%line)
            z.close()
            os.system("rm EVT_%s/%s/EVT_raw_data_%s.txt EVT_%s/%s/EVT_scaled_data_%s.txt;mkdir EVT_%s/%s/pdbfiles;mv EVT_%s/%s/*.pdb EVT_%s/%s/*.pml EVT_%s/%s/pdbfiles"%(pdb,label,label,pdb,label,label,pdb,label,pdb,label,pdb,label,pdb,label))
            return full_evt_matrix,full_evt_matrix_scaled

def full_evt(full_evt_matrix,full_evt_matrix_scaled):
        if EVT_centrality == 'y':
            import igraph
            weight_list = []
            for node in range(n):
                for node2 in range(n):
                    weight = full_evt_matrix_scaled[node][node2]
                    if weight != 0:
                        weight_list.append(weight)
            arr = np.asarray(weight_list)
            EVT_mean = np.mean(arr, axis=0)
            EVT_s = np.std(arr, axis=0)

            ulimit = 0
            m=nx.DiGraph()
            mod = igraph.Graph(n,directed=True)

            for node in range(n):
                for node2 in range(n):
                    vert1 = nodes_list[node]
                    vert2 = nodes_list[node2]
                    weight = full_evt_matrix_scaled[node][node2]
                    if weight > ulimit:
                        m.add_edge(vert1,vert2,{'w' : weight})
                        mod.add_edge(nodes_dict1[nodes_list[node][0]+str(nodes_list[node][1])],nodes_dict1[nodes_list[node2][0]+str(nodes_list[node2][1])], weights = weight)
            color_list =['red','blue','green','cyan','pink','orange','grey','yellow','white','black','purple','magenta',
                         'red','blue','green','cyan','pink','orange','grey','yellow','white','black','purple','magenta']

#           dir2undir():
            def getTopDictionaryKeys(dictionary,number):
                topList = []
                a = dict(dictionary)
                for i in range(number):
                    m = max(a, key=a.get)
                    topList.append([i,m,a[m]])
                    del a[m]
                return topList


            weight_total = 0
            weight_list = []
            UN = nx.Graph()
            for node in nodes_list:
                UN.add_node(node)
            for node1 in nodes_list:
                for node2 in nodes_list:
                    if node1 != node2:
                        weight1 = m.get_edge_data(node1,node2)['w']
                        weight2 = m.get_edge_data(node2,node1)['w']
                        weight3 = weight1+weight2
                        UN.add_edge(node1,node2,{'w': weight3})
                        weight_list.append(weight1)
                        weight_list.append(weight2)
            arr = np.asarray(weight_list)
            dir_mean = np.mean(arr, axis=0)
            dir_s = np.std(arr, axis=0)
            UN = nx.Graph()
            for node in nodes_list:
                    UN.add_node(node)
            for node1 in nodes_list:
                for node2 in nodes_list:
                    if node1 != node2:
                        weight1 = m.get_edge_data(node1,node2)['w']
                        weight2 = m.get_edge_data(node2,node1)['w']
                        weight3 = (weight1+weight2)/2
                        if weight3 > dir_mean+dir_s:
                            UN.add_edge(node1,node2,{'w': weight3})
            bet_cen_un = nx.betweenness_centrality(UN,weight = 'w')
            top_bet_cen_un = getTopDictionaryKeys(bet_cen_un,n)
def EVT_central():
            def getTopDictionaryKeys(dictionary,number):
                topList = []
                a = dict(dictionary)
                for i in range(number):
                    m = max(a, key=a.get)
                    topList.append([i,m,a[m]])
                    del a[m]
                return topList

            bet_cen = nx.betweenness_centrality(m,weight = 'w')
            top_bet_cent_i = getTopDictionaryKeys(bet_cen,n)
            a = open('EVT_%s/%s_bet_cent.out' %(pdb,pdb),'w')
            for line in top_bet_cent_i:
                a.write("%s\n" %line)
            a.close()

def network_analysis(g):
###############################################################
##################### NETWORK ANALYSIS ########################
###############################################################

    f = open('analysis_%s.out' %pdb,  'w')
    N,K = g.order(), g.size()
    avg_deg = float(K)/N
    avg_clust = nx.average_clustering(g)
    largest_k_1_clique = list(nx.k_clique_communities(g, 3))
    size_largest_k_1_clique = len(largest_k_1_clique[0])
    if no_phobe != 'y':
        top_two = len(largest_k_1_clique[0]) + len(largest_k_1_clique[1])
        top_three = len(largest_k_1_clique[0]) + len(largest_k_1_clique[1]) + len(largest_k_1_clique[2])
    size_largest_k_2_clique = list(k_2_communities.k_clique_communities_2(g,3))
    top_k_2 = len(size_largest_k_2_clique[0])
    s = []
    for line in size_largest_k_2_clique[0]:
            s.append(line)
    sub = nx.subgraph(g,s) # this creates the subgraph of g with the nodes in s.
    avg_clust_k_2 = nx.average_clustering(sub)
    largest_clique=nx.graph_clique_number(g)
    cliques_great_2 = []
    for clique in nx.find_cliques(g):
            if len(clique)>2:
                    cliques_great_2.append(clique)
    print >> f, 'Nodes:', N
    print >> f, 'Interactions:', K
    print >> f, 'Average Degree:', avg_deg
    print >> f, 'The largest k-1 clique:', size_largest_k_1_clique
    if no_phobe != 'y':
        print >> f, 'The sum of the top two largest (k-1) cliques:', top_two
        print >> f, 'The sum of the top three largest (k-1) cliques:', top_three
    print >> f, 'The size of the largest k-2 community is:',top_k_2
    print >> f, 'The average clustering coefficient is', avg_clust
    print >> f, 'The average clustering coefficient of the largest k-2 community is', avg_clust_k_2
    print >> f, 'The number of cliques greater than two are',len(cliques_great_2)
    print >> f,'\n'
    print >> f, 'The k(=3)-1-clique communities are', largest_k_1_clique
    print >>f, 'The cliques greater than two are', cliques_great_2
    print >>f, 'The degree historgram is', nx.degree_histogram(g)
###############################################################
##################### CENTRALITY###############################
###############################################################
    bet_cen = nx.betweenness_centrality(b, weight = 'w')
    clo_cen = nx.closeness_centrality(g)
    deg_cen = nx.degree_centrality(d)
############## top keys from a python dictionary ##############
    def getTopDictionaryKeys(dictionary,number):
        topList = []
        a = dict(dictionary)
        for i in range(number):
            m = max(a, key=a.get)
            topList.append([i,m,a[m]])
            del a[m]
        return topList
############## IDENTIFY TOP 30 and MAX CENTRALITIES ############
    bet_n = len(bet_cen)
    clo_n = len(clo_cen)
    deg_n = len(deg_cen)

    top_bet_cen_i = getTopDictionaryKeys(bet_cen,30)
    a = open('bet_cent/%s_bet_cent_%s.out' %(pdb, 30),'w')
    for line in top_bet_cen_i:
            a.write("%s\n" %line)
    a.close()
    top_bet_cen_i = getTopDictionaryKeys(bet_cen,bet_n)
    a = open('bet_cent/%s_bet_cent_%s.out' %(pdb, bet_n),'w')
    for line in top_bet_cen_i:
            a.write("%s\n" %line)
    a.close()

    top_clo_cen_i = getTopDictionaryKeys(clo_cen,30)
    a = open('clo_cent/%s_clo_cent_%s.out' %(pdb, 30),'w')
    for line in top_clo_cen_i:
            a.write("%s\n" %line)
    a.close()
    top_clo_cen_i = getTopDictionaryKeys(clo_cen,clo_n)
    a = open('clo_cent/%s_clo_cent_%s.out' %(pdb, clo_n),'w')
    for line in top_clo_cen_i:
            a.write("%s\n" %line)
    a.close()
    top_deg_cen_i = getTopDictionaryKeys(deg_cen,30)
    a = open('deg_cent/%s_deg_cent_%s.out' %(pdb, 30),'w')
    for line in top_deg_cen_i:
        a.write("%s\n" %line)
    a.close()
    top_deg_cen_i = getTopDictionaryKeys(deg_cen,deg_n)
    a = open('deg_cent/%s_deg_cent_%s.out' %(pdb, deg_n),'w')
    for line in top_deg_cen_i:
        a.write("%s\n" %line)
    a.close()

def quality_func():
            f = open('transition_profiles_%s.out' %pdb,  'w')
            print >> f, 'The Transition Profiles is as follows'
            trans_graph = open('%s/%s.trans_profile.largest_comp' %(pdb,pdb),'w') # TRANSITION OF LARGEST CONNECTED COMPONENT
            trans_graph_all = open('%s/%s.trans_profile.tot_comp' %(pdb,pdb),'w')    # TOTAL NUMBER OF NODES IN ALL CONNECTED COMPONENTS
            trans_K2_commun = open('%s/%s.trans_profile.K2' %(pdb,pdb),'w') #transition profile of K1 community
            trans_graph_degree = open('%s/%s.trans_profile.degree' %(pdb,pdb),'w') #TRANSITION OF TOTAL NUMBER OF EDGES
            trans_k1_commun = open('%s/%s.trans_profile.K1' %(pdb,pdb),'w') # transition profile of k2 community
            for HBD_cut in (float(-0.1),float(-0.25),float(-0.5),float(-0.75),float(-1.0),float(-1.5),float(-2.0),float(-2.5),clo_cutoff,float(-3.5),float(-4.0),float(-4.5),float(-5.0),float(-6.0),float(-7.0)) :
                    trans_set=open('hb_reslinks.%s' %HBD_cut)
                    hb_reslinks=[]
                    for line in trans_set:
                            hb_reslinks.append(eval(line.strip('\n')))
                    network_list_profile=phobe_reslinks+hb_reslinks
                    T = nx.Graph()
                    for l in network_list_profile:
                            T.add_edge(l[0],l[1],{'w': l[2]})
                    deg = T.size()
                    print >> trans_graph_degree, HBD_cut, deg
                    graphs=list(nx.connected_component_subgraphs(T))
                    lrg_comp = float(len(list(graphs[0])))
                    N = len(T.nodes())
                    prop = lrg_comp / N
                    num_nodes = float(T.order())
                    num_nodes_div_total = num_nodes / N
                    print >> f, 'The largest connected component at cutoff','%s' %HBD_cut,'is',lrg_comp ,'(which is',prop,' of total reisdue\
s)'
                    print >> f, 'The average clustering coefficient of the this component is', nx.average_clustering(T),'\n'
                    print >> trans_graph, HBD_cut, prop #PROPORTION OF NODES IN LARGEST SUBCOMPONENT, that is involved due to NC interaction
                    print >> trans_graph_all, HBD_cut,num_nodes_div_total # PROPORTION OF NODES
                    largest_k_1_clique = list(nx.k_clique_communities(T, 3))
                    size_largest_k_1_clique = len(largest_k_1_clique[0])
                    size_largest_k_2_clique = list(k_2_communities.k_clique_communities_2(T,3))
                    top_k_2 = len(size_largest_k_2_clique[0])
                    print >> trans_k1_commun,HBD_cut,size_largest_k_1_clique
                    print >> trans_K2_commun, HBD_cut,top_k_2
            trans_set.close()
            trans_graph.close()
            trans_graph_all.close()
            trans_k1_commun.close()
            trans_K2_commun.close()
            trans_graph_degree.close()
            f.close()

def shortest_path_func():
            f=open('%s_shortest_path.txt'%pdb)
            for line in f:
                a=line.split()
                if a[0] == 'begin':
                    begin = str(a[1]),a[2],str(a[3])
                if a[0] == 'end':
                    end = str(a[1]),a[2],str(a[3])

            for network in 'unBet','wClo':
                if network == 'unBet':
                    residue_count = []
                    tmp = open('%s.shortest_path_analysis.txt' %network, 'w')
                    z=b
                    shortest_path = nx.shortest_path(b,source=begin,target=end, weight='w')
                    weighted_path_length = nx.shortest_path_length(z,source=begin,target=end,weight = 'w')
                if network == 'wClo':
                    residue_count = []
                    tmp = open('%s.shortest_path_analysis.txt' %network,'w')
                    z=g
                    shortest_path = nx.shortest_path(z,source=begin,target=end,weight = 1)
                    weighted_path_length = nx.shortest_path_length(z,source=begin,target=end,weight = 1)
                num_residues_in_path = len(shortest_path)
                print >> tmp,'\n', 'THE ORIGINAL SHORTEST PATH IS:', shortest_path
                print >> tmp, 'this has weighted length', weighted_path_length
                for line in shortest_path:
                        residue_count.append(line)
                original_path = shortest_path
                count=0
                for residue in original_path :
                        if count < num_residues_in_path-2 :
                            if network == 'unBet' :
                                weight = z.get_edge_data(original_path[count],original_path[count+1])
                            if network =='wClo' :
                                weight = 1
                            z.remove_edge(original_path[count],original_path[count+1])
                            if network == 'unBet' :
                                shortest_path = nx.shortest_path(z,source=begin,target=end, weight='w')
                                weighted_path_length = nx.shortest_path_length(z,source=begin,target=end, weight='w')
                                z.add_edge(original_path[count],original_path[count+1],weight)
                            if network == 'wClo' :
                                shortest_path = nx.shortest_path(z,source=begin,target=end,weight = 1)
                                weighted_path_length = nx.shortest_path_length(z,source=begin,target=end, weight = 1)
                                z.add_edge(original_path[count],original_path[count+1])
                            print >> tmp, 'removed',original_path[count],original_path[count+1],'-->',shortest_path
                            print >> tmp, 'this has weighted length', weighted_path_length
                            z.add_edge(original_path[count],original_path[count+1])
                        for line in shortest_path:
                                residue_count.append(line)
                        count = count+1
                list_sorted = []
                list_sort = []
                for unique in set(residue_count):
                        number = residue_count.count(unique)
                        a = unique, number
                        list_sort.append(a)
                list_sorted = sorted(set(list_sort))
                sorted_appearences = sorted(list_sorted, key=lambda list_sorted : list_sorted[1])
                for line in sorted_appearences:
                        print >> tmp, line
                tmp.close()
                os.system("mv %s.shortest_path_analysis.txt %s" %(network,pdb))

def make_plot():
        import pylab as pl
        import math

        def centrality_scatter(dict1,dict2,path='',
        ylab='',xlab='',title='',line=False):
         # Create figure and drawing axis
         fig = plt.figure(figsize=(200,50))
         ax1 = fig.add_subplot(111)
         # Create items and extract centralities
         items1 = sorted(dict1.items())
         items2 = sorted(dict2.items())
         xdata=[b for a,b in items1]
         ydata=[b for a,b in items2]
         labels = []
         # Add each actor to the plot by ID
         print 'Constructing Closeness - Degree scatter plot with sequence number labels, if chain required modify at ax1.annotate()'
         for p in xrange(len(items1)):
             ax1.annotate(xy=(xdata[p], ydata[p]), s=str(items1[p][0][1]), size = 10, xytext = (100, 60), textcoords = 'offset points', ha = 'right', va = 'bottom',arrowprops = dict(arrowstyle = '->',fc="0.6", connectionstyle = 'arc3,rad=0.3'))
             labels.append(items1[p][0])
             if line:
           # use NumPy to calculate the best fit
                 slope, yint = plt.polyfit(xdata,ydata,1)
                 xline = plt.xticks()[0]
                 yline = map(lambda x: slope*x+yint,xline)
                 ax1.plot(xline,yline,ls='--',color='b')
            # Set new x- and y-axis limits
         plt.xlim((0.0,max(xdata)+(.15*max(xdata))))
         plt.ylim((0.0,max(ydata)+(.15*max(ydata))))
         # Add labels and save
         ax1.set_title(title)
         ax1.set_xlabel(xlab)
         ax1.set_ylabel(ylab)

        centrality_scatter(clo_cen, deg_cen,ylab='Degree Centrality',xlab='Closeness Centrality',title='Residues plotted by Closeness and Degree Centrality')
        plt.savefig('clo_deg_scatter.png')
        os.system("mv clo_deg_scatter.png %s" %pdb)

def subgraph_fun():
        print 'Subgraph According to Closeness Cutoff, change\
        g to b and add weight if you want betweeness'
        s = []
        resnum = {}
        for line in g.nodes():
            resnum[line[1]] = line
        for value in sublist:
                s.append(resnum[value])
        print 'Analysing Subgraph with Nodes:'
        for entry in s:
            print entry
#        k = nx.subgraph(b,s) CHANGE TO THIS IF YOU WANT BETWEENNESS
        k = nx.subgraph(g,s)
        f = open('%s/analysis_%s_SUBGRAPH.out' %(pdb,pdb),  'w')
        N,K = k.order(), k.size()
        n=len(k.nodes())
        avg_deg = float(K)/N
        avg_clust = nx.average_clustering(k)
        largest_k_1_clique = list(nx.k_clique_communities(k, 3))
        size_largest_k_1_clique = len(largest_k_1_clique[0])
        top_two = len(largest_k_1_clique[0]) + len(largest_k_1_clique[1])
        top_three = len(largest_k_1_clique[0]) + len(largest_k_1_clique[1]) + len(largest_k_1_clique[2])
        size_largest_k_2_clique = list(k_2_communities.k_clique_communities_2(k,3))
        largest_clique=nx.graph_clique_number(k)
        cliques_great_2 = []
        for clique in nx.find_cliques(k):
                if len(clique)>2:
                        cliques_great_2.append(clique)
        print >> f, 'Nodes:', N
        print >> f, 'Interactions:', K
        print >> f, 'Average Degree:', avg_deg
        print >> f, 'The sum of the top two largest (k-1) cliques:', top_two
        print >> f, 'The sum of the top three largest (k-1) cliques:', top_three
        print >> f, 'The average clustering coefficient is', avg_clust
        print >> f, 'The number of cliques greater than two are',len(cliques_great_2)
        print >> f,'\n'
        clo_cen = nx.closeness_centrality(g)
    ############## top keys from a python dictionary ##############
        def getTopDictionaryKeys(dictionary,number):
            topList = []
            a = dict(dictionary)
            for i in range(number):
                m = max(a, key=a.get)
                topList.append([i,m,a[m]])
                del a[m]
            return topList
        #Print top 10 centralities
        #NOTE: FOR THIS SCRIPT 100 has been replaced with MAX
        for i in (30,n):
            top_clo_cen_i = getTopDictionaryKeys(clo_cen,i)
            a = open('clo_cent/%s_clo_cent_%s_SUBGRAPH.out' %(pdb, i),'w')
            for line in top_clo_cen_i:
                    a.write("%s\n" %line)
            a.close()
            del line

def network_vis():
        col = []
        os.system("mkdir %s/Vis" %pdb)
        print 'Visualising weighted betweenness centrality network'
        for l in b.edges():
            col.append(b[l[0]][l[1]].values()[0])
        if node_col_opt == 'y':
            z=nx.Graph()
            for line in node_col_list:
                z.add_node(line)
            blue=z.nodes()

        pos=nx.spring_layout(b,k=1.0,iterations=50,weight='w')
        nx.draw_networkx_nodes(b,pos=pos)
        if node_col_opt == 'y':
            nx.draw_networkx_nodes(z,pos=pos,nodelist=blue,node_color=node_alt_col)
        nx.draw_networkx_edges(b,pos,edge_color=col)
        nx.draw_networkx_labels(b,pos,font_size=6,ax=None)
        plt.axis('off')
        plt.savefig('%s/Vis/spring.png' %pdb, dpi=300)
        nx.write_dot(b,'network.dot1')
        dot_input = 'graph  { \n size =  \"7,7\" ;\noverlap = false;\nsplines = true\nfixedsize=False\n node[fontsize = 30];'
        with open('dot_add.txt', 'w') as file1:
            file1.write(dot_input)
        os.system("sed -i 's/graph  {//g' network.dot1")
        os.system("cat dot_add.txt network.dot1 > network.dot ; rm dot_add.txt; rm network.dot1")
        os.system("neato -Tps network.dot -o network.eps")
        os.system("mv network.eps %s/Vis/"%pdb)
        os.system("mv network.dot %s/Vis/"%pdb)

def getVarFromFile(filename):
    import imp
    f = open(filename)
    global data
    import os
    data = imp.load_source('data', '', f)
    f.close()
    return data

def codons_fun():
    def chunkIt(seq, num):      #define function to break centrality list into 20 chunks
        avg = len(seq) / float(num)
        out = []
        last = 0.0
        while last < len(seq):
            out.append(seq[int(last):int(last + avg)])
            last += avg
        return out
    f=open('%s'%codons_path)
    codons_dict={}
    for line in f:                #create dictionary from list of codons
        l=line.strip().split()
        a=l[0][0:4]
        b=l[1]
        codons_dict['%s.%s,' %(a,b)]=line
        del a,b,line,l
    f.close()
    return

def codeml():
        for i in ('clo_cent', 'bet_cent', 'deg_cent'):
                os.chdir(i)
                f=open('%s_%s.out' %(pdb,i))
                b = []
                for line in f:
                     b.append(line.strip('\n'))
                f.close()
                del line

                res_cent_dict = {}
                for cent in b :
                    a=cent.split(' ')
                    res_cent_dict[a[0]] = a[1:]
                del a

                split=chunkIt(range(n),20)
                split_dict={}
                counter=0
                for bin in range(5,105,5):
                    split_dict[bin]=split[counter]
                    counter=counter+1
                    cent_bin = []
                    if len(split[0]) == 1:
                            for z in range(split_dict[bin][0],split_dict[bin][0]+len(split_dict[bin])):
                                    cent_bin.append(res_cent_dict['[%d,' %z])
                    else:
                            for z in range(split_dict[bin][1]-1,split_dict[bin][1]+len(split_dict[bin])-1):
                                    cent_bin.append(res_cent_dict['[%d,' %z])
                    list=[]
                    for line in cent_bin :
                            b=line
                            list.append(b)
                            del line
                    codon_file=open('full.%s.%s' %(pdb,bin),'w')
                    for item in list:
                        codon_file.write("%s\n" % item)
                    codon_file.close()
                del list
                os.chdir('%s' %work_dir)

def clean():
    os.system("cp -r *_cent %s;rm -rf *_cent" %pdb)
    os.system("mv analysis_%s.out %s" %(pdb,pdb))
    os.system("mv error_file.txt %s" %pdb)
    os.system("rm hb_cutoff* hb_reslinks*")
    os.system("mv hbonds.out hphobes.out %s" %pdb)
    os.system("mv %s.pdb %s_renumber.pdb %s"%(pdb,pdb,pdb))
    os.system("mv %s_renumber_results.txt %s/FIRST_results" %(pdb,pdb))
    os.system("rm %s_* %sH.pdb" %(pdb,pdb))

###############################################################
##################### RC IDENTIFICATION #######################
###############################################################
    if RC_identify == 'y':
        os.system("mkdir %s/Rigid_Clusters"%pdb)
        os.system("mv %s/%s_renumber.pdb ."%(pdb,pdb))
        os.system("FIRST %s_renumber.pdb  -E %f -L /u/tcmsf1/asf40/bin/FIRST-6.2.1-bin-64-gcc3.4.3-O3" %(pdb,RC_cutoff)) # identify noncovalent interaction network
        f=open('%s_renumber_RCD.pdb'%pdb)
        RC = {}
        for line in f:
            if line[0:3] != 'END' and line[13:15] == RC_atom:
                a1 = line[60:67]
                b = int(line[22:27])
                a = int(eval(a1))
                RC[b]=a
        RC_sets={}
        for rc in set(RC.values()): # for each RC
            values = []
            for node in RC.keys():     #looking at each node
                if int(RC[node]) == int(rc):
                    values.append(node)
            RC_sets[rc]=values
        for i in RC_sets.keys():
            if len(RC_sets[i]) > 1:
                f=open('%s/Rigid_Clusters/RC_%d'%(pdb,i),'w')
                for item in RC_sets[i]:
                    f.write("%s \n" % item)
                f.close()
        os.system("mv %s_renumber.pdb %s"%(pdb,pdb))
        os.system("rm %s_* "%pdb)
        print 'Finding residue %s atoms belonging to RCs'%RC_atom
        print 'Note: If an atom other than %s is in a RC then the node will'%RC_atom
        print '      not be identified as belonging to that RC'
    sec=time.time()-start
    HMS=str(datetime.timedelta(seconds=sec))
    print '##############################'
    print 'PGN took %.7s H:M:S to run' %HMS
    print '##############################'
    print '##############################'
    print '####CALCULATIONS CORRECT######'
    print '##############################'    


################################################
############ MAIN ##############################
################################################
data = getVarFromFile('pgn.in')
from data import *
if codons == 'y':
    codons = codons_path
    codons()
os.system("mkdir -p clo_cent bet_cent deg_cent")
pdb_files=open('%s'%pdb_list_file)
for pdb_file in pdb_files:
    pdb,ssbond,link,chain = format_pdb(pdb_file)
    hbond_cutoff(deg_cutoff,bet_cutoff,clo_cutoff)
    hb_reslinks = hbond_atom2res()
    phobe_reslinks = hphobe_atom2res()
    cov_links = cov_atom2res()
    if matrix_set == 'y':
        b,d,m,network_list_mat,g = make_network()        
    else:
        b,d,m,g = make_network()
    if matrix_set == 'y':
        n,nodes_list,nodes_dict,nodes_dict1,dist_matrix_full = interaction_matrix(b,d,m,network_list_mat)
#    if evt_matrix_restart == 'y':
        full_evt_matrix,full_evt_matrix_scaled = EVT(n,nodes_list,nodes_dict,nodes_dict1,dist_matrix_full)
    if EVT_centrality == 'y':
        full_evt(full_evt_matrix,full_evt_matrix_scaled)
#    if dir_to_un == 'y':
#        dir2undir()
    if EVT_centrality == 'y':
        EVT_central()
    network_analysis(g)
    if quality == 'y':
        quality_func()
    if shortest_path == 'y':
        shortest_path_func()
    if make_plot == 'y':
        make_plot()
    if subgraph == 'y':
        subgraph_fun()
    if net_vis == 'y':
        network_vis()
    if codons == 'y':
        codons_fun()
        codeml()
    clean()

##!/usr/bin/env python
## -*- coding: utf8 -*-
"""
Sandbox for functions in developpment

NOT FULLY TESTED! DON'T COMPLAIN IF YOUR COMPUTER EXPLODED!

"""
import copy, os, numpy
from NeuroTools.parameters import ParameterSet

def make_name(params_set,range_keys):
    range_keys.sort()
    name = ''
    chars_to_be_removed = ['[',']','/',' ',':','(',')','{','}',',','.',"'"]
    for key in range_keys:
        if key is not None:
            name = name + key + '_' + str(eval('params_set.'+key)) + '_'
    for char in chars_to_be_removed:
        name = name.replace(char,'')
    return name

def check_name(sim_name):
    if os.path.exists(sim_name+'running'):
        return False
    elif not os.path.exists(sim_name+'running'):
        os.system('touch '+sim_name+'running')
        return True

def string_table(tablestring):
    """Convert a table written as a multi-line string into a dict of dicts."""
    tabledict = {}
    rows = tablestring.strip().split('\n')
    column_headers = rows[0].split()
    for row in rows[1:]:
        row = row.split()
        row_header = row[0]
        tabledict[row_header] = {}
        for col_header,item in zip(column_headers[1:],row[1:]):
            tabledict[row_header][col_header] = float(item)
    return tabledict

def string_table_ParameterSet(tablestring):
    """Convert a table written as a multi-line string into a dict of dicts."""
    tabledict = ParameterSet({})
    rows = tablestring.strip().split('\n')
    column_headers = rows[0].split()
    for row in rows[1:]:
        row = row.split()
        row_header = row[0]
        tabledict[row_header] = ParameterSet({})
        for col_header,item in zip(column_headers[1:],row[1:]):
            tabledict[row_header][col_header] = float(item)
    return tabledict

def rUpdate(targetDict, itemDict):
    for key, val in itemDict.items():
        if type(val) == type({}):
            newTarget = targetDict.setdefault(key,{})
            rUpdate(newTarget, val)
        else:
            targetDict[key] = val
# ====================================================================== #
# creating experiments dict.
# ====================================================================== #

def get_experiment_list(params):
    """
    Takes params = dict with all parameters
    Calculates cross product of all and returns a list with all experiments.
    """
    f=lambda ss,row=[],level=0: len(ss)>1 \
       and reduce(lambda x,y:x+y,[f(ss[1:],row+[i],level+1) for i in ss[0]]) \
       or [row+[i] for i in ss[0]]

    tmplist=[]

    names = params.keys()
    for experiment in f(params.values()):
        tmptmpdict = {}
        for name , value in zip(names,experiment):
            tmptmpdict[name]=value
            tmplist.append(tmptmpdict)

    return tmplist

def get_experiment_dict(params):
    """
    Takes params = dict with all parameters
    Calculates cross product of all and returns a dict with all experiments.
    """
    f=lambda ss,row=[],level=0: len(ss)>1 \
       and reduce(lambda x,y:x+y,[f(ss[1:],row+[i],level+1) for i in ss[0]]) \
       or [row+[i] for i in ss[0]]

    count = 0
    tmpdict={}

    names = params.keys()
    for experiment in f(params.values()):
        exp_name = 'exp' + str(count)

        tmptmpdict = {}
        for name , value in zip(names,experiment):

            tmptmpdict[name]=value
            tmpdict[exp_name] = tmptmpdict
        count +=1
    return tmpdict

def make_experiments(parameters,parameters_template, use_name = True):
    experiments_tmp =  get_experiment_dict(parameters)
    experiments = {}
    for i, experiment_tmp in enumerate(experiments_tmp.values()):
        experiment = copy.deepcopy(parameters_template)
        experiment['run'] = copy.deepcopy(experiment_tmp)
        rUpdate(experiment,experiment_tmp)
        experiments[i] = experiment
    return experiments


# ====================================================================== #
def cross(*args):
    """
    Return the cross-product of a variable number of lists (e.g. of a list
    of lists).

    print cross(s1,s2,s3)
    OBSOLETE / LESS EFFICIENT than get_experiment_dict
    From:
    http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/159975
    """

    ans = [[]]
    for arg in args:
        ans = [x+[y] for x in ans for y in arg]
    return ans



def get_connectivity(params):
    """

    """
    a = params['a']
    radius = params['radius']
    radius_normalized = radius/a # when a is
    population = params['population']

    center = (int(population.dim[0]/2.),int(population.dim[1]/2.))
    offset = (int(round(population.dim[0]*radius_normalized)),int(round(population.dim[1]*radius_normalized)))

    targets={}
    targets_gid={}
    for n in range(center[0]-offset[0],center[0]+offset[0]+1):
        for m in range(center[1]-offset[1],center[1]+offset[1]+1):
            #print 'n: ',n, ' m: ',m
            gid = population[n,m]
            targets_tmp = pynest.getDict([gid])[0]['targets']
            targets_n_m = []
            for tgid in targets_tmp:
                targets_n_m.append(population.locate(tgid))

            targets[(n,m)]=targets_n_m
            targets_gid[gid]=targets_tmp
    return targets, targets_n_m




# TODO : this is not recommanded (SyntaxWarning)
from NeuroTools.benchmark import *
def run_simulations(model,url,tag):
    #lcn = LocalNetwork(0.1)
    #lcn = model
    print 'Simulations start'
    print '######################################'
    file = openHDF5File(url, "r")
    data_root = file.getStructure(nodepath = "/", structure = True)
    file.close()
    if data_root.has_key('benchmark_finished'):
        print 'Benchmark done. Data is in: ', url
        return

    params = data_root['params']



    experiments = data_root['run'].keys()
    #finished = False

    # create a tmp dict with the name of the tag
    root_dir = os.getcwd()
    if not(os.path.exists('results/'+tag)):
         os.mkdir('results/'+tag)

    # go there
    os.chdir('results/'+tag)

    # check if this node is the copy machine at the end
    # print 'check for copy node'
    # copy_node = False
    # if not(os.path.exists('final_copy_to_h5')):
     #   print 'I am the copy_node'
     #   copy_node = True
     #   os.system('touch final_copy_to_h5')
    #else:
    #    print 'NO, I am not the copy_node'

    if data_root['run'][experiments[0]].has_key('useHardware'):
        # check if have hardware
        haveHardware = False
        try:
            import pyNN.fhws1v2 as pyNN_
            haveHardware = True
        except ImportError:
            haveHardware = False

        if not haveHardware:
            print 'remove the hardware sim'
            # we have to sort out the sims that use the hardware
            experiments_tmp = []
            for param in data_root['run'].keys():
                if data_root['run'][param]['useHardware'] == False:
                    experiments_tmp.append(param)
            experiments = experiments_tmp


    #if data_root['run'][experiments[0]].has_key('analysed'):


    model.params = params
    for experiment in experiments:
        if data_root['run'][experiment].has_key('analysed'):
            continue

        name=''
        for key in data_root['run'][experiment].keys():
            name = name+'_'+key+'_'+str(data_root['run'][experiment][key])


        print 'checking name: ',name
        system_test_run=name+'_run'
        system_touch_run='touch '+name+'_run'
        system_touch_done='touch '+name+'_done'

        if os.path.exists(system_test_run):
            continue
        else:
            os.system(system_touch_run)

        # while not(finished):
        # experiment = experiments[int(numpy.floor(numpy.random.rand()*len(experiments)))] #TODO restrict to finished_experiments
        # if not(data_root['run'][experiment].has_key('begin')):#'end')):#
        # file = openHDF5File(url, "a")
        # file.setStructure({'begin' : datetime.datetime.now().isoformat()}, "/run/" + experiment, createparents = True)
        # file.close()

        #params = data_root['build']
        # update
        model.params.update(data_root['run'][experiment])
        print '\nExperiment:'
        print 'column_url', url
        #print 'retina_url', lcn.params['retina_url'], '\n'
        print 'Parameter to be simulated:'
        for key in data_root['run'][experiment].keys():
            print key,':',data_root['run'][experiment][key]


        model.params['name']=name


        #out_DATA = m7_1.run_m7_1(params)
        # run sim
        model.build_(model.params)
        # lcn.params['populations']=lcn.populations
        model.run_(model.params)
        #m7_1.run_m7_1(params)
        # touch the name_done such that the script knows that it has finished that sim
        os.system(system_touch_done)

    # Test if all sims have been simulated
    exps = experiments[:]
    #print 'exps ', exps
    all_done = False

    #for ex in exps:
    #    print ex


    for experiment in experiments:
        print 'experiment ',experiment
        name=''
        for key in data_root['run'][experiment].keys():
            name = name+'_'+key+'_'+str(data_root['run'][experiment][key])
        system_test_done=name+'_done'
        # print system_test_done
        if os.path.exists(system_test_done):
            exps.pop(exps.index(experiment))
        #print 'len exps', len(exps)
        #print exps
        if len(exps) == 0:
            all_done = True

    print 'all_done: ', all_done
    never_do=False
    #finished = False
    #if all_done: # until h5 is working properly, I will save data in gdf format, basta
    if never_do:
        exps = experiments[:]
        print 'I copy now'
        # while not(finished):
        for experiment in experiments:
            name=''
            for key in data_root['run'][experiment].keys():
                name = name+'_'+key+'_'+str(data_root['run'][experiment][key])

            #system_test_done=name+'_done'
            #    print system_test_done
            #    if os.path.exists(system_test_done):
            #experiments.pop(experiments.index(experiment))
            # print 'pop experiment'
            # update
            params.update(data_root['run'][experiment])
            params.update({'name': name})
            # read data
            #print lcn.populations
            print params
            out_DATA = model.get_data_(params)
            #out_DATA = m7_1.return_data(params)

            branch = '/run/' + experiment
            file = openHDF5File(url, "a")

            for pop in out_DATA:
                file.createSpikeList(branch, pop, rows = out_DATA[pop], dt = data_root['params']['dt'], spec = 'reltime_id', ref = None)

            # file.setStructure({'end' : datetime.datetime.now().isoformat()}, "/run/" + experiment, createparents = True)
            file.close()
            exps.pop(exps.index(experiment))


        if len(exps) == 0:
            #finished = True
            file = openHDF5File(url, "a")
            file.setStructure({'benchmark_finished':True}, "/", createparents = True)
            #file.setStructure({'benchmark_finished' : datetime.datetime.now().isoformat()}, "/run/", createparents = True)
            file.close()
            os.system('touch finished_all')

    if os.path.exists('finished_all'):
        print 'done all, data should be in: ', url
        # now del all in that dir
        files = os.listdir('../'+tag)
        for file in files:
            os.remove(file)
        os.chdir(root_dir)
        os.rmdir(tag)

    else:
        print 'not all done, still simulating, or copying.'

    os.chdir(root_dir)

import pysb.pathfinder

# Here set the path to BioNetGen
pysb.pathfinder.set_path('bng', '../BioNetGen-2.3.1-Linux/BioNetGen-2.3.1')
from pysb.importers.sbml import model_from_biomodels
from pysb.integrate import Solver
from pysb.tools import render_reactions, render_species
from collections import OrderedDict
from functools import partial
import pylab as pl
import numpy as np
import sys
import multiprocessing
from os import system
import pickle
import csv

# Choose the analysis and sampling methods
import SALib.analyze.sobol as analysis_method

an_m = 'sobol'
import SALib.sample.saltelli as sampling_method

smp_m = 'saltelli'

# Choose the model
# ACC_NO = '261'
# 259 is FEdeficient, 260 adequate and 261 loaded
ACC_NOs = {'deficient': '259', 'adequate': '260', 'loaded': '261'}
models = {key: model_from_biomodels(val) for key, val in ACC_NOs.items()}

# Set the simulation time
t = list(range(30))  # 30 minutes, as in publication.

# Choose species to be analyzed from the model
outputs = {key: list(models[key].observables.keys()) for key in
           ACC_NOs}  # all the observables! we need to run the analysis for all of them separately


# all three models have the same observables, they only differ in parameter values

# View the model
# print(model.monomers)
# print(models['adequate'].parameters)
# print(models['adequate'].rules)
# print(model.observables)
# print(model.initial_conditions)
# print(model.rules)
# print(model.compartments)

# Visualize the reactions and species in the model
def visualize():
    for ACC_NO in ACC_NOs:
        dot_r = render_reactions.run(models[ACC_NO])
        dotfile_r = open('model_reactions-' + ACC_NO + '.dot', 'w')
        dotfile_r.write(dot_r)
        dotfile_r.close()
        system('dot -T pdf model_reactions-' + ACC_NO + '.dot -o model_reactions-' + ACC_NO + '.pdf')
        dot_s = render_species.run(models[ACC_NO])
        dotfile_s = open('model_species-' + ACC_NO + '.dot', 'w')
        dotfile_s.write(dot_s)
        dotfile_s.close()
        system(
            'ccomps -x model_species-' + ACC_NO + '.dot | dot | gvpack -m0 | neato -n2 -T pdf -o model_species-' + ACC_NO + '.pdf')
        # Run the simulation
        solver = Solver(models[ACC_NO], t)
        solver.run()
        yout = solver.yobs
        # Plot simulation results for the observables in the model
        pl.ion()
        pl.figure()
        for observ in outputs[ACC_NO]:
            pl.plot(t, yout[observ], label=observ)
        pl.legend()
        pl.xlabel("Time (minutes)")
        pl.ylabel("Distribution of radioactive iron/body compartments")
        pl.savefig('model_results-' + ACC_NO + '.png')
        pl.clf()


visualize()

# Load the reference parameters - initials in the model, from publication
param_file = csv.reader(open('parameters.csv'), delimiter=' ')
headers = next(param_file)
param_bounds = {}
for row in param_file:
    param_bounds[row[0] + '_k1'] = [float(row[1]),
                                    float(row[2])]  # 0.9*min_bound,1.1*max_bound] # add some additional margin
for param in models['adequate'].parameters:
    if param.name not in param_bounds:
        param_bounds[param.name] = [0.99 * param.value, 1.01 * param.value]


def scatter_plot():
    n = 100
    mod = models['adequate']
    for param in param_bounds.keys():
        prob = {'num_vars': 1, 'names': [param], 'bounds': [param_bounds[param]]}
        par_sets = [p[0] for p in sampling_method.sample(prob, n)]
        solver = Solver(mod, t)
        res_sets = {out: [] for out in outputs['adequate']}
        for par_set in par_sets:
            pars = {p.name: p.value for p in mod.parameters}
            pars[param] = par_set
            solver.run(pars)
            for out in outputs['adequate']:
                res_sets[out].append(solver.yobs[out][-1])  # last result
        # Plot simulation results for the observables in the model
        pl.ion()
        pl.figure(figsize=(8, 8))
        for res in res_sets.keys():
            pl.scatter(par_sets, res_sets[res], label=res)
        pl.legend()
        pl.xlabel(param)
        pl.ylabel("Outputs")
        pl.savefig('scatterplots/' + param + '.png', dpi=300)
        pl.clf()


scatter_plot()

# Here the problem definition
N = 10

problem = {'num_vars': len(param_bounds), 'names': list(param_bounds.keys()),
           'bounds': list(param_bounds.values())}
# use this with Sobol
param_sets = sampling_method.sample(problem, N)

# use this with Morris
# param_sets = sampling_method.sample(problem, N, num_levels=4, grid_jump=2)
model = models['adequate']  # choose one, because essentialy they are all the same, if not for parameters and an error
solver = Solver(model, t)


def simulate(parameters):
    solver.run(parameters)
    return (solver.yobs)


p = multiprocessing.Pool(4)
func = partial(simulate)
simulations = p.map(func, param_sets)
results = {}
sa_results = {}
for output in outputs['adequate']:
    results[output] = [r[output][-1] for r in simulations]
    # Use this with Sobol
    sa_results[output] = analysis_method.analyze(problem, np.array(results[output]))

    # Use this with Morris
    # sa_results[output] = analysis_method.analyze(problem, param_sets, np.array(results[output]), conf_level=0.95,
    #                                             num_levels=4, grid_jump=2)

# Save results to a file for analysis
pickle.dump(sa_results, open("SA_results_" + str(N) + "_" + an_m + "_" + smp_m, "wb"))

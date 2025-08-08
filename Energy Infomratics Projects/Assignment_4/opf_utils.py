import pyomo.environ as pyo
from pyomo.opt.base import SolverFactory
import json

# Calculate the costs of power generation for a particular model.
#
# m = model
def calculate_costs(m):
    costs = 0
    for node, p in m.__getattribute__('x_p').items():
        costs += m.a[node] * p.value * p.value + m.b[node] * p.value + m.c[node]
    return costs

# Prints the results from a given solved model.
#
# m = model
# AC = alternative current system or not
def print_result(m, AC=True):
    def print_nodes_for_entity(name, entity_name):
        print('\nNode\t\t' + name)
        entity = m.__getattribute__(entity_name)
        for i, val in entity.items():
            print(i + 1, '\t\t', val.value)

    print('\n\n\nSolution\n')

    # Print Voltage Angle at Node
    print_nodes_for_entity('Voltage angle [deg]','x_v_angle')

    # Print Voltage Magnitude at Node
    if AC:
        print_nodes_for_entity('Voltage magnitude','x_v')

    # Print Generator Power at Node
    print_nodes_for_entity('Power','x_p')

    # Print Reactive Power at Node
    if AC:
        print_nodes_for_entity('Reactive Power','x_q')

    costs = calculate_costs(m)
    print('\nCost: ', costs)


# Solves a model and prints the result.
#
# m = model
# AC = alternative current system or not
def solve_and_print(m, AC=True):
    optim = pyo.SolverManagerFactory('neos')
    if AC:
        optim.solve(m, solver='conopt', tee=True)
    else:
        optim.solve(m, solver='cplex', tee=True)
    print_result(m, AC)

# Creates a json file for upload of optimization results.
# NOTE: This overrides the existing file.
#
# m = model
# team_name = your team name, ex. 'Group_A'
# AC = alternative current system or not
def create_results_json(m, team_name, AC=True):
    def get_entity_for_node(i, entity_name):
        entity = m.__getattribute__(entity_name)
        return entity[i].value

    results_dict = {}
    results_dict['name'] = team_name
    results_dict['task'] = 'AC' if AC else 'DC'
    results_dict['cost'] = calculate_costs(m)
    results_dict['nodes'] = []

    for i in m.nodes:
        node = {}
        node['id'] = i + 1
        node['voltage_angle'] = get_entity_for_node(i, 'x_v_angle')
        if AC:
            node['voltage_magnitude'] = get_entity_for_node(i, 'x_v')
        node['power'] = get_entity_for_node(i, 'x_p')
        if AC:
            node['reactive_power'] = get_entity_for_node(i, 'x_q')
        results_dict['nodes'].append(node)

    filename = '{team_name}_{task}.json'.format(team_name=results_dict['name'], task=results_dict['task'])
    with open(filename, 'w') as outfile:
        json.dump(results_dict, outfile, indent=2)
    print('Results written to ' + filename)
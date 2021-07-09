import itertools
import functools
import json
import networkx as nx
import matplotlib.pyplot as plt

import dimod

from dimod.generators.constraints import combinations
from dwave.system import LeapHybridSampler

def build_bqm(leaves, internals, root):

    non_root_nodes = leaves + internals
    non_leaf_nodes = internals + [root]
    all_nodes = non_root_nodes + [root]

    bqm = dimod.BinaryQuadraticModel({}, {}, 0.0, dimod.BINARY)

    def append_linear_safe(x, value):
        try:
            bqm.set_linear(x, bqm.get_linear(x) + value)
        except:
            bqm.set_linear(x, value)

    def append_quadratic_safe(x, value):
        try:
            bqm.quadratic[x] = bqm.quadratic[x] + value
        except:
            bqm.quadratic[x] = value

    # -------- H_A Hamiltonian --------
    # Unfortunately we need to introduce some auxiliary variables y_ab 
    # which is encoded as (a,b) where a -> b which means that b is a substring of a
    # Tree constraint: the join tree needs to have the root which contains the result

    root_bqm = combinations([root], 1, strength=1)
    bqm.update(root_bqm)

    # Tree constraint: the join tree needs to have all the single tables as leaves

    leaves_bqm = combinations(leaves, len(leaves), strength=1)
    bqm.update(leaves_bqm)


    # Tree constraint: every edge has the source and target vertex in the tree
    for x in all_nodes:

        c = 1
        nodes_lower = nodes_lower_than(x, all_nodes)
        
        for y in nodes_lower:
            
            y_vu_key = (y, x)
            key1 = (x, y_vu_key)
            key2 = (y, y_vu_key)
            key3 = (y, x)

            append_linear_safe(y_vu_key, 2*c)
            append_linear_safe(x, c)
            append_linear_safe(y, c)

            append_quadratic_safe(key1, -4*c)
            append_quadratic_safe(key2, -4*c)
            append_quadratic_safe(key3, 2*c)

    # Tree constraint: every non-root node needs to have exactly one parent
    for x in non_root_nodes:
        nodes_lower = nodes_lower_than(x, all_nodes)
        c = 1
        if len(nodes_lower) > 0:

            append_linear_safe(x, c)

            for y in nodes_lower:

                y_vu_key = (y, x)
                append_linear_safe(y_vu_key, c)

                xy_key = (x, y_vu_key)
                append_quadratic_safe(xy_key, -2*c)

            processed = []
            for y in nodes_lower:
                for z in nodes_lower:
                    if z != y and z not in processed:

                        processed.append(y)

                        yz_key = ((y, x), (z, x))
                        
                        append_quadratic_safe(yz_key, 2*c)


    # Tree constraint: every node expect the leaves need to have exactly two children 
    # i.e. the join is tree is a full binary tree
    # The opened bqm is (2x - y - z - ...)^2 = 4x + y + z + ... - 4xy - 4xz - ... + 2yz + ...
    for x in non_leaf_nodes:
        c = 32
        nodes_higher = nodes_higher_than(x, all_nodes)
        #print(nodes_higher)

        if len(nodes_higher) > 0:

            # 4x
            append_linear_safe(x, 4*c)

            for y in nodes_higher:
                
                # + y + z + ...
                y_vu_key = (x, y)
                append_linear_safe(y_vu_key, c)

                # - 4xy - 4xz - ...
                xy_key = (x, y_vu_key)
                append_quadratic_safe(xy_key, -4*c)

            # 2yz + ...
            processed = []
            for y in nodes_higher:
                for z in nodes_higher:
                    if z != y and z not in processed:
                        processed.append(y)

                        x_set = set([char for char in x])
                        y_set = set([char for char in y])
                        z_set = set([char for char in z])

                        if x_set == y_set.union(z_set):

                            yz_key = ((x, y), (x, z))
                            #print(yz_key)

                            append_quadratic_safe(yz_key, 2*c)


    # -------- H_B Hamiltonian --------
    # Cost constraints

    # In the worst case, any non-root table can be joined with any other non-root table except with itself 
    # for x in non_root_nodes:
    #     for y in non_root_nodes:
    #         if x != y:
    #             key = (x, y)
    #             bqm.quadratic[key] = 1 #cost_function(key)

    # print(bqm.linear)
    # print()
    # print(bqm.quadratic)
    # print()

    for elem in bqm.linear:
        print(elem, bqm.linear[elem])

    for elem in bqm.quadratic:
        print(elem, bqm.quadratic[elem])

    return bqm

def nodes_lower_than(x, all_nodes):
    x_set = set([char for char in x])
    nodes = []
    for y in all_nodes:
        y_set = set([char for char in y])
        if x_set.issubset(y_set) and y != x:
            nodes.append(y)
    return nodes

def nodes_higher_than(x, all_nodes):
    x_set = set([char for char in x])
    nodes = []
    for y in all_nodes:
        y_set = set([char for char in y])
        if y_set.issubset(x_set) and y != x:
            nodes.append(y)
    return nodes


def cost_function(tables):
    return 1

def solve_join_order(bqm):
    sampler = LeapHybridSampler()
    sampleset = sampler.sample(bqm, label='Join Order Selection')
    sample = sampleset.first.sample
    # energy = sampleset.first.energy
    for varname, value in sample.items():
        print(varname, value)

if __name__ == "__main__":

    tables = ["R", "S", "T"]
    table_string = functools.reduce(lambda a, b: a + b, tables, "")
    # The root of the join order plan is the final join of all the tables
    root = table_string
    leaves = tables
    internals = []

    for i in range(2, len(tables)):
        temp = itertools.combinations(table_string, i)
        for elem in temp:
            internals.append(functools.reduce(lambda a, b: a + b, elem, ""))

    #print(internals)
    #print(root)
    #print(leaves)

    # Solve BQM and update matrix
    bqm = build_bqm(leaves, internals, root)
    #g = dimod.to_networkx_graph(bqm)
    #nx.draw(g, with_labels = True)
    #plt.savefig("filename.png")
    solve_join_order(bqm)


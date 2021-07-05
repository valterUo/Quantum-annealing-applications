import itertools
import functools
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

    # -------- H_A Hamiltonian --------
    # Tree constraint: the join tree needs to have the root which contains the result
    root_bqm = combinations([root], 1, strength=1)
    bqm.update(root_bqm)

    # Tree constraint: the join tree needs to have all the single tables as leaves
    leaves_bqm = combinations(leaves, len(leaves), strength=1)
    bqm.update(leaves_bqm)

    # Tree constraint: every non-root node needs to have exactly one parent node
    for x in non_root_nodes:
        bqm.set_linear(x, bqm.get_linear(x) + 1)
        nodes_lower = nodes_lower_than(x, all_nodes)
        for y in nodes_lower:
            try:
                bqm.set_linear(y, bqm.get_linear(y) + 1)
            except:
                bqm.set_linear(y, 1)
            key = (x,y)
            try:
                bqm.quadratic[key] = bqm.quadratic[key] - 1
            except:
                bqm.quadratic[key] = -1
            for z in nodes_lower:
                if z != y:
                    key = (z,y)
                    try:
                        bqm.quadratic[key] = bqm.quadratic[key] + 1
                    except:
                        bqm.quadratic[key] = 1

    # Tree constraint: every node expect the leaves need to have exactly two children i.e. the join is tree is a full binary tree
    # Unfortunately we need to introduce some auxiliary variables y_ab which is encoded as (a,b) where a -> b which means that b is a substring of a
    # The opened bqm is (2x - y - z)^2 = 4x + y + z - 4xy - 4xz + 2yz

    for x in non_leaf_nodes:
        # x = x_v
        x_set = set([char for char in x])

        # the first linear term 4x

        try:
            bqm.set_linear(x, bqm.get_linear(x) + 4)
        except:
            bqm.set_linear(x, 4)

        for y in non_root_nodes:
            if y != x:
                y_set = set([char for char in y])
                for z in non_root_nodes:
                    if z != x and z != y:
                        z_set = set([char for char in z])
                        if y_set.issubset(x_set) and z_set.issubset(x_set):

                            y_vu_key = (x, y) # y = y_vu
                            y_vu2_key = (x, z) # z = y_vu2

                            # Linear terms y + z
                            
                            try:
                                bqm.set_linear(y_vu_key, bqm.get_linear(y_vu_key) + 1)
                            except:
                                bqm.set_linear(y_vu_key, 1)

                            try:
                                bqm.set_linear(y_vu2_key, bqm.get_linear(y_vu2_key) + 1)
                            except:
                                bqm.set_linear(y_vu2_key, 1)

                            # Quadratic terms i.e. the three last - 4xy - 4xz + 2yz

                            xy_key = (x, y_vu_key)
                            xz_key = (x, y_vu2_key)
                            yz_key = (y_vu_key, y_vu2_key)

                            try:
                                bqm.quadratic[xy_key] = bqm.quadratic[xy_key] - 4
                            except:
                                bqm.quadratic[xy_key] = -4

                            try:
                                bqm.quadratic[xz_key] = bqm.quadratic[xz_key] - 4
                            except:
                                bqm.quadratic[xz_key] = -4

                            try:
                                bqm.quadratic[yz_key] = bqm.quadratic[yz_key] + 2
                            except:
                                bqm.quadratic[yz_key] = 2


    # -------- H_B Hamiltonian --------
    # Cost constraints. The root i.e. the final result cannot be joined with any other table

    # In the worst case, any non-root table can be joined with any other non-root table except with itself 
    # for x in non_root_nodes:
    #     for y in non_root_nodes:
    #         if x != y:
    #             key = (x, y)
    #             bqm.quadratic[key] = 1 #cost_function(key)

    print(bqm.linear)
    print()
    print(bqm.quadratic)
    print()
    return bqm

def nodes_lower_than(x, all_nodes):
    nodes = []
    for y in all_nodes:
        if len(y) > len(x):
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


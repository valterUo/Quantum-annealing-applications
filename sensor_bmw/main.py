import dimod

from dimod.generators.constraints import combinations
from dwave.system import LeapHybridSampler

def build_bqm(leaves, internals, root):
    return None

def solve_join_order(bqm):
    sampler = LeapHybridSampler()
    sampleset = sampler.sample(bqm, label='BMW Sensor positioning')
    sample = sampleset.first.sample
    # energy = sampleset.first.energy
    for varname, value in sample.items():
        print(varname, value)

if __name__ == "__main__":

    # Solve BQM and update matrix
    bqm = build_bqm(leaves, internals, root)
    #g = dimod.to_networkx_graph(bqm)
    #nx.draw(g, with_labels = True)
    #plt.savefig("filename.png")
    solve_join_order(bqm)
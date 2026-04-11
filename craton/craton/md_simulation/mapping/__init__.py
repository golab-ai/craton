from .pairnetwork import PairNetwork

def pair_network_init(ligands,
                      topology="normal",
                      user_pair_list=None,
                      bias_nodes=None,
                      core=None,
                      nbunch=None):
    return PairNetwork.create_graph_from_molecules(ligands,
                                                   topology=topology,
                                                   user_pair_list=user_pair_list,
                                                   bias_nodes=bias_nodes,
                                                   core=core,
                                                   nbunch=nbunch)

def atom_mapping_calculate(gg):
    PairNetwork.calculate_atom_mapping(gg)

def molecule_similiarity_calculate(gg):
    PairNetwork.calculate_similarity(gg)

def pair_network_final(gg,topology="normal",nbunch=None,bias_nodes=None):
    if  topology== "normal":
        PairNetwork.reduce_normal_graph(gg,nbunch=nbunch,bias_nodes=bias_nodes)

def graph_attributes_report(gg,topology="normal"):
    PairNetwork.report_graph_attributes(gg,topology=topology)



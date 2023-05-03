import numpy as np
import pandas as pd
import numpy as np
import networkx as nx



class SynGraph:
    type_subs = {
        'Add': 'add',
        'Dissolve': 'add',
        'Extract': 'del',
        'Filter': 'del',
        'Partition': 'del',
        'Partiton': 'del',
        'partition': 'del',
        'Precipitate': 'del',
        'Purify': 'del',
        'Recover': 'del',
        'Remove': 'del',
        'Yield': 'del',
        'ApparatusAction': 'context_dependent',
        'Cool': 'trans',
        'Degass': 'trans',
        'Dry': 'del',
        'dry': 'del',
        'Heat': 'trans',
        'Stir': 'trans',
        'Synthesize': 'context_dependent',
        'synthesize': 'context_dependent',
        'Wait': 'trans',
        'wait': 'trans',
        'Wash': 'del',
        'wash': 'del',
        'Seal': 'trans',
        'Evacuate': 'trans',
        'Mix': 'add',
        'Suspend': 'add',
        'Pelletize': 'trans'
    }

    def __init__(self): pass


    def create_initial_network(self,sequence_df):
        nodes = {}
        for c, row in sequence_df.reset_index().iterrows():
            nodes[c] = {
                'name': row['Step type'],
                'steptype': self.type_subs[row['Step type']],
                'chemicals': row['chemicals and amount sued'],
                'time': row['Time'],
                'temperature': row['Temp'],
                'comments': row['comment']
            }

        nodes_with_attibutes = [(k, v) for k, v in nodes.items()]
        edges = []
        for c, _ in sequence_df.reset_index().iterrows():
            if c == 0:
                continue
            edges.append((c - 1, c))

        output = nx.DiGraph()
        output.add_nodes_from(nodes_with_attibutes)
        output.add_edges_from(edges)
        return output


    def create_coarse_grained_graph(self,initial_graph):
        from collections import defaultdict
        groups = defaultdict(list)
        undirected_graph = initial_graph.to_undirected()

        for type_string in ['add', 'del', 'trans', 'context_dependent']:
            indices = \
            np.where(np.vstack(list(nx.get_node_attributes(undirected_graph, 'steptype').values())) == type_string)[0]
            subgraph = nx.subgraph(undirected_graph, indices)
            groups[type_string].extend(list(nx.connected_components(subgraph)))

        new_groups = sorted([x for y in list(groups.values()) for x in y], key=lambda x: list(x)[0])
        #     print(new_groups)
        coarse_grained_graph = nx.DiGraph()

        for c, g in enumerate(new_groups):
            group_label = initial_graph.nodes.data()[list(g)[0]]['steptype']
            #         print(c, list(g)[0], group_label)
            coarse_grained_graph.add_node(c, steptype=group_label, original_nodes=g)
            if c > 0:
                coarse_grained_graph.add_edge(c - 1, c)

        for x in coarse_grained_graph.nodes.data():
            # print(x[0], x[1]['original_nodes'])
            aggregated_data = {
                'name': [],
                'steptypes': [],
                'chemicals': [],
                'times': [],
                'temps': [],
                'comments': []
            }
            for y in x[1]['original_nodes']:
                rawdata = initial_graph.nodes.data()[y]
                aggregated_data['name'].append(rawdata['name'])
                aggregated_data['steptypes'].append(rawdata['steptype'])
                aggregated_data['chemicals'].append(rawdata['chemicals'])
                aggregated_data['times'].append(rawdata['time'])
                aggregated_data['temps'].append(rawdata['temperature'])
                aggregated_data['comments'].append(rawdata['comments'])
                nx.set_node_attributes(coarse_grained_graph, {x[0]: aggregated_data})

        return coarse_grained_graph

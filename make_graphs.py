from pyvis.network import Network
import pandas as pd
import sys 
import os


folder = sys.argv[1]

file_list = os.listdir(folder)



for gene_file in file_list :
    
    # set the physics layout of the network
    reference = gene_file.split(".")[0]
    got_net = Network(height="75%", width="75%", font_color="black", heading=reference)
    got_net.barnes_hut()
    
    got_data = pd.read_csv(folder+"/"+gene_file,sep='\t')

    sources = got_data['Gene']
    targets = got_data['Others']
    gene_class = got_data['Gene_class']
    others_Class = got_data['Others_Class']
    weight_gene = got_data['Total number of interactions gene']
    weight_other = got_data['Total number of interactions others']

    edge_data = zip(sources, targets, weight_gene, weight_other, gene_class, others_Class)


    for e in edge_data:
        src = e[0]
        dst = e[1]
        weight_1 = e[2]
        weight_2 = e[3]
        gene_class = e[4]
        others_class = e[5]
        if gene_class == "Known" :
            if src == reference :
                got_net.add_node(src, src, title=src, group=1, shape = "dot", size = int(weight_1)*4+1000 , physics = True , borderWidth = 3 , borderWidthSelected = 7)
            else:
                got_net.add_node(src, src, title=src, group=2, shape = "dot", size = int(weight_1)*4 , physics = True , borderWidth = 3 , borderWidthSelected = 7)
        if gene_class == "Predicted" :
            if src == reference :
                got_net.add_node(src, src, title=src, group=1, shape = "dot", size = int(weight_1)*4+1000, physics = True , borderWidth = 3 , borderWidthSelected = 7)
            else:
                got_net.add_node(src, src, title=src, group=3, shape = "dot", size = int(weight_1)*4, physics = True , borderWidth = 3 , borderWidthSelected = 7)
        if others_class == "ANOTHERGEN" :
            got_net.add_node(dst, dst, title=dst, group=4, shape = "box", size = int(weight_2)+50, physics = True , level = 2)
        if others_class == "GOterm" :
            got_net.add_node(dst, dst, title=dst, group=5, shape = "box", size = int(weight_2)+50, physics = True , level = 2)
        if others_class == "disease" :
            got_net.add_node(dst, dst, title=dst, group=6, shape = "box", size = int(weight_2)+50, physics = True , level = 2)
        got_net.add_edge(src, dst, physics = True)

        neighbor_map = got_net.get_adj_list()

# add neighbor data to node hover data
    for node in got_net.nodes:
        node["title"] += " Neighbors:<br>" + "<br>".join(neighbor_map[node["id"]])
        node["value"] = len(neighbor_map[node["id"]])

    got_net.show(folder+gene_file.split(".")[0]+".html")

from pyvis.network import Network
import os
import pandas as pd

def net_visualization(folder, add_legend=False):
  #to change the node shape, change it in here
  #to change the node colour, change it in make_options
  #this is because got_net.add_node by default sets shape=dot if not provided, and this overrides the shape set in the groups options
	file_list = os.listdir(folder)
	file_list = [x for x in file_list if x.endswith('.txt')]
	for gene_file in file_list:
		# set the physics layout of the network
		reference = gene_file.split(".")[0]
		got_net = Network(height="75%", width="75%", font_color="black", heading=reference)
		got_net.barnes_hut()
		
		got_data = pd.read_csv(folder+"/"+gene_file,sep='\t')
		title_available = 'edge_title' in got_data.columns

		for index, row in got_data.iterrows():
			src = row['Gene']
			dst = row['Others']
			weight_1 = row['Total number of interactions gene']
			weight_2 = row['Total number of interactions others']
			gene_class = row['Gene_class']
			others_class = row['Others_Class']
			if gene_class == "Known" :
				if src == reference :
					got_net.add_node(src, src, title=src, group='Query', shape='triangle', size = int(weight_1)*4+1000 , physics = True , borderWidth = 3 , borderWidthSelected = 7)
				else:
					got_net.add_node(src, src, title=src, group=gene_class, shape='dot', size = int(weight_1)*4 , physics = True , borderWidth = 3 , borderWidthSelected = 7)
			if gene_class == "Predicted" :
				if src == reference :
					got_net.add_node(src, src, title=src, group='Query', shape='triangle',  size = int(weight_1)*4+1000, physics = True , borderWidth = 3 , borderWidthSelected = 7)
				else:
					got_net.add_node(src, src, title=src, group=gene_class, shape='dot', size = int(weight_1)*4, physics = True , borderWidth = 3 , borderWidthSelected = 7)
			if others_class == "ANOTHERGEN" :
				got_net.add_node(dst, dst, title=dst, group=others_class, shape="box",  size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "GOterm" :
				got_net.add_node(dst, dst, title=dst, group=others_class, shape="box",  size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "disease" :
				got_net.add_node(dst, dst, title=dst, group=others_class, shape="box",  size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "pubmedID" :
				got_net.add_node(dst, dst, title=dst, group=others_class, shape="box",  size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "expression" :
				got_net.add_node(dst, dst, title=dst, group=others_class, shape="box",  size = int(weight_2)+50, physics = True , level = 2)
			if others_class == "molecule" :
				got_net.add_node(dst, dst, title=dst, group=others_class, shape="box", size = int(weight_2)+50, physics = True , level = 2)
			
			if title_available:
				got_net.add_edge(src, dst, physics = True, title=row['edge_title'])
			else:
				got_net.add_edge(src, dst, physics = True)

		neighbor_map = got_net.get_adj_list() #this was inside the loop before and should not be


		# add neighbor data to node hover data
		for node in got_net.nodes:
			adj = neighbor_map[node["id"]]
			n_adj = len(adj)
			if n_adj > 10:
				node["title"] += f"<br>Neighbors:<br> {n_adj}"
			else:   
				node["title"] += "<br>Neighbors:<br>" + "<br>".join(adj)
			node["value"] = n_adj

		# add options object
		got_net.set_options(make_options())

		got_net.show(folder+'/graph__'+gene_file.split(".")[0]+".html")
		with open(folder+'/graph__'+gene_file.split(".")[0]+".html", 'r') as origin:
			modified = []
			for line in origin:
				if 'nodes = new vis.DataSet' in line:
					a = line.replace('        ','    ')
				elif 'edges = new vis.DataSet' in line:
					b = line.replace('        ','	')
				elif '// parsing and collecting nodes and edges from the python' in line:
					e = line
				else:
					modified.append(line)
		
		flag = 0

		groups = ["Known","Predicted","ANOTHERGEN","GOterm","disease","pubmedID","expression","molecule","Query"]
		if add_legend:
			legend = """
				// legend
				var mynetwork = document.getElementById("mynetwork");
				var x = -mynetwork.clientWidth * 1.5;
				var y = -mynetwork.clientHeight / 2 + 50;
				var step = 70;
				nodes.add([
				"""
			for g in range(len(groups)):
				grp = groups[g]
				s = f'{{x: x, y: y + step*{g}, label: "{grp}", group: "{grp}", fixed: true, physics: false, }},'
				legend += s
			legend += """
					]);
					// get a JSON object
					var allNodes = nodes.get({ returnType: "Object" });
				"""
		else:
			legend = """
			// get a JSON object
			var allNodes = nodes.get({ returnType: "Object" });
			"""

		modified.insert(modified.index('    var options, data;\n')+1,'    var highlightActive = false;\n')
		modified.insert(modified.index('    var options, data;\n')+2, e)
		modified.insert(modified.index('    var options, data;\n')+2, a)
		modified.insert(modified.index('    var options, data;\n')+3, b)
		modified.insert(modified.index('    function drawGraph() {\n')-2,legend)
		modified = modified[0:modified.index( '        network = new vis.Network(container, data, options);\n')+1]
		
		c = get_custom_js()


		d = '''
      </body>
      </html>
		'''

		os.remove(folder+'/graph__'+gene_file.split(".")[0]+".html")
		with open(folder+'/'+gene_file.split(".")[0]+".html", 'a') as newer:
			for i in range(len(modified)-2):
				if modified[i] == '<script type="text/javascript">\n':
					flag = 1
				if flag == 1:
					newer.write(modified[i])
			
			newer.write(c)
		with open(folder+'/graph_'+gene_file.split(".")[0]+".html", 'a') as highlighted:
			for i in range(len(modified)):
				highlighted.write(modified[i])
			highlighted.write(c)
			highlighted.write(d)

def make_options():
    
    options = """
    var options = {
    "configure": {
        "enabled": false
    },
    "edges": {
        "color": {
            "inherit": true
        },
        "smooth": {
            "enabled": false,
            "type": "continuous"
        }
    },
    "interaction": {
        "dragNodes": true,
        "hideEdgesOnDrag": false,
        "hideNodesOnDrag": false
    },
    "physics": {
        "barnesHut": {
            "avoidOverlap": 0,
            "centralGravity": 0.3,
            "damping": 0.09,
            "gravitationalConstant": -80000,
            "springConstant": 0.001,
            "springLength": 250
        },
        "enabled": true,
        "stabilization": {
            "enabled": true,
            "fit": true,
            "iterations": 500,
            "onlyDynamicEdges": false,
            "updateInterval": 50
        }
    },
    "layout": {
        "improvedLayout": true
    },
    "groups": {
      "Known": {
        "shape": "dot",
        "color": "#1f77b4"
      },
      "Predicted": {
        "shape": "dot",
        "color": "#aec7e8"
      },
      "ANOTHERGEN": {
        "shape": "box",
        "color": "#ff7f0e"
      },
      "GOterm": {
        "shape": "box",
        "color": "#ffbb78"
      },
      "disease": {
        "shape": "box",
        "color": "#2ca02c"
      },
      "pubmedID": {
        "shape": "box",
        "color": "#98df8a"
      },
      "expression": {
        "shape": "box",
        "color": "#d62728"
      },
      "molecule": {
        "shape": "box",
        "color": "#ff9896"
      },
      "Query": {
        "shape": "triangle",
        "color": "#9467bd"
      }
      
    }
}
    """
    return options

def get_custom_js():
  x = '''
						        
        network = new vis.Network(container, data, options);
        network.on("click", neighbourhoodHighlight);
        network.physics.physicsEnabled = false
		network.on("stabilizationProgress", function(params) {
      //check for the loading bar, otherwise small graphs that load quickly are buggy
      var lb = document.getElementById('loadingBar') !== null
      if(lb){
        document.getElementById('loadingBar').removeAttribute("style");
		var maxWidth = 496;
		var minWidth = 20;
		var widthFactor = params.iterations/params.total;
		var width = Math.max(minWidth,maxWidth * widthFactor);

		document.getElementById('bar').style.width = width + 'px';
		document.getElementById('text').innerHTML = Math.round(widthFactor*100) + '%';
      }
		
    });
    network.once("stabilizationIterationsDone", function() {
      var lb = document.getElementById('loadingBar') !== null
        if(lb){
      document.getElementById('text').innerHTML = '100%';
      document.getElementById('bar').style.width = '496px';
      document.getElementById('loadingBar').style.opacity = 0;
      // really clean the dom element
      setTimeout(function () {document.getElementById('loadingBar').style.display = 'none';}, 500);
      }
                      });
                      
          return network;

      }


      function neighbourhoodHighlight(params) {
        // if something is selected:
        if (params.nodes.length > 0) {
          highlightActive = true;
          var i, j;
          var selectedNode = params.nodes[0];
          var degrees = 2;

          // mark all nodes as hard to read.
          for (var nodeId in allNodes) {
            allNodes[nodeId].color = "rgba(200,200,200,0.5)";
            if (allNodes[nodeId].hiddenLabel === undefined) {
              allNodes[nodeId].hiddenLabel = allNodes[nodeId].label;
              allNodes[nodeId].label = undefined;
            }
          }
          var connectedNodes = network.getConnectedNodes(selectedNode);
          var allConnectedNodes = [];

          // get the second degree nodes
          for (i = 1; i < degrees; i++) {
            for (j = 0; j < connectedNodes.length; j++) {
              allConnectedNodes = allConnectedNodes.concat(
                network.getConnectedNodes(connectedNodes[j])
              );
            }
          }

          // all second degree nodes get a different color and their label back
          for (i = 0; i < allConnectedNodes.length; i++) {
            allNodes[allConnectedNodes[i]].color = "rgba(150,150,150,0.75)";
            if (allNodes[allConnectedNodes[i]].hiddenLabel !== undefined) {
              allNodes[allConnectedNodes[i]].label =
                allNodes[allConnectedNodes[i]].hiddenLabel;
              allNodes[allConnectedNodes[i]].hiddenLabel = undefined;
            }
          }

          // all first degree nodes get their own color and their label back
          for (i = 0; i < connectedNodes.length; i++) {
            allNodes[connectedNodes[i]].color = undefined;
            if (allNodes[connectedNodes[i]].hiddenLabel !== undefined) {
              allNodes[connectedNodes[i]].label =
                allNodes[connectedNodes[i]].hiddenLabel;
              allNodes[connectedNodes[i]].hiddenLabel = undefined;
            }
          }

          // the main node gets its own color and its label back.
          allNodes[selectedNode].color = undefined;
          if (allNodes[selectedNode].hiddenLabel !== undefined) {
            allNodes[selectedNode].label = allNodes[selectedNode].hiddenLabel;
            allNodes[selectedNode].hiddenLabel = undefined;
          }
        } else if (highlightActive === true) {
          // reset all nodes
          for (var nodeId in allNodes) {
            allNodes[nodeId].color = undefined;
            if (allNodes[nodeId].hiddenLabel !== undefined) {
              allNodes[nodeId].label = allNodes[nodeId].hiddenLabel;
              allNodes[nodeId].hiddenLabel = undefined;
            }
          }
          highlightActive = false;
        }

        // transform the object into an array
        var updateArray = [];
        for (nodeId in allNodes) {
          if (allNodes.hasOwnProperty(nodeId)) {
            updateArray.push(allNodes[nodeId]);
          }
        }
        nodes.update(updateArray);
      }
    

    var g = drawGraph();

    </script>'''
  return x

if __name__ == '__main__':
    folder = './tests'
    file_list = os.listdir(folder)
    html = [x for x in file_list if x.endswith('.html')]
    file_list = [x for x in file_list if x.endswith('.txt')]
    for fname in html:
        os.remove(folder+'/'+fname)
    net_visualization(folder, True)
# Network visualization with NetworkX and yED

For the visualization of the networks I'm using yED (Version 3.21.1) Ubuntu 22.04.3 LTS

The jupyter notebook __network_visualization_examples.ipynb__ contains exmples of how to use networkx's built-in functions (mainly relying on matplolib) to visualize graphs, and also how to add node properties that can be than used in yED when editing graphs manually. 

To map properties added to a networkx file (which should beexported to a .graphml) in yED, follow these steps: 

1. Load a .graphml into yED via File → Open or by dragging and dropping it through the screen. If the software gives warnings about the conversion of properties data just OK them. 

2. Open the Properties Mapper, by clicking Edit → Properties Mapper 

3. Load a configuration file (cnfx) or create your own mapping. 

5. Click Layout → Organic (here the settings can be changed) and OK

6. To remove directional arrows (all the networks are undirected) click on an edge and press Ctrl+A to select all edges. Then in the Properties View (bottom right) set the Target Arrow option to a line (just like the Source Arrow above it)

From this point every other change is manual fine-tuning to the author’s taste. One can move nodes around by first selecting them then moving them around (or by adjusting the layout parameters) and most properties (node size, label font, colors, etc) can be changed in the Properties View. For more systematic changes one can adjust the configuration in the Properties Mapper. You’re smart, you’ll figure it out. 

Bonus tip: add the color legend boxes as big nodes next to the network (or Copy-Paste them from the existing graphmls). To add a new node drag and drop the desired shape from the Palette panel on the top right. 

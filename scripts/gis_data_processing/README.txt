Scripts in this directory are used to process and generate the GIS and flow data used to initialised the model

Currently, to produce the required data the scripts need to be run in a certain order as they each do a different step of the processing. This file explains this process.

1. extractITNRouting.py

This gets the directoinality information for road links and road nodes for all of the ITN road link data and saves it

2. processRoadNetworkData

Selects the road network data that lies in the study area. Cleans the Open Roads data to simplify the lines to have minimal angular deviation along a line.

3. processOSTopographicData.py

Given a polygon covering the study area, this extracts the road links and nodes, pedestrian and vehicles topographic space, and lines obstructing pedestrian movement that should be included to model that study area.

It also processes and cleans this data, for example linking road links (both open road and ITN) with vehicle and pedestrian space

4. editORRoadsToMatchPedPolys.py

The processOSTopographicData.py script splits the pedestrian polygons so that they link with OR open road links. However, some links, due to being small or far from the pavement, do not get a pedestrian polygon.

This script removes and replaces these links so that the road network more accurately reflects the pavement network.

After running this script, the processOSTopographicData.py script needs to be re runs so that polygons are created using the updated versions of the OR road links.

5. makeITNdirectional.py

This uses the direction information extracted with the first script and edits the portion of the road network selected in the section so that it represents a directed road network.

6. checkDiNetworkRouting.py

NOTE: Requires manually placing shape files of vehicle OD nodes in the processed data directory.

Given the road nodes to use as vehicle ODs this procudes random flows between ODs, checking that a route is possible.

7. randomOD2.py

Very simple script to generate random flows between pedestrian ODs

8. createPedestrianNetwork.py

Creates a network approximating the paths available to a pedestrian in that separate side of the road have different nodes with links to connect them across junctions
Scripts in this directory are used to process and generate the GIS and flow data used to initialised the model

Currently, to produce the required data the scripts need to be run in a certain order as they each do a different step of the processing. This file explains this process.

1. extractITNRouting.py

This gets the directoinality information for road links and road nodes for all of the ITN road link data and saves it

2. processOSTopographicData.py

Given a polygon covering the study area, this extracts the road links and nodes, pedestrian and vehicles topographic space, and lines obstructing pedestrian movement that should be included to model that study area.

It also processes and cleans this data, for example linking road links with vehicle and pedestrian space.

3. editITN.py

This uses the direction information extracted with the first script and edits the portion of the road network selected in the section so that it represents a directed road network.

4. checkDiNetworkRouting.py

Given the road nodes to use as vehicle ODs this procudes random flows between ODs, checking that a route is possible.

5. randomOD2.py

Very simple script to generate random flows between pedestrian ODs
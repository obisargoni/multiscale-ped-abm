# repastInterSim

This is a Repast Simphony agent-based model (ABM) of pedestrian and vehicle movement in an urban road network.

The objective of this ABM is to model pedestrian road crossing behaviour that is compliant (pedestrians walk along pavements and crossing at established crossings) and non-compliant (pedestrians cross at locations and at times that minimise their journey time).

These behaviours are of interest for two main reasons. Firstly, by allowing pedestrian agents to have heterogeneous and spatiotemporally dependant road crossing behaviour better reflects the diversity in road crossing behaviour observed in urban areas and allows the models to better represent real pedestrian journeys. Secondly, models of diverse road crossing behaviour provide a good environment for testing models of pedestrian-vehicle interaction which are being developed for Autonomous Vehicles.

![A screenshot for the model, showing pedestrians in blue and vehicles in black. Note how a pedestrian agent is crossing the road at a location outside of the 'crossings' (the two small rectangle sections of the road). This behaviour is non-compliant.](screenshot_non_compliant.PNG)

*A screenshot for the model, showing pedestrians in blue and vehicles in black. Note how a pedestrian agent is crossing the road at a location outside of the 'crossings' (the two small rectangle sections of the road). This behaviour is non-compliant*

![Another screenshot, this time showing compliant pedestrian behaviour.](screenshot_compliant.PNG)

*Another screenshot, this time showing compliant pedestrian behaviour.*

## Acknowledgements

This project makes use of the respastCity project by Nick Malleson: https://github.com/nickmalleson/repastcity. 

You can read more about the repastCity project (including its license) here: https://github.com/nickmalleson/repastcity/blob/master/repastcity3/documentation/intro/repastcity3.txt
/*
©Copyright 2012 Nick Malleson
This file is part of RepastCity.

RepastCity is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RepastCity is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RepastCity.  If not, see <http://www.gnu.org/licenses/>.
*/

package repastInterSim.environment.contexts;

import repast.simphony.context.DefaultContext;
import repastInterSim.environment.OD;
import repastInterSim.main.GlobalVars;

public class VehicleDestinationContext extends DefaultContext<OD> {
	
	public VehicleDestinationContext() {
		super(GlobalVars.CONTEXT_NAMES.VEHICLE_DESTINATION_CONTEXT);
	}

}

package repastInterSim;

import exceptions.DuplicateIdentifierException;
import exceptions.NoIdentifierException;

/**
 * Interface for classes that can be identified. Useful for environment objects that must read an identifier
 * value from input data.
 * @author Nick Malleson
 * 
 * Not sure if this is useful for the interSim project but including for now
 * as an example of how interfaces work.
 *
 */
public interface Identified {
	
	String getIdentifier() throws NoIdentifierException;

	void setIdentifier(String id) throws DuplicateIdentifierException;

}

/*
©Copyright 2012 Nick Malleson
This file is part of RepastCity.
*/

package exceptions;


/**
 * Exception is thrown when an object tries to set a unique identifier that has already
 * been used.
 * @author Nick Malleson
 */
public class DuplicateIdentifierException extends Exception{

	private static final long serialVersionUID = 1L;

	public DuplicateIdentifierException(String message) {
		super(message);
	}

}

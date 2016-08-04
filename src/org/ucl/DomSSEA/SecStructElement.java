/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.ucl.DomSSEA;


public class SecStructElement
{
	public char type;
	public int length;
	public double hydrophob;

	public SecStructElement( char secType, int secLength, double secHydrophob )
	{
		type = secType;
		length = secLength;
		hydrophob = secHydrophob;
      	}

    @Override
	public String toString()
	{
		return String.valueOf( type ) + "_" + String.valueOf( length )+ "_" + String.valueOf( hydrophob );
	}

	public char getType()
	{
		return type;
	}

	public int getLength()
	{
		return length;
	}

	public double getHydrophob()
	{
		return hydrophob;
	}
}

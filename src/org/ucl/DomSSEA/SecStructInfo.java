package org.ucl.DomSSEA;

import java.lang.Math.*;
import java.util.*;
/*
*object containing struct info
*/
public class SecStructInfo
{
	public Vector structVect;
	public float length;
      	public String CATcode, structStr, primStr;

	public SecStructInfo( Vector secStructVect, float secLength, String secCATcode, String secStructStr, String secPrimStr )
	{
		structVect = secStructVect;
		length = secLength;
            	CATcode = secCATcode;
            	structStr = secStructStr;
		primStr = secPrimStr;
      	}

    @Override
	public String toString()
	{
		return String.valueOf( structVect ) + "_" + String.valueOf( length ) + "_" + CATcode + "\n" + structStr + "\n" + primStr;
	}

	public Vector getStructVect()
	{
		return structVect;
	}

	public float getLength()
	{
		return length;
	}

      	public String getCATcode()
	{
		return CATcode;
	}

      	public String getStructStr()
	{
		return structStr;
	}

	public String getPrimStr()
	{
		return primStr;
	}
}

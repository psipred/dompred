/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.ucl.dompred;

import java.io.*;

/**
Reads .horiz files from psipred output and outputs a fasta style files
*/
public class psipredReader2
{
	public psipredReader2( String inputfile, String outputfile)
	throws IOException
	{
		BufferedReader in;
		in = new BufferedReader( new FileReader( inputfile ) );

		String line = in.readLine();
                String pred = "";
                String conf = "";
                String AA = "";

		do
		{       if( line.startsWith( "Conf" ) )
                                conf = conf + line.substring( 6 );
                        if( line.startsWith( "Pred" ) )
                                pred = pred + line.substring( 6 );
                        if( line.startsWith( "  AA" ) )
                                AA = AA + line.substring( 6 );

                        line = in.readLine();
		}
                while (line != null);

                DataOutputStream output = new DataOutputStream( new FileOutputStream( outputfile ) );
                output.writeBytes( ">" + outputfile + "\n" + AA + "\n" + pred );
        }

        public static void main( String args[])
	{
                try
                {
                        psipredReader2 ppr2 = new psipredReader2( args[0], args[1]);
                }
                catch( IOException e )
                {
                        System.err.println( e );
                }
	}

}

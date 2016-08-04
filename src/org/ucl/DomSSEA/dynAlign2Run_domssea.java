package org.ucl.DomSSEA;

import java.io.*;
import java.util.*;
/*
hacked of dynAlign2Run_direct to include length filter i.e. don't try to align anything too big or too small
*/
public class dynAlign2Run_domssea
{
        public dynAlign2Run_domssea( String flatfile1, String targetfile, String outputfile )
        throws IOException
        {

                File realfasta = new File( flatfile1 );
                BufferedReader realin;

                //File predfasta = new File( targetsdirect );
                //String directory[] = predfasta.list();

                //File file = new File( targetsdirect );
                //String directory[] = file.list();

                SecStructInfo ssi1, ssi2;
                String filename = "", filename2 = "";
                Hashtable predhash = new Hashtable();

                //for( int i = 0; i < directory.length; i++ )
                //{
                        //filename = directory[i];
                        BufferedReader predin = new BufferedReader( new FileReader( targetfile ) );
                        //go through second file convert strings and put into hash
                        String predline = predin.readLine();
                        predline = predin.readLine();
                        predline = predin.readLine();
                        //System.out.println( filename + "\n" + predline );
                        ssi2 = convertStringToVect( predline );
                        predhash.put( targetfile, ssi2 );
                //}

                DataOutputStream output = new DataOutputStream( new FileOutputStream( outputfile ) );

                float toptenhits[] = new float[1000];



                //hashtable for storing filenames of toptenhits
                Hashtable toptenhash = new Hashtable();

                //go through first file convert strings and put into hash

//                 String realline = realin.readLine();
//                 Hashtable realhash = new Hashtable();
//                 do
//                 {
//                         if( realline.startsWith( ">" ) )
//                                 filename = realline.substring(1);
//                         else if( !realline.startsWith( ">" ) )
//                         {
//                                 //System2.out.println( realline );
//                                 ssi1 = convertStringToVect( realline );
//                                 realhash.put( filename, ssi1 );
//                                 System.out.println( realhash.size() );
//                         }
//
//                         realline = realin.readLine();
//
//                 }
//                 while( realline != null );



                System.out.println( predhash.size() );

                //enumerate through hashtables and do all v all alignments
                Vector realvect = new Vector(), predvect = new Vector(), checkvect = new Vector();
                for( Enumeration enumpred = predhash.keys(); enumpred.hasMoreElements(); )
                {
                        String filename1 = (String)enumpred.nextElement();
                        ssi2 = (SecStructInfo)predhash.get( filename1 );
                        predvect = ssi2.getStructVect();
                        float predlength = ssi2.getLength();

                        System.out.println( filename1 );
                        int filesread = 0;
                        int filesaligned = 0;

                        //reset topscore and cuts
                        float topscore = 0, normscore = 0;
                        String tophit = "";
                        StringBuffer cutbuf = new StringBuffer();
                        int numofdomains = 0;

                        //reset topscores - fill array with blanks
                        for( int i = 0; i < 1000; i++ )
                                toptenhits[i] = 0;

                        //reset toptenhash
                        toptenhash = new Hashtable();

                        ////hack////
                        realin = new BufferedReader( new FileReader( realfasta ) );
                        String realline = realin.readLine();
                        Hashtable realhash = new Hashtable();
                        Vector realVect = new Vector(), predVect = new Vector();

                        checkvect.addElement( filename1 );

                        int i = 0;
                        do
                        {
                                //i = 0;
                                if( realline.startsWith( ">" ) )
                                        filename2 = realline.substring(1).trim();

                                //if( filename1.equals( filename2.trim() ) ) //&& predhash.containsKey( filename2 ) )
                                //        System.out.println( "start"+ filename1 + " " + filename2 + "end");

                                if( !realline.startsWith( ">" ) && !filename1.equals( filename2 ) ) //&& predhash.containsKey( filename2 ) ) //&& !checkvect.contains(filename2) )
                                {

                                        if( realline.startsWith( "E" ) || realline.startsWith( "T" )
                                                || realline.startsWith( "H" ) || realline.startsWith( "B" )
                                                || realline.startsWith( "S" ) || realline.startsWith( "G" )
                                                || realline.startsWith( "I" ))
                                        //if( !realline.startsWith( "domains" ) )
                                        {
                                                float reallength = (new Integer( realline.length())).floatValue();

                                                //calculate lengthperc
                                                //float lengthperc = predlength*(lengthbuffer/20);

                                                filesread++;

                                                //if( reallength > predlength+lengthperc )
                                                //        System.out.println( "Skipping " + filename2 + " because it is over " + lengthbuffer +  "% longer than " + filename1 );
                                                //if( reallength < predlength-lengthperc )
                                                //        System.out.println( "Skipping " + filename2 + " because it is over " + lengthbuffer +  "% shorter than " + filename1 );

                                                //System.out.println( " cutoffperc: " + lengthperc + " predlength: " + predlength + " reallength: " + reallength + " pred-perc: " + (predlength-lengthperc) + " pred+perc: " + (predlength+lengthperc) );
                                                //System.out.println( );


                                                //System.out.println( lengthperc + " " + predlength + " " + reallength + " pred-perc: " + (predlength-lengthperc) + " pred+perc: " + (predlength+lengthperc) );


                                                //length filter make sure template is not too short and not too long
                                                //if( reallength < predlength+lengthperc && reallength > predlength-lengthperc )
                                                //{
                                                        //System.out.println( realline );
                                                        ssi1 = convertStringToVect( realline );
                                                        realvect = ssi1.getStructVect();
                                                        reallength = ssi1.getLength();

                                                        dynAlign3 da2 = new dynAlign3( predvect, realvect );

                                                        filesaligned++;

                                                        //find where cut is
                                                        predVect = da2.getSvect();
                                                        realVect = da2.getTvect();

                                                        float score = da2.getscore();
                                                        normscore = ( score/( ( reallength + predlength )/2 ) );
                                                        //System.out.println( filename1 + " " + filename2 +  "\n" + realVect + "\n" + predVect + "\n" + score );

//////////////////////////////////////////////////////////////////represent alignment as a string!!!
                                                        StringBuffer realstr = new StringBuffer();
                                                        StringBuffer predstr = new StringBuffer();
                                                        ;
                                                        realstr.append('|');
                                                        predstr.append('|');
                                                        for( int v = 0; v < realVect.size(); v++ )
                                                        {
                                                                SecStructElement sse1 = (SecStructElement)realVect.elementAt(v);
                                                                SecStructElement sse2 = (SecStructElement)predVect.elementAt(v);

                                                                char c1 = sse1.getType();
                                                                int i1 = sse1.getLength();
                                                                char c2 = sse2.getType();
                                                                int i2 = sse2.getLength();

                                                                for( int n = 1; n <= Math.max( i1, i2 ); n++ )
                                                                {
                                                                        if( n <= i1)
                                                                        {
                                                                                realstr.append( c1 );
                                                                        }

                                                                        if( n > i1)
                                                                        {
                                                                                realstr.append( '_' );
                                                                        }

                                                                        if( n <= i2)
                                                                        {
                                                                                predstr.append( c2 );
                                                                        }

                                                                        if( n > i2)
                                                                        {
                                                                                predstr.append( '_' );
                                                                        }
                                                                }

                                                                realstr.append('|');
                                                                predstr.append('|');

                                                        }

                                                        //get tophit
                                                        //if( normscore > topscore  )
                                                        //{
                                                                topscore = normscore;
                                                                tophit = filename2;

                                                                realline = realin.readLine();

                                                                //find cuts
                                                                StringTokenizer linetokens = new StringTokenizer( realline );
                                                                //System.out.println( filename2 + " " + realline );
                                                                numofdomains = (new Integer( linetokens.nextToken() )).intValue();
                                                                cutbuf = new StringBuffer();

                                                                while(  linetokens.hasMoreTokens() )
                                                                {
                                                                        int cut = (new Integer( linetokens.nextToken() )).intValue();
                                                                        int lengthpred = 0;
                                                                        int lengthreal = 0;
                                                                        boolean stop = true;

                                                                        for( i = 0; i < predVect.size() && stop; i++ )
                                                                        {
                                                                                SecStructElement ssepred = (SecStructElement)predVect.elementAt(i);
                                                                                SecStructElement ssereal = (SecStructElement)realVect.elementAt(i);

                                                                                lengthpred = lengthpred + ssepred.getLength();
                                                                                lengthreal = lengthreal + ssereal.getLength();

                                                                                if( lengthreal >= cut )
                                                                                {
                                                                                        cutbuf.append( lengthpred + " " );
                                                                                        //System.out.println( p );
                                                                                        stop = false;
                                                                                }

                                                                        }
                                                                        //System.out.println( filename1 + " " + filename2 + "\n" + realVect + "\n" + predVect + "\n" + normscore );
                                                                        //System.out.println( numofdomains + " " + cut + " " + cutbuf );

                                                                }
                                                        //}
                                                        //System.out.println( i++ );

                                                        //put top 20 scores into array list
                                                      if ( normscore > toptenhits[0] )
                                                      {
                                                                 //System.out.println( normscore + " " + toptenhits[0] );
                                                                 if( toptenhash.containsKey( new Float(normscore) ) )
                                                                 {
                                                                       String oldfilename2str = (String)toptenhash.get( new Float(normscore) );
                                                                       toptenhash.put( new Float(normscore), oldfilename2str+" "+filename2 + " " + numofdomains + " " + cutbuf.toString() + " \"" + predstr.toString() + realstr.toString() + "\"" );
                                                                 }
                                                                 else if( !toptenhash.containsKey( new Float(normscore) ) )
                                                                 {
                                                                        toptenhash.remove( new Float( toptenhits[0] ) );
                                                                        toptenhash.put( new Float(normscore), filename2 + " " + numofdomains + " " + cutbuf.toString() + " \"" + predstr.toString() + realstr.toString() + "\"" );
                                                                        toptenhits[0] = normscore;
                                                                 }
                                                      }
                                                      Arrays.sort( toptenhits );

                                                        //output.writeBytes( filename1 + " " + tophit + " " + topscore + " " + numofdomains + " " + cutbuf.toString() + "\n" );
                                                //}


                                        }

                                }

                                realline = realin.readLine();
                                //System.out.println( realline );
                        }
                        while( realline != null );
                        ////hack///

//                         for( Enumeration enumreal = realhash.keys(); enumreal.hasMoreElements(); )
//                         {
//                                 String filename2 = (String)enumreal.nextElement();
//                                 ssi1 = (SecStructInfo)realhash.get( filename2 );
//                                 realvect = ssi1.getStructVect();
//                                 float reallength = ssi1.getLength();
//
//                                 dynAlign2 da2 = new dynAlign2( predvect, realvect );
//
//                                 float score = da2.getscore();
//                                 float normscore = ( score/( ( reallength + predlength )/2 ) );
//
//                                 //put topscores into array list
//                                 if( normscore > topscore )
//                                 {
//                                    topscore = normscore;
//                                    tophit = filename2;
//                                 }
//                                 System.out.println( filename1 + " " + filename2 );
//                         }

                        //output.writeBytes( filename1 + " " + tophit + " " + topscore + " " + numofdomains + " " + cutbuf.toString() + "\n" );

                        System.out.println( filename1 + " - aligned: " + filesaligned + " out of: " + filesread );

			//toptenhits reverse
                        for( int j = 999; j >= 0; j-- )
                        {
//                               String info =  (String)toptenhash.get( new Float(toptenhits[j]) );
//                               StringTokenizer infotokens = new StringTokenizer( info );
//                               while( infotokens.hasMoreTokens() )
//                               {
//                                 String infotok = infotokens.nextToken();
//
//                                 if( infotok.startsWith("|") )
//                                 {
//                                         info = info + infotok;
//                                         output.writeBytes(  toptenhits[j] + " " + info  + "\n" );
//                                         info = info + infotok;
//                                 }
//
//
//                               }
                              output.writeBytes(  toptenhits[j] + " " + (String)toptenhash.get( new Float(toptenhits[j]) )  + "\n" );
                        }
                }

                //dynAlign2 da2 = new dynAlign2( structVect, errVect );

                //float score = da2.getscore();
                //float normscore = ( score/( ( length + predlength )/2 ) );

                //System.out.println( normscore + "\n" + length + "\n" + predlength +"\n" +da2.getSvect() + "\n" + da2.getTvect() );
                //System.out.println( normscore );
        }

        public SecStructInfo convertStringToVect( String secstruct )
        {
                String line = secstruct;

                float length = ( new Float( line.length() )).floatValue();

                Vector structVect;

                char currenttype;
                char structtype = ' ';
                int l = 0;
                float tothydrophob = 0;

                StringBuffer structbuf = new StringBuffer();
                structVect = new Vector();
                structVect.addElement( new SecStructElement( '_', 0, 0.000 ) );

                for( int m = 0; m <= line.length(); m++ )
                {
                        if( m  == line.length() )
                        {
                                currenttype = ' ';
                        }

                        else
                        currenttype = line.charAt( m );

                        if( currenttype != structtype || m == line.length() )
                        {

                                if( structtype == 'C' || structtype == 'H' && l >= 1 ||
                                        structtype == 'E' && l >= 1  )
                                {
                                        for( int n = 1; n <= l; n++ )
                                        {
                                                structbuf.append( structtype );
                                        }
                                }

                                if( structtype == 'H' && l < 1 || structtype == 'E' && l < 1 || structtype == 'T' ||
                                        structtype == 'S'|| structtype == ' ' || structtype == 'I' || structtype == 'B' || structtype == 'G')
                                {
                                        for( int n = 1; n <= l; n++ )
                                        {
                                                structbuf.append( 'C' );
                                        }
                                }

                                structtype = currenttype;
                                l = 0;
                        }

                        if( currenttype == structtype )
                        {
                                l++;
                        }

                }

                line = structbuf.toString();

                //System.out.println( line + "\n" + predline );

                l = 0;
                tothydrophob = 0;



                //convert real sec struct string to vector of SecStructElements
                for( int k = 0; k <= line.length() ; k++ )
                {
                        if( k  == line.length() )
                        {
                                currenttype = ' ' ;
                        }

                        else
                                currenttype = line.charAt( k );

                        if( currenttype != structtype )
                        {

                                if( structtype != ' ' )
                                {
                                                structVect.addElement( new SecStructElement( structtype, l, -0.4 ) );
                                }

                                structtype = currenttype;
                                l=0;
                                tothydrophob = 0;
                        }

                        if( currenttype == structtype )
                        {
                                l++;
                        }


                }

                //System.out.println(  structVect );

                //System.out.println( q3 );
                //chop the ends off if they are coil regions and adjust length accordingly
//                 int startcoillength = 0;
//                 int endcoillength = 0;
//
//                 if( ( ( SecStructElement )structVect.elementAt( 1 ) ).getType() == 'C' )
//                 {
//                         startcoillength = ( ( SecStructElement )structVect.elementAt( 1 ) ).getLength();
//                         structVect.removeElementAt( 1 );
//                 }
//
//                 if( ( ( SecStructElement )structVect.elementAt( structVect.size() - 1 ) ).getType() == 'C' )
//                 {
//                         endcoillength = ( ( SecStructElement )structVect.elementAt( structVect.size() - 1 ) ).getLength();
//                         structVect.removeElementAt( structVect.size() - 1 );
//                 }
//
//                 length = length - ( startcoillength + endcoillength );

                SecStructInfo ssi = new SecStructInfo( structVect, length, new String("1"), new String("2"), new String("3") );
                //SecStructInfo ssi = new SecStructInfo( structVect, 1, new String("1"), new String("2"), new String("3") );
                return ssi;
        }

        public static void main( String args[] )
        {
                if( args.length != 3 )
                {
                        System.out.println( "usage: java -jar dynAlign2Run_domssea.jar <databank file> <targets file> <output file>" );
                }

                else if( args.length == 3 )
                {
                        try
                        {
                                dynAlign2Run_domssea da2rd = new dynAlign2Run_domssea( args[0], args[1], args[2] );
                        }
                        catch( IOException e )
                        {
                                System.out.println( e );
                        }
                }
        }
}

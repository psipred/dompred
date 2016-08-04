/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.ucl.DomSSEA;


import java.util.*;
import java.lang.Math.*;
/**
* Carries out a pairwise alignment on elements within two vectors
*
* @author Liam James McGuffin
*
*/
public class dynAlign3
{
	public float covalue, val1, val2, val3, score;
      public int gap = 0, trace;
	public Vector S1, T1, S11, T11;

	public dynAlign3( Vector S, Vector T )
	{
		//System.out.println( S );
		//System.out.println( T );

		float coarray[][];
            int tracearray[][];

		coarray = new float[ S.size() ][ T.size() ];
		tracearray = new int[ S.size() ][ T.size() ];

		//System.out.println( coarray.size() );

		for( int i = 0; i < S.size(); i++ )
		{

			for( int j = 0; j < T.size(); j++ )
			{

				if( i == 0 && j == 0 )
				{
					coarray[ i ][ j ] = 0;

					tracearray[ i ][ j ] = 0;
				}

				if( i == 0 && j != 0)
				{
					coarray[ i ][ j ] = 0;
                                    //coarray[ i ][ j ] = -j;
				      tracearray[ i ][ j ] = 2;
				}

				if( j == 0  && i != 0)
				{
					coarray[ i ][ j ] = 0;
                                    //coarray[ i ][ j ] = -i;
				      tracearray[ i ][ j ] = 3;
			      }

				if( i != 0  && j != 0 )
				{
                                    float ls = ( ( SecStructElement )S.elementAt( i ) ).getLength();
                                    float lt = ( ( SecStructElement )T.elementAt( j ) ).getLength();
                                    char s = ( ( SecStructElement )S.elementAt( i ) ).getType();
                                    char t = ( ( SecStructElement )T.elementAt( j ) ).getType();

                 if( s == t )
					  {
                                          //original scoring
                                          if( ls <= lt )
                                          score = ls;

                                          else if( lt < ls )
                                          score = lt;


                                          //scoring as abs diff in length
                                          //if( Math.abs( ls - lt ) > 0 )
                                          //      score = 1-( Math.abs( ls - lt )/Math.max( ls, lt ) );

                                          //score = 1;

                   }


				     	else if( s != t )
                  {
                                         if( s != 'C' && t != 'C' )
                                         score = 0;

                                         else
                                          {
                                                //original scoring
                                                if( ls <= lt )
                                                score = ls/2;

                                                else if( lt < ls )
                                                score = lt/2;

                                                //scoring as abs diff in length
                                         //       if( Math.abs( ls - lt ) > 0 )
                                         //             score = ( 1-( Math.abs( ls - lt )/Math.max( ls, lt ) ) )/2;

                                         }

                                         //score = 0;
					}

				      //original method
					val1 = coarray[ i-1 ][ j-1 ] + score;
					val2 = coarray[ i-1 ][ j ] + gap;
					val3 = coarray[ i ][ j-1 ] + gap;

                                  //short helices and strands are ok
                                  //if( s == 'E' && ls > 4 )
                                  //val2 = coarray[ i-1 ][ j ] - ls;

                                  //if( s == 'H' && ls > 6 )
                                  //val2 = coarray[ i-1 ][ j ] - ls;

					//if( t == 'E' && ls > 4 )
                                   //val3 = coarray[ i ][ j-1 ] - lt;

                                  //if( t == 'H' && ls > 6 )
                                  //val3 = coarray[ i ][ j-1 ] - lt;

                                  //gap penalty = length of element opposite gap
                                  //val1 = coarray[ i-1 ][ j-1 ] + score;
                                  //val2 = coarray[ i-1 ][ j ] - ls;
					//val3 = coarray[ i ][ j-1 ] - lt;

				      	if( val1 >= val2 && val1 >= val3 )
					{
						covalue = val1;
				            trace = 1;
					}

					else if( val2 >= val1 && val2 >= val3 )
					{
					     	covalue = val2;
				            trace = 3;
				      }

					else if( val3 >= val1 && val3 >= val2 )
					{
					     	covalue = val3;
					      trace = 2;
					}

					coarray[ i ][ j ] = covalue;
				      tracearray[ i ][ j ] = trace;
                                    //System.out.println( "score: " + score + "covalue: " + covalue  + "ls: " + ls + "lt: " + lt + "max: " + Math.max( ls, lt ) + " " + Math.abs( ls - lt ) + " " + (1-( Math.abs( ls - lt )/Math.max( ls, lt ) ) ) );

      				}

			}

		}

            //baseline score
            score = coarray[ S.size() -1 ][ T.size() -1 ];
            //System.out.println( "score: " + score );

            	//for( int i = 0; i < coarray.length; i++ )
		//{
		//	for( int j = 0; j < coarray[ i ].length; j++ )
		//	{
		//		System.out.println( String.valueOf( coarray[ i ][ j ] ) );
                //
		//	}
		//}


            	//for( int i = 0; i < tracearray.length; i++ )
		//{
		//	for( int j = 0; j < tracearray[ i ].length; j++ )
		//	{
		//		System.out.println( String.valueOf( tracearray[ i ][ j ] ) );
		//
		//	}
            	//}

            	//trace back

		trace = tracearray[ S.size()-1 ][ T.size()-1 ];

		//System.out.println( trace );

		int i = S.size()-1, j = T.size()-1;

            S1 = new Vector();
            T1 = new Vector();
            S11 = new Vector();
            T11 = new Vector();
		while( trace != 0 )
		{

			if( trace == 1 )
		      {
			      S1.addElement( ( SecStructElement )( S.elementAt( i ) ) );
			      T1.addElement( ( SecStructElement )( T.elementAt( j ) ) );

			      trace = tracearray[ i-1 ][ j-1 ];
			      i--;
				j--;
		      }

		      else if( trace == 3 )
		      {
				S1.addElement( ( SecStructElement )( S.elementAt( i ) ) );
				T1.addElement( new SecStructElement( '-', 0, 0.000 ) );

				trace = tracearray[ i-1 ][ j ];
				i--;
		      }

		      else if( trace == 2 )
		      {
			      S1.addElement( new SecStructElement( '-', 0, 0.000  ) );
			      T1.addElement( ( SecStructElement )( T.elementAt( j ) ) );

			      trace = tracearray[ i ][ j-1 ];
				j--;

		      }

		      //System.out.println( trace );
		      //System.out.println( S1 );
		      //System.out.println( T1 );

		}

		for( int x = S1.size()-1; x >= 0; x-- )
		{
			S11.addElement( S1.elementAt( x ) );
			//System.out.println( Sbuf.charAt( x ) );
		}

		for( int x = T1.size()-1; x >= 0; x-- )
		{
			T11.addElement( T1.elementAt( x ) );
			//System.out.println( Tbuf.charAt( x ) );
		}

              //System.out.println( S1 );
              //System.out.println( T1 );

		//System.out.println( S11 );
		//System.out.println( T11 );

	}

	public float getscore()
	{
		return score;
	}

      public Vector getSvect()
      {
            return S11;
      }

      public Vector getTvect()
      {
            return T11;
      }


}

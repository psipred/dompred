import java.io.File;
import java.io.IOException;
import org.ucl.DomSSEA.dynAlign2Run_domssea;
import org.ucl.dompred.psipredReader2;

class DomSSEA
{
  public static void main(String args[])
  {
      if( args.length != 2 )
      {
            System.out.println( "usage: java DomSSEA <File_prefix> <All_cut3 file>" );
            System.exit(1);
      }

    File fHoriz = new File(args[0] + ".horiz");
    File fPred  = new File(args[0] + ".pred");
    File fDomssea = new File(args[0] + ".domssea");
    String strAllCut = args[1];

    System.out.println("   *** DOMSSEA COMPONENTS RUNNING ***");
    try
    {
      RunPsipredReader2(fHoriz, fPred);
    }
    catch (Exception ex)
    {
      System.out.println(ex.getMessage());
    }
    try
    {
      RunDynAlign2(strAllCut, fPred, fDomssea);
    }
    catch (Exception ex)
    {
      System.out.println(ex.getMessage());
    }
  }

  /*here we initialise the blast bits and pieces as per the job config*/
  protected static void RunPsipredReader2(File fHoriz, File fPred) throws Exception
   {
      try
      {
          System.out.println("    Running PsipredReader2\n    ");
          psipredReader2 pr2 = new psipredReader2( fHoriz.getCanonicalPath(), fPred.getCanonicalPath());

          if(fPred.exists() && fPred.length() > 0)
          {
          }
          else
          {
              throw (new Exception("no pred file was created"));
          }
      }
      catch (Exception ex)
      {
          throw (new Exception(ex.getMessage()));
      }
   }

  protected static void RunDynAlign2(String strAllCut, File fPred, File fDomssea) throws Exception
   {
      try
      {
          System.out.println("    Running: DynAlign2_domssea:\n     java dynAlign2Run_domssea "+ strAllCut + " " + fPred.getCanonicalPath() + " " + fDomssea.getCanonicalPath());
          dynAlign2Run_domssea da2rd = new dynAlign2Run_domssea(strAllCut, fPred.getCanonicalPath(), fDomssea.getCanonicalPath());

          if(fDomssea.exists() && fDomssea.length() > 0)
          {
          }
          else
          {
              throw (new Exception("no domssea file was created"));
          }
      }
      catch (Exception ex)
      {
        throw (new Exception(ex.getMessage()));
      }
   }

}

package tree.technique.nj;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.ParseException;
import java.util.NoSuchElementException;
import java.util.Scanner;
import tree.interfaces.DistanceMatrix;

/**
 * A class for reading data provided by PEx dmat files.
 * Details of PEx formats can be found in http://infoserver.lcad.icmc.usp.br/infovis2/PEXImage
 * 
 * @author Guilherme P. Telles and Jose Gustavo de Souza Paiva.
 */

public class DistanceMatrixReader {

  /**
   * Loads a PEx distance matrix file.  For a n x n matrix the 
   * file format is as follows.
   * 
   *  n
   *  L[0];L[1];...;L[n-1]
   *  C[0];C[1];...;C[n-1]
   *  M[1][0]
   *  M[2][0];M[2][1]
   *  M[3][0];M[3][1];M[3][2]
   *  ...
   *  M[n-1][0];M[n-1][1];M[n-1][2];...;M[n-1][n-2]
   * 
   * L[.] are object labels and C[.] are object classes. 
   * 
   * Details can be found in http://infoserver.lcad.icmc.usp.br/infovis2/PEXImage
   * 
   * @param file The file.
   * @return a PexMatrix record, where the matrix itself is a 
   * lower triangular matrix with an all zeros main diagonal
   * @throws ParseException If file is not properly formatted.
   * @throws IOException If an IO error occurs. 
   * 
   * @author Guilherme P. Telles 
   */
  public static PexMatrix loadPex(String file) throws IOException, ParseException {

      String st = "";
    try {
      Scanner s = new Scanner(new BufferedReader(new FileReader(file)));
      PexMatrix D = new PexMatrix();
  
      D.n = s.nextInt();
      s.nextLine();
      D.labels = s.nextLine().split("\\s*;\\s*");

      String[] cdata = s.nextLine().split(";");
      D.classes = new float[D.n];
      for (int i = 0; i < D.n; i++) {
        //Double n = Double.parseDouble(cdata[i]);
        Float n = Float.parseFloat(cdata[i]);
        D.classes[i] = n;//n.intValue();
      }

      s.useDelimiter("\\s*;\\s*|;?\\s*\\n");
      //s.useDelimiter("[;\\n\\s]*");
  
//      D.classes = new int[D.n];
//      for (int i = 0; i < D.n; i++)
//        D.classes[i] = (int) s.nextDouble();

      D.M = new double[D.n][];
      D.M[0] = new double[1];
      for (int i=1; i<D.n; i++) {
        D.M[i] = new double[i+1];
  
        for (int j=0; j<i; j++) {
          st = s.next();
          //D.M[i][j] = s.nextDouble();
          D.M[i][j] = Double.parseDouble(st);
          //System.out.println("OK : "+st);
        }
  
        D.M[i][i] = 0;
      }
      s.close();
  
      return D;
    }
    catch (NoSuchElementException e) {
        //System.out.println("Erro : "+st);
      throw new ParseException("Format mismatch in "+file,0);
    }
  }

  /**
   * Loads a PEx distance matrix object.
   * 
   * * Details can be found in http://infoserver.lcad.icmc.usp.br/infovis2/PEXImage
   * 
   * @param dmat the distance matrix object
   * 
   * @author Guilherme P. Telles
   */
  public static PexMatrix loadPex(DistanceMatrix dmat) {

      PexMatrix dm = new PexMatrix();

      dm.n = dmat.getElementCount();
      dm.ids = new int[dm.n];
      dm.labels = new String[dm.n];
      dm.classes = new float[dm.n];

      //populating ids data
      if (!dmat.getIds().isEmpty())
          for (int i=0;i<dm.n;i++)
              dm.ids[i] = dmat.getIds().get(i);

      //populating class data
      if (dmat.getClassData().length > 0)
          for (int i=0;i<dm.n;i++)
              dm.classes[i] = new Float(dmat.getClassData()[i]);//new Float(dmat.getClassData()[i]).intValue();

      //populating label data
      if (!dmat.getLabels().isEmpty()) {
          for (int i=0;i<dm.n;i++)
              dm.labels[i] = dmat.getLabels().get(i);
      }else {
          //In this case, matrix elements do not have labels. Creating labels based on its ids...
          if (!dmat.getIds().isEmpty()) {
              for (int i=0;i<dm.n;i++) {
                  dm.labels[i] = dmat.getIds().get(i).toString();
              }
          } else {
              //int this very rare case, matrix elements do not have ids. Creating ids...
              for (int i=0;i<dm.n;i++) {
                  dm.labels[i] = Integer.toString(i);
              }
          }
      }
      
      dm.M = new double[dm.n][];
      dm.M[0] = new double[1];
      dm.M[0][0] = 0.0;
      for (int i=0;i<dmat.getDistmatrix().length;i++) {
          dm.M[i+1] = new double[dmat.getDistmatrix()[i].length+1];
          for (int j=0;j<dmat.getDistmatrix()[i].length;j++) {
              dm.M[i+1][j] = dmat.getDistmatrix()[i][j];
          }
          dm.M[i+1][dmat.getDistmatrix()[i].length] = 0.0;
      }
      return dm;
    
  }

}




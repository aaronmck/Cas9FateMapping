package tree.technique.nj;

/**
 * A record for data provided by PEx dmat files: the matrix dimension n; 
 * the lower triangular, all zeros main diagonal matrix M; the labels of each row in 
 * the same order of the matrix rows; the class of each row in the same order of the 
 * matrix rows.
 * 
 * @author Guilherme P. Telles.
 */

public class PexMatrix {
  public int n;
  public double[][] M;
  public int[] ids;
  public String[] labels;
  public float[] classes;
}

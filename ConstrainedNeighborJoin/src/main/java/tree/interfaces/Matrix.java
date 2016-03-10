package tree.interfaces;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Description of the methods used by the algorithms to construct a NJ tree from a Matrix.
 * 
 * Format: Consult http://infoserver.lcad.icmc.usp.br/infovis2/PExImage, in MANUAL link
 * 
 *
 * @author Jose Gustavo de Souza Paiva
 */
public interface Matrix {

    public void load(String filename) throws IOException;

    public void save(String filename) throws IOException;

    public abstract Object clone() throws CloneNotSupportedException;

    public boolean contains(Vector vector);

    public void addRow(Vector vector);

    public void addRow(Vector vector, String label);

    public void setRow(int index, Vector vector);

    public void setRow(int index, Vector vector, String label);

    public Vector removeRow(int index);

    public int getRowCount();

    public int getDimensions();

    public ArrayList<Vector> getRows();

    public Vector getRow(int row);

    public String getLabel(int row);

    public void setLabels(ArrayList<String> labels);

    public void normalize();

    public float[][] toMatrix();

    public ArrayList<String> getAttributes();

    public void setAttributes(ArrayList<String> attributes);

    public ArrayList<Integer> getIds();

    public float[] getClassData();

    public ArrayList<String> getLabels();

}

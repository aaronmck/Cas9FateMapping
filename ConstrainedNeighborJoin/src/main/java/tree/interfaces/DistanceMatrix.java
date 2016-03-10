package tree.interfaces;

import java.io.IOException;
import java.util.ArrayList;

/**
 * Description of the methods used by the algorithms to construct a NJ tree from a Distance Matrix.
 *
 * Format: Consult http://infoserver.lcad.icmc.usp.br/infovis2/PExImage, in MANUAL link
 *
 * @author Jose Gustavo de Souza Paiva
 */
public interface DistanceMatrix {

    public float[][] getDistmatrix();

    public void setDistmatrix(float[][] distmatrix);

    public float getDistance(int indexA, int indexB);

    public void setDistance(int indexA, int indexB, float value);

    public float getMaxDistance();

    public void setMaxDistance(float maxDistance);

    public float getMinDistance();

    public void setMinDistance(float minDistance);

    public int getElementCount();

    public void setElementCount(int n);

    public void save(String filename) throws IOException;

    public void load(String filename) throws IOException;

    public float[] getClassData();

    public ArrayList<Integer> getIds();

    public void setClassData(float[] cdata);

    public void setIds(ArrayList<Integer> ids);

    public void setLabels(ArrayList<String> labels);

    public ArrayList<String> getLabels();

    public void removeElement(int id);

}

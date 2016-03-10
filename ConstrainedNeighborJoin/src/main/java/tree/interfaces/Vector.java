/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package tree.interfaces;

import java.io.BufferedWriter;
import java.io.IOException;

/**
 * Description of the methods used by the Matrix object to represent an instance of the data set.
 *
 * Format: Consult http://infoserver.lcad.icmc.usp.br/infovis2/PExImage, in MANUAL link
 *
 *
 * @author Jose Gustavo de Souza Paiva
 */
public interface Vector {

    public abstract void normalize();

    public abstract float dot(Vector vector);

    public abstract float[] toArray();

    public abstract float getValue(int index);

    public abstract void setValue(int index, float value);

    public abstract void write(BufferedWriter out) throws IOException;

    public abstract Object clone() throws CloneNotSupportedException;

    public void create(float[] vector, int id, float klass);

    public float norm();

    public void shouldUpdateNorm();

    public int size();

    public int getId();

    public void setId(int id);

    public float getKlass();

    public void setKlass(float klass);

    public boolean isNull();

    /**
     * This method it is only to be used if you know what you are doing.
     * @return The values stored on this vector.
     */
    public float[] getValues();

}

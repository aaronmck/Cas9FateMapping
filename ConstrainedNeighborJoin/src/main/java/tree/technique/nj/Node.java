package tree.technique.nj;

import java.util.ArrayList;

/**
 * A structure to represent a node. This structure is used by the techniques to construct the tree.
 *
 * @author Jose Gustavo de Souza Paiva
 */
public class Node {

    public int id;
    public ArrayList<Node> children = new ArrayList<Node>();
    public ArrayList<Double> distance = new ArrayList<Double>();

    public Node(int id) {
        this.id = id;
    }

    @Override
    public String toString() {
        String ret = Integer.toString(id);
        if (!children.isEmpty()) {
            for (int i=0;i<children.size();i++)
                ret += "\n - "+children.get(i).id+":"+distance.get(i);
        }
        return ret;
    }

}

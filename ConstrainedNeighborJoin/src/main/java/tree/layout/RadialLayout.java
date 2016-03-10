package tree.layout;

import java.util.ArrayList;
import tree.basics.Edge;
import tree.basics.Tree;

/**
 * Creation of a radial layout to the tree. This methods uses the algorithm described in:
 * M. Cuadros, F. V. Paulovich, R. Minghim, and G. P. Telles. 
 * Point Placement by Phylogenetic Trees and its Application for Visual Analysis of Document Collections. "
 * In Proceedings of IEEE Symposium on Visual Analytics Science and Technology (VAST.2007), pages 99.106, Sacramento, CA, USA, 2007.
 *
 * @author Jose Gustavo de Souza Paiva
 */
public class RadialLayout {

    public float[][] execute(Tree tree) {
        
        System.out.println("Radial Layout Procedure: M. Cuadros, F. V. Paulovich, R. Minghim, and G. P. Telles." 
            + " Point Placement by Phylogenetic Trees and its Application for Visual Analysis of Document Collections."
            + " In Proceedings of IEEE Symposium on Visual Analytics Science and Technology (VAST.2007),"
            + "pages 99.106, Sacramento, CA, USA, 2007.");
        
        //initializing the tree
        init(tree);

        for (int i=0;i<nodes.size();i++)
            nodes.get(i).dist = 1.0f/(1.0f+nodes.get(i).dist);

        //creating the layout
        postorderTraversal(root);

        root.x = 0.0f;
        root.y = 0.0f;
        root.w = (float) (2 * Math.PI);
        root.t = 0;

        preorderTraversal(root);

        int size = nodes.size();

        float[][] projection = new float[size][2];

        for (int i = 0; i < size; i++) {
            float[] vect = new float[2];
            vect[0] = nodes.get(i).x;
            vect[1] = nodes.get(i).y;
            projection[i] = vect;
        }

        return projection;
    }

    private void init(Tree tree) {
        
        ArrayList<Edge> edges = tree.getEdges();
        int nrinstances = tree.getSize();

        //creating the nodes
        nodes = new ArrayList<Node>();
        for (int i = 0; i < nrinstances; i++) {
            int id = tree.getNode(i).getId();
            Node node = new Node(id);
            nodes.add(node);
            if (tree.getRootId() == id)
                root = node;
        }

        for (int i = 0;i < edges.size();i++) {
            Edge edge = edges.get(i);
            Node parent = getNode(edge.getSource());
            Node ch = getNode(edge.getTarget());
            ch.parent = parent;
            ch.dist = edge.getWeight();
            parent.children.add(ch);
        }
    }

    private void postorderTraversal(RadialLayout.Node v) {
        if (v.children.isEmpty()) {
            v.l = 1;
        } else {
            v.l = 0;

            for (RadialLayout.Node w : v.children) {
                postorderTraversal(w);
                v.l = v.l + w.l;
            }
        }
    }

    private void preorderTraversal(RadialLayout.Node v) {
        if (v != root) {
            Node u = v.parent;
            v.x = u.x + v.dist * (float) Math.cos(v.t + v.w / 2);
            v.y = u.y + v.dist * (float) Math.sin(v.t + v.w / 2);
        }

        float n = v.t;

        for (RadialLayout.Node w : v.children) {
            w.w = w.l / root.l * (float) (2 * Math.PI);
            w.t = n;
            n = n + w.w;
            preorderTraversal(w);
        }
    }

    private Node getNode(int id) {
        for (int i=0;i<this.nodes.size();i++) {
            if (this.nodes.get(i).id == id) return this.nodes.get(i);
        }
        return null;
    }

    public class Node {

        public Node(int id) {
            this.id = id;
        }

        public int id;
        public float x;
        public float y;
        public float w;
        public float t;
        public float l;
        public Node parent;
        public float dist; //distance to parent
        public ArrayList<Node> children = new ArrayList<Node>();
    }

    private ArrayList<Node> nodes;
    private RadialLayout.Node root;
}

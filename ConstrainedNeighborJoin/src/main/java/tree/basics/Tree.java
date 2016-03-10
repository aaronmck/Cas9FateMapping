package tree.basics;

import java.io.*;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Representation of a NJ tree
 * 
 * @author Jose Gustavo de Souza Paiva
 */
public class Tree {

    private String type; //Algorithm used to generate the tree
    private ArrayList<ContentTree> nodes; //list of tree nodes
    private ArrayList<Edge> edges; //list of edge nodes

    public Tree() {
        nodes = new ArrayList<ContentTree>();
        edges = new ArrayList<Edge>();
    }

    /**
     * Constructor of a tree from a newick format file
     * 
     * @param newickFile newick file
     * @throws IOException 
     */
    public Tree(String newickFile) throws IOException {
        this.nodes = new ArrayList<ContentTree>();
        edges = new ArrayList<Edge>();
        BufferedReader in = null;
        try {
            in = new BufferedReader(new java.io.FileReader(newickFile));
            if (in != null) {
                Ref r = null;
                NewickTree nw = new NewickTree(r,in.readLine());
                convertToTree(nw);
            }
        } catch (IOException e) {
            throw new IOException(e.getMessage());
        } finally {
            if (in != null) {
                try {
                    in.close();
                } catch (IOException ex) {
                    Logger.getLogger(Tree.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
    }

    /**
     * Method used by constructor that builds a tree from a newick format file
     * 
     * @param nw 
     */
    private void convertToTree(NewickTree nw) {
        getTree(nw,nw);
    }

    /**
     * Recursive Method used by constructor that builds a tree from a newick format file
     * 
     * @param nw 
     */
    private ContentTree getTree(NewickTree nw,NewickTree root) {

        if (nw != null) {
            int id;
            try {
                id = Integer.parseInt(nw.name);
            } catch (NumberFormatException nfe) {
                id = nw.name.toLowerCase().hashCode();
            }
            ContentTree ct = new ContentTree(id);
            this.nodes.add(0,ct);
            if (nw.left != null) {
                ct.setValid(false);
                ContentTree left = getTree(nw.left,root);
                left.setParent(ct.getId());
                ct.setChildrenId(0,left.getId());
                ct.setDistChildren(0,root.getDistParent(nw.left.name,root));
                ContentTree right = getTree(nw.right,root);
                right.setParent(ct.getId());
                ct.setChildrenId(1,right.getId());
                ct.setDistChildren(1,root.getDistParent(nw.right.name,root));
                if (left.getLevel() > right.getLevel())
                    ct.setLevel(left.getLevel()+1);
                else
                    ct.setLevel(right.getLevel()+1);
            }else {
                ct.setLevel(0);
            }
            return ct;
        }else return null;

    }

    @Override
    public Tree clone() {
        Tree ntree = new Tree();
        ntree.type = this.type;

        ntree.nodes = new ArrayList<ContentTree>();
        
        for (int i=0;i<this.nodes.size();i++) {
            ContentTree n = (ContentTree)this.nodes.get(i);
            ContentTree p = n.clone();
            ntree.nodes.add(p);
        }

        ntree.edges = new ArrayList<Edge>();
        ntree.edges.addAll(this.edges);

        return ntree;
    }

    /**
     * Returns the algorithm type used to build the tree
     * 
     * @return the algorithm name
     */
    public String getType() {
        return type;
    }

    /**
     * Sets the algorithm type used to build the tree
     * 
     * @param type The algorithm name
     */
    public void setType(String type) {
        this.type = type;
    }

    /**
     * Adds a node to the tree
     * 
     * @param node Node to be added.
     */
    public void addNode(ContentTree node) {
        nodes.add(node);
    }

    /**
     * Returns the tree nodes list
     * 
     * @return the nodes list
     */
    public ArrayList<ContentTree> getNodes() {
        return nodes;
    }

    /**
     * Return the number of real (valid) nodes of the tree. 
     * 
     * @return the number of real (valid) nodes.
     */
    public int getNumberRealNodes() {
        int result = 0;
        for (int i=0;i<nodes.size();i++) {
            if (nodes.get(i).isValid())
                result++;
        }
        return result;
    }

    /**
     * Verifies if a node belongs to the tree.
     * 
     * @param node node to be verified
     * @return true if the node belongs to the tree, false if not.
     */
    public boolean isNode(ContentTree node) {
        for (int i=0;i<nodes.size();i++)
            if (nodes.get(i).equals(node)) return true;
        return false;
    }

    /**
     * Returns a node of a position index on the tree nodes list
     * 
     * @param index Index of the node on nodes list
     * @return The node on that position of the nodes list, null if this position does not exist, or if there is no node on this position
     */
    public ContentTree getNode(int index) {
        return nodes.get(index);
    }

    /**
     * Return the position, on the tree nodes list, of a node.
     * 
     * @param node The node to be searched on the tree nodes list
     * @return the position on the tree nodes list, -1 if the node is not on the list
     */
    public int getNodePosition(ContentTree node) {
        return nodes.indexOf(node);
    }

    /**
     * Returns a node of the tree, given its id
     * 
     * @param id The id of the node to be returned
     * @return the node that has the given id, null if there is no node with the given id.
     */
    public ContentTree getNodeById(Integer id) {
        for (int i=0;i<nodes.size();i++) {
            if (nodes.get(i).getId() == id) return nodes.get(i);
        }
        return null;
    }

    /**
     * Remove a node from the tree node list. This method does not make any adjustment on the tree after the exclusion of the node..
     * 
     * @param node Node to be excluded from the nodes list
     */
    public void removeNode(ContentTree node) {
        nodes.remove(node);
        //TODO verify if the node has children, and decide what to do in this case
    }

    /**
     * Returns the maximum level of the nodes of the tree. This values corresponds of the root level.
     * 
     * @return the maximum level of the nodes of the tree.
     */
    public int getMaxNivel() {
        int maxNivel = nodes.get(0).getLevel();
        for (int i=1;i<nodes.size();i++)
            if (nodes.get(i).getLevel() > maxNivel)
                maxNivel = nodes.get(i).getLevel();
        return maxNivel;
    }

    /**
     * Returns the total number of nodes on the tree.
     * 
     * @return The total number of nodes on the tree.
     */
    public int getSize() {
        return nodes.size();
    }

    /**
     * Returns the edges list of the tree.
     * 
     * @return The edges list of the tree.
     */
    public ArrayList<Edge> getEdges() {
        return edges;
    }

    /**
     * Set the edges list of the tree
     * 
     * @param edges The edges list of the tree
     */
    public void setEdges(ArrayList<Edge> edges) {
        this.edges = edges;
    }
    
    /**
     * Returns the root id of the tree
     * 
     * @return The root id of the tree
     */
    public int getRootId() {
        int id = 0;
        int maxNivel = nodes.get(0).getLevel();
        for (int i=1;i<nodes.size();i++)
            if (nodes.get(i).getLevel() >= maxNivel) {
                maxNivel = nodes.get(i).getLevel();
                id = nodes.get(i).getId();
            }
        return id;
    }

    /**
     * Verifies if there is a edge on the tree connecting two given nodes ids
     * 
     * @param sour the id of the source of the edge
     * @param targ the id of the target of the edge
     * @return true if there is an edge connecting the source id to the target id, false if not.
     */
    private boolean existEdge(int sour, int targ) {
        for (int i=0;i<edges.size();i++) {
            if (((edges.get(i).getSource()==sour)&&(edges.get(i).getTarget()==targ))||
                ((edges.get(i).getSource()==targ)&&(edges.get(i).getTarget()==sour)))
                return true;
        }
        return false;
    }

    /**
     * Generate the edges list, based on the hierarchy contained in the nodes list.
     * 
     */
    public void generateEdges() {
        int sour, targ;
        float dis;
        //generating the connectivity
        for (int i=0;i<this.getSize();i++) {
            ContentTree son = this.getNode(i);
            if (son != null) {
                ContentTree parent = this.getNodeById(son.getParent());
                if (parent != null) {
                    sour = parent.getId();
                    targ = son.getId();
                    int ind = parent.getChildrenIndex(son.getId());
                    dis = 0;
                    if (ind != -1)
                        dis = parent.getDistChildren(ind);
                    else
                        dis = 1000; //This distance is just a dummy distance, for the layout algorithm procedure, and is not related to any distance.
                    if (!existEdge(sour,targ)) {
                        edges.add(new Edge(sour, targ, dis));
                    }
                }
            }
        }
    }

    /**
     * Writes the list of nodes of the tree, with its attributes.
     */
    public void printNodes() {
        System.out.println("\n\n--- NODES ---");
        for (int i=0;i<nodes.size();i++) {
            System.out.print("ID:"+nodes.get(i).getId()+":"+
                             "CLASS:"+nodes.get(i).getKlass()+":"+
                             "NIVEL:"+nodes.get(i).getLevel()+":"+
                             "PARENT:"+nodes.get(i).getParent()+" | CHILDREN:");
            if (nodes.get(i).hasChildren())
                for (int j=0;j<nodes.get(i).getNumChildren();j++)
                    System.out.print(" "+nodes.get(i).getChildrenId(j)+"("+nodes.get(i).getDistChildren(j)+")");
            else
                System.out.print("NONE");
            System.out.println(" ");
        }
        System.out.println("-------------\n\n");
    }

    /**
     * Writes the list of edges of the tree, with its attributes.
     */
    public void printEdges() {
        System.out.println("\n\n--- EDGES ---");
        for (int i=0;i<edges.size();i++) {
            ContentTree source = this.getNodeById(edges.get(i).getSource());
            ContentTree target = this.getNodeById(edges.get(i).getTarget());
            System.out.println("SOURCE : "+source.getId()+" | TARGET : "+target.getId()+ " ("+edges.get(i).getWeight()+")");
        }
        System.out.println("-------------\n\n");
    }

    /**
     * Saves the tree hierarchy on a XML file
     * The format of the XML is as follows:
     * 
     * <tree>
     * <node label=XXXX distance=YYYY> Nodes with children (The distance parameter represents the distance from its parent, if applicable)
     *  <leaf label=XXXX distance=YYYY /> nodes without children (The distance parameter represents the distance from its parent, if applicable)
     * </node>
     * </tree>
     * 
     * @param output the XML file name that will contain the tree structure
     */
    public void saveXML(String output) {

        ContentTree lroot = getNode(getSize()-1);
        if (lroot == null) return;
        else {
            String ret = "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" ?>"+ System.getProperty("line.separator");
            ret = "<tree>" + System.getProperty("line.separator") + getXMLNode(lroot,-1) + "</tree>";
            if (output != null && !output.isEmpty()) {
                try{
                    if (!output.endsWith(".xml"))
                        output += ".xml";
                    BufferedWriter bw = new BufferedWriter(new FileWriter(new File(output)));
                    bw.write(ret);
                    bw.flush();
                    bw.close();
                } catch(IOException e){
                    e.printStackTrace();
                }
            }
        }

    }

    /**
     * Recursive method to return a node attributes
     * 
     * @param ct the node to be returned
     * @param distance the distance of this node to its parent, -1 if the node has no parent (root)
     * @return The string representation of the node
     */
    private String getXMLNode(ContentTree ct, float distance) {
        String ret = "";
        if (ct.hasChildren()) {
            if (distance != -1)
                ret += "<node label=\""+ct.getId()+"\" parentDistance=\""+distance+"\">"+System.getProperty("line.separator");
            else
                ret += "<node label=\""+ct.getId()+"\">"+System.getProperty("line.separator");
            for (int i=0;i<ct.getNumChildren();i++) {
                ContentTree n = getNodeById(ct.getChildrenId(i));
                if (n != null)
                    ret += getXMLNode(n,ct.getDistChildren(i));
                else
                    System.out.println("Node "+ct.getChildrenId(i)+" do not exist.");
            }
        ret += "</node>"+System.getProperty("line.separator");
        }else {
            if (distance != -1)
                ret += "<leaf label=\""+ct.getId()+"\" parentDistance=\""+distance+"\"/>"+System.getProperty("line.separator");
            else
                ret += "<leaf label=\""+ct.getId()+"\"/>"+System.getProperty("line.separator");
        }
        return ret;
    }
}
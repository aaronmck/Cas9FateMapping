package tree.basics;

public class ContentTree {

    private int id;
    private boolean valid;
    private int level; //0 leaf, > 0 non-leaf
    private int[] children;
    private float[] distChildren;
    private int parent;
    private float cdata;
    
    /**
     * Constructor of the ContentTree (a node of a tree).
     *
     * @param id the id of the ContentTree
     *
     * @author Jose Gustavo de Souza Paiva
     */
    public ContentTree(int id) {
        this.id = id;
        this.valid = true;
        this.level = 0;
        this.parent = -1;
        //this.distChildren = new float[2];
        //this.children = new int[2];
        //for (int i=0;i<this.children.length;i++)
        //    this.children[i] = -1;
    }

    @Override
    public ContentTree clone() {
        ContentTree t = new ContentTree(this.id);
        t.valid = this.valid;
        t.level = this.level;
        System.arraycopy(this.children, 0, t.children, 0, this.children.length);
        System.arraycopy(this.distChildren, 0, t.distChildren, 0, this.distChildren.length);
        t.parent = this.parent;
        return t;
    }

    /**
     * Returns the level, in the tree, of the node.
     * Level 0: leaf, > 0: non-leaf
     *
     * @return level of the node
     *
     */
    public int getLevel() {
        return level;
    }

    /**
     * Sets the level, in the tree, of a node
     * This method recursively increases, by one, the levels of the ancestors of this node
     *
     * @param l the level to be associated with the node
     *
     *
     */
    public void setLevel(int l) {
        this.level = l;
    }

    /**
     * Verify if the node is a real node of a virtual node
     *
     * @return true if the node is real, false if the node is virtual
     *
     */
    public boolean isValid() {
        return valid;
    }

    /**
     * Sets the validity of this node
     *
     * @param valid the value to be associated with the validity. True to real node, False to virtual node.
     *
     */
    public void setValid(boolean valid) {
        this.valid = valid;
    }

    /**
     * Returns the id of the node.
     *
     * @return the id of the node
     *
     */
    public int getId() {
        return id;
    }

    /**
     * Sets the id of a node
     *
     * @param id the id to be associated with the node
     *
     */
    public void setId(int id) {
        this.id = id;
    }

    /**
     * Returns the parent of the node.
     *
     * @return the id of the parent of the node
     *
     */
    public int getParent() {
        return parent;
    }
    
    /**
     * Sets the parent of a node
     *
     * @param parent the parent node to be associated with the node
     *
     */
    public void setParent(int parent) {
        this.parent = parent;
    }

    /**
     * Verify if the node is an ancestor of a descendant of a specific node
     *
     * @param node the node to verify
     * @return true if the node is related with the specified node, false if not
     *
     */
    public boolean isRelative(ContentTree node) {
        return false;
    }

    /**
     * Returns the root of the tree in which this node is inserted.
     *
     * @return the node that represents the root.
     *
     */
    public ContentTree getRoot() {
        return null;
    }

    /**
     * Verify if the node has children
     *
     *
     * @return true if the node has children, false if not
     *
     */
    public boolean hasChildren() {
        if (children != null) {
            for (int i=0;i<children.length;i++)
                if (children[i] != -1) return true;
            return false;
        }else
            return false;
    }

    /**
     * Returns the number of children of a node
     *
     *
     * @return the number of children
     *
     */
    public int getNumChildren() {
        int ct = 0;
        if (children != null)
            for (int i=0;i<children.length;i++)
                if (children[i] != -1)
                    ct++;
        return ct;
    }

    /**
     * Sets the distance of the ith children
     *
     * @param i the index of the children
     * @param d the distance value to be given
     *
     */
    public void setDistChildren(int i, float d) {
        if (i >= 0) {
            if (distChildren == null)
                distChildren = new float[2];
            if (i < distChildren.length)
                distChildren[i] = d;
            else { //root node. Temporary code.
                float[] temp = new float[distChildren.length];
                System.arraycopy(distChildren,0,temp,0,distChildren.length);
                distChildren = new float[3];
                System.arraycopy(temp,0,distChildren,0,temp.length);
                distChildren[2] = d;
            }
        }
    }

    /**
     * Returns the distance of this node to its ith children
     *
     * @param i the id of the children
     * @return the distance of the children
     *
     */
    public float getDistChildren(int i) {
        if (distChildren != null) {
            if (i >= 0 && i < distChildren.length)
                return distChildren[i];
            else
                return 0;
        }
        else return 0;
    }

    /**
     * Sets the id of the ith children
     *
     * @param i the index of the children
     * @param id the id to be given
     *
     */
    public void setChildrenId(int i, int id) {
        if (i >= 0) {
            if (children == null) {
                children = new int[2];
                for (int k=0;k<children.length;k++)
                    children[k] = -1;
            }
            if (i < children.length)
                children[i] = id;
            else { //root node. Temporary code.
                int[] temp = new int[children.length];
                System.arraycopy(children,0,temp,0,children.length);
                children = new int[3];
                System.arraycopy(temp,0,children,0,temp.length);
                children[2] = id;
            }
        }
    }

    /**
     * Returns the id of the ith children of a node
     *
     * @return the id of the ith children
     *
     */
    public int getChildrenId(int i) {
        if (children != null) {
            if (i >= 0 && i < children.length)
                return children[i];
            else
                return -1;
        }else return -1;
    }

    /**
     * Returns the index, on the list of children, of a given children id
     *
     * @param id the id of the children
     * @return the index of the children
     *
     */
    public int getChildrenIndex(int id) {
        if (children != null) {
            for (int i=0;i<children.length;i++)
                if (children[i] == id)
                    return i;
            return -1;
        }else return -1;
    }

    /**
     * Returns the class data of the node.
     *
     * @return the class data of the node
     *
     */
    public float getKlass() {
        return cdata;
    }

    /**
     * Sets the class data of a node
     *
     * @param c the class value of a node
     *
     */
    public void setKlass(float c) {
        cdata = c;
    }

    @Override
    public String toString() {
        return id+":("+(valid ? "valid":"non-valid")+")";
    }

}
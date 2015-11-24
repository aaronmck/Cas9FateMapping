/*
 Copyright (c) 2002 Compaq Computer Corporation

 SOFTWARE RELEASE

 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:

 - Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.

 - Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.

 - Neither the names of Compaq Research, Compaq Computer Corporation
 nor the names of its contributors may be used to endorse or promote
 products derived from this Software without specific prior written
 permission.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 IN NO EVENT SHALL COMPAQ COMPUTER CORPORATION BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

package main.java;

import java.io.IOException;
import java.util.*;

/**
 * A class representing a node of a (phylognenetic) tree. The tree that this
 * node belongs to is of type Tree. Nodes have fields that store a pre- and
 * post-ordering.
 *
 * A TreeNode has a list of children, a unique key, a leftmostleaf and a
 * rightmost leaf
 *
 * @author Tamara Munzner, Li Zhang, Yunhong Zhou
 * @version 2.2
 * @see Tree
 * @see GridCell
 */
public class TreeNode {

    /** Array of child nodes that are attached below this internal node.  Null if this is a leaf. */
    protected ArrayList children; // eventually turn this into an array (need
    // to change parser)

    /** key is unique for nodes in one tree.  Keys are pre-ordered (root = 0, depth-traversal ordering). */
    public int key;

    /** Height of font in font points used to draw the label. */
    private int fontSize;

    /** Score for a node in [0,1] that corresponds to the topological similarity between two tree drawers.
     @see TreePairs#getBestCorrNodeScore(Tree, TreeNode, Tree, int) */
    private Double bcnScore;

    // /**
    // * The offset of the point with respect to the cell. We only have
    // * this for the row offset as we assume that the vertical edges
    // * are all aligned. When computing the Y coordinate of a node, we
    // * add nodeOffsetR to the pointOffset[1], a fixed parameter set by
    // * AccordionDrawer.
    // */
    /**
     * The last frame that had a computed {@link #midYPosition}, for caching.
     */
    protected int computedFrame; // store frame midYPosition was last
    // calculated (needed to place parents)

    /** Cached location (world-space) of the mid point in the vertical of a cell where the horizontal tree edge is drawn.
     * This is (1/2 of cell size + minY) for leaves, midway between first and last child edge for internal nodes. */
    private double midYPosition;

    /**  Returns the minimum key value of nodes in the subtree rooted by this node.
     * @return The index of the smallest descendant node (which is the key for this node). */
    // this is the key for this node
    public int getMin() {
        return key;
    }

    /** Returns the maximum key value of nodes in the subtree rooted but this node.
     * @return The index of the smallest descendant node (which is the key for the rightmost leaf node). */
    // this is the key of the rightmost leaf
    public int getMax() {
        return rightmostLeaf.key;
    }

    /**
     * Returns the key for this node.
     * @return The value of {@link #key} for this node.
     */
    public int getKey() {
        return key;
    }

    /**
     * Returns the label for this node, which is {@link #name}.
     * @return The value of {@link #name} for this node.
     */
    public String getName() {
        return name;
    }

    /**
     * Tests to see if this node has a vertical or horizontal edge component.
     * @param xy 0/X for horizontal, 1/Y for vertical nodes.
     * @return True if this node has an edge in the chosen direction.  Only root nodes don't have a horizontal edge, and only leaves don't have vertical edges.
     */
    protected boolean getEdge(int xy) {
        if (xy == 0)
            return !isRoot();
        else
            return !isLeaf();
    }

    /** Implements Comparable interface - sorts on key field.
     * @param o The other object to compare this node to.
     * @return -1 if this is smaller than the object's key, +1 if greater, 0 if equal. */
    public int compareTo(Object o) {
        if (key == ((TreeNode) o).key)
            return 0;
        else if (key < ((TreeNode) o).key)
            return -1;
        else
            return 1;
    }

    /** The parent of this node.  This is null for the root node. */
    public TreeNode parent;

    /**
     * Node name with default "". Most internal nodes have no name and all leaf
     * nodes have a name.  This becomes the long version of the node name when fully
     * qualified names are used.
     */
    protected String name = ""; // the long form in fully qualified names

    /** The text that appears when the node is highlighted or has a name displayed. */
    public String label = ""; // always short form

    /** Distance from this node to the root node. The root is at height 1. */
    public int height;

    /** Weight is the horizontal edge length for the edge immediately above the node.  Edge lengths are not determined by this number currently; all edges are stretched to make leaves right aligned, with minimal integral lengths. */
    public float weight = 0.0f;

    /** Leftmost (minimum) leaf node under this internal node (or this node for leaves). */
    public TreeNode leftmostLeaf;
    /** Rightmost (maximum) leaf node under this internal node (or this node for leaves). */
    public TreeNode rightmostLeaf;

    /** The number of leaves under this internal node (or 1 for leaves). */
    public int numberLeaves;

    /** The next preorder node. */
    public TreeNode preorderNext = null;

    /** The next postorder node. */
    public TreeNode posorderNext = null;

    /**
     * Default tree node constructor.
     * Children list initially set to capacity 2 as in most case binary.
     * 	Used in 2 places: create the root when creating the tree;
     *  the parser uses this to create nodes attached to the root.
     */
    public TreeNode() {
        children = new ArrayList(2);
        bcnScore = new Double(0.0);
    }

    /**
     * Clean this node of children.
     */
    public void close() {
        children.clear();
    }

    /**
     * Destroy this node.  Runs {@link #close()}.
     */
    protected void finalize() throws Throwable {

        try {
            close();
        } finally {
            super.finalize();
            // System.out.println("finally clean treeNodes");
        }
    }

    /**
     * Set the name for this node, the name is usually the label drawn with this node.
     * @param s The new value of {@link #name}, the name for this node.
     */
    public void setName(String s) {
        name = s;
    }

    /**
     * Get the number of children under this node.
     * @return Number of nodes stored in the children array {@link #children}.
     */
    public int numberChildren() {
        return children.size();
    }

    /**
     * Get a given child for this node, with range checking and casting.
     * @param i The child index to get.
     * @return The i(th) child for this node.
     */
    public TreeNode getChild(int i) {
        if (i < children.size())
            return (TreeNode) children.get(i);
        else
            return null;
    }

    /**
     * Tests to determine if this node is a leaf.  Does not work for nodes not in the tree structure.
     * @return True if this node has no linked children, and therefore is a leaf node for the tree.
     */
    public boolean isLeaf() {
        return children.isEmpty();
    }

    /**
     * Tests to determine if this node is the root of its tree. Does not work for nodes not in the tree structure.
     * @return True if this node has no linked parent, and therefore is the root of the tree.
     */
    public boolean isRoot() {
        return (null == parent);
    }

    /**
     * Tests nodes for equality, based on the name of the node.
     * @param n Second node to test vs. this node.
     * @return True if the names of both nodes are the same, false otherwise.
     */
    public boolean equals(TreeNode n) {
        return (name.equals(n.name));
    }

    /**
     * Add a child to the end of the list of children.  Note there is no remove child method, this is permanent.
     * Additional processing for linking nodes (setting up pointers and leaf properties, for example) is done later.
     * @param n New child node for this node.
     */
    public void addChild(TreeNode n) {
        children.add(n);
        n.parent = this;
    }
    /**
     * Get the parent for this node.
     * @return Value of {@link #parent}.
     */
    public TreeNode parent() {
        return parent;
    }

    /**
     * Set the weight of this treenode, which encodes the length of the horizontal edge.
     * Edge weights are not implemented currently for drawing.
     * @param w New edge weight for this node, {@link #weight}.
     */
    public void setWeight(double w) {
        weight = (float) w;
    }

    /**
     * Get the weight of this treenode, which encodes the length of the horizontal edge.
     * Edge weights are not implemented currently for drawing.
     * @return Edge weight for this node, {@link #weight}.
     */
    public float getWeight() {
        return weight;
    }

    /** Get the first child of this node. Doesn't work with leaf nodes.
     * @return First child of this internal node.
     */
    protected TreeNode firstChild() {
        return (TreeNode) children.get(0);
    }

    /** Get the last child of this node. Doesn't work with leaf nodes.
     * @return Last child of this internal node.
     */
    public TreeNode lastChild() {
        return (TreeNode) children.get(children.size() - 1);
    }


    /**
     * Long form printing for a single node. Used in conjunction with
     * {@link #printSubtree()} to display a whole subtree.
     *
     */
    public void print() {
        if (name != null)
            System.out.print("node name: " + name + "\t");
        else
            System.out.print("node name null,\t");
        System.out.println("key: " + key);
    }

    /**
     * For debugging, prints the subtree contents, recursive.
     *
     */
    private void printSubtree() {
        print();
        for (int i = 0; i < children.size(); i++)
            getChild(i).printSubtree();
    }

    /**
     * Set the extreme leaves for this node.  This is done in leaf->root direction, so all linking can be done in O(n) time.
     *
     */
    public void setExtremeLeaves() {
        if (isLeaf()) {
            leftmostLeaf = this;
            rightmostLeaf = this;
            return;
        }
        leftmostLeaf = firstChild().leftmostLeaf;
        rightmostLeaf = lastChild().rightmostLeaf;
    }

    /** root->leaf traversal, depth first in direction of leftmost leaf. */
    public void linkNodesInPreorder() {
        if (isLeaf())
            return;
        preorderNext = firstChild();
        for (int i = 0; i < numberChildren() - 1; i++)
            getChild(i).rightmostLeaf.preorderNext = getChild(i + 1);
        // rightmostLeaf.preorderNext = null; // redundant
    }

    /** Leaf->root traversal, starting at leftmost leaf of tree. */
    public void linkNodesInPostorder() {
        if (isLeaf())
            return;
        // n.posorderNext = null; // redundant
        for (int i = 0; i < numberChildren() - 1; i++)
            getChild(i).posorderNext = getChild(i + 1).leftmostLeaf;
        lastChild().posorderNext = this;
    }

    /**
     * Sets the number of leaves, must be run on leaves first (pre-order)
     *
     * @return The number of leaves ({@link #numberLeaves}) including the
     *         current node (leaves = 1)
     */
    public int setNumberLeaves() {
        numberLeaves = 0;
        if (isLeaf())
            numberLeaves = 1;
        else
            for (int i = 0; i < children.size(); i++)
                numberLeaves += getChild(i).numberLeaves;
        return numberLeaves;
    }

    /**
     * String value of this node, name + key + tree height information.
     * @return String representation of this node.
     */
    public String toString() {
        // String edge[] = {edges[0]!=null?edges[0].toString():"X",
        // edges[1]!=null?edges[1].toString():"Y"};
        return name + "(" + key + " @ " + height + ")";
    }

    /**
     * Set the {@link #bcnScore} for this node.
     * @param n New value of {@link #bcnScore}.
     */
    public void setBcnScore(float n) {
        bcnScore = new Double(n);
    }

    /**
     * Get the BCN score for this treenode.
     * @return Value of {@link #bcnScore} for this node.
     */
    public Double getBcnScore() {
        return bcnScore;
    }
}
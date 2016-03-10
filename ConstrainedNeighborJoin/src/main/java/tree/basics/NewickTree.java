/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package tree.basics;

/**
 * Newick tree file parser
 * 
 * @author Jose Gustavo. Adapted from Eddie Yee Tak Ma, Recursive Descent newick tree parser,
 * found in http://eddiema.ca/2010/06/25/parsing-a-newick-tree-with-recursive-descent/.
 */
public class NewickTree {
    
    NewickTree left,right;
    String name;
    float branch_length;

    public NewickTree() {
    }

    protected NewickTree(Ref ref, String newick_string) {

        System.out.println("Recursive Descent Newick Tree Parser, adapted from Eddie Yee Tak Ma, in http://eddiema.ca/2010/06/25/parsing-a-newick-tree-with-recursive-descent/");
        
        try {
            if(ref == null) {
                ref = new Ref();
                ref.cursor = 0;
            }

            if (newick_string.charAt(ref.cursor) == '(') {
                ref.cursor++;
                //== 1 == Another open parenthesis...
                if (newick_string.charAt(ref.cursor) == '(') {
                    this.left = new NewickTree(ref, newick_string);
                }else {
                //== 2 == The name of a node on the LEFT side.
                    String lname = newick_string.substring(ref.cursor,newick_string.indexOf(":",ref.cursor));
                    if (lname != null && !lname.isEmpty()) {
                        this.left = new NewickTree();
                        this.left.name = lname;
                        ref.cursor += lname.length();
                    }
                }
            }
            String branch_length_string = "";
            if (newick_string.charAt(ref.cursor) == ':') {
                ref.cursor++;
                //The next value is a float.
                if (newick_string.indexOf(",",ref.cursor) != -1)
                    branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(",",ref.cursor));
                else if (newick_string.indexOf(")",ref.cursor) != -1)
                    branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(")",ref.cursor));
                try {
                    this.left.branch_length = Float.parseFloat(branch_length_string);
                } catch (NumberFormatException e) {
                    if (newick_string.indexOf(")",ref.cursor) != -1)
                        branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(")",ref.cursor));
                    else if (newick_string.indexOf(",",ref.cursor) != -1)
                        branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(",",ref.cursor));
                    try {
                        this.left.branch_length = Float.parseFloat(branch_length_string);
                    } catch (NumberFormatException e1) {
                        System.out.println("Could not parse string ["+branch_length_string+"]");
                        this.left.branch_length = 0.0f;
                    }
                }
                ref.cursor += branch_length_string.length();
            }

            if (newick_string.charAt(ref.cursor) == ',') {
                ref.cursor++;
                //== 1 == Another open parenthesis...
                if (newick_string.charAt(ref.cursor) == '(') {
                    this.right = new NewickTree(ref,newick_string);
                }else {
                //== 2 == The name of a node on the RIGHT side.
                    String rname = newick_string.substring(ref.cursor,newick_string.indexOf(":",ref.cursor));
                    if (rname != null && !rname.isEmpty()) {
                        this.right = new NewickTree();
                        this.right.name = rname;
                        ref.cursor += rname.length();
                    }
                }
            }

            if (newick_string.charAt(ref.cursor) == ':') {
                ref.cursor++;
                //The next value is a float.
                if (newick_string.indexOf(")",ref.cursor) != -1)
                    branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(")",ref.cursor));
                else if (newick_string.indexOf(",",ref.cursor) != -1)
                    branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(",",ref.cursor));
                try {
                    this.right.branch_length = Float.parseFloat(branch_length_string);
                } catch (NumberFormatException e) {
                    if (newick_string.indexOf(",",ref.cursor) != -1)
                        branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(",",ref.cursor));
                    else if (newick_string.indexOf(")",ref.cursor) != -1)
                        branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(")",ref.cursor));
                    try {
                        this.right.branch_length = Float.parseFloat(branch_length_string);
                    } catch (NumberFormatException e1) {
                        System.out.println("Could not parse string ["+branch_length_string+"]");
                        this.right.branch_length = 0.0f;
                    }
                }
                ref.cursor += branch_length_string.length();
            }

            if (newick_string.charAt(ref.cursor) == ')') {
                ref.cursor++;
            }

            String pname = "";
            if (newick_string.indexOf(":",ref.cursor) != -1)
                pname = newick_string.substring(ref.cursor,newick_string.indexOf(":",ref.cursor));
            else if (newick_string.indexOf(";",ref.cursor) != -1)
                pname = newick_string.substring(ref.cursor,newick_string.indexOf(";",ref.cursor));

            if (pname != null && !pname.isEmpty()) {
                this.name = pname;
                ref.cursor += pname.length();
            }

            if (newick_string.charAt(ref.cursor) == ':') {
                ref.cursor++;
                //The next value is a float.
                if (newick_string.indexOf(",",ref.cursor) != -1)
                    branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(",",ref.cursor));
                else if (newick_string.indexOf(")",ref.cursor) != -1)
                    branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(")",ref.cursor));
                try {
                    this.branch_length = Float.parseFloat(branch_length_string);
                } catch (NumberFormatException e) {
                    if (newick_string.indexOf(")",ref.cursor) != -1)
                        branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(")",ref.cursor));
                    else if (newick_string.indexOf(",",ref.cursor) != -1)
                        branch_length_string = newick_string.substring(ref.cursor,newick_string.indexOf(",",ref.cursor));
                    try {
                        this.branch_length = Float.parseFloat(branch_length_string);
                    } catch (NumberFormatException e2) {
                        System.out.println("Could not parse string ["+branch_length_string+"]");
                        this.branch_length = 0.0f;
                    }
                }
                ref.cursor += branch_length_string.length();
            }

            if (newick_string.charAt(ref.cursor) == ';') {
                ref.cursor++;
                if(this.right == null) {
                    if(this.left != null) {
                        this.left = this.left.left;
                        this.right = this.left.right;
                        this.branch_length = this.left.branch_length;
                    }
                }
            }
        } catch (Exception e1) {
            System.out.println("Could not parse string");
        }
    }

    public float getDistParent(String name, NewickTree n) {
        if (n.name.equalsIgnoreCase(name))
            return n.branch_length;
        else
            if (n.left != null) {
                float t = getDistParent(name,n.left);
                if (t != -1)
                    return t;
                else
                    return getDistParent(name,n.right);
            } else
                return -1;
    }
}
class Ref {
    int cursor;
}

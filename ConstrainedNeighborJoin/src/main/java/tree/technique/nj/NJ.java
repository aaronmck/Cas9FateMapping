package tree.technique.nj;

import java.io.IOException;
import java.text.ParseException;
import java.util.*;
import javax.swing.JOptionPane;
import tree.basics.ContentTree;
import tree.basics.Tree;
import tree.interfaces.DistanceMatrix;

/**
 * Implementations of NJ.
 * 
 * Every implementation takes a lower triangular matrix. For every i and j such that 0<=j<=i<n, Dij 
 * should be the distance between objects i an j.  It is unimportant whether or not Dii 
 * equals zero for every i<n.  
 * 
 * For instace, for n=4 a valid matrix is:
 * 
 * {{ 0 },
 *  { 7, 0 },
 *  { 8, 5, 0 },
 *  {11, 8, 5, 0 }}
 *  
 * @author Guilherme P. Telles (techniques implementation), Jose Gustavo S. Paiva (integration implementation).
*/
public class NJ {

    public enum NJAlgorithmType {

        ORIGINAL("Original Neighbor-Joining"),
        FAST("Fast Neighbor-Joining"),
        RAPID("Rapid Neighbor-Joining");

        private NJAlgorithmType(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        private final String name;
    }

    private boolean promotion = false;
    private NJAlgorithmType type;
    
    /**
    * Constructor of the NJ tree generator
    *
    * @param pnj True if the promotion leaf procedure has to be executed after the tree construction, false if not
    * @param type Algorithm to be used (Fast, Original or Rapid)
    *

    */
    public NJ(boolean pnj, NJAlgorithmType type) {
        this.promotion = pnj;
        this.type = type;
    }

    /**
    * NJ as proposed by Saitou and Nei.
    *
    * @param p A lower triangular matrix.
    * @return Hashmap containing the nodes of the constructed tree.
    */
    public static HashMap<Integer,Node> originalNJ(PexMatrix p) {
        
        System.out.println("ORIGINAL NJ: N. Saitou and M. Nei. "
            + "The Neighbor-Joining Method: A New Method for Reconstructing Phylogenetic Trees. "
            + "Molecular Biology and Evolution, 4(4):406–425, 1987");
        
        double Sum = 0;  // The distances sum.
        double[] sum = new double[p.n];  // The distances sum for node i:
        double[] joined = new double[p.n]; // Half of the distance between two nodes joined into a new node.
        String[] newick = new String[p.n];
        HashMap<Integer,Node> nodes = new HashMap<Integer,Node>();
        int nextId = p.n;
        int maxId = Integer.MIN_VALUE;

        for (int k=0; k<p.n; k++) {
            newick[k] = p.ids == null ? Integer.toString(k) : Integer.toString(p.ids[k]);
            sum[k] = 0;
            nodes.put(Integer.parseInt(newick[k]),new Node(Integer.parseInt(newick[k])));
            
            if (Integer.parseInt(newick[k]) > maxId)
                maxId = Integer.parseInt(newick[k]);

            for (int x=0; x<k; x++) {
                sum[k] += p.M[k][x];
                Sum += p.M[k][x];
            }
            for (int x=k+1; x<p.n; x++) {
                sum[k] += p.M[x][k];
                Sum += p.M[x][k];
            }
        }
        
        nextId = maxId+1;

        while (p.n > 2) {
            // Evals S[i][j] and gets minimum:
            double smin = Double.MAX_VALUE;
            int imin=0, jmin=0;

            for (int i=1; i<p.n; i++) {
                for (int j=0; j<i; j++) {

                    double s = (sum[i] - p.M[i][j] + sum[j] - p.M[i][j])/(2*(p.n-2)) +
                                p.M[i][j]/2 + (Sum - sum[i] - sum[j] + p.M[i][j])/(p.n-2);
                    if (s<smin) {
                        smin = s;
                        imin = i;
                        jmin = j;
                    }
                }
            }

            // Stores data on new node (imin,jmin) at jmin:
            double dikmin = (sum[imin] - p.M[imin][jmin])/(p.n-2);
            double djkmin = (sum[jmin] - p.M[imin][jmin])/(p.n-2);

            String oldvalue = newick[jmin];

            newick[jmin] = "(" +
                            newick[imin] + ":" +
                            new Formatter().format(Locale.US,"%1.2f",(p.M[imin][jmin]+dikmin-djkmin)/2-joined[imin]).toString() +
                            "," +
                            newick[jmin] + ":" +
                            new Formatter().format(Locale.US,"%1.2f",(p.M[imin][jmin]+djkmin-dikmin)/2-joined[jmin]).toString() +
                            ")" + Integer.toString(nextId);

            int id = nextId;
            nextId++;
            Node node = new Node(id);
            node.distance.add((p.M[imin][jmin]+dikmin-djkmin)/2-joined[imin]);
            node.distance.add((p.M[imin][jmin]+djkmin-dikmin)/2-joined[jmin]);
            if (newick[imin].lastIndexOf(")") != -1)
                node.children.add(nodes.get(Integer.parseInt(newick[imin].substring(newick[imin].lastIndexOf(")")+1))));
            else
                node.children.add(nodes.get(Integer.parseInt(newick[imin])));
            if (oldvalue.lastIndexOf(")") != -1)
                node.children.add(nodes.get(Integer.parseInt(oldvalue.substring(oldvalue.lastIndexOf(")")+1))));
            else
                node.children.add(nodes.get(Integer.parseInt(oldvalue)));
            nodes.put(id,node);

            joined[jmin] = p.M[imin][jmin]/2;

            sum[jmin] = 0;
            for (int k=0; k<jmin; k++)
                if (k != imin) {
                Sum -= p.M[jmin][k];
                sum[k] -= p.M[jmin][k];
                p.M[jmin][k] = (p.M[jmin][k] + (k<imin ? p.M[imin][k] : p.M[k][imin])) / 2;
                Sum += p.M[jmin][k];
                sum[k] += p.M[jmin][k];
                sum[jmin] += p.M[jmin][k];
                }

            for (int k=jmin+1; k<p.n; k++)
                if (k != imin) {
                Sum -= p.M[k][jmin];
                sum[k] -= p.M[k][jmin];
                p.M[k][jmin] = (p.M[k][jmin] + (k<imin ? p.M[imin][k] : p.M[k][imin])) / 2;
                Sum += p.M[k][jmin];
                sum[k] += p.M[k][jmin];
                sum[jmin] += p.M[k][jmin];
                }

            sum[jmin] += p.M[imin][jmin];

            // Moves n-1 to imin:
            for (int k=0; k<imin; k++) {
                Sum -= p.M[imin][k];
                sum[k] -= p.M[imin][k];
                p.M[imin][k] = p.M[p.n-1][k];
            }

            for (int k=imin+1; k<p.n-1; k++) {
                Sum -= p.M[k][imin];
                sum[k] -= p.M[k][imin];
                p.M[k][imin] = p.M[p.n-1][k];
            }

            Sum -= p.M[p.n-1][imin];

            newick[imin] = newick[p.n-1];
            joined[imin] = joined[p.n-1];
            sum[imin] = sum[p.n-1] - p.M[p.n-1][imin];
            p.n--;
        }

        // 3 points:
        double x = (p.M[1][0]+p.M[2][0]-p.M[2][1])/2 - joined[0];
        double y = (p.M[1][0]+p.M[2][1]-p.M[2][0])/2 - joined[1];
        double z = (p.M[2][0]+p.M[2][1]-p.M[1][0])/2 - joined[2];

        int id = nextId;
        Node node = new Node(id);
        node.distance.add(x);
        node.distance.add(y);
        node.distance.add(z);
        if (newick[0].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[0].substring(newick[0].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[0])));
        if (newick[1].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[1].substring(newick[1].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[1])));
        if (newick[2].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[2].substring(newick[2].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[2])));
        nodes.put(id,node);

        return nodes;

    }

    /**
    * An implementation of Fast NJ.
    *
    * I. Elias, J. Lagergren. Fast Neighbor Joining. Proc. of ICALP 2005.  *
    * @param p A lower triangular distance matrix.  D is destroyed during processing.
    *
    * @return Hashmap containing the nodes of the constructed tree.
    */
    public static HashMap<Integer,Node> fastNJ(PexMatrix p) {

        class closestPair {
            public int j;
            public double d;
        }

        System.out.println("FAST NJ: I. Elias and J. Lagergren. "
            + "Fast Neighbor Joining. In Proceedings of the 32nd International Colloquium on Automata, "
            + "Languages and Programming (ICALP’05), volume 3580, pages 1263–1274, 2005");
        
        int nextId = p.n;
        int maxId = Integer.MIN_VALUE;
        String[] newick = new String[p.n];
        HashMap<Integer,Node> nodes = new HashMap<Integer,Node>();
        for (int i=0; i<p.n; i++) {
            newick[i] = p.ids == null ? Integer.toString(i) : Integer.toString(p.ids[i]);
            nodes.put(Integer.parseInt(newick[i]),new Node(Integer.parseInt(newick[i])));
            if (Integer.parseInt(newick[i]) > maxId)
                maxId = Integer.parseInt(newick[i]);
        }
        nextId = maxId+1;

        // Evals the distances sum for node k:
        double[] sum = new double[p.n];
        for (int k=0; k<p.n; k++) {
            sum[k] = 0;
            for (int x=0; x<k; x++)
                sum[k] += p.M[k][x];
            for (int x=k+1; x<p.n; x++)
                sum[k] += p.M[x][k];
        }

        closestPair[] C = new closestPair[p.n+1];
        int c = p.n;
        // Builds the set of closest pairs:
        for (int i=0; i<p.n; i++) {
            C[i] = new closestPair();

            C[i].d = Double.MAX_VALUE;

            for (int j=0; j<i; j++) {
                double s = p.M[i][j] - (sum[i]+sum[j]) / (p.n-2);
                if (s<C[i].d) {
                C[i].d = s;
                C[i].j = j;
                }
            }

            for (int j=i+1; j<p.n; j++) {
                double s = p.M[j][i] - (sum[i]+sum[j]) / (p.n-2);
                if (s<C[i].d) {
                    C[i].d = s;
                    C[i].j = j;
                }
            }
        }

        while (p.n > 3) {

            // Gets minimum:
            int cmin = 0;
            for (int i=1; i<c; i++)
                if (C[i].d < C[cmin].d)
                cmin = i;

            int imin, jmin;

            if (cmin > C[cmin].j) {
                imin = cmin;
                jmin = C[cmin].j;
            }
            else {
                imin = C[cmin].j;
                jmin = cmin;
            }

            // Data on new node (imin,jmin) will be stored at jmin.
            // Evals branch lengths Lik and Ljk:
            double lik = 0.5 * (p.M[imin][jmin] + ((sum[imin]-sum[jmin])/(p.n-2)));
            double ljk = p.M[imin][jmin] - lik;

            // Updates tree:
            String oldvalue = newick[jmin];

            newick[jmin] = "(" +
                newick[imin] + ":" +
                new Formatter().format(Locale.US,"%.2f",lik).toString() +
                "," +
                newick[jmin] + ":" +
                new Formatter().format(Locale.US,"%.2f",ljk).toString() +
                ")" + Integer.toString(nextId);

            int id = nextId;
            nextId++;
            Node node = new Node(id);
            node.distance.add(Math.abs(lik));
            node.distance.add(Math.abs(ljk));
            if (newick[imin].lastIndexOf(")") != -1)
                node.children.add(nodes.get(Integer.parseInt(newick[imin].substring(newick[imin].lastIndexOf(")")+1))));
            else
                node.children.add(nodes.get(Integer.parseInt(newick[imin])));
            if (oldvalue.lastIndexOf(")") != -1)
                node.children.add(nodes.get(Integer.parseInt(oldvalue.substring(oldvalue.lastIndexOf(")")+1))));
            else
                node.children.add(nodes.get(Integer.parseInt(oldvalue)));
            nodes.put(id,node);

            // Updates D and sum:
            sum[jmin] = 0;
            for (int k=0; k<jmin; k++)
                if (k != imin) {
                    sum[k] -= p.M[jmin][k];
                    p.M[jmin][k] = (p.M[jmin][k] + (k<imin ? p.M[imin][k] : p.M[k][imin]) - p.M[imin][jmin]) / 2;
                    sum[k] += p.M[jmin][k];
                    sum[jmin] += p.M[jmin][k];
                }

            for (int k=jmin+1; k<p.n; k++)
                if (k != imin) {
                    sum[k] -= p.M[k][jmin];
                    p.M[k][jmin] = (p.M[k][jmin] + (k<imin ? p.M[imin][k] : p.M[k][imin]) - p.M[imin][jmin]) / 2;
                    sum[k] += p.M[k][jmin];
                    sum[jmin] += p.M[k][jmin];
                }

            sum[jmin] += p.M[imin][jmin];

            // Moves n-1 onto imin:

            for (int k=0; k<imin; k++) {
                sum[k] -= p.M[imin][k];
                p.M[imin][k] = p.M[p.n-1][k];
            }

            for (int k=imin+1; k<p.n-1; k++) {
                sum[k] -= p.M[k][imin];
                p.M[k][imin] = p.M[p.n-1][k];
            }

            newick[imin] = newick[p.n-1];
            sum[imin] = sum[p.n-1] - p.M[p.n-1][imin];

            p.n--;

            // Updates the closest pairs:

            // Moves n-1 onto imin:
            c--;
            closestPair a = C[imin];
            C[imin] = C[c];
            C[c] = a;

            for (int i=0; i<c; i++) {

                // The closest pairs for the new node, having j=imin or j=jmin must be rebuilt:
                if (i == jmin || C[i].j == imin || C[i].j == jmin) {
                    double s;
                    C[i].d = Double.MAX_VALUE;

                    for (int j=0; j<i; j++) {
                        s = p.M[i][j] - (sum[i]+sum[j]) / (p.n-2);
                        if (s<C[i].d) {
                            C[i].d = s;
                            C[i].j = j;
                        }
                    }

                    for (int j=i+1; j<p.n; j++) {
                        s = p.M[j][i] - (sum[i]+sum[j]) / (p.n-2);
                        if (s<C[i].d) {
                            C[i].d = s;
                            C[i].j = j;
                        }
                    }
                }

                // Any other closest pair must be updated:
                else {
                    if(C[i].j == p.n)
                        C[i].j = imin;
                    }

            }

        }

        // 3 points:
        double x = (p.M[1][0]+p.M[2][0]-p.M[2][1])/2;
        double y = (p.M[1][0]+p.M[2][1]-p.M[2][0])/2;
        double z = (p.M[2][0]+p.M[2][1]-p.M[1][0])/2;

        int id = nextId;
        Node node = new Node(id);
        node.distance.add(x);
        node.distance.add(y);
        node.distance.add(z);
        if (newick[0].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[0].substring(newick[0].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[0])));
        if (newick[1].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[1].substring(newick[1].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[1])));
        if (newick[2].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[2].substring(newick[2].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[2])));
        nodes.put(id,node);

        return nodes;
    }

    /**
    * An implementation of Rapid NJ.
    *
    * M. Simonsen, T. Mailund, C.N.S. Pedersen. Rapid Neighbour-Joining.
    * Proc. of WABI 2008.
    *
    * @param p A lower triangular distance matrix.  D is destroyed during processing.
    * 
    * @return Hashmap containing the nodes of the constructed tree.
    */
    public static HashMap<Integer,Node> rapidNJ(PexMatrix p) {

        System.out.println("RAPID NJ: M. Simonsen, T. Mailund, and C. N. Pedersen. "
            + "Rapid Neighbour-Joining. In Proceedings of WABI 2008, pages 113–122, Karlsruhe, Germany, September 2008");
        
        HashMap<Integer,Node> nodes = new HashMap<Integer,Node>();
        int nextId = p.n;
        int maxId = Integer.MIN_VALUE;

        // Makes room for a larger matrix:

        double[][] E = new double[2*p.n-3][];
        for (int i=0; i<p.n; i++)
            E[i] = p.M[i];
        p.M = E;

        // Labels:
        String[] newick = new String[2*p.n-3];
        for (int i=0; i<p.n; i++) {
            newick[i] = (p.ids == null ? Integer.toString(i) : Integer.toString(p.ids[i]));
            nodes.put(Integer.parseInt(newick[i]),new Node(Integer.parseInt(newick[i])));
            if (Integer.parseInt(newick[i]) > maxId)
                maxId = Integer.parseInt(newick[i]);
        }
        nextId = maxId+1;

        // The sorted matrix S and the indices I:
        double S[][] = new double[2*p.n-3][];
        int I[][] = new int[2*p.n-3][];   // I[i][j] has the position of S[i][j] in D[i].

        for (int i=0; i<p.n; i++) {
            S[i] = new double[i+1];
            I[i] = new int[i+1];

            for (int j=0; j<i; j++)
                S[i][j] = p.M[i][j];

            Sort.qsort(S[i], I[i], 0, i);
        }

        // Evals the distances sum for every node k and records the maximum sum:
        double[] sum = new double[2*p.n-3];
        double smax = Double.MIN_VALUE;

        for (int k=0; k<p.n; k++) {
        sum[k] = 0;
        for (int x=0; x<k; x++)
            sum[k] += p.M[k][x];
        for (int x=k+1; x<p.n; x++)
            sum[k] += p.M[x][k];

        if (sum[k]>smax)
            smax = sum[k];
        }

        // The last row in the matrix is l, the number of active rows is n:
        int l = p.n;

        while (p.n > 3) {

            // Gets minimum:
            double s, smin = Double.MAX_VALUE;
            int imin=0, jmin=0;

            for (int i=l-1; i>=0; i--)
                if (p.M[i] != null)
                    for (int j=0; j<i; j++)
                        if (p.M[I[i][j]] != null)
                            if (S[i][j] < Double.MAX_VALUE && S[i][j] - (sum[i]+smax) / (p.n-2) < smin) {
                                    s = S[i][j] - (sum[i]+sum[I[i][j]]) / (p.n-2);
                                    if (s < smin) {
                                    smin = s;
                                    imin = i;
                                    jmin = I[i][j];
                                }
                            }
                            else
                                break;

            // Data on new node (imin,jmin) will be stored at l.
            // Evals branch lengths Lik and Ljk:
            double lik = 0.5 * (p.M[imin][jmin] + ((sum[imin]-sum[jmin])/(p.n-2)));
            double ljk = p.M[imin][jmin] - lik;

            // Updates tree:
            newick[l] = "(" +
                newick[imin] + ":" +
                new Formatter().format(Locale.US,"%.2f",lik).toString() +
                "," +
                newick[jmin] + ":" +
                new Formatter().format(Locale.US,"%.2f",ljk).toString() +
                ")" + Integer.toString(nextId);//"A" + Integer.toString((l-n)/2);

            int id = nextId;
            nextId++;
            Node node = new Node(id);
            node.distance.add(lik);
            node.distance.add(ljk);
            if (newick[imin].lastIndexOf(")") != -1)
                node.children.add(nodes.get(Integer.parseInt(newick[imin].substring(newick[imin].lastIndexOf(")")+1))));
            else
                node.children.add(nodes.get(Integer.parseInt(newick[imin])));
            if (newick[jmin].lastIndexOf(")") != -1)
                node.children.add(nodes.get(Integer.parseInt(newick[jmin].substring(newick[jmin].lastIndexOf(")")+1))));
            else
                node.children.add(nodes.get(Integer.parseInt(newick[jmin])));
            nodes.put(id,node);

            newick[imin] = newick[jmin] = null;

            // Updates D, S, I, sum:
            sum[l] = 0;
            p.M[l] = new double[l+1];
            S[l] = new double[l+1];
            I[l] = new int[l+1];

            for (int k=0; k<jmin; k++)
                if (k != imin && p.M[k] != null) {
                    sum[k] -= p.M[jmin][k];
                    p.M[l][k] = (p.M[jmin][k] + (k<imin ? p.M[imin][k] : p.M[k][imin]) - p.M[imin][jmin]) / 2;
                    S[l][k] = p.M[l][k];
                    sum[l] += p.M[l][k];
                    sum[k] += p.M[l][k];
                }
                else  // Positions already removed at row l will be "sorted-out"
                    S[l][k] = Double.MAX_VALUE;

            for (int k=jmin+1; k<l; k++)
                if (k != imin && p.M[k] != null) {
                    sum[k] -= p.M[k][jmin];
                    p.M[l][k] = (p.M[k][jmin] + (k<imin ? p.M[imin][k] : p.M[k][imin]) - p.M[imin][jmin]) / 2;
                    S[l][k] = p.M[l][k];
                    sum[l] += p.M[l][k];
                    sum[k] += p.M[l][k];
                }
                else  // Positions already removed at row l will be "sorted-out"
                    S[l][k] = Double.MAX_VALUE;

            for (int k=0; k<imin; k++)
                if (k != jmin && p.M[k] != null)
                    sum[k] -= p.M[imin][k];

            for (int k=imin+1; k<l; k++)
                if (k != jmin && p.M[k] != null)
                    sum[k] -= p.M[k][imin];

            Sort.qsort(S[l],I[l],0,l);

            p.M[imin] = p.M[jmin] = null;
            S[imin] = S[jmin] = null;
            I[imin] = I[jmin] = null;

            // Updates the maximum sum:
            smax = Double.MIN_VALUE;
            for (int k=0; k<l; k++)
                if (p.M[k] != null && sum[k] > smax)
                    smax = sum[k];

            p.n--;
            l++;
        }

        // 3 points:
        // Finds the points out:
        int i,j,k;

        int a=0;
        for ( ; ; a++)
        if (p.M[a] != null) {
            i = a++;
            break;
        }

        for ( ; ; a++)
        if (p.M[a] != null) {
            j = a++;
            break;
        }

        for (; ; a++)
        if (p.M[a] != null) {
            k = a;
            break;
        }

        double x = (p.M[j][i]+p.M[k][i]-p.M[k][j])/2;
        double y = (p.M[j][i]+p.M[k][j]-p.M[k][i])/2;
        double z = (p.M[k][i]+p.M[k][j]-p.M[j][i])/2;

        int id = nextId;
        Node node = new Node(id);
        node.distance.add(x);
        node.distance.add(y);
        node.distance.add(z);
        if (newick[i].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[i].substring(newick[i].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[i])));
        if (newick[j].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[j].substring(newick[j].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[j])));
        if (newick[k].lastIndexOf(")") != -1)
            node.children.add(nodes.get(Integer.parseInt(newick[k].substring(newick[k].lastIndexOf(")")+1))));
        else
            node.children.add(nodes.get(Integer.parseInt(newick[k])));
        nodes.put(id,node);

        return nodes;

    }

    /**
    * Create the object tree, representing the structure created by the methods.
    *
    * @param nodes the list of the nodes, indexed by their identifier
    * @param n number of valid instances
    * @param ids the list of instances ids
    *
    * @return the constructed tree.

    */
    public Tree createTree(HashMap<Integer,Node> nodes, int n, int[] ids) {

        Tree tree = new Tree();
        tree.setType(this.type.toString());

        int maxId = Integer.MIN_VALUE;

        //Inserting valid nodes, from the distance matrix...
        for (int i=0;i<n;i++) {
            ContentTree ct = new ContentTree(ids[i]);
            tree.addNode(ct);
        }

        ids = null;

        ArrayList<Integer> virtualLabels = new ArrayList<Integer>();

        Iterator it = nodes.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry e = (Map.Entry) it.next();
            Integer id = (Integer)e.getKey();
            Node node = (Node)e.getValue();
            if (!node.children.isEmpty()) {//only virtual nodes has children, that is the condition to add them
                virtualLabels.add(id);
                if (id > maxId) maxId = id;
            }
        }
        //Ordering virtual keys so the access is done by creation order.
        Collections.sort(virtualLabels);

        //Inserting virtual nodes
        for (int i=0;i<virtualLabels.size();i++) {
            Node node = nodes.get(virtualLabels.get(i));
            if (node != null) {
                ContentTree ct = new ContentTree(virtualLabels.get(i));
                ct.setValid(false);
                tree.addNode(ct);
            }
        }

        virtualLabels = null;
        System.gc();

        for (int i=0;i<tree.getSize();i++) {
            Node node = nodes.get(tree.getNode(i).getId());
            if (node != null) {
                int maxNivel = -1;
                //Setting children and distances...
                for (int j=0;j<node.children.size();j++) {
                    ContentTree son = tree.getNodeById(node.children.get(j).id);
                    if (son != null) {
                        son.setParent(tree.getNode(i).getId());
                        if (son.getLevel() > maxNivel) maxNivel = son.getLevel();
                        tree.getNode(i).setChildrenId(j,son.getId());
                        tree.getNode(i).setDistChildren(j,node.distance.get(j).floatValue());
                    }
                }
                if (maxNivel != -1)
                    tree.getNode(i).setLevel(maxNivel+1);
            }
        }
        nodes = null;
        if (promotion) promoteLeafs(tree);
        tree.generateEdges();
        return tree;
    }

    /**
    * Executes the entire process to generate the tree from the distance matrix.
    *
    * @param dmat the data that came from the distance matrix file/object   *
    * @return the constructed tree.

    */
    public Tree execute(DistanceMatrix dmat) {

        PexMatrix pexMatrix = DistanceMatrixReader.loadPex(dmat);
        dmat = null;
        if (pexMatrix != null)
            return constructTree(pexMatrix);
        else return null;

    }

    /**
    * Executes the entire process to generate the tree from the distance matrix file.
    *
    * @param dmatFile the file containing the distance matrix
    * @return the constructed tree.

    */
    public Tree execute(String dmatFile) throws IOException, ParseException {

        PexMatrix pexMatrix = DistanceMatrixReader.loadPex(dmatFile);
        dmatFile = null;
        if (pexMatrix != null)
            return constructTree(pexMatrix);
        else return null;

    }

    private Tree constructTree(PexMatrix pexMatrix) {

        Tree tree = null;
        HashMap<Integer, Node> nodes = null;

        long linit, lend, diff, total;
        linit = System.currentTimeMillis();
        int n = pexMatrix.n;
        try {
            if (type.equals(NJAlgorithmType.ORIGINAL)) {
                nodes = originalNJ(pexMatrix);
            } else if (type.equals(NJAlgorithmType.FAST)) {
                nodes = fastNJ(pexMatrix);
            } else if (type.equals(NJAlgorithmType.RAPID)) {
                nodes = rapidNJ(pexMatrix);
            }
            
            if (nodes != null) {
                lend = System.currentTimeMillis();
                diff = lend - linit;
                total = diff;
                System.out.println("Time spent calculating ("+type+") -> " + (diff/1000.0f) + " seconds.");
                linit = System.currentTimeMillis();
                pexMatrix.M = null;
                tree = createTree(nodes,n,pexMatrix.ids);
                pexMatrix = null;
                lend = System.currentTimeMillis();
                diff = lend - linit;
                total += diff;
                System.out.println("- Time spent (Tree Creation) -> " + (diff/1000.0f) + " seconds.");
                System.out.println("TOTAL TIME PACKAGE NJ -> " + (total/1000.0f) + " seconds.");
            } else
                JOptionPane.showMessageDialog(null,"Error collecting tree nodes!");
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            return tree;
        }

    }

    /**
    * Promotion procedure.
    * 
    * This method promotes leafs, according to a pattern ocurrence on the tree. More information, please consult:
    * 
    * PAIVA, J. G., FLORIAN-CRUZ, L., PEDRINI, H., TELLES, G. P., MINGHIM, R.,
    * Improved Similarity Trees and their Application to Visual Data Classification,
    * IEEE Transactions on Visualization and Computer Graphics, 2011
    * 
    * @param t The constructed tree that will be submitted to promotion procedure.

    */
    public static void promoteLeafs(Tree t) {

        System.out.println("Promotion Procedure: PAIVA, J. G., FLORIAN-CRUZ, L., PEDRINI, H., TELLES, G. P., MINGHIM, R.,"
            + " Improved Similarity Trees and their Application to Visual Data Classification,"
            + " IEEE Transactions on Visualization and Computer Graphics, v. 17, p. 2459-2468, 2011");
        
        if (t == null) return;

        long linit, lend, diff;

        ContentTree virtualNode = null; //current virtual node
        ContentTree parentVirtualNode = null; //current virtual node parent
        ContentTree childrenParentVirtualNode = null; //current virtual node son

        linit = System.currentTimeMillis();

        for (int i=0;i<t.getSize();i++) {
            if ((!t.getNode(i).isValid())&&(t.getNode(i).getId() != t.getRootId())) { //only virtual nodes, and no roots (hihgest level)
                virtualNode = t.getNode(i);
                parentVirtualNode = t.getNodeById(virtualNode.getParent());
                if (parentVirtualNode != null) {
                    if (parentVirtualNode.hasChildren()) {
                        for (int j=0;j<parentVirtualNode.getNumChildren();j++) {
                            //analysing the other son...
                            if (parentVirtualNode.getChildrenId(j) != virtualNode.getId()) {
                                childrenParentVirtualNode = t.getNodeById(parentVirtualNode.getChildrenId(j));
                                if (childrenParentVirtualNode != null&&childrenParentVirtualNode.isValid()&&!childrenParentVirtualNode.hasChildren()) {
                                    //Pattern found, performing changing...
                                    if (virtualNode.hasChildren()) {
                                        for (int k=0;k<virtualNode.getNumChildren();k++) {
                                            if (t.getNodeById(virtualNode.getChildrenId(k)) != null) {
                                                t.getNodeById(virtualNode.getChildrenId(k)).setParent(childrenParentVirtualNode.getId());
                                                float d = virtualNode.getDistChildren(k) +
                                                        (parentVirtualNode.getDistChildren(parentVirtualNode.getChildrenIndex(virtualNode.getId()))/2) +
                                                        (parentVirtualNode.getDistChildren(parentVirtualNode.getChildrenIndex(childrenParentVirtualNode.getId()))/4);
                                                if (childrenParentVirtualNode.getChildrenId(0) == -1) {
                                                    childrenParentVirtualNode.setChildrenId(0,virtualNode.getChildrenId(k));
                                                    childrenParentVirtualNode.setDistChildren(0,d);
                                                } else {
                                                    childrenParentVirtualNode.setChildrenId(1,virtualNode.getChildrenId(k));
                                                    childrenParentVirtualNode.setDistChildren(1,d);
                                                }
                                            }
                                        }
                                    }
                                    //Setting promoted node father...
                                    childrenParentVirtualNode.setParent(parentVirtualNode.getParent());
                                    childrenParentVirtualNode.setLevel(parentVirtualNode.getLevel());
                                    if ((parentVirtualNode.getParent() != -1)&&(t.getNodeById(parentVirtualNode.getParent())!= null)&&
                                        (t.getNodeById(parentVirtualNode.getParent()).hasChildren())) {
                                        //In this case, a node does not appear as son of the other, and vice-versa, and I will not remove
                                        //it from the list of the other, nor add the other on the list of this node
                                        if ((t.getNodeById(parentVirtualNode.getParent()) != null)&&
                                            (t.getNodeById(t.getNodeById(parentVirtualNode.getParent()).getParent()) != null)&&
                                            (t.getNodeById(parentVirtualNode.getParent()).getParent() == parentVirtualNode.getId())) {
                                            t.getNodeById(childrenParentVirtualNode.getParent()).setParent(childrenParentVirtualNode.getId());
                                            //These are the last nodes (roots), setting their levels to equal values
                                            if (t.getNodeById(parentVirtualNode.getParent()).getLevel() > childrenParentVirtualNode.getLevel()) {
                                                childrenParentVirtualNode.setLevel(t.getNodeById(parentVirtualNode.getParent()).getLevel());
                                            } else {
                                                if (t.getNodeById(childrenParentVirtualNode.getParent()) != null)
                                                    t.getNodeById(childrenParentVirtualNode.getParent()).setLevel(childrenParentVirtualNode.getLevel());
                                            }
                                        }else {
                                            //collecting the distance between the promoted node and its parent...
                                            float dp = 0;
                                            if (t.getNodeById(parentVirtualNode.getParent()) != null) {
                                                int ind = t.getNodeById(parentVirtualNode.getParent()).getChildrenIndex(parentVirtualNode.getId());
                                                float k = t.getNodeById(parentVirtualNode.getParent()).getDistChildren(ind);
                                                float b = parentVirtualNode.getDistChildren(parentVirtualNode.getChildrenIndex(virtualNode.getId()));
                                                float a = parentVirtualNode.getDistChildren(parentVirtualNode.getChildrenIndex(childrenParentVirtualNode.getId()));
                                                dp = k + b/2 + a/2;
                                            }
                                            if (t.getNodeById(parentVirtualNode.getParent()) != null) {
                                                int index = t.getNodeById(parentVirtualNode.getParent()).getChildrenIndex(parentVirtualNode.getId());
                                                if (index != -1) {
                                                    t.getNodeById(parentVirtualNode.getParent()).setChildrenId(index,childrenParentVirtualNode.getId());
                                                    t.getNodeById(parentVirtualNode.getParent()).setDistChildren(index,dp);
                                                }
                                            }
                                        }
                                    }
                                    parentVirtualNode.setParent(-1);
                                    virtualNode.setParent(-1);
                                    //If the parentVirtualNode has 3 children, he is the tree root. I will exclude it, so 2 nodes will last,
                                    //and they will be joined, forming the new tree.
                                    if (parentVirtualNode.getNumChildren() == 3) {
                                        for (int x=0;x<parentVirtualNode.getNumChildren();x++) {
                                            if (parentVirtualNode.getChildrenId(x) != childrenParentVirtualNode.getId()&&
                                                parentVirtualNode.getChildrenId(x) != virtualNode.getId()) {
                                                childrenParentVirtualNode.setParent(parentVirtualNode.getChildrenId(x));
                                                if (t.getNodeById(parentVirtualNode.getChildrenId(x)) != null) {
                                                    t.getNodeById(parentVirtualNode.getChildrenId(x)).setParent(childrenParentVirtualNode.getId());
                                                    //In this case, the parentvirtualnode level was the highest one, so its children
                                                    //will be set with this level, because they will be the new roots.
                                                    t.getNodeById(parentVirtualNode.getChildrenId(x)).setLevel(parentVirtualNode.getLevel());
                                                }
                                                childrenParentVirtualNode.setLevel(parentVirtualNode.getLevel());
                                            }
                                        }
                                    }
                                    t.removeNode(parentVirtualNode);
                                    t.removeNode(virtualNode);
                                    virtualNode = null;
                                    i--;
                                    //This break is here because the root node is the only one that has 3 children, and the process
                                    //can reach this part before the last one. So it is for forcing the 'for' exit, because the change
                                    //was already done, and it is not necessary to iterate over the children.
                                    //sobre os filhos.
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
        lend = System.currentTimeMillis();
        diff = lend - linit;
        System.out.println("Time spent (Leaf promotion) -> " + (diff/1000.0f) + " seconds");
    }

}
package tree.technique.nj;

import java.util.Random;

/**
 * Sorting algorithms.
 * @author Guilherme P. Telles, Mar 2011. 
 */
public class Sort {

  /**
    A quicksort for v[l,r).

   * @param v
   * @param l
   * @param r
   */
  public static void qsort(int[] v, int l, int r) {

    Random rand = new Random();
    int p,i,j,tmp;
    r--;

    // A stack:
    int[] s = new int[20*(int) (Math.log(r-l)/Math.log(2)+1)];
    int ss = 0;

    s[ss++] = l;
    s[ss++] = r;

    while (ss>0) {
      // Pops a subproblem:
      r = j = s[--ss];
      l = i = s[--ss];
      p = l+rand.nextInt(r-l);
      tmp = v[i];  v[i] = v[p];  v[p] = tmp;
      p = i;
      
      // Partition:
      while (i < j) {
        while (i <= r && v[i] <= v[p])
          i++;

        while (v[j] > v[p])
          j--;

        if (i < j) {
          tmp = v[i];  v[i] = v[j];  v[j] = tmp;
        }
      }

      // j is the pivot position:
      tmp = v[p];  v[p] = v[j];  v[j] = tmp;

      // Pushes subproblems:
      if (j+1 < r) {
        s[ss++] = j+1;
        s[ss++] = r;
      }
      if (l < j-1) {
        s[ss++] = l;
        s[ss++] = j-1;
      }
    }
  }

  /**
    A quicksort for A[l,r) that also produces I[l,r) in-place, that
    indexes elements in A to their positions before sorting. 

   * @param A The array to sort.  Any position out of I[l,r) is left unchanged. 
   * @param I I[i] will hold the position of A[i] in A before sorting take place. Any 
   * position out of I[l,r) is left unchanged. 
   * @param l Left limit for the sorting interval (closed).
   * @param r Right limit for the sorting interval (open).
   */
  public static void qsort(double[] A, int[] I, int l, int r) {

    Random rand = new Random();
    int p,i,j,tmpi;
    double tmpd;
    r--;

    // Sets I to the identity:
    for (i=l; i<=r; i++)
      I[i] = i;

    if (r-l < 1)
      return;
      
    // A stack:
    int[] s = new int[20*(int) (Math.log(r-l)/Math.log(2)+1)];
    
    int ss = 0;

    s[ss++] = l;
    s[ss++] = r;

    while (ss>0) {
      // Pops a subproblem:
      r = j = s[--ss];
      l = i = s[--ss];
      p = l+rand.nextInt(r-l);
      
      tmpd = A[i];  A[i] = A[p];  A[p] = tmpd;
      tmpi = I[i];  I[i] = I[p];  I[p] = tmpi;
      p = i;
      
      // Partition:
      while (i < j) {
        while (i <= r && A[i] <= A[p])
          i++;

        while (A[j] > A[p])
          j--;

        if (i < j) {
          tmpd = A[i];  A[i] = A[j];  A[j] = tmpd;
          tmpi = I[i];  I[i] = I[j];  I[j] = tmpi;
        }
      }

      // j is the pivot position:
      tmpd = A[p];  A[p] = A[j];  A[j] = tmpd;
      tmpi = I[p];  I[p] = I[j];  I[j] = tmpi;

      // Pushes subproblems:
      if (j+1 < r) {
        s[ss++] = j+1;
        s[ss++] = r;
      }
      if (l < j-1) {
        s[ss++] = l;
        s[ss++] = j-1;
      }
    }
  }

  /**
   * An insertion sort for A[l,r).
   * 
   * @param A
   * @param l 
   * @param r
   */
  public static void isort(int[] A, int l, int r) {

    for (int j=l+1; j<r; j++) {
      int k = A[j];
      int i = j-1;
      while (i>=l && A[i]>k) {
        A[i+1] = A[i];
        i--;
      }
      A[i+1] = k;
    }
  }

  /**
   * A Shellsort for A[l,r).
   * 
   * @param A
   * @param l
   * @param r
   */
  public static void ssort(int[] A, int l, int r) {
  
    // Ciura increments:
    int[] incs = {960333, 436515, 198416, 90189, 40995, 18634, 8470, 3850,
        1750, 701, 301, 132, 57, 23, 10, 4, 1};
  
    // For each increment incs[t], element j is moved across positions pointed by i.
    for (int t=0; t<17; t++) {
      for (int j=l+incs[t]; j<r; j++) {
        int k = A[j];
        int i = j;
        while (i >= incs[t] && A[i-incs[t]] > k) {
          A[i] = A[i-incs[t]];
          i -= incs[t];
        }
        A[i] = k;
      }
    }
  }

  public static void main(String[] args) {

    int z[][] = {{49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0},
        {49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};
    
    
    double v[] = {6,9,6,5,0,8}; 
    
    int I[] = new int[50];

    int n = 6;

    int i;
    for (i=0; i<n; i++) {
      System.out.printf("%4d ", i);
    }
    System.out.printf("\n");

    for (i=0; i<n; i++) {
      System.out.printf("%4.1f ",v[i]);
    }
    System.out.printf("\n");

    Sort.qsort(v,I,0,n);

    for (i=0; i<n; i++) {
      System.out.printf("%4.1f ",v[i]);
    }
    System.out.printf("\n");

    for (i=0; i<n; i++) {
      System.out.printf("%4d ",I[i]);
    }
    System.out.printf("\n");
  }
 
}

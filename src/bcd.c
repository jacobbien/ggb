#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>

// given the lower triangle of a p by p matrix, this gets the jk-th element
// note: must have j > k
#define lt(j, k, p) p * k - k * (k + 1) / 2 + j - k - 1

void allocate_groups(int *dd, int p, int *M, int ****elem, int ***nelem, int ***ngelem);
void free_groups(int p, int *M, int ****elem, int ***nelem, int ***ngelem);
void allocate_v(int *M, int ***elem, int **nelem, int p, double ****v);
void free_v(int *M, int ***elem, int p, double ****v);
void compute_lamlist(double *rr, int p, int *M, int ***elem, int **nelem, int **ngelem, int nlam, double flmin, double *lamlist);
void compute_sumvv(double ***v, int p, int *M, int ***elem, int **nelem, int **ngelem, double *sumvv);
double objective(double ***v, double lam, int p, int *M, int ***elem, int **nelem, int **ngelem, double *rr0);
void print_matrix(double *mat, int nr, int nc);
void print_lowertri(double *mm, int nr, int nc);
void prox_bcd3(double *rr, int *M, double lam, int p, int ***elem, int **nelem,
               int **ngelem, double ***v, double *rr0, int maxiter,
               double tol, int verbose);

void pathwiseprox_bcd3(double *rr, int *dd, int *M, double *lamlist, int *nlam, double *flmin,
                       int *p, double *sumvv, double *obj, int *maxiter, double *tol, int *verbose) {
  /*
  Evaluate the ggb proximal operator with no psd constraint along a grid of lambda values.

   Returns sumvv which can be thought of as a cp2 by nlam matrix, each column being the
   lower triangle of sum_jb V_jb.

  argmin 0.5 * || R - sum_jb V_jb ||_F^2 + lam * sum_jb w_jb || V_jb ||_F
  where V_jb is 0 off of g_jb.

  rr = lower-tri(R)
  dd[lt(j, k, *p)] = num hops between j and k
  M[j] = maximal bandwidth for variable j

  g_jb = union(s_jb': 1 <= b' <= b) for b = 1, ..., M[j].

  The ith element of rr in s_jb will be rr[elem[j][b][i]] and
  the size of s_jb will be 2*nelem[j][b].

  Args:
  rr: length p * (p-1) / 2 array such that rr = lower-tri(R)
  lamlist: array of non-negative doubles
  p: dimension of R

   Note on w: for now, default weights will be calculated within: w[j][b] = |g_jb|^1/2
  */
  int i, ii, j, b;
  int cp2 = *p * (*p - 1) / 2;
  int l; // lambda index
  int ***elem; // elem[j][b][i] is ith element of s_jb
  int **nelem; // nelem[j][b] is |s_jb| / 2
  int **ngelem; // nelem[j][b] is |g_jb| / 2
  double ***v; // v[j][b][i] is ith element of (V_jb)_(g_jb)
  allocate_groups(dd, *p, M, &elem, &nelem, &ngelem); // also initializes
  allocate_v(M, elem, ngelem, *p, &v);

  // initialize each v[j][b] as all 0s
  for (j = 0; j < *p; j++) {
    if (M[j] == 0) continue;
    for (b = 0; b < M[j]; b++) for (ii = 0; ii < ngelem[j][b]; ii++) v[j][b][ii] = 0;
  }
  double *rr0 = malloc(cp2 * sizeof(double)); // an original copy of rr
  for (i = 0; i < cp2; i++) rr0[i] = rr[i];

  // get lamlist
  if (*flmin != -1) {
    // this means that lamlist was not passed and we need to compute it from flmin and nlam
    if (*verbose > 0)
      Rprintf("computing lamlist... with nlam = %d and flmin = %g\n", *nlam, *flmin);
    compute_lamlist(rr, *p, M, elem, nelem, ngelem, *nlam, *flmin, lamlist);
    if (*verbose > 0) print_matrix(lamlist, *nlam, 1);
  }

  for (l = 0; l < *nlam; l++) { // loop over lambda values
    if (*verbose > 0) Rprintf("lam = %g:\n", lamlist[l]);
    // update rr and v and sumvv for this lambda:
    prox_bcd3(rr, M, lamlist[l], *p, elem, nelem, ngelem, v, rr0, *maxiter, *tol, *verbose);
    // Form lower triangle of a p by p matrix sum_jb V_jb:
    compute_sumvv(v, *p, M, elem, nelem, ngelem, sumvv + (cp2 * l));
    obj[l] = objective(v, lamlist[l], *p, M, elem, nelem, ngelem, rr0);
    // note that sumvv + (cp2 * l) points to sumvv[cp2 *l], which is first element of the lth
    // array of length l.
  }
  free_v(M, elem, *p, &v);
  free_groups(*p, M, &elem, &nelem, &ngelem);
  free(rr0);
}


void prox_bcd3(double *rr, int *M, double lam, int p, int ***elem, int **nelem,
               int **ngelem, double ***v, double *rr0, int maxiter,
               double tol, int verbose) {
  /*
   Evaluate the ggb proximal operator with no psd constraint.

   argmin 0.5 * || R - sum_jb V_jb ||_F^2 + lam * sum_jb w_jb || V_jb ||_F
  where V_jb is 0 off of g_jb.

  rr = lower-tri(R)
  dd[lt(j, k, *p)] = num hops between j and k
  M[j] = maximal bandwidth for variable j

  g_jb = union(s_jb': 1 <= b' <= b) for b = 1, ..., M[j].

  The ith element of rr in s_jb will be rr[elem[j][b][i]] and
  the size of s_jb will be 2*nelem[j][b].

  Args:
  rr: length p * (p-1) / 2 array such that rr = lower-tri(R)
  lamlist: array of non-negative doubles
  p: dimension of R

  Note on w: for now, default weights will be calculated within: w[j][b] = |g_jb|^1/2
  */

  int i, ii, j, b, bb;
  int l; // cycle index
  double r, temp, del, maxdel;
  for (l = 0; l < maxiter; l++) {
    if (verbose > 1) {
      Rprintf(" objective is %.7g\n",
             objective(v, lam, p, M, elem, nelem, ngelem, rr0));
    }
    maxdel = 0; // largest change in any V[i] in this cycle
    // cycle over all groups
    for (j = 0; j < p; j++) {
      if (M[j] == 0) continue;
      for (b = 0; b < M[j]; b++) {
        // group jb
        ii = 0;
        for (bb = 0; bb <= b; bb++) {
          for (i = 0; i < nelem[j][bb]; i++) {
            rr[elem[j][bb][i]] += v[j][b][ii++];
          }
        }
        // now rr is jb-th partial residual R_jb
        // compute ||R_jb||:
        r = 0;
        for (bb = 0; bb <= b; bb++) {
          for (i = 0; i < nelem[j][bb]; i++) {
            r += rr[elem[j][bb][i]] * rr[elem[j][bb][i]];
          }
        }
        // now r is ||R_jb||^2 / 2
        r = sqrt(r / ngelem[j][b]) / lam; // using w_jb = 2*ngelem[j][b] here
        // now r = ||R_jb|| / (lam * w_jb)
        //printf("||R_jb|| / (lam * w_jb) = %g\n", r);
        if (r <= 1) {
          // zero out V_jb
          for (ii = 0; ii < ngelem[j][b]; ii++) {
            del = fabs(v[j][b][ii]);
            if (del > maxdel) maxdel = del;
            v[j][b][ii] = 0;
          }
        } else {
          // shrink V_jb
          r = 1 - 1 / r; // now r is 1 - lam * w_jb / ||R_jb||
          ii = 0;
          for (bb = 0; bb <= b; bb++) {
            for (i = 0; i < nelem[j][bb]; i++) {
              temp = r * rr[elem[j][bb][i]];
              del = fabs(v[j][b][ii] - temp);
              if (del > maxdel) maxdel = del;
              v[j][b][ii] = temp;
              ii++;
            }
          }
        }
        ii = 0;
        for (bb = 0; bb <= b; bb++) {
          for (i = 0; i < nelem[j][bb]; i++) {
            rr[elem[j][bb][i]] -= v[j][b][ii++];
          }
        }
        // now R is R (rather than R_jb)
      }
    }
    if (maxdel < tol) {
      if (verbose > 0) Rprintf("Converged after %d iterations.\n", l);
      break;
    }
  }
}

void allocate_v(int *M, int ***elem, int **nelem, int p, double ****v) {
  /*
   * Allocates (but does not initialize) three dimensional double array *v.
   * (*v)[j][b][i] is ith element of V_jb.
   */
  int j, b;
  *v = malloc(p * sizeof(double **));
  for (j = 0; j < p; j++) {
    if (M[j] == 0) continue;
    (*v)[j] = malloc(M[j] * sizeof(double *));
    for (b = 0; b < M[j]; b++) {
      (*v)[j][b] = malloc(nelem[j][b] * sizeof(double));
    }
  }
}

void free_v(int *M, int ***elem, int p, double ****v) {
  /*
  * Frees memory malloc'ed for three dimensional double array *v.
  */
  int j, b;
  for (j = 0; j < p; j++) {
    if (M[j] == 0) continue;
    for (b = 0; b < M[j]; b++) {
      free((*v)[j][b]); // since we did malloc(nelem[j][b] * sizeof(double));
    }
    free((*v)[j]); // since we did malloc(M[j] * sizeof(double *));
  }
  free(*v); // since we did malloc(p * sizeof(double **));
}


void allocate_groups(int *dd, int p, int *M, int ****elem, int ***nelem, int ***ngelem) {
  /*
   * Given the lower triangle of the p-by-p matrix of graph distances
   * and a p-vector M of maximal bandwidths, this function returns the
   * group structure g_jb in the form of elem, ends, ngroups.  In particular,
   * if rr is the lower triangle of a symmetric matrix R, then the elements
   * of rr that are in group s_jb are given by *elem[j][b-1]
   *
   * Explanation of int ****elem:
   * int *elem would point to an int
   * int **elem would point to an array of ints
   * int ***elem would point to an array of arrays of ints
   * but we want int ****elem since it points to an array of arrays of arrays of ints
   * in particular, (*elem)[j][b-1] is the array of indices in s_jb.
   * in words, elem points to a three-dimensional array of ints.  (*elem)[j][b-1][i]
   * is element i of group s_jb.
   *
   * recall that s_jb are {jk: j > k, d_G(j,k) == b}
   *
   * *nelem[j][b-1] gives |s_jb|/2, the length of *elem[j][b-1].
   * *ngelem[j][b-1] gives |g_jb|/2.
   *
   * The user inputs M: M[j] gives the number of arrays of the form *elem[j][b-1]
   * (1 <= b <= M[j]).
   *
   * Whenever we cycle over all groups, we first must check that M[j] > 0.
   * Likewise, when we cycle over elements within s_jb, we first check that
   * nelem[j][b-1] > 0.
   */
  int j, k, b;
  // allocate *nelem, an array of p int arrays

  *nelem = malloc(p * sizeof(int *));
  *ngelem = malloc(p * sizeof(int *));
  int **counter = malloc(p * sizeof(int *));

  for (j = 0; j < p; j++) {
    if (M[j] > 0) {
      // allocate *nelem[j] an array of M[j] ints.
      (*nelem)[j] = malloc(M[j] * sizeof(int));
      (*ngelem)[j] = malloc(M[j] * sizeof(int));
      counter[j] = malloc(M[j] * sizeof(int));
      // initialize as zero
      for (b = 0; b < M[j]; b++) {
        (*nelem)[j][b] = 0;
        (*ngelem)[j][b] = 0;
        counter[j][b] = 0;
      }
    }
  }
  // scan through dd, tallying num jk in each s_jb so *nelem[j][b] = |s_jb| / 2
  for (k = 0; k < p - 1; k++) {
    for (j = k + 1; j < p; j++) {
      b = dd[lt(j, k, p)];
      if (b <= M[j]) // if b > M[j], we don't add this to any group
        ((*nelem)[j][b - 1])++; // dd[lt(j, k, p)] is num hops b/w j and k
      if (b <= M[k]) // and same logic for k
        ((*nelem)[k][b - 1])++;
    }
  }
  // now initialize *ngelem = |g_jb| / 2:
  for (j = 0; j < p; j++) {
    if (M[j] == 0) continue;
    // |g_j0| = |s_j0|
    (*ngelem)[j][0] = (*nelem)[j][0];
    for (b = 1; b < M[j]; b++) {
      // |g_jb| = |g_j,b-1| + |s_jb|
      (*ngelem)[j][b] = (*ngelem)[j][b - 1] + (*nelem)[j][b];
    }
  }
  // allocate *elem, an array of p two-dimensional int arrays
  *elem = malloc(p * sizeof(int **));
  for (j = 0; j < p; j++) {
    if (M[j] > 0) {
      // allocate *elem[j], which is an array of M[j] one-dimensional int arrays.
      (*elem)[j] = malloc(M[j] * sizeof(int *));
      for (b = 0; b < M[j]; b++) {
        (*elem)[j][b] = malloc((*nelem)[j][b] * sizeof(int));
      }
    }
  }
  // scan through dd, giving indices contained in each s_jb.
  // counter[j][b] will keep track of how many of the elements of s_jb we have
  // filled in
  for (k = 0; k < p - 1; k++) {
    for (j = k + 1; j < p; j++) {
      // put index lt(j, k, p) in group g_jb where b is dd[lt(j, k, p)], the
      // num hops b/w j and k
      b = dd[lt(j, k, p)];
      if (b <= M[j]) // if b > M[j], we don't add this to any group
        (*elem)[j][b - 1][counter[j][b - 1]++] = lt(j, k, p);
      if (b <= M[k]) // same logic for k
        (*elem)[k][b - 1][counter[k][b - 1]++] = lt(j, k, p);
    }
  }
  for (j = 0; j < p; j++) {
    if (M[j] > 0) free(counter[j]); // since we had malloc(M[j] * sizeof(int));
  }
  free(counter); // since we had malloc(p * sizeof(int *));
}

void free_groups(int p, int *M, int ****elem, int ***nelem, int ***ngelem) {
  /*
   * This should be called after allocate_groups is called.
   */
  int j, b;

  for (j = 0; j < p; j++) {
    if (M[j] > 0) {
      for (b = 0; b < M[j]; b++) {
        free((*elem)[j][b]); // since we had called malloc((*nelem)[j][b] * sizeof(int));
      }
      free((*elem)[j]); // since we had called malloc(M[j] * sizeof(int *));
    }
  }
  free(*elem); // since we had called malloc(p * sizeof(int **));

  for (j = 0; j < p; j++) {
    if (M[j] == 0) continue;
    free((*nelem)[j]); // since before malloc(M[j] * sizeof(int));
    free((*ngelem)[j]); // since before malloc(M[j] * sizeof(int));
  }
  free(*nelem); // since before malloc(p * sizeof(int *));
  free(*ngelem); // since before malloc(p * sizeof(int *));
}

void compute_lamlist(double *rr, int p, int *M, int ***elem, int **nelem, int **ngelem,
                     int nlam, double flmin, double *lamlist) {
  int i, j, b;
  double tt;
  lamlist[0] = 0;
  // start by computing lamlist[0] = lammax:
  // lammax^2 = max_jb ||R_{jb}||_F^2 / |g_jb|
  for (j = 0; j < p; j++) {
    tt = 0; // tt will be ||R_g_jb||_F^2
    if (M[j] == 0) continue;
    for (b = 0; b < M[j]; b++) {
      // group jb
      // we use that ||R_g_jb||_F^2 = ||R_g_jb-1||_F^2 + ||R_s_jb||^2
      for (i = 0; i < nelem[j][b]; i++) {
          tt += pow(rr[elem[j][b][i]], 2);
      }
      if (tt / ngelem[j][b] > lamlist[0]) lamlist[0] = tt / ngelem[j][b];
    }
  }
  lamlist[0] = sqrt(lamlist[0]);
  tt = pow(flmin, 1.0 / (nlam - 1));
  for (i = 1; i < nlam; i++) {
    lamlist[i] = lamlist[i - 1] * tt;
  }
}

void compute_sumvv(double ***v, int p, int *M, int ***elem, int **nelem, int **ngelem, double *sumvv) {
  /*
   * This function takes v and produces the lower triangle of a p by p matrix sum(V_jb).  Assumes memory
   * for sumvv has already been allocated (as a cp2 length double array)
   */
  int i, ii, j, b, bb;
  int cp2 = p * (p - 1) / 2;
  for (i = 0; i < cp2; i++) sumvv[i] = 0;
  for (j = 0; j < p; j++) {
    if (M[j] == 0) continue;
    for (b = 0; b < M[j]; b++) {
      // group jb
      ii = 0;
      for (bb = 0; bb <= b; bb++) {
        for (i = 0; i < nelem[j][bb]; i++) {
          sumvv[elem[j][bb][i]] += v[j][b][ii++];
        }
      }
    }
  }
}

double objective(double ***v, double lam, int p, int *M, int ***elem, int **nelem, int **ngelem, double *rr0) {
  int i, ii, j, b;
  int cp2 = p * (p - 1) / 2;
  double obj = 0, tt;
  double *sumvv = malloc(cp2 * sizeof(double));
  compute_sumvv(v, p, M, elem, nelem, ngelem, sumvv);
  for (i = 0; i < cp2; i++) {
    tt = rr0[i] - sumvv[i];
    obj += tt * tt;
  }
  // obj is now 0.5 * || R - Sig ||^2_F
  for (j = 0; j < p; j++) {
    if (M[j] == 0) continue;
    for (b = 0; b < M[j]; b++) {
      tt = 0;
      for (ii = 0; ii < ngelem[j][b]; ii++) tt += v[j][b][ii] * v[j][b][ii];
      obj += lam * 2 * sqrt(ngelem[j][b] * tt); // using w_jb = |g_jb|^1/2
      // factor of 2 is b/c w_jb = sqrt(2 * ngelem[j][b]) and ||V_jb|| = sqrt(2*tt)
    }
  }
  free(sumvv);
  return obj;
}

void print_matrix(double *mat, int nr, int nc) {
  int i;
  int j;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      Rprintf("%.2g ", mat[i + nr * j]);
    }
    Rprintf("\n");
  }
}

void print_lowertri(double *mm, int nr, int nc) {
  int i, j, ii = 0;
  for (i = 0; i < nr; i++) {
    for (j = 0; j < nc; j++) {
      if (i > j)
        Rprintf("%.2g\t", mm[ii++]);
      else
        Rprintf("*\t");
    }
    Rprintf("\n");
  }
}


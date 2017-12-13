'use strict';

// n, /* A and L are n-by-n, where n >= 0 */
// Ap, /* input of size n + 1, not modified */
// Ai, /* input of size nz=Ap[n], not modified */
// Lp, /* output of size n + 1, not defined on input */
// Parent, /* output of size n, not defined on input */
// Lnz, /* output of size n, not defined on input */
// Flag /* workspace of size n, not defn. on input or output */

function ldlSymbolic(n, Ap, Ai, Lp, Parent, Lnz, Flag) {
    for (let i = 0; i < n; i++) {
        /* L(k,:) pattern: all nodes reachable in etree from nz in A(0:k-1,k) */
        Parent[i] = -1; /* parent of k is not yet known */
        Flag[i] = i; /* mark node k as visited */
        Lnz[i] = 0; /* count of nonzeros in column k of L */
        let p2 = Ap[i + 1];
        for (let j = Ap[i]; j < p2; j++) {
        /* A (row,k) is nonzero (original or permuted A) */
            let row = Ai[j];
            if (row < i) {
                /* follow path from row to root of etree, stop at flagged node */
                for (; Flag[row] !== i; row = Parent[row]) {
                    /* find parent of row if not yet determined */
                    if (Parent[row] === -1) {
                        Parent[row] = i;
                    }
                    Lnz[row]++; /* L (k,row) is nonzero */
                    Flag[row] = i; /* mark row as visited */
                }
            }
        }
    }
    /* construct Lp index array from Lnz column counts */
    Lp[0] = 0;
    for (let i = 0; i < n; i++) {
        Lp[i + 1] = Lp[i] + Lnz[i];
    }
}


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

// n, /* A and L are n-by-n, where n >= 0 */
// Ap, /* input of size n+1, not modified */
// Ai, /* input of size nz=Ap[n], not modified */
// Ax, /* input of size nz=Ap[n], not modified */
// Lp, /* input of size n+1, not modified */
// Parent, /* input of size n, not modified */
// Lnz, /* output of size n, not defn. on input */
// Li, /* output of size lnz=Lp[n], not defined on input */
// Lx, /* output of size lnz=Lp[n], not defined on input */
// D, /* output of size n, not defined on input */
// Y, /* workspace of size n, not defn. on input or output */
// Pattern, /* workspace of size n, not defn. on input or output */
// Flag /* workspace of size n, not defn. on input or output */

function ldl_numeric(n,Ap,Ai,Ax,Lp,Parent,Lnz,Li,Lx,D,Y,Pattern,Flag) {
    var yi, l_ki;
    var i, k, p, kk, p2, len, top;
    for (k = 0; k < n; k++) {
    /* compute nonzero Pattern of kth row of L, in topological order */
        Y[k] = 0.0; /* Y(0:k) is now all zero */
        top = n; /* stack for pattern is empty */
        Flag[k] = k; /* mark node k as visited */
        Lnz[k] = 0; /* count of nonzeros in column k of L */
        kk = (k); /* kth original, or permuted, column */
        p2 = Ap[kk + 1];
        for (p = Ap[kk]; p < p2; p++) {
            i = (Ai[p]); /* get A(i,k) */
            if (i <= k) {
                Y[i] += Ax[p]; /* scatter A(i,k) into Y (sum duplicates) */
                for (len = 0; Flag[i] !== k; i = Parent[i]) {
                    Pattern[len++] = i; /* L(k,i) is nonzero */
                    Flag[i] = k; /* mark i as visited */
                }
                while (len > 0) Pattern[--top] = Pattern[--len];
            }
        }
        /* compute numerical values kth row of L (a sparse triangular solve) */
        D[k] = Y[k]; /* get D(k,k) and clear Y(k) */
        Y[k] = 0.0;
        for (; top < n; top++) {
            i = Pattern[top]; /* Pattern[top:n-1] is pattern of L(:,k) */
            yi = Y[i]; /* get and clear Y(i) */
            Y[i] = 0.0;
            p2 = Lp[i] + Lnz[i];
            for (p = Lp[i]; p < p2; p++) {
                Y[Li[p]] -= Lx[p] * yi;
            }
            l_ki = yi / D[i]; /* the nonzero entry L(k,i) */
            D[k] -= l_ki * yi;
            Li[p] = k; /* store L(k,i) in column form of L */
            Lx[p] = l_ki;
            Lnz[i]++; /* increment count of nonzeros in col i */
        }

        if (D[k] === 0.0) return (k); /* failure, D(k,k) is zero */
    }

    return (n); /* success, diagonal of D is all nonzero */
}


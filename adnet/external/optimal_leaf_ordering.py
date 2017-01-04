# copied from: https://github.com/jamestwebber/jtw-utils/blob/master/optimal_leaf_ordering.py

# The MIT License (MIT)
#
# Copyright (c) 2015 James Webber
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.


# code initially adapted @markak on GitHub, algorithm is from
# Ziv Bar-Joseph et al., Bioinformatics 2001

import itertools

import numpy as np

from scipy.cluster import hierarchy
from scipy.spatial import distance


def order_tree(Z, rd, M):
    def swap_subtrees(n):
        n.right, n.left = n.left, n.right

    for v in xrange(Z.shape[0] * 2, Z.shape[0], -1):
        L,R = rd[v].left.pre_order(), rd[v].right.pre_order()
        u,w = min(itertools.product(L, R),
                  key=lambda (u,w): M[v, u, w])

        if rd[v].left.count > 1:
            LR = rd[v].left.right.pre_order()
            if u in LR:
                swap_subtrees(rd[v].left)

        if rd[v].right.count > 1:
            RL = rd[v].right.left.pre_order()
            if w in RL:
                swap_subtrees(rd[v].right)


def optimal_scores(Z, rd, dists):
    # Z - linkage matrix from scipy.cluster.hierarchy
    # rd - ClusterNode dictionary from to_tree
    # dists - distance matrix

    n_nodes = Z.shape[0] + 1

    M = {}

    # iterating through the linkage matrix guarantees
    # we never see a node before its children
    for i in xrange(Z.shape[0]):
        # linkage matrix starts at first non-leaf node
        v = n_nodes + i
        # the left and right nodes
        j,k = int(Z[i, 0]), int(Z[i, 1])

        if Z[i, 3] == 2:
            # both j and k are leaves, so there is no ordering to be done
            M[v, j, k] = M[v, k, j] = Z[i, 2]
        elif rd[j].is_leaf():
            # if j is a leaf, we calculate the distances to all
            # subtrees of k
            kwns = [kwn for kwn in M if kwn[0] == k]
            for k,w,n in kwns:
                M[v, j, n] = M[v, n, j] = M[k, w, n] + dists[j,w]
                M[v, j, w] = M[v, w, j] = M[k, w, n] + dists[j,n]
        elif rd[k].is_leaf():
            # symmetrically if k is a leaf
            jums = [jum for jum in M if jum[0] == j]
            for j,u,m in jums:
                M[v, m, k] = M[v, k, m] = M[j, u, m] + dists[k,u]
                M[v, u, k] = M[v, k, u] = M[j, u, m] + dists[k,m]
        else:
            # neither j nor k are leaves, so we consider combinations of subtrees
            LL,LR = rd[j].left.pre_order(), rd[j].right.pre_order()
            RL,RR = rd[k].left.pre_order(), rd[k].right.pre_order()

            for (this_L,that_L),(this_R,that_R) in itertools.product(((LL,LR), (LR,LL)),
                                                                     ((RL,RR), (RR,RL))):
                for u,w in itertools.product(this_L, this_R):
                    m_order = sorted(that_L, key=lambda m: M[j, u, m])
                    n_order = sorted(that_R, key=lambda n: M[k, w, n])
                    C = dists[np.ix_(m_order, n_order)].min()
                    Cmin = 1e10
                    for m,n in itertools.product(m_order, n_order):
                        if M[j, u, m] + M[k, w, n] + C >= Cmin:
                            break
                        C = M[j, u, m] + M[k, w, n] + dists[m,n]
                        if C < Cmin:
                            Cmin = C

                    M[v, u, w] = M[v, w, u] = Cmin

    return M


def optimal_ordering(Z, dists):
    # Z - linkage matrix
    # dists - the distance matrix

    # get the tree and a list of handles to its leaves
    tree,rd = hierarchy.to_tree(Z, True)

    # Generate scores
    M = optimal_scores(Z, rd, dists)
    # re-order the tree accordingly
    order_tree(Z, rd, M)

    # new leaf ordering
    row_reorder = tree.pre_order()

    return row_reorder


def plot_leaf_ordering(X, method, metric):
    dists = distance.squareform(distance.pdist(X, metric=metric))
    dists2 = distance.squareform(distance.pdist(X.T, metric=metric))

    Z = hierarchy.linkage(X, method=method, metric=metric)
    Z2 = hierarchy.linkage(X.T, method=method, metric=metric)

    t,rd = hierarchy.to_tree(Z, True)
    t2,rd2 = hierarchy.to_tree(Z2, True)

    M = optimal_scores(Z, rd, dists)
    order_tree(Z, rd, M)
    M2 = optimal_scores(Z2, rd2, dists2)
    order_tree(Z2, rd2, M2)

    rr = t.pre_order()
    rr2 = t2.pre_order()

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(8,8))
    gs = GridSpec(2, 2, top=0.95, bottom=0.05, left=0.05, right=0.95,
                  hspace=0.01, wspace=0.01,
                  width_ratios=(1,3), height_ratios=(1,3))

    ax01 = fig.add_subplot(gs[0,1])
    ax10 = fig.add_subplot(gs[1,0])
    ax11 = fig.add_subplot(gs[1,1])

    hierarchy.dendrogram(Z2, ax=ax01)
    ax01.set_axis_off()
    hierarchy.dendrogram(Z, orientation='right', ax=ax10)
    ax10.set_axis_off()

    ax11.matshow(X[np.ix_(rr,rr2)], cmap="Blues", aspect="auto")
    ax11.tick_params(**{s:'off' for s in ('top', 'bottom', 'right')})
    ax11.tick_params(labeltop='off', labelleft='off', labelright='on')

    ax11.set_xticks(np.arange(len(rr2)))
    ax11.set_xticklabels(rr2, fontsize=5.0)
    ax11.set_yticks(np.arange(len(rr)))
    ax11.set_yticklabels(rr, fontsize=5.0)

    plt.show()



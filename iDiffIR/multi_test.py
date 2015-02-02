#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#
#    Author: Mike Hamilton, Colorado State University, 2015
#
"""
Adjust for multiple testing
"""
import numpy as np
from scipy.interpolate import splrep, splev
import matplotlib
matplotlib.use('agg', warn=0) 
import matplotlib.pyplot as plt
from matplotlib import rc
rc('text', usetex=False)
def bonferroni( pvalues ):
    """
    Bonferroni adjustement: reject H_0 if p_i < alpha / N

    'probability of at least one type-1 error is at most alpha'
    """
    N = len(pvalues)
    return [ min( p*N, 1.0) for p in pvalues ]

def bh( pvalues ):
    """
    Benjamini-Hochberg adjustment:

    reject p_r, r=1...k, where k = max_i: p_i < i * alpha / m 

    'FDR <= # true null hypotheses / all tests * alpha'
    
    **preserves original order
    """
    order = np.argsort( np.array( pvalues ) )
    p = sorted(pvalues)
    m = len(pvalues)
    adj = [0]*m
    adj[m-1] = min(1.0,p[m-1] / (1.0+(m-1)) * m)

    for i in xrange(m-2,-1,-1):
        adj[i] = min( p[i] / (1.0+i) * m, adj[i+1] )
    rtn = [0]*m
    for i in xrange(m):
        rtn[ order[i]] = adj[i]
    return rtn

def qvalues( pvalues, splineDF=3, qplot=False ):
    """
    Compute $q$-values for the given set of $p$-values
    """
    m = len(pvalues)
    p = np.array(sorted(pvalues))
    order = np.argsort(pvalues)

    # Estimate proportion of truly null tests
    lams = np.arange(0, 0.95, 0.05)
    pi0s = [ np.mean(p >= li)/(1.0 - li) for li in lams ]
    spi0 = splrep( lams, pi0s, k = splineDF )
    preds = splev(np.arange(0, 1.05, 0.05), spi0  )

    pi0=min(1.0, splev(numpy.max(lams), spi0))

    qvals = [0]*m
    qvals[m-1] = min(pi0*p[m-1], 1.0)

    
    for i in xrange(m-2,-1,-1):
        qvals[i] = min( qvals[i+1], pi0*m*p[i]/float(i+1))
    if qplot:
        plt.figure()
        plt.plot(lams,pi0s,'.')
        plt.plot(lams,preds[:-2],'-')
        plt.savefig('lamdas.pdf')
    rtn = [0]*m
    for i in xrange(m):
        rtn[order[i]] = qvals[i]
    return rtn, preds,pi0

if __name__ == '__main__':
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
    from matplotlib.ticker import FormatStrFormatter
    from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
    from mpl_toolkits.axes_grid1.inset_locator import mark_inset

    # make some phony p-values
    pvals = np.hstack([np.random.uniform(size=400), np.random.uniform(0,0.001,30)])
    bhAdj = bh(pvals)
    bfAdj = bonferroni(pvals)
    qvAdj,preds,pi0 = qvalues( pvals, 4,qplot=1)
    import bisect
    k=bisect.bisect_right(sorted(pvals),0.0011)

    plt.figure(figsize=(14,7))

    ax = plt.subplot(121)
    ax.grid()
    ax.plot(pvals,bhAdj,'.', label='Benjamini-Hochberg')
    ax.plot(pvals,bfAdj,'.', label='Bonferroni')
    ax.plot(pvals,qvAdj,'.', label='q-values')
    ax.set_xlim(-.1,1.1)
    ax.set_ylim(-.1,1.1)
    ax.legend(loc=4,numpoints=1)
    ax.set_xlabel(r'$p$-value', size=18)
    ax.set_ylabel(r'adj $p$-value', size=18)
    axins = inset_axes(ax,
                   width="40%",
                   height="40%", 
                   loc=5)
    axins.set_ylim( 0,0.1)
    axins.set_xticks([0, sorted(pvals)[k-1]])

    
    majorFormatter = FormatStrFormatter(r'%0.1g')
    axins.xaxis.set_major_formatter(majorFormatter)
    axins.grid()

    axins.plot(sorted(pvals)[:k],sorted(bhAdj)[:k],'.', label='Benjamini-Hochberg')
    axins.plot(sorted(pvals)[:k],sorted(bfAdj)[:k],'.', label='Bonferonni')


    axins.plot(sorted(pvals)[:k],sorted(qvAdj)[:k],'.', label='q-values')
    axins.set_ylim( 0,0.1)
    mark_inset(ax, axins, loc1=2, loc2=3, fc=".01", ec=".5")
    #axins.set_xlim( 0,001)

    plt.xticks([0, sorted(pvals)[k-1]],rotation=45)


    plt.subplot(122)
    h = plt.hist(pvals,normed=1)
    plt.xlabel(r'$p$-value', size=18)
    plt.ylim(0,2)
    plt.plot(np.arange(0,1.05,0.05),preds,lw=2,color='k')
    plt.axhline(pi0,lw=2,ls='--',color='r')
    #ax.set_xticklabels([r'0.0', r'%0.4f'%sorted(pvals)[k-1]])
    plt.savefig('out.pdf')


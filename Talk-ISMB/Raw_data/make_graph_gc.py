import os
import sys
import numpy as np
from matplotlib import rc,use
use('pgf')
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from scipy import integrate

COLORS = ['r','b','g','y','k']

RAW_DATA = 'dump_gc_data.txt'
BIN_WIDTH = 0.005

def get_bin(x,width):
    x = int(x*10000)
    width = int(width*10000)
    bin = x/width
    return bin

def make_list_gc_per_round():
    raw_data = [x.strip() for x in open(RAW_DATA)]
    l_data = []
    for i,line in enumerate(raw_data):
        if i % 2 == 0:
            weight = float(line)
        else:
            gc_contents = line.split()
            l_data.append((weight, gc_contents))
    return l_data

def cluster_data(l_data,width):
    pass

def main():
    MAX_Y =1175 
    l_data = make_list_gc_per_round()
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title(r'Target 15% GC')
    ax.set_xlabel('GC content')
    ax.set_ylabel('Quantity')
    ax.set_ylim((0,MAX_Y))
    ax.set_xlim((0,1))
    ax.vlines(0.10,0,MAX_Y,color='k',linestyles='dashed')
    ax.vlines(0.20,0,MAX_Y,color='k',linestyles='dashed')
    #ax.get_yaxis().set_ticks([])
    total_generated = 0
    
    for i,data in enumerate(l_data):
        print i
        weight = data[0]
        a = np.array([float(x) for x in data[1]])
        sigma,mu = np.std(a),np.average(a)
        x = np.linspace(0,1,100)
        ax.plot(x,130*mlab.normpdf(x,mu,sigma),COLORS[i],
                label= '$\Pr(G\mid C) = %0.4f$' % (2*weight),linewidth=6)
        total_generated += 1000*integrate.quad(lambda x: mlab.normpdf(x,mu,sigma),
                                          0.1,0.2)[0]
        total_generated = int(total_generated)
        #total_generated += max([mlab.normpdf(x,mu,sigma) for x in 
        #                        np.linspace(0.1,0.2,1000)])
        #ax.hlines(total_generated,0.1,0.2,color=COLORS[i],linewidth=2,A
        #linestyle='dashed',label='Tot round %s' % i)
        ax.plot([0.1,0.2],[total_generated,total_generated],'%s--' % COLORS[i],
                linewidth=5,label='Tot round %s is %s' % ((i+1),total_generated))
        ax.fill_between([0.1,0.2],[total_generated]*2,[0]*2,facecolor=COLORS[i],
                       alpha=0.7)
        ax.legend()
        #fig.show()
        fig.savefig('dist_%s.pdf' % i, dpi=1200)
        print total_generated
        
    #a = np.array(l_data[0][])
    #print a
    #l_data_clustered = cluster_data(l_data,BIN_WIDTH)
    #print [x[0] for x in l_data]

if __name__ == '__main__':
    main()

###############################################################################
#Copyright (c) 2013, Vladimir Reinharz, Yann Ponty & Jerome Waldispuhl        #
#All rights reserved.                                                         #
#                                                                             #
#Redistribution and use in source and binary forms, with or without           #
#modification, are permitted provided that the following conditions are met:  #
#* Redistributions of source code must retain the above copyright             #
#notice, this list of conditions and the following disclaimer.                #
#* Redistributions in binary form must reproduce the above copyright          #
#notice, this list of conditions and the following disclaimer in the          #
#documentation and/or other materials provided with the distribution.         #
#* Neither the name of the <organization> nor the                             #
#names of its contributors may be used to endorse or promote products         #
#derived from this software without specific prior written permission.        #
#                                                                             #
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"  #
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE    #
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE   #
#ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY       #
#DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES   #
#(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; #
#LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND  #
#ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT   #
#(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS#
#SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE                  #
###############################################################################
import itertools
import random
import math
import sys
import os
import collections
from Energy import *
from pprint import pprint

def MPMATH_MISSING():
  print """The module `mpmath` was not found. This might impeed the
  processing of long RNAs.
  More information can be found on http://code.google.com/p/mpmath/
  We also recommand installation of `gmpy` (http://code.google.com/p/gmpy/),
  automatically leveraged by `mpmath` to increase the speed of computations"""

try: #For infinite precision
  from mpmath import mpf
except ImportError:
  MPMATH_MISSING()
  def mpf(n):
    return n

sys.setrecursionlimit(10000)

# now all the energy values from Turner 2004 are stored in Energy (stacking, internal, bulge, hairpin)
# moved isostericity values as well to Energy.py
nucleotides = ['A', 'C', 'G', 'U']
IUPACBASES = {
  'A':['A'],
  'C':['C'],
  'G':['G'],
  'U':['U'],
  'N':['A','C','G','U'],
  'R':['A','G'],
  'Y':['C','U'],
  'S':['G','C'],
  'W':['A','U'],
  'K':['G','U'],
  'M':['A','C'],
  'B':['C','G','U'],
  'D':['A','G','U'],
  'H':['A','C','U'],
  'V':['A','C','G']
  }
BASES = []

BOLTZMANN = 0.0019872041
global T
T = 310.15

class memoize(dict):
  """Generically memoizes a function results."""
  fun = None
  
  def __init__(self, f):
      self.fun = f
  
  def __call__(self,profile,ref_seq,struct,*args):
      nargs = (args)
      if nargs in self:
          return self[nargs]
      else:
          val = mpf(self.fun(profile,ref_seq,struct,*args))
          self[nargs] = val
          return val
  def resetCache(self):
      self.clear()

class memoize_iso(dict):
  """Generically memoizes a function results."""
  fun = None
  
  def __init__(self, f):
      self.fun = f
  
  def __call__(self,ref_seq,*args):
      nargs = (args)
      if nargs in self:
          return self[nargs]
      else:
          val = mpf(self.fun(ref_seq,*args))
          self[nargs] = val
          return val
  def resetCache(self):
      self.clear()

################################################################################# ss_decomposition
def nb_strand(nbs):
    nb = 1 
    for x in sorted(nbs)[:-1]:
        if not x+1 in nbs:
            nb +=1 
    return nb


def ss_to_bp(ss):
    bp = []
    l = []
    for i, x in enumerate(ss):
        if x == '(':
            l.append(i)
        elif x == ')':
            bp.append((l.pop(), i))
    return bp


def make_tree(ss):
    def make_sub_tree(start, end, bps):
        t = []
        x = start
        while x <= end:
            flag_aft_bp = False

            if x in bps:
                y, z = x, bps[x]
                while True:
                    if y+1 in bps and z-1 == bps[y+1]:
                        y += 1
                        z -= 1
                    else:
                        t.append((x,y,z,bps[x], make_sub_tree(y+1,z-1,bps)))
                        x = bps[x] + 1
                        flag_aft_bp = True
                        break
            if x > end:
                break
            if not flag_aft_bp:
                t.append(x)
                x += 1
        return t
    s = 0
    e = len(ss)-1
    bps = ss_to_bp(ss)
    d_bps = {x[0]:x[1] for x in bps}
    return [-1] + make_sub_tree(s,e,d_bps) + [len(ss)]


def make_components(tree):
    l_components = []
    this_comp = set()
    sub_comp = []
    flag = False
    if not tree:
        return []
    if tree[0] == -1:
        flag = True
        e = tree[-1]
        tree = tree[1:-1]

    for x in tree:
        if isinstance(x, tuple):
            if not flag or not (x[0] == 0 and x[3] == e-1):
                this_comp.add(x[0])
                this_comp.add(x[3])
            stack = [z for z in range(x[0], x[1]+1)]
            stack.extend(z for z in range(x[2], x[3]+1))
            l_components.append(set(stack))
            sub_comp.extend(make_components(x[-1]))
            sub_comp[-1].add(x[1])
            sub_comp[-1].add(x[2])
        else:
            this_comp.add(x)
    l_components.extend(sub_comp)
    if this_comp:
        l_components.append(this_comp)
    return l_components


def decompose(ss):
    #get the decomposition of the ss in elements
    dec = make_components(make_tree(ss))
    #a list of the bps
    bps = ss_to_bp(ss)
    #the dic needed for the forward recursion
    dic_forward = {}
    for el in dec: #for every element decomposed (a list of positions)
        el = sorted(el) #sort the elements
        if nb_strand(el) == 1: #nb strands (i.e. nb of stretch of consecutive nucleotides) if 1, must be hairpin
            dic_forward[el[0], el[-1]] = ('Hairpin', #type
                                         [(el[0],el[-1])],    #bps inside that element
                                         el)    #list of positions inside it
        elif nb_strand(el) > 2: #if more than 2 strands multiloop (we check for external at the end)
            dic_forward[el[0], el[-1]] = ('MultiLoop',
                                          [x for x in bps if x[0] in el and x[1] in el and ((el[0], el[-1]) != x)], #We want all base pairs inside that multiloop, except the most external one
                                         el)
         
        else: #else we must have 2 strands, of a stem or an interior loop / bulge
            if (len(el) % 2 == 0 and #if we have an even number of nucleotides
                all((el[i], el[-i-1]) in bps for i in range(len(el)/2))): #and they all form base pairs, we have a stem!
                dic_forward[el[0], el[-1]] = ('Stem',
                                             [(el[0], el[-1]), (el[-1+len(el)/2], el[len(el)/2])], # first bp and the last bp
                                             el)

            else:
                #if not a strand, we must be an interior loop, and then we need to record the internal base pair
		dic_forward[el[0], el[-1]] = ('InteriorLoop',
                                             [(el[0], el[-1])]+[x for x in bps if x[0] in el and x[1] in el and ((el[0], el[-1]) != x)],
                                             el)


        #We finally check for the external element, if their is one. The only ambiguous case is if it is a stem.
        if 0 in el and len(ss)-1 in el and dic_forward[el[0], el[-1]][0] != 'Stem':
            dic_forward[el[0], el[-1]] = ('External',
                                         [x for x in bps if x[0] in el and x[1] in el and ((el[0], el[-1]) != x)],
                                         el
                                        )
    return dic_forward

##ss = '(((...(.(....).)..((....))..)))...' # hard-coded at the moment to leave the original inputs as it is to compile and to test. After forward is checked, it would take from the user-provided txt file
ss = '(((...(.(....).)..((....))..)))...'
dic = decompose(ss)

################################################################# end of ss_decomposition

def stacking_energy((a,b),(a2,b2),alpha): # now called stacking energy
  #stacking energy of base pair (a,b) around base pair (a2,b2)
  E = STACKING_ENERGY[a,a2,b2,b]
  return  math.exp(-(alpha*E)/(BOLTZMANN*T))

def stacking_energy_short((a,b), alpha): # the short version of it with pre-computed values that depend on outer pair a,b
  E = STACKING_ENERGY_SHORT[a,b]
  return  math.exp(-(alpha*E)/(BOLTZMANN*T))

# new energy calculations
def hairpin_energy((a,b),(a2,b2),size,alpha):
	#stacking energy of base pair from the loop-closing pair of a-b and while a2-b2 are the first starting unpaired bases of the hairpin loop. (called terminal mismatch)
	#Then the initiation parameter (destablizing energies by loop size) is included
	if (size < 3) or (size > 30):
		print "The hairpin loop size has to be between 3 ~ 30"
		#return something like max_int to know the error in further stage
		#or maybe just 1 to exclude from affecting the energy model ('result' multiplication)
          	sys.exit(1)
	else:
		E = TERMINAL_MISMATCH[a,a2,b2,b] + HAIRPIN_INITIATION[size-3] # since the first entry starts with size 0
		return math.exp(-(alpha*E)/(BOLTZMANN*T))
def hairpin_energy_short((a,b), size, alpha): # the short pre-computed version so only the loop-closing pair a-b needed
	if (size < 3) or (size > 30):
		print "The hairpin loop size has to be between 3 ~ 30"
		sys.exit(1)
	else:
		E = TERMINAL_MISMATCH_SHORT[a,b] + HAIRPIN_INITIATION[size-3] # since the first entry starts with size 0
		return math.exp(-(alpha*E)/(BOLTZMANN*T))
def bulge_energy(size, alpha):
	if (size < 1) or (size > 30):
		print "The bulge size has to be between 1 ~ 30"
		sys.exit(1)
	E = BULGE_INITIATION[size-1]
	return  math.exp(-(alpha*E)/(BOLTZMANN*T))

def internal_loop_1x1_energy(a,b,c,d,e,f, alpha):
	# internal loop where a-f and c-d are the base pairs, and b,e form the 1x1 loop
	E = INTERNAL_ENERGY[a,b,c,d,e,f]
	return math.exp(-(alpha*E)/(BOLTZMANN*T))
def internal_loop_1x2_energy(a,b,c,d,e,f,g, alpha):
	E = INTERNAL_ENERGY[a,b,c,d,e,f,g]
	return math.exp(-(alpha*E)/(BOLTZMANN*T))
def internal_loop_2x2_energy(a,b,c,d,e,f,g,h, alpha):
	E = INTERNAL_ENERGY[a,b,c,d,e,f,g,h]
	return math.exp(-(alpha*E)/(BOLTZMANN*T))

def internal_loop_energy(size, alpha):
	if (size < 4) or (size > 30):
		print "The internal loop size has to be between 4 ~ 30"
		sys.exit(1)
	E = OTHER_INTERNAL_ENERGY[size - 4]
	return  math.exp(-(alpha*E)/(BOLTZMANN*T))

def dangling_energy_short((a,b), case, alpha):
	if case == 0: E = THREE_PRIME_DANGLING_END_SHORT[(a,b)]
	if case == 1: E = FIVE_PRIME_DANGLING_END_SHORT[(a,b)]
	if case == 2: E = TERMINAL_MISMATCH_SHORT[(a,b)]
	return math.exp(-(alpha*E)/(BOLTZMANN*T))

def multiloop_energy(asymmetry, num_helices, alpha): # cst parameters from Turner04
	E = 9.25 + 0.91*(min(2.0, asymmetry)) + (-0.63)*(num_helices)
	return math.exp(-(alpha*E)/(BOLTZMANN*T))

@memoize_iso
def isostericity(ref_seq,(i,j),(a,b), alpha):
  #isostericity of going from original base pair to (a,b)
  if not ref_seq:
    return 1
  iso = sum(ISO[(ref[i],ref[j]),(a,b)] for ref in ref_seq)/len(ref_seq)
  #iso_start = sum(ISO[(ref[i],ref[j]),(seq[i],seq[j])] for ref in ref_seq)
  #iso = mpf(iso_mut-iso_start)/len(ref_seq)
  return  math.exp(-((1-alpha)*iso)/(BOLTZMANN*T))

forward_dict = dic
@memoize
def forward(profile,ref_seq,struct,(i,j),(a,b),alpha): # currently only using i,j
  #alpha gives the weight energy vs isostericity
  result = 0.
  if i > j :
    result=1.
  else:
	k = struct[i]
	#l = struct[j]
	if k==-1: # if not paired, meaning it is '.'
	#pro = profile[i][a2]
		z = 0
		if -1 not in struct[i:len(struct)-1]:
			for tmp in range(i, len(struct)):
				if struct[i] != -1:
					z = tmp
					break
		f = forward(profile,ref_seq,struct,
                  (z,j), #new z goes here, where the next paired one starts
                  (a,b),
                  alpha)
      		result += f#*pro
      
 	else: # if paired
	    elem = []
	    for x in forward_dict:
		if (i,k) == forward_dict[x][1][0]:
			elem.append(forward_dict[x])
	    #print "ELEMENTS", elem, i,k
	    if elem:
		    for x in range(len(elem)):
			if elem[x][0] == 'Hairpin': # Hairpin
				for a in nucleotides:
					for b in nucleotides:
						if (a,b) in TERMINAL_MISMATCH_SHORT: # (a,b) is the loop-closing pair. So need to be paired (A:U, G:C, etc)
							result += hairpin_energy_short((a,b), len(elem[x][2])-2, alpha)

			elif elem[x][0] == 'Stem': # Stem
				p = 0
				q = len(elem[x][2])-1
				while ((q-p) > 1): # the list of the paired ones are in pairs by definitnion. So, multiple of 2. Index of the last paired differ by 1
					for a in nucleotides:
						for b in nucleotides:
							if (a,b) in STACKING_ENERGY_SHORT:
								result += stacking_energy_short((a,b), alpha)
					p = p+1
					q = q-1
				i1 = elem[x][1][1][0]
				j2 = elem[x][1][1][1]
				f = forward(profile,ref_seq,struct,(i1,j2),(a,b),alpha)
				result += f

			elif elem[x][0] == 'External': # External
				if (elem[x][1][0][0] == min(elem[x][2]) and elem[x][1][0][1] != max(elem[x][2])): # if dangling end at 3' prime side
					for a in nucleotides:
						for b in nucleotides:
							if (a,b) in THREE_PRIME_DANGLING_END_SHORT:
								result += dangling_energy_short((a,b), 0, alpha)
				elif (elem[x][1][0][0] != min(elem[x][2]) and elem[x][1][0][1] == max(elem[x][2])): # if dangling end at 5' prime side
					for a in nucleotides:
						for b in nucleotides:
							if (a,b) in FIVE_PRIME_DANGLING_END_SHORT:
								result += dangling_energy_short((a,b), 1, alpha)
				elif (elem[x][1][0][0] != min(elem[x][2]) and elem[x][1][0][1] != max(elem[x][2])): # if dangling end at both side, terminal mismatch 
					for a in nucleotides:
						for b in nucleotides:
							if (a,b) in TERMINAL_MISMATCH_SHORT:
								result += dangling_energy_short((a,b), 2, alpha)
				f = forward(profile,ref_seq,struct,(i+1,j),(a,b),alpha)
				result += f
			elif elem[x][0] == 'MultiLoop': # Multiloop
				pairs_list= elem[x][1]
				unpaired_size=0.0
				for n in range(len(elem[x][2])):
					if struct[elem[x][2][n]] == -1:
						unpaired_size += 1
				unpaired_size = unpaired_size/len(pairs_list)
				
				'''for n in range(0, len(pairs_list)):
					if n == 0: # if first nucleotide of the paired ones			
						unpaired_size.append(elem[x][2].index(pairs_list[n]))
					if n == (len(pairs_list)-1): # if last nucleotide of the paired ones
						unpaired_size.append(len(elem[x][2]) - elem[x][2].index(pairs_list[n]) - 1)
					elif ((n % 2) == 1): # we want to check the gap between the second nucleotide of the pair and the first nucleotide of next pair
						unpaired_size.append(elem[x][2].index(pairs_list[n+1]) - elem[x][2].index(pairs_list[n]) - 1)'''
				
				for m in range(len(pairs_list)):
					i1 = pairs_list[m][0]
					j1 = pairs_list[m][1]
					f = forward(profile,ref_seq,struct,(i1+1,j1+1),(a,b),alpha)
					result += f
				result += multiloop_energy(unpaired_size , len(elem[x][1])+1.0, alpha)
			elif elem[x][0] == 'InteriorLoop':  #Interior loops
				start_first = elem[x][1][0][0]
				start_second = elem[x][1][0][1]
				end_first = elem[x][1][1][0]
				end_second = elem[x][1][1][1]
				n1 = end_first - start_first -1
				n2 = end_second - start_second -1
				if (n1 > 0 and n2 > 0): # internal loop
					e = 0.6*abs(n1-n2)*internal_loop_energy((n1+n2),alpha) # TODO
				else: # bulges
					e = bulge_energy(max(n1,n2),alpha)
				f = forward(profile,ref_seq,struct,(end_first,end_second),(a,b),alpha)
				#back = backward()			
				result += e+f
  return result

flip = dic.items()
flip.reverse()
back_dict = collections.OrderedDict(flip)

@memoize
def backward(profile,ref_seq,struct,(i,j),(a,b),alpha):
  #print i,j,a,b, "backward"
  result = 0.
  if i<0: # same as before, may be adjusted when the input parameters are modified
    result=forward(profile,ref_seq,struct,
                   (j,len(struct)-1),
                   ('X','X'),
                   alpha)
  else: 
    k = struct[i]
    if k==-1:
	back = backward(profile,ref_seq,struct,
                        (i-1,j),
                        (a,b),
                        alpha)
	result += back

    else:
	elem = []
	for x in back_dict:
		if (i,k) == back_dict[x][1][-1]:
			elem.append(back_dict[x])
	if elem:
		for x in range(len(elem)):
			if elem[x][0] == 'Hairpin': # Hairpin
				f = 0.
				for a in nucleotides:
					for b in nucleotides:
						if (a,b) in TERMINAL_MISMATCH_SHORT: # (a,b) is the loop-closing pair. So need to be paired (A:U, G:C, etc)
							f += hairpin_energy_short((a,b), len(elem[x][2])-2, alpha)
				back = backward(profile,ref_seq,struct,(i-1,j+1),(a,b),alpha)
				result += back
			elif elem[x][0] == 'Stem': # Stem
				e = 0				
				p = 0
				q = len(elem[x][2])-1
				while ((q-p) > 1): # the list of the paired ones are in pairs by definitnion. So, multiple of 2. Index of the last paired differ by 1
					for a in nucleotides:
						for b in nucleotides:
							if (a,b) in STACKING_ENERGY_SHORT:
								e += stacking_energy_short((a,b), alpha)
					p = p+1
					q = q-1
				i1 = elem[x][1][0][0]
				j1 = elem[x][1][0][1]
				forw = forward(profile,ref_seq,struct,(i,k),(a,b),alpha)
				back = backward(profile,ref_seq,struct,(i1,j1),(a,b),alpha)
				
				result += forw+e+back

			elif elem[x][0] == 'External': # External
				e = 0
				if (elem[x][1][0][0] == min(elem[x][2]) and elem[x][1][0][1] != max(elem[x][2])): # if dangling end at 3' prime side
					for a in nucleotides:
						for b in nucleotides:
							if (a,b) in THREE_PRIME_DANGLING_END_SHORT:
								e += dangling_energy_short((a,b), 0, alpha)
				elif (elem[x][1][0][0] != min(elem[x][2]) and elem[x][1][0][1] == max(elem[x][2])): # if dangling end at 5' prime side
					for a in nucleotides:
						for b in nucleotides:
							if (a,b) in FIVE_PRIME_DANGLING_END_SHORT:
								e += dangling_energy_short((a,b), 1, alpha)
				elif (elem[x][1][0][0] != min(elem[x][2]) and elem[x][1][0][1] != max(elem[x][2])): # if dangling end at both side, terminal mismatch 
					for a in nucleotides:
						for b in nucleotides:
							if (a,b) in TERMINAL_MISMATCH_SHORT:
								e += dangling_energy_short((a,b), 2, alpha)
				forw = forward(profile,ref_seq,struct,(i+1,j-1),(a,b),alpha)
				# no backward				
				result += e+forw
			elif elem[x][0] == 'MultiLoop': # Multiloop
				pairs_list=elem[x][1]
				print pairs_list, " is the pair list"
				unpaired_size=[]
				'''
				#for pairs in elem[x][1]:
				#	pairs_list.append(pairs)
				for n in range(0, len(pairs_list)):
					if n == 0: # if first nucleotide of the paired ones			
						unpaired_size.append(elem[x][2].index(pairs_list[n]))
					if n == (len(pairs_list)-1): # if last nucleotide of the paired ones
						unpaired_size.append(len(elem[x][2]) - elem[x][2].index(pairs_list[n]) - 1)
					elif ((n % 2) == 1): # we want to check the gap between the second nucleotide of the pair and the first nucleotide of next pair
						unpaired_size.append(elem[x][2].index(pairs_list[n+1]) - elem[x][2].index(pairs_list[n]) - 1)'''
				for m in range(len(pairs_list)):
					i1 = pairs_list[m][0]
					j1 = pairs_list[m][1]
					f = forward(profile,ref_seq,struct,(i1,j1),(a,b),alpha)
					result += f
				back = backward(profile,ref_seq,struct,(i-1,j+1),(a,b),alpha)
				result += back
			elif elem[x][0] == 'InteriorLoop':  #Interior loops
				start_first = elem[x][1][0][0]
				start_second = elem[x][1][0][1]
				end_first = elem[x][1][1][0]
				end_second = elem[x][1][1][1]
				n1 = end_first - start_first -1
				n2 = end_second - start_second -1
				if (n1 > 0 and n2 > 0): # internal loop
					e = 0.6*abs(n1-n2)*internal_loop_energy((n1+n2),alpha) # TODO
				else: # bulges
					e = bulge_energy(max(n1,n2),alpha)
				f = forward(profile,ref_seq,struct,(end_first,end_second),(a,b),alpha)
				back = backward(profile,ref_seq,struct,(start_first-1,start_second+1),(a,b),alpha)		
				result += e+f+back

  return result


"""k = struct[i]
    
    if k==-1: # if the current nucleotides it not paired, meaning it is '.'
      print 'b'
      for a2 in BASES[i]:
        pro = profile[i][a2]
        back = backward(profile,ref_seq,struct,
                        (i-1,j),
                        (a2,b),
                        alpha)
        result += pro*back
    #BP to the left
    elif k<i:
      print 'c'
      for a2 in BASES[k]:
        for b2 in BASES[i]:
          pro = profile[k][a2]*profile[i][b2]
          back = backward(profile,ref_seq,struct,
                       (k-1,j),
                       (a2,b),
                       alpha)
          forw = forward(profile,ref_seq,struct,
                      (k+1,i-1),
                      (a2,b2),
                      alpha)
          iso = isostericity(ref_seq,
                             (k,i),
                             (a2,b2),
                             alpha)
          result += pro*back*forw*iso
    #BP to the right
    elif k>=j:
      print 'd'
      for a2 in BASES[i]:
        for b2 in BASES[k]:
          pro = profile[i][a2]*profile[k][b2]
          #No stack
          if not (j==k and struct[i+1]==j-1):
            back = backward(profile,ref_seq,struct,
                         (i-1,k+1),
                         (a2,b2),
                         alpha)
            forw = forward(profile,ref_seq,struct,
                        (j,k-1),
                        (b,b2),
                        alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result += pro*back*forw*iso
          #Stack
          else:
	    
            back = backward(profile,ref_seq,struct,
                         (i-1,k+1),
                         (a2,b2),
                         alpha)
            e = stacking_energy((a2,b2), # now called stacking energy
                       (a,b),
                       alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result += pro*back*e*iso"""

def parseStruct(dbn):
  p = []
  result = [-1 for c in dbn]
  for i in range(len(dbn)):
    c = dbn[i]
    if c=='(':
      p.append(i)
    elif c==')':
      j = p.pop()
      result[j] = i
      result[i] = j
  return result

def product_given_i(profile,ref_seq,struct,i,a,alpha):
  """Will compute the sum of boltzmann weights of structures 
  where the 'i-th' nucleotide is 'a'.
  """
  print "i is : ", i
  n = len(struct)
  tot = forward(profile,ref_seq,struct,(0,n-1),('X', 'X'),alpha)
  k = struct[i]
  result = mpf(0)
  if k == -1:
    pro = profile[i][a]
    result += pro*backward(profile,ref_seq,struct,(i-1,i+1),(a,a),alpha)
  elif k < i:
    for c in BASES[k]:
      pro = profile[k][c]*profile[i][a]
	  
      f = forward(profile,ref_seq,struct,(k+1,i-1),(c,a),alpha) 
      b = backward(profile,ref_seq,struct,(k-1,i+1),(c,a),alpha)
      
      iso = isostericity(ref_seq,(k,i),(c,a),alpha)
      result += pro*f*b*iso
  else:
    for c in BASES[k]:
      pro = profile[i][a]*profile[k][c]
      f = forward(profile,ref_seq,struct,(i+1,k-1),(a,c),alpha) 
      b = backward(profile,ref_seq,struct,(i-1,k+1),(a,c),alpha)
      iso = isostericity(ref_seq,(i,k),(a,c),alpha)
      result += pro*f*b*iso
  return result

def random_weighted_sampling(l_samples):
  tot = sum(x[1] for x in l_samples)
  scaled_weights = [x[1]/tot for x in l_samples]
  rand_nb = random.random()
  accumulation = 0
  for i,x in enumerate(scaled_weights):
    accumulation += x
    if accumulation > rand_nb:
      return l_samples[i][0]
  return l_samples[-1][0]


def backtrack(profile,ref_seq,struct,(i,j),(a,b),alpha):
  #alpha gives the weight energy vs isostericity
  result_list = []
  #max_seq_list = []
  max_seq = ''
  if i > j :
    return '' 
  else:
    k = struct[i]
    if k==-1:
      l_samples = []
      for a2 in BASES[i]:
        pro = profile[i][a2]
        f = forward(profile,ref_seq,struct,
                    (i+1,j),
                    (a2,b),
                    alpha)
        result =  pro*f
        l_samples.append((a2,result))
      a2 = random_weighted_sampling(l_samples)
      max_seq = a2 + backtrack(profile,ref_seq,struct,(i+1,j),(a2,b),alpha)
      
    elif i < k <= j: #If k > j we return 0
      l_samples = []
      for a2 in BASES[i]:
        for b2 in BASES[k]:
          pro = profile[i][a2]*profile[k][b2]
          #if not stacked outside (or border, then no stack possible)
          if i==0 or j==len(struct)-1 or not (j==k and struct[i-1]==j+1):
            f1 = forward(profile,ref_seq,struct,
                        (i+1,k-1),
                        (a2,b2),
                        alpha)
            f2 = forward(profile,ref_seq,struct,
                         (k+1,j),
                         (b2,b),
                         alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result = pro*f1*f2*iso
          #if stack, we add energy
          else :
            f = forward(profile,ref_seq,struct,
                       (i+1,k-1),
                       (a2,b2),
                       alpha)
            e = stacking_energy((a,b), # now called stacking energy
                       (a2,b2),
                       alpha)
            iso = isostericity(ref_seq,
                               (i,k),
                               (a2,b2),
                               alpha)
            result = pro*f*e*iso
          l_samples.append(((a2,b2),result))
      """
      for a2,b2 in result_list:
        for best_1 in backtrack(profile,ref_seq,struct,(i+1,k-1),(a2,b2),alpha):
          for best_2 in backtrack(profile,ref_seq,struct,(k+1,j),(b2,b),alpha): 
            max_seq_list.append(a2+best_1+b2+best_2)
      """
      a2,b2 = random_weighted_sampling(l_samples)
      best_1 = backtrack(profile,ref_seq,struct,(i+1,k-1),(a2,b2),alpha)
      best_2 = backtrack(profile,ref_seq,struct,(k+1,j),(b2,b),alpha)
      max_seq = a2+best_1+b2+best_2
  return max_seq

def probability_given_i(profile,ref_seq,struct,i,a,alpha):
  """Will compute the probability that the 'i-th' nucleotide
  is 'a' over all sequences at 'm' mutations from seq
  """
  n = len(struct)
  tot = forward(profile,ref_seq,struct,(0,n-1),('X', 'X'),alpha)
  result = product_given_i(profile,ref_seq,struct,i,a,alpha)
  if tot == 0:
    print """The total partition function is 0, you might want to increase
    the number of mutations allowed"""
    sys.exit(1)
  return result/tot

def testSingleSequence(profile,ref_seq,struct,alpha):
  forward.resetCache()
  backward.resetCache()
  n = len(struct)

  i = 29
  #print "Given i=%s and m=%s the probability of sequence is:" %(i,m)
  #print "\t", probability_given_i_most_m(ref_seq,struct,i,'C',m,alpha=1.0)
  print "  Forward: \t",forward(profile,ref_seq,struct,(0,n-1),('A','G'),alpha)
  i = n-1
  res = 0
  for j in BASES[i]:
    res += profile[i][j]*backward(profile,ref_seq,struct,(i-1,i+1),(j,j),alpha)
  print "  Backward of nuc %s:\t" % i, res

def test():
  
  seq = "UAUAUAUGUAAUACAACAAACAAUAAAAGGUG"
  dbn = "................................"
  ref_seq = seq
  alpha = 0.2 

  ref_seq = ref_seq*1
  seq = seq*1
  dbn = dbn*1
  struct = parseStruct(dbn)
  profile = tuple({'A':0.25,
                   'C':0.25,
                   'G':0.25,
                   'U':0.25} for x in range(len(seq)))



  testSingleSequence(profile,ref_seq,struct,alpha)

def parse_fasta(file_name):
  #Sequences in the MFE and the target secondary structure
  seq = []
  struct = []
  with open(file_name) as f:
    for line in f:
      line = line.strip()
      if not line:
        continue
      if all(x in 'AUGC' for x in line):
        seq.append(line)
        continue
      if all(x in '(.)' for x in line):
        struct.append(line)
        continue
  return seq, parseStruct(struct[0])

def parse_profile(file_name):
  #Parse a profile in the format of the output of RNApyro
  profile = []
  with open(file_name) as f:
    for l in f:
      l = l.strip().split()
      prob = {'A':mpf(l[0]),
              'C':mpf(l[1]),
              'G':mpf(l[2]),
              'U':mpf(l[3])}
      profile.append(prob)
  return tuple(profile)

def equiprob_profile(n):
      prob = {'A':0.25,
              'C':0.25,
              'G':0.25,
              'U':0.25}
      profile = tuple(prob for _ in xrange(n))
      return profile

def all_probabilities(profile,ref_seq, stuct, alpha):
  n = len(struct)
  results = []
  for i in range(n):
    results.append([])
    for a in IUPACBASES['N']:
      if a not in BASES[i]:
        results[-1].append(0)
      else:
        results[-1].append(
          probability_given_i(profile,ref_seq,struct,i,a,alpha))
  return results

def display_all_probabilities(results):
  for x in results:
    print x[0],x[1],x[2],x[3]

def gc_content(sequence,structure=None):
  if not structure:
    nb = len([x for x in sequence if x in 'GC'])
    return float(nb)/len(sequence)
  else:
    nb = len([x for i,x in enumerate(sequence) if x in 'GC' and 
             structure[i] > -1])
    nb_bp = len([x for x in structure if x > -1])
    return float(nb) / nb_bp

def update_profile(profile,max_bound,min_bound,increase=True):
  """Given a profile, will increase (resp. decrease) by half the probability
  of GC content at all positions. Will equaly distribute the remaining between
  AU
  """
  new_profile = []
  for position in profile:
    if increase:
      up_G = (max_bound-position['G'])/2
      up_C = (max_bound-position['C'])/2
    else:
      up_G = -(position['G']-min_bound)/2
      up_C = -(position['C']-min_bound) /2
    to_remove = up_G+up_C
    new_profile.append(
      {'G':position['G']+up_G,
       'C':position['C']+up_C,
       'A':position['A']-to_remove/2})
    #make sure sums to 1
    new_profile[-1]['U'] = 1.0 - sum(new_profile[-1][x] for x in 'GCA')
  return new_profile
     
def sample_gc_target(profile,ref_seq,struct,alpha,nb_gc_sample,gc_target,
                     file_gc_data,f_gc_only_structure,
                     max_err=0.1,sample_before_update=1000):
  """Will sample the required number of sequences, inside the gc_target

  """
  lower_gc = gc_target-max_err
  upper_gc = gc_target+max_err
  max_bound = 1.0
  min_bound = 0.0
  l_all_sample = []
  l_correct_gc = []

  if file_gc_data:
    l_contents = []
    f = open(file_gc_data,'w')

  while len(l_correct_gc) < nb_gc_sample:
    l_all_sample.append([backtrack(profile,ref_seq,struct,(0,n-1),('',''),alpha)
                         for _ in xrange(sample_before_update)])
    over_under = 0 
    if file_gc_data:
      l_contents[:] = []
    for sample in l_all_sample[-1]:
      if f_gc_only_structure:
        content = gc_content(sample,structure=struct)
      else:
        content = gc_content(sample)
      if file_gc_data:
        l_contents.append(content)
      over_under += (content - gc_target)
      if lower_gc <= content <= upper_gc: 
        l_correct_gc.append(sample)

    if file_gc_data:
      f.write('%s\n' % profile[0]['C']) 
      for c in l_contents[:-1]:
        f.write('%s\t' % c)
      f.write('%s\n' % l_contents[-1])
    if over_under < 0:
      min_bound = profile[0]['C']
      profile = update_profile(profile,max_bound,min_bound) 
      forward.clear()
    elif over_under > 0:
      max_bound = profile[0]['C']
      profile = update_profile(profile,max_bound,min_bound,increase=False)
      forward.clear()

  if file_gc_data:
    f.close()

  return l_correct_gc

def sub_seq_structure(sequence,struct):
  return ''.join(x for i,x in enumerate(sequence) if struct[i] > -1)

def diversity_seq(l_sequence,struct):
  n = len(l_sequence)
  n_set = float(len(set(l_sequence)))
  n_struct_set = float(len(set(sub_seq_structure(x,struct)
                               for x in l_sequence)))
  return n,n_set/n,n_struct_set/n



def help():
  print """
  Required:
    -d <file_path> 
      A file containing the target secondary structure and an optional MSA
    -a <float> 
      The value of alpha, between 0 and 1.
      1 takes only into account the secondary structure, 0 only the MSA

  Optional:
    -m <int> 
      The max penality for an invalid base pair, -1 for infinity
    -b <int> 
      print 'n' stochasticly backtracked sequences
    -no_profile <> 
      Doesn't output a profile
    -s_gc <target_gc> <nb_samples> 
      Sampling sequences with a 0<=target_gc<=1 and a given
      number of samples struct
    -gc_sec_struct <>
      Only the nucleotides with an interaction in the secondary structure
      will be considered for the GC content
    -gc_max_err <float>
      Max error from GC target allowed in sample, default 0.1
    -gc_data <file_name>
      Create a file with the given name, will print to it a first line
      the weight of 'C' and on the next the list of sampled GC content
    -t <float>
      The temperature (default 310.5K)
    -c <IUPAC sequence>
      An IUPAC sequence to constrain the outputed sequences
    -p <file_path> 
      A data file with  starting RNA profile (i.e. every lines contains
       nucleotides probability in order: 'ACGU')

  e.g.
    python IncaRNAtion -d data.txt -a 0.5 -m 20
    python IncaRNAtion -d data.txt -a 1 -m 20 -no_profile -s_gc 0.5 100
    python IncaRNAtion -d data.txt -a 0.5 -m 20 -b 5 -no_profile
    """

if __name__ == "__main__":
  opts = sys.argv
  l = len(opts)
  i = 1

  #Action Flags!!
  f_no_profile = False#Don't do profile
  f_backtrack = False#Do backtrack
  f_sample_gc = False#Sample with given gc content
  f_gc_only_structure = False#Only take into account GC with interactions
  f_gc_max_error = False#Custom max error for gc content, default 0.1
  file_gc_data=None#If a file given, print to it for each iter. weight C, GC content

  while i < l:
    cmd = opts[i]
    if not cmd.startswith('-'):
      i += 1
    else:
      if cmd == '-d':  
        file_name = opts[i+1]
        i += 2
        if not os.path.isfile(file_name):
          help()
          sys.exit(1)
      elif cmd == '-p':
        profile_path = opts[i+1]
        i += 2
        if not os.path.isfile(profile_path):
          help()
          sys.exit(1)
      #alpha
      elif cmd == '-a':
        alpha = opts[i+1]
        i += 2
        try:
          alpha = float(alpha)
        except ValueError:
          help()
          sys.exit(1)
        if not 0 <= alpha <= 1:
          help()
          sys.exit(1)
      #Max BP penality (def inf.)
      elif cmd == '-m':
        z = opts[i+1]
        i += 2
        try:
          z = float(z)
          if z >= 0:
            for x in STACKING_ENERGY:
              if STACKING_ENERGY[x] == sys.maxint:
                STACKING_ENERGY[x] = z
        except ValueError:
          pass
      #Don't output profile
      elif cmd == '-no_profile':
        f_no_profile = True
        i += 1
      #Backtrack N optimal sequences
      elif cmd == '-b':
        f_backtrack = True
        nb_backtrack = opts[i+1]
        i += 2
        try:
          nb_backtrack = int(nb_backtrack)
        except ValueError:
          print "Error -b"
          help()
          sys.exit(1)
        if nb_backtrack < 1:
          print "Error -b"
          help()
          sys.exit(1)
      #Sampling Sequences given GC target
      elif cmd == '-s_gc':
        f_sample_gc = True
        gc_target = opts[i+1]
        nb_gc_sample = opts[i+2]
        i += 3
        try:
          nb_gc_sample = int(nb_gc_sample)
          gc_target = float(gc_target)
        except ValueError:
          print "Error s_gc"
          help()
          sys.exit(1)
        if not (0 <= gc_target <= 1) or not nb_gc_sample > 0:
          print "Error s_gc"
          help()
          sys.exit(1)
      #Only check sec struct pos for GC content
      elif cmd == '-gc_sec_struct':
        f_gc_only_structure = True
        i += 1
      #Temperature of the system
      elif cmd == '-t':
        T = opts[i+1]
        i += 2
        try:
          T = float(T)
        except ValueError:
          help()
          sys.exit(1)
      #Max error from GC_Target allowedstruct
      elif cmd == '-gc_max_err':
        f_gc_max_error = True
        max_error = opts[i+1]
        i += 2
        try:
          max_error = float(max_error)
        except ValueError:
          print "Error gc_max_err"
          help()
          sys.exit(1)
        if not 0 <= max_error <= 1:
          print "Error gc_max_err"
          help()
          sys.exit(1)
      #If want data from GC sampling
      elif cmd == '-gc_data':
        file_gc_data = opts[i+1]
        i += 2
      elif cmd == '-c':
        iupac = opts[i+1]
        i+=2 
      else:
        print "Unrecognized arg"
        help()
        sys.exit(1)



  try:
    ref_seq,struct = parse_fasta(file_name)
  except NameError:
    help()
    sys.exit(1)

  n = len(struct)

  #iupac
  try:
    if len(iupac) != n:
      print "IUPAC code too short"
      help()
      sys.exit(1)
    for x in iupac.upper():
      try:
        BASES.append(IUPACBASES[x])
      except KeyError:
        print "Unrecognized IUPAC symbol"
        help()
        sys.exit(1)
  except NameError:
    for _ in range(n):
      BASES.append(IUPACBASES['N'])



  try:
    profile = parse_profile(profile_path)
  except NameError:
    profile = equiprob_profile(n)

  #Action!

  if not f_no_profile:
    results = all_probabilities(profile,ref_seq,struct,alpha)
    display_all_probabilities(results)

  if f_backtrack:
    for _ in xrange(nb_backtrack):
      res = backtrack(profile,ref_seq,struct,(0,n-1),('',''),alpha)
      print res, gc_content(res,struct)

  if f_sample_gc:
    if f_gc_max_error:
      res = sample_gc_target(profile,ref_seq,struct,alpha,nb_gc_sample,
                             gc_target,file_gc_data,f_gc_only_structure,
                             max_err=max_error)
    else:
      res = sample_gc_target(profile,ref_seq,struct,alpha,
                             nb_gc_sample,gc_target,file_gc_data,
                             f_gc_only_structure)
      print '\n'.join(res)

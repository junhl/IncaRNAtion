import os
import math

from ss_decomposition import *
from Energy import *

nucleotides = ['A', 'C', 'G', 'U']
BOLTZMANN = 0.0019872041
global T
T = 310.15

def stacking_energy((a,b),(a2,b2),alpha): # now called stacking energy
  #stacking energy of base pair (a,b) around base pair (a2,b2)
  E = STACKING_ENERGY[a,a2,b2,b]
  return E
  #return  math.exp(-(alpha*E)/(BOLTZMANN*T))

def stacking_energy_short((a,b), alpha): # the short version of it with pre-computed values that depend on outer pair a,b
  E = STACKING_ENERGY_SHORT[a,b]
  return E
  #return  math.exp(-(alpha*E)/(BOLTZMANN*T))

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
                return E
		#return math.exp(-(alpha*E)/(BOLTZMANN*T))
def hairpin_energy_short((a,b), size, alpha): # the short pre-computed version so only the loop-closing pair a-b needed
        if (size < 3) or (size > 30):
                print "The hairpin loop size has to be between 3 ~ 30"
                sys.exit(1)
        else:
                E = TERMINAL_MISMATCH_SHORT[a,b] + HAIRPIN_INITIATION[size-3] # since the first entry starts with size 0
                return E
		#return math.exp(-(alpha*E)/(BOLTZMANN*T))
def bulge_energy(size, alpha):
        if (size < 1) or (size > 30):
                print "The bulge size has to be between 1 ~ 30"
                sys.exit(1)
        E = BULGE_ENERGY[size-1]
        return E
	#return  math.exp(-(alpha*E)/(BOLTZMANN*T))

def internal_loop_1x1_energy(a,b,c,d,e,f, alpha):
        # internal loop where a-f and c-d are the base pairs, and b,e form the 1x1 loop
        E = INTERNAL_ENERGY[a,b,c,d,e,f]
        return E
	#return math.exp(-(alpha*E)/(BOLTZMANN*T))
def internal_loop_1x1_energy(a,b,c,d,e,f, alpha):
        # internal loop where a-f and c-d are the base pairs, and b,e form the 1x1 loop
        E = INTERNAL_ENERGY[a,b,c,d,e,f]
        return E
	#return math.exp(-(alpha*E)/(BOLTZMANN*T))
def internal_loop_1x2_energy(a,b,c,d,e,f,g, alpha):
        E = INTERNAL_ENERGY[a,b,c,d,e,f,g]
        return E
	#return math.exp(-(alpha*E)/(BOLTZMANN*T))
def internal_loop_2x2_energy(a,b,c,d,e,f,g,h, alpha):
        E = INTERNAL_ENERGY[a,b,c,d,e,f,g,h]
        return E
	#return math.exp(-(alpha*E)/(BOLTZMANN*T))

def internal_loop_energy(size, alpha):
        if (size < 4) or (size > 30):
                print "The internal loop size has to be between 4 ~ 30"
                sys.exit(1)
        E = OTHER_INTERNAL_ENERGY[size - 4]
        return E
	#return  math.exp(-(alpha*E)/(BOLTZMANN*T))

def dangling_energy_short((a,b), case, alpha):
        if case == 0: E = THREE_PRIME_DANGLING_END_SHORT[(a,b)]
        if case == 1: E = FIVE_PRIME_DANGLING_END_SHORT[(a,b)]
        if case == 2: E = TERMINAL_MISMATCH_SHORT[(a,b)]
        return E
	#return math.exp(-(alpha*E)/(BOLTZMANN*T))

def multiloop_energy(asymmetry, num_helices, alpha): # cst parameters from Turner04
        E = 9.25 + 0.91*(min(2.0, asymmetry)) + (-0.63)*(num_helices)
        return E
	#return math.exp(-(alpha*E)/(BOLTZMANN*T))

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
	print "elem : ", elem
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


ss = '((.....))'
d = decompose(ss)
print ss
print d
print backward(d,0)

ss1 = '((....))....((...))'
d1 = decompose(ss1)
print ss1
print d1
print backward(d1,0)

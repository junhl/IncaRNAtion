import re
import itertools
import sys

nucleotides = ['A', 'C', 'G', 'U']



#copy from original Incarnation
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
STACKING_ENERGY = {k:sys.maxint for k in itertools.product(
                  IUPACBASES['N'], repeat=4)}
#Adjust with turner04
#The order of the nucleotides is from 5' -> 3'
STACKING_ENERGY.update({('A', 'A', 'U', 'U'):-0.9,
                        ('A', 'C', 'G', 'U'):-2.2,
                        ('A', 'G', 'C', 'U'):-2.1,
                        ('A', 'G', 'U', 'U'):-0.6,
                        ('A', 'U', 'A', 'U'):-1.1,
                        ('A', 'U', 'G', 'U'):-1.4,
                        ('C', 'A', 'U', 'G'):-2.1,
                        ('C', 'G', 'U', 'G'):-1.4,
                        ('C', 'U', 'A', 'G'):-2.1,
                        ('C', 'U', 'G', 'G'):-2.1,
                        ('G', 'A', 'U', 'C'):-2.4,
                        ('G', 'C', 'G', 'C'):-3.4,
                        ('G', 'G', 'C', 'C'):-3.3,
                        ('G', 'G', 'U', 'C'):-1.5,
                        ('G', 'U', 'A', 'C'):-2.2,
                        ('G', 'U', 'G', 'C'):-2.5,
                        ('G', 'A', 'U', 'U'):-1.3,
                        ('G', 'C', 'G', 'U'):-2.5,
                        ('G', 'G', 'C', 'U'):-2.1,
                        ('G', 'G', 'U', 'U'):-0.5,
                        ('G', 'U', 'A', 'U'):-1.4,
                        ('G', 'U', 'G', 'U'):1.3,
                        ('U', 'A', 'U', 'A'):-1.3,
                        ('U', 'C', 'G', 'A'):-2.4,
                        ('U', 'G', 'C', 'A'):-2.1,
                        ('U', 'G', 'U', 'A'):-1.0,
                        ('U', 'U', 'A', 'A'):-0.9,
                        ('U', 'U', 'G', 'A'):-1.3,
                        ('U', 'A', 'U', 'G'):-1.0,
                        ('U', 'C', 'G', 'G'):-1.5,
                        ('U', 'G', 'C', 'G'):-1.4,
                        ('U', 'G', 'U', 'G'):0.3,
                        ('U', 'U', 'A', 'G'):-0.6,
                        ('U', 'U', 'G', 'G'):-0.5,
                        ('C', 'C', 'G', 'G'):-3.3,
                        ('C', 'G', 'C', 'G'):-2.4})



internal_1x1 = {}
internal_1x2 = {}
internal_2x2 = {}
internal = {}
other_internal = []

before_hairpin = {}
hairpin = []

bulge_initiation = []


# Internal loop energies
with open ('int11.txt') as int11:
	lines = int11.readlines()
	energy_values = []
	for i in range (0,len(lines)):
		if (i == 28 or i ==43 or i ==59 or i ==74 or i ==89 or i ==104): # the text file decided to troll with jumping by 15 OR 16 randomly
			five_prime = re.findall('[A,C,G,U]\s[A,C,G,U]', lines[i])
			#print five_prime
			three_prime = re.findall('[A,C,G,U]\s[A,C,G,U]', lines[i+1])			
			#print three_prime
			
			first_line_energy = re.findall('-?[0-9][\.][0-9]', lines[i+4])
			second_line_energy = re.findall('-?[0-9][\.][0-9]', lines[i+5])
			third_line_energy = re.findall('-?[0-9][\.][0-9]', lines[i+6])
			fourth_line_energy = re.findall('-?[0-9][\.][0-9]', lines[i+7])
			energy_values = [] + first_line_energy + second_line_energy + third_line_energy + fourth_line_energy
			#print energy_values
			#print len(energy_values)
			n = 0
			for x in nucleotides:
				for z in range (0,len(five_prime)):
					for y in nucleotides:
						#print (five_prime[0][0], x, five_prime[z][2], three_prime[z][2], y, three_prime[0][0])
						internal_1x1[five_prime[0][0], x, five_prime[z][2], three_prime[z][2], y, three_prime[0][0]] = float(energy_values[n])
						n = n+1
int11.close()

	
with open ('int21.txt') as int21:
	lines = int21.readlines()
	energy_values = []
	l = {26, 40, 54, 68, 82, 96, 110, 124, 138, 152, 167, 182, 196, 211, 226, 241, 256, 271, 286, 301, 316, 331, 346, 361}
	for  i in range (0, len(lines)):
		if i in l:
			five_prime = re.findall('[A,C,G,U]\s{2}[A,C,G,U]', lines[i])
			three_prime = re.findall('[A,C,G,U]\s{2}[A,C,G,U]', lines[i+1])
			nuc_2_side = re.search('[A,C,G,U]', lines[i+2]).group() #the nucleotide next to Y, which is on 2 side in 2x1
			
			first_line_energy = re.findall('-?[0-9][\.][0-9]', lines[i+4])
			second_line_energy = re.findall('-?[0-9][\.][0-9]', lines[i+5])
			third_line_energy = re.findall('-?[0-9][\.][0-9]', lines[i+6])
			fourth_line_energy = re.findall('-?[0-9][\.][0-9]', lines[i+7])
			energy_values = [] + first_line_energy + second_line_energy + third_line_energy + fourth_line_energy

			n = 0
			for x in nucleotides:
				for z in range(0,len(five_prime)):
					for y in nucleotides:
						internal_1x2[five_prime[0][0], x, five_prime[z][3], three_prime[z][3], nuc_2_side, y, three_prime[0][0]] = float(energy_values[n])
						n=n+1
int21.close()


with open ('int22.txt') as int22:
	lines = int22.readlines()
	energy_values = []
	l = []
	for j in range(38, len(lines)): # the spacing is regular for once..
		if (((j-38.0)%26.0) == 0.0):
			l.append(j)
	for i in range (0, len(lines)):
		if i in l:
			#print lines[i]
			five_prime = re.search('[A,C,G,U].*[A,C,G,U]', lines[i]).group()
			three_prime = re.search('[A,C,G,U].*[A,C,G,U]', lines[i+1]).group()
			energy_values = re.findall('-?[0-9][\.][0-9]', ''.join(lines[i+3:i+19]))

			n = 0
			for x1 in nucleotides:
				for x2 in nucleotides:
					for y1 in nucleotides:
						for y2 in nucleotides:
							internal_2x2[five_prime[0], x1, y1, five_prime[9], three_prime[9], y2, x2, three_prime[0]] = float(energy_values[n])
							n=n+1
int22.close()

# combining into one
internal.update(internal_1x1)
internal.update(internal_1x2)
internal.update(internal_2x2)

#if not 1x1, 1x2, 2x2,  it goes to length/size of the internal loop
other_internal = [1.1, 2.0, 2.0, 2.1, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 2.9, 3.0, 3.1, 3.1, 3.2, 3.3, 3.3, 3.4, 3.4, 3.5, 3.5, 3.5, 3.6, 3.6, 3.7, 3.7]
print len(other_internal)

# hairpin
with open ('tstack.txt') as hp:
	 # After the meeting, hairpin has been changed into size dependent list, instead of nucleotide depending to avoid 4^32 entries
	hairpin = [5.4, 5.6, 5.7, 5.4, 6.0, 5.5, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 6.9, 7.0, 7.1, 7.1, 7.2, 7.2, 7.3, 7.3, 7.4, 7.4, 7.5, 7.5, 7.5, 7.6, 7.6, 7.7]
	lines = hp.readlines()
	energy_values = []
	l = [23, 38, 53, 68]
	n = 0
	m = 0
	for i in range (0, len(lines)):
		if i in l:
			five_prime = ['AX', 'CX', 'GX', 'GX', 'UX', 'UX']
			three_prime = ['UY', 'GY', 'CY', 'UY', 'AY', 'GY']
			energy_values = re.findall('-?[0-9][\.][0-9]', ''.join(lines[i+3:i+7]))
			if (n <=1):
        	        	for x in nucleotides:
                                	for y in nucleotides:
                                              before_hairpin[five_prime[n][0], x,  y, three_prime[n][0]] = float(energy_values[m])                     
		                              m=m+1
                                n=n+1
                        else:
				for x in nucleotides:
                     		   	for y in nucleotides:
						before_hairpin[five_prime[n][0], x, y, three_prime[n][0]] = float(energy_values[m])
              		                        before_hairpin[five_prime[n+1][0], x, y, three_prime[n+1][0]] = float(energy_values[m+4])
               		        		m=m+1
					m=m+4
				n=n+2
			m=0
hp.close()

#bulge
# as of now, it is depending on nucleotides just like others, but might need to change to size only --> now changed
bulge_initiation = [3.8, 2.8, 3.2, 3.6, 4.0, 4.4, 4.6, 4.7, 4.8, 4.9, 5.0, 5.1, 5.2, 5.3, 5.4, 5.4, 5.5, 5.5, 5.6, 5.7, 5.7, 5.8, 5.8, 5.8, 5.9, 5.9, 6.0, 6.0, 6.0, 6.1]


with open ('Energy.py', 'w') as f:
	f.write('STACKING_ENERGY = ' + str(STACKING_ENERGY) + '\n')
	f.write('BULGE_INITIATION = ' + str(bulge_initiation) + '\n')
	f.write('INTERNAL_ENERGY = ' + str(internal) + '\n')
	f.write('OTHER_INTERNAL_ENERGY = ' + str(other_internal) + '\n')
	f.write('HAIRPIN_INITIATION = ' + str(hairpin) + '\n')
	f.write('HAIRPIN_TERMINAL_MISMATCH = ' + str(before_hairpin) + '\n')
f.close()


#testing purpose
#print internal_1x1['G', 'G', 'A', 'U', 'G', 'C']
#print internal_2x1['G', 'G', 'G', 'C', 'U', 'G', 'U']
#print internal_2x2['C', 'A', 'A', 'C', 'G', 'C', 'G', 'G']
#print before_hairpin['G', 'G', 'C', 'U']

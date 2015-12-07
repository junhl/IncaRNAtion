from pprint import pprint
import collections

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
    dic_forward = collections.OrderedDict()
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

if __name__ == '__main__':
    ss = '(((...(.(....).)..((....))..)))...'
    d = decompose(ss)
    pprint(d), #d.items()[0]
    flip = d.items()
    flip.reverse()
    print ""
    back_dict = collections.OrderedDict(flip)
    pprint(back_dict)

    print ss, len(ss)

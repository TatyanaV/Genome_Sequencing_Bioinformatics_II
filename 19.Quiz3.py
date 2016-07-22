#https://github.com/4DD8A19D69F5324F9D49D17EF78BBBCC/BioInfom_atics_Algorit_hms/blob/6fb1b85259368775e970fb5eb7b86de202bd327d/Week2/week2.py
from collections import Counter,deque
def str_to_array(strs):
    return [int(x) for x in strs.split()]

def cyclospectrum(pro,dict_=lambda x:mass[x],linear=False):
    ret =[0]
    for i in range(len(pro)):
        t = 0
        for j in range(i+1):
            t+=dict_(pro[j])
        ret.append(t)
        tmp = pro+pro[:i]
        if linear or i == len(pro)-1:
            tmp=pro
        for j in range(i+1,len(tmp)):
            t+=dict_(tmp[j])
            t-=dict_(tmp[j-(i+1)])
            ret.append(t)
    return sorted(ret)
peps = 'TMLA TALM IAMT ALTM MLAT TAIM'.split()
mass = {pro:int(mass) for pro,mass in (line.split() for line in open("integer_mass_table.txt"))}
spec1 = str_to_array('0 71 101 113 131 184 202 214 232 285 303 315 345 416')
print(list(filter(lambda x: cyclospectrum(x)==spec1,peps)))

def cycloscoring(seq,spec,dict_=lambda x:mass[x],linear=False):
    spec_real = cyclospectrum(seq,dict_,linear=linear)
    c1 = Counter(spec)
    c2 = Counter(spec_real)
    ret = 0
    for key in c1:
        if key in c2:
            ret+= min(c1[key],c2[key])
    return ret

spec4 = str_to_array('0 71 99 101 103 128 129 199 200 204 227 230 231 298 303 328 330 332 333')
peps2 = 'AVQ CTV TCE TCQ VAQ ETC'.split()
print(list(map(lambda x:cycloscoring(x,spec4,linear=True),peps2)))


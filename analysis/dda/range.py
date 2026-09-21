import numpy as np
PROTON=1.00727646
p='/home/user/csf/results/quandenser/consensus_spectra/Quandenser.consensus.part1.ms2'
z=mh=None; mz=[]; rows=[]
for line in open(p):
    c=line[0]
    if c=='S':
        if z is not None and mz: rows.append((z,mh,min(mz),max(mz),len(mz)))
        z=mh=None; mz=[]
    elif c=='Z':
        f=line.split(); z=int(f[1]); mh=float(f[2])
    elif c in 'HID': continue
    else: mz.append(float(line.split(None,1)[0]))
if z is not None and mz: rows.append((z,mh,min(mz),max(mz),len(mz)))
a=np.array(rows)
zz,mh,lo,hi,n = a[:,0],a[:,1],a[:,2],a[:,3],a[:,4]
S=mh+PROTON
print('%-4s %6s %9s %9s %9s %9s %9s'%('z','n','medS*','med maxpk','maxpk/S*','med minpk','frac maxpk>S*/2'))
for c in (1,2,3,4):
    m=zz==c
    if m.sum()<50: continue
    print('%-4d %6d %9.1f %9.1f %9.2f %9.1f %9.2f'%(
        c,m.sum(),np.median(S[m]),np.median(hi[m]),np.median(hi[m]/S[m]),
        np.median(lo[m]),np.mean(hi[m]>S[m]/2)))
print()
print('A complementary pair needs one peak above S*/2. Fraction of spectra')
print('whose highest peak is below S*/2 (no pair is even possible):')
for c in (1,2,3,4):
    m=zz==c
    if m.sum()<50: continue
    print('  z=%d  %5.1f%%'%(c,100*np.mean(hi[m]<=S[m]/2)))

import numpy as np
from pyteomics import mzml
PROTON=1.00727646
# MGF: index -> (charge, scans, mz, inten)
mg={}
idx=z=sc=None; mz=[]; it=[]
for line in open('/home/user/csf/calib/20150513_20_hTau_CSF_Frac9of12_rep1.mgf'):
    line=line.strip()
    if line=='BEGIN IONS': idx=z=sc=None; mz=[]; it=[]
    elif line=='END IONS':
        if idx is not None: mg[idx]=(z,sc,np.array(mz),np.array(it))
    elif line.startswith('TITLE='): idx=int(line.split('index=')[1])
    elif line.startswith('CHARGE='): z=int(line[7:].rstrip('+'))
    elif line.startswith('SCANS='): sc=int(line[6:])
    elif line and line[0].isdigit():
        a=line.split(); mz.append(float(a[0])); it.append(float(a[1]))
# mzML by scan number
mm={}
for s in mzml.read('/home/user/csf/mzml/20150513_20_hTau_CSF_Frac9of12_rep1.mzML'):
    if s.get('ms level')!=2: continue
    sn=int(s['id'].split('scan=')[1].split()[0])
    pre=s['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
    mm[sn]=(pre.get('charge state'), float(pre['selected ion m/z']),
            np.asarray(s['m/z array'],float))
print('MGF spectra %d   mzML MS2 %d'%(len(mg),len(mm)))
nmg=[];nmm=[];dz=0;dprec=[]
for i,(z,sc,m,t) in mg.items():
    if sc not in mm: continue
    z2,pmz,m2=mm[sc]
    nmg.append(len(m)); nmm.append(len(m2))
    if z2!=z: dz+=1
print('matched by scan number: %d'%len(nmg))
print('peaks per spectrum   MGF mean %.1f   mzML mean %.1f'%(np.mean(nmg),np.mean(nmm)))
print('charge disagreements: %d (%.1f%%)'%(dz,100*dz/len(nmg)))
# does the MGF contain peaks above precursor m/z (sign of no deconvolution
# is normal; deconvoluted spectra would have peaks up to the peptide mass)
above=[];frac=[]
for i,(z,sc,m,t) in list(mg.items())[:4000]:
    if sc not in mm: continue
    z2,pmz,m2=mm[sc]
    above.append(np.mean(m>pmz)); frac.append(np.mean(m2>pmz))
print('fraction of peaks above precursor m/z:  MGF %.3f   mzML %.3f'%(np.mean(above),np.mean(frac)))
# direct overlap of peak lists
ov=[]
for i,(z,sc,m,t) in list(mg.items())[:2000]:
    if sc not in mm: continue
    _,_,m2=mm[sc]; m2=np.sort(m2)
    hit=0
    for v in m:
        p=np.searchsorted(m2,v)
        if (p<len(m2) and abs(m2[p]-v)<=0.01) or (p>0 and abs(m2[p-1]-v)<=0.01): hit+=1
    ov.append(hit/max(len(m),1))
print('fraction of MGF peaks found in mzML within 0.01 Da: %.3f'%np.mean(ov))

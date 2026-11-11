"""Exact enumeration of two-periodic GF(2^k) code lengths.

Hypotheses: two information-error types (one substitution or one adjacent
swap), two untouched final check symbols, primitive alpha, nonzero beta !=1.
Uses the proved short-candidate restriction for m=2^(k-1)-c.
No external Python dependencies.
"""

from math import gcd
from time import perf_counter

def poly_mul(a,b,k,p):
 q=1<<k;v=0
 while b:
  if b&1:v^=a
  b>>=1;a<<=1
  if a&q:a^=p
 return v

def fastpow(a,n,k,p):
 y=1
 while n:
  if n&1:y=poly_mul(y,a,k,p)
  n>>=1
  if n:a=poly_mul(a,a,k,p)
 return y

def primes(n):
 f=[];p=2
 while p*p<=n:
  if n%p==0:
   f.append(p)
   while n%p==0:n//=p
  p+=1
 if n>1:f.append(n)
 return f

def primitive_polynomial(k):
 Q=(1<<k)-1
 fs=primes(Q)
 for poly in range((1<<k)+1,(1<<(k+1)),2):
  if fastpow(2,Q,k,poly)!=1:continue
  if all(fastpow(2,Q//f,k,poly)!=1 for f in fs):return poly
 raise ValueError('no primitive polynomial')

def candidates(Q,m):
 # exact list of (b,D) satisfying syndrome direction distinct and invertible systematic checks
 bs={}
 for b in range(2,Q):
  vv=sorted(set((j-(b if j%2==0 else 0))%Q for j in range(1,m+1)))
  if len(vv)!=m:continue
  if m%2 and b==Q-1:continue
  ds=[]
  for i in range(m):
   gap=(vv[(i+1)%m]-vv[i])%Q-1
   if gap>=m-1:
    start=(vv[i]+1)%Q
    for j in range(gap-m+2):ds.append((start+j-1)%Q)
  if ds: bs[b]=set(ds)
 return bs


def candidates_near(Q,m,c):
    s=m if m%2 else m-1
    loD=s; hiD=Q-m+1
    assert loD<=hiD
    out={}
    for B in range(-2*c-2,2*c+3):
        b=B%Q
        if b in (0,1):continue
        if m%2 and b==Q-1:continue
        V={(j-(b if j%2==0 else 0))%Q for j in range(1,m+1)}
        if len(V)!=m:continue
        Ds=[]
        for D in range(loD,hiD+1):
            if all((D+i)%Q not in V for i in range(1,m)):
                Ds.append(D)
        if Ds:out[b]=set(Ds)
    return out

def primitive_representatives(k,poly):
 Q=(1<<k)-1
 xp=[1]*Q
 for j in range(1,Q):xp[j]=poly_mul(xp[j-1],2,k,poly)
 lg={x:i for i,x in enumerate(xp)}
 assert len(lg)==Q
 # Frobenius orbits: r->2r modQ
 visited=set(); reps=[]
 for r in range(1,Q):
  if gcd(r,Q)>1 or r in visited:continue
  orbit=[(r*(1<<j))%Q for j in range(k)]
  assert len(set(orbit))==k
  visited.update(orbit)
  reps.append(r)
 return xp,lg,reps

def determine(k,max_c=None):
 t=perf_counter();q=1<<k;Q=q-1;h=q//2
 pp=primitive_polynomial(k)
 xp,lg,reps=primitive_representatives(k,pp)
 if max_c is None:max_c=min(h-2,100)
 for c in range(max_c+1):
  m=h-c
  if m<2:break
  bs=candidates_near(Q,m,c)
  for r in reps:
   inv=pow(r,-1,Q)
   v=lg[1^xp[r]]
   for b,Dset in bs.items():
    D=(inv*(v-lg[1^xp[(r*b)%Q]]))%Q
    if D in Dset:
     alpha=xp[r];beta=xp[(r*b)%Q]
     return {'k':k,'q':q,'best':m,'deficit':c,'alpha':alpha,'beta':beta,'b':b,'D':D,'r':r,'primitive_polynomial':hex(pp),'orbit_reps':len(reps),'formal_pairs':sum(len(ds) for ds in bs.values()),'duration':round(perf_counter()-t,3)}
 return {'k':k,'q':q,'failure':f'no result at c<={max_c}'}

if __name__ == '__main__':
    import sys
    maximum_k = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    expected = {2: 2, 3: 3, 4: 5, 5: 14, 6: 29, 7: 61,
                8: 122, 9: 250, 10: 508, 11: 1019, 12: 2033,
                13: 4090, 14: 8186, 15: 16373, 16: 32757,
                17: 65528, 18: 131066}
    for k in range(2, maximum_k + 1):
        record = determine(k)
        if k in expected:
            assert record['best'] == expected[k], record
        print('k={k:2d}, q={q:7d}, m_max={best:7d}, deficit={deficit:3d}, '
              'primitive_poly={primitive_polynomial}, alpha={alpha}, beta={beta}, '
              'b={b}, D={D}'.format(**record), flush=True)


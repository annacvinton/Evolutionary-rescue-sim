#!/usr/bin/env python3
# v3 landscapes: same spectral generator and seed scheme as v2, extended to 5 draws.
# Draws r0-r2 are IDENTICAL to landscapes_v2 (same seeds); r3, r4 are new.
import numpy as np, os

def field(ac, sd, n=100, seed=0):
    rng = np.random.default_rng(seed)
    kx = np.fft.fftfreq(n)[:, None]; ky = np.fft.fftfreq(n)[None, :]
    k = np.sqrt(kx**2 + ky**2); k[0, 0] = 1e-9
    amp = k**(-ac/2.0); amp[0, 0] = 0
    ph = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    f = np.real(np.fft.ifft2(amp*ph))
    f -= f.mean(); f *= sd/f.std() if sd > 0 else 0
    return f if sd > 0 else np.zeros((n, n))

os.makedirs("landscapes_v3", exist_ok=True)
count = 0
for sd in [0, 1, 2]:
    for ac in [0, 2, 4]:
        if sd == 0 and ac != 0: continue
        for r in range(5):
            f = field(ac, sd, seed=1000*ac + 10*sd + r)
            np.savetxt(f"landscapes_v3/L_ac{ac}_sd{sd}_r{r}.txt", f.ravel(), fmt="%.15f")
            count += 1
print(f"{count} landscapes written (7 conditions x 5 draws, minus flat duplicates)")
# verify r0-r2 match v2 exactly where v2 exists
import glob, filecmp
mism = 0
for f in glob.glob("landscapes_v3/*_r[012].txt"):
    v2 = f.replace("landscapes_v3", "landscapes_v2")
    if os.path.exists(v2) and not filecmp.cmp(f, v2, shallow=False):
        mism += 1; print("MISMATCH:", f)
print("r0-r2 identical to v2" if mism == 0 else f"{mism} mismatches vs v2")

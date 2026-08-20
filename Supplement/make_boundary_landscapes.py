#!/usr/bin/env python3
# Landscapes for the supplementary viability map: patch SD 0.5-3 x AC 0-4,
# 3 independent draws each, same spectral generator as the main design.
import numpy as np, os

def field(ac, sd, n=100, seed=0):
    rng = np.random.default_rng(seed)
    kx = np.fft.fftfreq(n)[:, None]; ky = np.fft.fftfreq(n)[None, :]
    k = np.sqrt(kx**2 + ky**2); k[0, 0] = 1e-9
    amp = k**(-ac/2.0); amp[0, 0] = 0
    ph = rng.normal(size=(n, n)) + 1j*rng.normal(size=(n, n))
    f = np.real(np.fft.ifft2(amp*ph))
    f -= f.mean(); f *= sd/f.std()
    return f

os.makedirs("landscapes_boundary", exist_ok=True)
SDS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
ACS = [0, 1, 2, 3, 4]
print(f"{'file':28s}{'SD':>5s}{'lag-1':>8s}")
for sd in SDS:
    for ac in ACS:
        for r in range(3):
            f = field(ac, sd, seed=int(10000*sd + 100*ac + r))
            name = f"landscapes_boundary/B_ac{ac}_sd{sd}_r{r}.txt"
            np.savetxt(name, f.ravel(), fmt="%.15f")
            if r == 0:
                l1 = np.corrcoef(f[:, :-1].ravel(), f[:, 1:].ravel())[0, 1]
                print(f"{name.split('/')[1]:28s}{f.std():5.2f}{l1:8.3f}")
print(f"\n{len(SDS)*len(ACS)*3} landscapes written")

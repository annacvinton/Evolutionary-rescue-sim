"""Spectral landscape generator, matching the Mathematica Random_Fields_Generator.

Power spectral density decays as k^-AC; white noise added per frequency; real part
of the inverse FFT; normalised to mean 0 and scaled to the target SD.

Writes several independent draws per (AC, SD) so that landscape identity is not
confounded with treatment level, which it was in the original runs.
"""
import numpy as np

N = 100
AC_LEVELS = [0, 2, 4]
SD_LEVELS = [0, 1, 3]          # from the pilot: viable window
N_DRAWS = 3                        # independent realisations per (AC, SD)


def make_field(ac, sd, n=N, seed=0):
    rng = np.random.default_rng(seed)
    kx = np.fft.fftfreq(n)[:, None]
    ky = np.fft.fftfreq(n)[None, :]
    k = np.sqrt(kx ** 2 + ky ** 2)
    k[0, 0] = 1e-9
    amp = k ** (-ac / 2.0)         # power ~ k^-ac  =>  amplitude ~ k^(-ac/2)
    amp[0, 0] = 0                  # zero mean
    ph = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    if sd == 0:
        return np.zeros((n, n))          # SD 0 is the pure gradient: exactly flat
    f = np.real(np.fft.ifft2(amp * ph))
    f -= f.mean()
    f *= sd / f.std()
    return f


if __name__ == "__main__":
    rows = []
    for ac in AC_LEVELS:
        for sd in SD_LEVELS:
            for r in range(N_DRAWS):
                f = make_field(ac, sd, seed=1000 * ac + 10 * sd + r)
                name = f"L_ac{ac}_sd{sd}_r{r}.txt"
                import os; os.makedirs("landscapes",exist_ok=True); np.savetxt("landscapes/"+name, f.ravel(), fmt="%.15f")
                l1 = (np.corrcoef(f[:, :-1].ravel(), f[:, 1:].ravel())[0, 1]
                      if sd > 0 else 0.0)
                rows.append((name, ac, sd, r, f.std(), l1))
    print(f"{'file':28s} {'AC':>3s} {'SD':>4s} {'draw':>5s} {'realised SD':>12s} {'lag-1 corr':>11s}")
    for n_, a, s, r, sdv, l1 in rows:
        print(f"{n_:28s} {a:>3d} {s:>4d} {r:>5d} {sdv:>12.3f} {l1:>11.3f}")
    print(f"\n{len(rows)} landscapes written")

import numpy as np, sympy as sp, matplotlib.pyplot as plt
from math import gcd, pi
import cmath as C


def _pow2_decompose(a, e):
    mod, tgt = 1 << e, a % (1 << e)
    for eps in (1, -1):
        order = 1 << (e - 2)
        for b in range(order):
            if (eps * pow(5, b, mod)) % mod == tgt:
                return eps, b


def _component(p, e, index):
    mod = p**e
    if p == 2 and e > 2:
        eps_n, a = _pow2_decompose(index % mod, e)
        denom = 1 << (e - 2)
        return (
            lambda m: 0j
            if m % 2 == 0
            else C.exp(
                2j
                * pi
                * (
                    ((1 - eps_n) * (1 - _pow2_decompose(m % mod, e)[0])) / 8
                    + (a * _pow2_decompose(m % mod, e)[1]) / denom
                )
            )
        )
    g = sp.primitive_root(mod)
    phi = (p - 1) * p ** (e - 1)
    val, table = 1 % mod, {1 % mod: 0}
    for b in range(1, phi):
        val = (val * g) % mod
        table[val] = b
    a = table[index % mod]
    return lambda m: 0j if m % p == 0 else C.exp(2j * pi * (a * table[m % mod]) / phi)


def make_dirichlet_character(q, index):
    comps = [_component(p, e, index % (p**e)) for p, e in sp.factorint(q).items()]
    return (
        lambda m: 0j
        if gcd(m, q) != 1
        else (lambda z: z)(np.prod([c(m) for c in comps], dtype=complex))
    )


def chi_interpolant_factory(chi, q):
    n = np.arange(q)
    chi_vals = np.array([chi(int(k)) for k in n], dtype=complex)
    pref = chi((-1) % q) * (np.conj(chi_vals) * np.exp(2j * np.pi * n / q)).sum() / q
    ks = n[:, None]
    return lambda x: pref * (
        np.conj(chi_vals)[:, None] * np.exp(2j * np.pi * ks * np.asarray(x)[None, :] / q)
    ).sum(0)


def gauss_sum(chi, q):
    return sum(chi(n) * C.exp(2j * pi * n / q) for n in range(q))


def plot_character(q, index, samples=120_000, cycles=1.0, linewidth=0.6):
    chi = make_dirichlet_character(q, index)
    f = chi_interpolant_factory(chi, q)

    X = np.linspace(0.0, q * cycles, int(samples), endpoint=False)
    Y = f(X)

    _, ax = plt.subplots(figsize=(7, 7), dpi=250)
    ax.plot(Y.real, Y.imag, lw=linewidth, alpha=0.95)

    P = np.array([chi(n) for n in range(q)])
    M = np.abs(P) > 0

    ax.scatter(P.real[M], P.imag[M], s=2, alpha=0.25)
    t = np.linspace(0, 2 * np.pi, 800)
    ax.plot(np.cos(t), np.sin(t), lw=0.3, alpha=0.25)
    ax.set_aspect("equal", "box")
    ax.axis("off")


if __name__ == "__main__":
    plot_character(q=22, index=5, samples=120_000, cycles=1.0, linewidth=0.5)
    plt.show()

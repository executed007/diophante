#!/usr/bin/env python3
"""diophante — équations diophantiennes linéaires et équation de Pell (stdlib seule, Python ≥ 3.8).

  python diophante.py "6x + 10y + 14z = 22"      solution générale (n variables)
  python diophante.py "3x + 5y = 100" -n         + solutions entières ≥ 0
  python diophante.py "x^2 - 61y^2 = 1" -k 3     Pell : x² − D·y² = ±1
"""
import argparse
import re
from itertools import islice
from math import gcd, isqrt


# ── Linéaire : a·x = c ──────────────────────────────────────────────────────
def solve_linear(a, c):
    """Résout a₁x₁ + … + aₙxₙ = c dans Z.

    Euclide « par colonnes » : on réduit (a₁…aₙ) à (g, 0, …, 0) par opérations
    unimodulaires, en les mémorisant dans `cols`. La colonne du pgcd donne une
    solution particulière, les colonnes nulles engendrent toutes les solutions
    du système homogène. Renvoie (particulière, base_du_noyau) ou None.
    """
    n = len(a)
    vals = list(a)
    cols = [[int(i == j) for j in range(n)] for i in range(n)]
    while True:
        nz = [i for i, v in enumerate(vals) if v]
        if len(nz) <= 1:
            break
        p = min(nz, key=lambda i: abs(vals[i]))
        for j in nz:
            if j != p:
                q = vals[j] // vals[p]
                vals[j] -= q * vals[p]
                cols[j] = [x - q * y for x, y in zip(cols[j], cols[p])]
    if not nz:
        return ([0] * n, cols) if c == 0 else None
    p = nz[0]
    if c % vals[p]:
        return None
    k = c // vals[p]
    return [k * x for x in cols[p]], [cols[i] for i in range(n) if i != p]


def nonneg(a, c):
    """Génère toutes les solutions x ≥ 0 de a·x = c (aᵢ > 0, c ≥ 0)."""
    if len(a) == 1:
        if c % a[0] == 0:
            yield (c // a[0],)
        return
    if len(a) == 2:                       # forme fermée : un pas de b/g à chaque solution
        x, y = a
        g = gcd(x, y)
        if c % g:
            return
        x, y, c = x // g, y // g, c // g
        u = c * pow(x, -1, y) % y
        while x * u <= c:
            yield (u, (c - x * u) // y)
            u += y
        return
    for u in range(c // a[0] + 1):
        for rest in nonneg(a[1:], c - a[0] * u):
            yield (u, *rest)


# ── Pell : x² − D·y² = ±1 ───────────────────────────────────────────────────
def pell(D, N=1):
    """Plus petite solution positive de x² − D·y² = N (N = ±1), ou None.

    Fraction continue de √D : la première réduite donnant ±1 tranche tout.
    Si c'est −1, la solution de +1 s'obtient en l'élevant au carré ;
    si c'est +1 d'abord, l'équation en −1 n'a pas de solution.
    """
    a0 = isqrt(D)
    if a0 * a0 == D:
        return None
    m, d, a = 0, 1, a0
    h0, h, k0, k = 1, a0, 0, 1
    while (v := h * h - D * k * k) not in (1, -1):
        m = d * a - m
        d = (D - m * m) // d
        a = (a0 + m) // d
        h0, h = h, a * h + h0
        k0, k = k, a * k + k0
    if v == N:
        return h, k
    return (h * h + D * k * k, 2 * h * k) if N == 1 else None


def pell_solutions(D, N=1):
    """Génère toutes les solutions positives, dans l'ordre croissant."""
    first = pell(D, N)
    if first:
        X, Y = pell(D, 1)
        x, y = first
        while True:
            yield x, y
            x, y = x * X + D * y * Y, x * Y + y * X


# ── Entrée / sortie ─────────────────────────────────────────────────────────
PELL = re.compile(r"x\^2([+-])(\d*)\*?y\^2=([+-]?1)")
LINEAR = re.compile(r"(?:[+-]?\d*\*?[a-z])+")
TERM = re.compile(r"([+-]?)(\d*)\*?([a-z])")


def parse_linear(s):
    lhs, rhs = s.split("=")
    if not LINEAR.fullmatch(lhs):
        raise ValueError("équation non reconnue")
    coeffs = {}
    for sign, num, var in TERM.findall(lhs):
        coeffs[var] = coeffs.get(var, 0) + (-1 if sign == "-" else 1) * int(num or 1)
    return list(coeffs), list(coeffs.values()), int(rhs)


def expr(p, coefs, params):
    s = str(p) if p else ""
    for c, t in zip(coefs, params):
        if c:
            s += f" {'-' if c < 0 else '+'} {'' if abs(c) == 1 else abs(c)}{t}"
    s = s.strip()
    if s[:2] in ("+ ", "- "):
        s = ("-" if s[0] == "-" else "") + s[2:]
    return s or "0"


def run(args):
    s = args.equation.replace(" ", "").replace("²", "^2").replace("**", "^").lower()

    try:
        if m := PELL.fullmatch(s):
            D, N = int(m[2] or 1) * (1 if m[1] == "-" else -1), int(m[3])
            sols = list(islice(pell_solutions(D, N), args.k)) if D > 0 else []
            if not sols:
                return print("Aucune solution non triviale.")
            for x, y in sols:
                print(f"(x, y) = ({x}, {y})")
            return

        names, a, c = parse_linear(s)
    except (ValueError, TypeError):
        return print("Équation non reconnue. Exemples : 6x+10y+14z=22  |  x^2-61y^2=1")

    res = solve_linear(a, c)
    if res is None:
        return print("Aucune solution entière.")

    p, basis = res
    params = ["t"] if len(basis) == 1 else [f"t{i + 1}" for i in range(len(basis))]
    print("Solution générale (paramètres entiers" + (f" {', '.join(params)}" if basis else "") + ") :")
    for i, name in enumerate(names):
        print(f"  {name} = {expr(p[i], [v[i] for v in basis], params)}")

    if args.nonneg:
        if min(a) <= 0 or c < 0:
            return print("\n--nonneg demande des coefficients > 0 et un second membre ≥ 0.")
        sols = list(islice(nonneg(a, c), args.limit + 1))
        print(f"\nSolutions ≥ 0 ({', '.join(names)}) :" if sols else "\nAucune solution ≥ 0.")
        for sol in sols[:args.limit]:
            print(" ", sol)
        if len(sols) > args.limit:
            print("  … (augmente --limit pour en voir plus)")


def main():
    ap = argparse.ArgumentParser(description="Résout des équations diophantiennes.")
    ap.add_argument("equation", nargs="?", help='ex. "6x+10y+14z=22" ou "x^2-61y^2=1" (sinon mode interactif)')
    ap.add_argument("-n", "--nonneg", action="store_true", help="lister les solutions ≥ 0")
    ap.add_argument("-l", "--limit", type=int, default=20, help="max de solutions listées (défaut 20)")
    ap.add_argument("-k", type=int, default=5, help="nombre de solutions de Pell (défaut 5)")
    args = ap.parse_args()

    if args.equation:
        return run(args)

    print("Mode interactif — ligne vide pour quitter. Ajoute « -n » en fin de ligne pour lister les solutions ≥ 0.")
    try:
        while eq := input("\nÉquation > ").strip():
            args.nonneg = eq.endswith(" -n")
            args.equation = eq[:-3] if args.nonneg else eq
            run(args)
    except (EOFError, KeyboardInterrupt):
        print()


if __name__ == "__main__":
    main()
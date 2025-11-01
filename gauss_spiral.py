import sys
import matplotlib.pyplot as plt


def prime_int(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    if any(n % p == 0 for p in small):
        return n in small
    d, s = n - 1, 0
    while d & 1 == 0:
        d >>= 1
        s += 1
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        x = pow(a % n, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def gauss_prime(a, b):
    if not (a or b) or (a, b) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        return False
    if a == 0 or b == 0:
        x = abs(a or b)
        return x % 4 == 3 and prime_int(x)
    return prime_int(a * a + b * b)


DIRS = {"E": (1, 0), "N": (0, 1), "W": (-1, 0), "S": (0, -1)}
ORDER = ["E", "N", "W", "S"]
turn_left = lambda d: ORDER[(ORDER.index(d) + 1) % 4]


def walk_cycle(a, b, d):
    prime_hits, seen, start = [], set(), (a, b, d)
    if gauss_prime(a, b):
        prime_hits.append((a, b))
        d = turn_left(d)
    while (a, b, d) not in seen:
        seen.add((a, b, d))
        dx, dy = DIRS[d]
        while True:
            a += dx
            b += dy
            if gauss_prime(a, b):
                prime_hits.append((a, b))
                d = turn_left(d)
                break
        if (a, b, d) == start:
            break
    return prime_hits


def plot_path(points, title):
    xs, ys = zip(*points)
    plt.figure(figsize=(8, 8))
    plt.plot(xs, ys, linewidth=0.5, color="black")
    plt.axis("equal")
    plt.axis("off")
    plt.title(title)
    plt.show()


def main():
    a0, b0, d0 = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3].upper()
    pts = walk_cycle(a0, b0, d0)
    plot_path(pts, f"{a0}+{b0}i, dir={d0} ({len(pts) - 1} turns)")


if __name__ == "__main__":
    main()

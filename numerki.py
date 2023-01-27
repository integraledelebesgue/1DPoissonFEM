def apply_f(f: list, a: list[float], xi: float) -> list[float]:
    return sum([fi(xi) * ai for fi, ai in zip(f, a)])


def q_form(f: list, a: list[float], sigma: list[float], x: list[float], y: list[float]) -> float:  # f - tablica funkcji anonimowych
    return sum([(apply_f(f, a, xi) - yi)**2 / sigmai**2 for xi, yi, sigmai in zip(x, y, sigma)])

# przykład:
def main():
    a = [1, -2, 3]
    f = [
        lambda x: x**1,
        lambda x: x**2,
        lambda x: x**3
    ]
    sigma = [10, 10, 100]
    x = [2, 4, 6, 8, -10, 12]
    y = [2, 3, 5, 7, 11, 13]

    print(q_form(f, a, sigma, x, y))

main()

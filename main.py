from Elliptic import Eleptic_curves

if __name__ == '__main__':
    print("E(a, b)(F_p)")
    key = [int(input("Введите значение a: ")), int(input("Введите значение b: ")), int(input("Введите значение p: "))]
    Elliptic = Eleptic_curves(key[0], key[1], key[2])
    Elliptic.find_the_generator()
    Elliptic.finding_points()
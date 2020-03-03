import numpy as np
import math as m

populatie_binar = []
populatie_decimal = []
f = []
probabilitati = []
intervale_selectie = []
cromozomi_selectati = []
cromozomi_selectati_incrucisare = []

with open("citire.txt", "r") as fin:
    linie = [int(x) for x in fin.readline().split()]
    dim_pop = linie[0]
    precizie = linie[1]
    nr_etape = linie[2]
    a_interv = linie[3]
    b_interv = linie[4]
    linie = [float(x) for x in fin.readline().split()]
    cross_prob = linie[0]
    mut_prob = linie[1]


def lung_cromozom():
    nr = (b_interv - a_interv) * (10 ** precizie)
    min = m.log(nr, 2)
    return int(min) + 1


def binary_to_interval(x, l_crom):
    nr = int(x, 2)
    nr_interv = (nr / (2 ** l_crom - 1)) * (b_interv - a_interv) + a_interv
    return nr_interv


def functie(x):
    return -x ** 2 + x + 2


def rand_gen(len_cromozom):   #generez aleator populatia initiala
    for i in range(dim_pop):
        bit_cromozom = np.random.randint(2, size=len_cromozom)
        populatie_binar.append(list(bit_cromozom))
    for elem in populatie_binar:
        stringlist = "".join(str(x) for x in elem)
        populatie_decimal.append(binary_to_interval(stringlist, len_cromozom))


def val_functie():
    s = 0
    for elem in populatie_decimal:
        valoare_fct = functie(elem)
        f.append(valoare_fct)
        s += valoare_fct     #performanta totala a populatiei
    return s


def calc_probab(s):
    for i in range(dim_pop):
        probabilitati.append(f[i] / s)   #asociem fiecarui individ probabilitatea de a fi selectat in functie de fitness


def det_interv_selec():
    intervale_selectie.append(probabilitati[0])
    for i in range(1, dim_pop):
        intervale_selectie.append(intervale_selectie[i - 1] + probabilitati[i])
    intervale_selectie.insert(0, 0)


def binary_search(x):
    low = 0
    high = len(intervale_selectie) - 1
    while low <= high:
        mid = (high + low) // 2
        if intervale_selectie[mid] == x:
            return mid
        elif x > intervale_selectie[mid]:
            low = mid + 1
        else:
            high = mid - 1
    return low


def selection():   #procesul de selectie
    global populatie_binar
    global populatie_decimal
    global f
    temp = []
    temp1 = []
    temp2 = []
    for i in range(dim_pop):
        nr_generat = np.random.rand()    #generez o variabila uniforma pe [0,1)
        crom = binary_search(nr_generat)  #determin intervalul caruia apartine
        cromozomi_selectati.append(crom)  #selectez cromozomii pentru etapa urmatoare(aleg cromozomul i+1 din interval)
        with open("evolutie.txt", "a") as fout:
            fout.write("u= " + str(nr_generat) + " selectam cromozomul " + str(crom) + "\n")
        temp.append(populatie_binar[crom - 1])    #updatez vectorii pentru cromozomii selectati
        temp1.append(populatie_decimal[crom - 1])
        temp2.append(f[crom - 1])
    populatie_binar = temp.copy()
    populatie_decimal = temp1.copy()
    f = temp2.copy()


def selectie_incrucisare():
    with open("evolutie.txt", "a") as fout:
        fout.write("Probabilitatea de incrucisare " + str(cross_prob) + "\n")
    for i in range(dim_pop):
        nr_generat = np.random.rand()
        with open("evolutie.txt", "a") as fout:
            fout.write(
                str(i + 1) + ": " + str(populatie_binar[i]) + " u= " + str(nr_generat) + "\n")
            if nr_generat <= cross_prob:      #se aleg cromozomii care participa la cross over
                cromozomi_selectati_incrucisare.append(i)
                fout.write("<" + str(cross_prob) + " participa\n")


def update_valori(l, i):
    stringlist = "".join(str(x) for x in populatie_binar[cromozomi_selectati_incrucisare[i]])
    populatie_decimal[cromozomi_selectati_incrucisare[i]] = binary_to_interval(stringlist, l)
    stringlist = "".join(str(x) for x in populatie_binar[cromozomi_selectati_incrucisare[i + 1]])
    populatie_decimal[cromozomi_selectati_incrucisare[i + 1]] = binary_to_interval(stringlist, l)
    f[cromozomi_selectati_incrucisare[i]] = functie(populatie_decimal[cromozomi_selectati_incrucisare[i]])
    f[cromozomi_selectati_incrucisare[i + 1]] = functie(populatie_decimal[cromozomi_selectati_incrucisare[i + 1]])


def update_valori_mutatie(l, i):
    stringlist = "".join(str(x) for x in populatie_binar[i])
    populatie_decimal[i] = binary_to_interval(stringlist, l)
    f[i] = functie(populatie_decimal[i])


def recombinare(l):
    if len(cromozomi_selectati_incrucisare) % 2:   #pentru lungime impara renunt la ultimul
        cromozomi_selectati_incrucisare.pop()
    lungime = len(cromozomi_selectati_incrucisare)
    for i in range(0, lungime, 2):   #fac perechi de cate 2 cromozomi
        if lungime == 0:
            break
        punct_rupere = np.random.randint(l)   #generez aleator punctul de rupere
        with open("evolutie.txt", "a") as fout:
            fout.write(
                "Recombinare dintre cromozomul " + str(
                    cromozomi_selectati_incrucisare[i] + 1) + " cu cromozomul " + str(
                    cromozomi_selectati_incrucisare[i + 1] + 1) + ":\n" + str(
                    populatie_binar[cromozomi_selectati_incrucisare[i]]) + " " + str(
                    populatie_binar[cromozomi_selectati_incrucisare[i + 1]]) + "punct " + str(
                    punct_rupere) + "\n")
        child1 = populatie_binar[cromozomi_selectati_incrucisare[i]][0:punct_rupere] + populatie_binar[
                                                                                           cromozomi_selectati_incrucisare[
                                                                                               i + 1]][
                                                                                       punct_rupere:]
        child2 = populatie_binar[cromozomi_selectati_incrucisare[i + 1]][0:punct_rupere] + populatie_binar[
                                                                                               cromozomi_selectati_incrucisare[
                                                                                                   i]][
                                                                                           punct_rupere:]
        populatie_binar[cromozomi_selectati_incrucisare[i]] = child1   #cu metoda ruletei
        populatie_binar[cromozomi_selectati_incrucisare[i + 1]] = child2
        update_valori(l, i)
        with open("evolutie.txt", "a") as fout:
            fout.write("Rezultat " + str(populatie_binar[cromozomi_selectati_incrucisare[i]]) + " " + str(
                populatie_binar[cromozomi_selectati_incrucisare[i + 1]]) + "\n")


def mutatie(l):
    cromozomi_selectati_mutatie = []
    for i in range(dim_pop):
        prob = np.random.rand()
        if prob <= mut_prob:
            cromozomi_selectati_mutatie.append(i)
    for elem in cromozomi_selectati_mutatie:
        with open("evolutie.txt", "a") as fout:
            fout.write(str(elem + 1) + "\n")
        gena = np.random.randint(l)
        populatie_binar[elem][gena] = 1 - populatie_binar[elem][gena]
        update_valori_mutatie(l, elem)


if __name__ == "__main__":
    maxim_functie = []
    valoare_medie=[]
    l_crom = lung_cromozom()
    for _ in range(nr_etape):
        rand_gen(l_crom)
        s = val_functie()
        with open("evolutie.txt", "a") as fout:
            fout.write("Populatia initiala:\n")
            for i in range(dim_pop):
                fout.write(
                    str(i + 1) + ": " + str(populatie_binar[i]) + " x= " + str(populatie_decimal[i]) + " f=" + str(
                        f[i]) + "\n")
        calc_probab(s)
        maxim_functie.append(max(f))
        valoare_medie.append(s/l_crom)
        with open("evolutie.txt", "a") as fout:
            fout.write("\n Probabilitati selectie:\n")
            for i in range(dim_pop):
                fout.write("cromozom " + str(i + 1) + " probabilitate " + str(probabilitati[i]) + "\n")
        det_interv_selec()
        with open("evolutie.txt", "a") as fout:
            fout.write("\n Intervale probabilitati selectie:\n")
            for i in range(dim_pop + 1):
                fout.write(str(intervale_selectie[i]) + "\n")
        selection()
        with open("evolutie.txt", "a") as fout:
            fout.write("\n Dupa selectie:\n")
            for i in range(dim_pop):
                fout.write(
                    str(i + 1) + ": " + str(populatie_binar[i]) + " x= " + str(populatie_decimal[i]) + " f=" + str(
                        f[i]) + "\n")
        selectie_incrucisare()
        recombinare(l_crom)
        with open("evolutie.txt", "a") as fout:
            fout.write("\n Dupa recombinare:\n")
            for i in range(dim_pop):
                fout.write(
                    str(i + 1) + ": " + str(populatie_binar[i]) + " x= " + str(populatie_decimal[i]) + " f=" + str(
                        f[i]) + "\n")
        with open("evolutie.txt", "a") as fout:
            fout.write("Probabilitate de mutatie pentru fiecare gena " + str(
                mut_prob) + "\n" + "Au fost modificati cromozomii: \n")
        mutatie(l_crom)
        with open("evolutie.txt", "a") as fout:
            fout.write("Dupa mutatie:\n")
            for i in range(dim_pop):
                fout.write(
                    str(i + 1) + ": " + str(populatie_binar[i]) + " x= " + str(populatie_decimal[i]) + " f=" + str(
                        f[i]) + "\n")
        intervale_selectie.clear()
        cromozomi_selectati_incrucisare.clear()
        cromozomi_selectati.clear()
        probabilitati.clear()
        f.clear()
        populatie_decimal.clear()
        populatie_binar.clear()
    print("Evolutia maximului: \n")
    print(maxim_functie)
    print(valoare_medie)

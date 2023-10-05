import matplotlib.pyplot as plt
import numpy as np

def readFile(interest, dt, filename):
    time = []
    t = 0
    
    data  = open(filename)
    next(data)
    if interest == 'energies':
        T = []
        V = []
        E = []
        for line in data:
            values = line.split()
            T.append(float(values[0]))
            V.append(float(values[1]))
            E.append(float(values[2]))
            time.append(t)
            t += dt
        data.close()
        return {"T":T, "V":V, "E":E, "t":time}
    elif interest == 'n_med':
        n_med = []
        for line in data:
            values = line.split()
            n_med.append(float(values[0]))
            time.append(t)
            t += dt
        data.close()
        return {"n":n_med, "t": time}
    elif interest == 'n_med_sphere':
        n_med = []
        for line in data:
            values = line.split()
            n_med.append(float(values[1]))
            time.append(t)
            t += dt
        data.close()
        return {"n":n_med, "t": time}
    

# Exercício 4.1
dt = 1*10**(+5)
energies = readFile('energies', dt, 'data1.txt')

plt.title("Energias cinética, potencial e mecânica em função do tempo")
plt.xlabel("t (anos)")
plt.ylabel(r'Energia $(pc^2 M_{solar} ano^{-2} )$')
plt.plot(energies["t"], energies["T"], label = "T (E. Cinética)")
plt.plot(energies["t"], energies["V"], label = "V (E. Potencial)")
plt.plot(energies["t"], energies["E"], label = "E (E. Mecânica)")
plt.legend()
plt.show()

# Exercício 4.2
dt = 1*(10**5)
n10  = readFile('n_med', dt, '10output.txt' )
n100 = readFile('n_med', dt, '100output.txt')
n500 = readFile('n_med', dt, '500output.txt')
    
plt.plot(n10["t"] , n10["n"] , label  = "Sistema de 10 estrelas")
plt.plot(n100["t"], n100["n"], label = "Sistema de 100 estrelas")
plt.plot(n500["t"], n500["n"], label = "Sistema de 500 estrelas")

plt.title(r'Fração de estrelas a uma distância inferior a $r_t=5$pc do centro do aglomerado em função do tempo')
plt.xlabel("t (anos)")
plt.ylabel("Fração de estrelas")
plt.grid(True)
plt.legend()
plt.show()

# Exercício 4.3
dt = 1*(10**5)
n = readFile('n_med_sphere', dt, 'sphere_output.txt')

plt.plot(n["t"], n["n"], marker='D', markerfacecolor = 'r', markeredgecolor = 'r', markersize=2.5)
plt.title(r'Média da fração de estrelas a uma distância inferior a $r_t=5$pc do centro do aglomerado')
plt.xlabel("t (anos)")
plt.ylabel("média da fração de estrelas")
#plt.yscale('log')
plt.grid(True)
plt.show()
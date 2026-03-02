# 1D Finite-Elemente-Methode (FEM) Löser

Dieses Projekt implementiert einen numerischen Löser basierend auf der **Finite-Elemente-Methode (FEM)** zur Lösung von eindimensionalen Differentialgleichungen (z. B. Poisson-Gleichung) in C++.

## 📖 Beschreibung
Das Programm berechnet die Verschiebung $u$ in einem Stab der Länge $L$, der einer konstanten Quellkraft $f$ ausgesetzt ist. [cite_start]Es nutzt lineare Ansatzfunktionen und eine Gauß-Integration zur Systemassemblierung[cite: 10, 11, 12].

## 🛠 Technische Merkmale
* **Symmetrische Matrix**: Effiziente Speicherung der Steifigkeitsmatrix $K$, bei der nur das untere Dreieck gespeichert wird, um Speicherplatz zu sparen.
* **Numerische Integration**: Verwendung der 2-Punkt Gauß-Quadratur (Gauß-Punkte: $\pm 0.57735$, Gewichte: $1.0$).
* **Gleichungslöser**: Implementierung einer **Cholesky-Zerlegung** zur Lösung des linearen Gleichungssystems $Ku = F$.
* **Randbedingungen**: Anwendung der Penalty-Methode zur Durchsetzung von Dirichlet-Randbedingungen ($u(0)=0, u(L)=0$).

## 🚀 Verwendung

### 1. Kompilierung (C++)
Sie benötigen einen C++ Compiler (wie GCC oder Clang). Kompilieren Sie den Quellcode mit folgendem Befehl:
```bash
g++ -O3 -o fem_loeser main.cpp
./fem_loeser

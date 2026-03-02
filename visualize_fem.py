import pandas as pd
import matplotlib.pyplot as plt

def visualize_fem():
    try:
        # Φόρτωση δεδομένων από το CSV
        data = pd.read_csv('results.csv')
        
        plt.figure(figsize=(10, 6))
        plt.plot(data['x'], data['u'], 'bo-', label='FEM Solution')
        
        # Προσθήκη θεωρητικής λύσης (αν f=10, L=2, u(0)=u(2)=0, τότε u = 5x(2-x))
        x_fine = data['x']
        u_exact = 5 * x_fine * (2 - x_fine)
        plt.plot(x_fine, u_exact, 'r--', label='Exact Solution')
        
        plt.title('1D FEM Results')
        plt.xlabel('Position (x)')
        plt.ylabel('Displacement (u)')
        plt.legend()
        plt.grid(True)
        plt.show()
    except FileNotFoundError:
        print("Σφάλμα: Δεν βρέθηκε το αρχείο results.csv. Τρέξε πρώτα τον C++ κώδικα!")

if __name__ == "__main__":
    visualize_fem()
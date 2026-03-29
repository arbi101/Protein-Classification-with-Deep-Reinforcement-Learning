import random
import time

def sequence_to_hp(sequence):
    """
    Converts an amino acid sequence to HP string
    H: Hydrophobic, P: Polar
    """
    hydrophobic = set('ACFILMVWY')  # Common hydrophobic amino acids
    hp_string = ''
    for aa in sequence:
        if aa in hydrophobic:
            hp_string += 'H'
        else:
            hp_string += 'P'
    return hp_string


def generate_2d_structure(hp_string, iterations=1000000):
    """
    Generates a 2D structure using Hill Climbing with Pivot Moves
    Returns list of positions (x,y) and energy
    """
    if not hp_string:
        return [], 0

    def calculate_energy(pos_list):
        e = 0
        pos_dict = {pos: i for i, pos in enumerate(pos_list)}
        for i in range(len(pos_list)):
            if hp_string[i] == 'H':
                x, y = pos_list[i]
                # Check non-consecutive neighbors
                for dx, dy in [(0,1), (1,0), (0,-1), (-1,0)]:
                    nx, ny = x + dx, y + dy
                    if (nx, ny) in pos_dict:
                        j = pos_dict[(nx, ny)]
                        if abs(i - j) > 1 and hp_string[j] == 'H':
                            e += 1
        return -(e // 2)

    # Initialize with a straight line (valid self-avoiding configuration)
    current_positions = [(i, 0) for i in range(len(hp_string))]
    current_energy = calculate_energy(current_positions)
    
    best_positions = list(current_positions)
    best_energy = current_energy

    for _ in range(iterations):
        if len(hp_string) <= 2:
            break
            
        # Scegli un punto di pivot casuale (tranne l'ultimo nodo)
        pivot_idx = random.randint(1, len(hp_string) - 2)
        # Scegli una rotazione: 90, -90, or 180 gradi
        angle = random.choice([90, -90, 180])
        
        # Applica rotazione ai nodi successivi al pivot
        cx, cy = current_positions[pivot_idx]
        new_positions = list(current_positions)
        
        # Pre-calcola seno e coseno
        if angle == 90:
            cos_a, sin_a = 0, 1
        elif angle == -90:
            cos_a, sin_a = 0, -1
        else: # 180
            cos_a, sin_a = -1, 0
            
        for i in range(pivot_idx + 1, len(hp_string)):
            x, y = current_positions[i]
            # Trasla rispetto al pivot
            tx, ty = x - cx, y - cy
            # Ruota
            rx = tx * cos_a - ty * sin_a
            ry = tx * sin_a + ty * cos_a
            # Ritorna alla posizione originale
            new_positions[i] = (rx + cx, ry + cy)
            
        # Verifica se è auto-evitante (nessuna sovrapposizione)
        if len(set(new_positions)) == len(new_positions):
            # Calcola nuova energia
            new_energy = calculate_energy(new_positions)
            
            # Acceptance condition: accept if new energy is non-increasing (allows moving in plateaus)
            if new_energy <= current_energy:
                current_positions = new_positions
                current_energy = new_energy
                
                # Update absolute best
                if current_energy < best_energy:
                    best_positions = list(current_positions)
                    best_energy = current_energy
                    
    return best_positions, best_energy


if __name__ == "__main__":
    # Define 10 example protein sequences to study
    p1 = "MKVILA"  # 6 AA
    p2 = "MKTIIALSYIFCLVF"  # 15 AA
    p3 = "MVHLTPEEKSAVTALWGK"  # 18 AA
    p4 = "MVHLTPEEKSAVTALWGKVN" # 20 AA
    p5 = "MVHLTPEEKSAVTALWGKVNVDEVG" # 25 AA
    p6 = "MVHLTPEEKSAVTALWGKVNVDEVGGEA" # 28 AA
    p7 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLL"  # 33 AA
    p8 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVY"  # 36 AA
    p9 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWT"  # 39 AA
    p10 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRF"  # 42 AA

    proteins = [
        ("Protein 1 (6 AA)", p1),
        ("Protein 2 (15 AA)", p2),
        ("Protein 3 (18 AA)", p3),
        ("Protein 4 (20 AA)", p4),
        ("Protein 5 (25 AA)", p5),
        ("Protein 6 (28 AA)", p6),
        ("Protein 7 (33 AA)", p7),
        ("Protein 8 (36 AA)", p8),
        ("Protein 9 (39 AA)", p9),
        ("Protein 10 (42 AA)", p10),
    ]

    # Different iterations to study the algorithm's convergence
    iterations_to_test = [1000, 10000, 50000, 100000]

    print("=== TEST HILL CLIMBING ALGORITHM WITH 10 PROTEINS ===\n")

    for name, seq in proteins:
        hp_str = sequence_to_hp(seq)
        
        # Count number of H
        h_count = hp_str.count('H')
        print(f"--- {name} ---")
        print(f"Sequence : {seq}")
        print(f"HP String: {hp_str} (H count: {h_count})\n")

        best_energies = []
        times_taken = []

        for iters in iterations_to_test:
            start_time = time.time()
            positions, energy = generate_2d_structure(hp_str, iterations=iters)
            end_time = time.time()
            
            elapsed_time = end_time - start_time
            best_energies.append(energy)
            times_taken.append(elapsed_time)

            print(f"  Iterations: {iters:<7} | Best Energy: {energy:<3} | Time: {elapsed_time:.4f} sec")
        
        print("\n" + "-"*50 + "\n")

import math
import random
import time

def sequence_to_hp(sequence):
    """
    Converts an amino acid sequence to HP string
    H: Hydrophobic, P: Polar
    """
    hydrophobic = set('ACFILMVWY')
    return ''.join(['H' if aa in hydrophobic else 'P' for aa in sequence])

def generate_2d_structure_mc(hp_string, iterations=100000, temperature=2.0):
    """
    Generates a 2D structure using Basic Monte Carlo with Pivot Moves at Constant Temperature.
    Returns list of positions (x,y) and energy.
    """
    if not hp_string:
        return [], 0

    def calculate_energy(pos_list):
        e = 0
        pos_dict = {pos: i for i, pos in enumerate(pos_list)}
        for i in range(len(pos_list)):
            if hp_string[i] == 'H':
                x, y = pos_list[i]
                for dx, dy in [(0,1), (1,0), (0,-1), (-1,0)]:
                    nx, ny = x + dx, y + dy
                    if (nx, ny) in pos_dict:
                        j = pos_dict[(nx, ny)]
                        if abs(i - j) > 1 and hp_string[j] == 'H':
                            e += 1
        return -(e // 2)

    current_positions = [(i, 0) for i in range(len(hp_string))]
    current_energy = calculate_energy(current_positions)
    
    best_positions = list(current_positions)
    best_energy = current_energy

    for _ in range(iterations):
        if len(hp_string) <= 2:
            break
            
        pivot_idx = random.randint(1, len(hp_string) - 2)
        angle = random.choice([90, -90, 180])
        
        cx, cy = current_positions[pivot_idx]
        new_positions = list(current_positions)
        
        if angle == 90:
            cos_a, sin_a = 0, 1
        elif angle == -90:
            cos_a, sin_a = 0, -1
        else:
            cos_a, sin_a = -1, 0
            
        for k in range(pivot_idx + 1, len(hp_string)):
            x, y = current_positions[k]
            tx, ty = x - cx, y - cy
            rx = tx * cos_a - ty * sin_a
            ry = tx * sin_a + ty * cos_a
            new_positions[k] = (rx + cx, ry + cy)
            
        if len(set(new_positions)) == len(new_positions):
            new_energy = calculate_energy(new_positions)
            
            # Acceptance condition
            if new_energy <= current_energy:
                current_positions = new_positions
                current_energy = new_energy
                
                if current_energy < best_energy:
                    best_positions = list(current_positions)
                    best_energy = current_energy
            else:
                # Basic Monte Carlo: Use constant temperature to accept bad moves
                delta_e = new_energy - current_energy
                probability = math.exp(-delta_e / temperature)
                
                if random.random() < probability:
                    current_positions = new_positions
                    current_energy = new_energy

    return best_positions, best_energy


if __name__ == "__main__":
    p1 = "MKVILA"
    p2 = "MKTIIALSYIFCLVF"
    p3 = "MVHLTPEEKSAVTALWGKVNVDEVG"
    
    proteins = [
        ("Protein 1 (6 AA)", p1),
        ("Protein 2 (15 AA)", p2),
        ("Protein 3 (25 AA)", p3),
    ]

    iterations_to_test = [10000, 50000]

    print("=== TEST BASIC MONTE CARLO ALGORITHM (CONSTANT TEMP) ===\n")

    for name, seq in proteins:
        hp_str = sequence_to_hp(seq)
        print(f"--- {name} ---")
        print(f"Sequence : {seq}")
        print(f"HP String: {hp_str}\n")

        for iters in iterations_to_test:
            start_time = time.time()
            positions, energy = generate_2d_structure_mc(hp_str, iterations=iters, temperature=2.0)
            elapsed = time.time() - start_time
            print(f"  Iterations: {iters:<7} | Best Energy: {energy:<3} | Time: {elapsed:.4f} sec")
        
        print("\n" + "-"*50 + "\n")

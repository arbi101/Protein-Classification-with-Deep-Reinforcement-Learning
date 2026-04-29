import math
import random
import time

def sequence_to_hp(sequence):
    """Convert an amino acid sequence into HP string (H: hydrophobic, P: polar)."""
    hydrophobic = set('ACFILMVWY')
    return ''.join(['H' if aa in hydrophobic else 'P' for aa in sequence])

def generate_2d_structure_remc(hp_string, iterations=300_000, num_replicas=10, t_min=1.0, t_max=30.0, swap_interval=200):
    """
    Optimized Replica Exchange Monte Carlo for real proteins (50-100 AA).
    """
    if not hp_string:
        return [], 0
    if len(hp_string) <= 2:
        return [(i, 0) for i in range(len(hp_string))], 0

    # Function to calculate the energy of the structure
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

    # Logarithmic spacing of temperatures across replicas
    if num_replicas > 1:
        temps = [t_min * (t_max / t_min) ** (i / (num_replicas-1)) for i in range(num_replicas)]
    else:
        temps = [t_min]

    # Initialize positions and energies for all replicas
    replicas_pos = [[(j, 0) for j in range(len(hp_string))] for _ in range(num_replicas)]
    replicas_en = [calculate_energy(r_pos) for r_pos in replicas_pos]
    
    global_best_pos = list(replicas_pos[0])
    global_best_energy = replicas_en[0]

    for step in range(1, iterations + 1):
        # Local Monte Carlo step for each replica
        for i in range(num_replicas):
            pivot_idx = random.randint(1, len(hp_string) - 2)
            angle = random.choice([90, -90])  # stable pivot
            cx, cy = replicas_pos[i][pivot_idx]
            new_positions = list(replicas_pos[i])
            
            cos_a, sin_a = (0, 1) if angle == 90 else (0, -1)
            
            for k in range(pivot_idx + 1, len(hp_string)):
                x, y = current = replicas_pos[i][k]
                tx, ty = x - cx, y - cy
                rx = tx * cos_a - ty * sin_a
                ry = tx * sin_a + ty * cos_a
                new_positions[k] = (rx + cx, ry + cy)
                
            # Accept move if no overlap
            if len(set(new_positions)) == len(new_positions):
                new_energy = calculate_energy(new_positions)
                # Metropolis acceptance criterion
                if new_energy <= replicas_en[i] or random.random() < math.exp(-(new_energy - replicas_en[i]) / temps[i]):
                    replicas_pos[i] = new_positions
                    replicas_en[i] = new_energy
                    # Update global best
                    if new_energy < global_best_energy:
                        global_best_pos = list(new_positions)
                        global_best_energy = new_energy

        # Replica Exchange (swap configurations)
        if step % swap_interval == 0 and num_replicas > 1:
            idx = random.randint(0, num_replicas - 2)
            j_idx = idx + 1
            delta_e = replicas_en[idx] - replicas_en[j_idx]
            delta_beta = (1.0 / temps[idx]) - (1.0 / temps[j_idx])
            exponent = delta_beta * delta_e
            # Metropolis criterion for swap
            if exponent >= 0 or random.random() < math.exp(exponent):
                # Swap positions and energies
                replicas_pos[idx], replicas_pos[j_idx] = replicas_pos[j_idx], replicas_pos[idx]
                replicas_en[idx], replicas_en[j_idx] = replicas_en[j_idx], replicas_en[idx]

    return global_best_pos, global_best_energy

if __name__ == "__main__":

    p = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLL"  # 33 AA

    # Test protein sequence (33 AA) and longer (70+ AA)
    p2 = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    hp_str = sequence_to_hp(p2)
    print("=== TEST OPTIMIZED REMC ===")
    print(f"Sequence : {p2}")
    print(f"HP String: {hp_str}")
    start = time.time()
    _, energy = generate_2d_structure_remc(hp_str, iterations=300_000)
    print(f"Best Energy found: {energy} (Time: {time.time()-start:.2f} s)")





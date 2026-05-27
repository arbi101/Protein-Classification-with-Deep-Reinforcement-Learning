"""
test_sa.py
----------
Implements the Simulated Annealing (SA) algorithm for 2D protein
structure prediction using the HP (Hydrophobic-Polar) lattice model.

SA extends Hill Climbing by occasionally accepting moves that WORSEN
the energy. This is controlled by a temperature parameter T that
decreases over time (cooling schedule), allowing the algorithm to
escape local minima early in the search.
"""

import math
import random
import time


def sequence_to_hp(sequence):
    """
    Converts a raw amino acid sequence into an HP binary string.

    Each residue is mapped to:
      - 'H' (Hydrophobic) if it belongs to {A, C, F, I, L, M, V, W, Y}
      - 'P' (Polar) otherwise

    Parameters
    ----------
    sequence : str  Raw one-letter amino acid sequence
    Returns  : str  HP binary string of the same length
    """
    hydrophobic = set('ACFILMVWY')  # Set of hydrophobic amino acid codes
    hp_string = ''
    for aa in sequence:
        if aa in hydrophobic:
            hp_string += 'H'
        else:
            hp_string += 'P'
    return hp_string


def generate_2d_structure_sa(hp_string, iterations=100000, initial_t=30.0, final_t=0.001,
                             return_trace=False, trace_interval=1000,
                             max_trace_frames=200):
    """
    Generates an approximate minimum-energy 2D protein conformation
    using Simulated Annealing with the full 4-move local move set.

    The algorithm starts from a straight-line conformation and at each
    step proposes a random move. It accepts improving moves always, and
    worsening moves with probability exp(-ΔE / T), where T decreases
    exponentially from initial_t to final_t.

    Parameters
    ----------
    hp_string   : str   Binary HP sequence (e.g. "HPHPH")
    iterations  : int   Total number of move attempts (default: 100,000)
    initial_t   : float Starting temperature T₀ (default: 30.0)
    final_t     : float Final temperature T_f (default: 0.001)

    Returns
    -------
    best_positions : list of (int, int)  Lattice coordinates of best conformation
    best_energy    : int                 Minimum energy found (≤ 0)
    """
    if not hp_string:
        # Empty sequence → nothing to fold
        if return_trace:
            return [], 0, []
        return [], 0

    def calculate_energy(pos_list):
        """
        Computes the HP energy of the given conformation.
        Energy = −(number of non-consecutive H-H topological contacts).
        Each H-H pair that is grid-adjacent but not sequence-adjacent
        contributes −1 to the total energy.
        """
        e = 0
        # Map each (x,y) position to its residue index for O(1) lookups
        pos_dict = {pos: i for i, pos in enumerate(pos_list)}
        for i in range(len(pos_list)):
            if hp_string[i] == 'H':       # Only H residues form HH contacts
                x, y = pos_list[i]
                # Iterate over the 4 orthogonal lattice neighbours
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if (nx, ny) in pos_dict:
                        j = pos_dict[(nx, ny)]
                        # Exclude consecutive residues (backbone bonds, |i-j|=1)
                        if abs(i - j) > 1 and hp_string[j] == 'H':
                            e += 1
        # Each contact counted twice → divide by 2; negate for convention
        return -(e // 2)

    # ── Initialization ────────────────────────────────────────────────────
    # Straight-line starting conformation: residue i at position (i, 0)
    current_positions = [(i, 0) for i in range(len(hp_string))]
    current_energy = calculate_energy(current_positions)

    # Global best tracker (separate from the random-walk current state)
    best_positions = list(current_positions)
    best_energy = current_energy
    trace_frames = []

    def add_trace_frame(step):
        if not return_trace:
            return
        if len(trace_frames) >= max_trace_frames:
            return
        trace_frames.append({
            "step": step,
            "energy": best_energy,
            "positions": list(best_positions),
        })

    trace_interval = max(1, int(trace_interval))
    max_trace_frames = max(2, int(max_trace_frames))
    add_trace_frame(0)

    # ── Main SA loop ──────────────────────────────────────────────────────
    for i in range(iterations):
        step = i + 1
        if len(hp_string) <= 2:
            break  # Chains too short to have any HH contacts

        # ── Exponential cooling schedule ──────────────────────────────────
        # Temperature decreases smoothly from initial_t to final_t:
        #   T_i = T₀ · (T_f / T₀)^(i / (iterations - 1))
        # At i=0: T = T₀ (maximum exploration)
        # At i=iterations-1: T = T_f (nearly greedy)
        if iterations > 1:
            fraction = i / float(iterations - 1)   # progress ratio [0, 1]
            T = initial_t * ((final_t / initial_t) ** fraction)
        else:
            T = final_t  # edge case: only 1 iteration

        # ── Move selection ─────────────────────────────────────────────────
        # Same weighted move set as HC: local moves are preferred
        move_types = ['end_flip', 'kink_jump', 'crankshaft', 'pivot']
        move = random.choices(move_types, weights=[0.2, 0.3, 0.3, 0.2])[0]

        new_positions = list(current_positions)  # working copy
        valid_move = False

        # ── END-FLIP ───────────────────────────────────────────────────────
        # Move one terminal residue to a free adjacent cell of its anchor
        if move == 'end_flip':
            idx = random.choice([0, len(hp_string) - 1])  # N- or C-terminus
            anchor_idx = 1 if idx == 0 else len(hp_string) - 2
            ax, ay = current_positions[anchor_idx]

            possible_moves = []
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nx, ny = ax + dx, ay + dy
                # Must be unoccupied AND different from the current terminal position
                if (nx, ny) not in current_positions and (nx, ny) != current_positions[idx]:
                    possible_moves.append((nx, ny))

            if possible_moves:
                new_positions[idx] = random.choice(possible_moves)
                valid_move = True

        # ── KINK-JUMP ──────────────────────────────────────────────────────
        # Flip an internal residue at a 90° bend to the opposite corner
        elif move == 'kink_jump':
            if len(hp_string) > 2:
                idx = random.randint(1, len(hp_string) - 2)
                p_prev = current_positions[idx - 1]
                p_curr = current_positions[idx]
                p_next = current_positions[idx + 1]

                # Diagonal neighbours → the residue is at a kink
                dx = abs(p_prev[0] - p_next[0])
                dy = abs(p_prev[1] - p_next[1])
                if dx == 1 and dy == 1:
                    # New position = reflection: p_prev + p_next - p_curr
                    nx = p_prev[0] + p_next[0] - p_curr[0]
                    ny = p_prev[1] + p_next[1] - p_curr[1]
                    if (nx, ny) not in current_positions:
                        new_positions[idx] = (nx, ny)
                        valid_move = True

        # ── CRANKSHAFT ─────────────────────────────────────────────────────
        # Rotate a U-shaped pair of consecutive residues 180°
        elif move == 'crankshaft':
            if len(hp_string) > 3:
                idx = random.randint(1, len(hp_string) - 3)
                p_prev  = current_positions[idx - 1]
                p_curr  = current_positions[idx]
                p_next  = current_positions[idx + 1]
                p_next2 = current_positions[idx + 2]

                # U-shape condition: p_prev and p_next2 must be adjacent
                dist_sq = (p_prev[0] - p_next2[0])**2 + (p_prev[1] - p_next2[1])**2
                if dist_sq == 1:
                    # Reflect both residues through the p_prev–p_next2 axis
                    nx_curr = p_prev[0] + p_next2[0] - p_next[0]
                    ny_curr = p_prev[1] + p_next2[1] - p_next[1]
                    nx_next = p_prev[0] + p_next2[0] - p_curr[0]
                    ny_next = p_prev[1] + p_next2[1] - p_curr[1]

                    # Both new cells must be free
                    if (nx_curr, ny_curr) not in current_positions and \
                       (nx_next, ny_next) not in current_positions:
                        new_positions[idx]     = (nx_curr, ny_curr)
                        new_positions[idx + 1] = (nx_next, ny_next)
                        valid_move = True

        # ── PIVOT ──────────────────────────────────────────────────────────
        # Rotate all residues after the pivot point by ±90° or 180°
        elif move == 'pivot':
            pivot_idx = random.randint(1, len(hp_string) - 2)
            angle = random.choice([90, -90, 180])
            cx, cy = current_positions[pivot_idx]  # pivot centre

            # 2D rotation matrix coefficients
            if angle == 90:
                cos_a, sin_a = 0, 1       # counter-clockwise 90°
            elif angle == -90:
                cos_a, sin_a = 0, -1      # clockwise 90°
            else:
                cos_a, sin_a = -1, 0      # 180° rotation

            for k in range(pivot_idx + 1, len(hp_string)):
                x, y = current_positions[k]
                tx, ty = x - cx, y - cy          # translate to pivot origin
                rx = tx * cos_a - ty * sin_a      # rotate x
                ry = tx * sin_a + ty * cos_a      # rotate y
                new_positions[k] = (rx + cx, ry + cy)  # translate back

            # Valid only if no two residues overlap (all positions unique)
            if len(set(new_positions)) == len(new_positions):
                valid_move = True

        # ── SA Acceptance criterion ────────────────────────────────────────
        if valid_move:
            new_energy = calculate_energy(new_positions)

            if new_energy <= current_energy:
                # Always accept improving or equal moves (same as HC)
                current_positions = new_positions
                current_energy = new_energy

                # Update global best only when strictly improving
                if current_energy < best_energy:
                    best_positions = list(current_positions)
                    best_energy = current_energy
            else:
                # ── Metropolis criterion ──────────────────────────────────
                # Accept a worsening move with probability exp(-ΔE / T).
                # When T is high → probability ≈ 1 (accept almost anything)
                # When T → 0    → probability → 0 (never accept worse moves)
                delta_e = new_energy - current_energy   # ΔE > 0 (worsening)
                probability = math.exp(-delta_e / T)    # Boltzmann factor

                if random.random() < probability:
                    # Accept the worse move to escape local minimum
                    current_positions = new_positions
                    current_energy = new_energy
                    # Note: global best is NOT updated here since energy worsened

        if step % trace_interval == 0:
            add_trace_frame(step)

    if return_trace and (not trace_frames or trace_frames[-1]["step"] != iterations):
        add_trace_frame(iterations)

    if return_trace:
        return best_positions, best_energy, trace_frames
    return best_positions, best_energy


# ── Script entry point ────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Protein sequences of increasing length for benchmarking
    p1  = "MKVILA"                                        #  6 AA
    p2  = "MKTIIALSYIFCLVF"                              # 15 AA
    p3  = "MVHLTPEEKSAVTALWGK"                           # 18 AA
    p4  = "MVHLTPEEKSAVTALWGKVN"                         # 20 AA
    p5  = "MVHLTPEEKSAVTALWGKVNVDEVG"                    # 25 AA
    p6  = "MVHLTPEEKSAVTALWGKVNVDEVGGEA"                 # 28 AA
    p7  = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLL"            # 33 AA
    p8  = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVY"         # 36 AA
    p9  = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWT"      # 39 AA
    p10 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRF"   # 42 AA

    # A longer realistic protein (64 AA) used in the final comparison tests
    p = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

    """
    proteins = [
        ("Protein 1 (6 AA)",  p1),  ...  ("Protein 10 (42 AA)", p10),
    ]
    """
    # For this run, only test the long protein
    proteins = [("Protein 11 (64 AA)", p)]

    iterations_to_test = [1000, 10000, 50000, 100000]

    print("=== TEST SIMULATED ANNEALING ALGORITHM WITH 10 PROTEINS ===\n")

    for name, seq in proteins:
        hp_str = sequence_to_hp(seq)
        h_count = hp_str.count('H')
        print(f"--- {name} ---")
        print(f"Sequence : {seq}")
        print(f"HP String: {hp_str} (H count: {h_count})\n")

        for iters in iterations_to_test:
            start_time = time.time()
            # Run SA with exponential cooling from 10.0 → 0.01
            positions, energy = generate_2d_structure_sa(hp_str, iterations=iters)
            end_time = time.time()
            elapsed_time = end_time - start_time
            print(f"  Iterations: {iters:<7} | Best Energy: {energy:<3} | Time: {elapsed_time:.4f} sec")

        print("\n" + "-"*50 + "\n")

"""
test_mc.py
----------
Implements the Basic Monte Carlo (MC) algorithm for 2D protein
structure prediction using the HP lattice model.

MC is similar to Simulated Annealing but uses a CONSTANT temperature
throughout the run. It applies the Metropolis acceptance criterion to
occasionally accept moves that worsen the energy, allowing limited
exploration around the current energy level — but without the
progressive cooling that makes SA so effective at minimization.
"""

import math
import random
import time


def sequence_to_hp(sequence):
    """
    Converts a raw amino acid sequence into an HP binary string.

    Each residue is classified as:
      - 'H' (Hydrophobic): belongs to {A, C, F, I, L, M, V, W, Y}
      - 'P' (Polar): all other residues

    Parameters
    ----------
    sequence : str  Raw amino acid sequence (one-letter codes)
    Returns  : str  HP binary string of the same length
    """
    hydrophobic = set('ACFILMVWY')  # Standard hydrophobic residue set
    # Build the HP string character by character using a list comprehension
    return ''.join(['H' if aa in hydrophobic else 'P' for aa in sequence])


def generate_2d_structure_mc(hp_string, iterations=100000, temperature=2.0,
                             return_trace=False, trace_interval=1000,
                             max_trace_frames=200):
    """
    Generates an approximate minimum-energy 2D conformation using
    Basic Monte Carlo at a CONSTANT temperature.

    Unlike Simulated Annealing, the temperature does not decrease over
    time. This allows thermodynamic sampling of the conformation ensemble
    at temperature T, but is less effective at finding the global minimum.

    Parameters
    ----------
    hp_string   : str   Binary HP sequence (e.g. "HPHHP")
    iterations  : int   Total number of move attempts (default: 100,000)
    temperature : float Constant Metropolis temperature (default: 2.0)

    Returns
    -------
    best_positions : list of (int, int)  Best lattice positions found
    best_energy    : int                 Minimum HP energy found (≤ 0)
    """
    if not hp_string:
        # Handle empty sequence gracefully
        if return_trace:
            return [], 0, []
        return [], 0

    def calculate_energy(pos_list):
        """
        Computes the HP energy of the current conformation.
        Energy = −(number of non-consecutive H-H topological contacts).
        Only pairs of H residues that are grid-adjacent (distance 1)
        AND non-consecutive in the chain (|i-j| > 1) count as contacts.
        """
        e = 0
        # Position → index dictionary for O(1) neighbour lookup
        pos_dict = {pos: i for i, pos in enumerate(pos_list)}
        for i in range(len(pos_list)):
            if hp_string[i] == 'H':   # Only H residues can form HH contacts
                x, y = pos_list[i]
                # Check all 4 orthogonal directions on the square lattice
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if (nx, ny) in pos_dict:
                        j = pos_dict[(nx, ny)]
                        # Ignore backbone bonds (|i-j| = 1); count only contacts
                        if abs(i - j) > 1 and hp_string[j] == 'H':
                            e += 1
        # Each contact is counted from both residues → divide by 2
        return -(e // 2)

    # ── Initialization ────────────────────────────────────────────────────
    # Start from a fully extended straight-line conformation: residue i at (i, 0)
    current_positions = [(i, 0) for i in range(len(hp_string))]
    current_energy = calculate_energy(current_positions)

    # Track the globally best solution encountered during the run
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

    # ── Main MC loop ──────────────────────────────────────────────────────
    for step in range(1, iterations + 1):
        if len(hp_string) <= 2:
            # Chains of length ≤ 2 cannot form any HH contacts
            break

        # ── Move selection ─────────────────────────────────────────────────
        # Sample one of 4 move types with weighted probabilities.
        # Local moves (kink-jump, crankshaft) are preferred because they
        # have higher acceptance rates in compact conformations.
        move_types = ['end_flip', 'kink_jump', 'crankshaft', 'pivot']
        move = random.choices(move_types, weights=[0.2, 0.3, 0.3, 0.2])[0]

        new_positions = list(current_positions)  # working copy
        valid_move = False

        # ── END-FLIP ───────────────────────────────────────────────────────
        # Relocate one terminal residue to a free neighbour of its anchor
        if move == 'end_flip':
            idx = random.choice([0, len(hp_string) - 1])   # pick N- or C-terminus
            anchor_idx = 1 if idx == 0 else len(hp_string) - 2  # bonded neighbour
            ax, ay = current_positions[anchor_idx]

            # Collect all unoccupied cells adjacent to the anchor
            possible_moves = []
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nx, ny = ax + dx, ay + dy
                # Cell must be free AND different from the current terminal cell
                if (nx, ny) not in current_positions and (nx, ny) != current_positions[idx]:
                    possible_moves.append((nx, ny))

            if possible_moves:
                new_positions[idx] = random.choice(possible_moves)
                valid_move = True

        # ── KINK-JUMP ──────────────────────────────────────────────────────
        # Move a residue at a 90° bend to the diagonally opposite corner
        elif move == 'kink_jump':
            if len(hp_string) > 2:
                idx = random.randint(1, len(hp_string) - 2)  # internal residue
                p_prev = current_positions[idx - 1]
                p_curr = current_positions[idx]
                p_next = current_positions[idx + 1]

                # Check if prev and next form a diagonal (kink condition)
                dx = abs(p_prev[0] - p_next[0])
                dy = abs(p_prev[1] - p_next[1])
                if dx == 1 and dy == 1:
                    # Reflect p_curr through the midpoint of (p_prev, p_next)
                    nx = p_prev[0] + p_next[0] - p_curr[0]
                    ny = p_prev[1] + p_next[1] - p_curr[1]
                    if (nx, ny) not in current_positions:
                        new_positions[idx] = (nx, ny)
                        valid_move = True

        # ── CRANKSHAFT ─────────────────────────────────────────────────────
        # Rotate a U-shaped pair of residues 180° around the U-axis
        elif move == 'crankshaft':
            if len(hp_string) > 3:
                idx = random.randint(1, len(hp_string) - 3)
                p_prev  = current_positions[idx - 1]
                p_curr  = current_positions[idx]
                p_next  = current_positions[idx + 1]
                p_next2 = current_positions[idx + 2]

                # U-shape: p_prev and p_next2 must be adjacent (squared dist = 1)
                dist_sq = (p_prev[0] - p_next2[0])**2 + (p_prev[1] - p_next2[1])**2
                if dist_sq == 1:
                    # Compute the reflected positions of both residues
                    nx_curr = p_prev[0] + p_next2[0] - p_next[0]
                    ny_curr = p_prev[1] + p_next2[1] - p_next[1]
                    nx_next = p_prev[0] + p_next2[0] - p_curr[0]
                    ny_next = p_prev[1] + p_next2[1] - p_curr[1]

                    # Both new positions must be unoccupied
                    if (nx_curr, ny_curr) not in current_positions and \
                       (nx_next, ny_next) not in current_positions:
                        new_positions[idx]     = (nx_curr, ny_curr)
                        new_positions[idx + 1] = (nx_next, ny_next)
                        valid_move = True

        # ── PIVOT ──────────────────────────────────────────────────────────
        # Rotate the tail of the chain around a pivot residue
        elif move == 'pivot':
            pivot_idx = random.randint(1, len(hp_string) - 2)
            angle = random.choice([90, -90, 180])
            cx, cy = current_positions[pivot_idx]  # pivot centre

            # 2D rotation matrix coefficients for the chosen angle
            if angle == 90:
                cos_a, sin_a = 0, 1    # 90° counter-clockwise
            elif angle == -90:
                cos_a, sin_a = 0, -1   # 90° clockwise
            else:
                cos_a, sin_a = -1, 0   # 180° rotation

            for k in range(pivot_idx + 1, len(hp_string)):
                x, y = current_positions[k]
                tx, ty = x - cx, y - cy          # translate to pivot origin
                rx = tx * cos_a - ty * sin_a      # rotate x-component
                ry = tx * sin_a + ty * cos_a      # rotate y-component
                new_positions[k] = (rx + cx, ry + cy)  # translate back

            # Valid only if no overlap in the resulting chain
            if len(set(new_positions)) == len(new_positions):
                valid_move = True

        # ── Acceptance criterion (constant-temperature Metropolis) ─────────
        if valid_move:
            new_energy = calculate_energy(new_positions)

            if new_energy <= current_energy:
                # Always accept improving or equal moves
                current_positions = new_positions
                current_energy = new_energy

                # Update global best only on strict improvement
                if current_energy < best_energy:
                    best_positions = list(current_positions)
                    best_energy = current_energy
            else:
                # ── Metropolis criterion at constant T ────────────────────
                # Accept a worsening move with probability exp(-ΔE / T).
                # Because T is CONSTANT (not decreasing), the acceptance rate
                # never drops → the algorithm keeps exploring at the same level,
                # which prevents convergence to deep minima.
                delta_e = new_energy - current_energy   # ΔE > 0 (energy worsened)
                probability = math.exp(-delta_e / temperature)  # Boltzmann factor

                if random.random() < probability:
                    # Accept the worse move (temperature-dependent randomness)
                    current_positions = new_positions
                    current_energy = new_energy
                    # Note: global best NOT updated here (energy worsened)

        if step % trace_interval == 0:
            add_trace_frame(step)

    if return_trace and (not trace_frames or trace_frames[-1]["step"] != iterations):
        add_trace_frame(iterations)

    if return_trace:
        return best_positions, best_energy, trace_frames
    return best_positions, best_energy


# ── Script entry point ────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Short test proteins for quick validation
    p1 = "MKVILA"                        #  6 AA
    p2 = "MKTIIALSYIFCLVF"              # 15 AA
    p3 = "MVHLTPEEKSAVTALWGKVNVDEVG"    # 25 AA

    proteins = [
        ("Protein 1 (6 AA)",  p1),
        ("Protein 2 (15 AA)", p2),
        ("Protein 3 (25 AA)", p3),
    ]

    # Test at two iteration budgets
    iterations_to_test = [10000, 50000]

    print("=== TEST BASIC MONTE CARLO ALGORITHM (CONSTANT TEMP) ===\n")

    for name, seq in proteins:
        hp_str = sequence_to_hp(seq)
        print(f"--- {name} ---")
        print(f"Sequence : {seq}")
        print(f"HP String: {hp_str}\n")

        for iters in iterations_to_test:
            start_time = time.time()
            # Run MC at constant temperature T=2.0
            positions, energy = generate_2d_structure_mc(hp_str, iterations=iters, temperature=2.0)
            elapsed = time.time() - start_time
            print(f"  Iterations: {iters:<7} | Best Energy: {energy:<3} | Time: {elapsed:.4f} sec")

        print("\n" + "-"*50 + "\n")

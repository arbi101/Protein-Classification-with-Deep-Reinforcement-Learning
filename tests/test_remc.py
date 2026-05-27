"""
test_remc.py
------------
Implements Replica Exchange Monte Carlo (REMC), also known as
Parallel Tempering, for 2D HP protein structure prediction.

REMC runs multiple independent Monte Carlo simulations (replicas) in
parallel, each at a different temperature. Periodically, adjacent
replicas attempt to SWAP their conformations. This allows low-temperature
replicas (which exploit good solutions) to benefit from high-temperature
configurations (which explore broadly), greatly improving the sampling
of the energy landscape and the quality of the final minimum found.
"""

import enum
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
    sequence : str  Raw one-letter amino acid sequence
    Returns  : str  HP binary string of the same length
    """
    hydrophobic = set('ACFILMVWY')  # Standard hydrophobic amino acid codes
    return ''.join(['H' if aa in hydrophobic else 'P' for aa in sequence])


def generate_2d_structure_remc(hp_string, iterations=300_000,
                                num_replicas=10, t_min=0.1, t_max=30.0,
                                swap_interval=200, return_trace=False,
                                trace_interval=1000, max_trace_frames=200):
    """
    Generates an approximate minimum-energy 2D protein conformation using
    Replica Exchange Monte Carlo (Parallel Tempering).

    REMC runs multiple independent MC simulations (replicas) in parallel,
    each at a different temperature. Periodically, adjacent replicas attempt
    to exchange their conformations based on a Metropolis criterion. Uses the
    same 4 move types as HC/SA/MC for fair comparison: end_flip, kink_jump,
    crankshaft, and pivot.

    Parameters
    ----------
    hp_string     : str   Binary HP sequence (e.g. "HPHHP")
    iterations    : int   Total MC steps per replica (default: 300,000)
    num_replicas  : int   Number of parallel replicas (default: 10)
    t_min         : float Lowest temperature replica (default: 1.0)
    t_max         : float Highest temperature replica (default: 30.0)
    swap_interval : int   Steps between replica-swap attempts (default: 200)

    Returns
    -------
    global_best_pos    : list of (int, int)  Best conformation found across all replicas
    global_best_energy : int                 Minimum energy found (≤ 0)
    """
    if not hp_string:
        # Empty sequence → nothing to fold
        if return_trace:
            return [], 0, []
        return [], 0
    if len(hp_string) <= 2:
        # Chains too short for any HH contacts → return straight line, energy 0
        positions = [(i, 0) for i in range(len(hp_string))]
        if return_trace:
            return positions, 0, [{"step": 0, "energy": 0, "positions": positions}]
        return positions, 0

    def calculate_energy(pos_list):
        """
        Computes the HP energy of a given conformation.
        Energy = −(number of non-consecutive H-H topological contacts).
        Each H residue checks its 4 grid neighbours; counts contact if
        the neighbour is also H and non-consecutive (|i-j| > 1).
        """
        e = 0
        # Map each position (x,y) → residue index for O(1) lookups
        pos_dict = {pos: i for i, pos in enumerate(pos_list)}
        for i in range(len(pos_list)):
            if hp_string[i] == 'H':   # Only H residues form HH contacts
                x, y = pos_list[i]
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if (nx, ny) in pos_dict:
                        j = pos_dict[(nx, ny)]
                        # Exclude covalently bonded neighbours (|i-j|=1)
                        if abs(i - j) > 1 and hp_string[j] == 'H':
                            e += 1
        # Each contact counted twice → divide by 2; negate for convention
        return -(e // 2)

    # ── Temperature ladder ────────────────────────────────────────────────
    # Distribute replica temperatures logarithmically between t_min and t_max.
    # Log spacing gives more replicas near the low end (where detail matters)
    # and fewer at high temperatures (where broad exploration is enough).
    if num_replicas > 1:
        temps = [t_min * (t_max / t_min) ** (i / (num_replicas - 1))
                 for i in range(num_replicas)]
    else:
        temps = [t_min]  # Single replica → just standard MC

    # ── Replica initialization ────────────────────────────────────────────
    # Each replica starts from the same straight-line conformation.
    # They will quickly diverge as each follows its own MC trajectory.
    replicas_pos = [[(j, 0) for j in range(len(hp_string))]
                    for _ in range(num_replicas)]
    # Pre-compute the initial energy for each replica
    replicas_en = [calculate_energy(r_pos) for r_pos in replicas_pos]

    # Global best tracks the single best conformation seen across ALL replicas
    global_best_pos    = list(replicas_pos[0])
    global_best_energy = replicas_en[0]
    trace_frames = []

    def add_trace_frame(step):
        if not return_trace:
            return
        if len(trace_frames) >= max_trace_frames:
            return
        trace_frames.append({
            "step": step,
            "energy": global_best_energy,
            "positions": list(global_best_pos),
        })

    trace_interval = max(1, int(trace_interval))
    max_trace_frames = max(2, int(max_trace_frames))
    add_trace_frame(0)

    # ── Main REMC loop ────────────────────────────────────────────────────
    for step in range(1, iterations + 1):

        # ── Local MC step for EACH replica ───────────────────────────────
        # Every replica performs one independent MC move at its own temperature.
        # Uses the same 4 move types as HC/SA/MC for fair comparison.
        for i in range(num_replicas):
            if len(hp_string) <= 2:
                # Chains too short → no optimization possible
                break

            # ── Move selection ─────────────────────────────────────────────
            # Choose one of 4 move types with weighted probabilities
            move_types = ['end_flip', 'kink_jump', 'crankshaft', 'pivot']
            move = random.choices(move_types, weights=[0.2, 0.3, 0.3, 0.2])[0]

            # Make a working copy of this replica's positions
            new_positions = list(replicas_pos[i])
            valid_move = False

            # ── END-FLIP move ──────────────────────────────────────────────
            if move == 'end_flip':
                idx = random.choice([0, len(hp_string) - 1])
                anchor_idx = 1 if idx == 0 else len(hp_string) - 2
                ax, ay = replicas_pos[i][anchor_idx]

                possible_moves = []
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = ax + dx, ay + dy
                    if (nx, ny) not in replicas_pos[i] and (nx, ny) != replicas_pos[i][idx]:
                        possible_moves.append((nx, ny))

                if possible_moves:
                    new_positions[idx] = random.choice(possible_moves)
                    valid_move = True

            # ── KINK-JUMP move ────────────────────────────────────────────
            elif move == 'kink_jump':
                if len(hp_string) > 2:
                    idx = random.randint(1, len(hp_string) - 2)
                    p_prev = replicas_pos[i][idx - 1]
                    p_curr = replicas_pos[i][idx]
                    p_next = replicas_pos[i][idx + 1]

                    dx = abs(p_prev[0] - p_next[0])
                    dy = abs(p_prev[1] - p_next[1])
                    if dx == 1 and dy == 1:
                        nx = p_prev[0] + p_next[0] - p_curr[0]
                        ny = p_prev[1] + p_next[1] - p_curr[1]
                        if (nx, ny) not in replicas_pos[i]:
                            new_positions[idx] = (nx, ny)
                            valid_move = True

            # ── CRANKSHAFT move ───────────────────────────────────────────
            elif move == 'crankshaft':
                if len(hp_string) > 3:
                    idx = random.randint(1, len(hp_string) - 3)
                    p_prev  = replicas_pos[i][idx - 1]
                    p_curr  = replicas_pos[i][idx]
                    p_next  = replicas_pos[i][idx + 1]
                    p_next2 = replicas_pos[i][idx + 2]

                    dist_sq = (p_prev[0] - p_next2[0])**2 + (p_prev[1] - p_next2[1])**2
                    if dist_sq == 1:
                        nx_curr = p_prev[0] + p_next2[0] - p_next[0]
                        ny_curr = p_prev[1] + p_next2[1] - p_next[1]
                        nx_next = p_prev[0] + p_next2[0] - p_curr[0]
                        ny_next = p_prev[1] + p_next2[1] - p_curr[1]

                        if (nx_curr, ny_curr) not in replicas_pos[i] and \
                           (nx_next, ny_next) not in replicas_pos[i]:
                            new_positions[idx]     = (nx_curr, ny_curr)
                            new_positions[idx + 1] = (nx_next, ny_next)
                            valid_move = True

            # ── PIVOT move ─────────────────────────────────────────────────
            elif move == 'pivot':
                pivot_idx = random.randint(1, len(hp_string) - 2)
                angle = random.choice([90, -90, 180])
                cx, cy = replicas_pos[i][pivot_idx]

                if angle == 90:
                    cos_a, sin_a = 0, 1
                elif angle == -90:
                    cos_a, sin_a = 0, -1
                else:
                    cos_a, sin_a = -1, 0

                for k in range(pivot_idx + 1, len(hp_string)):
                    x, y = replicas_pos[i][k]
                    tx, ty = x - cx, y - cy
                    rx = tx * cos_a - ty * sin_a
                    ry = tx * sin_a + ty * cos_a
                    new_positions[k] = (rx + cx, ry + cy)

                if len(set(new_positions)) == len(new_positions):
                    valid_move = True

            # ── Metropolis acceptance criterion at this replica's temperature ─
            if valid_move:
                new_energy = calculate_energy(new_positions)

                # Accept if energy improves or with probability exp(-ΔE / T)
                delta_e = new_energy - replicas_en[i]
                if delta_e <= 0 or random.random() < math.exp(-delta_e / temps[i]):
                    replicas_pos[i] = new_positions   # update replica conformation
                    replicas_en[i]  = new_energy      # update replica energy

                    # Check if this is the new global best across all replicas
                    if new_energy < global_best_energy:
                        global_best_pos    = list(new_positions)
                        global_best_energy = new_energy

        # ── Replica Exchange step ─────────────────────────────────────────
        # Every `swap_interval` steps, attempt to swap two adjacent replicas.
        # This is the key REMC mechanism: allows low-T replicas to receive
        # conformations discovered at high T and vice versa.
        if step % swap_interval == 0 and num_replicas > 1:
            # Pick a random adjacent pair of replicas (idx, idx+1)
            idx   = random.randint(0, num_replicas - 2)
            j_idx = idx + 1

            # ── Swap acceptance probability ───────────────────────────────
            # Derived from detailed balance:
            # P(swap) = min(1, exp((β_k − β_{k+1}) · (E_k − E_{k+1})))
            # where β_k = 1/T_k (inverse temperature).
            delta_e    = replicas_en[idx] - replicas_en[j_idx]   # ΔE = E_k - E_{k+1}
            delta_beta = (1.0 / temps[idx]) - (1.0 / temps[j_idx])  # Δβ = β_k - β_{k+1}
            exponent   = delta_beta * delta_e  # Combined exponent for the swap criterion

            # Accept swap always if exponent ≥ 0, otherwise with probability exp(exponent)
            if exponent >= 0 or random.random() < math.exp(exponent):
                # Swap BOTH positions and energies between the two replicas
                replicas_pos[idx], replicas_pos[j_idx] = \
                    replicas_pos[j_idx], replicas_pos[idx]
                replicas_en[idx], replicas_en[j_idx] = \
                    replicas_en[j_idx], replicas_en[idx]

        if step % trace_interval == 0:
            add_trace_frame(step)

    if return_trace and (not trace_frames or trace_frames[-1]["step"] != iterations):
        add_trace_frame(iterations)

    # Return the globally best conformation seen across all replicas and all steps
    if return_trace:
        return global_best_pos, global_best_energy, trace_frames
    return global_best_pos, global_best_energy


# ── Script entry point ────────────────────────────────────────────────────────
if __name__ == "__main__":
    # Short test protein (33 AA) — for quick validation
    p = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLL"

    # Longer realistic protein (74 AA) — for full REMC benchmark
    p2 = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

    # Convert to HP and run REMC
    hp_str = sequence_to_hp(p2)
    print("=== TEST OPTIMIZED REMC ===")
    print(f"Sequence : {p2}")
    print(f"HP String: {hp_str}")

    start = time.time()
    # Run REMC with 300,000 steps — default 20 replicas between T=0.1 and T=30.0
    _, energy = generate_2d_structure_remc(hp_str, iterations=300_000, num_replicas=20, t_min=0.1)
    print(f"Best Energy found: {energy} (Time: {time.time()-start:.2f} s)")

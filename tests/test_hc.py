"""
test_hc.py
----------
Implements the Hill Climbing (HC) algorithm for 2D protein structure
prediction using the HP (Hydrophobic-Polar) lattice model.

HC is a greedy local search: it only accepts moves that improve or
maintain the current energy (no temperature, no randomness beyond
the move selection itself).
"""

import random
import time


def sequence_to_hp(sequence):
    """
    Converts a raw amino acid sequence into an HP binary string.

    Each amino acid is classified as:
      - 'H' (Hydrophobic) if it belongs to the set {A, C, F, I, L, M, V, W, Y}
      - 'P' (Polar) for all other residues

    Parameters
    ----------
    sequence : str
        One-letter amino acid sequence (e.g. "MKVILA")

    Returns
    -------
    str
        HP string of the same length (e.g. "HHHHHH")
    """
    # Define the set of hydrophobic amino acid characters
    hydrophobic = set('ACFILMVWY')
    hp_string = ''
    for aa in sequence:
        # Map each amino acid to 'H' or 'P'
        if aa in hydrophobic:
            hp_string += 'H'
        else:
            hp_string += 'P'
    return hp_string


def generate_2d_structure(hp_string, iterations=1000000, return_trace=False,
                          trace_interval=1000, max_trace_frames=200):
    """
    Generates an approximate minimum-energy 2D protein conformation
    using Hill Climbing with a shared set of local moves.

    The protein is embedded on a 2D integer lattice as a self-avoiding walk.
    The algorithm starts from a straight-line conformation and iteratively
    proposes random moves, accepting only those that do not increase energy.

    Parameters
    ----------
    hp_string : str
        Binary HP sequence (e.g. "HPHPH")
    iterations : int
        Number of move attempts (default: 1,000,000)

    Returns
    -------
    best_positions : list of (int, int)
        The lattice coordinates of each residue in the best conformation found
    best_energy : int
        The minimum HP energy found (negative; more negative = more stable)
    """
    if not hp_string:
        # Edge case: empty sequence → return immediately
        if return_trace:
            return [], 0, []
        return [], 0

    def calculate_energy(pos_list):
        """
        Computes the HP energy of a given conformation.

        The energy is the negative count of non-consecutive H-H contacts:
        for every pair of H residues that are adjacent on the grid but
        NOT adjacent in the sequence (|i-j| > 1), we count -1.

        Parameters
        ----------
        pos_list : list of (int, int)
            Current lattice positions of all residues

        Returns
        -------
        int
            Energy value (≤ 0; lower is more stable)
        """
        e = 0
        # Build a lookup dict: position → residue index for O(1) neighbour checks
        pos_dict = {pos: i for i, pos in enumerate(pos_list)}
        for i in range(len(pos_list)):
            if hp_string[i] == 'H':      # Only H residues can form HH contacts
                x, y = pos_list[i]
                # Check all 4 orthogonal neighbours on the 2D lattice
                for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                    nx, ny = x + dx, y + dy
                    if (nx, ny) in pos_dict:
                        j = pos_dict[(nx, ny)]
                        # Count only non-consecutive pairs (|i-j|>1 avoids counting
                        # backbone bonds as contacts)
                        if abs(i - j) > 1 and hp_string[j] == 'H':
                            e += 1
        # Each contact is counted twice (once from each side) → divide by 2
        # Negate because more contacts = lower (better) energy
        return -(e // 2)

    # ── Initialization ──────────────────────────────────────────────────────
    # Start from a fully extended straight-line conformation: residue i at (i, 0)
    # This is always a valid self-avoiding walk and has energy = 0
    current_positions = [(i, 0) for i in range(len(hp_string))]
    current_energy = calculate_energy(current_positions)

    # Keep track of the globally best conformation seen during the search
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

    # ── Main Hill Climbing loop ──────────────────────────────────────────────
    for step in range(1, iterations + 1):
        if len(hp_string) <= 2:
            # Chains of 1 or 2 residues cannot form contacts → nothing to optimize
            break

        # ── Move selection ───────────────────────────────────────────────────
        # Choose one of 4 move types with weighted probabilities.
        # Local moves (kink, crankshaft) are favoured over global ones (pivot)
        # because pivot moves have a high rejection rate in compact structures.
        move_types = ['end_flip', 'kink_jump', 'crankshaft', 'pivot']
        move = random.choices(move_types, weights=[0.2, 0.3, 0.3, 0.2])[0]

        # Make a working copy of the current positions to modify
        new_positions = list(current_positions)
        valid_move = False   # Flag: becomes True only if the proposed move is legal

        # ── END-FLIP move ─────────────────────────────────────────────────────
        # Move the first OR last terminal residue to a free adjacent site of
        # its only bonded neighbour.
        if move == 'end_flip':
            # Randomly pick which terminus to move: 0 (N-term) or n-1 (C-term)
            idx = random.choice([0, len(hp_string) - 1])
            # The anchor is the residue bonded to the chosen terminus
            anchor_idx = 1 if idx == 0 else len(hp_string) - 2
            ax, ay = current_positions[anchor_idx]

            # Collect all free neighbouring cells of the anchor
            possible_moves = []
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nx, ny = ax + dx, ay + dy
                # The target cell must be: (a) unoccupied AND (b) different from current position
                if (nx, ny) not in current_positions and (nx, ny) != current_positions[idx]:
                    possible_moves.append((nx, ny))

            if possible_moves:
                # Place the terminal residue at a randomly chosen free site
                new_positions[idx] = random.choice(possible_moves)
                valid_move = True

        # ── KINK-JUMP move ────────────────────────────────────────────────────
        # Applicable when an internal residue forms a 90° corner ("kink").
        # Condition: the previous and next residues differ by 1 in BOTH x and y
        # (they form a diagonal), meaning the current residue is at a corner.
        # The move flips the residue to the diagonally opposite corner of the
        # 2×2 square defined by prev, curr, and next.
        elif move == 'kink_jump':
            if len(hp_string) > 2:
                # Pick a random internal residue (not the first or last)
                idx = random.randint(1, len(hp_string) - 2)
                p_prev = current_positions[idx - 1]  # residue before
                p_curr = current_positions[idx]      # residue to move
                p_next = current_positions[idx + 1]  # residue after

                # Check if it's a kink: both coordinates differ by exactly 1
                dx = abs(p_prev[0] - p_next[0])
                dy = abs(p_prev[1] - p_next[1])
                if dx == 1 and dy == 1:
                    # The new position is the reflection of p_curr through the midpoint
                    # of (p_prev, p_next): new = p_prev + p_next - p_curr
                    nx = p_prev[0] + p_next[0] - p_curr[0]
                    ny = p_prev[1] + p_next[1] - p_curr[1]
                    # Accept only if the target cell is free (no overlap)
                    if (nx, ny) not in current_positions:
                        new_positions[idx] = (nx, ny)
                        valid_move = True

        # ── CRANKSHAFT move ───────────────────────────────────────────────────
        # Applicable when two consecutive internal residues form a U-shape:
        # residue i-1 and i+2 are adjacent (distance 1), while i and i+1
        # are in the "belly" of the U. Rotating the U 180° gives a new valid
        # conformation.
        elif move == 'crankshaft':
            if len(hp_string) > 3:
                # Pick a random starting index for the pair (i, i+1)
                idx = random.randint(1, len(hp_string) - 3)
                p_prev  = current_positions[idx - 1]
                p_curr  = current_positions[idx]
                p_next  = current_positions[idx + 1]
                p_next2 = current_positions[idx + 2]

                # Check U-shape: squared distance between p_prev and p_next2 must be 1
                dist_sq = (p_prev[0] - p_next2[0])**2 + (p_prev[1] - p_next2[1])**2
                if dist_sq == 1:
                    # Compute new positions by reflecting both residues through the
                    # axis defined by p_prev and p_next2
                    nx_curr = p_prev[0] + p_next2[0] - p_next[0]
                    ny_curr = p_prev[1] + p_next2[1] - p_next[1]
                    nx_next = p_prev[0] + p_next2[0] - p_curr[0]
                    ny_next = p_prev[1] + p_next2[1] - p_curr[1]

                    # Both new cells must be free (no self-overlap)
                    if (nx_curr, ny_curr) not in current_positions and \
                       (nx_next, ny_next) not in current_positions:
                        new_positions[idx]     = (nx_curr, ny_curr)
                        new_positions[idx + 1] = (nx_next, ny_next)
                        valid_move = True

        # ── PIVOT move ────────────────────────────────────────────────────────
        # Select a random pivot residue k. Rotate all residues from k+1 to n-1
        # by ±90° or 180° around the pivot point p_k.
        # The rotation is applied using a 2D rotation matrix.
        elif move == 'pivot':
            # Pick a random pivot (not the first or last residue)
            pivot_idx = random.randint(1, len(hp_string) - 2)
            angle = random.choice([90, -90, 180])
            cx, cy = current_positions[pivot_idx]  # pivot coordinates

            # Set rotation matrix coefficients for the chosen angle
            if angle == 90:
                cos_a, sin_a = 0, 1      # 90° counter-clockwise
            elif angle == -90:
                cos_a, sin_a = 0, -1     # 90° clockwise
            else:
                cos_a, sin_a = -1, 0     # 180° rotation

            # Apply rotation to all residues after the pivot point
            for i in range(pivot_idx + 1, len(hp_string)):
                x, y = current_positions[i]
                # Translate to pivot-centred coordinate system
                tx, ty = x - cx, y - cy
                # Apply 2D rotation: r = R · t
                rx = tx * cos_a - ty * sin_a
                ry = tx * sin_a + ty * cos_a
                # Translate back to global coordinates
                new_positions[i] = (rx + cx, ry + cy)

            # The pivot move is valid only if the resulting chain has no overlaps
            # (all n positions must be distinct)
            if len(set(new_positions)) == len(new_positions):
                valid_move = True

        # ── Acceptance criterion (Hill Climbing) ─────────────────────────────
        # Only accept the move if it leads to a non-increasing energy.
        # This is the key difference from SA: HC NEVER accepts worse moves.
        if valid_move:
            new_energy = calculate_energy(new_positions)

            # Accept if the new energy is equal or better (lower)
            if new_energy <= current_energy:
                current_positions = new_positions
                current_energy = new_energy

                # Update the global best solution if we found a new minimum
                if current_energy < best_energy:
                    best_positions = list(current_positions)
                    best_energy = current_energy

        if step % trace_interval == 0:
            add_trace_frame(step)

    if return_trace and (not trace_frames or trace_frames[-1]["step"] != iterations):
        add_trace_frame(iterations)

    if return_trace:
        return best_positions, best_energy, trace_frames
    return best_positions, best_energy


# ── Script entry point ───────────────────────────────────────────────────────
if __name__ == "__main__":
    # Define 10 protein sequences of increasing length for benchmarking
    p1 = "MKVILA"                                        #  6 AA
    p2 = "MKTIIALSYIFCLVF"                              # 15 AA
    p3 = "MVHLTPEEKSAVTALWGK"                           # 18 AA
    p4 = "MVHLTPEEKSAVTALWGKVN"                         # 20 AA
    p5 = "MVHLTPEEKSAVTALWGKVNVDEVG"                    # 25 AA
    p6 = "MVHLTPEEKSAVTALWGKVNVDEVGGEA"                 # 28 AA
    p7 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLL"            # 33 AA
    p8 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVY"         # 36 AA
    p9 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWT"      # 39 AA
    p10 = "MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQRF"  # 42 AA

    proteins = [
        ("Protein 1 (6 AA)",  p1),
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

    # Test each protein at 4 different iteration budgets to study convergence
    iterations_to_test = [1000, 10000, 50000, 100000]

    print("=== TEST HILL CLIMBING ALGORITHM WITH 10 PROTEINS ===\n")

    for name, seq in proteins:
        # Convert the amino acid sequence to HP binary string
        hp_str = sequence_to_hp(seq)

        # Count H residues (determines the theoretical minimum energy)
        h_count = hp_str.count('H')
        print(f"--- {name} ---")
        print(f"Sequence : {seq}")
        print(f"HP String: {hp_str} (H count: {h_count})\n")

        best_energies = []
        times_taken = []

        for iters in iterations_to_test:
            start_time = time.time()
            # Run Hill Climbing with the given iteration budget
            positions, energy = generate_2d_structure(hp_str, iterations=iters)
            end_time = time.time()

            elapsed_time = end_time - start_time
            best_energies.append(energy)
            times_taken.append(elapsed_time)

            print(f"  Iterations: {iters:<7} | Best Energy: {energy:<3} | Time: {elapsed_time:.4f} sec")

        print("\n" + "-"*50 + "\n")

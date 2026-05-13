"""
test_ql.py
----------
Implements Tabular Q-Learning for 2D HP protein structure prediction.

Q-Learning is a model-free reinforcement learning algorithm. The agent
maintains a table of Q-values Q(state, action) and updates them using
the Bellman equation after each step. It uses an ε-greedy policy to
balance exploration (trying random moves) and exploitation (using the
best known move).

This approach is effective for short sequences (≤ ~25 AA). For longer
chains, the state space grows exponentially, making the Q-table too
large to store in memory (the "curse of dimensionality").
"""

import random
import math
import time
from collections import defaultdict


# ─────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ─────────────────────────────────────────────────────────────────────────────

def sequence_to_hp(sequence):
    """
    Converts a raw amino acid sequence into an HP binary string.

    Classification:
      - 'H' (Hydrophobic): belongs to {A, C, F, I, L, M, V, W, Y}
      - 'P' (Polar): all other residues

    Parameters
    ----------
    sequence : str  One-letter amino acid sequence (e.g. "MKVILA")
    Returns  : str  HP binary string (e.g. "HHHHHH")
    """
    hydrophobic = set('ACFILMVWY')  # Standard hydrophobic residue set
    # Map each amino acid character to 'H' or 'P'
    return ''.join('H' if aa in hydrophobic else 'P' for aa in sequence)


def calculate_energy(pos_list, hp_string):
    """
    Computes the HP energy of a given 2D lattice conformation.

    The energy equals the negative count of non-consecutive H-H contacts:
    for each pair (i, j) where both residues are 'H', are adjacent on the
    grid (Manhattan distance = 1), and are NOT bonded in the chain
    (|i-j| > 1), we count one contact (energy contribution: -1).

    Parameters
    ----------
    pos_list  : list of (int, int)  Current lattice positions of all residues
    hp_string : str                 HP binary string

    Returns
    -------
    int  HP energy (≤ 0; more negative = more stable)
    """
    e = 0
    # Build a position → index dictionary for O(1) neighbour lookup
    pos_dict = {pos: i for i, pos in enumerate(pos_list)}
    for i, (x, y) in enumerate(pos_list):
        if hp_string[i] == 'H':   # Only H residues form HH contacts
            # Check all 4 orthogonal directions on the 2D grid
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nb = (x + dx, y + dy)    # neighbouring cell
                if nb in pos_dict:
                    j = pos_dict[nb]
                    # Count only non-consecutive pairs (exclude backbone bonds)
                    if abs(i - j) > 1 and hp_string[j] == 'H':
                        e += 1
    # Each contact is counted from BOTH residues → divide by 2
    return -(e // 2)


# ─────────────────────────────────────────────────────────────────────────────
# State representation
# ─────────────────────────────────────────────────────────────────────────────

def positions_to_state(pos_list):
    """
    Encodes the current protein conformation as a canonical, hashable tuple
    suitable for use as a Q-table key.

    The encoding removes two types of symmetry to reduce state space size:
      1. TRANSLATION symmetry: the chain is translated so that the first
         residue is always at the origin (0, 0).
      2. ROTATION symmetry: the chain is rotated so that the second residue
         always points in the +x direction.

    This means two conformations that are geometrically equivalent (just
    shifted or rotated) will map to the SAME state key, effectively
    reducing the Q-table size by up to 4×.

    Parameters
    ----------
    pos_list : list of (int, int)  Raw lattice positions

    Returns
    -------
    tuple  Canonical position tuple (hashable, usable as dict key)
    """
    if len(pos_list) < 1:
        return ()  # Edge case: empty chain

    # Step 1: Translate so the first residue is at the origin
    ox, oy = pos_list[0]
    translated = [(x - ox, y - oy) for x, y in pos_list]

    # Step 2: Rotate so the second residue is always at (+1, 0) (positive x-axis)
    if len(translated) > 1:
        dx, dy = translated[1]
        if dx == 0 and dy == 1:       # currently pointing UP → rotate -90°
            translated = [(y, -x) for x, y in translated]
        elif dx == -1 and dy == 0:    # currently pointing LEFT → rotate 180°
            translated = [(-x, -y) for x, y in translated]
        elif dx == 0 and dy == -1:    # currently pointing DOWN → rotate +90°
            translated = [(-y, x) for x, y in translated]
        # If already at (+1, 0): no rotation needed

    return tuple(translated)  # Return as immutable tuple (hashable)


# ─────────────────────────────────────────────────────────────────────────────
# Move set (same moves as HC / SA / MC for consistency)
# ─────────────────────────────────────────────────────────────────────────────

# The 4 available actions the RL agent can choose from
MOVE_TYPES = ['end_flip', 'kink_jump', 'crankshaft', 'pivot']


def apply_move(current_positions, move_type, hp_string):
    """
    Applies one of the 4 move types to the current protein conformation.

    Each move modifies the chain while preserving connectivity
    (consecutive residues remain adjacent) and self-avoidance
    (no two residues at the same position).

    Parameters
    ----------
    current_positions : list of (int, int)  Current lattice positions
    move_type         : str                 One of MOVE_TYPES
    hp_string         : str                 HP binary string (for length info)

    Returns
    -------
    list of (int, int) if the move is valid, or None if invalid/no valid position
    """
    n = len(hp_string)
    new_positions = list(current_positions)  # work on a copy

    # ── END-FLIP ─────────────────────────────────────────────────────────────
    # Move the first or last residue to a free adjacent cell of its anchor
    if move_type == 'end_flip':
        idx = random.choice([0, n - 1])            # pick N- or C-terminus
        anchor_idx = 1 if idx == 0 else n - 2      # the bonded anchor residue
        ax, ay = current_positions[anchor_idx]

        # Collect all unoccupied cells adjacent to the anchor
        candidates = []
        pos_set = set(current_positions)           # for O(1) membership test
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = ax + dx, ay + dy
            # Must be unoccupied AND different from current terminal position
            if (nx, ny) not in pos_set and (nx, ny) != current_positions[idx]:
                candidates.append((nx, ny))

        if not candidates:
            return None  # No valid position available → move is impossible
        new_positions[idx] = random.choice(candidates)
        return new_positions

    # ── KINK-JUMP ─────────────────────────────────────────────────────────────
    # Flip a residue at a 90° corner to the diagonally opposite position
    elif move_type == 'kink_jump':
        if n <= 2:
            return None  # Need at least 3 residues for a kink

        idx = random.randint(1, n - 2)             # random internal residue
        p_prev = current_positions[idx - 1]
        p_curr = current_positions[idx]
        p_next = current_positions[idx + 1]

        # Kink condition: prev and next differ by 1 in BOTH x and y (diagonal)
        if abs(p_prev[0] - p_next[0]) == 1 and abs(p_prev[1] - p_next[1]) == 1:
            # New position = reflection of p_curr through midpoint of (p_prev, p_next)
            nx = p_prev[0] + p_next[0] - p_curr[0]
            ny = p_prev[1] + p_next[1] - p_curr[1]
            if (nx, ny) not in set(current_positions):
                new_positions[idx] = (nx, ny)
                return new_positions
        return None  # Not a kink position, or target cell occupied

    # ── CRANKSHAFT ─────────────────────────────────────────────────────────────
    # Rotate a U-shaped pair of residues 180° around the U-axis
    elif move_type == 'crankshaft':
        if n <= 3:
            return None  # Need at least 4 residues for a crankshaft

        idx = random.randint(1, n - 3)             # pick start of the pair
        p_prev  = current_positions[idx - 1]
        p_curr  = current_positions[idx]
        p_next  = current_positions[idx + 1]
        p_next2 = current_positions[idx + 2]

        # U-shape condition: p_prev and p_next2 must be adjacent (distance = 1)
        dist_sq = (p_prev[0] - p_next2[0])**2 + (p_prev[1] - p_next2[1])**2
        if dist_sq == 1:
            # Reflect both residues through the axis defined by p_prev and p_next2
            nx_curr = p_prev[0] + p_next2[0] - p_next[0]
            ny_curr = p_prev[1] + p_next2[1] - p_next[1]
            nx_next = p_prev[0] + p_next2[0] - p_curr[0]
            ny_next = p_prev[1] + p_next2[1] - p_curr[1]

            pos_set = set(current_positions)
            # Both new positions must be unoccupied
            if (nx_curr, ny_curr) not in pos_set and (nx_next, ny_next) not in pos_set:
                new_positions[idx]     = (nx_curr, ny_curr)
                new_positions[idx + 1] = (nx_next, ny_next)
                return new_positions
        return None

    # ── PIVOT ──────────────────────────────────────────────────────────────────
    # Rotate all residues after a pivot point by ±90° or 180°
    elif move_type == 'pivot':
        if n <= 2:
            return None

        pivot_idx = random.randint(1, n - 2)       # random pivot residue
        angle = random.choice([90, -90, 180])
        cx, cy = current_positions[pivot_idx]       # pivot centre

        # 2D rotation matrix coefficients for the chosen angle
        if angle == 90:
            cos_a, sin_a = 0, 1    # 90° counter-clockwise
        elif angle == -90:
            cos_a, sin_a = 0, -1   # 90° clockwise
        else:
            cos_a, sin_a = -1, 0   # 180° rotation

        # Apply rotation to all residues after the pivot
        for i in range(pivot_idx + 1, n):
            x, y = current_positions[i]
            tx, ty = x - cx, y - cy          # translate to pivot origin
            new_positions[i] = (tx * cos_a - ty * sin_a + cx,
                                tx * sin_a + ty * cos_a + cy)

        # Check that the resulting chain has no self-intersections
        if len(set(new_positions)) == len(new_positions):
            return new_positions
        return None  # Overlap detected → invalid move

    return None  # Unknown move type


# ─────────────────────────────────────────────────────────────────────────────
# Tabular Q-Learning — core algorithm
# ─────────────────────────────────────────────────────────────────────────────

def generate_2d_structure_ql(
    hp_string,
    episodes=500,
    max_steps_per_episode=200,
    alpha=0.1,           # learning rate α: controls how fast Q-values are updated
    gamma=0.95,          # discount factor γ: importance of future rewards
    epsilon_start=1.0,   # initial exploration rate (100% random actions)
    epsilon_end=0.05,    # final exploration rate (5% random actions)
    epsilon_decay=None,  # number of steps over which ε decays (default: auto)
):
    """
    Tabular Q-Learning agent for 2D HP protein folding.

    The agent starts each episode from a straight-line conformation and
    learns which moves lead to lower-energy structures. It uses an
    ε-greedy policy (explore randomly with prob ε, exploit best-known
    action with prob 1−ε) and updates Q-values using the Bellman equation:

        Q(s,a) ← Q(s,a) + α [ R(s,a) + γ · max_a' Q(s',a') − Q(s,a) ]

    Reward design:
      - R = E_old − E_new  (positive if energy decreased = improved)
      - R = −1             (negative if the move was invalid/collision)

    Parameters
    ----------
    hp_string              : str   Binary HP sequence
    episodes               : int   Number of training episodes (default: 500)
    max_steps_per_episode  : int   Max moves per episode (default: 200)
    alpha                  : float Learning rate (default: 0.1)
    gamma                  : float Discount factor (default: 0.95)
    epsilon_start          : float Initial ε for ε-greedy (default: 1.0)
    epsilon_end            : float Final ε for ε-greedy (default: 0.05)
    epsilon_decay          : int   Steps over which ε decays (default: auto)

    Returns
    -------
    best_positions : list of (int, int)  Best conformation found
    best_energy    : int                 Minimum energy found
    q_table        : defaultdict         Learned Q-table {state: {action: value}}
    """
    if not hp_string or len(hp_string) <= 1:
        # Edge case: sequences too short to fold
        return [(0, 0)] * len(hp_string), 0, {}

    n = len(hp_string)
    actions = MOVE_TYPES  # The 4 possible actions the agent can take

    # ── Q-table initialization ────────────────────────────────────────────
    # defaultdict automatically creates a new entry {action: 0.0} for any
    # state key that hasn't been seen before (avoids KeyError exceptions)
    q_table = defaultdict(lambda: {a: 0.0 for a in actions})

    # ── Epsilon decay schedule ────────────────────────────────────────────
    # ε decays LINEARLY from epsilon_start to epsilon_end over the first
    # half of training. After that, ε stays at epsilon_end.
    total_steps = episodes * max_steps_per_episode
    if epsilon_decay is None:
        epsilon_decay = total_steps // 2   # decay over the first half of training

    step_counter = 0  # global step counter (across all episodes)

    # ── Global best tracker ───────────────────────────────────────────────
    # Initialize the best solution as the straight-line conformation
    best_positions = [(i, 0) for i in range(n)]
    best_energy    = calculate_energy(best_positions, hp_string)

    # ── Training loop ─────────────────────────────────────────────────────
    for episode in range(episodes):
        # Reset to straight-line conformation at the START of every episode.
        # The agent must re-learn to fold from scratch each time, but the
        # Q-table carries over all learned knowledge from previous episodes.
        current_positions = [(i, 0) for i in range(n)]
        current_energy    = calculate_energy(current_positions, hp_string)
        current_state     = positions_to_state(current_positions)

        for _ in range(max_steps_per_episode):
            # ── ε-greedy action selection ─────────────────────────────────
            # Compute current ε: linearly decays from epsilon_start to epsilon_end
            epsilon = max(
                epsilon_end,
                epsilon_start - (epsilon_start - epsilon_end) * step_counter / epsilon_decay
            )
            step_counter += 1  # increment global step count

            if random.random() < epsilon:
                # EXPLORE: pick a random action (ignores Q-values)
                action = random.choice(actions)
            else:
                # EXPLOIT: pick the action with the highest Q-value
                action = max(actions, key=lambda a: q_table[current_state][a])

            # ── Apply the chosen move ──────────────────────────────────────
            new_positions = apply_move(current_positions, action, hp_string)

            if new_positions is None:
                # INVALID move (collision or geometrically impossible):
                # Give a negative reward and stay in the same state
                reward      = -1
                next_state  = current_state
                next_energy = current_energy
            else:
                # VALID move: compute new energy and reward
                next_energy = calculate_energy(new_positions, hp_string)
                # Reward = energy improvement (positive = better structure)
                # E.g. if E went from -5 to -7: reward = -5 - (-7) = +2
                reward      = float(current_energy - next_energy)
                next_state  = positions_to_state(new_positions)

            # ── Q-Learning (Bellman) update ────────────────────────────────
            # TD target = R + γ · max_a' Q(s', a')
            # TD error  = TD target - Q(s, a)
            # New Q     = old Q + α · TD error
            best_next_q = max(q_table[next_state].values())  # greedy future value
            old_q       = q_table[current_state][action]
            q_table[current_state][action] = (
                old_q + alpha * (reward + gamma * best_next_q - old_q)
            )

            # ── State transition ──────────────────────────────────────────
            if new_positions is not None:
                # Move to the new state only if the move was valid
                current_positions = new_positions
                current_energy    = next_energy
                current_state     = next_state

                # Update global best if this conformation is better than all previous
                if current_energy < best_energy:
                    best_energy    = current_energy
                    best_positions = list(current_positions)

    return best_positions, best_energy, q_table


# ─────────────────────────────────────────────────────────────────────────────
# Greedy exploitation pass using the learned Q-table
# ─────────────────────────────────────────────────────────────────────────────

def exploit_q_table(hp_string, q_table, steps=500):
    """
    Runs a purely greedy pass (ε = 0) using the TRAINED Q-table.
    At each step, the agent always picks the highest Q-value action.
    This extracts the best folding trajectory the agent has learned.

    Parameters
    ----------
    hp_string : str          HP binary sequence
    q_table   : defaultdict  Trained Q-table from generate_2d_structure_ql
    steps     : int          Number of greedy steps (default: 500)

    Returns
    -------
    best_positions : list of (int, int)  Best conformation found in greedy pass
    best_energy    : int                 Minimum energy in greedy pass
    """
    n = len(hp_string)
    actions = MOVE_TYPES

    # Start from the same straight-line conformation as training
    current_positions = [(i, 0) for i in range(n)]
    current_energy    = calculate_energy(current_positions, hp_string)
    best_positions    = list(current_positions)
    best_energy       = current_energy

    for _ in range(steps):
        # Encode current conformation as a state key
        state  = positions_to_state(current_positions)
        # Greedy action: always pick the action with the highest Q-value
        action = max(actions, key=lambda a: q_table[state][a])

        new_positions = apply_move(current_positions, action, hp_string)
        if new_positions is None:
            continue  # Skip invalid moves silently (no state change)

        new_energy = calculate_energy(new_positions, hp_string)
        current_positions = new_positions
        current_energy    = new_energy

        # Update best if this step improved the conformation
        if current_energy < best_energy:
            best_energy    = current_energy
            best_positions = list(current_positions)

    return best_positions, best_energy


# ─────────────────────────────────────────────────────────────────────────────
# Main test – small proteins only (≤ 25 AA recommended for tabular QL)
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Small proteins well-suited for Tabular Q-Learning
    # (larger sequences cause the Q-table to grow too large in memory)
    p1 = "MKVILA"                     #  6 AA
    p2 = "MKTIIALSYIFCLVF"           # 15 AA
    p3 = "MVHLTPEEKSAVTALWGK"        # 18 AA
    p4 = "MVHLTPEEKSAVTALWGKVN"      # 20 AA
    p5 = "MVHLTPEEKSAVTALWGKVNVDEVG" # 25 AA

    proteins = [
        ("Protein 1 (6 AA)",  p1),
        ("Protein 2 (15 AA)", p2),
        ("Protein 3 (18 AA)", p3),
        ("Protein 4 (20 AA)", p4),
        ("Protein 5 (25 AA)", p5),
    ]

    # Different episode × step configurations to study convergence
    # total_steps = episodes × max_steps
    episode_configs = [
        {"episodes": 200,  "max_steps": 100},   # 20,000 total steps
        {"episodes": 500,  "max_steps": 200},   # 100,000 total steps
        {"episodes": 1000, "max_steps": 300},   # 300,000 total steps
        {"episodes": 2000, "max_steps": 400},   # 800,000 total steps
    ]

    print("=" * 65)
    print("  TABULAR Q-LEARNING  –  HP Model 2-D Protein Folding")
    print("  alpha=0.1  |  gamma=0.95  |  epsilon: 1.0 → 0.05")
    print("=" * 65)
    print()

    for name, seq in proteins:
        # Convert amino acid sequence to binary HP string
        hp_str  = sequence_to_hp(seq)
        h_count = hp_str.count('H')   # number of H residues (max possible contacts)
        print(f"--- {name} ---")
        print(f"  Sequence  : {seq}")
        print(f"  HP String : {hp_str}  (H: {h_count})")
        print()

        for cfg in episode_configs:
            eps  = cfg["episodes"]
            stps = cfg["max_steps"]

            start = time.time()
            # Train the Q-Learning agent with these episode/step settings
            positions, energy, q_table = generate_2d_structure_ql(
                hp_str,
                episodes=eps,
                max_steps_per_episode=stps,
                alpha=0.1,
                gamma=0.95,
            )
            elapsed = time.time() - start

            # Run a purely greedy exploitation pass after training to squeeze
            # out the best policy the agent learned
            ex_positions, ex_energy = exploit_q_table(hp_str, q_table, steps=1000)
            # Take the best result between training run and exploitation pass
            final_energy = min(energy, ex_energy)

            total_steps = eps * stps  # total number of environment steps
            print(
                f"  episodes={eps:<5} steps/ep={stps:<4} "
                f"| total_steps={total_steps:<8} "
                f"| Best E={final_energy:<4} "
                f"| Q-states={len(q_table):<6} "   # number of unique states visited
                f"| Time={elapsed:.3f}s"
            )

        print()
        print("-" * 65)
        print()

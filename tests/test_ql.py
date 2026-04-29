import random
import math
import time
from collections import defaultdict

# ─────────────────────────────────────────────────────────────────────────────
# Utility helpers
# ─────────────────────────────────────────────────────────────────────────────

def sequence_to_hp(sequence):
    """
    Converts an amino acid sequence to an HP string.
    H: Hydrophobic  |  P: Polar
    """
    hydrophobic = set('ACFILMVWY')
    return ''.join('H' if aa in hydrophobic else 'P' for aa in sequence)


def calculate_energy(pos_list, hp_string):
    """
    HP model energy: counts non-consecutive H-H contacts on the 2-D lattice.
    Each contact contributes -1 to the energy.
    """
    e = 0
    pos_dict = {pos: i for i, pos in enumerate(pos_list)}
    for i, (x, y) in enumerate(pos_list):
        if hp_string[i] == 'H':
            for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
                nb = (x + dx, y + dy)
                if nb in pos_dict:
                    j = pos_dict[nb]
                    if abs(i - j) > 1 and hp_string[j] == 'H':
                        e += 1
    return -(e // 2)


# ─────────────────────────────────────────────────────────────────────────────
# State representation
# ─────────────────────────────────────────────────────────────────────────────

def positions_to_state(pos_list):
    """
    Encodes the current conformation as a canonical (hashable) tuple.

    We translate the chain so that the first residue is at the origin and
    the second residue is always to the right (+x), removing translation and
    one rotation symmetry so that symmetrically equivalent conformations map
    to the same state.
    """
    if len(pos_list) < 1:
        return ()
    ox, oy = pos_list[0]
    translated = [(x - ox, y - oy) for x, y in pos_list]
    # Enforce canonical orientation: second residue on positive x-axis
    if len(translated) > 1:
        dx, dy = translated[1]
        if dx == 0 and dy == 1:          # rotated 90° → rotate -90°
            translated = [(y, -x) for x, y in translated]
        elif dx == -1 and dy == 0:       # rotated 180°
            translated = [(-x, -y) for x, y in translated]
        elif dx == 0 and dy == -1:       # rotated 270° → rotate 90°
            translated = [(-y, x) for x, y in translated]
    return tuple(translated)


# ─────────────────────────────────────────────────────────────────────────────
# Move set  (same moves as HC / SA / MC for consistency)
# ─────────────────────────────────────────────────────────────────────────────

MOVE_TYPES = ['end_flip', 'kink_jump', 'crankshaft', 'pivot']

def apply_move(current_positions, move_type, hp_string):
    """
    Applies one of the 4 move types to the chain.
    Returns new_positions if the move is valid, or None otherwise.
    """
    n = len(hp_string)
    new_positions = list(current_positions)

    if move_type == 'end_flip':
        idx = random.choice([0, n - 1])
        anchor_idx = 1 if idx == 0 else n - 2
        ax, ay = current_positions[anchor_idx]
        candidates = []
        pos_set = set(current_positions)
        for dx, dy in [(0, 1), (1, 0), (0, -1), (-1, 0)]:
            nx, ny = ax + dx, ay + dy
            if (nx, ny) not in pos_set and (nx, ny) != current_positions[idx]:
                candidates.append((nx, ny))
        if not candidates:
            return None
        new_positions[idx] = random.choice(candidates)
        return new_positions

    elif move_type == 'kink_jump':
        if n <= 2:
            return None
        idx = random.randint(1, n - 2)
        p_prev = current_positions[idx - 1]
        p_curr = current_positions[idx]
        p_next = current_positions[idx + 1]
        if abs(p_prev[0] - p_next[0]) == 1 and abs(p_prev[1] - p_next[1]) == 1:
            nx = p_prev[0] + p_next[0] - p_curr[0]
            ny = p_prev[1] + p_next[1] - p_curr[1]
            if (nx, ny) not in set(current_positions):
                new_positions[idx] = (nx, ny)
                return new_positions
        return None

    elif move_type == 'crankshaft':
        if n <= 3:
            return None
        idx = random.randint(1, n - 3)
        p_prev  = current_positions[idx - 1]
        p_curr  = current_positions[idx]
        p_next  = current_positions[idx + 1]
        p_next2 = current_positions[idx + 2]
        dist_sq = (p_prev[0] - p_next2[0])**2 + (p_prev[1] - p_next2[1])**2
        if dist_sq == 1:
            nx_curr = p_prev[0] + p_next2[0] - p_next[0]
            ny_curr = p_prev[1] + p_next2[1] - p_next[1]
            nx_next = p_prev[0] + p_next2[0] - p_curr[0]
            ny_next = p_prev[1] + p_next2[1] - p_curr[1]
            pos_set = set(current_positions)
            if (nx_curr, ny_curr) not in pos_set and (nx_next, ny_next) not in pos_set:
                new_positions[idx]     = (nx_curr, ny_curr)
                new_positions[idx + 1] = (nx_next, ny_next)
                return new_positions
        return None

    elif move_type == 'pivot':
        if n <= 2:
            return None
        pivot_idx = random.randint(1, n - 2)
        angle = random.choice([90, -90, 180])
        cx, cy = current_positions[pivot_idx]
        if angle == 90:
            cos_a, sin_a = 0, 1
        elif angle == -90:
            cos_a, sin_a = 0, -1
        else:
            cos_a, sin_a = -1, 0
        for i in range(pivot_idx + 1, n):
            x, y = current_positions[i]
            tx, ty = x - cx, y - cy
            new_positions[i] = (tx * cos_a - ty * sin_a + cx,
                                tx * sin_a + ty * cos_a + cy)
        if len(set(new_positions)) == len(new_positions):
            return new_positions
        return None

    return None


# ─────────────────────────────────────────────────────────────────────────────
# Tabular Q-Learning  –  core algorithm
# ─────────────────────────────────────────────────────────────────────────────

def generate_2d_structure_ql(
    hp_string,
    episodes=500,
    max_steps_per_episode=200,
    alpha=0.1,          # learning rate
    gamma=0.95,         # discount factor
    epsilon_start=1.0,  # initial exploration rate
    epsilon_end=0.05,   # final exploration rate
    epsilon_decay=None, # steps over which epsilon decays (default: auto)
):
    """
    Tabular Q-Learning for 2-D HP protein folding.

    The agent starts from a straight-line conformation and iteratively
    improves it by learning which moves lead to lower-energy structures.

    State  : canonical encoding of the current conformation
             (positions_to_state).
    Action : one of 4 move types (end_flip, kink_jump, crankshaft, pivot).
    Reward : ΔE = E_old − E_new  (positive when energy is reduced).

    Returns
    -------
    best_positions : list of (x, y) tuples
    best_energy    : int (negative → more HH contacts)
    q_table        : the learned Q-table (for inspection / reuse)
    """
    if not hp_string or len(hp_string) <= 1:
        return [(0, 0)] * len(hp_string), 0, {}

    n = len(hp_string)
    actions = MOVE_TYPES

    # Q-table: state → {action: Q-value}
    q_table = defaultdict(lambda: {a: 0.0 for a in actions})

    # Epsilon decay schedule
    total_steps = episodes * max_steps_per_episode
    if epsilon_decay is None:
        epsilon_decay = total_steps // 2   # decay over first half of training

    step_counter = 0

    # Track global best
    best_positions = [(i, 0) for i in range(n)]
    best_energy    = calculate_energy(best_positions, hp_string)

    for episode in range(episodes):
        # Reset to straight-line conformation each episode
        current_positions = [(i, 0) for i in range(n)]
        current_energy    = calculate_energy(current_positions, hp_string)
        current_state     = positions_to_state(current_positions)

        for _ in range(max_steps_per_episode):
            # ε-greedy action selection
            epsilon = max(
                epsilon_end,
                epsilon_start - (epsilon_start - epsilon_end) * step_counter / epsilon_decay
            )
            step_counter += 1

            if random.random() < epsilon:
                action = random.choice(actions)   # explore
            else:
                action = max(actions, key=lambda a: q_table[current_state][a])  # exploit

            # Apply the chosen move
            new_positions = apply_move(current_positions, action, hp_string)

            if new_positions is None:
                # Invalid move: small negative reward, stay in same state
                reward      = -0.1
                next_state  = current_state
                next_energy = current_energy
            else:
                next_energy = calculate_energy(new_positions, hp_string)
                # Reward = energy improvement (lower energy is better)
                reward      = float(current_energy - next_energy)
                next_state  = positions_to_state(new_positions)

            # ── Q-Learning update (off-policy, Bellman equation) ──
            best_next_q = max(q_table[next_state].values())
            old_q       = q_table[current_state][action]
            q_table[current_state][action] = (
                old_q + alpha * (reward + gamma * best_next_q - old_q)
            )

            # Transition
            if new_positions is not None:
                current_positions = new_positions
                current_energy    = next_energy
                current_state     = next_state

                # Update global best
                if current_energy < best_energy:
                    best_energy    = current_energy
                    best_positions = list(current_positions)

    return best_positions, best_energy, q_table


# ─────────────────────────────────────────────────────────────────────────────
# Greedy exploitation pass using the learned Q-table
# ─────────────────────────────────────────────────────────────────────────────

def exploit_q_table(hp_string, q_table, steps=500):
    """
    Runs a purely greedy pass (ε = 0) using the trained Q-table.
    Useful to extract the best policy found by the agent.

    Returns
    -------
    best_positions, best_energy
    """
    n = len(hp_string)
    actions = MOVE_TYPES

    current_positions = [(i, 0) for i in range(n)]
    current_energy    = calculate_energy(current_positions, hp_string)
    best_positions    = list(current_positions)
    best_energy       = current_energy

    for _ in range(steps):
        state  = positions_to_state(current_positions)
        action = max(actions, key=lambda a: q_table[state][a])

        new_positions = apply_move(current_positions, action, hp_string)
        if new_positions is None:
            continue

        new_energy = calculate_energy(new_positions, hp_string)
        current_positions = new_positions
        current_energy    = new_energy

        if current_energy < best_energy:
            best_energy    = current_energy
            best_positions = list(current_positions)

    return best_positions, best_energy


# ─────────────────────────────────────────────────────────────────────────────
# Main test – small proteins only (≤ 25 AA recommended for tabular QL)
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Small proteins well-suited for Tabular Q-Learning
    p1  = "MKVILA"                    #  6 AA
    p2  = "MKTIIALSYIFCLVF"          # 15 AA
    p3  = "MVHLTPEEKSAVTALWGK"       # 18 AA
    p4  = "MVHLTPEEKSAVTALWGKVN"     # 20 AA
    p5  = "MVHLTPEEKSAVTALWGKVNVDEVG" # 25 AA

    proteins = [
        ("Protein 1 (6 AA)",  p1),
        ("Protein 2 (15 AA)", p2),
        ("Protein 3 (18 AA)", p3),
        ("Protein 4 (20 AA)", p4),
        ("Protein 5 (25 AA)", p5),
    ]

    # Episode configurations to study convergence
    episode_configs = [
        {"episodes": 200,  "max_steps": 100},
        {"episodes": 500,  "max_steps": 200},
        {"episodes": 1000, "max_steps": 300},
        {"episodes": 2000, "max_steps": 400},
    ]

    print("=" * 65)
    print("  TABULAR Q-LEARNING  –  HP Model 2-D Protein Folding")
    print("  alpha=0.1  |  gamma=0.95  |  epsilon: 1.0 → 0.05")
    print("=" * 65)
    print()

    for name, seq in proteins:
        hp_str  = sequence_to_hp(seq)
        h_count = hp_str.count('H')
        print(f"--- {name} ---")
        print(f"  Sequence  : {seq}")
        print(f"  HP String : {hp_str}  (H: {h_count})")
        print()

        for cfg in episode_configs:
            eps  = cfg["episodes"]
            stps = cfg["max_steps"]

            start = time.time()
            positions, energy, q_table = generate_2d_structure_ql(
                hp_str,
                episodes=eps,
                max_steps_per_episode=stps,
                alpha=0.1,
                gamma=0.95,
            )
            elapsed = time.time() - start

            # Exploitation pass with the trained Q-table
            ex_positions, ex_energy = exploit_q_table(hp_str, q_table, steps=1000)
            final_energy = min(energy, ex_energy)

            total_steps = eps * stps
            print(
                f"  episodes={eps:<5} steps/ep={stps:<4} "
                f"| total_steps={total_steps:<8} "
                f"| Best E={final_energy:<4} "
                f"| Q-states={len(q_table):<6} "
                f"| Time={elapsed:.3f}s"
            )

        print()
        print("-" * 65)
        print()

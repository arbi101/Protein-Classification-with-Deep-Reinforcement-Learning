"""
Inference utilities for the Colab-trained advanced DQN model.

This model is the constructive DQN exported from colab_dqn_advanced.py. It is
not compatible with the older DQNAgent in agent_dqn.py: the advanced model uses
a 2-channel 21x21 local grid and 3 relative actions.
"""

import os
import random

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class ConstructiveHPEnv:
    def __init__(self, hp_string, grid_size=21):
        self.hp_string = hp_string
        self.n = len(hp_string)
        self.grid_size = grid_size
        self.half_grid = grid_size // 2
        self.reset()

    def reset(self):
        if self.n == 0:
            self.positions = []
            self.current_step = 0
            self.current_dir = (1, 0)
            self.done = True
            self.total_energy = 0
            return self._get_state()

        self.positions = [(0, 0)]
        if self.n > 1:
            self.positions.append((1, 0))

        self.current_step = len(self.positions)
        self.current_dir = (1, 0)
        self.done = self.current_step >= self.n
        self.total_energy = calculate_hp_energy(self.positions, self.hp_string)
        return self._get_state()

    def step(self, action):
        if self.done:
            return self._get_state(), 0.0, True

        dx, dy = self.current_dir
        if action == 0:
            new_dir = (dx, dy)
        elif action == 1:
            new_dir = (-dy, dx)
        elif action == 2:
            new_dir = (dy, -dx)
        else:
            raise ValueError(f"Unknown action index: {action}")

        new_pos = (
            self.positions[-1][0] + new_dir[0],
            self.positions[-1][1] + new_dir[1],
        )

        if new_pos in self.positions:
            self.done = True
            return self._get_state(), -1.0, True

        before_energy = self.total_energy
        self.positions.append(new_pos)
        self.current_dir = new_dir
        self.current_step += 1
        self.total_energy = calculate_hp_energy(self.positions, self.hp_string)

        reward = float(before_energy - self.total_energy)
        if self.current_step == self.n:
            self.done = True
            reward += 0.5

        return self._get_state(), reward, self.done

    def _get_state(self):
        state = np.zeros((2, self.grid_size, self.grid_size), dtype=np.float32)
        if not self.positions:
            return state

        head_x, head_y = self.positions[-1]
        dx, dy = self.current_dir

        for i, (x, y) in enumerate(self.positions):
            rel_x = x - head_x
            rel_y = y - head_y
            rot_x = rel_x * dy - rel_y * dx
            rot_y = rel_x * dx + rel_y * dy
            grid_x = rot_x + self.half_grid
            grid_y = rot_y + self.half_grid

            if 0 <= grid_x < self.grid_size and 0 <= grid_y < self.grid_size:
                channel = 0 if self.hp_string[i] == "H" else 1
                state[channel, grid_y, grid_x] = 1.0

        return state


class DQN_CNN(nn.Module):
    def __init__(self, grid_size=21):
        super().__init__()
        self.conv1 = nn.Conv2d(2, 16, kernel_size=3, stride=1, padding=1)
        self.conv2 = nn.Conv2d(16, 32, kernel_size=3, stride=1, padding=1)
        self.flat_size = 32 * grid_size * grid_size
        self.fc1 = nn.Linear(self.flat_size, 128)
        self.fc2 = nn.Linear(128, 3)

    def forward(self, x):
        x = F.relu(self.conv1(x))
        x = F.relu(self.conv2(x))
        x = x.view(x.size(0), -1)
        x = F.relu(self.fc1(x))
        return self.fc2(x)


_MODEL_CACHE = {}


def _torch_load_weights(weights_path, device):
    try:
        return torch.load(weights_path, map_location=device, weights_only=True)
    except TypeError:
        return torch.load(weights_path, map_location=device)


def load_advanced_policy(weights_path):
    abs_path = os.path.abspath(weights_path)
    try:
        stat = os.stat(abs_path)
        cache_key = (abs_path, stat.st_mtime_ns, stat.st_size)
    except OSError:
        cache_key = (abs_path, None, None)

    if cache_key in _MODEL_CACHE:
        return _MODEL_CACHE[cache_key]

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    policy_net = DQN_CNN(grid_size=21).to(device)
    state_dict = _torch_load_weights(weights_path, device)
    policy_net.load_state_dict(state_dict)
    policy_net.eval()

    stale_keys = [
        key for key in _MODEL_CACHE
        if (key[0] if isinstance(key, tuple) else key) == abs_path
    ]
    for key in stale_keys:
        del _MODEL_CACHE[key]

    _MODEL_CACHE[cache_key] = (policy_net, device)
    return policy_net, device


def calculate_hp_energy(positions, hp_string):
    energy = 0
    limit = min(len(positions), len(hp_string))

    for i in range(limit):
        if hp_string[i] != "H":
            continue
        x1, y1 = positions[i]
        for j in range(i + 2, limit):
            if hp_string[j] != "H":
                continue
            x2, y2 = positions[j]
            if abs(x1 - x2) + abs(y1 - y2) == 1:
                energy -= 1

    return energy


def straight_line_fallback(hp_string):
    positions = [(i, 0) for i in range(len(hp_string))]
    return positions, calculate_hp_energy(positions, hp_string)


def _rank_candidate(candidate):
    positions, energy = candidate[:2]
    return (-len(positions), energy)


def _run_constructive_rollout(policy_net, device, hp_string, epsilon,
                              return_trace=False, trace_interval=1,
                              max_trace_frames=200):
    env = ConstructiveHPEnv(hp_string, grid_size=21)
    state = env.reset()
    trace_frames = []

    def add_trace_frame(step):
        if not return_trace:
            return
        if len(trace_frames) >= max_trace_frames:
            return
        trace_frames.append({
            "step": step,
            "energy": calculate_hp_energy(env.positions, hp_string),
            "positions": list(env.positions),
        })

    trace_interval = max(1, int(trace_interval))
    max_trace_frames = max(2, int(max_trace_frames))
    add_trace_frame(len(env.positions))

    while not env.done:
        if random.random() < epsilon:
            action = random.randrange(3)
        else:
            with torch.no_grad():
                state_tensor = torch.as_tensor(
                    state, dtype=torch.float32, device=device
                ).unsqueeze(0)
                action = policy_net(state_tensor).argmax(dim=1).item()

        state, _, _ = env.step(action)
        if len(env.positions) % trace_interval == 0:
            add_trace_frame(len(env.positions))

    final_energy = calculate_hp_energy(env.positions, hp_string)
    if return_trace and (not trace_frames or trace_frames[-1]["step"] != len(env.positions)):
        add_trace_frame(len(env.positions))

    if return_trace:
        return list(env.positions), final_energy, trace_frames
    return list(env.positions), final_energy


def run_dqn_advanced_inference(hp_string, weights_path, max_rollouts=64,
                               return_trace=False, trace_interval=1,
                               max_trace_frames=200):
    """
    Fold an HP string using the advanced constructive DQN weights.

    max_rollouts controls how many greedy/stochastic constructions are sampled.
    The returned structure prefers complete folds first, then lower energy.
    """
    if not hp_string:
        if return_trace:
            return [], 0, []
        return [], 0

    if not os.path.exists(weights_path):
        positions, energy = straight_line_fallback(hp_string)
        if return_trace:
            trace_frames = [
                {"step": idx + 1, "energy": calculate_hp_energy(positions[:idx + 1], hp_string), "positions": positions[:idx + 1]}
                for idx in range(0, len(positions), max(1, int(trace_interval)))
            ]
            if not trace_frames or trace_frames[-1]["step"] != len(positions):
                trace_frames.append({"step": len(positions), "energy": energy, "positions": positions})
            return positions, energy, trace_frames[:max_trace_frames]
        return positions, energy

    try:
        policy_net, device = load_advanced_policy(weights_path)
    except Exception as exc:
        print(f"Failed to load advanced DQN weights: {exc}")
        positions, energy = straight_line_fallback(hp_string)
        if return_trace:
            return positions, energy, [{"step": len(positions), "energy": energy, "positions": positions}]
        return positions, energy

    rollout_count = max(1, min(int(max_rollouts), 300))
    epsilons = [0.0, 0.05, 0.10, 0.20]
    candidates = []

    for idx in range(rollout_count):
        epsilon = epsilons[idx % len(epsilons)]
        candidates.append(_run_constructive_rollout(
            policy_net,
            device,
            hp_string,
            epsilon,
            return_trace=return_trace,
            trace_interval=trace_interval,
            max_trace_frames=max_trace_frames,
        ))

    complete = [candidate for candidate in candidates if len(candidate[0]) == len(hp_string)]
    selected = min(complete, key=lambda candidate: candidate[1]) if complete else min(candidates, key=_rank_candidate)

    if return_trace:
        positions, energy, trace_frames = selected
        return positions, energy, trace_frames
    return selected

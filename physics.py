# physics.py

import numpy as np
from scipy.spatial import Delaunay
from numba import njit

class PhysicsEngine:
    def __init__(self, positions, velocities, params):
        self.positions = positions
        self.velocities = velocities
        self.params = params
        self.N = positions.shape[0]
        self.domain_size = params.domain_size

        # Initialize neighbor lists and rest lengths
        self.neighbors = self.compute_neighbors()
        self.rest_lengths = self.compute_rest_lengths()

    def compute_neighbors(self):
        # Compute Delaunay triangulation
        delaunay = Delaunay(self.positions)
        neighbors = [set() for _ in range(self.N)]
        for simplex in delaunay.simplices:
            for idx in simplex:
                others = simplex[simplex != idx]
                neighbors[idx].update(others)
        return neighbors

    def compute_rest_lengths(self):
        rest_lengths = {}
        for i, nbrs in enumerate(self.neighbors):
            for j in nbrs:
                key = tuple(sorted((i, j)))
                if key not in rest_lengths:
                    delta = self.positions[i] - self.positions[j]
                    delta = self.apply_periodic_boundary(delta)
                    rest_lengths[key] = np.linalg.norm(delta)
        return rest_lengths

    def apply_periodic_boundary(self, delta):
        L = self.domain_size
        delta = delta - np.round(delta / L) * L
        return delta

    def compute_forces(self):
        forces = np.zeros_like(self.positions)
        positions = self.positions
        velocities = self.velocities
        k = self.params.spring_constant
        b = self.params.damping_coefficient
        L = self.domain_size

        for i, nbrs in enumerate(self.neighbors):
            for j in nbrs:
                key = tuple(sorted((i, j)))
                delta_pos = positions[i] - positions[j]
                delta_pos = self.apply_periodic_boundary(delta_pos)
                delta_vel = velocities[i] - velocities[j]
                distance = np.linalg.norm(delta_pos)
                if distance == 0:
                    continue  # Avoid division by zero
                direction = delta_pos / distance
                rest_length = self.rest_lengths[key]
                # Spring force
                spring_force = -k * (distance - rest_length) * direction
                # Damping force
                damping_force = -b * delta_vel
                # Total force
                forces[i] += spring_force + damping_force
        return forces

    def apply_brownian_motion(self, forces):
        random_forces = np.random.normal(scale=self.params.noise_scale, size=forces.shape)
        forces += random_forces
        return forces


    def update(self):
        # Compute forces
        forces = self.compute_forces()
        # Apply Brownian motion
        forces = self.apply_brownian_motion(forces)
        # Update velocities and positions using Velocity-Verlet integration
        dt = self.params.time_step
        self.velocities += forces * dt
        self.positions += self.velocities * dt
        # Enforce reflecting boundaries
        self.enforce_reflecting_boundaries()
        # Recompute neighbors and rest lengths if necessary
        if self.params.dynamic_neighbors:
            self.neighbors = self.compute_neighbors()
            self.rest_lengths = self.compute_rest_lengths()

    def enforce_reflecting_boundaries(self):
        for i in range(self.N):
            for dim in range(2):  # x and y dimensions
                if self.positions[i, dim] < 0:
                    self.positions[i, dim] = -self.positions[i, dim]
                    self.velocities[i, dim] = -self.velocities[i, dim] * self.params.boundary_damping
                elif self.positions[i, dim] > self.domain_size:
                    self.positions[i, dim] = 2 * self.domain_size - self.positions[i, dim]
                    self.velocities[i, dim] = -self.velocities[i, dim] * self.params.boundary_damping

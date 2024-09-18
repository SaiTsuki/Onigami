# parameters.py

class Parameters:
    def __init__(self):
        self.N = 8
        self.domain_size = 1.0  # Domain size (0 to domain_size in both x and y)
        self.spring_constant = 1.0
        self.damping_coefficient = 0.1
        self.noise_scale = 0.01
        self.time_step = 0.01
        self.dynamic_neighbors = False  # Set to True if neighbors are updated dynamically
        self.boundary_damping = 1.0  # Coefficient of restitution for boundary collisions

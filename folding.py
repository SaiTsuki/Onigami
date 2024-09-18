# folding.py

import numpy as np
from shapely.geometry import Point, Polygon
from scipy.spatial import Voronoi

class FoldingMechanism:
    def __init__(self, params):
        self.params = params
        self.folding_states = None  # To be initialized after first Voronoi computation
        self.voronoi = None

    def update_voronoi(self, positions):
        self.voronoi = Voronoi(positions)
        if self.folding_states is None:
            self.folding_states = [False] * positions.shape[0]

    def update_folding_states(self, positions):
        for i, pos in enumerate(positions):
            region_index = self.voronoi.point_region[i]
            region = self.voronoi.regions[region_index]
            if -1 in region or len(region) == 0:
                continue  # Skip regions with infinite vertices
            vertices = self.voronoi.vertices[region]
            cell_polygon = Polygon(vertices)
            nucleus_point = Point(pos)
            if not cell_polygon.contains(nucleus_point):
                self.folding_states[i] = not self.folding_states[i]  # Toggle state

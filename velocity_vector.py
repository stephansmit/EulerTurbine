import numpy as np

class VelocityVector():
    def __init__(self, omega=0, radius=0):
        self.r = 0
        self.theta = omega * radius
        self.vec = np.array([self.r, self.theta])

    @classmethod
    def blade_velocity(cls, omega, radius):
        return cls( omega, radius)

    def set_components_with_vector(self):
        self.r = self.vec[0]
        self.theta = self.vec[1]

    def set_vector_with_components(self):
        self.vec = np.array([self.r, self.theta])

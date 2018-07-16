import numpy as np

class VelocityVector():
    def __init__(self, omega=0, radius=0):
        self.r = 0
        self.theta = omega * radius
        self.vec = np.array([self.r, self.theta])
        self.mag = np.linalg.norm(self.vec)
        self.angle = 0
    @classmethod
    def blade_velocity(cls, omega, radius):
        return cls( omega, radius)

    def set_components_with_vector(self):
        self.r = self.vec[0]
        self.theta = self.vec[1]
        self.mag = np.sqrt(self.vec[0]**2 + self.vec[1]**2)
        self.angle = np.degrees(np.arctan(self.theta / self.r))

    def set_vector_with_components(self):
        self.vec = np.array([self.r, self.theta])
        self.mag = np.linalg.norm(self.vec)
        self.angle = np.degrees(np.arctan(self.theta/self.r))

    def set_vector_with_thetacomponent_mag(self, theta, mag):
        self.theta = theta
        self.r = np.sqrt(mag**2 - self.theta**2 )
        self.set_vector_with_components()

    def set_vector_with_rcomponent_mag(self, r, mag):
        self.r = r
        self.theta = -np.sqrt(mag ** 2 - self.r ** 2)
        self.set_vector_with_components()

    def get_velocity_info(self):
        return dict({'r': self.r, 'theta': self.theta, 'mag': self.mag})

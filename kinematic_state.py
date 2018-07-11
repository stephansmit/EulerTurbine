from velocity_vector import VelocityVector
import numpy as np

class KinematicState():
    def __init__(self, omega=0, r=0):
        self.w = VelocityVector()
        self.c = VelocityVector()
        self.u = VelocityVector(omega, r)
        self.alpha = 0
        self.beta = 0
        self.r = r
        self.omega = omega

    def set_state_alpha_cmag(self, alpha, cmag):
        self.c.r = cmag * np.cos(np.radians(alpha))
        self.c.theta = cmag * np.sin(np.radians(alpha))
        self.c.vec = np.append(self.c.r, self.c.theta)
        self.set_w_with_c_u()

    def set_w_with_c_u(self):
        self.w.vec = np.subtract(self.c.vec, self.u.vec)
        self.w.set_components_with_vector()

    def set_state_wtheta_cr(self, wtheta, cr):
        self.w.theta = wtheta
        self.w.r = cr
        self.w.set_vector_with_components()
        self.set_c_with_w_u()

    def set_state_beta_wr(self, beta, wr):
        self.beta = beta
        self.w.r = wr
        self.w.theta = self.w.r * np.tan(np.radians(beta))
        self.w.set_vector_with_components()
        self.set_c_with_w_u()

    def set_c_with_w_u(self):
        self.c.vec = np.add(self.w.vec, self.u.vec)
        self.c.set_components_with_vector()

    def set_mach_numbers(self, A):
        self.c.mach = np.linalg.norm(self.c.vec)/A
        self.w.mach = np.linalg.norm(self.w.vec)/A

    #functions for plotting
    def draw_with_c(self, ax, origin):
        ax.quiver(origin[0], origin[1],
                  self.c.r, self.c.theta,
                  color='r', scale=1, angles='xy', scale_units="xy")
        ax.quiver(origin[0] + self.c.r, origin[1] + self.c.theta,
                  -self.u.r, -self.u.theta,
                  color='g', scale=1, angles='xy', scale_units="xy")
        ax.quiver(origin[0], origin[1],
                  self.w.r, self.w.theta,
                  color='b', scale=1, angles='xy', scale_units="xy")

    def draw_with_w(self, ax, origin):
        ax.quiver(origin[0], origin[1],
                  self.w.r, self.w.theta,
                  color='b', scale=1, angles='xy', scale_units="xy")
        ax.quiver(origin[0] + self.w.r, origin[1] + self.w.theta,
                  self.u.r, self.u.theta,
                  color='g', scale=1, angles='xy', scale_units="xy")
        ax.quiver(origin[0], origin[1],
                  self.c.r, self.c.theta,
                  color='r', scale=1, angles='xy', scale_units="xy")

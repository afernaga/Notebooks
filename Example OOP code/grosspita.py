from typing import Any

import numpy as np
from tqdm import tqdm
from typing_extensions import SupportsIndex


class GrossPitaevskiiProblem:
    """Class that generates a Gross-Pitaevskii problem instance,
    in harmonic oscillator units, for a 3D harmonic trap.
    """

    def __init__(
        self,
        particle_number,
        grid_size,
        grid_step,
        scattering_length,
        sigma,
        time_step,
        iterations,
        interacting_system=True,
        thomas_fermi=False,
    ):
        # Put default values in the declaration above
        self.particle_number = particle_number
        self.grid_size = grid_size
        self.grid_step = grid_step
        self.scattering_length = scattering_length
        self.sigma = sigma
        self.time_step = time_step
        self.iterations = iterations
        self.interacting_system = interacting_system
        self.thomas_fermi = thomas_fermi
        self.discreted_r = np.arange(0, self.grid_size, self.grid_step)
        self.evolved_psi = np.zeros(len(self.discreted_r))
        self.evolved_ddpsi = np.zeros(len(self.discreted_r))
        self.mu = np.zeros(len(self.discreted_r))
        self.interaction = (
            self.scattering_length * self.particle_number
            if self.interacting_system
            else 0
        )
        self._has_evolved = False

    def evolution(self):
        psi = self.ansatz()
        mu = np.zeros(len(self.discreted_r))
        psi0 = psi
        for _ in tqdm(range(self.iterations)):
            normalization = self.simpson_integral(psi0 ** 2, self.grid_step)
            psi0 = psi0 / np.sqrt(normalization)
            is_normal = self.simpson_integral(psi0 ** 2, self.grid_step)
            if abs(is_normal) > 1.1:
                raise ValueError(
                    f"Wave function is not normalized, Normalization: {is_normal}"
                )
            ddpsi = self.second_derivative(psi0, self.grid_step)
            mu = self._calculate_mu(
                r_vector=self.discreted_r,
                psi=psi0,
                dpsi=ddpsi,
                interaction=self.interaction,
            )

            for j, psi_term in enumerate(psi0):
                if j == 0:
                    continue
                psi[j] = psi_term - self.time_step * mu[j] * psi_term
            psi0 = psi
        normalization = self.simpson_integral(psi0 ** 2, self.grid_step)
        psi0 = psi0 / np.sqrt(normalization)
        self.evolved_psi = psi0
        self.evolved_ddpsi = self.second_derivative(self.evolved_psi, self.grid_step)
        self.mu = mu
        self._has_evolved = True
        return self.simpson_integral(mu * psi0 ** 2, self.grid_step)

    def _calculate_mu(self, r_vector, psi, dpsi, interaction):
        """Method to calculate the chemical potential

        Args:
            r_vector (array like): values of the r-vector
            psi (array like): values of the wave function
            dpsi (array like): values of the second derivative of the wave function

        Returns:
            float: chemical potential
        """
        mu = np.zeros(len(r_vector))
        kinetic_coefficient = 0.000000000000001 if self.thomas_fermi else 0.5
        for i, x in enumerate(r_vector):
            if i == 0:
                mu[i] = 0
                continue
            mu[i] = (
                -kinetic_coefficient * dpsi[i] / psi[i]
                + 0.5 * x ** 2
                + interaction * (psi[i] / r_vector[i]) ** 2
            )
        return mu

    @staticmethod
    def second_derivative(f, h):
        derivative = np.zeros(len(f))
        for i, _ in enumerate(f):
            if i == 0:
                derivative[i] = 0
            elif i == len(f) - 1:
                derivative[i] = (f[i - 1] - 2 * f[i]) / (h ** 2)
            else:
                derivative[i] = (f[i + 1] - 2 * f[i] + f[i - 1]) / (h ** 2)
        return derivative

    @staticmethod
    def simpson_integral(f, h):
        """Method to calculate integrals using Simpson's rule

        Args:
            f (array like): values of the function to integrate, separated by h
            h (float): step size

        Returns:
            float: integral value
        """
        n = len(f) - 1
        integral = sum(2 * f[i] if i % 2 == 0 else 4 * f[i] for i in range(1, n))
        integral += f[0] + f[n]
        return integral * h / 3

    @staticmethod
    def trapezoidal_integral(f, h):
        """Method to calculate integrals using the trapezoidal rule

        Args:
            f (array like): values of the function to integrate, separated by h
            h (float): step size

        Returns:
            float: integral value
        """
        n = len(f)
        integral = sum(2 * f[i] for i in range(n))
        return integral * h

    def ansatz(self):
        r = np.arange(0, self.grid_size, self.grid_step)
        psi = np.zeros(len(r))
        cvar = 2 * np.sqrt(self.sigma) ** 3 / np.sqrt(np.sqrt(np.pi))
        for i, x in enumerate(r):
            psi[i] = cvar * x * np.exp(-0.5 * self.sigma ** 2 * x ** 2)
        return psi

    @property
    def kinetic_term(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        kinetic_coefficient = 0.000000000000001 if self.thomas_fermi else 0.5
        return -kinetic_coefficient * self.simpson_integral(
            self.evolved_ddpsi[1:] * self.evolved_psi[1:], self.grid_step
        )

    @property
    def trap_term(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        return 0.5 * self.simpson_integral(
            self.discreted_r[1:] ** 2 * self.evolved_psi[1:] ** 2, self.grid_step
        )

    @property
    def interaction_term(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        return (
            0.5
            * self.interaction
            * self.simpson_integral(
                self.evolved_psi[1:] ** 4 / self.discreted_r[1:] ** 2, self.grid_step
            )
        )

    @property
    def potential_term(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        return self.trap_term + self.interaction_term

    @property
    def energy(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        return self.kinetic_term + self.potential_term

    @property
    def density(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        return (self.evolved_psi[1:] / self.discreted_r[1:]) ** 2 * (1 / (4 * np.pi))

    @property
    def virial(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        return 2 * self.kinetic_term - 2 * self.trap_term + 3 * self.interaction_term

    @property
    def radius(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        return np.sqrt(
            self.simpson_integral(
                self.evolved_psi[1:] ** 2 * self.discreted_r[1:] ** 2, self.grid_step
            )
        )

    def check_density_normalization(self):
        if not self._has_evolved:
            raise ValueError("System has not evolved yet")
        return self.simpson_integral(
            self.density * 4 * np.pi * self.discreted_r[1:] ** 2, self.grid_step
        )

    def __str__(self) -> str:
        return f"{self.particle_number} Bosons in a spherical trap \n r-grid in {self.grid_size} points, r-step {self.grid_step} \n A0={self.scattering_length}, sigma={self.sigma} \n time={self.time_step}, number-iter={self.iterations}"

if __name__ == "__main__":
    # Create a problem instance
    problem = GrossPitaevskiiProblem(100, 100, 0.1, 1, 1, 0.1, 100)
    print(problem)

"""Range-energy relation for an alpha in He:CO2 90/10.

The table is a list of (distance_mm, energy_MeV) for an alpha that starts at 11.40 MeV
and is tabulated at 150 torr.  Ranges scale inversely with density, hence with pressure,
so a length measured at P torr is converted to the table's scale by multiplying by
P/150 before looking it up.  Energies are unaffected by that scaling.
"""
import numpy as np


class RangeTable:
    def __init__(self, path, pressure, p_table=150.0):
        d, e = [], []
        for line in open(path):
            if line.startswith("#") or not line.strip():
                continue
            a, b = line.split()
            d.append(float(a)); e.append(float(b))
        self.d = np.asarray(d)          # distance travelled, mm at p_table, increasing
        self.e = np.asarray(e)          # energy remaining, MeV, decreasing
        self.scale = pressure / p_table
        self.total_mm = self.d[-1] / self.scale      # full range at the real pressure

    def energy_of_range(self, R_mm):
        """Energy of an alpha whose REMAINING range is R_mm (i.e. a track of that length)."""
        R = np.asarray(R_mm, dtype=float)
        depth = self.d[-1] - np.clip(R * self.scale, 0, self.d[-1])
        return np.interp(depth, self.d, self.e)

    def energy_after(self, path_mm):
        """Energy of the beam after travelling path_mm from the entrance window."""
        p = np.asarray(path_mm, dtype=float)
        return np.interp(p * self.scale, self.d, self.e)

    def range_of_energy(self, E_MeV):
        """Range in mm, at the real pressure, of an alpha of energy E."""
        E = np.asarray(E_MeV, dtype=float)
        depth = np.interp(E, self.e[::-1], self.d[::-1])     # both must increase
        return (self.d[-1] - depth) / self.scale

    def depth_of_energy(self, E_MeV):
        """How far the beam must travel from the window to be left with energy E."""
        E = np.asarray(E_MeV, dtype=float)
        return np.interp(E, self.e[::-1], self.d[::-1]) / self.scale

    def dEdx(self, E_MeV):
        """Stopping power in MeV/mm at the real pressure."""
        g = -np.gradient(self.e, self.d) * self.scale        # MeV/mm at real pressure
        return np.interp(np.asarray(E_MeV, dtype=float), self.e[::-1], g[::-1])

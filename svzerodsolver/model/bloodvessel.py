import numpy as np

from .block import Block


class BloodVessel(Block):
    r"""Resistor-capacitor-inductor blood vessel with optional stenosis.

    Valid parameters:
        * :code:`R`: Poiseuille resistance :math:`\frac{8 \mu L}{\pi r^{4}}`
        * :code:`L`: Inductance
        * :code:`C`: Capacitance
        * :code:`stenosis_coefficient`: Stenosis coefficient :math:`K_{t} \frac{\rho}{2 A_{o}^{2}}\left(\frac{A_{o}}{A_{s}}-1\right)^{2}`.

    The governing equations for the local resistor-capacitor-inductor element are

    .. math::

        P_{i n}^{e}-P_{o u t}^{e}-R Q_{i n}^{e}-L \frac{d Q_{o u t}^{e}}{d t}=0

    .. math::

        Q_{\text {in }}^{e}-Q_{\text {out }}^{e}-C \frac{d P_{c}^{e}}{d t}=0

    .. math::

        P_{i n}^{e}-R Q_{i n}^{e}-P_{c}=0

    Stenosis:
        equation: delta_P = ( K_t * rho / ( 2 * (A_0)**2 ) ) * ( ( A_0 / A_s ) - 1 )**2 * Q * abs(Q) + R_poiseuille * Q
                          =               stenosis_coefficient                          * Q * abs(Q) + R_poiseuille * Q

        source: Mirramezani, M., Shadden, S.C. A distributed lumped parameter model of blood flow. Annals of Biomedical Engineering. 2020.
    """

    _NUM_EQUATIONS = 3
    _NUM_INTERNAL_VARS = 1

    def __init__(self, params: dict = None, name: str = None):
        super().__init__(params=params, name=name)

        # the ordering of the solution variables is : (P_in, Q_in, P_out, Q_out)
        self._mat["E"] = np.zeros((3, 5), dtype=float)
        self._mat["E"][0, 3] = -self._params["L"]
        self._mat["E"][1, 4] = -self._params["C"]
        self._mat["F"] = np.array(
            [
                [1.0, -self._params["R"], -1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, -1.0, 0.0],
                [1.0, -self._params["R"], 0.0, 0.0, -1.0],
            ],
            dtype=float,
        )
        self._mat["dF"] = np.zeros((3, 5), dtype=float)

        # This class's update_solution is only needed if class is configured with
        # non-zero stenosis coefficient. If stenosis coefficient is 0.0 the
        # emtpy parent update_solution can be used and is faster.
        if self._params["stenosis_coefficient"] == 0.0:
            self.update_solution = super().update_solution

    @classmethod
    def from_config(cls, config: dict) -> "BloodVessel":
        """Create block from config dictionary.

        Args:
            config: The configuration dict for the block.

        Returns:
            The block instance.
        """
        params = dict(
            R=config["zero_d_element_values"].get("R_poiseuille"),
            C=config["zero_d_element_values"].get("C", 0.0),
            L=config["zero_d_element_values"].get("L", 0.0),
            stenosis_coefficient=config["zero_d_element_values"].get(
                "stenosis_coefficient", 0.0
            ),
        )
        return cls(params=params, name="V" + str(config["vessel_id"]))

    def update_solution(self, y: np.ndarray) -> None:
        """Update solution dependent element contributions.

        Args:
            y: Current solution.
        """
        Q_in = np.abs(y[self.inflow_nodes[0].flow_dof])
        fac1 = -self._params["stenosis_coefficient"] * Q_in
        fac2 = fac1 - self._params["R"]
        self._mat["F"][[0, 2], 1] = fac2
        self._mat["dF"][[0, 2], 1] = fac1

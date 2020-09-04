import numpy as np
from typing import List, Sized, Optional, Sequence, Container
from decay_chain import DecayChain


class DecayChainSimulation(DecayChain):
    istate: List[int]
    state: List[int]

    def __init__(self, λ: List[float], μ: List[str], istate: Optional[List[int]] = None, n: int = 0):
        super().__init__(λ, μ)
        if istate is None:
            self.istate = [0] * len(λ)
            self.istate[0] = n
        else:
            assert len(istate) == len(λ)
            self.istate = istate
        self.state = self.istate

    def setup_simulation(self, n: Optional[int] = None, istate: Optional[List[int]] = None):
        if istate is None:
            assert n is not None
            self.istate = [0] * len(self)
            self.istate[0] = n
        else:
            assert len(istate) == len(self)
            self.istate = istate
        self.state = self.istate

    def reset_simulation(self):
        self.state = self.istate

    def simulate_count(self, dt: float, counted: Optional[List[str]] = None) -> int:
        if counted is None:
            counted = []
        total_decays = np.array([0] * len(self))
        rng = np.random.default_rng()
        for i in range(len(self)):
            if self.state[i] > 0:
                decays = rng.exponential([1./r for r in self.λ[i:]], (self.state[i], len(self)-i))
                total_decays[i:] += np.sum(np.cumsum(decays, 1) < dt, 0)
        self.state = [self.state[i] - total_decays[i] + (int(i>0) and total_decays[i-1]) for i in range(len(self))]
        return sum([n for i, n in enumerate(total_decays) if self.μ[i] in counted])

    def simulate(self, dt: float) -> None:
        self.simulate_count(dt)

    def simulate_counts(self, ts: float, p: int, tf: float = 0., counted: Optional[List[str]] = None) -> List[int]:
        if counted is None:
            counted = []
        sc = [0] * p
        self.simulate(tf)
        for i in range(p):
            sc[i] = self.simulate_count(ts, counted)
        return sc

    def get_state(self) -> List[int]:
        return self.state

import numpy as np
from typing import List, Sized, Optional, Sequence, Container


class DecayChain(Sized):
    λ: List[float]
    μ: List[str]

    def __init__(self, λ: List[float], μ: List[str]):
        assert (len(λ) > 0) and (len(λ) == len(μ))
        self.λ = λ
        self.μ = μ

    def __len__(self):
        return len(self.λ)

    def expected_count(self, s1: float, s2: float, counted: Optional[List[str]] = None) -> float:
        assert s1 <= s2
        if s1 == s2:
            return 0.
        if counted is None:
            counted = []
        ec = 0
        for j in range(len(self)):
            if self.μ[j] in counted:
                for r in range(j+1):
                    tmp = np.exp(-self.λ[r] * s1) - np.exp(-self.λ[r] * s2)
                    for q in range(j+1):
                        if q != r:
                            tmp *= self.λ[q] / (self.λ[q] - self.λ[r])
                    ec += tmp
        return ec

    def expected_counts(self, ts: Optional[float] = None, p: Optional[int] = None, tf: Optional[float] = None,
                        s: Optional[Sequence[float]] = None, counted: Optional[Container[str]] = None) -> List[float]:
        if s is None:
            assert (ts is not None) and (p is not None)
            if tf is None:
                tf = 0.
            s = [tf + ts * i for i in range(p + 1)]
        ec = [0.] * (len(s)-1)
        for i in range(len(s) - 1):
            ec[i] = self.expected_count(s[i],s[i+1],counted)
        return ec

    def expected_remaining(self, t: float) -> List[float]:
        er = [0.] * len(self)
        for j in range(len(self)):
            for r in range(j+1):
                tmp = self.λ[r] * np.exp(-self.λ[r] * t)
                for q in range(j+1):
                    if q != r:
                        tmp *= self.λ[q] / (self.λ[q] - self.λ[r])
                er[j] += tmp
            er[j] /= self.λ[j]
        return er



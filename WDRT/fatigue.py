# Copyright 2016 Sandia Corporation and the National Renewable Energy
# Laboratory
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
from pylab import find


def EqLoad(F, N, m):
    '''
    Calculates an equivalent fatigue load based on Miner’s Rule,
    using the rainflow counting method presented in,
    Downing & Socie, 1982, “Simple rainflow counting algorithms”

    Parameters
    ----------
        F : (np.array)
            Force or stress history
        N : (float)
            Number of lifetime cycles
        m : (float)
            S-N curve slope

    Returns
    -------
        Feq : (float)
            Equivalent force or stress
    '''

    # Extrema
    dF = np.diff(F)
    A = F[np.append(np.append(
        0, np.array(np.nonzero((dF[0:-1] * dF[1:]) < 0.0)).ravel() + 1),
        F.size - 1)]

    # Reorder
    Imax = A.argmax(0)
    if A[0] == A[A.size - 1]:
        A = np.append(A[Imax:], A[1:Imax + 1])
    else:
        A = np.append(A[Imax:], A[0:Imax + 1])

    # Rainfow Count
    i = -1
    j = -1
    E = A * 0
    R = A * 0
    M = A * 0
    while A.any():
        j = j + 1
        E[j] = A[0]
        A = np.delete(A, 0)
        while j >= 2:
            X = abs(E[j] - E[j - 1])
            Y = abs(E[j - 1] - E[j - 2])
            if X < Y:
                break
            i = i + 1
            R[i] = Y
            M[i] = (E[j - 1] + E[j - 2]) / 2
            j = j - 2
            E[j] = E[j + 2]

    # Equivalent Load
    Feq = (sum(R[0:i + 1]**m) / N)**(1 / m)

    return Feq

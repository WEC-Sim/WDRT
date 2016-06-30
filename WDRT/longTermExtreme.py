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
import scipy.interpolate as interp
import scipy.optimize as optim


'''WDRT longTermExtreme module

This module contains tools for analysis of long-term extremes.

'''


class fullLongTermSurvival():
    '''Class for full sea state long-term approach
    '''

    def __init__(self, Fr, fs):
        '''
        Parameters
        ----------
            Fr : list
                List where each element is a CCDF for a given sea state.
            fs : list
                Weighting for each sea state.
        '''
        self.Fr = Fr
        self.fs = fs

    def __call__(self, x):
        '''Survival function (S(x) = 1 - F(x))

        Parameters
        ----------
            x : np array
                Points at which to evaluate survival function.

        Returns
        -------
            S : np array
                Survival function value(s) at x.
        '''
        S = 0
        for jj in range(len(self.Fr)):
            S += self.Fr[jj](x) * self.fs[jj]
        return S

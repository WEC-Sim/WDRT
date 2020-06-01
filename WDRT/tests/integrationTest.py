import unittest
from os.path import abspath, dirname, join, isfile
import os
import numpy as np
import pandas as pd
import WDRT.fatigue as fatigue
import json, codecs

class TestFatigue(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.diameter = 1
        self.height = 2
        self.width = 3
        self.diameters = [1,2,3,4]

    @classmethod
    def tearDownClass(self):
        pass

    def test_fatigue(self):
        # Reduced size joint probability distribution
        # Wave energy periods
        Te = [[6.0, 10.0, 14.0], 
              [6.0, 10.0, 14.0], 
              [6.0, 10.0, 14.0]]
        # Significant wave heights
        Hs = [[1.25, 1.25, 1.25], 
              [2.75, 2.75, 2.75], 
              [4.25, 4.25, 4.25]]
        # Probability
        P = np.array([[0.2345, 0.2478, 0.0164],
                      [0.0918, 0.2821, 0.0411],
                      [0.0005, 0.05  , 0.0234]])
        [h, t] = np.shape(P)
        N1Hour=0
        for i in range(h):
            for j in range(t):
                # Average N in 1 hour (Tavg = 0.82476*Te)
                N1Hour += 1 * 60 * 60 * P[i][j] / (Te[i][j] * .82476)
        N1Year = N1Hour *24*365
        
        # Assume an S-N curve slope (m) of 6 (representative of cast iron)
        m = 6.
        Fequivalent1Hour = np.zeros((h, t))
        Fequivalent1Year = 0
       
        # Read pre-calculated PTO force histories for each sea state
        timeseriesForcePTO="FPTO.json"
        #json.dump(FPTO,
        #          codecs.open(file_path, 'w', encoding='utf-8'),
        #          separators=(',', ':'), sort_keys=True, indent=4)
        obj_text = codecs.open(timeseriesForcePTO, 'r', encoding='utf-8').read()
        FPTO = json.loads(obj_text)


        for i in range(h):
            for j in range(t):
                Hval=str(int(Hs[i][j]))
                Tval= str(int(Te[i][j]))
                dictKey=f'H{Hval}T{Tval}'
                Fpto = np.array(FPTO[dictKey])
                # Equivalent fatigue load for a 1 hour timeframe
                Fequivalent1Hour[i][j] = fatigue.EqLoad(Fpto, N1Hour, m)
                Fequivalent1Year += (Fequivalent1Hour[i][j]**m) * N1Hour * P[i][j]


        # Equivalent fatigue load for a 1 year timeframe
        Fequivalent1Year = (Fequivalent1Year / N1Year)**(1 / m)
        
        print('1 hour equivalent fatigue loads: (in Newtons)')
        print(Fequivalent1Hour)
        print('1 year equivalent fatigue load: (in Newtons)')
        print(Fequivalent1Year)



if __name__ == '__main__':
      unittest.main() 


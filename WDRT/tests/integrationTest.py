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
        pass
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
        FEquivalent1Hour = np.zeros((h, t))
        FEquivalent1Year = 0
       
        # Read pre-calculated PTO force histories for each sea state
        timeseriesForcePTO=os.path.join("data","FPTO.json")
        forcesFile = codecs.open(timeseriesForcePTO, 'r', encoding='utf-8')
        forcesFileData = forcesFile.read()
        FPTO = json.loads(forcesFileData)
        forcesFile.close()

        for i in range(h):
            for j in range(t):
                Hval=str(int(Hs[i][j]))
                Tval= str(int(Te[i][j]))
                dictKey=f'H{Hval}T{Tval}'
                Fpto = np.array(FPTO[dictKey])
                # Equivalent fatigue load for a 1 hour timeframe
                FEquivalent1Hour[i][j] = fatigue.EqLoad(Fpto, N1Hour, m)
                FEquivalent1Year += (FEquivalent1Hour[i][j]**m) * N1Hour * P[i][j]


        # Equivalent fatigue load for a 1 year timeframe
        FEquivalent1Year = (FEquivalent1Year / N1Year)**(1 / m)
       

        expected1HourSolution=[[ 846265.91946806, 2126737.82456905, 2129283.43688851],
                               [1974807.39144152, 4315273.86105664, 5138299.64239418],
                               [2961883.66545167, 5930901.74577105, 8090902.26582692]]
        expected1YearSolution=1044099.6832087268

        self.assertTrue(np.allclose(FEquivalent1Hour, np.array(expected1HourSolution)))
        self.assertEqual(FEquivalent1Year, expected1YearSolution)



if __name__ == '__main__':
      unittest.main() 


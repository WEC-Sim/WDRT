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

import requests
import bs4
import urllib2
import re
import numpy as np
from datetime import datetime, date
import os


def fetchFromWeb(buoyNum, savePath='./'):
    '''Searches ndbc.noaa.gov for the historical spectral wave density
    data of a given device and writes the annual files from the website
    to a single .txt file, and stores the values in an array

    Parameters
    ----------
    buoyNum : string
        Buoy number to download data for. Can find the appropriate buoy at
        http://www.ndbc.noaa.gov. Only buoys with historical data (i.e. have a
        http://www.ndbc.noaa.gov/station_history.php?station=##### page) can be
        used for environmental contour analyses.
    savePath : string
        Relative path to place directory with data files.

    Returns
    -------

    swdList : list
        List that contains numpy arrays of the spectral wave density data,
        separated by year.
    freqList: list
        List that contains numpy arrays that contain the frequency values
        for each year
    dateList : list
        List that contains numpy arrays of the date values for each line of
        spectral data, separated by year.

    Example
    -------
    To download data from the NDBC website for the first time

    >>> import NDBCdata
    >>> swdList, freqList, dateList = NDBCdata.fetchFromWeb(46022)
    '''
    numLines = 0
    numCols = 0
    numDates = 0
    dateVals = []
    spectralVals = []
    dateList = []
    freqList = []
    swdList = []
    url = "http://www.ndbc.noaa.gov/station_history.php?station=%s" % (buoyNum)
    ndbcURL = requests.get(url)
    ndbcURL.raise_for_status()
    ndbcHTML = bs4.BeautifulSoup(ndbcURL.text, "lxml")
    b = ndbcHTML.findAll("b", text="Spectral wave density data: ")

    if len(b) == 2:
        b = b[1]
    else:
        b = b[0]

    links = [a["href"] for a in b.find_next_siblings("a", href=True)]
    # Grab the device number so the filename is more specific

    saveDir = os.path.join(savePath, 'NDBC%s' % (buoyNum))

    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    for link in links:
        dataLink = "http://ndbc.noaa.gov" + link
        year = int(re.findall("[0-9]+", link)[1])

        #certain years have multiple files marked with the letter 'b'
        if ('b' + str(year)) not in link:
            swdFile = open(os.path.join(saveDir, "SWD-%s-%d.txt" %
                           (buoyNum, year)), 'w')
        else:
            swdFile = open(os.path.join(saveDir, "SWD-%s-%s.txt" %
                           (buoyNum, str(year) + 'b')), 'w')

        fileName = dataLink.replace('download_data', 'view_text_file')
        data = urllib2.urlopen(fileName)
        print "Reading from:", fileName



        # dates after 2004 contain a time-value for minutes
        if (year > 2004):
            numDates = 5
        else:
            numDates = 4

        #First Line of every file contains the frequency data
        frequency = data.readline()
        swdFile.write(frequency)
        frequency = np.array(frequency.split()[numDates:], dtype = np.float)


        for line in data:
            swdFile.write(line)
            currentLine = line.split()
            numCols = len(currentLine)

            if float(currentLine[numDates+1]) < 999:
                numLines += 1
                for j in range(numDates):
                    dateVals.append(currentLine[j])
                for j in range(numCols - numDates):
                    spectralVals.append(currentLine[j + numDates])

        dateValues = np.array(dateVals, dtype=np.int)
        spectralValues = np.array(spectralVals, dtype=np.float)
        dateValues = np.reshape(dateValues, (numLines, numDates))
        spectralValues = np.reshape(spectralValues, (numLines,
                                                     (numCols - numDates)))
        del dateVals[:]
        del spectralVals[:]

        numLines = 0
        numCols = 0
        dateList.append(dateValues)
        swdList.append(spectralValues)
        freqList.append(frequency)

    swdFile.close()
    return swdList, freqList, dateList


def loadFromText(dirPath):
    '''Loads NDBC data previously downloaded to a series of text files in the
    specified directory.

    Parameters
    ----------
        dirPath : string
            Relative path to directory containing NDBC text files (created by
            NBDCdata.fetchFromWeb)

    Returns
    -------
        swdList : list
            List that contains numpy arrays of the spectral wave density data,
            separated by year.
        freqList: list
            List that contains numpy arrays that contain the frequency values
            for each year
        dateList : list
            List that contains numpy arrays of the date values for each line of
            spectral data, separated by year.

    Example
    -------
    To load data from previously downloaded files

    >>> import NDBCdata
    >>> swdList, freqList, dateList = NDBCdata.loadFromText('./NDBC460022')
    '''
    dateList = []
    swdList = []
    freqList = []
    dateVals = []
    spectralVals = []
    numLines = 0

    for fileName in os.listdir(dirPath):
        fn = os.path.join(dirPath, fileName)
        print 'Reading from: %s' % (fn)
        f = open(fn, 'r')
        frequency = f.readline().split()
        numCols = len(frequency)
        if frequency[4] == 'mm':
            frequency = np.array(frequency[5:], dtype=np.float)
            numTimeVals = 5

        else:
            frequency = np.array(frequency[4:], dtype=np.float)
            numTimeVals = 4

        for line in f:
            currentLine = line.split()
            if float(currentLine[numTimeVals + 1]) < 999:
                numLines += 1
                for i in range(numTimeVals):
                    dateVals.append(currentLine[i])
                for i in range(numCols - numTimeVals):
                    spectralVals.append(currentLine[i + numTimeVals])

        dateValues = np.array(dateVals, dtype=np.int)
        spectralValues = np.array(spectralVals, dtype=np.double)
        dateValues = np.reshape(dateValues, (numLines, numTimeVals))
        spectralValues = np.reshape(
            spectralValues, (numLines, (numCols - numTimeVals)))

        del dateVals[:]
        del spectralVals[:]

        numLines = 0
        numCols = 0
        dateList.append(dateValues)
        swdList.append(spectralValues)
        freqList.append(frequency)

    return swdList, freqList, dateList


def prepData(swdList, freqList, dateList):
    '''Runs getStats and getDataNums for full set of data, then removes any
    NaNs.

    Parameters
    ----------
        swdList : list
            List that contains numpy arrays of the spectral wave density data,
            separated by year.
        freqList: list
            List that contains numpy arrays that contain the frequency values
            for each year.
        dateList : list
            List of date vals from NDBC.fetchFromWeb

    Returns
    -------
        Hm0 : list
            Significant wave height.
        Te : list
            Energy period.
        dateNum : list
            List of datetime objects.

    Example
    -------
    To process data pulled from NDBC website

    >>> import WDRT.NDBCdata as NDBCdata
    >>> swdList, freqList, dateVals = NDBCdata.fetchFromWeb(46089, savePath='data')
    >>> Hs, T, DateNum = NDBCdata.prepData(swdList, freqList, dateVals)
    '''
    n = len(swdList)
    Hs = []
    T = []
    DateNum = []
    for ii in range(n):
        tmp1, tmp2 = getStats(swdList[ii], freqList[ii])
        Hs.extend(tmp1)
        T.extend(tmp2)
        DateNum.extend(getDateNums(dateList[ii]))
    Hs = np.array(Hs, dtype=np.float)
    T = np.array(T, dtype=np.float)
    DateNum = np.array(DateNum, dtype=np.float)

    # Removing NaN data, assigning T label depending on input (Te or Tp)
    Nanrem = np.logical_not(np.isnan(T) | np.isnan(Hs))
    # Find NaN data in Hs or T
    DateNum = DateNum[Nanrem]  # Remove any NaN data from DateNum
    Hs = Hs[Nanrem]  # Remove any NaN data from Hs
    T = T[Nanrem]  # Remove any NaN data from T
    return Hs, T, DateNum

def getStats(swdArr, freqArr):
    '''Significant wave height and energy period

    Parameters
    ----------
        swdArr : np.array
            Numpy array of the spectral wave density data for a specific year
        freqArr: np.array
            Numpy array that contains the frequency values for a specific year

    Returns
    -------
        Hm0 : list
            Significant wave height.
        Te : list
            Energy period.
    '''
    Hm0 = []
    Te = []

    for line in swdArr:
        m_1 = np.trapz(line * freqArr ** (-1), freqArr)
        m0 = np.trapz(line, freqArr)
        Hm0.append(4.004 * m0 ** 0.5)
        np.seterr(all='ignore')
        Te.append(m_1 / m0)
    return Hm0, Te


def getDateNums(dateArr):
    '''datetime objects

    Parameters
    ----------
        dateArr : np.array
            Array of a specific years date vals from NDBC.fetchFromWeb

    Returns
    -------
        dateNum : np.array
            Array of datetime objects.
    '''
    dateNum = []
    for times in dateArr:
        if  times[0] < 1900:
            times[0] = 1900 + times[0]
        if times[0] < 2005:
            dateNum.append(date.toordinal(datetime(times[0], times[1],
                                                   times[2], times[3])))
        else:
            dateNum.append(date.toordinal(datetime(times[0], times[1],
                                                   times[2], times[3],
                                                   times[4])))
    return dateNum

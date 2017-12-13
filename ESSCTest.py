import numpy as np
import WDRT.ESSC as ESSC
import matplotlib.pyplot as plt
import collections

#Parameters for testing
n_size = 40. 
bin_1_limit = 1. 
bin_step = 0.25
Time_SS = 1.  # Sea state duration (hrs)
Time_R = 100  # Return periods (yrs) of interest

tol = 1e-5
def test():
	testFiles = loadTestData()
	testData = generateTestData()
	differences = getDiffereces(testFiles,testData)
	for contourType in differences:
		if differences[contourType] > tol:
			print "TEST FAILED - " + contourType
			print contourType + " - ", differences[contourType]
		else:
			print "TEST PASSED - " + contourType



def getDiffereces(testFiles, testData):
	differencePairs = [("PCA-HS" , np.sum(testFiles["PCA-Hs"] - testData["PCA-Hs"])),
				   ("PCA-T" , np.sum(testFiles["PCA-T"] - testData["PCA-T"])),
				   ("Gumbel-Hs" , np.sum(testFiles["Gumbel-Hs"] - testData["Gumbel-Hs"])),
				   ("Gumbel-T" , np.sum(testFiles["Gumbel-T"] - testData["Gumbel-T"])),
				   ("Gauss-Hs" , np.sum(testFiles["Gauss-Hs"] - testData["Gauss-Hs"])),
				   ("Gauss-T" , np.sum(testFiles["Gauss-T"] - testData["Gauss-T"])),
				   ("CC-Hs" , np.sum(testFiles["CC-Hs"] - testData["CC-Hs"])),
				   ("CC-T" , np.sum(testFiles["CC-T"] - testData["CC-T"])),
				   ("NanaParaGauss-Hs", np.sum(testFiles["NonParaGauss-Hs"] - testData["NonParaGauss-Hs"])),
				   ("NanaParaGauss-T", np.sum(testFiles["NonParaGauss-T"] - testData["NonParaGauss-T"])),
				   ("NanaParaGum-Hs", np.sum(testFiles["NonParaGum-Hs"] - testData["NonParaGum-Hs"])),
				   ("NanaParaGum-T", np.sum(testFiles["NonParaGum-T"] - testData["NonParaGum-T"])),
				   ("NanaParaCC-Hs", np.sum(testFiles["NonParaCC-Hs"] - testData["NonParaCC-Hs"])),				   
				   ("NanaParaCC-T", np.sum(testFiles["NonParaCC-T"] - testData["NonParaCC-T"])),
				   ("Rosen-Hs" , np.sum(testFiles["Rosen-Hs"] - testData["Rosen-Hs"])),
				   ("Rosen-T" , np.sum(testFiles["Rosen-T"] - testData["Rosen-T"]))]

	return collections.OrderedDict(differencePairs)





def loadTestData():
	testFiles = {}
	testFiles["PCA-Hs"] = np.load("TestFiles/pca46022-Hs.npy")
	testFiles["PCA-T"] = np.load("TestFiles/pca46022-T.npy")
	testFiles["Gauss-Hs"] = np.load("TestFiles/gauss46022-Hs.npy")
	testFiles["Gauss-T"] = np.load("TestFiles/gauss46022-T.npy")
	testFiles["Gumbel-Hs"] = np.load("TestFiles/gumbel46022-Hs.npy")
	testFiles["Gumbel-T"] = np.load("TestFiles/gumbel46022-T.npy")
	testFiles["CC-Hs"] = np.load("TestFiles/cc46022-Hs.npy")
	testFiles["CC-T"] = np.load("TestFiles/cc46022-T.npy")	
	testFiles["NonParaGauss-Hs"] = np.load("TestFiles/Gauss-NonPara46022-Hs.npy")
	testFiles["NonParaGauss-T"] = np.load("TestFiles/Gauss-NonPara46022-T.npy")
	testFiles["NonParaGum-Hs"] = np.load("TestFiles/Gum-NonPara46022-Hs.npy")
	testFiles["NonParaGum-T"] = np.load("TestFiles/Gum-NonPara46022-T.npy")
	testFiles["NonParaCC-Hs"] = np.load("TestFiles/CC-NonPara46022-Hs.npy")
	testFiles["NonParaCC-T"] = np.load("TestFiles/CC-NonPara46022-T.npy")
	testFiles["Rosen-Hs"] = np.load("TestFiles/rosen46022-Hs.npy")
	testFiles["Rosen-T"] = np.load("TestFiles/rosen46022-T.npy")
	print "--> Test Data loaded"
	return testFiles


def generateTestData():
	buoy46022 = ESSC.Buoy('46022','NDBC')

	buoy46022.fetchFromWeb()
	compareLoadMethods(buoy46022)

	pca46022 = ESSC.PCA(buoy46022)
	Gauss46022 = ESSC.GaussianCopula(buoy46022, n_size = n_size, bin_1_limit = bin_1_limit, bin_step = bin_step)
	Gumbel46022 = ESSC.GumbelCopula(buoy46022, n_size = n_size, bin_1_limit = bin_1_limit, bin_step = bin_step)
	Clayton46022 = ESSC.ClaytonCopula(buoy46022, n_size = n_size, bin_1_limit = bin_1_limit, bin_step = bin_step)
	rosen46022 = ESSC.Rosenblatt(buoy46022, n_size = n_size, bin_1_limit = bin_1_limit, bin_step = bin_step)
	NonParaGauss46022 = ESSC.NonParaGaussianCopula(buoy46022)
	NonParaClay46022 = ESSC.NonParaClaytonCopula(buoy46022)
	NonParaGum46022 = ESSC.NonParaGumbelCopula(buoy46022)


	pca_Hs_Return, pca_T_Return = pca46022.getContours(Time_SS, Time_R)
	Gauss_Hs_Return, Gauss_T_Return = Gauss46022.getContours(Time_SS, Time_R)
	Gumbel_Hs_Return, Gumbel_T_Return = Gumbel46022.getContours(Time_SS, Time_R)
	Clayton_Hs_Return, Clayton_T_Return = Clayton46022.getContours(Time_SS, Time_R)
	rosen_Hs_Return, rosen_T_Return = rosen46022.getContours(Time_SS, Time_R)
	NonParaGau_Hs_Return, NonParaGau_T_Return = NonParaGauss46022.getContours(Time_SS, Time_R)
	NonParaClay_Hs_Return, NonParaClay_T_Return = NonParaClay46022.getContours(Time_SS, Time_R)
	NonParaGum_Hs_Return, NonParaGum_T_Return = NonParaGum46022.getContours(Time_SS, Time_R)


	testData = {"PCA-Hs" : pca_Hs_Return,
				"PCA-T" : pca_T_Return,
				"Gauss-Hs" : Gauss_Hs_Return,
				"Gauss-T" : Gauss_T_Return,
				"Gumbel-Hs" : Gumbel_Hs_Return,
				"Gumbel-T" : Gumbel_T_Return,
				"CC-Hs" : Clayton_Hs_Return,
				"CC-T" : Clayton_T_Return,
				"NonParaGauss-Hs" : NonParaGau_Hs_Return,
				"NonParaGauss-T" : NonParaGau_T_Return,
				"NonParaGum-Hs" : NonParaGum_Hs_Return,
				"NonParaGum-T" : NonParaGum_T_Return,
				"NonParaCC-Hs" : NonParaClay_Hs_Return,
				"NonParaCC-T" : NonParaClay_T_Return,
				"Rosen-Hs" : rosen_Hs_Return,
				"Rosen-T" :rosen_T_Return}
	return testData


def compareLoadMethods(buoy):
	buoy.saveAsTxt('.\TestTxt')
	buoy.saveAsH5()
	txtBuoy = ESSC.Buoy('46022', 'NDBC')
	txtBuoy.loadFromText('.\TestTxt\NDBC46022\\')
	h5Buoy = ESSC.Buoy('46022', 'NDBC')
	h5Buoy.loadFromH5()

	error = False
	for i in range(len(buoy.Hs)):
		if buoy.Hs[i] - txtBuoy.Hs[i] > tol:
			error = True
			print "TEST FAILED = .txt Files"
			break
	if not error:
		print "TEST PASSED - .txt Files"
	error = False
	for i in range(len(buoy.Hs)):
		if buoy.Hs[i] - h5Buoy.Hs[i] > tol:
			error = True
			print "TEST FAILED = .h5 Files"
			break
	if not error:
		print "TEST PASSED - .h5 Files"

test()
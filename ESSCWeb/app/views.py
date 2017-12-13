from flask import render_template, flash, redirect, request
from app import app
from .forms import LoginForm, ESSCForm
import WDRT.ESSC as ESSC
import matplotlib.pyplot as plt

#def calculateContour(contourTypes):


@app.route('/')
@app.route('/ESSC_Form', methods=['GET', 'POST'])
def ESSC_Form():
	form = ESSCForm()
	if form.validate_on_submit():
		flash('NDBC Buoy Number =%s' % 
			(form.buoyNum.data))
		return render_template('ESSC_Results.html',
			                    title=('Contours for NDBC#' + form.buoyNum.data))
	return render_template('ESSC_Form.html',
		                    form=form)

@app.route('/ESSC_Results', methods=['GET','POST'])
def ESSC_Results():
	contours = {}
	buoyNum = request.form['buoyNum']
	buoyType = request.form.getlist('Buoy Type')[0]

	buoy = ESSC.Buoy(buoyNum, buoyType)
	buoy.fetchFromWeb()

	print request.form.getlist('Contour')

	if "PCA" in request.form.getlist('Contour'):
		PCA = ESSC.PCA(buoy)
		PCA.getContours(1.,100)
		contours["PCA"] = PCA
		print "pca"

	if "Gauss" in request.form.getlist('Contour'):
		gauss = ESSC.GaussianCopula(buoy)
		gauss.getContours(1.,100)
		contours["Gaussian Copula"] = gauss
		print "gauss"

	if "CC" in request.form.getlist('Contour'):
		CC = ESSC.ClaytonCopula(buoy)
		CC.getContours(1.,100)
		contours["Clayton Copula"] = CC
		print "cc"

	if "Gumbel" in request.form.getlist('Contour'):
		gumbel = ESSC.GumbelCopula(buoy)
		gumbel.getContours(1.,100)
		contours["Gumbel Copula"] = gumbel
		print "gumbel"

	if "Rosen" in request.form.getlist('Contour'):
		rosen = ESSC.GaussianCopula(buoy)
		rosen.getContours(1.,100)
		contours["Rosenblatt"] = rosen
		print "rosen"


	plt.figure()
	plt.plot(buoy.T, buoy.Hs, 'bo', alpha=0.1, label='NDBC data')
	for contour in contours:
		print contour
		plt.plot(contours[contour].T_ReturnContours, contours[contour].Hs_ReturnContours, '-', label=contour)
	#plt.plot(self.T_SampleFSS, self.Hs_SampleFSS, 'ro', label='full sea state samples')
	#plt.plot(self.T_SampleCA, self.Hs_SampleCA, 'y^', label='contour approach samples')
	plt.legend(loc='lower right', fontsize='small')
	plt.grid(True)
	plt.xlabel('Energy period, $T_e$ [s]')
	plt.ylabel('Sig. wave height, $H_s$ [m]')
	filePath =  buoyNum + '-'
	fileName = filePath
	for contour in request.form.getlist('Contour'):
		fileName += contour
	fileName += ".png"
	plt.savefig("./app/static/images/" + fileName)
	return render_template('ESSC_Results.html',
							title = "ESSC Results",
		                    imgsrc= "../static/images/" + fileName)


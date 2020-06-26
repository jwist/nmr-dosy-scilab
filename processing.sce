
// Procesamiento (Doneshot45) (20180327)
clear

// to close all the open window
//mclose('all')

// =============================================================================
// automatic extraction of parameters and spectra
// =============================================================================
// this section of the script is case specific and Bruker specific
// it allows to grab the 1D traces for each gradients and to extract
// the important parameters
//

// =============================================================================
// modify these parameters according to your needs
// path to experiment (expname)
//path = "/home/jul/git/jwist/nmr-dosy-scilab/20180912-coaxial-Ciclohexane_convection/";
//path = "/home/jul/git/jwist/nmr-dosy-scilab/20181217-coaxial-ciclohexane_D20_H2O_G/";
path = "/home/jul/git/jwist/nmr-dosy-scilab/gradient_calibration/"; // not well shimmed
//path='/run/media/jul/44B6-1E16/data/summer-school-test/'

// experiment number of the peudo 2D
expno = "3";

// spectral window where signal is to be observed (point index, not ppm)
integrationLowerLimit = 9500;
integrationUpperLimit = 10500;
// =============================================================================
// extraction of the information

// open and reads acqus bruker file
idr=mopen(path + expno +"/acqus", "r");
acqus = mgetl(idr,-1); mclose(idr);

// extracts delays
delayIndex = grep(acqus, "##$D=");
delayText  = tokens(acqus(delayIndex + 1) + acqus(delayIndex + 2));
D          = msscanf(-1, delayText, "%g")';

// extracts pulse lengths
pulseLengthIndex = grep(acqus,"##$P=");
pulseLengthText  = tokens(acqus(pulseLengthIndex + 1) + acqus(pulseLengthIndex + 2) + acqus(pulseLengthIndex + 3));
P                = msscanf(-1, pulseLengthText,"%g")';

// extracts TD
TDIndex           = grep(acqus,"##$TD=");
TDText            = tokens(acqus(TDIndex));
TD                = strtod(TDText(2));

// extracts constants
constantsIndex = grep(acqus,"##$CNST=");
constantsText  = tokens(acqus(constantsIndex + 1) + acqus(constantsIndex + 2));
CNST           = msscanf(-1, constantsText,"%g")';

// load gradient calibration file
trueGradientList = fscanfMat(path + expno + "/difflist");
// load gradient calibration file
//trueGradientList = fscanfMat("/home/jul/git/jwist/nmr-dosy-scilab/gm1.txt");

// number of different gradients acquired
numberOfGradients = length(trueGradientList);

// create the experiment number of the 1D trace after "splitser" is performed in topspin
// the 1D traces are created consecutively.
expnos1D = linspace(strtod(expno) + 1, strtod(expno) + numberOfGradients, numberOfGradients)

for i = 1:length(expnos1D)
  // opens the processed binary file 1r and store it in a buffer
  spect = mopen(path + string(expnos1D(1,i)) + '/pdata/1/1r', 'rb')
  Spectrum(:, i) = mtlb_fread(spect, TD,'int32')
  mclose(spect)

  // opens and read the procs file to... 
  idr = mopen(path + string(expnos1D(1, i)) + '/pdata/1/procs', 'r')
  procs = mgetl(idr, -1)
  mclose(idr)

  // ...obtain necessary processing parameters
  [w, r] = grep(procs,'##$NC_proc')
  [a, b, c, d] = regexp(procs(w), '/(?P<name>\w+)=\s+(?P<digit>\D*\d+)/')
  fa = strtod(d(2))

  Spectrum(:, i) = Spectrum(:, i) .* (2.^fa)
end

// =============================================================================
// definition of the region of interest (the signal you wish to analyze) and
// extraction of the signal heigth/intensity
// =============================================================================

// if not using Bruker equipment you must provide the following two information.

// 1. a vector of true(calibrated) magnitudes for gradients in increading order
// trueGradientList = [1, 5, 10, 15, 20, 30, 50, 80];

// 2. a matrix of spectra, each spectrum is a column
// Spectrum = [[spect1];[spect2]];

// =============================================================================
// here we display the spectra to identify the region of interest
// stack plot of the extracted spectra
spectrumLength = size(Spectrum, 1);
// control the position of the trace in the stack plot
xOffset = 1000;
yOffset = 2.5e6;
xAxis = linspace(1, spectrumLength, spectrumLength);

scf(0)
for k = 1:length(expnos1D)
  subplot(1,3,1);
  plot2d(xOffset * k + xAxis, k * yOffset + Spectrum(:, k), style = 2);
  e = gce();
  e.children(1).thickness = 2;
end

// define the region of interest
// when run for the first time use the following loop to better
// determine the range of point that contains the signal to be fitted
scf(1)
  
for k = 1:8
  subplot(1,2,1);
  plot2d(Spectrum(:,k));
  title("whole window", "fontsize", 3, "color", "blue");
  xlabel("ppm", "fontsize", 3, "color", "blue");
  ylabel("relative intensity", "fontsize", 3, "color", "blue");

  a = gca();
  a.font_size = 2;
  e = gce();
  e.children(1).thickness = 2;
end

// this allows to verify that the region is correctly selected
for k = 1:8
  subplot(1,2,2);
  plot2d(Spectrum(integrationLowerLimit:integrationUpperLimit,k));
  title("selected region", "fontsize", 3, "color", "blue");
  xlabel("ppm", "fontsize", 3, "color", "blue");
  ylabel("relative intensity", "fontsize", 3, "color", "blue");
  
  a = gca();
  a.font_size = 2;
  e = gce();
  e.children(1).thickness = 2;
end

// =============================================================================
// extraction of peak height/intensity from the region of interest

// plain integration by summing the discrete point of the spectra. More complexe
// procedure are to be tried before concluding that such a simplistic approach
// can indeed be used.
// finding max of the peak. Beware this is a very simplistic approach doesn't 
// corresponds to "true peak" maximum position that shoulittleDelta be obtained 
// by fitting the signal by an ad-hoc function.
for i = 1:size(expnos1D, 2)
  regionOfInterest = real(Spectrum(integrationLowerLimit:integrationUpperLimit, i))
  observedIntensities(1, i) = sum(regionOfInterest);
  observedPositions(1, i) = max(regionOfInterest);
end

// =============================================================================
// extraction of the diffusion coefficients
// =============================================================================

// =============================================================================
// if not using Bruker instruments, please provide the following information
// manually:

// big delta [seg]
bigDelta = D(21); 

// little delta [seg] 
littleDelta = 2 * P(30) * 1e-6;

// gyromagnetic constant [rad/s]
r = 2 * %pi * 4.258e+03

// alpha (specific to oneshot sequence, to balance gradients)
alpha = CNST(15);

// tau [???]
tau = CNST(18);

// =============================================================================
// definition of function for common experiments

// function oneshot-45
// Assuming half-sine gradient pulses [most common on Bruker systems]
function y = yth(Diff, trueGradientList, r, littleDelta, bigDelta, alpha, tau)
  y = Diff(1) * exp(-Diff(2) * 1e+4 * trueGradientList.^2 * (r * littleDelta)^2 * (bigDelta - littleDelta * ((5 - 3 * alpha^2) / 16) - tau * ((1 - alpha^2) / 2)));
endfunction

// function Bipolar pulse pair (BPP) LED 
function e = errorFunction(Diff, trueGradientList, r, littleDelta, bigDelta, alpha, tau, ym, wm)
  e = wm .* (yth(Diff, trueGradientList, r, littleDelta, bigDelta, alpha, tau) - ym)
endfunction


// =============================================================================
// fitting of the observed data

// normalization to the biggest intensity in the decay (not necessary)
normalizedObservedIntensities = (observedIntensities ./ max(observedIntensities))' .* 100;
normalizedObservedIntensities = (observedIntensities ./ max(observedIntensities))' .* 100;

// weights to ponderate the experimental points if necessary
weights = ones(size(trueGradientList, 1), 1);

// initial values for fitting (should be close enough... whatever this means)
Diff0=[0.3 9.919e-10]

options = list(errorFunction, ...
trueGradientList, ...
r, ...
littleDelta, ...
bigDelta, ...
alpha, ...
tau, ...
normalizedObservedIntensities, ...
weights);

// minimization of the error function
[f,xopt, gopt] = leastsq(options, Diff0);

// computation of the fitted values
fittedIntensities = yth(xopt, trueGradientList, r, littleDelta, bigDelta, alpha, tau);

// =============================================================================
// graph the results

scf(0);
subplot(1, 3, 2);
plot2d(trueGradientList, [normalizedObservedIntensities, fittedIntensities]);

xlabel("(gauss/cm)", "fontsize", 3,"color", "blue");
ylabel("S/S(0)", "fontsize", 3, "color", "blue");

e = gce();
e.children(1).thickness = 3;
e.children(1).foreground = 5;

e.children(2).foreground = 1;
e.children(2).mark_style = 9;
e.children(2).line_mode="off"
hl=captions(e.children,["fitted intensities"; "observed intensities"]);

logScale = trueGradientList.^2 .* ((r*littleDelta)^2 * ...
((bigDelta - littleDelta * ((5 - 3 * alpha^2) / 16) - ...
tau * ((1 - alpha^2) / 2))) * 1e-4);

subplot(1, 3, 3);
plot2d(logScale, [log(fittedIntensities), log(normalizedObservedIntensities)]);

//plot2d(logScale, log(normalizedObservedIntensities),style=-5);
xlabel("q2((bigDelta-littleDelta*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))", "fontsize", 3,"color", "blue");
ylabel("Ln(S/S(0))", "fontsize", 3, "color", "blue");
e = gce();
e.children(2).thickness = 3;
e.children(2).foreground = 5;

e.children(1).foreground = 1;
e.children(1).mark_style = 9;
e.children(1).line_mode="off"
hl=captions(e.children,["fitted intensities"; "observed intensities"]);

// =============================================================================
// display the diffusion coefficient

result = "intersection: " + string(xopt(1)) + " | diffusion coeff: " + string(xopt(2));
disp(result);

q=sqrt(log(normalizedObservedIntensities)*(-1/(2.299D-09))*(1/(r*littleDelta))*(1/(r*littleDelta))*(1/(bigDelta-littleDelta*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2))))*0.01;
gcal=sqrt(log(fittedIntensities./xopt(1,1))*(1/(2.29e-9*(r)^2*(littleDelta)^2*1e+4*(bigDelta-(1*littleDelta/3)-(tau/2)))))



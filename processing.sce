// Procesamiento (Doneshot45) (20180327)
clear

// to close all the open window
//mclose('all')

// path to experiment (expname)
path = '/home/jul/git/selfdiffusion_tutorial/20180912-coaxial-Ciclohexane_convection/' 

// experiment number of the peudo 2D
expno = '39'

// load gradient calibration file
gm = fscanfMat('/home/jul/git/selfdiffusion_tutorial/G_Cal_5%_80%_ref_100%max_H2O_end.txt')

// number of different gradients acquired
numberOfGradients = 8 

// create the experiment number of the 1D trace after "splitser" is performed in topspin
// the 1D traces are created consecutively.
expnos1D = linspace(strtod(expno) + 1, strtod(expno) + 8, numberOfGradients)

for i = 1:length(expnos1D)
    // opens the processed binary file 1r and store it in a buffer
    spect = mopen(path + string(expnos1D(1,i)) + '/pdata/1/1r', 'rb')
    Spectrum(:, i) = mtlb_fread(spect, 2*16384,'int32')
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

// stack plot of the extracted spectra
spectrumLength = size(Spectrum, 1)
xOffset = 1000
yOffset = 2.5e6
xAxis = linspace(1, spectrumLength, spectrumLength)

for k = 1:length(expnos1D)
    subplot(1,3,1)
    plot2d(xOffset * k + xAxis, k * yOffset + Spectrum(:, k), style = 2) 
end

// calculate max
[n, m] = max(Spectrum(:, 1))

// define the region of interest
// when run for the first time use the following loop to better
// determine the range of point that contains the signal to be fitted
if 1 == 2 then
    for k=1:8
        plot2d(Spectrum(:,k))
    end
end

// when finished set the range here
signal_1_int_reg_low = 14900
signal_1_int_reg_upp = 15300

for i = 1:size(expnos1D, 2)
I1(1,i)=sum(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i)));
end

for i=1:size(expnos1D,2)
I2(1,i)=max(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i)));
end


//path='/Users/christianpantoja/Documents/opt/20180320_dioxane_splitser/7/'
//path='/Users/christianpantoja/Documents/opt/20180322_dioxane_splitser_9/9/'

//idr=mopen('/Users/christianpantoja/Documents/opt/20180320_dioxane_splitser/7/acqus','r');
idr=mopen(path + expno +'/acqus','r');
//idr=mopen('/Users/christianpantoja/Documents/opt/20180322_dioxane_splitser_9/9/acqus','r');

acqus = mgetl(idr,-1); mclose(idr);

[w,r]=grep(acqus,'##$D=');

mputl(acqus(w+1),path + expno +'/D.txt')
mputl(acqus(w+2),path + expno +'/D2.txt')
//mputl(acqus(w+3),path + expno +'/D3.txt')
d1=fscanfMat(path + expno +'/D.txt')
d2=fscanfMat(path + expno +'/D2.txt')
//d3=fscanfMat(path + expno +'/D3.txt')

D=cat(2,d1,d2)

[w,r]=grep(acqus,'##$P=');

mputl(acqus(w+1),path + expno +'/P.txt')
mputl(acqus(w+2),path + expno +'/P2.txt')
mputl(acqus(w+3),path + expno +'/P3.txt')
P1=fscanfMat(path + expno +'/P.txt')
P2=fscanfMat(path + expno +'/P2.txt')
P3=fscanfMat(path + expno +'/P3.txt')

P=cat(2,P1,P2,P3)

[w,r]=grep(acqus,'##$CNST=');
mputl(acqus(w+1),path + expno +'/CNST.txt')
mputl(acqus(w+2),path + expno +'/CNST2.txt')
CNST=fscanfMat(path + expno +'/CNST.txt')
CNST2=fscanfMat(path + expno +'/CNST2.txt')

CNST=cat(2,CNST,CNST2)

// Parametros 

BD=D(1,21) //Big delta (seg)
LD=2*P(1,31)*1e-6 // little delta (seg) 
r=2*%pi*4.258e+03
alpha=CNST(1,15)// 
tau=CNST(1,18)//
//tau=0.001223

// function Doneshot45 // Assuming half-sine gradient pulses [most common on Bruker systems]

function y=yth(Diff,g,r,LD,BD,alpha,tau)
   y = Diff(1)*exp(-Diff(2)*1e+4*g.^2*(r*LD)*(r*LD)*(BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))
endfunction


//function y=yth(Diff,g,r,LD,BD,alpha,tau)
   //y = 1*exp(-Diff(1)*1e+4*g^2*(r*LD)*(r*LD)*(BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))
//endfunction
//ym=(I1./I1(1,1))'
ym=(I1./1.4e+8)'
//Señal muy baja////ym=(I1./1.4e+7)'
//ym=(I1./0.4e+8)'
//ym=(I2./10e+6)'
//ym=fscanfMat('/Users/christianpantoja/Desktop/I1.txt')'

//ym(4,1)=1.4
wm = ones(size(gm,1),1);

Diff0=[0.3 9.919e-10]

//Diff0=[0.01 1.5e-9]
//Diff0=[9.919e-10]
//error function 

function e=myfun(Diff,gm,r,LD,BD,alpha,tau,ym,wm)
   e = wm.*(yth(Diff,gm,r,LD,BD,alpha,tau) - ym)
endfunction

// function Bipolar pulse pair (BPP) LED 
[f,xopt, gopt] = leastsq(list(myfun,gm,r,LD,BD,alpha,tau,ym,wm),Diff0)


yy = yth(xopt,gm,r,LD,BD,alpha,tau);

subplot(1,3,2)
plot2d(gm,ym,style=-9)
xlabel("(chemical shift (ppm))", "fontsize", 3,"color", "blue");
ylabel("Relative intensity", "fontsize", 3, "color", "blue");

subplot(1,3,2)
plot2d(gm,yy,style=5)
xlabel("(gauss/cm)", "fontsize", 3,"color", "blue");
ylabel("S/S(0)", "fontsize", 3, "color", "blue");

subplot(1,3,3)
plot2d(gm.^2*(r*LD)*(r*LD)*((BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))*1e-4,log(yy),style=5)
plot2d(gm.^2*(r*LD)*(r*LD)*((BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))*1e-4,log(ym),style=-5)
xlabel("q2((BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))", "fontsize", 3,"color", "blue");
ylabel("Ln(S/S(0))", "fontsize", 3, "color", "blue");
a = gca()
//a.data_bounds = [0.1,-2.8;15,0.0]


disp(xopt)

//yy = yth(Diff0,gm,r,LD,BD,alpha,tau);

//subplot(1,3,2)
//plot2d(gm,yy,style=3)``


//q2=sqrt(log(yy/xopt(1,1))*(1/(BD-(LD/3)-(tau/2)))*(-1/(1.902D-09)))*(1/r)*(1/LD)*0.01

// water

//q1=sqrt(log(ym/xopt(1,1))*(-1/(2.299D-09))*(1/(r*LD))*(1/(r*LD))*(1/(BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2))))*0.01

q=sqrt(log(ym)*(-1/(2.299D-09))*(1/(r*LD))*(1/(r*LD))*(1/(BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2))))*0.01
// cyclohexane

//q=sqrt(log(ym/xopt(1,1))*(-1/(1.4021D-09))*(1/(r*LD))*(1/(r*LD))*(1/(BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2))))*0.01

// dodecane

//q=sqrt(log(ym/xopt(1,1))*(-1/(0.80329D-09))*(1/(r*LD))*(1/(r*LD))*(1/(BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2))))*0.01

//dioxane

//q=sqrt(log(ym/xopt(1,1))*(-1/(1.07446D-09))*(1/(r*LD))*(1/(r*LD))*(1/(BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2))))*0.01


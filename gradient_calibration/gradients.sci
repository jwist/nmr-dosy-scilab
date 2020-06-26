// Procesamiento (Doneshot45) (20180327)
clear
//mclose('all')


//path='/home/jul/git/jwist/nmr-dosy-scilab/gradient_calibration/'
path='/run/media/jul/44B6-1E16/data/summer-school-test/'
path2=path + '60/'

n1=60

// Gradients values 80% max, ramp (5%-80%) 8-pts
// Calibration H2O_doneshot_2d_laura  (20180611)

//Gradientes arrojados en DIFFLIST
gm = fscanfMat(path2+ 'difflist');

////"Al corregir cambiar gm por valor arrojado de la ultima ecuacion, q1"
////............Lista corregidos..............
//gm=[1.0487698  
//    5.8394707  
//    9.3746098  
//    13.390563  
//    17.062744  
//    20.758365  
//    24.734899  
//    28.772867]:
//gm=fscanfMat('C:\Users\Lorraine\Desktop\NEW_CODE\G_Cal_5%_80%_ref_100%max_H2O_end.txt')
//


N=linspace(n1+1,n1+11,11)

//   for i=1:size(N,2)
//        idr=mopen(path+string(N(1,i))+'/pdata/1/1r','rb');
//        Spectrum(i,:) = mtlb_fread(idr, 32768/2,'int32'); mclose(idr);
//        idr = mopen(path+string(N(1,i))+'/pdata/1/procs','r');
//        procs = mgetl(idr,-1); mclose(idr);
//        [w,r]=grep(procs,'##$NC_proc');
//        [a,b,c,d]=regexp(procs(w),'/(?P<name>\w+)=\s+(?P<digit>\D*\d+)/');
//        fa = strtod(d(2));
//        Spectrum(:,i)= Spectrum(:,i).*(2.^fa);
//   end
      for i=1:size(N,2)
        idr=mopen(path+string(N(1,i))+'/pdata/1/1r','rb');
        Spectrum(:,i) = mtlb_fread(idr, 16384,'int32'); mclose(idr);
        idr = mopen(path+string(N(1,i))+'/pdata/1/procs','r');
        procs = mgetl(idr,-1); mclose(idr);
        [w,r]=grep(procs,'##$NC_proc');
        [a,b,c,d]=regexp(procs(w),'/(?P<name>\w+)=\s+(?P<digit>\D*\d+)/');
        fa = strtod(d(2));
        Spectrum(:,i)= Spectrum(:,i).*(2^fa);
   end
   
   
   for k=1:size(N,2)
       subplot(1,3,1)
       plot2d((1000*k)+linspace(1,size(Spectrum,1),size(Spectrum,1)),k*(1.5e+6)+Spectrum(:,k),style=2) 
   end

// calculate max
[n,m]=max(Spectrum(:,1))

//signal_1_int_reg_low=m-250
//signal_1_int_reg_upp=m+310

//H2O 298_K LB= 0.3 
signal_1_int_reg_low=10000
signal_1_int_reg_upp=10600
//

for i=1:size(N,2)
I1(1,i)=sum(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i)));
end
//subplot(1,3,2)
//plot2d(gm,I1./0.7e+8,style=-9)

for i=1:size(N,2)
I2(1,i)=max(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i)));
end

idr=mopen(path2+'/acqus','r');

acqus = mgetl(idr,-1); mclose(idr);

[w,r]=grep(acqus,'##$D=');

mputl(acqus(w+1),path2+'/D.txt')
mputl(acqus(w+2),path2+'/D2.txt')
//mputl(acqus(w+3),path2+'/D3.txt')
d1=fscanfMat(path2+'/D.txt')
d2=fscanfMat(path2+'/D2.txt')
//d3=fscanfMat(path2+'/D3.txt')

D=cat(2,d1,d2)

[w,r]=grep(acqus,'##$P=');

mputl(acqus(w+1),path2+'/P.txt')
mputl(acqus(w+2),path2+'/P2.txt')
mputl(acqus(w+3),path2+'/P3.txt')
P1=fscanfMat(path2+'/P.txt')
P2=fscanfMat(path2+'/P2.txt')
P3=fscanfMat(path2+'/P3.txt')

P=cat(2,P1,P2,P3)

[w,r]=grep(acqus,'##$CNST=');
mputl(acqus(w+1),path2+'/CNST.txt')
mputl(acqus(w+2),path2+'/CNST2.txt')
CNST=fscanfMat(path2+'/CNST.txt')
CNST2=fscanfMat(path2+'/CNST2.txt')

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

//ym=(I1./I1(1,1))'
ym=(I1./1.4e+8)'
//Señal muy baja////ym=(I1./1.4e+7)'
//ym=(I1./0.4e+8)'
//ym=(I2./10e+6)'

//ym(4,1)=1.4
wm = ones(size(gm,1),1);

Diff0=[0.3 9.919e-10]

//error function 

function e=myfun(Diff,gm,r,LD,BD,alpha,tau,ym,wm)
   e = wm.*(yth(Diff,gm,r,LD,BD,alpha,tau) - ym)
endfunction

// function Bipolar pulse pair (BPP) LED 
[f,xopt, gopt] = leastsq(list(myfun,gm,r,LD,BD,alpha,tau,ym,wm),Diff0)

yy = yth(xopt,gm,r,LD,BD,alpha,tau);

subplot(1,2,1)
plot2d(gm,ym,style=-9)
xlabel("(chemical shift (ppm))", "fontsize", 3,"color", "blue");
ylabel("Relative intensity", "fontsize", 3, "color", "blue");

subplot(1,2,1)
plot2d(gm,yy,style=5)
xlabel("(gauss/cm)", "fontsize", 3,"color", "blue");
ylabel("S/S(0)", "fontsize", 3, "color", "blue");

subplot(1,2,2)
plot2d(gm.^2*(r*LD)*(r*LD)*((BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))*1e-4,log(yy),style=5)
plot2d(gm.^2*(r*LD)*(r*LD)*((BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))*1e-4,log(ym),style=-5)
xlabel("q2((BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2)))", "fontsize", 3,"color", "blue");
ylabel("Ln(S/S(0))", "fontsize", 3, "color", "blue");
a = gca()
//a.data_bounds = [0.1,-2.8;15,0.0]


disp(xopt)


//q2=sqrt(log(yy/xopt(1,1))*(1/(BD-(LD/3)-(tau/2)))*(-1/(1.902D-09)))*(1/r)*(1/LD)*0.01

//correccion de gradientes  HOLZ 2000 (Diff=2.299D-09)
//El valor q1 arroja los gradientes corregidos
///
// water
q1=sqrt(log(ym/xopt(1,1))*(-1/(2.299D-09))*(1/(r*LD))*(1/(r*LD))*(1/(BD-LD*((5-3*alpha*alpha)/16)-tau*((1-alpha*alpha)/2))))*0.01




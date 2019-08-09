clear
mclose('all') 

// BPP-LED

path='/run/media/jul/44B6-1E16/data/summer-school-test/'
path2='/run/media/jul/44B6-1E16/data/summer-school-test/200/'


//path='/opt/topspin3.5pl7/data/600cryo/20190318_ref_dataset/' // 600 cryo
//path2='/opt/topspin3.5pl7/data/600cryo/20190318_ref_dataset/13' // 600 cryo


// N: Vector con experimentos 

// n1: Folder inicial 
n1=201
// nt: Numero total de folders 
nt=11

// N 
N=linspace(n1,n1+(nt-1),nt)

// Grandient list (cm2/s)

//gm=fscanfMat('/Users/cpantoj/Documents/Calibration_z_Gradient/difflist_G_800cryo.txt')
//gm=fscanfMat('/Users/cpantoj/Documents/Calibration_z_Gradient/difflist_G_800cryo_c.txt')
//gm=fscanfMat('/Users/cpantoj/Documents/Calibration_z_Gradient/difflist_G_800cryo_zprofile.txt')

//gm=fscanfMat('/Users/cpantoj/Documents/Calibration_z_Gradient/difflist_G_800cryo_bbp_led.txt')
//gm=fscanfMat('/Users/cpantoj/Documents/Calibration_z_Gradient/difflist_G_600cryo_bbp_led_cal.txt')

//gm=fscanfMat('/Users/cpantoj/Documents/Calibration_z_Gradient/difflist_G_800cryo_bbp_led_c.txt')
//gm=fscanfMat('/Users/cpantoj/Documents/Calibration_z_Gradient/difflist_G_800cryo_bbp_led_c_7_5mm.txt')

//gm=fscanfMat('/Users/cpantoj/Documents/Calibration_z_Gradient/difflist_G_800cryo_zprofile_bbpled.txt')// z_profile_calibration 
gm = fscanfMat(path2+ 'difflist');




//gm=gm(1:7,1);
// ****** Matrix espectral *****

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
       
              //plot2d((0*k)+linspace(1,size(Spectrum,1),size(Spectrum,1)),k*(10e+1)+Spectrum(:,k),style=5) 
       plot2d(linspace(9300,9500,401)'+k*10, k*(10e+4)+Spectrum(9200:9600,k),style=18) // 800cryo
              //plot2d(k*(10e+1)+Spectrum(8174:8208,k),style=2) // 600cryo
              //plot2d(k*(10e+1)+Spectrum(20840:20920,k),style=14)
//              //plot2d((1*k)+linspace(1,size(Spectrum(5000:25000,1)),size(Spectrum(5000:25000,1)),k*(10e+1)+Spectrum(5000:25000,k),style=5) 
   end


   
   
   //         // *** Integration area 
   signal_1_int_reg_low=9350
   signal_1_int_reg_upp=9450
   


//   //         // *** Integration area  //600cryo
//   signal_1_int_reg_low=8174
//   signal_1_int_reg_upp=8208
   
   
   // Loop integration ***
  for i=1:size(N,2)
  I1(1,i)=sum(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i))); 
    I1(1,i)=max(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i))); 
  end


// Extracción de parametros 

idr=mopen(path2+'/acqus','r');
acqus = mgetl(idr,-1); mclose(idr);

[w,r]=grep(acqus,'##$D=');

mputl(acqus(w+1),path2+'/D.txt')
mputl(acqus(w+2),path2+'/D2.txt')


d1=fscanfMat(path2+'/D.txt')
d2=fscanfMat(path2+'/D2.txt')

// **** // ****

D=cat(2,d1,d2)

[w,r]=grep(acqus,'##$P=');
// **** // ****
mputl(acqus(w+1),path2+'/P.txt')
mputl(acqus(w+2),path2+'/P2.txt')
//mputl(acqus(w+3),path2+'/P3.txt')
P1=fscanfMat(path2+'/P.txt')
P2=fscanfMat(path2+'/P2.txt')
//P3=fscanfMat(path2+'/P3.txt')

// **** // ****
//P=cat(2,P1,P2,P3)
P=cat(2,P1,P2)
//
//// Parametros 

BD=D(1,21) //**** Diffusion_delta (seg)
tau=D(1,17) // *****  tau (seg)
LD=2*P(1,31)*1e-6 //**** Little_delta_1 (seg) 
//p28=1*P(1,29)*1e-6 // 
//r=42.58e+6
//r=267.53803e+6
//r=26752.219
r=26753.803// rad s-1 G-1
//LD2=2*50e-6+2*d19*9+2*p28*0.087+2*p28*0.206+2*p28*0.413+2*p28*0.778+2*p28*1.491 //**** Little_delta_2 (seg)
//GPZ2=9.6299
//GPZ2=12.96
//

function y=yth(Diff,g,r,LD,BD,tau)
   //y = Diff(1)*exp(-Diff(2)*(r*LD)*(r*LD)*(BD-(4*LD/3)-2*LD2)*((1e+4*GPZ2)-(1e+4*g))*((1e+4*GPZ2)-(1e+4*g)))
   // y = Diff(1)*exp(-Diff(2)*(r)^2*1e+4*(LD)^2*(g)^2 * (BD-(1*LD/3)-(tau/2)))
   y = Diff(1)*exp(-Diff(2)*((r)*(LD)*(g)*100)^2 * (BD-(1*LD/3)-(tau/2)))
   
      //y = Diff(1)*exp(-Diff(2)*(r*LD)*(r*LD)*(BD-(4*LD/3)-(LD2))*1e+4*(GPZ2-(g))^2)
   
endfunction

//****

//ym=(I1./1.4e+5)'
ym=(I1./max(I1))'
//
wm = ones(size(gm,1),1);
//
Diff0=[0.3 9.919e-9]
//
function e=myfun(Diff,g,r,LD,BD,tau,ym,wm)
   e = wm.*(yth(Diff,gm,r,LD,BD,tau) - ym)
endfunction
//
//// function Bipolar pulse pair (BPP) LED 
[f,xopt, gopt] = leastsq(list(myfun,gm,r,LD,BD,tau,ym,wm),Diff0)
//
yy = yth(xopt,gm,r,LD,BD,tau);
//
//
scf(1)
subplot(1,2,1)
plot2d(gm,ym,style=-9)
plot2d(gm,yy,style=5)
xlabel("(gauss/cm)", "fontsize", 3,"color", "blue");
ylabel("S/S(0)", "fontsize", 3, "color", "blue");
//
subplot(1,2,2)
plot2d((r)^2*(LD)^2*1e+4*(gm)^2*(BD-(1*LD/3)-(tau/2)),log(yy),style=5)
plot2d((r)^2*(LD)^2*1e+4*(gm)^2*(BD-(1*LD/3)-(tau/2)),log(ym),style=-5)
xlabel("q2((r)^2*(LD)^2*1e+4*(g)^2*(BD-(1*LD/3)-(tau/2)))", "fontsize", 3,"color", "blue");
ylabel("Ln(S/S(0))", "fontsize", 3, "color", "blue");
//a = gca()
//
disp(xopt)
//
//


gcal=sqrt(log(yy/xopt(1,1))*(1/(-2.29e-9*((r)*(LD)*100)^2*(BD-(1*LD/3)-(tau/2)))))







// setdiffparm
// xf2
// phase corr
// dosy2d

// gradpar

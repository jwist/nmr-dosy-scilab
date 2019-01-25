clear
mclose('all')

path='/opt/topspin3.5pl7/data/chpa800cryo/20181228_dropplet_milk/'
path2='/opt/topspin3.5pl7/data/chpa800cryo/20181228_dropplet_milk/53'

// N: Vector con experimentos 

// n1: Folder inicial 
n1=54

// nt: Numero total de folders 
nt=12

N=linspace(n1,n1+(nt-1),nt)

N=N(1,1:12);

// Grandient list (cm2/s)
//gm=fscanfMat('/Users/christianpantoja/Desktop/gm800cryo.txt')
//gm=fscanfMat('/Users/christianpantoja/Desktop/gm_5%_80%_8pts_800cryo.txt')
gm=fscanfMat('/Users/christianpantoja/Desktop/G_13C_gm.txt')
//gm=gm(1:7,1);
// ****** Matrix espectral *****

   for i=1:size(N,2)
        idr=mopen(path+string(N(1,i))+'/pdata/1/1r','rb');
        Spectrum(:,i) = mtlb_fread(idr, 1*16384,'int32'); mclose(idr);
        idr = mopen(path+string(N(1,i))+'/pdata/1/procs','r');
        procs = mgetl(idr,-1); mclose(idr);
        [w,r]=grep(procs,'##$NC_proc');
        [a,b,c,d]=regexp(procs(w),'/(?P<name>\w+)=\s+(?P<digit>\D*\d+)/');
        fa = strtod(d(2));
        Spectrum(:,i)= Spectrum(:,i).*(2^fa);
   end
   
      for k=1:size(N,2)
       subplot(1,3,1)
       plot2d((100*k)+linspace(1,size(Spectrum,1),size(Spectrum,1)),k*(1e+9)+Spectrum(:,k),style=5) 
   end

   
   
// *** Integration area 
   signal_1_int_reg_low=7050
   signal_1_int_reg_upp=7080
   
   // Loop integration ***
  for i=1:size(N,2)
  I1(1,i)=sum(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i))); 
  end


// Extracción de parametros 

idr=mopen(path2+'/acqus','r');
acqus = mgetl(idr,-1); mclose(idr);

[w,r]=grep(acqus,'##$D=');

mputl(acqus(w+1),path2+'/D.txt')
mputl(acqus(w+2),path2+'/D2.txt')

d1=fscanfMat(path2+'/D.txt')
d2=fscanfMat(path2+'/D2.txt')

D=cat(2,d1,d2)

[w,r]=grep(acqus,'##$P=');

mputl(acqus(w+1),path2+'/P.txt')
mputl(acqus(w+2),path2+'/P2.txt')
mputl(acqus(w+3),path2+'/P3.txt')
P1=fscanfMat(path2+'/P.txt')
P2=fscanfMat(path2+'/P2.txt')
P3=fscanfMat(path2+'/P3.txt')

P=cat(2,P1,P2,P3)

// Parametros 
BD=D(1,21) //Big delta (seg)
LD=1*P(1,31)*1e-6 // little delta (seg) 
r=2*%pi*1.07084e+03


function y=yth(Diff,g,r,LD,BD)
   y = Diff(1)*exp(-Diff(2)*1e+4*g^2*(r*LD)*(r*LD)*BD)
endfunction


ym=(I1./1.4e+6)'

wm = ones(size(gm,1),1);

Diff0=[0.3 9.919e-10]

function e=myfun(Diff,gm,r,LD,BD,ym,wm)
   e = wm.*(yth(Diff,gm,r,LD,BD) - ym)
endfunction

// function Bipolar pulse pair (BPP) LED 
[f,xopt, gopt] = leastsq(list(myfun,gm,r,LD,BD,ym,wm),Diff0)

yy = yth(xopt,gm,r,LD,BD);


subplot(1,3,2)
plot2d(gm,ym,style=-9)
xlabel("(chemical shift (ppm))", "fontsize", 3,"color", "blue");
ylabel("Relative intensity", "fontsize", 3, "color", "blue");

subplot(1,3,2)
plot2d(gm,yy,style=5)
xlabel("(gauss/cm)", "fontsize", 3,"color", "blue");
ylabel("S/S(0)", "fontsize", 3, "color", "blue");

subplot(1,3,3)
plot2d(gm^2*(r*LD)*(r*LD)*BD*1e-4,log(yy),style=5)
plot2d(gm^2*(r*LD)*(r*LD)*BD*1e-4,log(ym),style=-5)
xlabel("q2(gm^2*(r*LD)*(r*LD)*BD)", "fontsize", 3,"color", "blue");
ylabel("Ln(S/S(0))", "fontsize", 3, "color", "blue");
a = gca()

disp(xopt)









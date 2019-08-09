// Initial condition Toluene_Cyclohexane_294K 
clear

//path='/Users/christianpantoja/Documents/opt/20181019-Toluene_Cyclohexane_data/'
path='/home/jul/data/diffusion/2019-07-18-ciclohexano_tolueno/'
//path='/run/media/jul/44B6-1E16/data/jwist/nmr/ref_ledbpgp2s_1d_2d/'

// folder inicial 
n1=62

// Number of total folders
nt=17

// incremental time dv
//dv=fscanfMat('/Users/christianpantoja/Documents/Backup_20180520_data/Backup_June_2018/pulse_program_T1/vdlist1ss.txt')

// Matrix processing 
N=linspace(n1,n1+nt-1,nt)

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

mclose('all')
plot2d(Spectrum(9800:10600,:))
signal_1_int_reg_low=9800
signal_1_int_reg_upp=10600

for i=1:size(N,2)
I1(1,i)=sum(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i)));
end

plot2d(Spectrum(9000:9800,:))
signal_2_int_reg_low=9000
signal_2_int_reg_upp=9800

plot2d(Spectrum(3700:4600,:))
signal_2_int_reg_low=3700
signal_2_int_reg_upp=4600

for i=1:size(N,2)
I2(1,i)=sum(real(Spectrum(signal_2_int_reg_low:signal_2_int_reg_upp,i)));
end


sp=linspace(-36,36,nt)


OMEG_low=-36000 // lutidine_2,6
OMEG_upp=+36000 //

OMEGA=linspace(OMEG_upp,OMEG_low,nt)
r=42.57747892*10.^2
z= OMEGA./(r*10.6104) + 2.0

z=z(1,$:-1:1)


//for i=1:19
//    X_T(1,i)=I2(1,i)*1/5/((I2(1,i))*1/5+(I1(1,i))*1/12);
//    X_C(1,i)=I1(1,i)*1/12/((I2(1,i))*1/5+(I1(1,i))*1/12);
//end

// Using methyl group of Toluene 


for i=1:nt
    X_T(1,i)=I2(1,i)*1/3/((I2(1,i))*1/3+(I1(1,i))*1/12);
    X_C(1,i)=I1(1,i)*1/12/((I2(1,i))*1/3+(I1(1,i))*1/12);
end

X_T0=cat(1,z,X_T)
X_C0=cat(1,z,X_C)

X_T0=X_T0'
X_C0=X_C0'


plot2d(X_T0(:,1),X_T0(:,2),style=-5)
plot2d(X_C0(:,1),X_C0(:,2),style=-9)




//fprintfMat('/Users/christianpantoja/Desktop/backup_oct_2018/Cyclohexane_Toluene_system/X_T0.txt',X_T0) 
//fprintfMat('/Users/christianpantoja/Desktop/backup_oct_2018/Cyclohexane_Toluene_system/X_C0.txt',X_C0) 


for i=1:nt
    plot2d([1:16384]+100*i, Spectrum(:,i)+1e5*i)
end


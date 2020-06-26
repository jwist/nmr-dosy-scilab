//Condición inicial Ciclohexano_Tolueno

clear
clc

//Dirección de la carpeta donde se encuentran los experimentos

path='/Users/Jarvis/Desktop/Datosnmr/2019-07-18-ciclohexano_tolueno/'

//Carpeta inicial

nb=170

//Carpeta final 

nf=849

//Numero de puntos

nt=17

//Numero de mapeos

nm=40

//BUCLE DE MAPEOS//

M=linspace(nb,nb+nm-1,nm)

//Carpeta que hace el seguimineto a las carpetas iniciales de cada mapeo

n1=nb

for j=1:size(M,2)

if n1<nf then 

//Matriz de procesamiento

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

//Limites de integracion

//------------------------------------------------------------------------------

//Ciclohexano_Para_Aromatico
//plot2d(Spectrum(9900:11200,:))
signal_1_int_reg_low=9900
signal_1_int_reg_upp=11200

for i=1:size(N,2)
I1(1,i)=sum(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,i)));
end

//Tolueno_Aromatico
//plot2d(Spectrum(3800:4900,:))
signal_2_int_reg_low=3800
signal_2_int_reg_upp=4900

for i=1:size(N,2)
I2(1,i)=sum(real(Spectrum(signal_2_int_reg_low:signal_2_int_reg_upp,i)));
end

//------------------------------------------------------------------------------

//Ciclohexano_Para_Metilo
//plot2d(Spectrum(9900:11200,:))
signal_4_int_reg_low=9900
signal_4_int_reg_upp=11200

for i=1:size(N,2)
I4(1,i)=sum(real(Spectrum(signal_4_int_reg_low:signal_4_int_reg_upp,i)));
end

//Tolueno_Metilo
//plot2d(Spectrum(9000:9900,:))
signal_3_int_reg_low=9000
signal_3_int_reg_upp=9900

for i=1:size(N,2)
I3(1,i)=sum(real(Spectrum(signal_3_int_reg_low:signal_3_int_reg_upp,i)));
end

//------------------------------------------------------------------------------

//Tratamiento espacial a la coordenada z

OMEG_low=-36000
OMEG_upp=+36000

OMEGA=linspace(OMEG_upp,OMEG_low,nt)
r=42.57747892*10.^2
z= OMEGA./(r*10.6104) + 2.0

z=z(1,$:-1:1)

//Tramiento para encontrar la fraccion molar a partir de las integrales

for i=1:nt
    X_CA(1,i)=I1(1,i)*1/12/((I2(1,i))*1/5+(I1(1,i))*1/12);
    X_A(1,i)=I2(1,i)*1/5/((I2(1,i))*1/5+(I1(1,i))*1/12);
    X_M(1,i)=I3(1,i)*1/3/((I3(1,i))*1/3+(I4(1,i))*1/12);
    X_CM(1,i)=I4(1,i)*1/12/((I3(1,i))*1/3+(I4(1,i))*1/12);
    X_CT(1,i)=(((I1(1,i)*1/12/((I2(1,i))*1/5+(I1(1,i))*1/12))+(I4(1,i)*1/12/((I3(1,i))*1/3+(I4(1,i))*1/12)))/2);
    X_TT(1,i)=(((I2(1,i)*1/5/((I2(1,i))*1/5+(I1(1,i))*1/12))+(I3(1,i)*1/3/((I3(1,i))*1/3+(I4(1,i))*1/12)))/2);
end

X_CA0=cat(1,z,X_CA)
X_A0=cat(1,z,X_A)
X_CM0=cat(1,z,X_CM)
X_M0=cat(1,z,X_M)
X_CT0=cat(1,z,X_CT)
X_TT0=cat(1,z,X_TT)

X_CA0=X_CA0'
X_A0=X_A0'
X_CM0=X_CM0'
X_M0=X_M0'
X_CT0=X_CT0'
X_TT0=X_TT0'

//Graficas

//------------------------------------------------------------------------------

//plot2d(X_CA0(:,1),X_CA0(:,2),style=-11)//Ciclohexano_Para_Aromatico (Cuadro)
//plot2d(X_A0(:,1),X_A0(:,2),style=-5)//Tolueno_Aromatico (Diamante)

//------------------------------------------------------------------------------

//plot2d(X_CM0(:,1),X_CM0(:,2),style=-9)//Ciclohexano_Para_Metilo (Circulo)
//plot2d(X_M0(:,1),X_M0(:,2),style=-3)//Tolueno_Metilo(Circulo_Cruz)

//------------------------------------------------------------------------------

plot2d(X_CT0(:,1),X_CT0(:,2),style=-12)//Ciclohexano (Triangulo)
plot2d(X_TT0(:,1),X_TT0(:,2),style=-8)//Tolueno (Rombo cruz)

//------------------------------------------------------------------------------

//Solapamiento de espectros

//plot2d(Spectrum(:,j))
//for i=1:nt
//    plot2d([1:16384], Spectrum(:,i))
//end
//Angulo
//+100*i
//+1e5*i

//------------------------------------------------------------------------------

n1=n1+nt

end

end

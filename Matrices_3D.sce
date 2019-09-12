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

//Hora del primer experimento___Mirar acqus

Dia=8
Horas=18
Minutos=20
Segundos=39
Milisegundos=017

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

//Tramiento para encontrar la fraccion molar a partir de las integrales

for i=1:nt
    X_CA(1,i)=I1(1,i)*1/12/((I2(1,i))*1/5+(I1(1,i))*1/12);
    X_A(1,i)=I2(1,i)*1/5/((I2(1,i))*1/5+(I1(1,i))*1/12);
    X_M(1,i)=I3(1,i)*1/3/((I3(1,i))*1/3+(I4(1,i))*1/12);
    X_CM(1,i)=I4(1,i)*1/12/((I3(1,i))*1/3+(I4(1,i))*1/12);
    X_CT(1,i)=(((I1(1,i)*1/12/((I2(1,i))*1/5+(I1(1,i))*1/12))+(I4(1,i)*1/12/((I3(1,i))*1/3+(I4(1,i))*1/12)))/2);
    X_TT(1,i)=(((I2(1,i)*1/5/((I2(1,i))*1/5+(I1(1,i))*1/12))+(I3(1,i)*1/3/((I3(1,i))*1/3+(I4(1,i))*1/12)))/2);
end

//Tratamiento espacial a la coordenada z

OMEG_low=-36000
OMEG_upp=+36000

OMEGA=linspace(OMEG_upp,OMEG_low,nt)
r=42.57747892*10.^2
z= OMEGA./(r*10.6104) + 2.0

z=z(1,$:-1:1)

//Tratamiento para el tiempo

   for i=1:size(N,2)
        idr = mopen(path+string(N(1,i))+'/acqus','r');
        procs = mgetl(idr,-1); mclose(idr);
        [w,r]=grep(procs,'$$ 2019-08-');
        [a,b,c,d]=regexp(procs(w),'/(?P<name>\w+)\:(?P<digit>\D*\d+)/');
        [e,f,g,h]=regexp(procs(w),'/(?P<name>\w+)(?P<digit>\D*\d+)/');
        [k,l,m,nnn]=regexp(h(7),'/\:(?P<digit>\D*\d+)/');
        fa = strtod(h(6));
        fe = strtod(h(3));
        fi = strtod(nnn(1));
        fo = strtod(h(4));
        fu = strtod(h(2));
        if fu==Dia then
            time(:,i)=((fa*3600)+(fe*60)+fi+(fo/1000)-((Horas*3600)+(Minutos*60)+Segundos+(Milisegundos/1000)))
        else
            time(:,i)=((fa*3600)+(fe*60)+fi+(fo/1000)+((24*3600)-((Horas*3600)+(Minutos*60)+Segundos+(Milisegundos/1000))))
        end
   end

//------------------------------------------------------------------------------

X_CA0=cat(1,z,time,X_CA)
X_A0=cat(1,z,time,X_A)
X_CM0=cat(1,z,time,X_CM)
X_M0=cat(1,z,time,X_M)
X_CT0=cat(1,z,time,X_CT)
X_TT0=cat(1,z,time,X_TT)

X_CA0=X_CA0'
X_A0=X_A0'
X_CM0=X_CM0'
X_M0=X_M0'
X_CT0=X_CT0'
X_TT0=X_TT0'

n1=n1+nt

end

//------------------------------------------------------------------------------


Y_CA0z(:,j)=X_CA0(:,1)
Y_CA0t(:,j)=X_CA0(:,2)
Y_CA0x(:,j)=X_CA0(:,3)

Y_A0z(:,j)=X_A0(:,1)
Y_A0t(:,j)=X_A0(:,2)
Y_A0x(:,j)=X_A0(:,3)

//------------------------------------------------------------------------------

Y_CM0z(:,j)=X_CM0(:,1)
Y_CM0t(:,j)=X_CM0(:,2)
Y_CM0x(:,j)=X_CM0(:,3)

Y_M0z(:,j)=X_M0(:,1)
Y_M0t(:,j)=X_M0(:,2)
Y_M0x(:,j)=X_M0(:,3)

//------------------------------------------------------------------------------

Y_CT0z(:,j)=X_CT0(:,1)
Y_CT0t(:,j)=X_CT0(:,2)
Y_CT0x(:,j)=X_CT0(:,3)

Y_TT0z(:,j)=X_TT0(:,1)
Y_TT0t(:,j)=X_TT0(:,2)
Y_TT0x(:,j)=X_TT0(:,3)

//------------------------------------------------------------------------------

end

//------------------------------------------------------------------------------

Cadenaz_CA0=[]
Cadenat_CA0=[]
Cadenax_CA0=[]

Cadenaz_A0=[]
Cadenat_A0=[]
Cadenax_A0=[]

//------------------------------------------------------------------------------

Cadenaz_CM0=[]
Cadenat_CM0=[]
Cadenax_CM0=[]

Cadenaz_M0=[]
Cadenat_M0=[]
Cadenax_M0=[]

//------------------------------------------------------------------------------


Cadenaz_CT0=[]
Cadenat_CT0=[]
Cadenax_CT0=[]

Cadenaz_TT0=[]
Cadenat_TT0=[]
Cadenax_TT0=[]

//------------------------------------------------------------------------------


for j=1:size(M,2)
Cadenaz_CA0=cat(1,Cadenaz_CA0,Y_CA0z(:,j));
Cadenat_CA0=cat(1,Cadenat_CA0,Y_CA0t(:,j));
Cadenax_CA0=cat(1,Cadenax_CA0,Y_CA0x(:,j));
Cadenaz_A0=cat(1,Cadenaz_A0,Y_A0z(:,j));
Cadenat_A0=cat(1,Cadenat_A0,Y_A0t(:,j));
Cadenax_A0=cat(1,Cadenax_A0,Y_A0x(:,j));
end

//------------------------------------------------------------------------------

for j=1:size(M,2)
Cadenaz_CM0=cat(1,Cadenaz_CM0,Y_CM0z(:,j));
Cadenat_CM0=cat(1,Cadenat_CM0,Y_CM0t(:,j));
Cadenax_CM0=cat(1,Cadenax_CM0,Y_CM0x(:,j));
Cadenaz_M0=cat(1,Cadenaz_M0,Y_M0z(:,j));
Cadenat_M0=cat(1,Cadenat_M0,Y_M0t(:,j));
Cadenax_M0=cat(1,Cadenax_M0,Y_M0x(:,j));
end

//------------------------------------------------------------------------------

for j=1:size(M,2)
Cadenaz_CT0=cat(1,Cadenaz_CT0,Y_CT0z(:,j));
Cadenat_CT0=cat(1,Cadenat_CT0,Y_CT0t(:,j));
Cadenax_CT0=cat(1,Cadenax_CT0,Y_CT0x(:,j));
Cadenaz_TT0=cat(1,Cadenaz_TT0,Y_TT0z(:,j));
Cadenat_TT0=cat(1,Cadenat_TT0,Y_TT0t(:,j));
Cadenax_TT0=cat(1,Cadenax_TT0,Y_TT0x(:,j));
end

//------------------------------------------------------------------------------

Matriz_CA0=cat(2,Cadenaz_CA0,Cadenat_CA0,Cadenax_CA0)
Matriz_A0=cat(2,Cadenaz_A0,Cadenat_A0,Cadenax_A0)

//------------------------------------------------------------------------------

Matriz_CM0=cat(2,Cadenaz_CM0,Cadenat_CM0,Cadenax_CM0)
Matriz_M0=cat(2,Cadenaz_M0,Cadenat_M0,Cadenax_M0)

//------------------------------------------------------------------------------

Matriz_CT0=cat(2,Cadenaz_CT0,Cadenat_CT0,Cadenax_CT0)
Matriz_TT0=cat(2,Cadenaz_TT0,Cadenat_TT0,Cadenax_TT0)

//------------------------------------------------------------------------------

//plot3d(Matriz_CA0(:,3),Matriz_CA0(:,1),Matriz_CA0(:,2))
//plot3d(Matriz_A0(:,3),Matriz_A0(:,1),Matriz_A0(:,2))

//plot3d(Matriz_CM0(:,3),Matriz_CM0(:,1),Matriz_CM0(:,2))
//plot3d(Matriz_M0(:,3),Matriz_M0(:,1),Matriz_M0(:,2))

plot3d(Matriz_CT0(:,3),Matriz_CT0(:,1),Matriz_CT0(:,2))
plot3d(Matriz_TT0(:,3),Matriz_TT0(:,1),Matriz_TT0(:,2))

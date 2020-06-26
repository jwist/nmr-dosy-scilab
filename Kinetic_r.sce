// 06 Nov / 2018 

// Cyclohexane-Toluene *** Cinetica_2_Replica ***

clear
mclose('all')

//<< El archivo fun_read.sce contiene las herramientas necesarias para extraer 
//los valores de spoffs0 como los respectivos espectros procesados a partir de 
//Topspin >>  
exec("/Users/christianpantoja/Documents/new_processing_nmr/fun_read.sce");

// << Se carga la dirección correspondiente a los archivos destinados a procesar >>
//path='/Users/christianpantoja/Documents/opt/20181018-Toluene_Cyclohexane_data/'
path='/Users/christianpantoja/Documents/opt/20181104-Toluene_Cyclohexane_L_R1/'
// Establecimiento de los rangos de procesamiento ***

// Numero de mapeos espaciales $ 
np=30

// Numero de puntos para cada mapeo $
np_m=19

// folder inicial $
fold_in=43

// folder final $
fold_fi=612

//  Limites de integración inferior **int_reg_low** y superior **int_reg_upp** signal_1 $
//signal_1_int_reg_low=10000
//signal_1_int_reg_upp=11000

signal_1_int_reg_low=10400
signal_1_int_reg_upp=10700

//  Limites de integración inferior **int_reg_low** y superior **int_reg_upp** signal_2 $
//signal_2_int_reg_low=3800
//signal_2_int_reg_upp=4600

//signal_2_int_reg_low=9650
//signal_2_int_reg_upp=9900

signal_2_int_reg_low=9500
signal_2_int_reg_upp=10000

// El vector xs debe contener los valores iniciales asociados al punto inicial de
// cada uno de los mapeos espaciales 

xs=linspace(fold_in,fold_fi-(np_m-1),np)

// El vector xs debe contener los valores finales asociados al punto final de
// cada uno de los mapeos espaciales 
 
ys=linspace(fold_in+(np_m-1),fold_fi,np) 


//<< Este ciclo permite cargar e integrar en la regíon previamente seleccionada
// para cada unas de las señales correspondientes >>

for qw=1:np 

folder_ini = xs(1,qw)
folder_end = ys(1,qw)
jump_folder=1

GH=(folder_ini:jump_folder:folder_end);//list the folders


// ****** Acquisiton and Processing parameters ******

list_par = ['BF1','SW_h','O1','TD']; 
idr=mopen(path+string(GH(1,1))+'/acqus','r');
acqus = mgetl(idr,-1); mclose(idr);

for j = 1:size(list_par,2)
    [w,r]=grep(acqus,'##$'+list_par(j)+'=');
    [a,b,c,d]=regexp(acqus(w),'/(?P<name>\w+)=\s+(?P<digit>\D*\d+\.*\d*)/');
    execstr(d(1)+'='+d(2));
end

O1P = O1/BF1;
DW = 1/(2*SW_h);
AQ = DW*TD;
SW = SW_h/BF1;
ppm = -1*linspace(-SW/2,SW/2,TD/2)+O1P;

[Spectrum,GH] = R_data(path,GH);

// Ordenar el set espectral en un modo particular, (-1) o (+1)
GH = Order(GH,-1);
 

// << Esta ciclo se encarga de integrar en las regiones previamente seleccionadas >>

for i=1:np_m
I1(1,i)=sum(real(Spectrum(signal_1_int_reg_low:signal_1_int_reg_upp,GH(3,i))))// Water 
I2(1,i)=sum(real(Spectrum(signal_2_int_reg_low:signal_2_int_reg_upp,GH(3,i))))// TEA_CH2_Signal
//I3(1,i)=sum(real(Spectrum(signal_3_int_reg_low:signal_3_int_reg_upp,GH(3,i))))// TEA_CH3_Signal
end

// << Almacenamiento de la integral obtenida >> 

A(qw,:)=I1
B(qw,:)=I2
//C(qw,:)=I3

// << La dirección correspondiente al lugar donde se guardaran los datos debe ser especificada >>

//save('/Users/christianpantoja/Documents/Maestria_Datos_Tesis/cinetica_6_278.15_rep/integrales/I1_1','A')
//load('/Users/christianpantoja/Documents/Maestria_Datos_Tesis/cinetica_6_278.15_rep/integrales/I1_1','A')
//save('/Users/christianpantoja/Documents/Maestria_Datos_Tesis/cinetica_6_278.15_rep/integrales/I2','B')
//load('/Users/christianpantoja/Documents/Maestria_Datos_Tesis/cinetica_6_278.15_rep/integrales/I2','B')
//save('/Users/christianpantoja/Documents/Maestria_Datos_Tesis/cinetica_6_278.15_rep/integrales/I3','C')
//load('/Users/christianpantoja/Documents/Maestria_Datos_Tesis/cinetica_6_278.15_rep/integrales/I3','C')

end


////<< En la siguientes lineas la intesidad extraida de los datos se convierte a fracción molar usando 
//// las relaciones de area y numero de protones correspondiente a cada integral >>




//for j=1:np
//for i=1:np_m
//    X_T(j,i)=B(j,i)*1/5/((B(j,i))*1/5+(A(j,i))*1/12);
//    X_C(j,i)=A(j,i)*1/12/((B(j,i))*1/5+(A(j,i))*1/12);
//end
//end


for j=1:np
for i=1:np_m
    X_T(j,i)=B(j,i)*1/3/((B(j,i))*1/3+(A(j,i))*1/12);
    X_C(j,i)=A(j,i)*1/12/((B(j,i))*1/3+(A(j,i))*1/12);
end
end
//
//
////<En las siguientes lineas la coordenada espacial es calculada usando la intensidad del pulse field 
////magnetic gradient G = 10.6104 gauss/cm2 
//
//
//// $ Los valores de spoffs0 deben ser suministrados en esta parte 
//
OMEG_low=-36000 // 
OMEG_upp=+36000 //
//
//
OMEGA=linspace(OMEG_upp,OMEG_low,np_m)
r=42.57747892*10.^2
z= OMEGA./(r*10.6104) + 2.0
//
//
////<< Estos vectores se encargan de la unificación de los datos en la nueva matrix data >>
//
x=[1:np_m:((fold_fi-fold_in)+1)]
y=[np_m:np_m:((fold_fi-fold_in)+1)]
//
//
//// << La matrix de tiempo es cargada y normalizada en las siguientes lineas >>
//
//load('/Users/christianpantoja/Desktop/Cyclohexane_Toluene_system/time_toluene_cyclohexane','G1')  
//load('/Users/christianpantoja/Desktop/backup_oct_2018/Cyclohexane_Toluene_system/time_toluene_cyclohexane_2','G1') 
//load('/Users/christianpantoja/Desktop/backup_oct_2018/Cyclohexane_Toluene_system/time_toluene_cyclohexane_cinetica_2_L','G1') 

//load('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/time_toluene_cyclohexane_cinetica_2_L_10_pts','G1') // 10_pts
//load('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/time_toluene_cyclohexane_cinetica_2_L_15_pts','G1') // 15_pts
//load('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/time_toluene_cyclohexane_cinetica_2_L_20_pts','G1') // 20_pts
//load('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/time_toluene_cyclohexane_cinetica_2_L_25_pts','G1') // 15_pts
//load('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/time_toluene_cyclohexane_cinetica_2_L_30_pts','G1') // 30_pts
load('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/R1_20181104_TC/time_toluene_cyclohexane_cinetica_2_R_30_pts','G1') // 30_pts
 
//
G1=G1-(min(G1)-1);
G1=G1+1403
// 
for j=1:np
    M1(1,x(1,j):y(1,j))=z(1,:);
    M2(1,x(1,j):y(1,j))=G1(j,:);
    M3(1,x(1,j):y(1,j))=X_T(j,:)
end
//
data=cat(1,M1,M2,M3)
//
data=data'
//
//fprintfMat('/Users/christianpantoja/Desktop/backup_oct_2018/Cyclohexane_Toluene_system/toluene_cyclohexane_cinetica_2.txt',data)
//fprintfMat('/Users/christianpantoja/Desktop/backup_oct_2018/Cyclohexane_Toluene_system/toluene_cyclohexane_cinetica_2_LO.txt',data) // Cinetica_"2"_Lorraine 1.0---0.8 (Tolueno)

//fprintfMat('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/toluene_cyclohexane_cinetica_2_L_10pts.txt',data) // 10_pts_Cinetica
//fprintfMat('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/toluene_cyclohexane_cinetica_2_L_15pts.txt',data) // 15_pts_Cinetica
//fprintfMat('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/toluene_cyclohexane_cinetica_2_L_20pts.txt',data) // 20_pts_Cinetica
//fprintfMat('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/toluene_cyclohexane_cinetica_2_L_25pts.txt',data) // 25_pts_Cinetica
//fprintfMat('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/toluene_cyclohexane_cinetica_2_L_30pts.txt',data) // 30_pts_Cinetica



// *** Kinetic_2_replica *** 

fprintfMat('/Users/christianpantoja/Desktop/Tolueno_Ciclohexano_cinetica_2_294K/R1_20181104_TC/toluene_cyclohexane_cinetica_2_R_30pts.txt',data) // 25_pts_Cinetica


//// Calculate initial conditions and save file txt.
//
//// New Matrix with initial condition
//
//M(:,1)=data(1:23,1)  
//M(:,2)=data(1:23,3) 
//
//// save txt. initial conditions TEA
//fprintfMat('/Users/christianpantoja/Documents/Maestria_Datos_Tesis/cinetica_6_278.15_rep/dataset_XT0.txt',M) 
//
//M(:,2)=1-data(1:23,3) 
//
//// save txt. initial conditions H2O
//fprintfMat('/Users/christianpantoja/Documents/Maestria_Datos_Tesis/cinetica_6_278.15_rep/dataset_XT0.txt',M) 


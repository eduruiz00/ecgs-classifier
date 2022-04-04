function [sentido, quad_peaks_arriba,quad_peaks_abajo]=sentido_picos(vect) %se introduce el vector y se enytrega el sentido de los picos del ECG
if ~isvector(vect)
    error('Input must be a vector')
end
L=length(vect);
vect_quad=vect.^2; %se cuadratiza el ECG para incrementar su desviaci�n est�ndar e incrementar los picos.
count=1;
count2=1;
quad_peaks_arriba=[];
quad_peaks_abajo=[];
wind=round(length(vect)*2/60); %se decide la longitud de enventanado para la funci�n window_cov_std m�s abajo explicada
st_quad_fin=window_cov_std(vect_quad,wind); %se aplica la funci�n al vector cuadratizado
 %para comenzar a detectar el sentido, comenzaremos a detectar hacia el signo de la suma de todos los elementos del vector.
for m=2:L-1
    if vect(m)-vect(m-1)>=0 && vect(m+1)-vect(m)<=0 && vect_quad(m)>2*st_quad_fin(m) %c�lculo de picos positivos
        quad_peaks_arriba(count)=m;
        count=count+1;
    end
    if vect(m)-vect(m-1)<=0 && vect(m+1)-vect(m)>=0 && vect_quad(m)>2*st_quad_fin(m) %c�lculo de los picos negativos
        quad_peaks_abajo(count2)=m;
        count2=count2+1;
    end
    
end
if length(quad_peaks_arriba)>length(quad_peaks_abajo) %hacia donde haya m�s picos, ser� el sentido.
        sentido="arriba";
    else
        sentido="abajo";
    end
end
function st_cov_fin=window_cov_std(vect_cov,wind) %Inputs: vector y longitud de ventana. Outputs:vector enventanado
    
st_cov=zeros(1,length(vect_cov)+wind); 
st_cov(1,fix(wind/2):length(vect_cov)+fix(wind/2)-1)=vect_cov;
for n=1:length(vect_cov)
    st_cov_fin(n)=std(st_cov(n:n+wind)); %se calcula la desiaci�n est�ndar para cada punto del vector en el cual est� centrada la ventana.
end
end
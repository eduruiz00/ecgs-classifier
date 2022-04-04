%Outputs: Quad_peaks:vector con los índices de los diferentes picos en
%orden, new_d_peaks: distancia entre cada pico (usada para el clasificador,
%el detector no debería llevar este output.; Inputs: vect: vector del que
%se quieren obtener los picos, sentido: dirección de los picos, se obtiene
%de la función sentido. wind: patrón de correlación, forma caracterísitica
%del pico.

function [quad_peaks,new_d_peaks]=detect_picos_cov(vect,sentido,wind)

if ~isvector(vect)%se comprueba que el input sea un vector. De lo contrario, se lanza un error.
    error('Input must be a vector')
end

L=length(vect);
vect_corr=v_cov(vect, wind); %con la función v_cov, comentada abajo, se obtiene el valor de la correlación con la ventana para cada punto del vector.
count=1;
quad_peaks=[];
threshold=1.6*std(vect_corr); %se establece un threshold mínimo para considerar un pico.
if string(sentido) == "arriba" %el sentido determina el signo de las derivadas alrededor del pico: "Arriba":[+,-], "Abajo":[-,+].
    for m=2:L-1 %para cada punto del vector se evalúa el cambio de tendencia y que este tenga lugar por encima del threshold.
        if vect(m)-vect(m-1)>=0 && vect(m+1)-vect(m)<=0 && vect_corr(m)>threshold
            quad_peaks(count)=m;
            count=count+1; %la variable count se actualiza cada vez que se detecta un pico.
        end
    end
end
if string(sentido) == "abajo"
    for m=2:L-1
        if vect(m)-vect(m-1)<=0 && vect(m+1)-vect(m)>=0 && vect_corr(m)<-threshold
            quad_peaks(count)=m;
            count=count+1;
        end
    end
end

d_peaks=[0 0];
if length(quad_peaks)>2 %se evalúa que haya  más de un pico picos calcuar las distancias entre ellos.
    for g=1:length(quad_peaks)-1
        d_peaks(g)=quad_peaks(g+1)-quad_peaks(g);
    end
    aver=mean(d_peaks);
    count2=0;
    for h=1:length(d_peaks)
        if d_peaks(h)<aver/5 %se eliminan los picos que estén más cerca que 1/5 de la distancia media entre picos (eliminación de outliers demasiado cercanos). 
            quad_peaks(h-count2)=[];
            count2=count2+1;
        else
            quad_peaks(h-count2)=quad_peaks(h-count2);
        end
    end
end
new_d_peaks=[]; %se crea el nuevo vector de distancias habiendo eliminado los picos cercanos
for m=1:length(quad_peaks)-1
        new_d_peaks(m)=quad_peaks(m+1)-quad_peaks(m);
end
end

function vect_cov=v_cov(vect,wind_cov) %Inputs: vect:Vector sobre el que se va a desplazar la ventana, wind_cov: ventana con el patrón que se desea correlacionar (en realidad se calcula la covarianza) punto por punto.

L=length(vect);
vect_ad=zeros(1,length(vect)+length(wind_cov)-1);
vect_ad(1,fix(length(wind_cov)/2)+1:L+fix(length(wind_cov)/2))=vect; %se añaden ceros a ambos lados del vector para calcular la covarianza desde el primer punto sin salirnos de los márgenes
vect_cov=zeros(1,L);
L_c=length(wind_cov);

for n=1:L
    mat_c=[ transpose(wind_cov)  transpose(vect_ad(1,n:n+L_c-1)) ];
    c=cov(mat_c);
    vect_cov(n)=c(1,2); %se tiene en cuenta tan solo la covarianza entre el tramo del vector vect y el patrón (las autocovarianzas no).  
end

end

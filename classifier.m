%% generador clasificadores

%los documentos deben estar en la carpeta de training2017 para que cargue

%N = 100; %determine the number of samples
N = 8528; %use this in case you want to charge all of them, it takes
%almost 30 min

%generamos la señal patron para la correlacion
ecg1 = load('A00001').val;
patron = ecg1(4130:4150);
wind_cov=[300 200 200 300];
mat_vect_cov=zeros(200,150);

%list of classes
reference_total = readtable('REFERENCE.csv');
reference = reference_total(:,2);
list_reference = table2array(reference);
list_reference = list_reference(1:N);


list = 1:N;
cad = string(zeros(1, length(list)));
for i = 1:length(list)
    if i<10
        cad(i)="A0000"+string(i);
    else if i>=10 && i<100
        cad(i)="A000"+string(i);
    else if i>=100 && i<1000
        cad(i)="A00"+string(i);
    else if i>=1000
        cad(i)="A0"+string(i);
    end
    end
    end
    end
end

%generate empty vectors
poinc = zeros(1, N);
t = zeros(1, N);
suma_dif = zeros(1, N);
sd = zeros(1, N);
rmssd = zeros(1, N);

for j = 1:length(cad)
    if rem(j, 100) == 0
        disp(j);
    end
    ecg = load(cad(j));
    signal = ecg.val;
    %calculo de picos para los clasificadores
    [sentido] = sentido_picos(signal);
    [quad_peaks,d_peaks]=detect_picos_cov(signal, sentido, patron);
    
    rmssd(j) = rms(d_peaks); %descriptro rmssd
    sd(j) = std(d_peaks); %descriptor de la desviacion estándar
    t(j)=length(signal)/300; %tiempo
    frec(j,1)=mean(d_peaks); % descriptor frecuencia picos
    
    
    %descriptor diferencia filtro mediana (al final no lo utilizamos)
    m_wind = 20;
    new_signal = [];
    
    for p = 1:(length(signal)-m_wind) % aplicamos filtro de mediana
        mediana = median(signal(p:p+m_wind));
        new_signal = [new_signal, mediana];
    end
    
    dif_total = 0; %buscamos la diferencia ente ambos, el area
    sig_int = signal((m_wind/2)+1:end-(m_wind/2));
    [sentido2] = sentido_picos(sig_int);
    [quad_peaks2,d_peaks2]=detect_picos_cov(sig_int, sentido2, patron);
    longitud = 0;
    for h = 1:length(quad_peaks2)-1
        tramo = sig_int((quad_peaks2(h)+8):(quad_peaks2(h+1)-8));
        longitud = longitud + length(tramo);
        new_sig_tram = new_signal((quad_peaks2(h)+8):(quad_peaks2(h+1)-8));
        dif = abs(tramo-new_sig_tram);
        dif_total = dif_total + sum(dif);
    end
    suma_dif(j) = dif_total/longitud;
    
    
    %descriptor dispersion poincare
    ejex = [];
    ejey = [];
    for i = 1:length(d_peaks)-1
        ejex = [ejex, d_peaks(i)];
        ejey = [ejey, d_peaks(i+1)];
    end
    mx = mean(ejex);
    my = mean(ejey);
    distances = sqrt((ejex-mx).^2 + (ejey-my).^2);
    av_dist = mean(distances);
    poinc(j) = av_dist;
    
    
    
    %descriptor count_zer y count_der
    L=length(d_peaks);
    vect_ad=zeros(1,length(d_peaks)+length(wind_cov)-1);
    vect_ad(1,fix(length(wind_cov)/2)+1:L+fix(length(wind_cov)/2))=d_peaks;
    vect_cov=zeros(1,L);
    L_c=length(wind_cov);
    for n=1:L
        mat_c=[ transpose(wind_cov)  transpose(vect_ad(1,n:n+L_c-1)) ];
        c=cov(mat_c);
        vect_cov(n)=c(1,2);
    end
    vect_d_cov=diff(vect_cov);
    count_der(j,1)=length(find(abs(vect_d_cov)<250));
    count_zer(j,1)=length(find(abs(vect_cov)<800));

end
disp('DONE!')

%% classifier
num_training = 200; %indicar el numero de muestras en training
%tener cuidado, en funcion de la señal puede sobrepasar el valor, mirar
%idx_noise

idx_noise = [];
idx_normal = [];
idx_arr = [];
idx_other = [];
int_reference = [];

%calcular los indices de cada clase
for i = 1:length(list_reference)
    if strcmpi(list_reference(i),'N')
        idx_normal = [idx_normal, i];
        int_reference = [int_reference, 1];
    elseif strcmpi(list_reference(i), '~')
        idx_noise = [idx_noise, i];
        int_reference = [int_reference, 2];
    elseif strcmpi(list_reference(i), 'A')
        idx_arr = [idx_arr, i];
        int_reference = [int_reference, 3];
    else
        idx_other = [idx_other, i];
        int_reference = [int_reference, 4];
    end
end

%realizamos una elección aleatoria de los indices
normal_train = randsample(length(idx_normal),num_training);
normal_test = setdiff(1:length(idx_normal), normal_train);
noise_train = randsample(length(idx_noise),num_training);
noise_test = setdiff(1:length(idx_noise), noise_train);
arr_train = randsample(length(idx_arr),num_training);
arr_test = setdiff(1:length(idx_arr), arr_train);
other_train = randsample(length(idx_other),num_training);
other_test = setdiff(1:length(idx_other), other_train);

%determinamos los valores de training y test
training = [idx_normal(normal_train), idx_noise(noise_train), idx_arr(arr_train)];
training = [training, idx_other(other_train)];
test = [idx_normal(normal_test), idx_noise(noise_test), idx_arr(arr_test)];
test = [test, idx_other(other_test)];
ref_train = int_reference(training);
ref_test = int_reference(test);


%calculo de los descriptores para training y test
suma_train = suma_dif(training);
suma_test = suma_dif(test);
poinc_train = poinc(training);
poinc_test = poinc(test);
rmssd_train = rmssd(training);
rmssd_test = rmssd(test);
count_der_train = count_der(training)';
count_der_test = count_der(test)';
count_zer_train = (count_zer(training).^2)';
count_zer_test = (count_zer(test).^2)';
frec_train = frec(training)';
frec_test = frec(test)';

%lista descriptores
descriptors = [count_der_train; poinc_train; frec_train; count_zer_train; rmssd_train];

num_nei = 6; %numero de vecinos

pre = [count_der_test; poinc_test; frec_test; count_zer_test; rmssd_test];



%hacemos el knn
mdl = fitcknn(descriptors',ref_train,'NumNeighbors',num_nei);
knn = predict(mdl, pre');

%hacemos el lda
MdlLinear = fitcdiscr(descriptors',ref_train);
lda = predict(MdlLinear,pre');


%hacemos la confusion matrix
[C,order] = confusionmat(ref_test, knn);
confusionchart(C, order);
diag = trace(C);
num = sum(sum(C));
disp(diag/num); %calculamos la eficacia

clear; clc;

% executar para todos os possíveis
files = dir('../outputs/*.out');
for i=1:size(files,1)
    filename = strcat('../outputs/',files(i).name(1:end-4),'.out');
    % recebe dados gerados pelo PCA
    A = csvread(filename);
    % converte para uint8
    B = uint8(A);
    % escreve como imagem
    img = strcat('../results/',files(i).name(1:end-4),'.jpg');
    imwrite(B, img);
    disp(['Arquivo ' img ' gerado']);
    % deletar output do PCA
    %delete(filename);
end
 
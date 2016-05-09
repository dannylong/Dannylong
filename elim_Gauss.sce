clear;
// Geração de números aleatórios de 1 a 6
//fix(6*rand()+1);
A = [6 2 -1; 2 4 1; 3 2 8];
b = [7 7 13];

function Ab = gauss(A, b)
    Ab = [A b']; //Matriz aumentada Ab
    [linhas colunas] = size(A);
    for j = 1:colunas
        pivo  = Ab(j, j);        
        // PIVOTAMENTO PARCIAL
        //Pivotamento
        maiorLinha = j;
        for i = (j+1):linhas
            //procurar a linha do maior elemento
            if (abs(Ab(i,j)) > abs(Ab(maiorLinha, j))) then
                maiorLinha = i;
            end
        end
        //Permutar a linha do pivo(j) com a linha maiorLinha
        aux = Ab(j,:);
        Ab(j,:) = Ab(maiorLinha, :);
        Ab(maiorLinha, :) = aux;
        
        //Atualizar o pivo após o pivotamento
        pivo = Ab(j,j)        
        for i = (j+1):linhas
            //zerar todos os termos abaixo do pivo
            //zerar o elemento Ab(i,j) 
            Ab(i,:) = Ab(i, :) - (Ab(i,j)/pivo)*Ab(j,:);
        end
    end
endfunction

// Tratamento da Matriz aumentada "triangular superior"
function x = resolMatTrSupAumentada(Ab)
    [linhas colunas] = size(Ab);
    x = (1:linhas)'; // Vetor solução, inicializado com uma sequência de 1 à ordem da Matriz A
    for i = linhas:-1:1
        pivo = Ab(i, i);
        x(i) = Ab(i, colunas)/pivo;
        for j = 1:(colunas-1);
            if (j<>i) then
                x(i) = x(i) - x(j)*Ab(i, j)/pivo;
            end
        end
    end    
endfunction

// Tratamento da Matriz aumentada "triangular inferior"
function x = resolMatTrInfAumentada(Ab)
    [linhas colunas] = size(Ab);
    x = (1:linhas)'; // Vetor solução, inicializado com uma sequência de 1 à ordem da Matriz A
    for i = 1:linhas    
        pivo = Ab(i, i);
        x(i) = Ab(i, colunas)/pivo;
        for j = 1:(colunas-1);
            if (j<>i) then
                x(i) = x(i) - x(j)*Ab(i, j)/pivo;
            end
        end
    end    
endfunction

function x = resolverSistemaLinear(A, b)
    x = resolMatTrSupAumentada(gauss(A,b));
endfunction




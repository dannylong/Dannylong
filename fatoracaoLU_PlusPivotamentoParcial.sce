clear;
// Fatoração LU COM pivotamento parcial
function [L, U, P] = fatoracaoLU(A)
    
    [linhas colunas] = size(A);
    L = eye(linhas, colunas);
    
    // Matriz para as permutações de linhas de A e b -> A' = PA = LU; Ly = Pb; Ux = y
    P = eye(linhas, colunas);
    
    for j = 1:colunas
        pivo = A(j,j);
        //Pivotamento
        maiorLinha = j;
        for i = (j+1):linhas
            //procurar a linha do maior elemento
            if (abs(A(i,j)) > abs(A(maiorLinha, j))) then
                maiorLinha = i;
            end
        end
        //Permutar a linha do pivo(j) com a linha maiorLinha
        aux = A(j,:);
        A(j,:) = A(maiorLinha, :);
        A(maiorLinha, :) = aux;
        
        P_n                 = eye(linhas, colunas);
        P_aux               = P_n(j, :);
        P_n(j, :)           = P_n(maiorLinha, :);
        P_n(maiorLinha, :)  = P_aux;
        
        // Atualizar a Matriz P, ao final de cada permutação parcial
        P = P_n * P;
        
        // Atualizar a Matriz L, ao final de cada permutação parcial
        L_aux                   = L(j, 1:j-1);
        L(j, 1:j-1)             = L(maiorLinha, 1:j-1);
        L(maiorLinha, 1:j-1)    = L_aux;
        
        //Atualizar o pivo após o pivotamento
        pivo = A(j,j)        
        for i = (j+1):linhas            
            lambda = A(i,j)/pivo;
            //zerar todos os termos abaixo do pivo (elemento genérico A(i, j))
            A(i,:) = A(i, :) - lambda * A(j,:);
            L(i, j) = lambda;
        end

    end
    // No final(para um SL de ordem 3 sem pivotamento, e. g.), L = [1 0 0; lambda21 1 0; lambda31 lambda32 1];
    
    // ----------------------------------------------------------------------------------------------------------------------------    
    
    // U recebe todos os termos da matriz A modificada (matriz triangular superior), mas deve ser corrigida posteriormente, pois...
    // erros de arredondamente podem deixar termos abaixo dos pivôs muito pequenos, mas ainda assim não nulos
    U = A;
    // Garantir que a parte à esquerda(ou abaixo, sei lá xD) da diagonal principal tenha todos os termos nulos
    for j = 1:colunas
        for i = (j+1):linhas
            U(i, j) = 0;
        end
    end 
    
endfunction

A = [3 -4 1; 1 2 2; 4 0 -3];
B = [2 1 1 0;4 3 3 1;8 7 9 5;6 7 9 8];
C = [2 3;1 7];
// Exemplo de uso da função acima no console: [L U P] = fatoracaoLU(B)

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

// 5ª questão da lista de Sistemas Lineares
function Minv = inversaComLUPP(A)
    [l c]       = size(A);
    I           = eye(l, c);
    [mL mU mP]  = fatoracaoLU(A);
    Minv        = A;
    
    for (j = 1:c)
        e_n             = I(:, j);
        y_n             = resolMatTrInfAumentada([mL mP*e_n]);
        x_n             = resolMatTrSupAumentada([mU y_n]);
        Minv(:, j)      = x_n;
    end
endfunction
// #FIM da 5ª questão da lista de Sistemas Lineares

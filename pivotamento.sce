clear;
function Ab = gauss(A, b)
    Ab = [A b];
    [linhas colunas] = size(A);
    for j = 1:colunas
        pivo = Ab(j,j);
        
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
        
        //Atualizar o pivo apos o pivotamento
        pivo = Ab(j,j)        
        for i = (j+1):linhas
            //zerar todos os termos abaixo do pivo
            //zerar o elemento Ab(i,j) 
            Ab(i,:) = Ab(i, :) - (Ab(i,j)/pivo)*Ab(j,:);
        end
    end
endfunction

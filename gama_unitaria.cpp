
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
    DATA_VECTOR(Y);         //observações
    DATA_MATRIX(X);         //matriz de efeito fixo
    PARAMETER_VECTOR(beta); //vetor de parâmetros
    PARAMETER(logphi);      //parâmetro de precisão
    
    Type phi = exp(logphi);
    
    // preditor linear para média
     vector<Type> mu = exp(X*beta)/(1 + exp(X*beta)); // função logito
    //vector<Type> mu = (1-exp(-exp(X*beta)));           // função c. log-log
    vector<Type> tau = pow(mu,1/phi)/(1-pow(mu,1/phi));


    // log-verossimilhança negativa
    Type nll = 0;
    for(int i=0; i < Y.size(); i++)
   
         nll -= phi*log(tau[i]) - lgamma(phi) + (tau[i] - 1)*log(Y[i])+    
                (phi -1)*log(-log(Y[i]));

    // método delta
     ADREPORT(phi);
    
     return nll;    
        
}

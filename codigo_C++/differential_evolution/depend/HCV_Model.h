#ifndef _HCV_Model_H_
#define _HCV_Model_H_

#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>

using namespace std;

const int    TIME  = 50; //tempo de simulacao
const int    AGE   = 500; // idade da infecao (dias)
const int    buffer   = 2;

class HCV_Model{

  private:

    double I[buffer][AGE]; //celulas infectadas

    double Rp[buffer][AGE]; // RNA intracelular positivo
    double Rn[buffer][AGE]; // RNA intracelular negativo
    double Rt[buffer][AGE]; // RNA positivo traduzido

    double T; // celulas alvo
    double V; // virus

    int simCase;
    int days;
    int points;
    double deltaT;
    double deltaA;
    int iterPerDay;
    double tol;
    int vardelta;
    int varrho;

    double  N; //valores iniciais


    double s; // taxa constante de producao de celulas alvo
    double d; // decaimento natural de celulas alvo
    double beta; // taxa de infecao
    double delta; // decaimento de celulas infectadas
    double rho; // taxa de exporta��o de RNA positivo
    double c; // taxa de eliminacaoo do virus pelo SI

    double alpha; // taxa de replica��o de RNA positivo

    double r; // taxa de replica��o de RNA negativo / complexo de replica��o (Rn)

    double k; // coeficiente da fun��o exponencial de atraso na exporta��o de RNA positivo
    double tau; // tempo de atraso para a exporta��o de RNA positivo
    double n; //atraso de delta
    double Rmax; // n�mero m�ximo de RNA negativo / complexo de replica��o (Rn)
    double sigma; // taxa de forma��o de complexos de replica��o
    double mu_t; // decaimento natural de Rt
    double theta; // taxa de disponibilidade para tradu��o
    double mu_c; // decaimento natural de Rc e Rn

    double epsilon_s; // efetividade da terapia em diminuir ou bloquear a exporta��o de RNA positivo
    double epsilon_alpha; // efetividade da terapia em diminuir ou bloquear a replica��o de RNA positivo
    double epsilon_r; // efetividade da terapia em diminuir ou bloquear a replica��o de RNA negativo

    double kappa_t; // fator para aumentar a degrada��o de RNA positivo dispon�vel para tradu��o
    double kappa_c; // fator para aumentar a degrada��o de RNA positivo e negativo no complexo de replica��o

    int saveFiles;
    char *dir;
    FILE* dataInfected;
    FILE* dataVirus;
    FILE* dataTarget;
    FILE* dataRNA_Positivo;
    FILE* dataRNA_Negativo;
    FILE* dataRNA_Traduzido;
    FILE* datadelta;
    FILE* datarho;

    std::string Header();
    std::string Footer(long int t);
    int checkFile(FILE* theFile);
    void initialize();
    double calcIntegral(double vec1[][AGE], double vec2[][AGE], double vec3[][AGE]);
    double calcIntegral2(double a, double b, double vec1[][AGE], double vec2[][AGE], double delta, double rho, double deltaA);
    void update(double vec[][AGE]);

public:
    HCV_Model();
    int solve();

};
/**
* Constructor
*/
HCV_Model::HCV_Model(){

}

/**
* Set conditions and parameter values
*/
void HCV_Model::initialize(){
    fstream param;
    param.open("parametros_DE.txt");
    std::string aux_string;

    getline(param, aux_string, ',');
    double V0 = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    delta = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    mu_t = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    r = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    mu_c = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    epsilon_alpha = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    epsilon_r = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    sigma = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    theta = atof(aux_string.c_str());
    
    getline(param, aux_string, ',');
    rho = atof(aux_string.c_str());

    getline(param, aux_string, ',');
    alpha = atof(aux_string.c_str());

    param.close();

    // cout << "V0: " << V0 << " delta: " << delta << " mu_t: " << mu_t << " r: " << r << " mu_c: " << mu_c << " ep_alpha: " << epsilon_alpha << " epsilon_r: " << epsilon_r << endl;
    /**
    * number of days simulated
    */
    days      = 2;
    /**
    * number of points saved
    */
    points    = 10;
    /**
    * output directory
    */
    //dir       = (char *) "saida/";
    /**
    * simulation
    */
    deltaT     = 0.01; //pow(10,-4);
    /**
    * age
    */
    deltaA     = 0.1;
    /**
    * number of iterations per day
    */
    iterPerDay = 100;

    tol        = pow(10,-4);

    /**
     * Vary delta 1, fix delta 0
     */
    vardelta = 0;//1; // 0;
    /**
     * Vary rho 1, fix delta 0
     */
    varrho = 1;//1; // 0;

    /**
    * baseline parameters
    */
    d     = 0.010;
    s     = 130000;
    // delta = 0.58; //0.47; //PAT8 = 0.58; PAT42 = 0.64; PAT68 = 0.1; PAT69 = 0.47; PAT83 = 0.62; //0.01 ~1.8 // valor estimado 0.14
    beta  = 5*pow(10,-8);
    c     = 22.30;//19
    // rho   = 8.180;//12
    // alpha = 30.0;
    Rmax  = 50.0;
    // r     = 1.49; //PAT8 = 1.49; PAT42 = 1.1; PAT68 = 5.08; PAT69 = 2.24; PAT83 = 1.61;
    tau   = 0.50;
    n     = 1.00;
    k     = 0.80;
    // mu_t  = 0.89; //PAT68 = 0.88; //PAT8 = PAT42 = PAT69 = PAT83 = 0.89;
    // mu_c  = 2.55; //PAT8 = 2.55; PAT42 = 1.72; PAT68 = 3.38; PAT69 = 3.15; PAT83 = 2.39;
    // sigma = 1.30;
    // theta = 1.20; //ou 1.2;
    /**
    * therapy parameters
    */
    // epsilon_alpha = 0.928;  //PAT8 = 0.928; PAT42 = 0.909; PAT68 = 0.992; PAT69 = 0.936; PAT83 = 0.924;
    // epsilon_r     = 0.47; //PAT8 = 0.47; PAT42 = 0.12; PAT68 = 0.61; PAT69 = 0.36; PAT83 = 0.29;
    epsilon_s     = 0.998;
    kappa_t       = 1.00;
    kappa_c       = 1.00;
    
    /**
    * Initial Conditions
    */
    for(int j=0; j<2;j++){
        for(int i=0; i<AGE;i++){
            I[j][i] = 0.0;
            Rt[j][i] = 0.0;
            Rp[j][i] = 0.0;
            Rn[j][i] = 0.0;
        }
    }
    Rt[0][0] = 1.00;


    double rho1;
    if (varrho) rho1 = 0.00;// = rho;
    else rho1 = rho;

    double soma=0;
    soma+=Rp[0][0] + Rt[0][0];

    for(int a=1; a<AGE;a++){
        if(varrho){
            if((deltaA*(double) a)<tau){
                rho1=0.0;
            }else{
                rho1=(1-exp(-k*((deltaA*(double)a)-tau)))*rho;
            }
        }else rho1 = rho;
        Rn[0][a] = (r*Rp[0][a-1] - r*Rp[0][a-1]*(Rn[0][a-1]/Rmax) -
                mu_c*Rn[0][a-1])*deltaA + Rn[0][a-1];

        Rp[0][a] = (alpha*Rn[0][a-1] + sigma*Rt[0][a-1] -
                theta*Rp[0][a-1]- rho1*Rp[0][a-1] -
                mu_c*Rp[0][a-1])*deltaA + Rp[0][a-1];

        Rt[0][a] = (theta*Rp[0][a-1] - sigma*Rt[0][a-1] -
                rho1*Rt[0][a-1] - mu_t*Rt[0][a-1])*deltaA + Rt[0][a-1];

        //if (a<250) soma+=(Rp[0][a] + Rt[0][a]);

        //printf("a = %d Rt = %.4lf Rp = %.4lf Rn = %.4lf \n", a, Rt[0][a], Rp[0][a], Rn[0][a]);
    }
    
    T        = 1.3*pow(10,5);
    V        = V0;//1.07*pow(10,6)//7.15*pow(10,6);///5.64*pow(10,5); //PAT8 = 5.64*pow(10,5), PAT42 = 5.65, PAT68 = 7.5*pow(10,6), PAT69 = 6.14, PAT83 = 5.45
    I[0][0]  = beta*T*V;

    double delta1;
    if (vardelta) delta1 = 0.01;
    else delta1 = delta;
    for(int a=1; a<AGE;a++){
        //delta variavel
        if(vardelta) delta1=(1-exp(-k*(a*deltaA)))*delta;
        else delta1 = delta;
        I[0][a] = (beta*V*T*exp(-delta1*deltaA*(double)a))*deltaA + I[0][a-1];
       // printf("a = %d I0 = %lf \n", a, I[0][a]);
    }
}

/**
 * Updates the current results to position 1
 */
void  HCV_Model::update(double vec[][AGE]){
    for(int a = 0; a < AGE; a++) {
        vec[0][a] = vec[1][a];
	}
}

double  HCV_Model::calcIntegral(double vec1[][AGE], double vec2[][AGE],double vec3[][AGE]){
    int a;
    double sum=0.0;
    for(a = 0; a < AGE; a++){
        sum += (vec1[0][a]+vec2[0][a])*vec3[0][a];
    }
    sum = sum;
    return sum/(2.0*AGE);
}

/******************************************************************************
* Solve model equations
*******************************************************************************/
int HCV_Model::solve(){

    long int t = 0;
    long int a = 0;

    //set initial conditions
    initialize();

    fstream saida;
    saida.open("saida.txt", std::ios::out);

    /**
    * begin time loop
    */
    do{

        // if (t == 0) cout << "Calculating...\n";

        int value = ((int)iterPerDay*days)/points;
        float time_save = (float) t/(float)iterPerDay;
        if(true) {
            // cout << "Saving files : iteration ..."<< time_save << "\n";
            
            saida << time_save << "," << V << endl;

        }
        /**
        * ODEs
        */
        T = (s - d*T - beta*V*T)*deltaT + T;
        V = ((1-epsilon_s)*rho*calcIntegral(I,Rp,Rt) - c*V)*deltaT + V;
        
        //printf("t = %ld V = %lf T = %lf \n", t,V,T);
        /**
        * begin age loop
        */
        double rho1;// = rho;
        double delta1;// = delta;
        if(varrho)
            rho1=0.0;
        else
            rho1=rho;
        a=0;
        Rp[1][a]  = 0.0;
        Rn[1][a]  = 0.0;
        Rt[1][a]  = 1.0;
        I[1][a]   = beta*V*T;

        if (vardelta) delta1=0.01;
        else delta1 = delta;

        for(a = 1; a < AGE-1; a++){
            if(varrho){
                if(((double)a*deltaA)<tau){
                    rho1=0;
                }else{
                    rho1=(1-exp(-k*(((double)a*deltaA)-tau)))*rho;
                }
            }else{
                rho1 = rho;
            }
            //delta variavel
            if(vardelta)
                delta1=delta*(1-exp(-k*(a*deltaA)));
            else
                delta1 = delta;

            // delta variavel
            I[1][a] = (-delta1*I[0][a]-(I[0][a]-I[0][a-1])/(deltaA))*deltaT + I[0][a];

		    Rn[1][a] = ((1-epsilon_r)*r*Rp[0][a] - (1-epsilon_r)*r*Rp[0][a]*(Rn[0][a]/Rmax) -
                        kappa_c*mu_c*Rn[0][a] - (Rn[0][a]-Rn[0][a-1])/(deltaA))*deltaT + Rn[0][a];

		    Rp[1][a] = ((1-epsilon_alpha)*alpha*Rn[0][a]+sigma*Rt[0][a]-theta*Rp[0][a] -
                        (1-epsilon_s)*rho1*Rp[0][a]-kappa_c*mu_c*Rp[0][a] -
                        (Rp[0][a]-Rp[0][a-1])/(deltaA))*deltaT + Rp[0][a];

            Rt[1][a] = (theta*Rp[0][a]-sigma*Rt[0][a]-(1-epsilon_s)*rho1*Rt[0][a] -
                        kappa_t*mu_t*Rt[0][a] - ((Rt[0][a]-Rt[0][a-1])/(deltaA)))*deltaT + Rt[0][a];
	    }

	    update(I);
	    update(Rn);
	    update(Rp);
	    update(Rt);

        t++;

    }while (t < (iterPerDay*days)) ;

    saida.close();

    return 0;
 }

#endif

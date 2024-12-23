//The following ifndef/define/endif pattern is called a 
//scope guard, and prevents the C++ compiler (actually, preprocessor)
//from including a header file more than once.
#ifndef _Measure_H_
#define _Measure_H_

#include <string>
#include <iostream>
#include "itensor/all.h"

using namespace itensor;
class Measure
    {
    private:

    double t,J;
    int Nx,Ny;
    int holespin,holenum,sigmatj;//sigmatj=-1: sigma cdot t-J model
    std::vector<std::string> PST_t;
    int pssize;
    int Ifbasis,Ifsym;
    int IfHt,IfHJ,IfA, IfZi, Ifnh;
    //std::vector<int> irowHt, icolHt, irowHJ, icolHJ, irowA, icolA , irowZi , icolZi;
    //std::vector<std::complex<double>> valHt,valHJ, valA, valZi;
    std::vector<double> theta_ij;
    std::vector<itensor::MPS> ctilde_array_psl;
    MPS phi0;
    tJ sitesnew; //the sites of phi0 with QN
    tJ sitesold; //the sites of basis without QN

    public:
    
    Measure(double t_,double J_,int Nx_, int Ny_, int holespin_,int holenum_,int sigmatj_,std::vector<std::string> PST_t_,int pssize_,int Ifbasis_,int Ifsym_,int IfHt_,int IfHJ_, int IfA_, int IfZi_,int Ifnh_)
      : t(t_), J(J_),Nx(Nx_), Ny(Ny_), holespin(holespin_), holenum(holenum_),sigmatj(sigmatj_),PST_t(PST_t_),pssize(pssize_),
      Ifbasis(Ifbasis_),Ifsym(Ifsym_),IfHt(IfHt_),IfHJ(IfHJ_),IfA(IfA_), IfZi(IfZi_), Ifnh(Ifnh_)
        {  
            generate_theta();
            generate_lattice();
            if(Ifbasis==1) load_basis();
            else generate_basis();

        }

    void generate_theta();
    
    double get_theta(const int & h, const int & l);

    int get_pssign(std::string pstt, unsigned i);

    void ctilde(int site_i, int pssign, const itensor::MPS& phi, itensor::MPS& psi);

    void generate_basis();

    void load_basis();

    void generate_lattice();

    void generalmeasure();

    void measureAred();

    void measureHt();

    void measureHJ();

    void measureZi();

    void measurenh();

    unsigned sym_ope(unsigned h, unsigned sym_index);

    bool Isconj(unsigned sym_index);

    std::complex<double> sym_phase(int m, int n, unsigned sym_index, unsigned hs);

    void coutbasis();

    //void coutdata();

    void write_sparse_mat(std::string filename, std::vector<int> &row, std::vector<int> &col, std::vector<std::complex<double>> &data);

    };



    void Measure::generate_theta()
    {
        if(theta_ij.size() == 0){ 
            for(int iy = 0; iy < 2*Ny - 1; iy++){ 
                for(int ix = 0; ix < 2*Nx - 1; ix++){ 
                    double dx = (double)ix - (double)(Nx - 1), dy =  ((double)iy - (double)(Ny - 1)); 
                    theta_ij.push_back(std::atan2(dy, dx)); 
                }
            }
        }
    }

    double Measure::get_theta(const int & h, const int & l){
        int dy = (h - 1)% Ny - (l - 1)% Ny;
        int dx =  (h - 1)/Ny - (l - 1)/ Ny;
        return theta_ij[(2*Nx - 1)*(dy + (Ny - 1))+ dx + Nx - 1]; //calculate dh- dl;
    } 

    int Measure::get_pssign(std::string pstt, unsigned i){
        if (pstt[i] == '+')
            return 1;
        else if (pstt[i] == '-')
            return -1;
        else if (pstt[i] == '0')
            return 0;
        else
            return 0;
    }

    void Measure::ctilde(int site_i, int pssign, const itensor::MPS& phi, itensor::MPS& psi) {
        auto N = Nx * Ny;
        auto lattice = squareLattice(Nx,Ny);
        psi = phi;
        
        if(pssign != 0){
        for(int j = 1; j <= N; j++)
        {
            if(j==site_i) continue;
            psi.position(j);
            double theta = 0.0;
            theta = pssign * get_theta(site_i, j);
            auto psop = std::complex<double>(std::cos(theta) - 1, std::sin(theta))*op(sitesold,"Ndn",j);
            psop += op(sitesold,"Id",j);
            auto newpsi = psop * psi(j);
            newpsi.noPrime();
            psi.set(j,newpsi);
        }
        
        if(holespin==1)
        {
            for(int j = 1; j <= N; j++)
            {
            if(j==site_i) continue;
            psi.position(j);
            double theta = 0.0;
            theta =  pssign * get_theta(j,site_i);
            auto psop =std::complex<double>(1-std::cos(theta), (-1)*std::sin(theta))*op(sitesold,"Ntot",j);
            psop += std::complex<double>(std::cos(theta), std::sin(theta))*op(sitesold,"Id",j);
            auto newpsi = psop * psi(j);
            newpsi.noPrime();
            psi.set(j,newpsi);
            }
        }

        }
        psi = removeQNs(psi);
        psi.position(site_i);

        if (site_i==1) sitesnew = siteInds(psi);
        auto ania = op(sitesnew,"Aup",site_i);
        if(holespin == 1) ania = op(sitesnew,"Adn",site_i);
        auto newpsi = ania * psi(site_i);
        newpsi.noPrime();
        psi.set(site_i,newpsi);            

        for (int i = site_i - 1 + holespin; i >=1;i--)
        {
            psi.position(i);
            auto jordan = op(sitesnew,"F",i);
            auto newpsi = jordan * psi(i);
            newpsi.noPrime();
            psi.set(i,newpsi);   
        }

        std::cout<<"finish"<<site_i<<std::endl;

    }

    void Measure::generate_basis()
    {
        readFromFile("sites_file",sitesold);
        MPS phi0temp(sitesold);
        readFromFile("psi_file",phi0temp);
        phi0 = phi0temp;
        int N = Nx * Ny;
        for (int i = 1; i <= N ; i++)
        {
            for (int psl = 0; psl < pssize; psl ++)
            {
            MPS psil(sitesold);
            ctilde(i,get_pssign(PST_t[psl],0),phi0,psil);
            ctilde_array_psl.push_back(psil);
            }
        }
    }

    void Measure::load_basis()
    {
        readFromFile("sites_file",sitesold);
        MPS phi0temp(sitesold);
        readFromFile("psi_file",phi0temp);
        phi0 = phi0temp;
        auto f = h5_open("basis.h5",'r');
        std::string basisname = "basis";
        for (int i = 0; i< Nx*Ny ;i++)
        {
            std::string basisnumnew = basisname + std::to_string(i);
            MPS psil= h5_read<MPS>(f,basisnumnew);
            ctilde_array_psl.push_back(psil);
        }
        close(f);
        sitesnew = siteInds(ctilde_array_psl[0]);
        std::cout<<"success_load"<<std::endl;
    }

    void Measure::generate_lattice()
    {
        int N = Nx * Ny;
        std::vector<int> script;
        script.resize(N);
        std::vector<int> index;
        index.resize(N);
        for (int i=0; i< N; i++) script[i]= N;
        int num = 0;
        for (int i=0; i < N; i++)
        {
            int iy = i%Ny, ix = i/Ny;
            script[i] = num;
            index.push_back(i);
            num++;
        }
        int N_h = num;
        std::cout << "hdim = " << N_h << ";" << std::endl;
        std::ofstream scriptfile;
        scriptfile.open("script.dat");
        for (auto & iter:script)
            scriptfile << iter << "\n";
        scriptfile << std::endl;
        scriptfile.close();
            
        std::ofstream indexfile;
        indexfile.open("index.dat");
        for (auto & iter:index)
            indexfile << iter << "\n";
        indexfile << std::endl;
        indexfile.close();
    }

    void Measure::generalmeasure()
    {
        if(IfA == 1) measureAred();
        if(IfHt == 1) measureHt();
        if(IfHJ == 1) measureHJ();
        if(IfZi == 1) measureZi();
        if(Ifnh == 1) measurenh(); 
    }

    void Measure::measureAred()
    {
        auto N = Nx * Ny;
        int sym_num = (Nx == Ny)? (7):(3);
        if (Ifsym==0) sym_num = 0 ;
        std::vector<int>  irowA, icolA;
        std::vector<std::complex<double>>  valA;      
        for (int psr = 0; psr < pssize; psr++){
            for (int psl = 0; psl < pssize; psl++){

                std::vector<std::vector<bool>> label_visited(N,std::vector<bool>(N));
                for (int j = 0; j < N;j++)
                for(int i = 0; i < N; i++)
                    label_visited[i][j] = 0;

                for (int i = 1; i <= N; i++){
                    for (int j = i; j <= N; j++)
                    {
                        if ((label_visited[i-1][j-1])||(label_visited[j-1][i-1])) continue;
                        std::complex<double> Avalue = innerC(ctilde_array_psl[(i-1)*pssize+psl],ctilde_array_psl[(j-1)*pssize+psr]);
                        if(std::abs(Avalue)<1e-15) Avalue=0;
                        for (int sym_index = 0; sym_index < sym_num + 1; sym_index++){
                            int S_i=sym_ope(i,sym_index),S_j=sym_ope(j,sym_index);
                            if ((label_visited[S_i-1][S_j-1])||(label_visited[S_j-1][S_i-1])) continue;
                            int Iscon = Isconj(sym_index);
                            if(S_i==S_j)
                            {
                                irowA.push_back(S_i * pssize + psl); icolA.push_back(S_j * pssize + psr);
                                valA.push_back(std::real(Avalue));
                            }
                            else
                            {
                                irowA.push_back(S_i * pssize + psl); icolA.push_back(S_j * pssize + psr);
                                std::complex<double> Avaluesym =(Iscon)?(std::conj(Avalue)):(Avalue);
                                Avaluesym = Avaluesym * sym_phase(get_pssign(PST_t[psl], 0), get_pssign(PST_t[psr], 0), sym_index, holespin);
                                valA.push_back(Avaluesym);
                                irowA.push_back(S_j * pssize + psr);icolA.push_back(S_i * pssize + psl);
                                valA.push_back(std::conj(Avaluesym));
                            }
                            label_visited[S_i-1][S_j-1] = 1;
                            label_visited[S_j-1][S_i-1] = 1;
                        }
                    }
                    std::cout<<"finishA"<<i<<std::endl;
                }
            }
        }
        write_sparse_mat("A.dat", irowA, icolA, valA);
    }

    void Measure::measureHt()
    {
        auto N = Nx * Ny;
        auto lattice = squareLattice(Nx,Ny);
        int sym_num = (Nx == Ny)? (7):(3);
        if (Ifsym==0) sym_num = 0 ;
        std::vector<int>  irowHt, icolHt;
        std::vector<std::complex<double>>  valHt;  
        auto ampo = AutoMPO(sitesnew);    
        for(auto bnd : lattice)
        {
            ampo += -t,"Cdagup",bnd.s1,"Cup",bnd.s2;
            ampo += -t,"Cdagup",bnd.s2,"Cup",bnd.s1;
            ampo += -sigmatj*t,"Cdagdn",bnd.s1,"Cdn",bnd.s2;
            ampo += -sigmatj*t,"Cdagdn",bnd.s2,"Cdn",bnd.s1;
        }
        auto Ht = toMPO(ampo);
        for (int psr = 0; psr < pssize; psr++){
            for (int psl = 0; psl < pssize; psl++){

                std::vector<std::vector<bool>> label_visited(N,std::vector<bool>(N));
                for (int j = 0; j < N;j++)
                for(int i = 0; i < N; i++)
                    label_visited[i][j] = 0;

                for (int i = 1; i <= N; i++){
                    for (int j = i; j <= N; j++)
                    {
                        if ((label_visited[i-1][j-1])||(label_visited[j-1][i-1])) continue;
                        std::complex<double> Htvalue = innerC(ctilde_array_psl[(i-1)*pssize+psl],Ht,ctilde_array_psl[(j-1)*pssize+psr]);
                        if(std::abs(Htvalue)<1e-15) Htvalue=0;
                        for (int sym_index = 0; sym_index < sym_num + 1; sym_index++){
                            int S_i=sym_ope(i,sym_index),S_j=sym_ope(j,sym_index);
                            if ((label_visited[S_i-1][S_j-1])||(label_visited[S_j-1][S_i-1])) continue;
                            int Iscon = Isconj(sym_index);
                            if(S_i==S_j)
                            {
                                irowHt.push_back(S_i * pssize + psl); icolHt.push_back(S_j * pssize + psr);
                                valHt.push_back(std::real(Htvalue));
                            }
                            else
                            {
                                irowHt.push_back(S_i * pssize + psl); icolHt.push_back(S_j * pssize + psr);
                                std::complex<double> Htvaluesym =(Iscon)?(std::conj(Htvalue)):(Htvalue);
                                Htvaluesym = Htvaluesym * sym_phase(get_pssign(PST_t[psl], 0), get_pssign(PST_t[psr], 0), sym_index, holespin);
                                valHt.push_back(Htvaluesym);
                                irowHt.push_back(S_j * pssize + psr);icolHt.push_back(S_i * pssize + psl);
                                valHt.push_back(std::conj(Htvaluesym));
                            }
                            label_visited[S_i-1][S_j-1] = 1;
                            label_visited[S_j-1][S_i-1] = 1;
                        }
                    }
                    std::cout<<"finishHt"<<i<<std::endl;
                }
            }
        }
        write_sparse_mat("Hteff.dat", irowHt, icolHt, valHt);
    }

    void Measure::measureHJ()
    {
        auto N = Nx * Ny;
        int sym_num = (Nx == Ny)? (7):(3);
        std::vector<int>  irowHJ, icolHJ;
        std::vector<std::complex<double>>  valHJ;
        if (Ifsym==0) sym_num = 0 ;
        auto lattice = squareLattice(Nx,Ny);
        auto ampo = AutoMPO(sitesnew);        
        for(auto bnd : lattice)
        {
            ampo +=  0.5*J,"S+",bnd.s1,"S-",bnd.s2;
            ampo +=  0.5*J,"S-",bnd.s1,"S+",bnd.s2;
            ampo +=      J,"Sz",bnd.s1,"Sz",bnd.s2;
            ampo +=-0.25*J,"Ntot",bnd.s1,"Ntot",bnd.s2;
        }
        auto HJ = toMPO(ampo);
        for (int psr = 0; psr < pssize; psr++){
            for (int psl = 0; psl < pssize; psl++){

                std::vector<std::vector<bool>> label_visited(N,std::vector<bool>(N));
                for (int j = 0; j < N;j++)
                for(int i = 0; i < N; i++)
                    label_visited[i][j] = 0;

                for (int i = 1; i <= N; i++){
                    for (int j = i; j <= N; j++)
                    {
                        if ((label_visited[i-1][j-1])||(label_visited[j-1][i-1])) continue;
                        std::complex<double> HJvalue = innerC(ctilde_array_psl[(i-1)*pssize+psl],HJ,ctilde_array_psl[(j-1)*pssize+psr]);
                        if(std::abs(HJvalue)<1e-15) HJvalue=0;
                        for (int sym_index = 0; sym_index < sym_num + 1; sym_index++){
                            int S_i=sym_ope(i,sym_index),S_j=sym_ope(j,sym_index);
                            if ((label_visited[S_i-1][S_j-1])||(label_visited[S_j-1][S_i-1])) continue;
                            int Iscon = Isconj(sym_index);
                            if(S_i==S_j)
                            {
                                irowHJ.push_back(S_i * pssize + psl); icolHJ.push_back(S_j * pssize + psr);
                                valHJ.push_back(std::real(HJvalue));
                            }
                            else
                            {
                                irowHJ.push_back(S_i * pssize + psl); icolHJ.push_back(S_j * pssize + psr);
                                std::complex<double> HJvaluesym =(Iscon)?(std::conj(HJvalue)):(HJvalue);
                                HJvaluesym = HJvaluesym * sym_phase(get_pssign(PST_t[psl], 0), get_pssign(PST_t[psr], 0), sym_index, holespin);
                                valHJ.push_back(HJvaluesym);
                                irowHJ.push_back(S_j * pssize + psr);icolHJ.push_back(S_i * pssize + psl);
                                valHJ.push_back(std::conj(HJvaluesym));
                            }
                            label_visited[S_i-1][S_j-1] = 1;
                            label_visited[S_j-1][S_i-1] = 1;
                        }
                    }
                    std::cout<<"finishHJ"<<i<<std::endl;
                }
            }
        }
        write_sparse_mat("HJeff.dat", irowHJ, icolHJ, valHJ);
    }

    void Measure::measureZi()
    {
        auto N = Nx * Ny;
        int sym_num = (Nx == Ny)? (7):(3);
        std::vector<int>  irowZi, icolZi;
        std::vector<std::complex<double>>  valZi;
        if (Ifsym==0) sym_num = 0 ;
        //sym_num = 0;
        std::vector<std::string> PST_Zi;
        PST_Zi.push_back("0");
        std::vector<itensor::MPS> barec_array;
        for (int i = 1; i <= N ; i++)
        {
            for (int psl = 0; psl < (int)PST_Zi.size(); psl ++)
            {
            MPS psil(sitesold);
            ctilde(i,get_pssign(PST_Zi[psl],0),phi0,psil);
            barec_array.push_back(psil);
            }
        }

        for (int psl = 0; psl < (int)PST_Zi.size(); psl++){
        for (int psr = 0 ;psr < pssize; psr++){

            std::vector<std::vector<bool>> label_visited(N,std::vector<bool>(N));
            for (int j = 0; j < N;j++)
                for(int i = 0; i < N; i++)
                    label_visited[i][j] = 0;

            for (int i = 1; i <= N; i++)
            {   
                for (int j = 1; j <= N; j++)
                {   
                    if (label_visited[i-1][j-1]) continue;
                    std::complex<double> Zivalue = innerC(barec_array[(i-1)*(int)PST_Zi.size()+psl],ctilde_array_psl[(j-1)*pssize+psr]);
                    if(std::abs(Zivalue)<1e-15) Zivalue=0;    
                    for (int sym_index = 0; sym_index < sym_num + 1; sym_index++){
                        int S_i=sym_ope(i,sym_index),S_j=sym_ope(j,sym_index);
                        if (label_visited[S_i-1][S_j-1]) continue;
                        int Iscon = Isconj(sym_index);
                        std::complex<double> Zivaluesym =(Iscon)?(std::conj(Zivalue)):(Zivalue);
                        Zivaluesym = Zivaluesym * sym_phase(get_pssign(PST_Zi[psl], 0), get_pssign(PST_t[psr], 0), sym_index, holespin);
                        irowZi.push_back(S_i * (int)PST_Zi.size() + psl); icolZi.push_back(S_j * pssize + psr);
                        valZi.push_back(Zivaluesym);
                        label_visited[S_i-1][S_j-1] = 1;                        
                    }
                }
                std::cout<<"finishZi"<<i<<std::endl;
            }
        } 
        }
        write_sparse_mat("Zi.dat", irowZi, icolZi, valZi);  
    }

    void Measure::measurenh()
    {
        auto N = Nx * Ny;
        int sym_num = (Nx == Ny)? (7):(3);
        std::vector<int>  irownh, icolnh , inumhole;
        std::vector<std::complex<double>>  valnh;
        if (Ifsym==0) sym_num = 0 ;
        //sym_num = 0;
        for (int psl = 0; psl < pssize; psl++){
        for (int psr = 0 ;psr < pssize; psr++){

            std::vector<std::vector<std::vector<bool>>> label_visited(N,std::vector<std::vector<bool>>(N,std::vector<bool>(N)));
            for (int j = 0; j < N;j++)
                for(int i = 0; i < N; i++)
                    for(int k = 0; k < N; k++)
                        label_visited[i][j][k] = 0;
            for (int h = 1 ;h <= N ;h++)
            {
                auto nhole_h = (-1)*op(sitesnew,"Ntot",h);
                nhole_h += op(sitesnew,"Id",h);

            for (int i = 1; i <= N; i++)
            {
                for (int j = i; j <= N; j++)
                {   
                    if ((label_visited[i-1][j-1][h-1])||(label_visited[j-1][i-1][h-1])) continue;
                    MPS nhstate = ctilde_array_psl[(j-1)*pssize+psr];
                    auto newstate = nhole_h * nhstate(h);
                    newstate.noPrime();
                    nhstate.set(h,newstate);
                    std::complex<double> nhvalue = innerC(ctilde_array_psl[(i-1)*pssize+psl],nhstate);
                    if(std::abs(nhvalue)<1e-15) nhvalue=0;    
                    for (int sym_index = 0; sym_index < sym_num + 1; sym_index++){
                        int S_i=sym_ope(i,sym_index),S_j=sym_ope(j,sym_index),S_h=sym_ope(h,sym_index);
                        if ((label_visited[S_i-1][S_j-1][S_h-1])||(label_visited[S_j-1][S_i-1][S_h-1])) continue;
                        int Iscon = Isconj(sym_index);
                        std::complex<double> nhvaluesym =(Iscon)?(std::conj(nhvalue)):(nhvalue);
                        nhvaluesym = nhvaluesym * sym_phase(get_pssign(PST_t[psl], 0), get_pssign(PST_t[psr], 0), sym_index, holespin);
                        irownh.push_back(S_i * pssize + psl); icolnh.push_back(S_j * pssize + psr); inumhole.push_back(S_h);
                        if(S_i == S_j) valnh.push_back(std::real(nhvaluesym));
                        else
                        {
                            valnh.push_back(nhvaluesym);
                            irownh.push_back(S_j * pssize + psr); icolnh.push_back(S_i * pssize + psl); inumhole.push_back(S_h);
                            valnh.push_back(std::conj(nhvaluesym));
                        }
                        label_visited[S_i-1][S_j-1][S_h-1] = 1;label_visited[S_j-1][S_i-1][S_h-1] = 1;                       
                    }
                }
            }
            std::cout<<"finishnh"<<h<<std::endl;
            }
        } 
        }
        std::ofstream nhfile;
        nhfile.open("nh.dat");
        for (int i = 0; i< (int) irownh.size();i++)
        nhfile << irownh[i]  << "\t" << icolnh[i]  << "\t" << inumhole[i]  << "\t" << valnh[i].real() << std::showpos << valnh[i].imag() << "i"<<std::endl;
        nhfile << std::endl;
        nhfile.close();  
    }

    //sym_index: 0 I; 1 Px; 2 Py; 3 Parity; 4 C41; 5 C43; 6 Pxy; 7 Pyx;  Px means two points are symmetric along y axis.
    unsigned Measure::sym_ope(unsigned h, unsigned sym_index){
        static int N = Nx * Ny;
        switch(sym_index){
            case 0 :
                return h;
            case 1 :
                return N - Ny +2 + 2*((h - 1)% Ny) - h;
            case 2 :
                return Ny - 1 -2 * ((h - 1)% Ny) + h;
            case 3 :
                return N - h + 1;
            case 4 :
                return N - (-(unsigned)(( h - 1) /Ny) + Ny* ( ( h - 1) % Ny) + Ny) + 1; 
            case 5 :
                return -(unsigned)(( h - 1) /Ny) + Ny* ( ( h - 1) % Ny) + Ny;
            case 6 :
                return Ny * ((h - 1)% Ny) + (unsigned)(( h - 1) /Ny) + 1;
            case 7 :
                return N - (Ny * ((h - 1)% Ny) + (unsigned)(( h - 1) /Ny)  ) ;
            default :
                return h;
        }
    }


    bool Measure::Isconj(unsigned sym_index){
        switch(sym_index){
            case 0 :
                return 0;
            case 1 :
                return 1;
            case 2 :
                return 1;
            case 3 :
                return 0;
            case 4 :
                return 0;
            case 5 :
                return 0;
            case 6 :
                return 1;
            case 7 :
                return 1;
            default :
                return 0;
        }
    }    

    //sym_index: 0 I; 1 Px; 2 Py; 3 Parity; 4 C41; 5 C43; 6 Pxy; 7 Pyx;
    std::complex<double> Measure::sym_phase(int m, int n, unsigned sym_index, unsigned hs){
        static int N = Nx * Ny; 
        static const double pi = 3.141592653589793 ;
        double theta;
        int psl = m, psr= n;
        int dnspinnum;
        if((hs==0)&&(holenum%2==1)) dnspinnum=(holenum+1)/2;
        else dnspinnum=holenum/2;
        switch(sym_index){
            case 0 :
                theta = 0;
                break;
            case 1 :
                theta = (psr - psl) * pi * (N / 2 - hs -dnspinnum) ;
                break;
            case 2 :
                theta = 0;
                break;
            case 3 :
                theta = - (psr - psl) * pi * (N / 2 - hs -dnspinnum) ;
                break;
            case 4 :
                theta = (psr - psl) * pi / 2 * (N / 2 - hs -dnspinnum) ;
                break;
            case 5 :
                theta = - (psr - psl) * pi / 2 * (N / 2 - hs -dnspinnum) ;
                break;
            case 6 :
                theta = (psr- psl) * (N / 2 - hs -dnspinnum) * pi / 2 ;
                break;
            case 7 :
                theta = - (psr- psl) * (N / 2 - hs -dnspinnum) * pi / 2 ;
                break;
            default :
                theta = 0;
                break;
        }
        return std::complex<double>(std::cos(theta),std::sin(theta));
    }

    void Measure::coutbasis()
    {
        auto f = h5_open("basis.h5",'w');
        std::string basisname = "basis";
        for (int i = 0; i< (int)ctilde_array_psl.size() ;i++)
        {
            std::string basisnumnew = basisname + std::to_string(i);
            h5_write(f,basisnumnew,ctilde_array_psl[i]);
        }
        close(f);
    }

    /*
    void Measure::coutdata()
    {
        if(IfA == 1) write_sparse_mat("A.dat", irowA, icolA, valA);
        if(IfHt == 1) write_sparse_mat("Hteff.dat", irowHt, icolHt, valHt);
        if(IfHJ == 1) write_sparse_mat("HJeff.dat", irowHJ, icolHJ, valHJ); 
        if(IfZi == 1) write_sparse_mat("Zi.dat", irowZi, icolZi, valZi);       
    }*/

    void Measure::write_sparse_mat(std::string filename, std::vector<int> &row, std::vector<int> &col, std::vector<std::complex<double>> &data)
    {
        std::ofstream afile;
        afile.open(filename);
        for (size_t i = 0; i < col.size(); i++)
            afile << row[i]  << "\t" << col[i]  << "\t" << data[i].real() << std::showpos << data[i].imag() << "i"<<std::endl;
        afile.close();
    }


#endif //_Measure_H_

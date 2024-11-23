#include "itensor/all.h"
#include <time.h>
#include "Measure.h"

using namespace itensor;

int 
main(int argc, char* argv[])
    {
    if(argc < 2) 
    { 
    printfln("Usage: %s input_file",argv[0]); 
    return 0; 
    }
    auto input = InputGroup(argv[1],"input");
    auto Nx = input.getInt("Nx");  auto Ny = input.getInt("Ny");
    auto t = input.getReal("t",3.); auto J = input.getReal("J",1.);
    auto N = Nx * Ny;
    int Ifphi0 = input.getInt("Ifphi0"); int holenum = input.getInt("holenum");
    int sweepnum = input.getInt("sweepnum"); int sweepmaxdim = input.getInt("sweepmaxdim"); auto sweepcutoff = input.getReal("sweepcutoff",1E-16);//sweepcutoff = 1E-16;
    int Ifbasis = input.getInt("Ifbasis"); int Ifmeasure = input.getInt("Ifmeasure");
    int holespin = input.getInt("holespin");int measurenum = input.getInt("measurenum");int sigmatj = input.getInt("sigmatj");
    int Phasestring = input.getInt("Phasestring");
    int Ifsym = input.getInt("Ifsym"); int IfHt = input.getInt("IfHt"); int IfHJ = input.getInt("IfHJ"); int IfA = input.getInt("IfA");
    int IfZi = input.getInt("IfZi"); int Ifnh = input.getInt("Ifnh");

    if((Ifphi0==1)&&(Ifbasis==1)) std::abort();

    auto lattice = squareLattice(Nx,Ny);
    auto sites = tJ(N,{"ConserveQNs=",true});

    std::time_t start_time ,  time2 , timegenerate;
    std::time(&start_time);
    /*
    printfln("\nt = %.2f\nJ = %.2f\nNx = %.2f\nNy = %.2f\nNtot = %.2f",t,J,Nx,Ny,N);
    for(auto& bnd : lattice) // bnd is of type LatticeBond
    {
    printfln("Bond from site %d -> %d",bnd.s1,bnd.s2);
    printfln("  Connecting points (%s,%s) -> (%s,%s)",bnd.x1,bnd.y1,bnd.x2,bnd.y2);
    printfln("  This bond is of type \"%s\"",bnd.type);
    }
    */

    if (Ifphi0 != 0)
    {
    auto ampo = AutoMPO(sites);
    for(auto bnd : lattice)
    {
        ampo += -t,"Cdagup",bnd.s1,"Cup",bnd.s2;
        ampo += -t,"Cdagup",bnd.s2,"Cup",bnd.s1;
        ampo += -t,"Cdagdn",bnd.s1,"Cdn",bnd.s2;
        ampo += -t,"Cdagdn",bnd.s2,"Cdn",bnd.s1;
    }
    auto Ht = toMPO(ampo);
    for(auto bnd : lattice)
    {
        ampo +=  0.5*J,"S+",bnd.s1,"S-",bnd.s2;
        ampo +=  0.5*J,"S-",bnd.s1,"S+",bnd.s2;
        ampo +=      J,"Sz",bnd.s1,"Sz",bnd.s2;
        ampo +=-0.25*J,"Ntot",bnd.s1,"Ntot",bnd.s2;
    }
    auto H = toMPO(ampo);

    auto state = InitState(sites);
    std::vector<bool> initvisit(N);
    for (int i = 0 ;i < N; i++ ) initvisit[i] = 0;
    for (int i = 0; i < holenum; i++)
    {
        if(i%2 ==1)
        {state.set(Ny*(i-1)/2 + Ny*(Nx-1)/2,"Emp"); initvisit[Ny*(i-1)/2 + Ny*(Nx-1)/2-1] = 1;}
        if(i%2 ==0)
        {state.set(Ny*(i)/2 + Ny*(Nx-1)/2+1,"Emp"); initvisit[Ny*(i)/2 + Ny*(Nx-1)/2] = 1;}
    }
    for(int i = 1; i <= N; ++i) 
        {
        if(initvisit[i-1]) continue;
        if(((i-1)%Ny+(int)((i-1)/Ny))%2 == 0)
            state.set(i,"Up");
        else
            state.set(i,"Dn");
    }
    auto psi0 = MPS(state);
    println("\nTotal QN of phi0 = ",totalQN(psi0));

    /*
    std::vector<std::string> PST_t;
    PST_t.push_back("+");
    Measure measure_halffilling(t,J,Nx,Ny,0,1000,1,PST_t,1,0,0,0,1);
    for (int j= 2 ;j<=6; j++){
    auto amponew = AutoMPO(sites);
    double theta = measure_halffilling.get_theta(1,j);
    amponew += std::complex<double>(std::cos(theta) - 1, std::sin(theta)), "Ndn", j;
    amponew += 1, "Id", j;
    auto psop = toMPO(amponew);
    psi0test = applyMPO(psop,psi0test,{"Method=","DensityMatrix","MaxDim=",300,"Cutoff=",1E-16});
    psi0test.noPrime();
    }
    */
    auto sweeps = Sweeps(sweepnum); //number of sweeps is 5
    //sweeps.maxdim() = 500,600,700,800,900,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000,1000;
    sweeps.maxdim() = sweepmaxdim;
    sweeps.cutoff() = sweepcutoff;
    println(sweeps);

    auto [energy,psi] = dmrg(H,psi0,sweeps,"Quiet");
    printfln("\nGround State Energy = %.10f",energy);
    double Etvalue = inner(psi,Ht,psi);
    printfln("\nE_t = %.10f\nE_J = %.10f", Etvalue , energy-Etvalue );
    println("\nTotal QN of Ground State = ",totalQN(psi));
    writeToFile("sites_file",sites);
    writeToFile("psi_file",psi);
    }


    //to test the Lz of the phi0 state    
    
    tJ sitesold;
    readFromFile("sites_file",sitesold);
    MPS phi0test(sitesold);
    readFromFile("psi_file",phi0test);
    //MPS phih0rot(sitesold);
    auto teststate1=InitState(sitesold);
    auto teststate2=InitState(sitesold);
    for(int i = 1; i <= N; ++i) 
        {
        int roti=N - (-(int)(( i - 1) /Ny) + Ny* ( ( i - 1) % Ny) + Ny) + 1;
        if((i<holenum + Ny*(Nx-1)/2)&&(i>= + Ny*(Nx-1)/2)) {teststate1.set(i,"Emp");teststate2.set(roti,"Emp");continue;}
        if(((i-1)%Ny+(int)((i-1)/Ny))%2 == 1)
            {teststate1.set(i,"Up");teststate2.set(roti,"Up");}
        else
            {teststate1.set(i,"Dn");teststate2.set(roti,"Dn");}
    }
    auto psitest1 = MPS(teststate1);
    auto psitest2 = MPS(teststate2);
    println("\nTotal QN of Ground State = ",totalQN(psitest1));
    println("\nTotal QN of Ground State = ",totalQN(phi0test));
    std::cout<<innerC(psitest1,phi0test)<<std::endl;
    std::cout<<innerC(psitest2,phi0test)<<std::endl;
    std::cout<<innerC(psitest1,phi0test)/innerC(psitest2,phi0test)<<std::endl; //exp(iLz*pi/2)
    

    if(Ifmeasure!=0){

    std::vector<std::string> PST_t;
    if (Phasestring==1) PST_t.push_back("+");
    if (Phasestring==-1) PST_t.push_back("-");
    if (Phasestring==0) PST_t.push_back("0");
    Measure measure_halffilling(t,J,Nx,Ny,holespin,holenum,sigmatj,PST_t,PST_t.size(),Ifbasis,Ifsym,IfHt,IfHJ,IfA,IfZi,Ifnh);
    std::time(&timegenerate);
    std::time_t elapsedgenerate = timegenerate - start_time;
    std::printf("\ntimegenerate = %lds",elapsedgenerate);
    std::printf("\ntimegenerate = %.4fmin",double(elapsedgenerate)/(60));
    measure_halffilling.generalmeasure();
    //measure_halffilling.coutdata();
    //measure_halffilling.coutbasis();
    }

    std::time(&time2);
    std::time_t elapsed = time2 - start_time;
    std::printf("\ntime = %lds",elapsed);
    std::printf("\ntime = %.4fmin",double(elapsed)/(60));
    return 0;
    }



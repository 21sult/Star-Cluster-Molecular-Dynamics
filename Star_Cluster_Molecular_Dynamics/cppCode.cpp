#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>

using namespace std;

// Constantes globais~
double kg_in_sm = (1.9891 * pow(10,-30)); // 1 kg em massas solares
double km_in_pc = (3.24078 * pow(10,-14)); // 1 km em parsecs
double s_in_yr = (3.17098 * pow(10,-8)); // 1 segundo em anos 
double r_t = 5.0; // pc
double G = 4.49353036*pow(10,-15); // pc+3 Msolar-1 ano-2 
bool Print = false; // printar resultados

class Body
{
  public:
    vector<double> r = vector<double>(3,0);
    vector<double> v = vector<double>(3,0);
    vector<double> a = vector<double>(3,0);
    
    vector<double> r_prev = vector<double>(3,0);
    vector<double> v_prev = vector<double>(3,0);
    vector<double> a_prev = vector<double>(3,0);
    double m;
};

void printFrame(int frameNo, int nBodies, vector<Body> bodies, vector<double> T, vector<double> V, vector<double> E)
{
    cout << "###### FRAME " << frameNo+1 << " ###### " << endl;
    for (int i = 0; i < nBodies; i++)
    {
        cout << "Body #" << i+1 << ":" << endl;
        cout << "r_x = " << bodies[i].r[0] << endl;
        cout << "r_y = " << bodies[i].r[1] << endl;
        cout << "r_z = " << bodies[i].r[2] << endl;
        cout << "v_x = " << bodies[i].v[0] << endl;
        cout << "v_y = " << bodies[i].v[1] << endl;
        cout << "v_z = " << bodies[i].v[2] << endl;
        cout << "a_x = " << bodies[i].a[0] << endl;
        cout << "a_y = " << bodies[i].a[1] << endl;
        cout << "a_z = " << bodies[i].a[2] << endl;
    }
    cout << "T = " << T[frameNo] << " V = " << V[frameNo] << " E = " << E[frameNo] << endl << endl;
}

// Funções globais
vector<double> vecXscalar(vector<double> v, double k)
// Devolve todos os elementos de um vetor v vecXscalariplicados por um escalar k
{
    int n = sizeof(v);
    vector<double> vNew(n);
    for (int i = 0; i < n; i++) vNew[i] = v[i]*k;

    return vNew;
}

vector<double> vecSum(vector<double> v1, vector<double> v2)
// Devolve um vetor v = v1 + v2
{
    int dim = 3;
    int i = 0;
    vector<double> v(dim);
    for (i = 0; i < dim; i++) v[i] = v1[i] + v2[i];
    
    return v;
}

double mag(vector<double> v)
// Devolve a magnitude de um vetor tridimensional
{
    double mag = abs(sqrt( pow(v[0],2) + pow(v[1],2) + pow(v[2],2) ));
    return mag;
}

vector<double> dist(Body body1, Body body2)
// Devolve a distância do corpo 1 ao corpo 2 na forma de um vetor
{
    double x = body1.r[0] - body2.r[0]; // se < 0: corpo 1 à esquerda do corpo 2. se > 0: corpo 1 à direita do corpo 2.
    double y = body1.r[1] - body2.r[1]; // se < 0: corpo 1 abaixo do corpo 2. se > 0: corpo 1 acima do corpo 2.
    double z = body1.r[2] - body2.r[2];
    
    vector<double> d;
    d.push_back(x);
    d.push_back(y);
    d.push_back(z);
    
    return d;
}

vector<double> position(Body body, double dt)
// Devolve o vetor posição (x,y,z) de um corpo
{
    vector<double> adt2 = vecXscalar(body.a, pow(dt, 2));
    double rx = 2*body.r[0] - body.r_prev[0] + adt2[0];
    double ry = 2*body.r[1] - body.r_prev[1] + adt2[1];
    double rz = 2*body.r[2] - body.r_prev[2] + adt2[2];
    
    vector<double> r;
    r.push_back(rx);
    r.push_back(ry);
    r.push_back(rz);

    return r;
}

vector<double> velocity(Body body, double dt)
// Devolve o vetor velocidade (vx, vy, vz) de um corpo
{
    vector<double> adt = vecXscalar(body.a, dt);
    double vx = body.v[0] + adt[0];
    double vy = body.v[1] + adt[1];
    double vz = body.v[2] + adt[2];
    
    vector<double> v;
    v.push_back(vx);
    v.push_back(vy);
    v.push_back(vz);
    
    return v;
}

vector<double> acceleration(Body body1, Body body2)
// Devolve a aceleração de um corpo devido a força gravitica de outro corpo
{
    double r = mag(dist(body1, body2));
    //if (r<0.1) r=0.1;
    double ax = abs(G*body2.m/pow(r,2)) * (body2.r[0]-body1.r[0])/r;
    double ay = abs(G*body2.m/pow(r,2)) * (body2.r[1]-body1.r[1])/r;
    double az = abs(G*body2.m/pow(r,2)) * (body2.r[2]-body1.r[2])/r;
    
    vector<double> a;
    a.push_back(ax);
    a.push_back(ay);
    a.push_back(az);
    
    return a;
}

double kinetic(Body body) // OK
// Calcula energia cinética de um corpo
{
    double T = body.m*pow(mag(body.v),2)/2;
    return T;
}

double potential(Body body1, Body body2) // OK
// Calcula energia potencial de um sistema binário
{
    double d = mag(dist(body1, body2));
    double V = -G*body1.m*body2.m/d;
    return V;
}

vector<double> centerOfMass(int nBodies, vector<Body> bodies)
// Calcula o centro de massa de uma lista de corpos
{
    vector<double> CM = vector<double>(3,0);
    vector<double> sum = vector<double>(3,0);
    double mTotal = 0;
    for (int i = 0; i < nBodies; i++)
    {
        sum = vecSum(sum, bodies[i].r);
        mTotal += bodies[i].m;
    }
    CM = vecSum(CM, sum);
    CM[0] = CM[0]/mTotal;
    CM[1] = CM[1]/mTotal;
    CM[2] = CM[2]/mTotal;
    
    return CM;
}

vector<double> distToPoint(Body body, vector<double> point)
// Calcula a distancia de um corpo a um ponto no espaço
{
    vector<double> dist;
    dist.push_back(body.r[0] - point[0]);
    dist.push_back(body.r[1] - point[1]);
    dist.push_back(body.r[2] - point[2]);
    
    return dist;
}

vector<double> prevDistToPoint(Body body, vector<double> point)
// Calcula a distancia previa de um corpo a um ponto no espaço
{
    vector<double> dist;
    dist.push_back(body.r_prev[0] - point[0]);
    dist.push_back(body.r_prev[1] - point[1]);
    dist.push_back(body.r_prev[2] - point[2]);
    
    return dist;
}

double lowerThanRT(vector<Body> bodies, vector<double> CM, int nBodies)
// Devolve fração de corpos que têm a distância ao centro de massa menor que r_t.
{
    double lowerThanRT = 0;
    for (int i = 0; i < nBodies; i++)
    {
        double magDist = mag(distToPoint(bodies[i], CM));
        if (magDist < r_t) lowerThanRT += 1;
    }
    
    return (double)lowerThanRT/(double)nBodies;
}

vector<double> updateEnergies(int nBodies, vector<Body> bodies)
// Devolve T, V, E do aglomerado fornecido
{
    double Ti = 0;
    double Vi = 0;
    double Ei = 0;
    // Calcular energias
    for (int i = 0; i < nBodies; i++)
    {
        Ti += kinetic(bodies[i]);
        for (int j = i+1; j < nBodies; j++)
        {
            Vi += potential(bodies[i], bodies[j]);
        }
    }
    Ei = Ti + Vi;
    
    vector<double> energies;
    energies.push_back(Ti);
    energies.push_back(Vi);
    energies.push_back(Ei);
    
    return energies;
}

vector<Body> readList(int nBodies, string filename)
{
    vector<Body> bodyList;
    string text;
    ifstream infile(filename);
    for (int i = 0; i < nBodies; i++)
    {
        getline(infile, text);
        string xStr;
        string yStr;
        string zStr;
        string vxStr;
        string vyStr;
        string vzStr;
        string mStr;
        stringstream s(text);
        s >> xStr >> yStr >> zStr >> vxStr >> vyStr >> vzStr >> mStr;
        
        Body newBody;
        newBody.r[0] = stod(xStr);
        newBody.r[1] = stod(yStr);
        newBody.r[2] = stod(zStr);
        newBody.v[0] = stod(vxStr)*km_in_pc/s_in_yr;
        newBody.v[1] = stod(vyStr)*km_in_pc/s_in_yr;
        newBody.v[2] = stod(vzStr)*km_in_pc/s_in_yr;
        newBody.m = stod(mStr);
        
        bodyList.push_back(newBody);
    }
    infile.close();
    
    return bodyList;
}

double kingDistribution(double r)
// Distribuição de King
{
    double king = (((1/sqrt(1+pow(r,2))) - (1/sqrt(1+pow(5,2)))));
    return king;
}

double kroupaDistribution(double m, double alpha)
{
    double kroupa = pow(m, -alpha);
    return kroupa;
}

double generateMass()
// Devolve uma massa aleatória entre 0 e 1.
{
    bool ok = false;
    double alpha;
    double m;
    while (ok == false)
    {
        double u = drand48();
        m = 0.01 + (100-0.1)*drand48();
        if (m >= 0.01 && m < 0.08)
        {
            alpha = 0.3;
            if (u <= kroupaDistribution(m, alpha)/0.1869450902893) ok = true;
        }
        else if (m >= 0.08 && m < 0.5)
        {
            alpha = 1.3;
            if (u <= kroupaDistribution(m, alpha)/3.0075302995944) ok = true;
        }
        else if (m >= 0.5 && m <= 100)
        {
            alpha = 2.3;
            if (u <= kroupaDistribution(m, alpha)/1.8921361078910) ok = true; 
        }
    }
    return m;
}

double generatePosition(int n0)
// Gera uma posição aleatória segundo perfil de King
{
    bool ok = false;
    double u;
    double x;
    while (ok == false)
    {
        u = drand48();
        x = rand() % 5 + 1;
        if (u <= kingDistribution(x)/1.331860) ok = true;
    }
    
    return u;
}

vector<double> distributeInSphere(double p)
// Distribui aleatoriamente (uniforme) um ponto numa esfera
{
    double theta = rand()% 2*M_PI;
    double phi = acos(2*drand48() - 1);
    
    double x = p*cos(theta)*sin(phi);
    double y = p*sin(theta)*sin(phi);
    double z = p*cos(phi);
    
    vector<double> coords = {x,y,z};
	return coords;
}

vector<double> generateSphereCoords(double n0)
// Gera um ponto aleatório (dist. uniforme) numa esfera
{
    vector<double> coords = distributeInSphere(n0);

	return coords;
}

double generateVelocity()
// Gera uma velocidade aleatória
{
	double conv = pow(10,3)*(3.085678 * pow(10, -16))*3.16887646 * pow(10,+8); //pc / ano
	double y1 = drand48();
	double y2 = drand48();

	double v1 = sqrt(-2*log(y1))*sin(2*M_PI*y2)*conv;
	double v2 = sqrt(-2*log(y1))*cos(2*M_PI*y2)*conv;
    
	double V = (v1+v2)/2;
	return V;
}

void verletSim(string interest, vector<Body> bodies, int nBodies, int nIterations, string filename, double dt)
{
    // Listas e variáveis:
    vector<double> T = vector<double>(nIterations, 0); // lista das energias cinéticas do sistema
    vector<double> V = vector<double>(nIterations, 0); // lista das energias potenciais do sistema
    vector<double> E = vector<double>(nIterations, 0); // lista das energias mecânicas do sistema
    vector<double> energies;
    vector<double> lowerRtList;
    vector<double> CM;
    double d;
    
    ofstream data(filename);
    if (interest == "energies") data << "T" << "\t\t" << "V" << "\t\t" << "E" <<  endl;
    else if (interest == "n_med")  data << "lowerThanRT" << endl;
    
    // Como o algoritmo utiliza t+dt, t e t-dt, precisamos calcular os primeiros 2 frames antes do loop principal.
    // > frame 0:
    // > com as massas e posições fornecidas, podemos calcular as acelerações de cada corpo
    for (int i = 0; i < nBodies; i++)
    {
        for (int j = i+1; j < nBodies; j++)
        {
            vector<double> aNew = acceleration(bodies[i], bodies[j]);
            for (int k = 0; k < 3; k++)
            {
                bodies[i].a[k] += aNew[k];
                bodies[j].a[k] -= aNew[k];
            }
        }
    }
    
    // Calcular energias
    if (interest == "energies")
    {
        energies = updateEnergies(nBodies, bodies);
        T[0] = energies[0];
        V[0] = energies[1];
        E[0] = energies[2];
        data << T[0] << "\t" << V[0] << "\t" << E[0] <<  endl;
    }
    
    // Anotar quantos corpos estão a uma distancia inferior a rt do centro do aglomerado
    if (interest == "n_med")
    {
        CM = centerOfMass(nBodies, bodies);
        lowerRtList.push_back(lowerThanRT(bodies, {0,0,0}, nBodies));
        data << lowerRtList[0] << endl;
    }
    
    // > frame 1:
    for (int i = 0; i < nBodies; i++)
    {
        vector<double> rprev = bodies[i].r;

        bodies[i].r = vecSum( vecSum(bodies[i].r, vecXscalar(bodies[i].v, dt)) , vecXscalar(bodies[i].a, (1/2)*pow(dt,2)) );
        bodies[i].v = velocity(bodies[i], dt);
        
        vector<double> a_total(3,0);
        for (int j = 0; j < nBodies; j++)
        {
            if (i != j) a_total = vecSum(a_total, acceleration(bodies[i], bodies[j]));
        }
        bodies[i].a = a_total;
        
        bodies[i].r_prev = rprev;
    }
    
    // Calcular energias
    if (interest == "energies")
    {
        energies = updateEnergies(nBodies, bodies);
        T[1] = energies[0];
        V[1] = energies[1];
        E[1] = energies[2];
        data << T[1] << "\t" << V[1] << "\t" << E[1] <<  endl;
    }
    
    // Anotar quantos corpos estão a uma distancia inferior a rt do centro do aglomerado
    else if (interest == "n_med")
    {
        CM = centerOfMass(nBodies, bodies);
        lowerRtList.push_back(lowerThanRT(bodies, {0,0,0}, nBodies));
        data << lowerRtList[1] << endl;
    }
    
    // Loop principal
    for (int i = 2; i < nIterations; i++)
    {
        // > calcular prox a(t), v(t), r(t)
        for (int i = 0; i < nBodies; i++)
        {
            vector<double> rprev = bodies[i].r;
    
            bodies[i].r = position(bodies[i], dt);
            bodies[i].v = velocity(bodies[i], dt);
            
            vector<double> a_total(3,0);
            for (int j = 0; j < nBodies; j++)
            {
                if (i != j) a_total = vecSum(a_total, acceleration(bodies[i], bodies[j]));
            }
            bodies[i].a = a_total;
            bodies[i].r_prev = rprev;
        }
        
        // Calcular energias
        if (interest == "energies")
        {
            energies = updateEnergies(nBodies, bodies);
            T[i] = energies[0];
            V[i] = energies[1];
            E[i] = energies[2];
            data << T[i] << "\t" << V[i] << "\t" << E[i] << endl;
        }
        
        // Anotar quantos corpos estão a uma distancia inferior a rt do centro do aglomerado
        else if (interest == "n_med")
        {
            CM = centerOfMass(nBodies, bodies);
            lowerRtList.push_back(lowerThanRT(bodies, {0,0,0}, nBodies));
            data << lowerRtList[i] << endl;
        }
    }
    data.close();
    
    cout << "Simulation ended successfully." << endl;
    // Fim da simulação
}

void sphereSim(vector<Body> bodies, int nBodies, int nIterations, int nConfigs, string filename, double rt, double dt)
{
    vector<double> avgList = vector<double>(nIterations, 0);
    ofstream data(filename);
    data << "t" << "\t" << "N_medio" << endl;
    
    for (int config = 0; config < nConfigs; config++)
    {
        cout << "Config: " << config << endl;
        for (int i = 0; i < nBodies; i++)
        {
            // Gerar primeira posição, velocidade e massa
            double r = generatePosition(nBodies);
            double v = generateVelocity();
            bodies[i].r_prev = generateSphereCoords(r);
            bodies[i].v = generateSphereCoords(v);
            bodies[i].m = generateMass();
        }
        // Gerar segunda posição
        for (int i = 0; i < nBodies; i++) bodies[i].r = vecSum( bodies[i].r_prev, vecXscalar(bodies[i].v, dt) );
        
        // Quantas estrelas começam abaixo de r_t: dois primeiros frames
        int c1 = 0;
        int c2 = 0;
        
        for (int i = 0; i < nBodies; i++)
        {
            vector<double> CM = centerOfMass(nBodies, bodies);
            double r1 = mag(distToPoint(bodies[i], {0,0,0}));
            double r2 = mag(prevDistToPoint(bodies[i], {0,0,0}));
            if (r1 < rt) c1 += 1;
            if (r2 < rt) c2 += 1;
        }
        
        avgList[0] += (double)c1/(double)nBodies;
        avgList[1] += (double)c2/(double)nBodies;
        
        // Evolução do sistema
        for (int it = 2; it < nIterations; it++)
        {
            int inside = 0; // estrelas dentro do raio
            
            // Atualizar acelerações
            for (int i = 0; i < nBodies; i++)
            {
                for (int j = 0; j < nBodies; j++)
                {
                    if (i != j) bodies[i].a = vecSum( bodies[i].a, acceleration(bodies[i], bodies[j]));
                }
            }
            
            // Atualizar posições
            for (int i = 0; i < nBodies; i++)
            {
                vector<double> r_next = position(bodies[i], dt);
                bodies[i].r_prev = bodies[i].r;
                bodies[i].r = r_next;
                
                vector<double> CM = centerOfMass(nBodies, bodies);
                double r = mag(distToPoint(bodies[i], {0,0,0}));
                if (r < rt) inside += 1;
            }
            avgList[it] += (double)inside/(double)nBodies;
        }
    }
    for (int it = 0; it < nIterations; it++) data << dt*it << "\t" << (double)avgList[it]/(double)nConfigs << endl;
    data.close();
}


int main()
{
    srand(777);
    srand48(777);
    
    // Configurações
    bool part1 = true;
    bool part2 = true;
    bool part3 = false;
    double dt; // anos
    int nBodies;
    int nIterations;
    
    // 4.1
    if (part1 == true)
    {
        nIterations = 1000;
        dt = 1*pow(10,+4);
        Body body1;
        Body body2;
        body1.r = {-1, 0, 0}; // pc
        body2.r = {+1, 0, 0};
        body1.v = {0, +2.04542*pow(10,-8), 0}; // pc/ano
        body2.v = {0, -2.04542*pow(10,-8), 0};
        body1.m = 1; // massas solares
        body2.m = 1;
        vector<Body> bodiesBinary = {body1, body2};
        
        verletSim("energies", bodiesBinary, 2, nIterations, "data1.txt" , dt);
    }
    
    // 4.2
    if (part2 == true)
    {
        nIterations = 100;
        dt = 1*pow(10,+5);
        
        // > Lista de 10 corpos
        nBodies = 10;
        vector<Body> bodies10 = readList(nBodies, "10.dat");
        verletSim("n_med", bodies10, nBodies, nIterations, "10output.txt" , dt);
        
        // > Ler lista de 100 corpos
        nBodies = 100;
        vector<Body> bodies100 = readList(nBodies, "100.dat");
        verletSim("n_med", bodies100, nBodies, nIterations, "100output.txt", dt);
        
        // Ler lista de 500 corpos
        nBodies = 500;
        vector<Body> bodies500 = readList(nBodies, "500.dat");
        verletSim("n_med", bodies500, nBodies, nIterations, "500output.txt", dt);
    }
    
    // 4.3
    if (part3 == true)
    {
        dt = 1*pow(10, +5);
        nIterations = 100;
        nBodies = 100;
        int nConfigs = 20;
        double rt = 5; //pc
        vector<Body> bodiesSphere;
        for (int i = 0; i < nBodies; i++)
        {
            Body newBody;
            bodiesSphere.push_back(newBody);
        }
        sphereSim(bodiesSphere, nBodies, nIterations, nConfigs, "sphere_output.txt", rt, dt);
    }
    
    return 0;
}




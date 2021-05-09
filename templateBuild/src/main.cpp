#include <iostream>
#include <fstream>
#include "structs.h"
#include <random>
#include<string>
#include <boost/thread/thread.hpp>
#include <chrono>

// 'biology' parameters & in the parameter sheet
vector<int> ctld = {12500,25000}; // known ctl densities on day 7/8.

static const int buffer =5; //

float tt = static_cast<float>(0.);
float dt = static_cast<float>(1. / (24. * 60.));
float end_t = static_cast<float>(14.);
int frame = 0;

// configThreads populates global vectors with start and end coordinates for
// threads to iterate over. Then returns the number of threads to use.
vector<int> iMinThread; int nThreads;
vector<int> iMaxThread;
int configThreads(const population &p, const grid &g, const parms &par){
    iMinThread.clear();
    iMaxThread.clear();
    // diffusion operates in this range (of i):
    int rmin=g.centre-p.rad-g.buffer;
    int rmax = g.centre + p.rad + g.buffer;
    int range = rmax-rmin;

    // how many threads needed @ 10 slices per thread?
    int minThreads = ceil(float(range)/10.);

    // use user requested threads unless that would be less than 10 slices
    // per thread:
    int nThreads = min(par.Nthreads,minThreads);

    // calculate how many sites each thread is responsible for:
    int sitesPerThread = round(float(range)/(float)(nThreads));

    // calculate start and end i values.

    for(int i = 0; i<nThreads; i++){
        iMinThread.push_back(rmin+sitesPerThread*i);
        iMaxThread.push_back(rmin+sitesPerThread*(1+i));
    }

    return(nThreads);

}

// search neighboring coordinates exhaustively for a free site
int nFreeNeighbours(Coordinate c, grid &g, neighbours &n){
    int nFree=0;
    Coordinate trial;
    for(int i =0; i<26; i++){
        trial.x=c.x+n.n[i].x;
        trial.y=c.y+n.n[i].y;
        trial.z=c.z+n.n[i].z;
        if(!g.isoccupied(trial)) nFree++;
    }
    return(nFree);
};

// here we will write all tumour cells, all CTL and all IFN to file.
// write ti specifically to a special folder in data.
// this function should now take arguments to determine what it writes
void rAll( vector<vector<vector<float>>>& ifn, population &p, tumour &t,int frame, grid&g, parms par, string& filename){

    // CELL ID CODES (VAR) (4 CODES -> 2 BITS...):
    // TUMOUR CELL 0/1
    // CTL 2
    // IFN 3

    // VAL REPRESENTS...
    // IFN (FLOAT)
    // CTL REMAINING DAMAGE (INT/FLOAT...)
    // TUMOUR CELL CONJUGATE (BOOL, INT, 1/0)

    neighbours n; n.initalise();

    list<tumour_cell>::iterator it;
    list<CTL>::iterator it2;
    list<Coordinate>::iterator it3;

    string dir = filename + "output";
    string save = dir + "/" + "slice_" + to_string(frame) + ".csv";

    ofstream s_sum(save.c_str(),std::ofstream::out | std::ofstream::trunc);


    // write out the living tumour cells
    for (it = t.cells.begin(); it != t.cells.end(); ++it){
    if(it->location.z<6 & it ->location.z >(-6)) { // sample grid
        s_sum << it->location.x << ',' << it->location.y << ',' << it->location.z << ',' <<
        0 << ',' << it->health <<  '\n';
    }
    }

    // if CTLs exists write them out along with anything they killed
    if(par.ctl==true){
        for (it3 = t.deadcells.begin(); it3 != t.deadcells.end(); ++it3){
            if(it3->z<6 & it3 ->z >(-6)) { // sample grid
                s_sum <<it3->x << ',' << it3->y   << ',' << it3->z   << ","
                << 1 << ',' << 0   << ',' <<
                 '\n';
            }
        }

        // conjuagted ctl
        for(it2 = p.con.begin(); it2 != p.con.end(); it2++) {
            if (it2->location.z<6 & it2->location.z>(-6)) {
                s_sum << it2->location.x << ',' << it2->location.y << ','
                << it2->location.z <<"," << 2 << ',' << 0 << ',' << '\n';
            }
        }
        //migrating ctl
        for(it2 = p.mig.begin(); it2 != p.mig.end(); it2++){
            if (it2->location.z<6 & it2->location.z>(-6)){
                s_sum << it2->location.x << ',' << it2->location.y   << ','
                << it2->location.z <<  "," << 2 << ',' << 1 << ','  << '\n';
            }

        }

    }

    if(par.ifn==true){
        // try to restrict the slowdown damage by only writing within the known tumour expanse...
        int s=g.centre-ceil(t.radius);
        int e = g.centre + ceil(t.radius);
        for(int i=s; i<e;i++) {
            for (int j = s; j < e; j++) {
                for (int k = s; k < e; k++) {
                    if(k <g.centre+6 & k >(g.centre-6)){
                        if (ifn[i][j][k] > 0.001) { // skip writing pointlessly small floats.
                            s_sum  << i-g.centre << ',' << j-g.centre << ',' << k  - g.centre<< ","
                            << 3 << ',' << (float)ifn[i][j][k]  << '\n';

                        }
                    }

                }
            }
        }
    }

    s_sum.close();

};

// find free site
Coordinate find_free(grid &g, float rad){
    Coordinate trial{};
    bool placed = false;
    int r = static_cast<int>(floor(rad));
    int attempts = 0;

    // randomly generate cells within radius until free site found.
    do{
        trial.x=round(r*dis2(gen));//rand() % (2 *r + 1) - r -1;
        trial.y=round(r*dis2(gen));
        trial.z=round(r*dis2(gen));

        attempts++;

        if(trial.displacement()<r & !g.isoccupied(trial) ) placed = true;

        if(attempts>10000000){
            cout << "warning likely stuck in find_free loop!" << endl;
            placed=true;
        }

    }while(!placed);

    return(trial);

};

// find an occupied grid site
Coordinate find_occupied(grid &g, float rad) {

    Coordinate index{};
    bool occupied;
    do {

        index.x = static_cast<int>((dis(gen) - 0.5) * 2 * rad); // 0.5 + disgen?? what are we searching
        index.y = static_cast<int>((dis(gen) - 0.5) * 2 * rad);  // what about negative values you dolt
        index.z = static_cast<int>((dis(gen) - 0.5) * 2 * rad);

        occupied = g.isoccupied(index);




    } while (!occupied);

    return(index);

};

// add a cell
void add_cell(Coordinate &c,char type, tumour &t, grid &g){
    if((c.displacement()+5)< (g.centre-buffer)){ // keep cells away from the border
        list<tumour_cell>::iterator p;
        p = t.add_cell(c,type);
        g.add(p);

    }
};

// initialise tumour on grid
void startcancer(parms &par, grid &g, tumour &t, ofstream &s_growth){
    char type;
    for(int i=0;i<par.ni;i++){
        Coordinate c = find_free(g,par.initrad);
        if(dis(gen)<0.5 & par.mix==true) type = '4';
        else type='7';
        float tt = 0.;
        add_cell(c,type,t,g);
      //  s_growth<<tt << "," << c.x << "," << c.y << "," << c.z << "," <<type << '\n';
    }
};

// pick a direction and keep going. Return first free grid site found
// WARNING ROUGH!! THINK MORE LATER
Coordinate force_divide(Coordinate c, grid &g){
    Coordinate trial;
    trial=c;

    int dx, dy, dz;
    do{dx=dsite(gen), dy = dsite(gen), dz=dsite(gen);}while(0 == dx & dy == 0 & dz == 0);
    do{
        trial.x+=dx;trial.y+=dy;trial.z+=dz;
        if(trial.displacement()>g.max_permitted_radius) throw "tumour too confined";
    }while(g.isoccupied(ref(trial)));

    return(trial);

};

// perform dipersal attempt
Coordinate trial_disperse(Coordinate c, float numsteps,  grid &g){
    Coordinate cold = c;

    int steps=0;
    auto nsteps = static_cast<int>(floor(numsteps));

    for(;;){

        Coordinate trial_c=c;

        trial_c.x+=dsite(gen);
        trial_c.y+=dsite(gen);
        trial_c.z+=dsite(gen);

         c=trial_c;
         steps++;


        if(nsteps == steps) {
            if(c.displacement()>(g.centre-buffer)) return(cold) ;
            else return(c);
        }
    }
}

// search neighboring coordinates exhaustively for a free site
Coordinate search_neighbours(Coordinate c, grid &g, neighbours &n){
    n.shuffle();
    Coordinate trial{};
    for(int i =0; i<26; i++){
        trial.x=c.x+n.n[i].x;
        trial.y=c.y+n.n[i].y;
        trial.z=c.z+n.n[i].z;

        if(!g.isoccupied(trial)) return(trial);

    }

    return(c);
};



// do all growths in a time step
void growth(parms& par,  float dt, grid &g, tumour &t , neighbours &n, vector<vector<vector<float>>>& ifn1, ofstream &s_growth ) {

    // randomly generate the number of replications and dispersals
    std::poisson_distribution<int> reps(t.cell_count*par.rg*dt);
    bool block =false;
    bool dispersed;

    int nrep = reps(gen);

    Coordinate c,trial;
    for (int i = 0; i < nrep; i++) {

        dispersed = false;
        c = find_occupied(g,t.radius);
        if(par.ifn==true){
            block = ifn1[c.x+g.centre][c.y+g.centre][c.z+g.centre]>par.tifn;
        }

        trial = search_neighbours(c,g,n);


        if (c.x==trial.x & c.y==trial.y & c.z==trial.z) { // this means there were no free neighbouring sites

            if (par.force & !block ){
                trial = force_divide(c,ref(g));
                char type = g.type(c);
                add_cell(trial, type, t, g);
               // s_growth<<tt << "," << trial.x << "," << trial.y << "," << trial.z << "," <<type << '\n';
            }

            if(!par.force & !block & dis(gen)<par.rd){ // try to perform a dispersal attempt
                dispersed = true;
                float r = pow((float) t.cell_count * 3. / (4. * 3.14), (1. / 3.));
                Coordinate trial = trial_disperse(c, static_cast<float>(pow(r, 2.) * par.dist), g);
                if (!g.isoccupied(trial)) {
                    char type = g.type(c);
                    add_cell(trial, type, t, g);
                  //  s_growth<<tt << "," << trial.x << "," << trial.y << "," << trial.z << "," <<type << '\n';

                }
            }else dispersed = false; // only disperse if there were no free neighbours
        }

        if(!dispersed & !g.isoccupied(trial)& !block){
            char type = g.type(c);
            add_cell(trial,type,t,g);
          //  s_growth<<tt << "," << trial.x << "," << trial.y << "," << trial.z << "," <<type << '\n';
        }

    }


}

float convolve(int i, int j, int k, vector<vector<vector<float>>> &kern, vector<vector<vector<float>>> &ifn1){

    float cin=0.;
    float cout=0.;
    float kernval,ifna,ifnb;
    int x,y,z;
    for(int l=-1; l<2;l++) {
        for (int m = -1; m < 2; m++) {
            for (int n = -1; n < 2; n++) {
                kernval=kern[l+1][m+1][n+1];
                x=(i+l);
                y=(j+m);
                z=(k+n);
                ifna=ifn1[x][y][z];
                ifnb=ifn1[i][j][k];
                cin+= (ifna*kernval);
                cout+= (ifnb*kernval);
            }
        }
    }
    float net = ifn1[i][j][k]+cin-cout;
    return(net);
}

void diffuse(grid& g, population&p , vector<vector<vector<float>>> kern, vector<vector<vector<float>>> &ifn1, vector<vector<vector<float>>> &ifn2, int si, int ei){
    // we will split i evenly between threads. Therefore do all j and k iters but only some of i per thread
    // note - believe order is important here, do i{j{k{stuff;}}} never k{j{i{stuff;}}}.
    int sjk=g.centre-p.rad-g.buffer;
    int ejk = g.centre + p.rad + g.buffer;

    float net;
    //cout << si << "," << ei << "," << sjk << "," << ejk << ","
      //   << sjk << "," << ejk << "," << '\n';
     for(int i=si; i<ei;i++){
        for(int j=sjk; j<ejk;j++){
            for(int k=sjk; k<ejk;k++){
                net = convolve(i,j,k,kern,ifn1);
                if(g.grid[i][j][k].alive==true)net=net*g.rc;
                if(g.grid[i][j][k].alive==false)net=net*g.rd;
                if(net<0) net=0.;
                ifn2[i][j][k] = net;
            }
        }
    }

}

// set up difusion kernel
void setupkern(float s, vector<vector<vector<float>>>& kern){

    float kijk; float count = 0;

    for (int i = 0; i < (3); ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {

                kijk= 1./pow(s*sqrt(2.*3.14),(float)3)*exp(-
                                                                   (pow((float)(i-1),(float)2)+
                                                                    pow((float)(j-1),(float)2)+
                                                                    pow((float)(k-1),(float)2))/
                                                           (2.*pow(s,(float)2)));
                count += kijk;
                kern[i][j][k] = kijk;
            }
        }
    }

    for (int i = 0; i < (3); ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                kern[i][j][k] = kern[i][j][k]/count;
            }
        }
    }
}

void erase(grid& g, tumour& t,  float& tt){

    list<Coordinate>::iterator dummy,it1;
    it1=t.deadcells.begin();
    int nerr; //ninc
    nerr = round(t.dead_count*t.dt*t.deg);
    for(int i=0;i<nerr;i++){
       dummy = it1; dummy--;
       g.remov(*it1);
       t.deadcells.erase(it1);
       it1=dummy;
       t.dead_count--;
       it1++;
    }

}




int main(int argc, char **argv) {

    try {
    // process argument directory  and load parameters
    string dir = argv[1];
    string pardir = dir + "/";
    string outdir = dir + "/output/";
    parms par=parms(dir);

    par.gsize;

    // derived parameters
    par.rdeg=1-par.rdeg; // ifn degredation
    par.rc=1-par.rc; // ifn conspumtipn
    float sigma = 2 * sqrt(par.D*dt);
    // initialise saving stats stream
    string statfile = outdir + "stats_" + to_string(time(NULL))+".csv";
    ofstream s_stats(statfile.c_str(),std::ofstream::out | std::ofstream::trunc);
    s_stats << "time" << "," << "N"<< "," <<"dead"<< "," <<  "CTL" <<  "," <<  "con" <<
                ","  <<  "tumrad" <<  "," <<  "ctlrad" <<   '\n';



    // Yet another 'temporary' solution:
    // will need to throw time, x, y, z, (maybe) what; onto these streams.
    // We will worry about capturing locations for dead cells separately in the visualisations.
    string growthfile = outdir + "growth_" + to_string(time(NULL))+".csv"; // time, x,y,z, conjugate (1/0).
    ofstream s_growth(growthfile.c_str(),std::ofstream::out | std::ofstream::trunc);

    string ctlfile = outdir + "ctl_" + to_string(time(NULL))+".csv"; // time, x,y,z, conjugate (1/0).
    ofstream s_ctl(ctlfile.c_str(),std::ofstream::out | std::ofstream::trunc);

    string killfile = outdir + "kill_" + to_string(time(NULL))+".csv"; // time, x,y,z.
    ofstream s_kill(killfile.c_str(),std::ofstream::out | std::ofstream::trunc);



        // initialise spatial grid and tumour
    tumour t; neighbours n;
    grid g; //= grid(501);
    g.initialise(par.gsize);
    site s;
    vector<vector<vector<site>>> temp(par.gsize,vector<vector<site>>(par.gsize,vector<site>(par.gsize,s)));
    g.grid=temp;
    g.buffer=buffer;
    t.maxhp=par.hits; t.deg=par.deg; t.dt=dt;
    startcancer(par,g,t,ref(s_growth));
    n.initalise();

    //ctls
    population p; p.nmig = 0,p.ncon = 0,p.nctl = 0; p.buffer=buffer;
    p.f=par.f; p.pdet=par.pdet;
    if (par.ctl) p.pcon=par.pcon; p.pdet=par.pdet;
    g.tifn=par.tifn; g.rc=par.rc; g.rd=par.rdeg;

    vector<vector<vector<float>>> ifn1;
    vector<vector<vector<float>>> ifn2;



    if(par.ifn==true) {
        vector<vector<vector<float>>> tmp1(par.gsize, vector<vector<float>>(par.gsize, vector<float>(par.gsize, 0)));
        vector<vector<vector<float>>> tmp2(par.gsize, vector<vector<float>>(par.gsize, vector<float>(par.gsize, 0)));
        ifn1=tmp1; ifn2=tmp2;

    }
    vector<vector<vector<float>>> kern(3, vector<vector<float>>(3, vector<float>(3, 0)));
    setupkern(sigma,kern);



    auto start = chrono::system_clock::now();
    auto end = chrono::system_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);

    for (;;){


            if (tt > 5 & par.ctl == true) {
                erase(g, t, tt); // erase any dead cells
                if (par.ifn == true) {
                    for (int i = 0; i < par.dsteps; i++) {
                        vector<vector<vector<float>>> tmp;
                        p.diffuse(ifn1, g);

                        std::vector<boost::thread *> t;
                        nThreads=configThreads(ref(p),ref(g),ref(par));
                        for (int  j = 0; j < nThreads; j++)

                            t.push_back(new boost::thread(diffuse,ref(g), ref(p), kern, ref(ifn1), ref(ifn2),
                                    iMinThread[j],iMaxThread[j]));

                        for (int j = 0; j < nThreads; j++)
                            t[j]->join();

                        tmp = ifn1;
                        ifn1 = ifn2;
                        ifn2 = tmp;

                    }
                }
                p.kills(g, t, par.k * dt, tt,ref(s_ctl),ref(s_kill));
                p.update(ctld, t, tt, par);
                p.migrate(g, t, tt,ref(s_ctl));


            };

            growth(par, dt, g, t, n, ifn1,ref(s_growth));
            //dispersal(par,dt,g,t,ifn1);

            s_stats << tt << "," << t.cell_count << "," << t.dead_count << "," << p.nctl << "," <<
                    p.nmig << "," << t.radius << "," << p.rad << '\n'; // pipe to statistics fle

            if (par.video == true &(( tt>6 & tt < 6.001)|( tt>7 & tt < 7.001)
            |( tt>8 & tt < 8.001)) |( tt>10 & tt < 10.001) ) {

                rAll(ifn1, p, t, frame, g, par, pardir);

            }
            if (frame % (60 * 2) == 0) {
                cout << "ABM SUMMARY:" << '\n';
                cout << "time:                  " << tt << '\n';
                cout << "tumour cell count:     " << t.cell_count << '\n';
                cout << "dead cells:            " << t.dead_count << '\n';
                cout << "Conjugate CTL:         " << p.ncon << '\n';
                cout << "Migratory CTL:         " << p.nmig << '\n';
                cout << "tumour radius:         " << t.radius << '\n';
                cout << "Max conjugate radius:  " << p.rad << '\n';
                cout << "Max migrated radius:   " << p.mrad << '\n';
                cout << "With number of threads:   " << nThreads << '\n';
                end = chrono::system_clock::now();
                elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
                cout << "time for last period in ms: " << elapsed.count() << '\n';
                start = chrono::system_clock::now();;
            }

            frame++;
            tt += dt;

            if (tt > end_t | t.cell_count < 1) {
                s_stats.close();

                return 0;
            }
        }

    }catch (const std::bad_alloc& e) {
        cout << "Allocation failed! Not enough memory for this!! ....  " << e.what() << endl;
    } catch(...){
        cout << "unknown failure! bad design! Bjarn is angry!" << endl;
    }

    return 0;
}
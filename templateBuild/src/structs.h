#include<list>
#include<string>
#include<math.h>
#include<vector>
#include <algorithm>
#include <array>
#include <random>
#include <ctime>
#include "functions.h"

using namespace std;

// setup

std::mt19937 gen(0); //gen(rd());
std::uniform_real_distribution<double> dis(0.0,1.0);
std::uniform_real_distribution<double> dis2(-1.0,1.0);

// this is flat out wrong Richard. You should use here a gaussian, I think. <BUG>
std::uniform_int_distribution<int> dsite(-1,1);
struct parms{
    string dir;
    vector<string> parnames {"ri", "ni", "rg","rd","dd","k","pcon","pdet","hits","deg","ctl",
                             "ifn","dsteps","f","tifn","D","rdeg","rc","mix","edge","gsize",
                             "video","Nthreads", "force"};
    vector<float> parvals{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    // parameters for growth only case
    int initrad=10; // initial radius of tumour
    int ni=100; // number of cells at t=0 to disperse within rini from grid centre
    float rg=1.8; // tumour cell growth rate
    float rd=0.08; // tumour cell dispersal probability
    float dist=0.01;// range of dispersal

    // parameters for CTL infiltration
    float k=4; // killing rate
    float pcon=1;
    float pdet=0;
    int hits=1;
    float deg=1;
    bool ctl=false;
    int Nthreads = 1;

    // ifn parameters
    bool ifn = false;
    int dsteps=3;
    float f = 1200./(float)dsteps;
    float tifn=1;
    float D = 750;  // controlling range of ifn diffusion
    float rdeg = 0.3;
    float rc = 0.1;
    bool mix = false;
    bool edge = false;
    bool force = true;

    int gsize=301;

    bool video = 0;

    parms(string dir){
        this->dir=dir;
        initialise();
    }
    parms();
    int match(string);
    void initialise();


};

int parms::match(string name){

    for(int i = 0;i<parnames.size();i++){
        if(parnames[i]==name) return(i) ;
    }
    return(-1);

};

void parms::initialise(){

    // read parameters
    string parms = dir + "/parms.csv";

    ifstream file;
    file.open(parms);
    if ( ! file ) {
        cout << "Error: Can't open the parameter file.\n";
        cout << "tried to open: \n" << endl;
        cout << parms << endl;
        exit(1);
    }
    string name;
    for (int i = 0; i < parvals.size(); i++) {
        file >> name;
        int index = match(name);
        file >> parvals[index];
    }

    initrad = round(parvals[0]); // initial radius of tumour
    ni = round (parvals[1]); // number of cells at t=0 to disperse within rini from grid centre
    rg = parvals[2]; // tumour cell growth rate
    rd = parvals[3]; // tumour cell dispersal probability
    dist = parvals[4]; // range of dispersal
    k =parvals[5];
    pcon =parvals[6];
    pdet =parvals[7];
    hits = parvals[8];
    deg =parvals[9];
    ctl =(bool)(round(parvals[10]));
    ifn = (bool)(round(parvals[11]));
    dsteps=round(parvals[12]);
    f = parvals[13]/(float)dsteps;
    tifn=parvals[14];
    D = parvals[15];  // controlling range of ifn diffusion
    rdeg = parvals[16];
    rc = parvals[17];
    mix = (bool)(round(parvals[18]));
    edge= (bool)(round(parvals[19]));
    gsize=round(parvals[20]);
    video = (bool)(round(parvals[21]));
    Nthreads = parvals[22];
    force = (bool)parvals[23];
    file.close();
};

struct Coordinate {
    int x, y, z;
    float displacement();
    void randinit(float r){
        do{x=round(r*dis2(gen));
        y=round(r*dis2(gen));
        z=round(r*dis2(gen));}while(displacement()>r);
    };
    void step(){
        x+=dsite(gen);
        y+=dsite(gen);
        z+=dsite(gen);
    };

    void randinit(float,float);
};

float Coordinate::displacement(){
    return(sqrt(pow((double)x,2.)+pow((double)y,2.)+pow((double)z,2.)));
};

void Coordinate::randinit(float r1, float r2){
do{
x=round(r2*dis2(gen));
y=round(r2*dis2(gen));
z=round(r2*dis2(gen));}
while(displacement()<r1 | displacement()>r2);
};

struct neighbours{
    vector<Coordinate> n;
    void initalise();
    void shuffle();
};

void neighbours::initalise() {
    for(int i=-1; i<2; i++){
        for(int j=-1; j<2; j++) {
            for (int k = -1; k < 2; k++) {
                Coordinate c;
                c.x=i;
                c.y=j;
                c.z=k;
                if(i!=0 | j!=0 | k!=0 )n.push_back(c);
            }
        }
    }
}

void neighbours::shuffle() {
    random_shuffle(n.begin(), n.end());
}

struct tumour_cell{

    Coordinate location; // could point to a grid site
    short int health; // could be 1 byte
    char type; // could be 1 bit
    short unsigned int nNeighbours;

};

struct tumour{
    list<Coordinate> deadcells;
    list<Coordinate>::iterator it = deadcells.begin();
    list<tumour_cell> cells;
    float radius=1;
    int cell_count=0;
    int dead_count = 0;
    int itnum = 50;
    float deg;
    float dt;
;
    short int maxhp;
    list<tumour_cell>::iterator add_cell(Coordinate,char);
    //list<tumour_cell*> initialise(int, int);
    void kill(list<tumour_cell>::iterator);
    void check_radius();


};

void tumour::check_radius(){
    radius = 0;
    for (list<tumour_cell>::iterator i = cells.begin(); i!=cells.end();i++){
        if(i->location.displacement()>radius) radius = i->location.displacement();
    }
};

void tumour::kill(list<tumour_cell>::iterator p) {
    Coordinate d=p->location;
    cells.erase(p);
//    deadcells.push_back(d);
    cell_count--;

if((radius-d.displacement())<0.1) check_radius();

}

list<tumour_cell>::iterator tumour::add_cell(Coordinate coord, char type) {
    if(coord.displacement()>radius) this->radius = coord.displacement();
    tumour_cell c; // create a tumour cell
    c.location=coord; //
    c.health=maxhp;
    c.type=type;
    list<tumour_cell>::iterator p; // create a pointer to it
    cells.push_back(c);

    // is this not odd?
    p= cells.end(); p--;
    cell_count++;
    return(p); // give the memory location of tumour cell to the grid

}

struct site{
    bool occupied=false;
    bool alive=false;
    list<tumour_cell>::iterator p; // point to a tumour cell
};

struct grid{
    vector<vector<vector<site>>> grid;
    int size,centre; // of the grid
    float tifn,rd,rc; // threshold for ifn effect,ifn degredation & consumption rates
    float maxr; // max radius of tumour cell
    int buffer = 5; // how far we keep cells away from grid edge.
    float max_permitted_radius;

    void initialise(int);
    void add(list<tumour_cell>::iterator);
    bool isoccupied(Coordinate&);
    bool isalive(Coordinate&);
    Coordinate map(Coordinate&);
    char type(Coordinate&);
    list<tumour_cell>::iterator kill(Coordinate&);
    list<tumour_cell>::iterator point(Coordinate&);
    void remov(Coordinate&);

};



list<tumour_cell>::iterator grid::point(Coordinate& c) {
    Coordinate index = map(c);
    return(grid[index.x][index.y][index.z].p);
}

list<tumour_cell>::iterator grid::kill(Coordinate& c){
    Coordinate index = map(c);
    grid[index.x][index.y][index.z].alive=false;
    return(grid[index.x][index.y][index.z].p);
}

void grid::remov(Coordinate& c){
    Coordinate index = map(c);
    grid[index.x][index.y][index.z].occupied=false;

};

char grid::type(Coordinate& c){

    Coordinate i=map(c);
    return(grid[i.x][i.y][i.z].p->type);

};

Coordinate grid::map(Coordinate& c) {

    Coordinate index;
    index.x = c.x+centre, index.y = c.y+centre, index.z=c.z+centre;

    return(index);

}

void grid::initialise(int gsize) {
    size = gsize;
    centre = (int)round((gsize-1)/2);
    maxr = (int)round((gsize-1)/2)-5;
    max_permitted_radius = ((float)gsize/2.)-(float)buffer;


}

void grid::add(list<tumour_cell>::iterator p) {


    Coordinate index = map(p->location);
    if(grid[index.x][index.y][index.z].occupied == true) cout << "err added 2nd cell to gridsite" <<endl;
    grid[index.x][index.y][index.z].occupied= true;
    grid[index.x][index.y][index.z].alive= true;
    grid[index.x][index.y][index.z].p= p;

}

bool grid::isoccupied(Coordinate& c) {
    Coordinate index = map(c);
    if(grid[index.x][index.y][index.z].occupied == true) return(true);
    else return(false);
}

bool grid::isalive(Coordinate& c) {
    Coordinate index = map(c);
    if(grid[index.x][index.y][index.z].alive == true) return(true);
    else return(false);
}

struct CTL{
    Coordinate location;
};

struct population{

    int buffer; // how far we simulate outside the furthest tumour cell
    float f; // ifn production molecules/min
    float pcon=1; // probability of conjugate formation
    float pdet=0; // ... .... .. of conjugate detachment
    float rad=0; // radius of outermost conjugated CTL
    float mrad=0;
    float p_replace=0.1;
    int densityMethod = 0;
    int nctl;
    list<CTL> mig;
    list<CTL> con;
    int nmig,ncon;
    void migrate(grid&, tumour&, float&, ofstream&);//migration function
    void kills(grid&, tumour&, float, float&, ofstream&, ofstream&);
    void update(vector<int>&, tumour&, float, parms&);
    void check_radius();
    void diffuse(vector<vector<vector<float>>>&,grid&);

};

void population::check_radius(){
    rad = 0;
    for (list<CTL>::iterator i = con.begin(); i!=con.end();i++){
        if(i->location.displacement()>rad) rad = i->location.displacement();
    }
};

void population::update(vector<int> &ctld, tumour &t, float tt, parms &par){

int ntar=(t.cell_count+t.dead_count)*getETRatio(tt,densityMethod);
nmig=mig.size();
ncon=con.size();
nctl = nmig + ncon;

int add=ntar-(nctl);

if(add>0){
    for(int i=0; i<add;i++){
        CTL c;

        if(par.edge==false)  c.location.randinit(t.radius+(float)buffer);
        if(par.edge==true) c.location.randinit(t.radius*0.5,t.radius);
        mig.push_back(c);
    }
}


}

void population::migrate(grid& g, tumour &t, float &tt, ofstream & s_ctl){ // assuming conjugates always formed;
    mrad=0;bool erased;
    list<CTL>::iterator dummy;

    for(list<CTL>::iterator i = mig.begin(); i != mig.end(); i++ ){
        //s_ctl << tt << ","
        //<< i->location.x << "," << i->location.y << "," << i->location.z << ","
        //<< 0 << "," << '\n';

        erased=false;
        i->location.step(); // step location

        // check if we are furthest CTL in sim
        if(i->location.displacement()>mrad) mrad= i->location.displacement();

        // test if still in tumour:
        if(i->location.displacement()>(t.radius+buffer)) {
            dummy = i; dummy--;
            mig.erase(i);
            i=dummy;
            erased=true;
        }

        // also remove with probability
        if(dis(gen)<p_replace &!erased) {
            dummy = i; dummy--;
            mig.erase(i);
            i=dummy;
            erased=true;
        }

        // attempt to form conjugate:
        if(!erased) {
            if (g.isalive(i->location)){
                if (g.type(i->location) == '7') {
                    if ((dis(gen) < pcon)) {
                        dummy = i;
                        dummy--;
                        con.push_back(*i);
                        if (i->location.displacement() > rad) rad = i->location.displacement();
                        mig.erase(i);
                        i = dummy;
                    }
                }
            }
        }
    }

};

void population::kills(grid &g, tumour &t, float p,  float &tt, ofstream &s_ctl, ofstream &s_kill){
    // here one of 4 things happens
    // we find no target and we become migrating
    // we find a target and we kill it
    // we detach
    // nothing
    // they are mutually exclusive

    bool check,alive;
    list<CTL>::iterator dummy,i;
    i= con.begin();
    int cm =con.size();
    check=false;
    for(int c = 0; c<cm;c++){


       // s_ctl << tt << ","
       //       << i->location.x << "," << i->location.y << "," << i->location.z << ","
       //       << 1 << "," << '\n';

        dummy=i;

        // we find no living target and we become migrating
        alive = g.isalive(i->location);
        if(alive!=true){
            if((rad-i->location.displacement())<0.1) check =true;
            mig.push_back(*i);
            dummy++;
            con.erase(i);
        }

        // we find a target and we kill it
        // only place in code for tumour cells to die
        bool kill = dis(gen)<p;
        if(kill==true & alive){

            //  point to a tumour cell (residing in our tumour cell list)
            // how do we know its still there??
            list<tumour_cell>::iterator tar = g.point(i->location);
            tar->health--;

            if(tar->health<1) {// kill tumour cell
                if((rad-i->location.displacement())<0.1) check =true;
             //   s_kill << tt << "," << i->location.x << "," << i->location.y << "," << i->location.z << '\n';
                list<tumour_cell>::iterator it = g.kill(i->location);
                t.kill(it);
                t.deadcells.push_back(i->location); t.dead_count++;
                mig.push_back(*i);
                dummy++;
                con.erase(i);
                // what does i now point to
            }
        };

        // we detach
        if((kill!=true) & dis(gen)<pdet){
            if((rad-i->location.displacement())<0.1) check =true;
            mig.push_back(*i);
            dummy++;
            con.erase(i);
        };
        if(i==dummy) dummy++;
        i = dummy;
    };

    if (check==true){
        check_radius();
    }
};

void population::diffuse(vector<vector<vector<float>>>& ifn,grid &g){
    for(list<CTL>::iterator p = con.begin(); p != con.end(); p++){
        Coordinate i = g.map(p->location);
        ifn[i.x][i.y][i.z]+=f;
    }
}


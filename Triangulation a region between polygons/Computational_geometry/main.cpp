#include <math.h>
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <queue>
#include <set>
#include <string>
#include "splay.h"
#include <sys/time.h>
#include <time.h>
#include <stack>
#include <vector>
using namespace std;
long   int    l_id=0,p_id=0;
using namespace std;

#define sqr(t)  (t)*(t)

const double    PI=3.141592653589793238462643383279502884197169399375105820974944592308;
enum  Type      { UNKNOWN, INPUT, INSERT, START, END, MERGE, SPLIT, REGULAR_UP, REGULAR_DOWN};

template <class T, class KeyType>      class    SplayTree;
//base class for points;
class Pointbase
{
public:
    //constructors and destructor
    Pointbase() {}
    Pointbase(const Pointbase& pb);
    
    Pointbase(double xx, double yy)
    :id(0), x(xx), y(yy), type(UNKNOWN) { }
    
    Pointbase(int idd, double xx, double yy)
    :id(idd), x(xx), y(yy), type(UNKNOWN) { }
    
    Pointbase(double xx, double yy, Type ttype)
    :id(0), x(xx), y(yy), type(ttype) { }
    
    Pointbase(int idd, double xx, double yy, Type ttype)
    :id(idd),x(xx), y(yy), type(ttype) { }
    
    //rotate a point by angle theta, not used;
    void rotate(double theta)
    {
        double cosa=cos(theta),sina=sin(theta),newx,newy;
        
        newx=x*cosa-y*sina;
        newy=x*sina+y*cosa;
        x=newx;
        y=newy;
    }
    
    //operator overloading
    friend  bool operator==(const Pointbase&, const Pointbase&);
    friend  bool operator>(const Pointbase&, const Pointbase&);
    friend  bool operator<(const Pointbase&, const Pointbase&);
    friend  bool operator!=(const Pointbase&, const Pointbase&);
    
    //public data
    unsigned int    id;              //id of point;
    double          x, y;            //coordinates;
    Type            type;            //type of points;
    bool            left;            //left chain or not;
    int my_id;
    
};


//base class for polygon boundary
//Linebase class is a directed line segment with start/end point
class Linebase
{
public:
    //constructors and destructor
    Linebase();
    Linebase(Pointbase* ep1, Pointbase* ep2, Type type);
    Linebase(const Linebase& line);
    ~Linebase() {};
    
    unsigned int id() const { return _id; }
    
    //two end points
    Pointbase*   endPoint(int i) const { return _endp[i]; }
    Type         type() const { return _type; }
    double       keyValue() const { return _key; }
    void         setKeyValue(double y);
    //slightly increased the key to avoid duplicated key for searching tree.
    void         increaseKeyValue(const double diff) { _key+=diff; }
    //reverse a directed line segment; reversable only for inserted diagonals
    void         reverse();
    
    //set and return helper of a directed line segment;
    void         setHelper(unsigned int i) { _helper=i; }
    unsigned int helper() { return _helper; }
    
    //operator overloading
    
protected:
    unsigned int _id;           //id of a line segment;
    Pointbase*   _endp[2];      //two end points;
    
    Type         _type;         //type of a line segement, input/insert
    double       _key;          //key of a line segment for splay tree searching
    unsigned int _helper;       //helper of a line segemnt
};



class Polygon
{
public:
    
    map<unsigned int, Pointbase*>& points()    { return _points; }
    map<unsigned int, Linebase*> &      edges()     { return _edges; }
    
    unsigned int            _ncontours;   //number of contours
    vector<unsigned int>    _nVertices;   //
    map<unsigned int, Pointbase*>            _points;      //all vertices
    map<unsigned int, Linebase*>                  _edges;       //all edges
    double                  _xmin,_xmax, _ymin,_ymax; //boundary box for polygon
    string                  _prefix;     //prefix of associated polygon bdm file;
    //constructor and destructor
    Polygon(vector<vector<pair<double, double > > > a);
    ~Polygon();
    
    // main member function for polygon triangulation;
    void         partition2Monotone();
    void         searchMonotones();
    void         triangulation();
    void         delTriangulation() { cout<<"not available for this version\n";}
    
    //return all triangles
    list<vector<unsigned int> >    triangles() { return _triangles; }
    
    //output file format;
    void         setDebugOption(bool debug);
    void         saveAsShowme();
    void         saveAsTecplot();
    void         saveAsMetaPost();
    
    //private member functions.
private:
    //rotate input polygon by angle theta, not used;
    void         rotate(double theta);
    void         initializate();
    
    //prev or next point/edge id for a given ith point/edge;
    unsigned int prev(unsigned int i);
    unsigned int next(unsigned int i);
    
    //handle event vertext according to vertex type;
    void         handleStartVertex(unsigned int);
    void         handleEndVertex(unsigned int);
    void         handleSplitVertex(unsigned int);
    void         handleMergeVertex(unsigned int);
    void         handleRegularVertexUp(unsigned int);
    void         handleRegularVertexDown(unsigned int);
    
    //add diagonal between two vertices;
    void         addDiagonal(unsigned int i, unsigned int j);
    
    
    //angle ABC for three given points, for monotone polygon searching purpose;
    double       angleCosb(double *A, double *B, double *C);
    //find the next edge, for monotone polygon searching purpose;
    unsigned int selectNextEdge(Linebase* edge);
    
    //triangulate a monotone polygon piece;
    void         triangulateMonotone(list<unsigned int> & mpoly);
    
    //private data memebers
    priority_queue<Pointbase>       _qpoints;                            //priority queue for event points
    SplayTree<Linebase*, double>       _edgebst;                            //edge binary searching tree (splaytree)
    list<list<unsigned int> >   _mpolys;                             //all monotone polygon piece list;
    list<vector<unsigned int> >   _triangles;                          //all triangle list;
    
    //data for monotone piece searching purpose;
    map<unsigned int, set<unsigned int> >   _startAdjEdgeMap;                    //all edges starting from given points (map)
    map<unsigned int, Linebase*>      _diagonals;                          //added diagonals to partition polygon to
    //monotont pieces, not all diagonals of
    //given polygon
    
    bool        _debug;                              //debug option;                            //log file for debug purpose;
    
};

//Jonathan schewchuk's exact arithmetic code, see predicates.cc for detais;
extern double orient2d(double* pa, double* pb, double* pc){
    return pa[0]*(pb[1]-pc[1])+pb[0]*(pc[1]-pa[1])+pc[0]*(pa[1]-pb[1]);
}

//----------------------------------------------------------------------------
//square of the distance of two points;
//----------------------------------------------------------------------------
double dist_sqr(const Pointbase& sp, const Pointbase& ep)
{
    return sqr(sp.x-ep.x)+sqr(sp.y-ep.y);
}

//----------------------------------------------------------------------------
//square of the distance of two points;
//----------------------------------------------------------------------------
double dist_sqr(double *pa, double *pb)
{
    return sqr(pa[0]-pb[0])+sqr(pa[1]-pb[1]);
}

void UpdateKey(BTreeNode<Linebase*,double>* node, double y)
{
    node->data()->setKeyValue(y);
}

void CheckTime(struct timeval tv0, struct timeval tv1)
{
    double runtime=1000l*(tv1.tv_sec - tv0.tv_sec)+(tv1.tv_usec - tv0.tv_usec) / 1000l;
//    cout<<"Triangulation Time::"<<runtime/1000.0<<" ( s )\n";
}

//----------------------------------------------------------------------------
//copy constructor
//----------------------------------------------------------------------------
Pointbase::Pointbase(const Pointbase& pb)
{
    this->id=pb.id;
    this->x=pb.x;
    this->y=pb.y;
    this->type=pb.type;
    this->left=pb.left;
}

//----------------------------------------------------------------------------
//operator ( ==, >, < and != ) overloading for pointbase class
//----------------------------------------------------------------------------
bool operator==(const Pointbase& pa, const Pointbase& pb)
{
    return (pa.x==pb.x && pa.y==pb.y);
}

//----------------------------------------------------------------------------
bool operator>(const Pointbase& pa, const Pointbase& pb)
{
    return( (pa.y > pb.y) || ( (pa.y==pb.y) && (pa.x < pb.x)) );
}

//----------------------------------------------------------------------------
bool operator<(const Pointbase& pa, const Pointbase& pb)
{
    return( (pa.y < pb.y) || ( (pa.y==pb.y) && (pa.x > pb.x)) );
}

//----------------------------------------------------------------------------
bool operator!=(const Pointbase& pa, const Pointbase& pb)
{
    return !(pa.x==pb.x && pa.y==pb.y);
}



//----------------------------------------------------------------------------
//Linebase construct
//----------------------------------------------------------------------------
Linebase::Linebase():_type(UNKNOWN)
{
    for(int i=0; i<2; i++) _endp[i]=0;
    _id=0;
}

//-----------------------------------------------------------------------------
//Linebase construct
//-----------------------------------------------------------------------------
Linebase::Linebase(Pointbase* sp, Pointbase* ep, Type type):_type(type)
{
    _endp[0]=sp;
    _endp[1]=ep;
    //_key=_endp[0]->x < _endp[1]->x ? _endp[0]->x:_endp[1]->x;
    _id=++l_id;
}

//----------------------------------------------------------------------------
//copy constructor
//----------------------------------------------------------------------------
Linebase::Linebase(const Linebase& line)
{
    this->_id=line._id;
    this->_endp[0]=line._endp[0];
    this->_endp[1]=line._endp[1];
    this->_key=line._key;
    this->_helper=line._helper;
}


//----------------------------------------------------------------------------
//reverse a directed line segment, reverseable only for insert diagonals
//----------------------------------------------------------------------------
void Linebase::reverse()
{
    assert(_type==INSERT);
    Pointbase* tmp=_endp[0];
    _endp[0]=_endp[1];
    _endp[1]=tmp;
}

void Linebase::setKeyValue(double y)
{
    if( _endp[1]->y==_endp[0]->y )
        _key=_endp[0]->x < _endp[1]->x ? _endp[0]->x:_endp[1]->x;
    else    _key=( y - _endp[0]->y ) * ( _endp[1]->x - _endp[0]->x ) / (_endp[1]->y - _endp[0]->y ) + _endp[0]->x;
}




//----------------------------------------------------------------------------
//polygon class constructor
//----------------------------------------------------------------------------
struct pt {
    double x, y;
};

bool cmp (pt a, pt b) {
    return a.x < b.x || a.x == b.x && a.y < b.y;
}

bool cw (pt a, pt b, pt c) {
    return a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y) < 0;
}

bool ccw (pt a, pt b, pt c) {
    return a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y) > 0;
}

void convex_hull (vector<pt> & a) {
    if (a.size() == 1)  return;
    sort (a.begin(), a.end(), &cmp);
    pt p1 = a[0],  p2 = a.back();
    vector<pt> up, down;
    up.push_back (p1);
    down.push_back (p1);
    for (size_t i=1; i<a.size(); ++i) {
        if (i==a.size()-1 || !ccw (p1, a[i], p2)) {
            while (up.size()>=2 && ccw (up[up.size()-2], up[up.size()-1], a[i]))
                up.pop_back();
            up.push_back (a[i]);
        }
        if (i==a.size()-1 || !cw (p1, a[i], p2)) {
            while (down.size()>=2 && cw (down[down.size()-2], down[down.size()-1], a[i]))
                down.pop_back();
            down.push_back (a[i]);
        }
    }
    a.clear();
    for (size_t i=0; i<up.size(); ++i)
        a.push_back (up[i]);
    for (size_t i=down.size()-2; i>0; --i)
        a.push_back (down[i]);
}

pair<pair<double, double>, pair<double, double> > tr(int i, vector<pair<double, double> > & aa){
    int n = (int)aa.size();
    return {{aa[i].first, aa[i].second}, {aa[(i + 1) % n].first, aa[(i + 1) % n].second}};
}
bool eq(pair<pair<double, double>, pair<double, double> > a, pair<pair<double, double>, pair<double, double> > b){
    double eps = 1e-20;
    if(fabs(a.first.first - b.first.first) > eps) return false;
    if(fabs(a.first.second - b.first.second) > eps) return false;
    if(fabs(a.second.first - b.second.first) > eps) return false;
    if(fabs(a.second.second - b.second.second) > eps) return false;
    return true;
}
vector<vector<pair<double, double>  > > main_filename;
Polygon::Polygon(vector<vector<pair<double, double>  > >  filename)
{
    
    vector<pt > aaa;
    for(auto it : filename){
        for(auto itt : it){
            aaa.push_back({itt.first, itt.second});
        }
    }
    convex_hull(aaa);
    reverse(aaa.begin(), aaa.end());
    vector<pair<double, double > > aa;
    for(auto it : aaa){
        aa.push_back({it.x, it.y});
    }
    // we need homotety to make convex hull a little bit bigger
    double xx = 0;
    double yy = 0;
    int nn = (int)aa.size();
    for(auto it : aa){
        xx += it.first;
        yy += it.second;
    }
    xx /= nn;
    yy /= nn;
    for(auto &it: aa){
        it.first = xx + (it.first - xx) * 1.0001;
        it.second = yy + (it.second - yy) * 1.0001;
    }
    unsigned int i = 1, first, last, num;
    double x,y;
    Type type;
    _ncontours=0;
    vector<pair<double, double>  > current_vector;
    reverse(filename.begin(), filename.end());
    filename.push_back(aa);
    reverse(filename.begin(), filename.end());
    for(int qq = 0; qq < (int)filename.size(); ++qq)
    {
        current_vector.clear();
        if(qq >= 1){
            current_vector = filename[qq];
            reverse(current_vector.begin(), current_vector.end());
            
        }
        else{
            current_vector = filename[qq];
        }
        num = (int)current_vector.size();
        
        _nVertices.push_back( num );
        
        first = i;
        last = first + _nVertices[_ncontours] - 1;
        for (unsigned int j = 0; j < _nVertices[_ncontours]; j++, i++)
        {
            x = current_vector[j].first;
            y = current_vector[j].second;
            aa.push_back({x, y});
            type=INPUT;
            
            Pointbase* point=new Pointbase(i,x,y,type);
            if(x > _xmax ) _xmax=x;
            if(x < _xmin ) _xmin=x;
            if(y > _ymax ) _ymax=y;
            if(y < _ymin ) _ymin=y;
            //point->rotate(PI/2.0);
            point -> my_id = 0;
            if(qq) point -> my_id = 0;
            _points[i]=point;
        }
        _ncontours++;
    }
    
    int sid,eid;
    num=0;
    
    for(unsigned int j=0; j<_ncontours; j++)
    {
        for(i=1; i<=_nVertices[j]; i++)
        {
            sid=num+i;
            eid=(i==_nVertices[j])?num+1:num+i+1;
            type = INPUT;
            Linebase* line=new Linebase(_points[sid], _points[eid], type);
            _edges[l_id]=line;
        }
        num+=_nVertices[j];
    }
    
    int sum=0;
    for(unsigned int i=0; i<_ncontours; i++)
    {
        sum+= _nVertices[i];
        _nVertices[i]=sum;
    }
    
    p_id=num;
    
    initializate();
    _debug=false;
    main_filename = filename;
}

//----------------------------------------------------------------------------
//polygon destructor
//----------------------------------------------------------------------------
Polygon::~Polygon()
{
    
    //clear all dynamic allocated memory
    map<unsigned int, Pointbase*>::iterator itp=_points.begin();
    for(; itp!=_points.end(); itp++)
    {
        delete itp->second;
    }
    
    map<unsigned int, Linebase*> ::iterator itl=_edges.begin();
    for(; itl!=_edges.end(); itl++)
    {
        delete itl->second;
    }
    
}

//----------------------------------------------------------------------------
//return the previous point (or edge) id for a given ith point (or edge);
//----------------------------------------------------------------------------
unsigned int Polygon::prev(unsigned int i)
{
    unsigned int j(0),prevLoop(0),currentLoop(0);
    
    while ( i > _nVertices[currentLoop] )
    {
        prevLoop=currentLoop;
        currentLoop++;
    }
    
    if( i==1 || (i==_nVertices[prevLoop]+1) ) j=_nVertices[currentLoop];
    else if( i <= _nVertices[currentLoop] ) j=i-1;
    
    return j;
}

//----------------------------------------------------------------------------
//return the next point (or edge) id for a given ith point (or edge);
//----------------------------------------------------------------------------
unsigned int Polygon::next(unsigned int i)
{
    unsigned int j(0),prevLoop(0),currentLoop(0);
    
    while ( i > _nVertices[currentLoop] )
    {
        prevLoop=currentLoop;
        currentLoop++;
    }
    
    if( i < _nVertices[currentLoop] ) j=i+1;
    else if ( i==_nVertices[currentLoop] )
    {
        if( currentLoop==0) j=1;
        else j=_nVertices[prevLoop]+1;
    }
    
    return j;
}


//----------------------------------------------------------------------------
//rotate a polygon by angle theta, reference point (0,0), not used;
//----------------------------------------------------------------------------
void Polygon::rotate(double theta)
{
    map<unsigned int, Pointbase*>::iterator it=_points.begin();
    for(; it!=_points.end(); it++)
        it->second->rotate(theta);
}


//----------------------------------------------------------------------------
//polygon initialization;
//to find types of all polygon vertices;
//create a priority queue for all vertices;
//construct an edge set for each vertex (the set holds all edges starting from
//the vertex, only for loop searching purpose).
//----------------------------------------------------------------------------
void Polygon::initializate()
{
    map<unsigned int, Pointbase*>::iterator it=_points.begin();
    for(; it!=_points.end(); it++)
    {
        int id=it->first;
        int idp=prev(id);
        int idn=next(id);
        Pointbase p=*_points[id], pnext=*_points[idn], pprev=*_points[idp];
        
        if( p > pnext && pprev > p )
            _points[id]->type=REGULAR_DOWN;
        else if (p > pprev && pnext > p)
            _points[id]->type=REGULAR_UP;
        else
        {
            double pa[2], pb[2], pc[2];
            
            pa[0]=_points[idp]->x;
            pa[1]=_points[idp]->y;
            
            pb[0]=_points[id]->x;
            pb[1]=_points[id]->y;
            
            pc[0]=_points[idn]->x;
            pc[1]=_points[idn]->y;
            
            double area=orient2d(pa,pb,pc);
            
            if( pprev > p && pnext > p ) _points[id]->type=((area >0) ^ _points[id] -> my_id) ? END: MERGE ;
            if( pprev < p && pnext < p ) _points[id]->type=((area >0) ^ _points[id] -> my_id) ? START : SPLIT;
            
        }
        
        _qpoints.push(*(it->second));
        
        _startAdjEdgeMap[id].insert(id);
        
    }
}

//----------------------------------------------------------------------------
//Add a diagonal from point id i to j
//----------------------------------------------------------------------------
void Polygon::addDiagonal(unsigned int i, unsigned int j)
{
    Type type=INSERT;
    Linebase* diag=new Linebase(_points[i], _points[j], type);
    _edges[diag->id()]=diag;
    
    _startAdjEdgeMap[i].insert(diag->id());
    _startAdjEdgeMap[j].insert(diag->id());
    
    _diagonals[diag->id()]=diag;
    
    if(_debug) cout<<"Add Diagonal from "<<i<<" to "<<j<<'\n';
}

//----------------------------------------------------------------------------
//Handle start vertex
//----------------------------------------------------------------------------
void Polygon::handleStartVertex(unsigned int i)
{
    double y=_points[i]->y;
    _edgebst.InOrder(UpdateKey, y);
    
    _edges[i]->setHelper(i);
    _edges[i]->setKeyValue(y);
    _edgebst.Insert(_edges[i]);
    
    if(_debug)
    {
        cout<<"set e"<<i<<" helper to "<<i<<'\n';
        cout<<"Insert e"<<i<<" to splay tree\n";
        cout<<"key:"<<_edges[i]->keyValue()<<'\n';
    }
}

//----------------------------------------------------------------------------
//Handle end vertex
//----------------------------------------------------------------------------
void Polygon::handleEndVertex(unsigned int i)
{
    double y=_points[i]->y;
    _edgebst.InOrder(UpdateKey, y);
    
    unsigned int previ=prev(i);
    Linebase* edge=_edges[previ];
    unsigned int helper=_edges[previ]->helper();
    
    
    if(_points[helper]->type==MERGE) addDiagonal(i, helper);
    _edgebst.Delete(edge->keyValue());
    
    if(_debug)
    {
        cout<<"Remove e"<<previ<<" from splay tree\n";
        cout<<"key:"<<edge->keyValue()<<'\n';
    }
}

//----------------------------------------------------------------------------
//Handle split vertex
//----------------------------------------------------------------------------
void Polygon::handleSplitVertex(unsigned int i)
{
    double x=_points[i]->x, y=_points[i]->y;
    _edgebst.InOrder(UpdateKey, y);
    
    BTreeNode<Linebase*, double>*  leftnode;
    _edgebst.FindMaxSmallerThan(x, leftnode);
    Linebase* leftedge=leftnode->data();
    
    unsigned int helper=leftedge->helper();
    addDiagonal(i, helper);
    
    if(_debug)
    {
        cout<<"Search key:"<<x<<" edge key:"<<leftedge->keyValue()<<'\n';
        cout<<"e"<<leftedge->id()<<" is directly left to v"<<i<<'\n';
        cout<<"Set e"<<leftedge->id()<<" helper to "<<i<<'\n';
        cout<<"set e"<<i<<" helper to "<<i<<'\n';
        cout<<"Insert e"<<i<<" to splay tree\n";
        cout<<"Insert key:"<<_edges[i]->keyValue()<<'\n';
    }
    
    leftedge->setHelper(i);
    _edges[i]->setHelper(i);
    _edges[i]->setKeyValue(y);
    _edgebst.Insert(_edges[i]);
}


//----------------------------------------------------------------------------
//Handle merge vertex
//----------------------------------------------------------------------------
void Polygon::handleMergeVertex(unsigned int i)
{
    double x=_points[i]->x, y=_points[i]->y;
    _edgebst.InOrder(UpdateKey, y);
    
    unsigned int previ=prev(i);
    unsigned int helper=_edges[previ]->helper();
    if (_points[helper]->type==MERGE) addDiagonal(i, helper);
    _edgebst.Delete(_edges[previ]->keyValue());
    if(_debug)
    {
        cout<<"e"<<previ<<" helper is "<<helper<<'\n';
        cout<<"Remove e"<<previ<<" from splay tree.\n";
    }
    
    BTreeNode<Linebase*, double>*  leftnode;
    _edgebst.FindMaxSmallerThan(x, leftnode);
    Linebase* leftedge=leftnode->data();
    
    helper=leftedge->helper();
    if(_points[helper]->type==MERGE) addDiagonal(i, helper);
    
    leftedge->setHelper(i);
    
    if(_debug)
    {
        cout<<"Search key:"<<x<<" found:"<<leftedge->keyValue()<<'\n';
        cout<<"e"<<leftedge->id()<<" is directly left to v"<<i<<'\n';
        cout<<"Set e"<<leftedge->id()<<" helper to "<<i<<'\n';
    }
}

//----------------------------------------------------------------------------
//Handle regular down vertex
//----------------------------------------------------------------------------
void Polygon::handleRegularVertexDown(unsigned int i)
{
    double y=_points[i]->y;
    _edgebst.InOrder(UpdateKey, y);
    
    unsigned int previ=prev(i);
    unsigned int helper=_edges[previ]->helper();
    if(_points[helper]->type==MERGE) addDiagonal(i, helper);
    
    _edgebst.Delete(_edges[previ]->keyValue());
    _edges[i]->setHelper(i);
    _edges[i]->setKeyValue(y);
    _edgebst.Insert(_edges[i]);
    
    if(_debug)
    {
        cout<<"e"<<previ<<" helper is "<<helper<<'\n';
        cout<<"Remove e"<<previ<<" from splay tree.\n";
        cout<<"Set e"<<i<<" helper to "<<i<<'\n';
        cout<<"Insert e"<<i<<" to splay tree\n";
        cout<<"Insert key:"<<_edges[i]->keyValue()<<'\n';
    }
}


//----------------------------------------------------------------------------
////Handle regular up vertex
//----------------------------------------------------------------------------
void Polygon::handleRegularVertexUp(unsigned int i)
{
    double x=_points[i]->x, y=_points[i]->y;
    _edgebst.InOrder(UpdateKey, y);
    
    BTreeNode<Linebase*, double>*  leftnode;
    _edgebst.FindMaxSmallerThan(x, leftnode);
    
    Linebase* leftedge=leftnode->data();
    
    unsigned int helper=leftedge->helper();
    if(_points[helper]->type==MERGE) addDiagonal(i, helper);
    leftedge->setHelper(i);
    
    if(_debug)
    {
        cout<<"Search key:"<<x<<" found:"<<leftedge->keyValue()<<'\n';
        cout<<"e"<<leftedge->id()<<" is directly left to v"<<i<<" and its helper is:"<<helper<<'\n';
        cout<<"Set e"<<leftedge->id()<<" helper to "<<i<<'\n';
    }
}

//----------------------------------------------------------------------------
//partition polygon to monotone polygon pieces
//----------------------------------------------------------------------------
void Polygon::partition2Monotone()
{
    if(_qpoints.top().type!=START)
    {
        cout<<"Please check your input polygon:\n1)orientations?\n2)duplicated points?\n";
        cout<<"stopped.\n";
        exit(1);
    }
    
    while(!_qpoints.empty())
    {
        Pointbase vertex=_qpoints.top();
        _qpoints.pop();
        unsigned int id=vertex.id;
        
        if(_debug)
        {
            string stype;
            switch (vertex.type)
            {
                case START:        stype="START";       break;
                case END:          stype="END";         break;
                case MERGE:        stype="MERGE";       break;
                case SPLIT:        stype="SPLIT";       break;
                case REGULAR_UP:   stype="REGULAR_UP";  break;
                case REGULAR_DOWN: stype="REGULAR_DOWN";break;
                default:
                    cout<<"No duplicated points please! stopped\n";
                    exit(1); break;
            }
            
            cout<<"\n\nHandle vertex:"<<vertex.id<<" type:"<<stype<<'\n';
        }
        
        
        switch(vertex.type)
        {
            case START:        handleStartVertex(id);       break;
            case END:          handleEndVertex(id);         break;
            case MERGE:        handleMergeVertex(id);       break;
            case SPLIT:        handleSplitVertex(id);       break;
            case REGULAR_UP:   handleRegularVertexUp(id);   break;
            case REGULAR_DOWN: handleRegularVertexDown(id); break;
            default:
                cout<<"No duplicated points please! stopped\n";
                exit(1); break;
        }
    }
}


//----------------------------------------------------------------------------
//two Auxiliary functions to find monotone polygon pieces
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
//calculate angle B for A, B, C three given points
//----------------------------------------------------------------------------
double Polygon::angleCosb(double *pa, double *pb, double *pc)
{
    double dxab = pa[0] - pb[0];
    double dyab = pa[1] - pb[1];
    
    double dxcb = pc[0] - pb[0];
    double dycb = pc[1] - pb[1];
    
    double dxab2 = dxab * dxab;
    double dyab2 = dyab * dyab;
    double dxcb2 = dxcb * dxcb;
    double dycb2 = dycb * dycb;
    double ab = dxab2 + dyab2;
    double cb = dxcb2 + dycb2;
    
    double cosb = dxab * dxcb + dyab * dycb;
    double denom = sqrt( ab * cb);
    
    cosb/=denom;
    
    return cosb;
}

//----------------------------------------------------------------------------
//for any given edge, find the next edge we should choose when searching for
//monotone polygon pieces;
//----------------------------------------------------------------------------
unsigned int Polygon::selectNextEdge(Linebase* edge)
{
    
    unsigned int eid= edge->endPoint(1)->id;
    set<unsigned int> edges=_startAdjEdgeMap[eid];
    assert(!edges.empty());
    
    unsigned int nexte=0;
    if( edges.size() == 1 )  nexte=*(edges.begin());
    else if( edges.size() > 1 )
    {
        unsigned int nexte_ccw(0), nexte_cw(0);
        double max=-2.0,min=2.0;
        
        
        set<unsigned int>::iterator it=edges.begin();
        for(; it!=edges.end(); it++)
        {
            if(*it==edge->id()) continue;
            double A[2], B[2], C[2];
            A[0]=edge->endPoint(0)->x;        A[1]=edge->endPoint(0)->y;
            B[0]=edge->endPoint(1)->x;        B[1]=edge->endPoint(1)->y;
            
            if(edge->endPoint(1)!=_edges[*it]->endPoint(0)) _edges[*it]->reverse();
            C[0]=_edges[*it]->endPoint(1)->x; C[1]=_edges[*it]->endPoint(1)->y;
            
            double area=orient2d(A, B, C);
            double cosb=angleCosb(A, B, C);
            
            if( area > 0 && max < cosb ) { nexte_ccw=*it; max=cosb; }
            else if( min > cosb ) { nexte_cw=*it; min=cosb; }
        }
        
        nexte = (nexte_ccw!=0) ? nexte_ccw : nexte_cw;
    }
    
    return nexte;
}

//----------------------------------------------------------------------------
//searching all monotone pieces;
//----------------------------------------------------------------------------
void Polygon::searchMonotones()
{
    int loop=0;
    
    map<unsigned int, Linebase*>  edges=_edges;
    
    while( edges.size() > _diagonals.size() )
    {
        loop++;
        list<unsigned int>  poly;
        map<unsigned int, Linebase*> ::iterator it=edges.begin();
        Pointbase* startp=startp=it->second->endPoint(0);
        Pointbase* endp=0;
        Linebase*  next=it->second;
        
        poly.push_back(startp->id);
        
        if(_debug)
        {
            cout<<"Searching for loops:"<<loop<<'\n';
            cout<<"vertex index:"<<startp->id<<" ";
        }
        
        for(;;)
        {
            endp=next->endPoint(1);
            if(next->type()!=INSERT)
            {
                edges.erase(next->id());
                _startAdjEdgeMap[next->endPoint(0)->id].erase(next->id());
            }
            if(endp==startp) break;
            poly.push_back(endp->id);
            
            if(_debug) cout<<endp->id<<" ";
            
            unsigned int nexte=selectNextEdge(next);
            
            if(nexte==0)
            {
                cout<<"Please check your input polygon:\n";
                cout<<"1)orientations?\n2)with duplicated points?\n3)is a simple one?\n";
                cout<<"stopped.\n";
                exit(1);
            }
            //assert( nexte > 0);
            next=edges[nexte];
            if(next->endPoint(0) !=endp ) next->reverse();
        }
        
        if(_debug) cout<<"\nloop closed!\n\n";
        
        _mpolys.push_back(poly);
    }
}


//----------------------------------------------------------------------------
//triangulate a monotone polygon;
//----------------------------------------------------------------------------
void  Polygon::triangulateMonotone(list<unsigned int> & mpoly)
{
    
    priority_queue<Pointbase>   qvertex;
    list<unsigned int> ::iterator it=mpoly.begin(), itnext;
    for(; itnext=it, it!=mpoly.end(); it++)
    {
        itnext++;
        if(itnext==mpoly.end()) itnext=mpoly.begin();
        Pointbase point=*_points[*it], pointnext=*_points[*itnext];
        point.left=(point > pointnext)? true:false;
        qvertex.push(point);
    }
    
    stack<Pointbase> spoint;
    for(int i=0; i<2; i++) { spoint.push(qvertex.top()); qvertex.pop(); }
    
    while ( qvertex.size() > 1 )
    {
        Pointbase topQueuePoint=qvertex.top();
        Pointbase topStackPoint=spoint.top();
        
        if(topQueuePoint.left!=topStackPoint.left)
        {
            while ( spoint.size()  > 1 )
            {
                Pointbase p1=spoint.top();
                spoint.pop();
                Pointbase p2=spoint.top();
                vector<unsigned int> v(3);
                v[0]=topQueuePoint.id;
                v[1]=p1.id;
                v[2]=p2.id;
                _triangles.push_back(v);
                
                if(_debug) cout<<"Add triangle:"<<v[0]<<" "<<v[1]<<" "<<v[2]<<'\n';
                
            }
            spoint.pop();
            spoint.push(topStackPoint);
            spoint.push(topQueuePoint);
        }
        else
        {
            while( spoint.size() > 1 )
            {
                Pointbase stack1Point=spoint.top();
                spoint.pop();
                Pointbase stack2Point=spoint.top();
                spoint.push(stack1Point);
                double pa[2], pb[2], pc[2];
                pa[0]=topQueuePoint.x; pa[1]=topQueuePoint.y;
                pb[0]=stack2Point.x;   pb[1]=stack2Point.y;
                pc[0]=stack1Point.x;   pc[1]=stack1Point.y;
                
                if(_debug)
                {
                    cout<<"current top queue vertex index="<<topQueuePoint.id<<'\n';
                    cout<<"Current top stack vertex index="<<stack1Point.id<<'\n';
                    cout<<"Second stack vertex index="<<stack2Point.id<<'\n';
                }
                
                double area=orient2d(pa,pb,pc);
                bool   left=stack1Point.left;
                if( (area > 0 && left) || (area < 0 && !left ) )
                {
                    vector<unsigned int>  v(3);
                    v[0]=topQueuePoint.id;
                    v[1]=stack2Point.id;
                    v[2]=stack1Point.id;
                    _triangles.push_back(v);
                    if(_debug) cout<<"Add triangle:"<<v[0]<<" "<<v[1]<<" "<<v[2]<<'\n';
                    spoint.pop();
                } else break;
            }
            
            spoint.push(topQueuePoint);
            
        }
        
        qvertex.pop();
        
    }
    
    Pointbase lastQueuePoint=qvertex.top();
    while( spoint.size() !=1 )
    {
        Pointbase topPoint=spoint.top();
        spoint.pop();
        Pointbase top2Point=spoint.top();
        
        vector<unsigned int>  v(3);
        v[0]=lastQueuePoint.id;
        v[1]=topPoint.id;
        v[2]=top2Point.id;
        _triangles.push_back(v);
        
        if(_debug) cout<<"Add triangle:"<<v[0]<<" "<<v[1]<<" "<<v[2]<<'\n';
    }
}

//----------------------------------------------------------------------------
//main triangulation function;
////----------------------------------------------------------------------------
void Polygon::triangulation()
{
    
    struct timeval tv0, tv1;
    struct timezone tz;
    gettimeofday(&tv0, &tz);
    partition2Monotone();
    searchMonotones();
    list<list<unsigned int> >::iterator it=_mpolys.begin();
    for(; it!=_mpolys.end(); it++)
        triangulateMonotone(*it);
    gettimeofday(&tv1, &tz);
    CheckTime(tv0, tv1);
//    cout<<"Total number of triangles:"<<_triangles.size()<<'\n';
    
}


void Polygon::saveAsShowme()
{

    cout <<(int)_triangles.size() + main_filename.size()<< " ";
    cout <<(int)main_filename.size() << endl;
    for(auto itt : main_filename){
        for(auto tmp : itt){
            cout << tmp.first<< " " << tmp.second << " ";
        }
        cout << endl;
    }
    auto itt=_triangles.begin();
    for(; itt!=_triangles.end(); itt++){
        cout<<_points[(*itt)[0]] -> x << " " << _points[(*itt)[0]] -> y << " ";
        cout<<_points[(*itt)[1]] -> x << " " << _points[(*itt)[1]] -> y << " ";
        cout<<_points[(*itt)[2]] -> x << " " << _points[(*itt)[2]] -> y << endl;
    }
    
}

int main()
{
    freopen("/Users/stetsyk/Desktop/Triangulation/input.txt", "r", stdin);
    freopen("/Users/stetsyk/Desktop/Triangulation/output.txt", "w" , stdout);
    vector<vector<pair<double, double > > > a;
    int m;
    cin >> m;
    for(int iterator = 0; iterator < m; ++iterator){
        int n;
        cin >> n;
        vector<pair<double, double > > d;
        for(int i = 0; i < n; ++i){
            double x, y;
            cin >> x >> y;
            d.push_back({x, y});
        }
        a.push_back(d);
    }
    Polygon poly(a);
	poly.triangulation();                    //main triangulation function
    poly.saveAsShowme();
    
	return 0;
}






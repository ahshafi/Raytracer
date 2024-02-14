#ifndef DEBUG_H
#define DEBUG_H
void err(istream_iterator<string> it) {cout<<endl;}
 template<typename T, typename... Args>
 void err(istream_iterator<string> it, T a, Args... args){
  cout << *it << " = " << a << ", ";err(++it, args...);
 }
 template<class T1, class T2>
 ostream &operator <<(ostream &os, pair<T1,T2>&p) {
  os<<"{"<<p.first<<", "<<p.second<<"} ";
  return os;
 }
 #define dbg1(args...) { string _s = #args; replace(_s.begin(), _s.end(), ',', ' '); stringstream _ss(_s); istream_iterator<string> _it(_ss); cout<<"line "<<__LINE__<<": "; err(_it, args); }

#define Gene template< class
#define Rics printer& operator,
Gene c> struct rge{c b, e;};
Gene c> rge<c> range(c i, c j){ return {i, j};}
struct printer{
    ~printer(){cerr<<endl;}
    Gene c >Rics(c x){ cerr<<boolalpha<<x; return *this;}
    Rics(string x){cerr<<x;return *this;}
    Gene c, class d >Rics(pair<c, d> x){ return *this,"(",x.first,", ",x.second,")";}
    Gene ... d, Gene ...> class c >Rics(c<d...> x){ return *this, range(begin(x), end(x));}
    Gene c >Rics(rge<c> x){
        *this,"["; for(auto it = x.b; it != x.e; ++it)
            *this,(it==x.b?"":", "),*it; return *this,"]";}
};

#define debug() cerr<<"LINE "<<__LINE__<<" >> ", printer()
#define dbg2(x) debug(), "[",#x,": ",(x),"] "
#define stop getchar()
#endif
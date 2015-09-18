#ifndef VECX
#define VECX
// Vector (Math) class
#include <iostream>
#include <vector>
#include <map>
#include <initializer_list>
#include <algorithm>
#include <math.h>

using namespace std;
//class Vec3;
typedef unsigned long long uint_16;  
typedef unsigned long uint_8; 
typedef unsigned int uint_4; 
typedef unsigned short uint_2;
typedef long long int_16; 
typedef long int_8; 
typedef int int_4; 
typedef short int_2; 


template <typename Key, typename T>
ostream &operator<<(ostream &out, const map<Key, T> &a){
  typename map<Key,T>::const_iterator it;
  for ( it=a.begin() ; it != a.end(); it++ )
    cout << (*it).first << ":" << (*it).second << " ";
  return out;
}


template <typename T> class VecX {
 public: 
  //----------------------
  // Properties
  //----------------------
  typedef initializer_list<int_8> ilist1; 
  typedef initializer_list<T> ilist2; 
  map<int_8, T> px; 
  vector<T> data;  
  int_8 rank; 
  
  typedef typename vector<T>::iterator it_vec; 
  typedef typename vector<T>::const_iterator it_cvec;
  typedef typename map<int_8, T>::iterator it_map;   
  typedef typename map<int_8, T>::const_iterator it_cmap;   
  typedef typename ilist1::iterator it_list1; 
  typedef typename ilist1::const_iterator it_clist1; 
  typedef typename ilist2::iterator it_list2; 
  typedef typename ilist2::const_iterator it_clist2; 

  //----------------------
  // Constructors
  //----------------------
  // C1: Empty
  VecX() {rank=0;}
  // C2: (long int) - Sets rank and contains - uncompressed empty
  VecX(int_8 const &a) {rank=a; data.assign(a, T(0));}
  // C3: (int, T*) - Sets rank and contains - uncompressed
  VecX(int_8 const &a, T const *b) {
    rank = a; 
    //data.resize(rank); 
    data.assign(b, b+rank);
  }; 
  VecX(int_8 const &a, T b) {
    rank = a; 
    //    data.resize(a); 
    data.assign(a, b); 
  }

  // C4: ({list}) - Assigns list values - uncompressed
  VecX(ilist2 a) {
    rank = a.size(); 
    data.assign(a.begin(), a.end()); 
   };
  // C5: ({int list}, {list}) - index-value couples - compressed
  VecX(ilist1 b, ilist2 c, int_8 const &a = 0) {
    if (b.size() != c.size()) {
      cout<< "Init VecX Error: list sizes are not the same"<<endl; 
      exit(1); 
    }    
    it_list2 xit=c.begin(); 
    for(it_list1 it=b.begin(); it!=b.end(); ++it, ++xit) 
      px[*it] = *xit; 
    data.clear(); 
    rank = *max_element(b.begin(), b.end()) + 1; 
    rank = (a > rank) ? a : rank; 
   };

  //----------------------
  // Methods
  //----------------------
  int_8 size() {if (px.empty()) return data.size(); else return rank; }
  void reserve(int_8 const &a) {if (!px.empty()) data.reserve(a);}
  void resize(int_8 const &a) {
    //rank = (rank<a)?a:rank; 
    rank = a; 
    if (px.empty()) 
      data.resize(rank);
  }
  void assign(int_8 const &a, T b) { 
    if (px.empty())
      data.assign(a, b); 
  }
  void assign(VecX<T> a) {
    if (px.empty() && a.px.empty()) {
      data.assign(a.data.begin(), a.data.end());
      rank = a.data.size(); 
    } else if (!a.px.empty()) {
      cout << "DATA is compressed!" <<endl; 
    }
    rank = a.rank; 
  }

  T min() { 
    if (px.empty()) 
      return *min_element(data.begin(), data.end()); 
    else
      return 0; //*min_element(px.begin(), px.end()); 
  }

  T max() { 
    if (px.empty()) 
      return *max_element(data.begin(), data.end()); 
    else
      return 0;//*max_element(px.begin(), px.end()); 
  }

  void empty() {px.empty(); data.empty();}
  void insert(ilist2 const &a, int_8 const &b=0) {
    if (px.empty()) {
      if (data.size() < b + a.size()) resize(b + a.size()); 
      it_vec dit=data.begin()+b; 
      for (it_list2 it=a.begin(); it!=a.end(); ++it, ++dit)
	*dit = *it; 
    } else {
      int_8 i=b; 
      for(it_list2 it=a.begin(); it!=a.end(); ++it,++i) 
	if (*it != 0) px[i] = *it; 
      rank = rank<b+a.size() ? b+a.size() : rank; 
    }
  }
  // Removes zeros and converts the data structure into 
  // index-> value pairs
  void compress() {
    if (px.empty()) {
      for (int_8 i = 0; i<data.size(); i++) {
	if (data[i] == 0) continue;
      	px[i]=data[i]; 
      }
    }
    if (!data.empty()) data.clear(); 
  }
  void uncompress() {
    data.resize(rank); 
    for (it_cmap it=px.begin(); it!=px.end(); ++it) 
      data[it->first] = it->second; 
    px.clear();
  }
  // Push_pair(int, T) : add index->value to the list 
  // priority is on compressed if empty
  void push_pair(int_8 const &a, T const &b) {    
    if (!data.empty()) {
      if (a+1>data.size()) {
	data.resize(a+1); 
	rank = a+1;
      }
      data[a] = b; 
    } else {
      if (b!=0) px[a] = b; 
    } 
  }
  // Push_pair
  void push_pair(vector<int_8> const a, vector<T> const b) {
    for(vector<int_8>::size_type i=0; i<a.size(); ++i)  
      push_pair(a[i], b[i]); 
  }
  void push_pair(ilist1 const a, ilist2 const b) {
    it_list2 ix=b.begin(); 
    for(it_list1 it=a.begin(); it!=a.end(); ++it,++ix)  
      push_pair(*it, *ix); 
  }
  void push_back(T b) {
    if (px.empty()) {
      data.push_back(b); ++rank;  
    } else
      cout << "Uncompress the data before pushing back a new value"<<endl; 
  }

  void clear() {data.clear(); px.clear(); }
  double abs() {
    double sum = 0; 
    if (!px.empty()) {
      for (it_cmap it=px.begin(); it!=px.end(); ++it) 
	sum += (it->second)*(it->second); 
    } else {
      for (int_8 i = 0; i<rank; ++i) 
	sum += data[i]*data[i]; 
    }
    return sqrt(sum);       
  }

  VecX<double> comp(int k=0); 
  // Get value at a linear index
  // T get(uint_8 a) {
  //   if (a > rank) {
  //     cout << "Error: index greater than the rank" << endl; 
  //     exit(1); 
  //   }
  //   if (!px.empty()) return px[a]; // if out of bounds, map returns zero; 
  //   else if (!data.empty() && a < data.size()) return data[a];
  //   else {
  //     cout << "Error: not a good request!"<<endl; 
  //     exit(1); 
  //   }
  // }

  //----------------------
  // Operators
  //----------------------
  // O1: Assignment ! WITH VECTOR
  VecX<T>& operator =(VecX<T> const &a) {
    if (this != &a) {
      this->data.assign(a.data.begin(), a.data.end()); 
      // std::copy(a.data.begin(), a.data.end(), data.begin()); 
      this->px = a.px; 
      this->rank = a.rank; 
    }
    return *this; 
  }					      
  // O2: Assignment | WITH Init-List
  VecX<T>& operator =(const ilist2 a){
    rank = (a.size() > rank) ? a.size() : rank;
    if (!px.empty()) {
      int_8 i=0; 
      for (it_clist2 it=a.begin(); it!=a.end(); ++it, ++i) 
	px[i] = *it; 
    } else if (!data.empty()) 
      data.assign(a.begin(), a.end());
    return *this; 
  };
  // O3: Bracket    | ACCESS to either data or px element
  //                  Return a new position if element don't exists
  T& operator[] (int_8 const &index) {
    if (!px.empty()) {
      if (index > rank) rank = index+1; 
      return px[index];       
    } else if (!data.empty()) {
      if (index > rank) {
	rank = index+1; 
	data.resize(index+1); 
      } 
      return data[index]; 
    } else {
      return px[index]; 
    }    
  };
  // const T& operator[] (uint_8 const &index) const {
  //    T result = 0; 
  //    if (!px.empty()) {
  //      result = px[index];       
  // } else if (!data.empty()) {
  //      result = data[index];
  //    } 
  //    return result;    
  // };
  // O4: Multiplication | DOT PRODUCT
  friend T operator*(VecX<T> const &a, VecX<T> const &b) {
    // cout << "&&&&" << endl; 
    // a.info(); b.info(); 
    if (a.rank != b.rank) {
      cout<<"Error: Cannot multiply two vectors having different ranks: "<<a.rank<< " " <<b.rank<<endl;
      exit(1); 
    }
    T sum = 0; 
    if (!a.px.empty()) {
      if (!b.px.empty()) {
	for (it_cmap it=a.px.begin(); it!=a.px.end(); ++it) {
	  it_cmap bt = b.px.find(it->first);
	  sum += it->second * bt->second; 
	} 
      } else {
	for (it_cmap it=a.px.begin(); it!=a.px.end(); ++it) {
	  sum += it->second * b.data[it->first]; 
	}
      }
    } else if (!b.px.empty()) {
      for (it_cmap it=b.px.begin(); it!=b.px.end(); ++it) 
	sum += it->second * a.data[it->first]; 
    } else {
      for (int i=0; i<a.rank; ++i) 
      	sum += a.data[i]*b.data[i]; 
    }
    return sum;
  };
  // 05a: Multiplication : WITH a real number
  // friend VecX<T> operator*(VecX<T> const &a, T const &b) {
  //   VecX<T> result(a.rank); 
  //   if (b == 0) return result;
  //   if (!a.px.empty()) {
  //     for (it_cmap it=a.px.begin(); it!=a.px.end(); ++it)
  // 	result.px[it->first] = it->second*b;
  //   } else if (!a.data.empty()) {
  //     result.resize(a.rank);
  //     for (int i=0; i<a.rank; i++)
  // 	result.data[i] = a.data[i]*b; 
  //   }
  //   return result;
  // }; 
  // // 05b: Multipication : Variant of A
  // friend VecX<T> operator*(T const &b, VecX<T> const &a) {
  //   return a*b; 
  // };

  // 05a: Multiplication : WITH a real number
  friend VecX<T> operator*(VecX<T> const &a, double const &b) {
    VecX<T> result(a.rank); result.compress();  
    if (b == 0) return result;
    if (!a.px.empty()) {
      for (it_cmap it=a.px.begin(); it!=a.px.end(); ++it)
	result.px[it->first] = it->second*b;
    } else if (!a.data.empty()) {
      result.resize(a.rank);
      for (int i=0; i<a.rank; i++)
	result.data[i] = a.data[i]*b; 
    }
    return result;
  }; 
  // 05b: Multipication : Variant of A
  friend VecX<T> operator*(double const &b, VecX<T> const &a) {
    return a*b; 
  };

  VecX<T>& operator+=(double const &a) {
    if (a == 0.0) return *this; 
    if (!px.empty()) this->uncompress();
    if (data.size() < rank) data.resize(rank); 
    for (it_vec it=data.begin(); it!=data.end(); ++it)
      //data[i] = data[i]+a; 
      *it = *it + a; 
    return *this; 
  };
  VecX<T>& operator=(double const &a) { this->clear(); *this += a; return *this; }

  VecX<T>& operator+=(VecX<T> const &a) {
    if (a.data.empty() && a.px.empty()) return *this;
    rank = rank < a.rank ? a.rank : rank;  
     if (!px.empty() && !a.px.empty()) {
       for (it_cmap it=a.px.begin(); it!=a.px.end(); ++it) {
	 if (px[it->first] + it->second == 0) continue;
	 px[it->first] = px[it->first] + it->second;
       }
    } else if(!px.empty()) {
      data = a.data; 
      for (it_cmap it=px.begin(); it!=px.end(); ++it) 
	data[it->first] = data[it->first] + it->second;
      px.clear(); 
    } else if(!a.px.empty()) {
      if (!data.empty()) {
	for (it_cmap it=a.px.begin(); it!=a.px.end(); ++it)
	  data[it->first] = data[it->first] + it->second; 
      } else {
	px = a.px; 
      }	
    } else {
      if (!data.empty()) {
	for (int i=0; i<rank; ++i)
	  data[i] = data[i] + a.data[i];      
      } else {
	data = a.data; 
      }
    }
    return *this; 
  }

  VecX<T>& operator/=(VecX<T> const &a) {
    if (a.data.empty() && a.px.empty()) return *this; 
    if (a.rank != rank) return *this; 
    // a should not be compressed --> leads to div by zero
    if (a.data.empty()) { cout << "Error: div by zero! compressed"<< endl; exit(1);}
    if (!px.empty()) {
      for (it_cmap it=px.begin(); it!=px.end(); ++it) {
	if (a.data[it->first] == 0) {cout<< "Error: div by zero @"<<it->first<<endl; exit(1);}
	px[it->first] /= a.data[it->first];
      }
    } else {
      if (!data.empty()) 
	for (int i=0; i<rank; ++i) {
	  if (a.data[i] == 0) {cout<< "Error: div by zero @"<<i<<endl; exit(1);}
	  data[i] /= a.data[i];     
	} 
    }
    return *this; 
  }
  VecX<T>& operator/=(T const &a) {
    if (a == 0) {cout << "Error: div by zero! constant"<<endl; exit(1);}
    if (!px.empty()) {
     for (it_cmap it=px.begin(); it!=px.end(); ++it) 
	px[it->first] /= a;
    } else {
      if (!data.empty()) 
	for (int i=0; i<rank; ++i)
	  data[i] /= a;      
    }
    return *this; 
  }
   

  VecX<T> operator-() const {VecX<T> c; c = -1.0*(*this); return c; };

  VecX<T>& operator-=(double const &a) {*this += -a; return *this;}
  VecX<T>& operator-=(VecX<T> const &a) {*this += -a; return *this; }

  //friend operator Vec3 () {Vec3 a(*this); return a; }; 

  // 06: Addition : VecX and double
  friend VecX<T> operator+(VecX<T> const &a, VecX<T> const &b) {VecX<T> c=a; c += b; return c;}
  friend VecX<T> operator+(VecX<T> const &a, T const &b) {VecX<T> c=a; c += b; return c;};
  friend VecX<T> operator+(T const &b, VecX<T> const &a) {VecX<T> c=a; c += b; return c;};
  //friend VecX<T> operator-(VecX<T> const &a, VecX<T> b) {return -b + a;};
  //friend VecX<T> operator-(VecX<T> const &a, T b) {return -b + a;};
  friend VecX<T> operator-(VecX<T> const &a, VecX<T> const &b) {VecX<T> c=a; c -= b; return c;}
  friend VecX<T> operator-(VecX<T> const &a, T const &b) {VecX<T> c=a; c -= b; return c;};
  friend VecX<T> operator-(T const &b, VecX<T> const &a) {VecX<T> c = -a; c += b; return c;};
  
  friend VecX<T> operator/(VecX<T> const &a, VecX<T> const &b) {VecX<T> c=a; c /= b; return c;}
  friend VecX<T> operator/(VecX<T> const &a, T const &b) {VecX<T> c=a; c /= b; return c;}
  friend VecX<T> operator/(T const &b, VecX<T> const &a) {VecX<T> c=a; c /= b; return c;}

  //----------------------
  // Output options
  //----------------------
  friend ostream &operator<<(ostream &out, const VecX<T> &a){
    if (!a.px.empty()) {
      out << (a.px); 
    } else {
      for (typename vector<T>::size_type i = 0; i<a.data.size(); i++)
	out<< a.data[i] << " " ; 
    }
    return out;
  }
  void info() {
    cout<<endl<<"--------------------------------"<<endl; 
    cout<<"      Rank:\t"<<rank<<endl; 
    cout<<"  Contains:\t"<<(!px.empty()?px.size():data.size())<<endl; 
    cout<<"Size(1 el):\t"<<sizeof(T)<<" bytes"<<endl;
    cout<<"Total size:\t~"<<data.size()*sizeof(typename vector<T>::value_type) + px.size()*sizeof(typename map<int_8,T>::value_type)<<" bytes"<<endl;
    cout<<"Compressed:\t"<<(!px.empty()?"Yes":"No")<<endl; 
    if (!px.empty()) cout<<"     Dense:\t"<<(double)px.size()/rank*100<<"%"<<endl; 
    cout<<"--------------------------------"<<endl; 
  }
  it_cvec begin() { if (!px.empty()) return data.end(); else return data.begin(); } 
  it_cvec end() { return data.end(); }
};

 class Vec3: public VecX<double> {
 public:
   Vec3 () {resize(3);}
   Vec3 (double const &a) {resize(3); data[0]=a; };
   Vec3 (double const & a, double const & b) {resize(3); data[0]=a; data[1]=b;};
   Vec3 (double const & a, double const & b, double const & c) {
     resize(3); data[0]=a; data[1]=b; data[2]=c;  
   };
   Vec3 (Vec3 const &a) {    
     data.assign(a.data.begin(), a.data.end()); 
     resize(3); 
   }
   Vec3 (VecX<double> const &a) {data.assign(a.data.begin(), a.data.begin() + std::min(3, (int)a.data.size())); resize(3);}
   Vec3 (initializer_list<double> a) {     
     data.assign(a.begin(), a.end()); 
     resize(3); 
   };
   
   void set(Vec3 const &a) {
     data.assign(a.data.begin(), a.data.end()); 
     //resize(3); 
   }
   
   double x() {return data[0];}
   double y() {return data[1];}
   double z() {return data[2];}

   friend bool operator==(Vec3 const &a, Vec3 const &b) {
     bool isit = true;
     if (a.rank != b.rank) {
       isit = false; 
     } else {
       int_8 icnt =0; auto jt = b.data.begin();
       for (auto it = a.data.begin(); it != a.data.end(); ++it, ++icnt) {
	 jt = b.data.begin() + icnt; 
	 if (*jt == *it) {continue;}
	 isit = false; 
	 break; 
       }
     }
     return isit;
   }
   
};

// // + needs to be outside function
// inline Vec3D operator+(Vec3D lhs, Vec3D const &rhs ) { lhs += rhs; return lhs; };
// inline Vec3D operator+(Vec3D lhs, double const &rhs ) { lhs += rhs; return lhs; };
// inline Vec3D operator+(double const &rhs,Vec3D lhs ) { lhs += rhs; return lhs; };
// inline Vec3D operator-(Vec3D lhs, Vec3D const &rhs ) { lhs -= rhs; return lhs; };
// inline Vec3D operator-(Vec3D lhs, double const &rhs ) { lhs -= rhs; return lhs; }
// inline Vec3D operator-(double const &rhs,Vec3D lhs ) { lhs -= rhs; return lhs; };
// inline double operator*(Vec3D const &lhs, Vec3D const &rhs ) {return lhs.x*rhs.x + lhs.y*rhs.y + lhs.z*rhs.z; };


inline Vec3 operator^(Vec3 const& lhs, Vec3 const& rhs ) {
  Vec3 a;
  a.data[0] = (lhs.data[1]*rhs.data[2] - lhs.data[2]*rhs.data[1]);
  a.data[1] = (lhs.data[2]*rhs.data[0] - lhs.data[0]*rhs.data[2]);
  a.data[2] = (lhs.data[0]*rhs.data[1] - lhs.data[1]*rhs.data[0]);
  return a;
};
inline double abs(Vec3 const& a) {return sqrt(a*a);};

template <class Vec3>
VecX<double> VecX<Vec3>::comp(int k) {
  VecX<double> result(rank); 
  if (px.empty()) result.uncompress(); 
  for (int i=0; i<rank; ++i) { 
    result[i] = data[i][k]; 
  }
  return result; 
}


#endif

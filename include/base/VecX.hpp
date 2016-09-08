/*
 ***************************************************************************
 *   Copyright (C) 2015, Eray Uzgoren                                      *
 *                                                                         *
 *   This program is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation, either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>. *
 ***************************************************************************
*/
#ifndef VECX
#define VECX

// template <typename T>
// ostream &operator<<(ostream &out, const vector<T> &a){
//   for (auto it=a.begin() ; it != a.end(); it++ )
//     cout << (*it) << " ";
//   return out;
// }


template <typename T> class VecX {
 public: 
  //----------------------
  // Properties
  //----------------------
  typedef initializer_list<int_8> ilist_i8;
  typedef initializer_list<T> ilist_T; 
  vector<T> data;  
  
  //----------------------
  // Constructors
  //----------------------
  // C1: Empty
  VecX() {};
  // C2: (long int) - Sets rank and contains - uncompressed empty
  VecX(int_8 const &a) {if ( a > 0 ) data.assign(a, T(0));}  
  // C3: (int, T*) - Sets rank and contains - uncompressed
  //  VecX(int_8 const &a, T const *b) {data.assign(b, b+a);}; 
  VecX(int_8 const &a, T b) {data.assign(a, b);}

  // C4: ({list}) - Assigns list values - uncompressed
  VecX(ilist_T a) {data.assign(a.begin(), a.end());};

  typename vector<T>::iterator begin() { return data.begin(); }
  typename vector<T>::iterator end()   { return data.end(); }
  
  //----------------------
  // Methods
  //----------------------
  int_8 size() {return data.size();}
  void reserve(int_8 const &a) {data.reserve(a);}
  int_8 capacity() {return data.capacity();}
  void resize(int_8 const &a) {data.resize(a);}
  void assign(int_8 const &a, T b=T(0)) {vector<T>(a,T(b)).swap(data);}
  void assign(VecX<T> &a) {data.assign(a.data.begin(), a.data.end());}
  void assign(ilist_T const &a, int_8 const &b=0) {
    if (data.size() < b + a.size()) resize(b + a.size());
    auto it = a.begin();
    for (auto i=0; i < a.size(); ++i, ++it) 
      data[b+i] = *it;   
  }
  void swap(VecX<T> &a) {a.data.swap(data);}
  
  T min() {return *min_element(data.begin(), data.end());}
  T max() {return *max_element(data.begin(), data.end());}

  void empty() {std::vector<T>().swap(data);} // capacity is zero (no delete)
  void clear() {data.clear();} // capacity is unchanged

  void zeros(int_8 n) {vector<T>(n, T(0)).swap(data);}
  void ones(int_8 n) {vector<T>(n, T(1)).swap(data);}
  void linear(double const &x0, double const &x1, int_8 const &n) {
    auto dx = (x1-x0)/(double)(n-1); data.resize(n); 
    for (auto i = 0; i < n; ++i)
      data[i] = x0 + i*dx;
  }
  void linear(double const &x0, double const &dx, double const &x1) {
    data.clear(); 
    for (auto x = x0; x < x1 || almost_equal(x, x1, 10); x+=dx)
      data.emplace_back(x);
  }
  
  // Push_pair(int, T) : add index->value to the list 
  // priority is on compressed if empty
  void push_pair(int_8 const &a, T const &b) {
    if (a >= data.size()) data.resize(a+1);     
    data[a] = b;  
  }
  // Push_pair
  void push_pair(vector<int_8> const &a, vector<T> const &b) {
    for(auto i=0; i < a.size(); ++i)  
      push_pair(a[i], b[i]); 
  }
  void push_pair(ilist_i8 const a, ilist_T const b) {
    push_pair(vector<int_8>(a), vector<T>(b)); 
  }  
  void push_back(T const &b) {data.push_back(b);}

  double abs() {
    double sum = 0; 
    for (auto i = 0; i<size(); ++i) 
      sum += data[i]*data[i]; 
    return sqrt(sum);       
  }
  
  VecX<double> comp(int k);

  //----------------------
  // Operators
  //----------------------
  // O1: Assignment ! WITH VECTOR
  VecX<T>& operator =(VecX<T> const &a) {
    if (this != &a) {
      this->data.assign(a.data.begin(), a.data.end()); 
    }
    return *this; 
  }					      
  // O2: Assignment | WITH Init-List
  VecX<T>& operator =(const ilist_T a){
    data.assign(a.begin(), a.end());
    return *this; 
  };

  // O3: Bracket    | ACCESS to either data or px element
  //                  Return a new position if element don't exists
  T& operator[] (int_8 const &index) {
    if (index >= data.size()) {
      cout << "Error: VecX out of bounds! "<< endl;
      exit(1);
    }
    return data[index]; 
  };

  // O4: Multiplication | DOT PRODUCT
  friend T operator*(VecX<T> const &a, VecX<T> const &b) {
    if (a.data.size() != b.data.size()) {
      cout<<"Error: sizes are not equal!!! A: "<<a.data.size() << " B: " <<b.data.size()<<endl;
      exit(1); 
    }
    T sum = 0; 
    for (auto i=0; i < a.data.size(); ++i) 
      sum += a.data[i]*b.data[i]; 
    
    return sum;
  };

  // 05a: Multiplication : WITH a real number
  friend VecX<T> operator*(VecX<T> const &a, double const &b) {
    VecX<T> result(a.data.size()); 
    if (b == 0) return result;
    for (auto i=0; i < a.data.size(); i++)
      result[i] = a.data[i]*b; 
    return result;
  };
  
  // 05b: Multipication : Variant of A
  friend VecX<T> operator*(double const &b, VecX<T> const &a) {
    return a*b; 
  };

  VecX<T>& operator+=(double const &a) {
    if (a == 0.0) return *this; 
    for (auto it=data.begin(); it!=data.end(); ++it)
      *it = *it + a; 
    return *this; 
  };
  
  VecX<T>& operator=(double const &a) {
    assign(data.size(), T(a)); 
  }

  VecX<T>& operator+=(VecX<T> const &a) {
    for (auto i=0; i < size(); ++i)
      data[i] = data[i] + a.data[i];      
    return *this; 
  }

  VecX<T>& operator/=(VecX<T> const &a) {
    for (int i=0; i<size(); ++i) {
      if (a.data[i] == 0) {cout<< "Error: div by zero @"<<i<<endl; exit(1);}
      data[i] /= a.data[i];     
    }
    return *this; 
  }
  
  VecX<T>& operator/=(T const &a) {
    if (a == 0) {cout << "Error: div by zero! constant"<<endl; exit(1);}
    for (int i=0; i<size(); ++i)
      data[i] /= a;      
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
    out << " [";
    for (auto i = 0; i<a.data.size(); i++)
      out<< a.data[i] << " " ;
    out << "] " ; 
    return out;
  }  
};

template <typename T>
VecX<T> abs(VecX<T> const &a) {
  VecX<T> result;
  for (auto el : a.data) 
    result.data.push_back(abs(el));
  return result; 
}



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
     if (a.data.size() != b.data.size()) {
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
  VecX<double> result(size()); 
  for (auto i=0; i<size(); ++i) { 
    result[i] = data[i][k]; 
  }
  return result; 
}


#endif
